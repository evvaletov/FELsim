#!/usr/bin/env python3
"""C3: FR3+MGE Stage 11-only CMA-ES optimization.

Stages 1-10 are fixed to their FIT-converged values (from the FR3→FR3+MGE
chain run). Only Stage 11's 4 variables (Ic, I, I2, I3) are optimized
externally via CMA-ES, with each evaluation running a single COSY lattice
propagation (~6-8s).

Author: Eremey Valetov
"""

import os
import sys
import re
import json
import math
import subprocess
import time
import argparse
import numpy as np

backend_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, backend_dir)

# ── Beam parameters (must match UHM_beamline_opt_cosy.py) ───────────────────
ENERGY = 40
EPSILON_N = 8
X_STD = 0.8

E0 = 0.511
gamma = (ENERGY + E0) / E0
beta_rel = np.sqrt(1 - 1 / gamma**2)
norm = gamma * beta_rel
EPSILON = EPSILON_N / norm

K = 1.2
LAMBDA_U = 2.3e-2
BETA_YM = gamma * LAMBDA_U / (2 * np.pi * K)
ALPHA_YM = 0.0
BETA_XM = 1.4
ALPHA_XM = 0.47
BETA_0 = X_STD**2 / EPSILON


def generate_fox_template(fox_path):
    """Convert a FIT-enabled FOX file into a no-FIT evaluation template.

    Strips all FIT...ENDFIT blocks, post-FIT verification sections, and
    redundant LATTICE calls. Replaces Stage 11 variable assignments with
    placeholders. Moves the LATTICE call after all variable assignments.
    """
    with open(fox_path) as f:
        lines = f.readlines()

    output = []
    in_fit = False
    in_postfit = False

    for line in lines:
        stripped = line.strip()

        # Skip FIT...ENDFIT blocks (including internal LATTICE and OBJ lines)
        if re.match(r'\s*FIT\s+S\d', stripped):
            in_fit = True
            continue
        if in_fit:
            if re.match(r'\s*ENDFIT\s', stripped):
                in_fit = False
                in_postfit = True
            continue

        # Skip post-FIT verification sections
        if in_postfit:
            if 'POST-FIT' in stripped or '{ Post-FIT verification }' in stripped:
                continue
            if (stripped.startswith('LATTICE') or
                stripped.startswith('LOOP IDX') or
                stripped.startswith('B0(') or
                stripped.startswith('A0(') or
                stripped.startswith('G0(') or
                stripped.startswith('WRITE 6')):
                continue
            in_postfit = False

        # Skip standalone POST-FIT lines
        if 'POST-FIT' in stripped:
            continue

        # Skip ALL standalone LATTICE calls — we'll add one at the right place
        if stripped == 'LATTICE ;':
            continue

        # Replace Stage 11 variable assignments with placeholders
        m = re.match(r'S11_(Ic|I2|I3|I)\s*:=\s*[\d.e+-]+\s*;', stripped)
        if m:
            indent = line[:len(line) - len(line.lstrip())]
            var = 'S11_' + m.group(1)
            output.append(f"{indent}{var} := ___{var}___ ;\n")
            continue

        # Insert LATTICE call just before the PM command (after all assignments)
        if stripped.startswith('CO 1 ; PM 99'):
            output.append("    LATTICE ;\n")

        output.append(line)

    return ''.join(output)


def create_evaluation_template(fox_path):
    """Create a minimal FOX template for single-evaluation CMA-ES."""
    template = generate_fox_template(fox_path)

    # Verify placeholders were inserted
    for var in ['S11_Ic', 'S11_I', 'S11_I2', 'S11_I3']:
        placeholder = f'___{var}___'
        if placeholder not in template:
            raise ValueError(f"Placeholder {placeholder} not found in template")

    # Verify no FIT blocks remain
    if re.search(r'\bFIT\b.*\bENDFIT\b', template, re.DOTALL):
        raise ValueError("FIT blocks still present in template")

    return template


def read_transfer_map(results_dir):
    """Read 6x6 transfer map from COSY fort.99."""
    fort99 = os.path.join(results_dir, 'fort.99')
    M = np.eye(6)

    num_re = re.compile(r'[+-]?\d+\.?\d*(?:[Ee][+-]?\d+)?')
    index_to_col = {
        '100000': 0, '010000': 1, '001000': 2,
        '000100': 3, '000010': 4, '000001': 5,
    }

    with open(fort99) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('--'):
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            index = parts[-1]
            if index not in index_to_col:
                continue
            col = index_to_col[index]
            data_part = stripped[:stripped.rfind(index)]
            values = [float(x) for x in num_re.findall(data_part)]
            for row in range(min(5, len(values))):
                M[row, col] = values[row]

    return M


def compute_objective(M):
    """Stability-aware Twiss MSE from transfer map."""
    cos_mu_x = (M[0, 0] + M[1, 1]) / 2
    cos_mu_y = (M[2, 2] + M[3, 3]) / 2
    instability = max(abs(cos_mu_x), abs(cos_mu_y))

    if instability > 1:
        # Smooth sigmoid transition instead of hard step
        excess = instability - 1
        penalty = 1e3 * (1 + np.log1p(excess))
        return penalty

    b0, g0 = BETA_0, 1.0 / BETA_0
    bx = M[0, 0]**2 * b0 + M[0, 1]**2 * g0
    ax = -(M[0, 0] * M[1, 0] * b0 + M[0, 1] * M[1, 1] * g0)
    by = M[2, 2]**2 * b0 + M[2, 3]**2 * g0
    ay = -(M[2, 2] * M[3, 2] * b0 + M[2, 3] * M[3, 3] * g0)

    mse = ((bx - BETA_XM)**2 + (by - BETA_YM)**2 +
           (ax - ALPHA_XM)**2 + (ay - ALPHA_YM)**2) / 4

    return mse if math.isfinite(mse) else 1e6


def evaluate(fox_template, cosy_bin, run_dir, s11_ic, s11_i, s11_i2, s11_i3):
    """Run one COSY evaluation with given Stage 11 currents."""
    fox = fox_template
    fox = fox.replace('___S11_Ic___', str(s11_ic))
    fox = fox.replace('___S11_I___', str(s11_i))
    fox = fox.replace('___S11_I2___', str(s11_i2))
    fox = fox.replace('___S11_I3___', str(s11_i3))

    with open(os.path.join(run_dir, 'input.fox'), 'w') as f:
        f.write(fox)

    try:
        result = subprocess.run(
            [cosy_bin, 'input.fox'], cwd=run_dir,
            capture_output=True, text=True, timeout=120,
        )
    except subprocess.TimeoutExpired:
        return 1e6

    if result.returncode != 0:
        return 1e6

    output = result.stdout + result.stderr
    fatal = ['COMMAND PLACEMENT', 'NOT DECLARED', 'ARRAY INDEX', 'VARIABLE EXHAUSTED']
    if any(m in output for m in fatal):
        return 1e6

    try:
        M = read_transfer_map(run_dir)
    except Exception:
        return 1e6

    return compute_objective(M)


def main():
    try:
        import cma
    except ImportError:
        sys.exit("ERROR: 'cma' package required. pip install cma")

    parser = argparse.ArgumentParser(description='C3: Stage 11-only CMA-ES for FR3+MGE')
    parser.add_argument('--fox-source', type=str, default=None,
                        help='Path to FIT-enabled FOX file (default: generate from adapter)')
    parser.add_argument('--sigma', type=float, default=0.5,
                        help='CMA-ES initial step size')
    parser.add_argument('--max-eval', type=int, default=2000,
                        help='Maximum function evaluations')
    parser.add_argument('--popsize', type=int, default=20,
                        help='Population size (0=CMA default)')
    parser.add_argument('--restarts', type=int, default=3,
                        help='Number of BIPOP restarts')
    parser.add_argument('--run-dir', type=str, default=None,
                        help='COSY run directory (default: backend/results)')
    parser.add_argument('--save', type=str,
                        default='results/cosy_s1_fr3_mge_cmaes_s11.json',
                        help='Output JSON file')
    args = parser.parse_args()

    # Generate FOX template
    if args.fox_source:
        fox_source = args.fox_source
    else:
        # Generate from the adapter
        from pathlib import Path
        from cosyAdapter import COSYAdapter
        from cosyOptHelper import add_stages, parse_beamline_felsim_indexed

        file_path = Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx'

        # Load Stage 1-10 converged currents from chain run
        chain_path = Path(__file__).resolve().parent / 'results' / 'cosy_s1_fr3_mge_warm_chain.json'
        with open(chain_path) as f:
            chain = json.load(f)

        # Import stage definitions
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from UHM_beamline_opt_cosy import build_stages, compute_targets, apply_warm_start

        targets = compute_targets()
        stages = build_stages(targets)
        apply_warm_start(stages, chain['currents'])

        config = {'simulation': {'KE': ENERGY, 'order': 3, 'dimensions': 3}}
        adapter = COSYAdapter(
            lattice_path=str(file_path), mode='transfer_matrix',
            config=config, fringe_field_order=3, use_mge_for_dipoles=True,
            debug=False
        )
        sim = adapter.get_native_simulator()
        sim.beamline = parse_beamline_felsim_indexed(str(file_path))[:118]
        sim.set_geometric_emittance(targets['epsilon'])
        sim.set_initial_twiss(beta_x=BETA_0, alpha_x=0.0, beta_y=BETA_0, alpha_y=0.0)
        sim.fit_nmax = 1000
        sim.fit_eps = 1e-8
        sim.fit_nalgorithm = 1
        sim.fit_combined_mse = True
        add_stages(sim, stages)

        fox_source = sim.generate_input(output_dir='results')
        print(f"Generated FOX source: {fox_source}")

    template = create_evaluation_template(fox_source)

    # Run directory setup
    if args.run_dir:
        run_dir = os.path.abspath(args.run_dir)
        os.makedirs(run_dir, exist_ok=True)
    else:
        run_dir = os.path.join(backend_dir, 'results')
    cosy_bin = os.path.join(backend_dir, 'results', 'cosy')
    if not os.path.isfile(cosy_bin):
        cosy_bin = 'cosy'

    # Starting point: FR3 converged values for Stage 11
    fr3_path = os.path.join(os.path.dirname(__file__), 'results', 'cosy_s1_fr3_warm.json')
    with open(fr3_path) as f:
        fr3 = json.load(f)
    x0 = [
        fr3['currents']['87'],  # S11_Ic
        fr3['currents']['93'],  # S11_I
        fr3['currents']['95'],  # S11_I2
        fr3['currents']['97'],  # S11_I3
    ]
    var_names = ['S11_Ic', 'S11_I', 'S11_I2', 'S11_I3']

    print(f"C3: Stage 11-only CMA-ES for FR3+MGE")
    print(f"  Variables: {var_names}")
    print(f"  Start: {[f'{v:.3f}' for v in x0]}")
    print(f"  Sigma: {args.sigma}, Max eval: {args.max_eval}, Popsize: {args.popsize}")
    print(f"  Targets: bx={BETA_XM}, ax={ALPHA_XM}, by={BETA_YM:.4f}, ay={ALPHA_YM}")

    # Test evaluation
    t0 = time.time()
    test_mse = evaluate(template, cosy_bin, run_dir, *x0)
    dt = time.time() - t0
    print(f"  Test eval: RMS={math.sqrt(test_mse):.4e} ({dt:.1f}s)")
    print(f"  Est. time: {dt * args.max_eval / 3600:.1f}h for {args.max_eval} evals")

    n_eval = [0]
    best_mse = [float('inf')]
    best_x = [None]
    t_start = time.time()

    def objective(x):
        mse = evaluate(template, cosy_bin, run_dir, x[0], x[1], x[2], x[3])
        n_eval[0] += 1
        if mse < best_mse[0]:
            best_mse[0] = mse
            best_x[0] = list(x)
            elapsed = time.time() - t_start
            stable = "STABLE" if mse < 1000 else "unstable"
            print(f"  [{n_eval[0]:5d}] RMS={math.sqrt(mse):.6e} ({stable}) "
                  f"Ic={x[0]:.3f} I={x[1]:.3f} I2={x[2]:.3f} I3={x[3]:.3f} "
                  f"t={elapsed:.0f}s")
        return mse

    opts = {
        'maxfevals': args.max_eval,
        'bounds': [[-10, -10, -10, -10], [10, 10, 10, 10]],
        'seed': 42,
        'verb_disp': 50,
        'verb_log': 0,
        'tolfun': 1e-10,
    }
    if args.popsize > 0:
        opts['popsize'] = args.popsize

    if args.restarts > 0:
        xopt, es = cma.fmin2(objective, x0, args.sigma, opts,
                             restarts=args.restarts, bipop=True)
    else:
        es = cma.CMAEvolutionStrategy(x0, args.sigma, opts)
        es.optimize(objective)
        xopt = es.result.xbest

    print(f"\nDone: {n_eval[0]} evaluations, best RMS={math.sqrt(best_mse[0]):.6e}")
    print(f"  S11_Ic = {best_x[0][0]:.6f}")
    print(f"  S11_I  = {best_x[0][1]:.6f}")
    print(f"  S11_I2 = {best_x[0][2]:.6f}")
    print(f"  S11_I3 = {best_x[0][3]:.6f}")

    # Compute final Twiss
    evaluate(template, cosy_bin, run_dir, *best_x[0])
    M = read_transfer_map(run_dir)
    b0, g0 = BETA_0, 1.0 / BETA_0
    bx = M[0, 0]**2 * b0 + M[0, 1]**2 * g0
    ax = -(M[0, 0] * M[1, 0] * b0 + M[0, 1] * M[1, 1] * g0)
    by = M[2, 2]**2 * b0 + M[2, 3]**2 * g0
    ay = -(M[2, 2] * M[3, 2] * b0 + M[2, 3] * M[3, 3] * g0)

    print(f"\nFinal Twiss at undulator:")
    print(f"  beta_x  = {bx:.6f}  (target {BETA_XM})")
    print(f"  beta_y  = {by:.6f}  (target {BETA_YM:.6f})")
    print(f"  alpha_x = {ax:.6f}  (target {ALPHA_XM})")
    print(f"  alpha_y = {ay:.6f}  (target {ALPHA_YM})")

    # Load chain data for full current set
    chain_path = os.path.join(os.path.dirname(__file__), 'results', 'cosy_s1_fr3_mge_warm_chain.json')
    with open(chain_path) as f:
        chain = json.load(f)

    # Merge Stage 11 optimized values into the chain currents
    currents = dict(chain['currents'])
    currents['87'] = best_x[0][0]
    currents['93'] = best_x[0][1]
    currents['95'] = best_x[0][2]
    currents['97'] = best_x[0][3]

    data = {
        'config': 'S1_2ps',
        'energy_MeV': ENERGY,
        'epsilon_n': EPSILON_N,
        'fringe_field_order': 3,
        'mge': True,
        'optimizer': 'cma-es-s11',
        'n_evaluations': n_eval[0],
        'targets': {
            'beta_xm': BETA_XM, 'alpha_xm': ALPHA_XM,
            'beta_ym': BETA_YM, 'alpha_ym': ALPHA_YM,
            'epsilon': EPSILON, 'beta_0': BETA_0,
        },
        'mse': float(best_mse[0]),
        'twiss_undulator': {
            'beta_x': float(bx), 'alpha_x': float(ax),
            'beta_y': float(by), 'alpha_y': float(ay),
        },
        'currents': {k: float(v) for k, v in sorted(currents.items(), key=lambda x: int(x[0]))},
    }
    with open(args.save, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"\nResults saved to {args.save}")


if __name__ == '__main__':
    main()
