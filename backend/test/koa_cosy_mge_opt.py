#!/usr/bin/env python3
"""Self-contained COSY FR3+MGE optimization using Python CMA-ES.

Designed for deployment on HPC clusters (Koa) without glyfada.
All COSY interaction is via FOX file generation + subprocess.

Usage:
    python3 koa_cosy_mge_opt.py [--sigma 0.2] [--max-eval 5000] [--popsize 30]

Requires: numpy, cma, COSY binary (./cosy) and COSY.bin in working directory.

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
import atexit
import numpy as np

# ── Beam parameters ──────────────────────────────────────────────────────────
ENERGY = 40  # MeV
EPSILON_N = 8  # pi.mm.mrad (normalized)
X_STD = 0.8  # mm

E0 = 0.511  # MeV, electron rest mass
gamma = (ENERGY + E0) / E0
beta_rel = np.sqrt(1 - 1 / gamma**2)
norm = gamma * beta_rel
EPSILON = EPSILON_N / norm  # geometric emittance

K = 1.2
LAMBDA_U = 2.3e-2  # m
BETA_YM = gamma * LAMBDA_U / (2 * np.pi * K)
ALPHA_YM = 0.0
BETA_XM = 1.4
ALPHA_XM = 0.47
BETA_0 = X_STD**2 / EPSILON

# ── Variable definitions ─────────────────────────────────────────────────────
# (var_name, default_value, lo_bound, hi_bound)
# Defaults from FR3 warm-converged result (MSE=4.39e-9 without MGE).
VARIABLES = [
    ('S1_I',    0.842235,  -10, 10),
    ('S1_I2',   1.057249,  -10, 10),
    ('S2_I',    3.979245,  -10, 10),
    ('S3_I',    3.441128,  -10, 10),
    ('S3_I2',   4.693946,  -10, 10),
    ('S3_I3',   1.218056,  -10, 10),
    ('S4_I',    4.719248,  -10, 10),
    ('S5_I3',  -2.755024,  -10, 10),
    ('S5_I2',  -3.727535,  -10, 10),
    ('S5_I',   -1.986620,  -10, 10),
    ('S6_I',    4.723484,  -10, 10),
    ('S7_I',    3.227473,  -10, 10),
    ('S7_I2',   3.361248,  -10, 10),
    ('S8_I',    5.212214,  -10, 10),
    ('S8_I2',   3.892062,  -10, 10),
    ('S9_I',    4.711847,  -10, 10),
    ('S10_I',   3.800455,  -10, 10),
    ('S10_I2',  5.464208,  -10, 10),
    ('S10_I3',  2.078793,  -10, 10),
    ('S11_Ic',  2.057687,  -10, 10),
    ('S11_I',   3.178241,  -10, 10),
    ('S11_I2',  3.536943,  -10, 10),
    ('S11_I3',  3.590107,  -10, 10),
]


VAR_NAMES = {v[0] for v in VARIABLES}


def generate_fox_template(fox_path):
    """Convert a FIT-enabled FOX file into a CMA-ES evaluation template.

    Strips all FIT...ENDFIT blocks, post-FIT verification sections, and
    standalone LATTICE calls. Replaces all stage variable assignments with
    __var_name__ placeholders. Inserts LATTICE call before PM 99.
    """
    with open(fox_path) as f:
        lines = f.readlines()

    output = []
    in_fit = False
    in_postfit = False

    for line in lines:
        stripped = line.strip()

        if re.match(r'\s*FIT\s+S\d', stripped):
            in_fit = True
            continue
        if in_fit:
            if re.match(r'\s*ENDFIT\s', stripped):
                in_fit = False
                in_postfit = True
            continue

        if in_postfit:
            if ('POST-FIT' in stripped or
                stripped.startswith('LATTICE') or
                stripped.startswith('LOOP IDX') or
                stripped.startswith('B0(') or
                stripped.startswith('A0(') or
                stripped.startswith('G0(') or
                stripped.startswith('WRITE 6')):
                continue
            in_postfit = False

        if 'POST-FIT' in stripped:
            continue
        if stripped == 'LATTICE ;':
            continue

        # Replace ANY stage variable assignment with placeholder
        m = re.match(r'(S\d+_\w+)\s*:=\s*[\d.e+-]+\s*;', stripped)
        if m and m.group(1) in VAR_NAMES:
            indent = line[:len(line) - len(line.lstrip())]
            var = m.group(1)
            output.append(f"{indent}{var} := __{var}__ ;\n")
            continue

        if stripped.startswith('CO 1 ; PM 99'):
            output.append("    LATTICE ;\n")

        output.append(line)

    template = ''.join(output)

    # Verify all placeholders present
    for name in VAR_NAMES:
        if f'__{name}__' not in template:
            raise ValueError(f"Placeholder __{name}__ not found in template")

    return template


def read_transfer_map(results_dir):
    """Read 6x6 linear transfer map from COSY fort.99 (PM command output).

    COSY PM outputs 5 values per line (x, x', y, y', l) — the δK row
    is omitted (identity in magnetic lattices). Values can concatenate
    without whitespace in Fortran format when they overflow column width.
    """
    import re
    fort99 = os.path.join(results_dir, 'fort.99')
    M = np.eye(6)  # identity handles M(5,5) = 1

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
    """Compute objective from 6x6 transfer map.

    Returns stability-aware objective:
    - Unstable: log-scale penalty from half-trace
    - Stable: Twiss MSE at undulator entrance
    """
    cos_mu_x = (M[0, 0] + M[1, 1]) / 2
    cos_mu_y = (M[2, 2] + M[3, 3]) / 2
    instability = max(abs(cos_mu_x), abs(cos_mu_y))

    if instability > 1:
        return 1e3 * (1 + math.log(instability))

    b0, g0 = BETA_0, 1.0 / BETA_0
    bx = M[0, 0]**2 * b0 + M[0, 1]**2 * g0
    ax = -(M[0, 0] * M[1, 0] * b0 + M[0, 1] * M[1, 1] * g0)
    by = M[2, 2]**2 * b0 + M[2, 3]**2 * g0
    ay = -(M[2, 2] * M[3, 2] * b0 + M[2, 3] * M[3, 3] * g0)

    mse = ((bx - BETA_XM)**2 + (by - BETA_YM)**2 +
           (ax - ALPHA_XM)**2 + (ay - ALPHA_YM)**2) / 4

    return mse if math.isfinite(mse) else 1e6


def setup_cosy_binary(run_dir):
    """Locate COSY binary: $COSY_BIN env var, or copy to /tmp if noexec.

    Returns the path to the executable cosy binary.
    """
    # Allow SLURM job script to pre-stage the binary
    env_bin = os.environ.get('COSY_BIN')
    if env_bin and os.path.isfile(env_bin):
        print(f"  COSY binary from $COSY_BIN: {env_bin}")
        return env_bin

    cosy_src = os.path.join(run_dir, 'cosy')
    if not os.path.isfile(cosy_src):
        sys.exit("ERROR: 'cosy' binary not found in working directory")

    # Try executing from run_dir first
    try:
        subprocess.run([cosy_src], capture_output=True, timeout=5)
        return cosy_src
    except PermissionError:
        pass
    except (subprocess.TimeoutExpired, Exception):
        return cosy_src  # executable, just needs input

    # noexec filesystem — copy to /tmp
    import shutil
    import tempfile
    tmp_cosy = os.path.join(tempfile.gettempdir(), f'cosy_{os.getpid()}')
    shutil.copy2(cosy_src, tmp_cosy)
    os.chmod(tmp_cosy, 0o755)
    atexit.register(lambda: os.path.isfile(tmp_cosy) and os.remove(tmp_cosy))
    print(f"  COSY binary copied to {tmp_cosy} (noexec workaround)")
    return tmp_cosy


COSY_BIN = None  # set by main()


def evaluate(fox_template, var_names, run_dir, currents):
    """Run one COSY evaluation with given quad currents."""
    fox = fox_template
    for name, val in zip(var_names, currents):
        fox = fox.replace(f"__{name}__", str(val))

    with open(os.path.join(run_dir, 'input.fox'), 'w') as f:
        f.write(fox)

    try:
        result = subprocess.run(
            [COSY_BIN, 'input.fox'], cwd=run_dir,
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
        sys.exit("ERROR: 'cma' package required. Install: python3 -m pip install --user cma")

    parser = argparse.ArgumentParser(description='COSY FR3+MGE CMA-ES optimization')
    parser.add_argument('--sigma', type=float, default=0.5,
                        help='CMA-ES initial step size')
    parser.add_argument('--max-eval', type=int, default=10000,
                        help='Maximum function evaluations')
    parser.add_argument('--popsize', type=int, default=0,
                        help='Population size (0=CMA default)')
    parser.add_argument('--restarts', type=int, default=0,
                        help='Number of BIPOP restarts (0=none)')
    parser.add_argument('--warm-start', type=str, default=None,
                        help='JSON result file with currents for warm start')
    parser.add_argument('--save', type=str, default='result_cosy_mge.json',
                        help='Output JSON file')
    parser.add_argument('--fox-template', type=str, default=None,
                        help='Pre-built FOX template file')
    parser.add_argument('--fox-source', type=str, default=None,
                        help='FIT-enabled FOX file to auto-generate template from')
    parser.add_argument('--run-dir', type=str, default=None,
                        help='COSY run directory (default: cwd)')
    args = parser.parse_args()

    var_names = [v[0] for v in VARIABLES]
    defaults = [v[1] for v in VARIABLES]
    lo = [v[2] for v in VARIABLES]
    hi = [v[3] for v in VARIABLES]

    # Warm start from result JSON (maps stage var names to values)
    if args.warm_start:
        with open(args.warm_start) as f:
            warm = json.load(f)
        warm_currents = warm.get('currents', warm)
        # Accept either {var_name: val} or {elem_idx: val} format
        if any(k.startswith('S') for k in warm_currents):
            for i, name in enumerate(var_names):
                if name in warm_currents:
                    defaults[i] = warm_currents[name]
        # No automatic element-index mapping — use var names only

    # FOX template: auto-generate from source or read pre-built
    if args.fox_source:
        fox_template = generate_fox_template(args.fox_source)
    elif args.fox_template:
        with open(args.fox_template) as f:
            fox_template = f.read()
    else:
        sys.exit("ERROR: Provide --fox-source or --fox-template")

    if args.run_dir:
        run_dir = os.path.abspath(args.run_dir)
        os.makedirs(run_dir, exist_ok=True)
    else:
        run_dir = os.getcwd()

    global COSY_BIN
    COSY_BIN = setup_cosy_binary(run_dir)

    print(f"COSY FR3+MGE CMA-ES Optimization")
    print(f"  Variables: {len(var_names)}")
    print(f"  Sigma: {args.sigma}, Max eval: {args.max_eval}")
    print(f"  Targets: bx={BETA_XM}, ax={ALPHA_XM}, by={BETA_YM:.4f}, ay={ALPHA_YM}")
    print(f"  Initial beta_0: {BETA_0:.4f} m")

    # Test one evaluation
    t0 = time.time()
    test_mse = evaluate(fox_template, var_names, run_dir, defaults)
    dt = time.time() - t0
    print(f"  Test evaluation: RMS={math.sqrt(test_mse):.4e} ({dt:.1f}s)")
    est_hours = dt * args.max_eval / 3600
    print(f"  Estimated time: {est_hours:.1f} hours for {args.max_eval} evaluations")

    n_eval = [0]
    best_mse = [test_mse]
    best_x = [list(defaults)]
    t_start = time.time()

    def objective_wrapper(x):
        mse = evaluate(fox_template, var_names, run_dir, list(x))
        n_eval[0] += 1
        if mse < best_mse[0]:
            best_mse[0] = mse
            best_x[0] = list(x)
            elapsed = time.time() - t_start
            stable = "STABLE" if mse < 1000 else "unstable"
            print(f"  [{n_eval[0]:5d}] RMS={math.sqrt(mse):.6e} ({stable}) t={elapsed:.0f}s")
        return mse

    opts = {
        'maxfevals': args.max_eval,
        'bounds': [lo, hi],
        'seed': 42,
        'verb_disp': 500,
        'verb_log': 0,
        'tolfun': 1e-10,
        'tolx': 1e-4,
    }
    if args.popsize > 0:
        opts['popsize'] = args.popsize

    if args.restarts > 0:
        xopt, es = cma.fmin2(objective_wrapper, defaults, args.sigma, opts,
                             restarts=args.restarts, bipop=True)
    else:
        es = cma.CMAEvolutionStrategy(defaults, args.sigma, opts)
        es.optimize(objective_wrapper)

    result = es.result
    print(f"\nOptimization complete: {n_eval[0]} evaluations")
    print(f"  Best RMS: {math.sqrt(best_mse[0]):.6e}")
    print(f"  Stable: {'YES' if best_mse[0] < 1000 else 'NO'}")

    print("\nOptimal currents:")
    for name, val in zip(var_names, best_x[0]):
        print(f"  {name:12s} = {val:.6f}")

    # Final evaluation to get Twiss
    evaluate(fox_template, var_names, run_dir, best_x[0])
    M = read_transfer_map(run_dir)

    b0, g0 = BETA_0, 1.0 / BETA_0
    bx = M[0, 0]**2 * b0 + M[0, 1]**2 * g0
    ax = -(M[0, 0] * M[1, 0] * b0 + M[0, 1] * M[1, 1] * g0)
    by = M[2, 2]**2 * b0 + M[2, 3]**2 * g0
    ay = -(M[2, 2] * M[3, 2] * b0 + M[2, 3] * M[3, 3] * g0)

    print(f"\nFinal Twiss:")
    print(f"  beta_x  = {bx:.6f}  (target {BETA_XM})")
    print(f"  beta_y  = {by:.6f}  (target {BETA_YM:.6f})")
    print(f"  alpha_x = {ax:.6f}  (target {ALPHA_XM})")
    print(f"  alpha_y = {ay:.6f}  (target {ALPHA_YM})")

    data = {
        'config': 'S1_2ps',
        'energy_MeV': ENERGY,
        'epsilon_n': EPSILON_N,
        'fringe_field_order': 3,
        'mge': True,
        'optimizer': 'cma-es',
        'n_evaluations': n_eval[0],
        'mse': float(best_mse[0]),
        'twiss_undulator': {
            'beta_x': float(bx), 'alpha_x': float(ax),
            'beta_y': float(by), 'alpha_y': float(ay),
        },
        'currents': {name: float(val) for name, val in zip(var_names, best_x[0])},
    }
    with open(args.save, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"\nResults saved to {args.save}")


if __name__ == '__main__':
    main()
