#!/usr/bin/env python3
"""C3v4: COSY FR3+MGE optimization with stability-guaranteed objective.

Key fix over v1-v3: the original objective allows unstable solutions
(penalty ~1000) to beat barely-stable solutions with bad Twiss (MSE >> 1000).
CMA-ES then *prefers* the unstable side of the boundary.

This version caps stable-region MSE at 999, guaranteeing that ANY stable
solution ranks better than ANY unstable one. Once CMA-ES crosses the
stability boundary, it can then refine Twiss within the stable basin.

Approaches tested (in order of expected success):
  A) Capped objective + cold start from [0]*23 + large sigma + BIPOP
  B) Capped objective + warm start from v1 result + BIPOP
  C) Two-phase: stability-only objective for first half, then Twiss MSE

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
# Defaults from v1 best result (MSE=1029.6, unstable but near boundary).
VARIABLES_V1_WARM = [
    ('S1_I',    1.843955,  -10, 10),
    ('S1_I2',   1.342348,  -10, 10),
    ('S2_I',    3.543533,  -10, 10),
    ('S3_I',    1.185194,  -10, 10),
    ('S3_I2',   4.089469,  -10, 10),
    ('S3_I3',   3.161799,  -10, 10),
    ('S4_I',    3.646672,  -10, 10),
    ('S5_I3',   0.007356,  -10, 10),
    ('S5_I2',  -2.059263,  -10, 10),
    ('S5_I',   -1.903788,  -10, 10),
    ('S6_I',    4.834105,  -10, 10),
    ('S7_I',    3.015729,  -10, 10),
    ('S7_I2',   3.427726,  -10, 10),
    ('S8_I',    4.839480,  -10, 10),
    ('S8_I2',   3.843532,  -10, 10),
    ('S9_I',    4.598188,  -10, 10),
    ('S10_I',   3.653179,  -10, 10),
    ('S10_I2',  4.227379,  -10, 10),
    ('S10_I3',  0.155622,  -10, 10),
    ('S11_Ic',  1.628827,  -10, 10),
    ('S11_I',   1.609510,  -10, 10),
    ('S11_I2',  3.618329,  -10, 10),
    ('S11_I3',  3.583319,  -10, 10),
]

# Cold start: center of bounds
VARIABLES_COLD = [(name, 0.0, lo, hi) for name, _, lo, hi in VARIABLES_V1_WARM]

VAR_NAMES = {v[0] for v in VARIABLES_V1_WARM}


def generate_fox_template(fox_path):
    """Convert a FIT-enabled FOX file into a CMA-ES evaluation template."""
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
    for name in VAR_NAMES:
        if f'__{name}__' not in template:
            raise ValueError(f"Placeholder __{name}__ not found in template")
    return template


def read_transfer_map(results_dir):
    """Read 6x6 linear transfer map from COSY fort.99."""
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


def compute_objective_capped(M):
    """Stability-guaranteed objective: stable solutions ALWAYS beat unstable.

    Unstable: 1000 + 1000*log(instability), range [1000, ~13000]
    Stable:   min(Twiss_MSE, 999), range [0, 999]

    The cap at 999 ensures CMA-ES always prefers crossing the stability
    boundary. Without it, barely-stable solutions with bad Twiss (MSE >> 1000)
    rank worse than the unstable penalty, trapping the optimizer.
    """
    cos_mu_x = (M[0, 0] + M[1, 1]) / 2
    cos_mu_y = (M[2, 2] + M[3, 3]) / 2
    instability = max(abs(cos_mu_x), abs(cos_mu_y))

    if instability > 1:
        return 1000 + 1000 * math.log(instability)

    b0, g0 = BETA_0, 1.0 / BETA_0
    bx = M[0, 0]**2 * b0 + M[0, 1]**2 * g0
    ax = -(M[0, 0] * M[1, 0] * b0 + M[0, 1] * M[1, 1] * g0)
    by = M[2, 2]**2 * b0 + M[2, 3]**2 * g0
    ay = -(M[2, 2] * M[3, 2] * b0 + M[2, 3] * M[3, 3] * g0)

    mse = ((bx - BETA_XM)**2 + (by - BETA_YM)**2 +
           (ax - ALPHA_XM)**2 + (ay - ALPHA_YM)**2) / 4

    if not math.isfinite(mse):
        return 999.0

    return min(mse, 999.0)


def compute_objective_stability_only(M):
    """Phase 1 objective: minimize instability measure, ignore Twiss entirely."""
    cos_mu_x = (M[0, 0] + M[1, 1]) / 2
    cos_mu_y = (M[2, 2] + M[3, 3]) / 2
    instability = max(abs(cos_mu_x), abs(cos_mu_y))

    if instability <= 1:
        # Stable — return negative reward proportional to stability margin
        return -1.0 + instability  # range [-1, 0) for stable, closer to -1 = more stable
    return instability  # range (1, inf) for unstable


def setup_cosy_binary(run_dir):
    """Locate COSY binary: $COSY_BIN env var, or copy to /tmp if noexec."""
    env_bin = os.environ.get('COSY_BIN')
    if env_bin and os.path.isfile(env_bin):
        print(f"  COSY binary from $COSY_BIN: {env_bin}")
        return env_bin

    cosy_src = os.path.join(run_dir, 'cosy')
    if not os.path.isfile(cosy_src):
        sys.exit("ERROR: 'cosy' binary not found in working directory")

    try:
        subprocess.run([cosy_src], capture_output=True, timeout=5)
        return cosy_src
    except PermissionError:
        pass
    except (subprocess.TimeoutExpired, Exception):
        return cosy_src

    import shutil
    import tempfile
    tmp_cosy = os.path.join(tempfile.gettempdir(), f'cosy_{os.getpid()}')
    shutil.copy2(cosy_src, tmp_cosy)
    os.chmod(tmp_cosy, 0o755)
    atexit.register(lambda: os.path.isfile(tmp_cosy) and os.remove(tmp_cosy))
    print(f"  COSY binary copied to {tmp_cosy} (noexec workaround)")
    return tmp_cosy


COSY_BIN = None


def evaluate(fox_template, var_names, run_dir, currents, objective_fn):
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

    return objective_fn(M)


def main():
    try:
        import cma
    except ImportError:
        sys.exit("ERROR: 'cma' package required. Install: pip install cma")

    parser = argparse.ArgumentParser(description='C3v4: COSY FR3+MGE with stability-guaranteed objective')
    parser.add_argument('--approach', choices=['A', 'B', 'C'], default='A',
                        help='A=capped+cold, B=capped+warm, C=two-phase')
    parser.add_argument('--sigma', type=float, default=2.0,
                        help='CMA-ES initial step size (default 2.0 for broad search)')
    parser.add_argument('--max-eval', type=int, default=5000,
                        help='Max evaluations per BIPOP restart')
    parser.add_argument('--restarts', type=int, default=15,
                        help='Number of BIPOP restarts')
    parser.add_argument('--warm-start', type=str, default=None,
                        help='JSON result file for warm start (approach B)')
    parser.add_argument('--save', type=str, default='result_cosy_mge_v4.json')
    parser.add_argument('--fox-template', type=str, default=None)
    parser.add_argument('--fox-source', type=str, default=None)
    parser.add_argument('--run-dir', type=str, default=None)
    args = parser.parse_args()

    # Select variable set based on approach
    if args.approach == 'A':
        variables = VARIABLES_COLD
        desc = "Cold start [0]*23, capped objective, sigma=2.0"
    elif args.approach == 'B':
        variables = VARIABLES_V1_WARM
        desc = "Warm start from v1, capped objective"
    else:
        variables = VARIABLES_COLD
        desc = "Two-phase: stability-only → Twiss MSE"

    var_names = [v[0] for v in variables]
    defaults = [v[1] for v in variables]
    lo = [v[2] for v in variables]
    hi = [v[3] for v in variables]

    if args.warm_start:
        with open(args.warm_start) as f:
            warm = json.load(f)
        warm_currents = warm.get('currents', warm)
        if any(k.startswith('S') for k in warm_currents):
            for i, name in enumerate(var_names):
                if name in warm_currents:
                    defaults[i] = warm_currents[name]

    if args.fox_source:
        fox_template = generate_fox_template(args.fox_source)
    elif args.fox_template:
        with open(args.fox_template) as f:
            fox_template = f.read()
    else:
        sys.exit("ERROR: Provide --fox-source or --fox-template")

    run_dir = os.path.abspath(args.run_dir) if args.run_dir else os.getcwd()
    os.makedirs(run_dir, exist_ok=True)

    global COSY_BIN
    COSY_BIN = setup_cosy_binary(run_dir)

    print(f"C3v4: COSY FR3+MGE Optimization")
    print(f"  Approach {args.approach}: {desc}")
    print(f"  Variables: {len(var_names)}")
    print(f"  Sigma: {args.sigma}, Max eval/restart: {args.max_eval}, Restarts: {args.restarts}")
    print(f"  Max total evals: ~{args.max_eval * (args.restarts + 1)}")
    print(f"  Targets: bx={BETA_XM}, ax={ALPHA_XM}, by={BETA_YM:.4f}, ay={ALPHA_YM}")

    # Select objective function
    if args.approach == 'C':
        # Phase 1: stability-only for first half of evals
        phase1_budget = args.max_eval * (args.restarts + 1) // 2
        objective_fn = compute_objective_stability_only
        phase_desc = "Phase 1 (stability-only)"
    else:
        objective_fn = compute_objective_capped
        phase_desc = "Capped objective"
    print(f"  Objective: {phase_desc}")

    # Test evaluation
    t0 = time.time()
    test_mse = evaluate(fox_template, var_names, run_dir, defaults, objective_fn)
    dt = time.time() - t0
    est_hours = dt * args.max_eval * (args.restarts + 1) / 3600
    print(f"  Test eval: obj={test_mse:.4e} ({dt:.1f}s)")
    print(f"  Est. max time: {est_hours:.1f}h")

    n_eval = [0]
    best_obj = [float('inf')]
    best_x = [list(defaults)]
    first_stable = [None]
    t_start = time.time()

    def objective_wrapper(x):
        nonlocal objective_fn

        # Two-phase: switch objective after budget
        if args.approach == 'C' and n_eval[0] >= phase1_budget:
            if objective_fn is compute_objective_stability_only:
                objective_fn = compute_objective_capped
                print(f"\n  === Phase 2: switching to capped Twiss MSE objective ===\n")

        obj = evaluate(fox_template, var_names, run_dir, list(x), objective_fn)
        n_eval[0] += 1

        # Track stability crossings
        if obj < 0 or (0 <= obj < 1000):
            if first_stable[0] is None:
                first_stable[0] = n_eval[0]
                elapsed = time.time() - t_start
                print(f"  *** FIRST STABLE SOLUTION at eval {n_eval[0]} "
                      f"(obj={obj:.6e}, t={elapsed:.0f}s) ***")

        if obj < best_obj[0]:
            best_obj[0] = obj
            best_x[0] = list(x)
            elapsed = time.time() - t_start
            stable = "STABLE" if obj < 1000 else "unstable"
            print(f"  [{n_eval[0]:5d}] obj={obj:.6e} ({stable}) t={elapsed:.0f}s")

        return obj

    opts = {
        'maxfevals': args.max_eval,
        'bounds': [lo, hi],
        'seed': 42,
        'verb_disp': 500,
        'verb_log': 0,
        'tolfun': 1e-12,
        'tolx': 1e-5,
    }

    xopt, es = cma.fmin2(objective_wrapper, defaults, args.sigma, opts,
                         restarts=args.restarts, bipop=True)

    print(f"\nOptimization complete: {n_eval[0]} evaluations")
    print(f"  Best objective: {best_obj[0]:.6e}")
    stable = best_obj[0] < 1000 and best_obj[0] >= 0
    print(f"  Stable: {'YES' if stable else 'NO'}")
    if first_stable[0]:
        print(f"  First stable solution at eval {first_stable[0]}")

    # Final evaluation with uncapped objective for true MSE
    from koa_cosy_mge_opt import compute_objective as compute_objective_uncapped
    final_mse = evaluate(fox_template, var_names, run_dir, best_x[0],
                         compute_objective_uncapped)
    print(f"  True (uncapped) MSE: {final_mse:.6e}")

    # Read final Twiss
    M = read_transfer_map(run_dir)
    b0, g0 = BETA_0, 1.0 / BETA_0
    bx = M[0, 0]**2 * b0 + M[0, 1]**2 * g0
    ax = -(M[0, 0] * M[1, 0] * b0 + M[0, 1] * M[1, 1] * g0)
    by = M[2, 2]**2 * b0 + M[2, 3]**2 * g0
    ay = -(M[2, 2] * M[3, 2] * b0 + M[2, 3] * M[3, 3] * g0)

    cos_mu_x = (M[0, 0] + M[1, 1]) / 2
    cos_mu_y = (M[2, 2] + M[3, 3]) / 2

    print(f"\nFinal state:")
    print(f"  cos_mu_x = {cos_mu_x:.8f}  (|.| <= 1 for stable)")
    print(f"  cos_mu_y = {cos_mu_y:.8f}")
    print(f"  beta_x   = {bx:.6f}  (target {BETA_XM})")
    print(f"  beta_y   = {by:.6f}  (target {BETA_YM:.6f})")
    print(f"  alpha_x  = {ax:.6f}  (target {ALPHA_XM})")
    print(f"  alpha_y  = {ay:.6f}  (target {ALPHA_YM})")

    print("\nOptimal currents:")
    for name, val in zip(var_names, best_x[0]):
        print(f"  {name:12s} = {val:.6f}")

    data = {
        'config': 'S1_2ps',
        'energy_MeV': ENERGY,
        'epsilon_n': EPSILON_N,
        'fringe_field_order': 3,
        'mge': True,
        'optimizer': f'cma-es-v4-{args.approach}',
        'approach': args.approach,
        'approach_desc': desc,
        'n_evaluations': n_eval[0],
        'first_stable_eval': first_stable[0],
        'capped_objective': float(best_obj[0]),
        'mse': float(final_mse),
        'stability': {
            'cos_mu_x': float(cos_mu_x),
            'cos_mu_y': float(cos_mu_y),
            'stable': stable,
        },
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
