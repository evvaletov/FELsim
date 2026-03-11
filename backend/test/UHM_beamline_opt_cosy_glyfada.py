"""COSY FR3+MGE beamline optimization via glyfada evolutionary optimizer.

Treats COSY as a black-box evaluator: each evaluation generates a FOX file
with fixed quad currents (no internal FIT), runs COSY, reads the transfer
map, and computes Twiss MSE. Glyfada optimizes the quad currents externally.

Author: Eremey Valetov
"""

import sys
import os
import json
import math
import subprocess
import shutil
import argparse
from pathlib import Path
import numpy as np

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

from cosyAdapter import COSYAdapter
from cosyOptHelper import parse_beamline_felsim_indexed
from cosyResultsReader import COSYResultsReader
from glyfadaAdapter import GlyfadaOptimizer
from UHM_beamline_opt_cosy import (
    Energy, epsilon_n, compute_targets, build_stages,
    apply_warm_start,
)


def extract_variables(stages):
    """Extract ordered (felsim_idx, prefixed_name) for all independent quad variables."""
    variables = []
    for stage_num, stage in enumerate(stages, 1):
        prefix = f"S{stage_num}_"
        for idx, name in sorted(stage['variables'].items()):
            variables.append((idx, prefix + name))
    return variables


def setup_run_directory(sim, run_dir):
    """Copy COSY binary and support files to the run directory."""
    os.makedirs(run_dir, exist_ok=True)
    for fname in ['cosy', 'COSY.bin', 'cosy.fox', 'SYSCA.DAT']:
        src = sim._find_file(fname)
        if src:
            dst = os.path.join(run_dir, fname)
            shutil.copy(src, dst)
            if fname == 'cosy':
                os.chmod(dst, 0o755)


def build_fox_template(file_path, targets, stages, run_dir):
    """Generate FOX template with symbolic quad variables and placeholder assignments.

    Returns (fox_template, var_names) where fox_template has ``__VAR_NAME__``
    placeholders that must be substituted with numeric current values before
    running COSY.
    """
    config = {'simulation': {'KE': Energy, 'order': 3, 'dimensions': 3}}
    adapter = COSYAdapter(
        lattice_path=str(file_path), mode='transfer_matrix',
        config=config, fringe_field_order=3,
        use_mge_for_dipoles=True, debug=False,
    )
    sim = adapter.get_native_simulator()
    sim.beamline = parse_beamline_felsim_indexed(str(file_path))[:118]
    sim.set_geometric_emittance(targets['epsilon'])
    beta_0 = targets['beta_0']
    sim.set_initial_twiss(beta_x=beta_0, alpha_x=0.0, beta_y=beta_0, alpha_y=0.0)

    # Apply symbolic variable mapping for all quads (including mirrors)
    for stage_num, stage in enumerate(stages, 1):
        prefix = f"S{stage_num}_"
        var_mapping = {}
        for idx, name in stage['variables'].items():
            var_mapping[idx] = {"current": prefix + name}
        for target_idx, source_idx in stage.get('mirror', {}).items():
            if source_idx in stage['variables']:
                var_mapping[target_idx] = {"current": prefix + stage['variables'][source_idx]}
        sim.apply_variable_mapping(var_mapping, validation=False)

    setup_run_directory(sim, run_dir)
    sim.generate_input(output_dir=run_dir)

    fox_path = os.path.join(run_dir, 'input.fox')
    with open(fox_path) as f:
        fox = f.read()

    # Insert placeholder variable assignments before the OV line
    variables = extract_variables(stages)
    var_names = [name for _, name in variables]
    assignment_block = '\n'.join(
        f"    {name} := __{name}__ ;" for name in var_names
    ) + '\n'

    ov_marker = f"    OV {sim.order}"
    if ov_marker not in fox:
        raise RuntimeError(f"'{ov_marker}' not found in generated FOX")
    fox_template = fox.replace(ov_marker, assignment_block + ov_marker, 1)

    return fox_template, var_names


def make_objective(fox_template, var_names, run_dir, targets, beta_0):
    """Create a picklable objective function for glyfada.

    The returned function takes a list of quad currents and returns MSE.
    All state is captured in the closure for cloudpickle serialization.
    """
    _fox = str(fox_template)
    _names = list(var_names)
    _dir = str(run_dir)
    _bxm = float(targets['beta_xm'])
    _axm = float(targets['alpha_xm'])
    _bym = float(targets['beta_ym'])
    _aym = float(targets['alpha_ym'])
    _b0 = float(beta_0)

    def objective(currents):
        import os
        import subprocess as sp
        import math

        fox = _fox
        for name, val in zip(_names, currents):
            fox = fox.replace(f"__{name}__", str(val))

        with open(os.path.join(_dir, 'input.fox'), 'w') as f:
            f.write(fox)

        try:
            result = sp.run(
                ['./cosy', 'input.fox'], cwd=_dir,
                capture_output=True, text=True, timeout=120,
            )
        except sp.TimeoutExpired:
            return 1e6

        if result.returncode != 0:
            return 1e6
        # Only reject fatal COSY errors, not GT "unrepresentable map"
        output = result.stdout + result.stderr
        fatal = ['COMMAND PLACEMENT', 'NOT DECLARED', 'ARRAY INDEX',
                 'VARIABLE EXHAUSTED']
        if any(m in output for m in fatal):
            return 1e6

        try:
            from cosyResultsReader import COSYResultsReader as Reader
            reader = Reader(_dir)
            M = reader.read_linear_transfer_map()
        except (FileNotFoundError, ValueError, KeyError, IndexError):
            return 1e6

        # Stability check: |Tr(M)/2| must be <= 1 for stable optics
        cos_mu_x = (M[0, 0] + M[1, 1]) / 2
        cos_mu_y = (M[2, 2] + M[3, 3]) / 2
        instability = max(abs(cos_mu_x), abs(cos_mu_y))

        if instability > 1:
            # Unstable: smooth log-scale penalty gives gradient info to CMA-ES
            # |Tr/2|=1.01 → 1001, |Tr/2|=10 → 3300, |Tr/2|=180000 → 13100
            return 1e3 * (1 + math.log(instability))

        # Stable: compute Twiss MSE from transfer map (alpha0 = 0)
        b0, g0 = _b0, 1.0 / _b0
        bx = M[0, 0]**2 * b0 + M[0, 1]**2 * g0
        ax = -(M[0, 0] * M[1, 0] * b0 + M[0, 1] * M[1, 1] * g0)
        by = M[2, 2]**2 * b0 + M[2, 3]**2 * g0
        ay = -(M[2, 2] * M[3, 2] * b0 + M[2, 3] * M[3, 3] * g0)

        mse = ((bx - _bxm)**2 + (by - _bym)**2 +
               (ax - _axm)**2 + (ay - _aym)**2) / 4

        return mse if math.isfinite(mse) else 1e6

    return objective


def main():
    parser = argparse.ArgumentParser(
        description='COSY FR3+MGE beamline optimization via glyfada')
    parser.add_argument('--warm-start', type=str, default=None,
                        help='JSON with previous results for initial values')
    parser.add_argument('--pop-size', type=int, default=30)
    parser.add_argument('--max-gen', type=int, default=300)
    parser.add_argument('--sigma', type=float, default=0.1)
    parser.add_argument('--algorithm', type=str, default='CMA_ES',
                        choices=['ULS', 'CMA_ES', 'auto'])
    parser.add_argument('--save', type=str, default=None)
    parser.add_argument('--dry-run', action='store_true',
                        help='Test one evaluation without running glyfada')
    args = parser.parse_args()

    file_path = (Path(__file__).resolve().parent.parent.parent
                 / 'beam_excel' / 'Beamline_elements.xlsx')

    targets = compute_targets()
    stages = build_stages(targets)

    if args.warm_start:
        with open(args.warm_start) as f:
            warm_data = json.load(f)
        apply_warm_start(stages, warm_data['currents'])
        mse_val = warm_data.get('mse')
        mse_str = f"{mse_val:.2e}" if mse_val is not None else "?"
        rms_str = f"{math.sqrt(mse_val):.2e}" if mse_val is not None else "?"
        print(f"Warm-started from {args.warm_start} "
              f"(FR {warm_data.get('fringe_field_order', '?')}, "
              f"RMS {rms_str})")

    variables = extract_variables(stages)
    var_names = [name for _, name in variables]

    defaults, bounds = [], []
    for stage_num, stage in enumerate(stages, 1):
        for idx, name in sorted(stage['variables'].items()):
            defaults.append(stage['start_point'][name]['start'])
            bounds.append(stage['start_point'][name]['bounds'])

    run_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), 'results', 'glyfada_cosy_mge'))

    print("Building FOX template (FR3+MGE)...")
    fox_template, var_names = build_fox_template(
        file_path, targets, stages, run_dir)
    print(f"  {len(fox_template)} chars, {len(var_names)} variables")
    print(f"  Run directory: {run_dir}")

    beta_0 = targets['beta_0']
    objective = make_objective(fox_template, var_names, run_dir, targets, beta_0)

    if args.dry_run:
        print(f"\nDry run: single evaluation with default currents...")
        for name, val in zip(var_names, defaults):
            print(f"  {name:12s} = {val:.4f}")
        mse = objective(defaults)
        print(f"\n  RMS = {math.sqrt(mse):.6e}")

        if mse < 1e5:
            reader = COSYResultsReader(run_dir)
            twiss = reader.get_twiss_from_transfer_map(
                initial_twiss_x={'beta': beta_0, 'alpha': 0.0},
                initial_twiss_y={'beta': beta_0, 'alpha': 0.0},
            )
            print(f"  beta_x  = {twiss['beta_x']:.4f}  (target {targets['beta_xm']:.4f})")
            print(f"  beta_y  = {twiss['beta_y']:.4f}  (target {targets['beta_ym']:.4f})")
            print(f"  alpha_x = {twiss['alpha_x']:.4f}  (target {targets['alpha_xm']:.4f})")
            print(f"  alpha_y = {twiss['alpha_y']:.4f}  (target {targets['alpha_ym']:.4f})")
        return

    print(f"\nGlyfada optimization: {args.algorithm}, pop={args.pop_size}, "
          f"gen={args.max_gen}, sigma={args.sigma}")
    print(f"  {len(var_names)} variables, bounds [0, 10] A")

    optimizer = GlyfadaOptimizer(
        objective_func=objective,
        variable_names=var_names,
        bounds=bounds,
        default_values=defaults,
        pop_size=args.pop_size,
        max_gen=args.max_gen,
        sigma=args.sigma,
        n_processes=1,
        timeout_minutes=5,
        algorithm=args.algorithm,
        extra_config={'cma_es': {'initial_sigma': args.sigma}},
    )

    result = optimizer.optimize()

    print(f"\nResult: RMS = {math.sqrt(result.fun):.6e}")
    print("Optimal currents:")
    for name, val in zip(var_names, result.x):
        print(f"  {name:12s} = {val:.6f}")

    if args.save:
        currents = {}
        for (idx, _), val in zip(variables, result.x):
            currents[idx] = float(val)
        for stage in stages:
            for target_idx, source_idx in stage.get('mirror', {}).items():
                if source_idx in currents:
                    currents[target_idx] = currents[source_idx]

        objective(list(result.x))
        reader = COSYResultsReader(run_dir)
        twiss = reader.get_twiss_from_transfer_map(
            initial_twiss_x={'beta': beta_0, 'alpha': 0.0},
            initial_twiss_y={'beta': beta_0, 'alpha': 0.0},
        )

        data = {
            'config': 'S1_2ps',
            'energy_MeV': Energy,
            'epsilon_n': epsilon_n,
            'fringe_field_order': 3,
            'mge': True,
            'optimizer': 'glyfada',
            'algorithm': args.algorithm,
            'targets': {k: v for k, v in targets.items()
                        if isinstance(v, (int, float))},
            'mse': result.fun,
            'twiss_undulator': {k: float(v) for k, v in twiss.items()},
            'currents': {str(k): float(v) for k, v in sorted(currents.items())},
        }
        with open(args.save, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"Results saved to {args.save}")


if __name__ == '__main__':
    main()
