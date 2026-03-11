"""S7: Verification Runs at Key Transition Points

Re-runs S4 transition points with 1000 and 2000 particles to confirm
that 500-particle results are statistically robust.

Key points selected from S4:
  - Emittance: ε_n = 1, 3, 5, 8, 14, 16, 20 (spans full transition range)
  - Energy spread: σ_E = 0.4, 0.55, 0.7% (around the Acceptable dip)

Author: Eremey Valetov
"""

import sys
import time
import csv
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from UHM_beamline_opt_05ps_params import (
    run_optimization, BASELINE, QUAD_INDICES, MSE_THRESHOLDS,
)

OUTDIR = Path(__file__).resolve().parent / 'results' / 'S7'

# Key points from S4 transitions
EMITTANCE_POINTS = [1, 3, 5, 8, 14, 16, 20]
ENERGY_SPREAD_POINTS = [0.4, 0.55, 0.7]
PARTICLE_COUNTS = [500, 1000, 2000]


def _print(msg):
    print(msg, flush=True)


def classify(mse):
    if mse < MSE_THRESHOLDS['Excellent']:
        return 'Excellent'
    elif mse < MSE_THRESHOLDS['Acceptable']:
        return 'Acceptable'
    elif mse < MSE_THRESHOLDS['Marginal']:
        return 'Marginal'
    return 'Failed'


def run_verification(param_name, param_values, particle_counts, seed=42):
    """Run verification grid: param_values × particle_counts."""
    csv_path = OUTDIR / f'verify_{param_name}.csv'

    # Load checkpoint
    completed = set()
    if csv_path.exists():
        with open(csv_path) as f:
            for r in csv.DictReader(f):
                completed.add((float(r['param_value']), int(float(r['nb_particles']))))

    total = len(param_values) * len(particle_counts)
    _print(f"\n  {param_name}: {len(param_values)} values × {len(particle_counts)} particle counts = {total}")
    _print(f"  Already completed: {len(completed)}")

    header = ['param_value', 'nb_particles', 'mse', 'alpha_x', 'alpha_y',
              'beta_x', 'beta_y', 'disp_resid', 'time_s', 'nfev']
    if not csv_path.exists():
        with open(csv_path, 'w', newline='') as f:
            csv.writer(f).writerow(header)

    done = len(completed)
    t_start = time.time()

    for pval in param_values:
        for npart in particle_counts:
            if (pval, npart) in completed:
                continue
            done += 1

            _print(f"  [{done}/{total}] {param_name}={pval}, N={npart}")

            kwargs = dict(BASELINE)
            kwargs[param_name] = pval
            kwargs['nb_particles'] = npart
            kwargs['seed'] = seed

            try:
                t0 = time.time()
                res = run_optimization(**kwargs)
                dt = time.time() - t0

                row = [pval, npart, res['mse'], res['alpha_x'], res['alpha_y'],
                       res['beta_x'], res['beta_y'], res['disp_resid'],
                       res['time_s'], res.get('nfev', 0)]

                with open(csv_path, 'a', newline='') as f:
                    csv.writer(f).writerow(row)

                _print(f"           RMS={math.sqrt(res['mse']):.2e} ({classify(res['mse'])}), "
                       f"t={res['time_s']:.1f}s")

            except Exception as e:
                _print(f"           FAILED: {e}")
                row = [pval, npart] + [float('nan')] * (len(header) - 2)
                with open(csv_path, 'a', newline='') as f:
                    csv.writer(f).writerow(row)


def plot_verification():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, param_name, xlabel in [
        (axes[0], 'epsilon_n', r'$\varepsilon_n$ ($\pi$·mm·mrad)'),
        (axes[1], 'energy_std_percent', r'$\sigma_E$ (%)'),
    ]:
        csv_path = OUTDIR / f'verify_{param_name}.csv'
        if not csv_path.exists():
            ax.set_title(f'{param_name}: no data')
            continue

        with open(csv_path) as f:
            rows = list(csv.DictReader(f))

        for npart in PARTICLE_COUNTS:
            subset = [r for r in rows if int(float(r['nb_particles'])) == npart]
            if not subset:
                continue
            x = [float(r['param_value']) for r in subset]
            y = [math.sqrt(float(r['mse'])) for r in subset]
            ax.semilogy(x, y, 'o-', label=f'N={npart}', markersize=5)

        for name, thresh in MSE_THRESHOLDS.items():
            rms_thresh = math.sqrt(thresh)
            ax.axhline(rms_thresh, ls=':', color='gray', alpha=0.5)
            ax.text(ax.get_xlim()[0], rms_thresh * 1.5, name, fontsize=7, color='gray')

        ax.set_xlabel(xlabel)
        ax.set_ylabel('RMS Twiss Mismatch')
        ax.set_title(f'S7: Verification — {param_name}')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'S7_verification.{ext}', dpi=150)
    _print(f"  Saved: S7_verification.{{eps,png}}")
    plt.close(fig)

    # Print summary table
    _print("\n── S7 Verification Summary ──")
    _print(f"{'Param':>12} {'Value':>8} {'N=500':>12} {'N=1000':>12} {'N=2000':>12} {'Consistent?':>12}")
    _print("-" * 70)

    for param_name in ['epsilon_n', 'energy_std_percent']:
        csv_path = OUTDIR / f'verify_{param_name}.csv'
        if not csv_path.exists():
            continue
        with open(csv_path) as f:
            rows = list(csv.DictReader(f))

        values = sorted(set(float(r['param_value']) for r in rows))
        for val in values:
            mses = {}
            for r in rows:
                if float(r['param_value']) == val:
                    mses[int(float(r['nb_particles']))] = float(r['mse'])

            cats = {n: classify(m) for n, m in mses.items()}
            consistent = len(set(cats.values())) == 1
            cols = [f"{math.sqrt(mses.get(n, float('nan'))):.2e}" for n in PARTICLE_COUNTS]
            _print(f"{param_name:>12} {val:8.3g} {cols[0]:>12} {cols[1]:>12} {cols[2]:>12} "
                   f"{'Yes' if consistent else 'NO':>12}")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    _print("S7: Verification Runs at Key Transition Points")

    _print("\n── Emittance verification ──")
    run_verification('epsilon_n', EMITTANCE_POINTS, PARTICLE_COUNTS)

    _print("\n── Energy spread verification ──")
    run_verification('energy_std_percent', ENERGY_SPREAD_POINTS, PARTICLE_COUNTS)

    plot_verification()

    _print("\n" + "=" * 60)
    _print("  S7 Complete")
    _print("=" * 60)


if __name__ == "__main__":
    main()
