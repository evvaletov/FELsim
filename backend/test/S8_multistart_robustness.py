"""S8: Multi-Start Robustness Study

Runs 10 random starts at 5 extreme emittance points to characterize the
optimizer landscape. Addresses S7's finding that single-start results at
extreme emittances depend on the particle realization.

Key question: Is the NM landscape convex enough for single-start, or are
there multiple local minima with different quality classifications?

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

OUTDIR = Path(__file__).resolve().parent / 'results' / 'S8'

EMITTANCE_POINTS = [1, 3, 5, 14, 16]
N_STARTS = 10
STAGE11_QUADS = [87, 93, 95, 97]


def _print(msg):
    print(msg, flush=True)


def classify(mse):
    if math.isnan(mse):
        return 'Failed'
    if mse < MSE_THRESHOLDS['Excellent']:
        return 'Excellent'
    elif mse < MSE_THRESHOLDS['Acceptable']:
        return 'Acceptable'
    elif mse < MSE_THRESHOLDS['Marginal']:
        return 'Marginal'
    return 'Failed'


def generate_starts(n_starts, seed=42):
    """Generate n_starts random Stage 11 starting points."""
    rng = np.random.RandomState(seed)
    starts = []
    for _ in range(n_starts):
        starts.append({
            "Ic": {"bounds": (0, 10), "start": rng.uniform(0, 10)},
            "I":  {"bounds": (0, 10), "start": rng.uniform(0, 10)},
            "I2": {"bounds": (0, 10), "start": rng.uniform(0, 10)},
            "I3": {"bounds": (0, 10), "start": rng.uniform(0, 10)},
        })
    return starts


def run_multistart(nb_particles=500, seed=42):
    csv_path = OUTDIR / 'multistart_emittance.csv'

    header = ['epsilon_n', 'start_idx', 'mse', 'alpha_x', 'alpha_y',
              'beta_x', 'beta_y', 'disp_resid',
              'q87', 'q93', 'q95', 'q97', 'time_s', 'nfev']

    # Load checkpoint
    completed = set()
    if csv_path.exists():
        with open(csv_path) as f:
            for r in csv.DictReader(f):
                completed.add((float(r['epsilon_n']), int(r['start_idx'])))

    if not csv_path.exists():
        with open(csv_path, 'w', newline='') as f:
            csv.writer(f).writerow(header)

    total = len(EMITTANCE_POINTS) * N_STARTS
    done = len(completed)
    starts = generate_starts(N_STARTS, seed=seed + 7777)

    _print(f"S8: Multi-Start Robustness Study")
    _print(f"  {len(EMITTANCE_POINTS)} emittance points × {N_STARTS} starts = {total}")
    _print(f"  Already completed: {done}")

    for eps_n in EMITTANCE_POINTS:
        _print(f"\n── ε_n = {eps_n} ──")
        for si, sp in enumerate(starts):
            if (eps_n, si) in completed:
                continue
            done += 1
            _print(f"  [{done}/{total}] ε_n={eps_n}, start {si}")

            kwargs = dict(BASELINE)
            kwargs['epsilon_n'] = eps_n
            kwargs['nb_particles'] = nb_particles
            kwargs['seed'] = seed
            kwargs['stage11_startPoint'] = sp

            try:
                res = run_optimization(**kwargs)
                row = [eps_n, si, res['mse'], res['alpha_x'], res['alpha_y'],
                       res['beta_x'], res['beta_y'], res['disp_resid'],
                       res['quad_currents'][87], res['quad_currents'][93],
                       res['quad_currents'][95], res['quad_currents'][97],
                       res['time_s'], res.get('nfev', 0)]
                _print(f"           MSE={res['mse']:.2e} ({classify(res['mse'])}), "
                       f"t={res['time_s']:.1f}s")
            except Exception as e:
                _print(f"           FAILED: {e}")
                row = [eps_n, si] + [float('nan')] * (len(header) - 2)

            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow(row)


def plot_results():
    csv_path = OUTDIR / 'multistart_emittance.csv'
    if not csv_path.exists():
        _print("No data to plot.")
        return

    with open(csv_path) as f:
        rows = list(csv.DictReader(f))

    eps_vals = sorted(set(float(r['epsilon_n']) for r in rows))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel 1: MSE box plot
    mse_data = []
    labels = []
    for eps in eps_vals:
        mses = [float(r['mse']) for r in rows
                if float(r['epsilon_n']) == eps and not math.isnan(float(r['mse']))]
        if mses:
            mse_data.append(mses)
            labels.append(f'{eps:.0f}')

    # Convert MSE → RMS for display
    rms_data = [[math.sqrt(m) for m in group] for group in mse_data]
    bp = ax1.boxplot(rms_data, labels=labels, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('#3498db')
        patch.set_alpha(0.6)

    ax1.set_yscale('log')
    for name, thresh in MSE_THRESHOLDS.items():
        rms_thresh = math.sqrt(thresh)
        ax1.axhline(rms_thresh, ls=':', color='gray', alpha=0.5)
        ax1.text(len(labels) + 0.3, rms_thresh, name, fontsize=7, color='gray',
                 va='center')
    ax1.set_xlabel(r'$\varepsilon_n$ ($\pi$·mm·mrad)')
    ax1.set_ylabel('RMS Twiss Mismatch')
    ax1.set_title('S8: Multi-Start RMS Distribution')
    ax1.grid(True, alpha=0.3, axis='y')

    # Panel 2: Quality classification stacked bar
    categories = ['Excellent', 'Acceptable', 'Marginal', 'Failed']
    colors = {'Excellent': '#2ecc71', 'Acceptable': '#f1c40f',
              'Marginal': '#e67e22', 'Failed': '#e74c3c'}
    x = np.arange(len(eps_vals))
    bottom = np.zeros(len(eps_vals))

    for cat in categories:
        counts = []
        for eps in eps_vals:
            n = sum(1 for r in rows
                    if float(r['epsilon_n']) == eps
                    and classify(float(r['mse'])) == cat)
            counts.append(n)
        ax2.bar(x, counts, bottom=bottom, color=colors[cat], label=cat, width=0.6)
        bottom += np.array(counts)

    ax2.set_xticks(x)
    ax2.set_xticklabels([f'{e:.0f}' for e in eps_vals])
    ax2.set_xlabel(r'$\varepsilon_n$ ($\pi$·mm·mrad)')
    ax2.set_ylabel('Count (out of 10)')
    ax2.set_title('S8: Quality Classification by Start Point')
    ax2.legend(fontsize=9)
    ax2.set_ylim(0, N_STARTS + 1)

    plt.tight_layout()
    for ext in ['eps', 'png', 'pdf']:
        fig.savefig(OUTDIR / f'S8_multistart_robustness.{ext}', dpi=150,
                    bbox_inches='tight')
    _print(f"  Saved: S8_multistart_robustness.{{eps,png,pdf}}")
    plt.close(fig)

    # Summary table (display RMS = sqrt(MSE))
    _print(f"\n── S8 Summary (RMS = √MSE) ──")
    _print(f"{'ε_n':>6} {'Best RMS':>12} {'Median':>12} {'Worst':>12} "
           f"{'Exc':>4} {'Acc':>4} {'Mar':>4} {'Fail':>4}")
    _print("-" * 70)

    for eps in eps_vals:
        mses = [float(r['mse']) for r in rows
                if float(r['epsilon_n']) == eps and not math.isnan(float(r['mse']))]
        if not mses:
            continue
        rms = [math.sqrt(m) for m in mses]
        cats = {c: 0 for c in categories}
        for m in mses:
            cats[classify(m)] += 1
        _print(f"{eps:6.0f} {min(rms):12.2e} {np.median(rms):12.2e} "
               f"{max(rms):12.2e} {cats['Excellent']:4d} {cats['Acceptable']:4d} "
               f"{cats['Marginal']:4d} {cats['Failed']:4d}")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    _print("S8: Multi-Start Robustness Study")
    _print(f"  Emittance points: {EMITTANCE_POINTS}")
    _print(f"  Starts per point: {N_STARTS}")

    run_multistart()
    plot_results()

    _print(f"\n{'=' * 60}")
    _print("  S8 Complete")
    _print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
