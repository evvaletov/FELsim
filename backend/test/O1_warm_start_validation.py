"""O1: Warm-Start from Neighbors — Validation

Re-runs the S4 emittance sweep (ε_n = 1–20) with warm_start=True and
compares against the cold-start baseline. Measures improvement in MSE,
convergence speed, and current continuity across the sweep.

Author: Eremey Valetov
"""

import sys
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
    run_scan, BASELINE, QUAD_INDICES, MSE_THRESHOLDS, read_csv,
)

OUTDIR = Path(__file__).resolve().parent / 'results' / 'O1'
S4_DIR = Path(__file__).resolve().parent / 'results' / 'params_05ps'

EMITTANCE_RANGE = np.linspace(1, 20, 20)


def _print(msg):
    print(msg, flush=True)


def classify(mse):
    if math.isnan(mse):
        return 'Failed'
    for label in ['Excellent', 'Acceptable', 'Marginal']:
        if mse < MSE_THRESHOLDS[label]:
            return label
    return 'Failed'


def run_warm_sweep():
    _print("O1: Warm-start emittance sweep (ε_n = 1–20, 20 points)")
    run_scan(
        scan_name='o1_emittance_warm',
        param_name='epsilon_n',
        param_values=EMITTANCE_RANGE,
        outdir=OUTDIR,
        nb_particles=500,
        seed=42,
        warm_start=True,
    )


def run_cold_sweep():
    _print("O1: Cold-start emittance sweep (baseline comparison)")
    run_scan(
        scan_name='o1_emittance_cold',
        param_name='epsilon_n',
        param_values=EMITTANCE_RANGE,
        outdir=OUTDIR,
        nb_particles=500,
        seed=42,
        warm_start=False,
    )


def plot_comparison():
    warm_path = OUTDIR / 'scan_o1_emittance_warm.csv'
    cold_path = OUTDIR / 'scan_o1_emittance_cold.csv'
    s4_path = S4_DIR / 'scan_emittance.csv'

    datasets = []
    for path, label in [(cold_path, 'Cold-start'), (warm_path, 'Warm-start'),
                        (s4_path, 'S4 original')]:
        if path.exists():
            datasets.append((read_csv(path), label))

    if not datasets:
        _print("No data to plot.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel 1: RMS comparison
    ax = axes[0, 0]
    markers = ['o-', 's--', '^:']
    for (rows, label), marker in zip(datasets, markers):
        x = [r['param_value'] for r in rows]
        y = [math.sqrt(r['mse']) for r in rows]
        ax.semilogy(x, y, marker, label=label, markersize=5)
    for name, thresh in MSE_THRESHOLDS.items():
        rms_thresh = math.sqrt(thresh)
        ax.axhline(rms_thresh, ls=':', color='gray', alpha=0.5)
        ax.text(0.5, rms_thresh * 1.3, name, fontsize=7, color='gray')
    ax.set_xlabel(r'$\varepsilon_n$ ($\pi$·mm·mrad)')
    ax.set_ylabel('RMS Twiss Mismatch')
    ax.set_title('RMS: Cold vs Warm Start')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel 2: RMS ratio (warm / cold)
    ax = axes[0, 1]
    if len(datasets) >= 2:
        cold_rows, warm_rows = datasets[0][0], datasets[1][0]
        if len(cold_rows) == len(warm_rows):
            x = [r['param_value'] for r in cold_rows]
            ratio = []
            for cr, wr in zip(cold_rows, warm_rows):
                if cr['mse'] > 0 and wr['mse'] > 0:
                    ratio.append(math.sqrt(wr['mse']) / math.sqrt(cr['mse']))
                else:
                    ratio.append(float('nan'))
            ax.semilogy(x, ratio, 'ko-', markersize=5)
            ax.axhline(1.0, ls='--', color='red', alpha=0.5, label='Equal')
            ax.set_xlabel(r'$\varepsilon_n$ ($\pi$·mm·mrad)')
            ax.set_ylabel('RMS ratio (warm / cold)')
            ax.set_title('Warm-Start Improvement Factor')
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)

    # Panel 3: β_x comparison
    ax = axes[1, 0]
    for (rows, label), marker in zip(datasets, markers):
        x = [r['param_value'] for r in rows]
        ax.plot(x, [r['beta_x'] for r in rows], marker, label=label + r' $\beta_x$',
                markersize=4)
    ax.axhline(1.4, ls=':', color='blue', alpha=0.4, label=r'$\beta_x$ target')
    ax.set_xlabel(r'$\varepsilon_n$ ($\pi$·mm·mrad)')
    ax.set_ylabel(r'$\beta_x$ (m)')
    ax.set_title(r'$\beta_x$ at Undulator')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 4: Stage 11 quad current continuity
    ax = axes[1, 1]
    for qi, color in zip([87, 93, 95, 97], ['#3498db', '#e74c3c', '#2ecc71', '#f39c12']):
        for (rows, label), ls in zip(datasets[:2], ['-', '--']):
            x = [r['param_value'] for r in rows]
            y = [r[f'quad_{qi}'] for r in rows]
            ax.plot(x, y, ls, color=color, markersize=3,
                    label=f'q{qi} ({label})' if ls == '-' else f'q{qi} ({label})')
    ax.set_xlabel(r'$\varepsilon_n$ ($\pi$·mm·mrad)')
    ax.set_ylabel('Current (A)')
    ax.set_title('Stage 11 Quad Current Continuity')
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    for ext in ['eps', 'png', 'pdf']:
        fig.savefig(OUTDIR / f'O1_warm_start_comparison.{ext}', dpi=150,
                    bbox_inches='tight')
    _print(f"  Saved: O1_warm_start_comparison.{{eps,png,pdf}}")
    plt.close(fig)

    # Summary table
    _print(f"\n── O1 Summary ──")
    _print(f"{'ε_n':>6} {'Cold RMS':>12} {'Warm RMS':>12} {'Ratio':>8} "
           f"{'Cold':>10} {'Warm':>10}")
    _print("-" * 62)

    if len(datasets) >= 2:
        cold_rows, warm_rows = datasets[0][0], datasets[1][0]
        for cr, wr in zip(cold_rows, warm_rows):
            eps = cr['param_value']
            ratio = wr['mse'] / cr['mse'] if cr['mse'] > 0 else float('nan')
            _print(f"{eps:6.1f} {math.sqrt(cr['mse']):12.2e} {math.sqrt(wr['mse']):12.2e} {ratio:8.2f} "
                   f"{classify(cr['mse']):>10} {classify(wr['mse']):>10}")

    # Win/loss count
    if len(datasets) >= 2:
        wins = sum(1 for cr, wr in zip(cold_rows, warm_rows)
                   if wr['mse'] < cr['mse'])
        _print(f"\n  Warm wins: {wins}/{len(cold_rows)}")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    _print("O1: Warm-Start from Neighbors — Validation")
    _print(f"  Emittance range: {EMITTANCE_RANGE[0]:.1f}–{EMITTANCE_RANGE[-1]:.1f}, "
           f"{len(EMITTANCE_RANGE)} points")

    run_cold_sweep()
    run_warm_sweep()
    plot_comparison()

    _print(f"\n{'=' * 60}")
    _print("  O1 Complete")
    _print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
