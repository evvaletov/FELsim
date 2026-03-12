"""P12: Multi-Seed Robustness Study

S7 tested particle count sensitivity (fixed seed=42), S8 tested Stage 11
optimizer landscape (fixed seed, varied starting point). Neither varied the
random beam realization itself. P12 fills this gap: does the 11-stage
optimizer produce consistent results across different random seeds?

This answers whether single-seed results (used in all S4/S5/W2 studies) are
trustworthy or just lucky draws.

Author: Eremey Valetov
"""

import sys
import csv
import math
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use(str(Path(__file__).resolve().parent / 'felsim.mplstyle'))

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from UHM_beamline_opt_05ps_params import (
    run_optimization, BASELINE, QUAD_INDICES, MSE_THRESHOLDS,
)

OUTDIR = Path(__file__).resolve().parent / 'results' / 'P12'
CSV_PATH = OUTDIR / 'multi_seed_robustness.csv'

DEFAULT_EMITTANCES = [5, 8, 14]
DEFAULT_PARTICLES = 500
DEFAULT_N_SEEDS = 20

QUALITY_COLORS = {
    'Excellent': '#2ecc71', 'Acceptable': '#f1c40f',
    'Marginal': '#e67e22', 'Failed': '#e74c3c',
}
QUALITY_ORDER = ['Excellent', 'Acceptable', 'Marginal', 'Failed']


def _print(msg):
    print(msg, flush=True)


def classify(mse):
    if math.isnan(mse) or math.isinf(mse):
        return 'Failed'
    if mse < MSE_THRESHOLDS['Excellent']:
        return 'Excellent'
    elif mse < MSE_THRESHOLDS['Acceptable']:
        return 'Acceptable'
    elif mse < MSE_THRESHOLDS['Marginal']:
        return 'Marginal'
    return 'Failed'


def csv_header():
    cols = ['seed', 'epsilon_n', 'mse', 'alpha_x', 'alpha_y',
            'beta_x', 'beta_y', 'disp_resid']
    cols += [f'quad_{idx}' for idx in QUAD_INDICES]
    cols += ['time_s', 'nfev', 'converged']
    return cols


def load_csv(path):
    rows = []
    if not path.exists():
        return rows
    with open(path) as f:
        for r in csv.DictReader(f):
            row = {}
            for k, v in r.items():
                try:
                    row[k] = float(v)
                except (ValueError, TypeError):
                    row[k] = v
            rows.append(row)
    return rows


def run_seed_study(seeds, emittances, n_particles):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    header = csv_header()

    # Load completed pairs
    completed = set()
    if CSV_PATH.exists():
        with open(CSV_PATH) as f:
            for r in csv.DictReader(f):
                completed.add((int(float(r['seed'])), float(r['epsilon_n'])))
    else:
        with open(CSV_PATH, 'w', newline='') as f:
            csv.writer(f).writerow(header)

    total = len(seeds) * len(emittances)
    done = len(completed)

    _print(f"P12: Multi-Seed Robustness Study")
    _print(f"  {len(emittances)} emittances × {len(seeds)} seeds = {total} runs")
    _print(f"  {n_particles} particles, cold-start (no warm-starting)")
    _print(f"  Already completed: {done}")

    for eps_n in emittances:
        _print(f"\n── ε_n = {eps_n} ──")
        for seed in seeds:
            if (seed, float(eps_n)) in completed:
                continue
            done += 1
            _print(f"  [{done}/{total}] ε_n={eps_n}, seed={seed}")

            kwargs = dict(BASELINE)
            kwargs['epsilon_n'] = eps_n
            kwargs['nb_particles'] = n_particles
            kwargs['seed'] = seed

            try:
                res = run_optimization(**kwargs)
                rms = math.sqrt(res['mse'])
                row = [seed, eps_n, res['mse'], res['alpha_x'], res['alpha_y'],
                       res['beta_x'], res['beta_y'], res['disp_resid']]
                row += [res['quad_currents'][idx] for idx in QUAD_INDICES]
                row += [res['time_s'], res.get('nfev') or 0,
                        int(res.get('converged') or False)]
                _print(f"           RMS={rms:.2e} ({classify(res['mse'])}), "
                       f"t={res['time_s']:.1f}s")
            except Exception as e:
                _print(f"           FAILED: {e}")
                row = [seed, eps_n] + [float('nan')] * (len(header) - 2)

            with open(CSV_PATH, 'a', newline='') as f:
                csv.writer(f).writerow(row)


def plot_mse_boxplot(rows, outdir):
    """Box plot of RMS per emittance (log y-axis, threshold lines)."""
    eps_vals = sorted(set(r['epsilon_n'] for r in rows))

    fig, ax = plt.subplots(figsize=(8, 5))

    rms_data, labels = [], []
    for eps in eps_vals:
        mses = [r['mse'] for r in rows
                if r['epsilon_n'] == eps and not math.isnan(r['mse'])]
        if mses:
            rms_data.append([math.sqrt(m) for m in mses])
            labels.append(f'{eps:.0f}')

    bp = ax.boxplot(rms_data, tick_labels=labels, patch_artist=True,
                    widths=0.5, showfliers=True)
    for patch in bp['boxes']:
        patch.set_facecolor('#3498db')
        patch.set_alpha(0.6)

    # Overlay individual points (seeded jitter for reproducible plots)
    jitter_rng = np.random.RandomState(0)
    for i, data in enumerate(rms_data):
        x = jitter_rng.normal(i + 1, 0.04, len(data))
        ax.scatter(x, data, alpha=0.4, s=20, c='#2c3e50', zorder=3)

    ax.set_yscale('log')
    for name, thresh in MSE_THRESHOLDS.items():
        rms_thresh = math.sqrt(thresh)
        ax.axhline(rms_thresh, ls=':', color='gray', alpha=0.5)
        ax.text(len(labels) + 0.4, rms_thresh, name, fontsize=8, color='gray',
                va='center')

    ax.set_xlabel(r'$\varepsilon_n$ ($\pi$·mm·mrad)')
    ax.set_ylabel('RMS Twiss Mismatch')
    ax.set_title('P12: RMS Distribution Across Random Seeds')

    fig.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(outdir / f'rms_boxplot.{ext}')
    plt.close(fig)
    _print(f"  Saved rms_boxplot.{{pdf,png}}")


def plot_quad_spread(rows, outdir):
    """CV% of each quad current per emittance (grouped bar chart)."""
    eps_vals = sorted(set(r['epsilon_n'] for r in rows))
    quad_keys = [f'quad_{idx}' for idx in QUAD_INDICES]

    fig, ax = plt.subplots(figsize=(14, 5))

    n_eps = len(eps_vals)
    n_quads = len(quad_keys)
    x = np.arange(n_quads)
    width = 0.8 / n_eps
    colors = ['#3498db', '#e74c3c', '#2ecc71', '#9b59b6']

    for i, eps in enumerate(eps_vals):
        cvs = []
        for qk in quad_keys:
            vals = [r[qk] for r in rows
                    if r['epsilon_n'] == eps and not math.isnan(r.get(qk, float('nan')))]
            if vals and np.mean(vals) != 0:
                cvs.append(np.std(vals) / abs(np.mean(vals)) * 100)
            else:
                cvs.append(0)
        offset = (i - n_eps / 2 + 0.5) * width
        ax.bar(x + offset, cvs, width, label=f'ε_n={eps:.0f}',
               color=colors[i % len(colors)], alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels([f'q{idx}' for idx in QUAD_INDICES],
                       rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('CV (%)')
    ax.set_title('P12: Quad Current Variability Across Seeds')
    ax.legend(fontsize=9)

    fig.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(outdir / f'quad_spread.{ext}')
    plt.close(fig)
    _print(f"  Saved quad_spread.{{pdf,png}}")


def plot_quality_histogram(rows, outdir):
    """Stacked bar: quality categories per emittance."""
    eps_vals = sorted(set(r['epsilon_n'] for r in rows))

    fig, ax = plt.subplots(figsize=(7, 5))

    x = np.arange(len(eps_vals))
    bottom = np.zeros(len(eps_vals))

    for cat in QUALITY_ORDER:
        counts = []
        for eps in eps_vals:
            n = sum(1 for r in rows
                    if r['epsilon_n'] == eps and classify(r['mse']) == cat)
            counts.append(n)
        ax.bar(x, counts, bottom=bottom, color=QUALITY_COLORS[cat],
               label=cat, width=0.6)
        bottom += np.array(counts)

    ax.set_xticks(x)
    ax.set_xticklabels([f'{e:.0f}' for e in eps_vals])
    ax.set_xlabel(r'$\varepsilon_n$ ($\pi$·mm·mrad)')
    ax.set_ylabel('Count')
    ax.set_title('P12: Quality Classification by Seed')
    ax.legend(fontsize=9)
    n_seeds = max(sum(1 for r in rows if r['epsilon_n'] == eps) for eps in eps_vals)
    ax.set_ylim(0, n_seeds + 1)

    fig.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(outdir / f'quality_histogram.{ext}')
    plt.close(fig)
    _print(f"  Saved quality_histogram.{{pdf,png}}")


def plot_twiss_scatter(rows, outdir):
    """2×2 Twiss params scatter, color-coded by quality."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    params = [
        ('beta_x', r'$\beta_x$ (m)'),
        ('alpha_x', r'$\alpha_x$'),
        ('beta_y', r'$\beta_y$ (m)'),
        ('alpha_y', r'$\alpha_y$'),
    ]
    eps_vals = sorted(set(r['epsilon_n'] for r in rows))
    markers = {5: 'o', 8: 's', 14: '^'}

    for ax, (key, label) in zip(axes.flat, params):
        for eps in eps_vals:
            subset = [r for r in rows if r['epsilon_n'] == eps
                      and not math.isnan(r['mse'])]
            vals = [r[key] for r in subset]
            rms = [math.sqrt(r['mse']) for r in subset]
            colors = [QUALITY_COLORS[classify(r['mse'])] for r in subset]
            marker = markers.get(eps, 'o')
            ax.scatter(vals, rms, c=colors, marker=marker, s=40,
                       alpha=0.7, edgecolors='#2c3e50', linewidths=0.5,
                       label=f'ε_n={eps:.0f}')
        ax.set_xlabel(label)
        ax.set_ylabel('RMS')
        ax.set_yscale('log')

    # Legend for emittance markers (top-right panel)
    handles = [plt.Line2D([0], [0], marker=markers.get(e, 'o'), color='gray',
               linestyle='', markersize=6, label=f'ε_n={e:.0f}')
               for e in eps_vals]
    axes[0, 1].legend(handles=handles, fontsize=8, loc='upper right')

    fig.suptitle('P12: Twiss Parameters vs RMS', fontsize=13)
    fig.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(outdir / f'twiss_scatter.{ext}')
    plt.close(fig)
    _print(f"  Saved twiss_scatter.{{pdf,png}}")


def print_summary(rows):
    eps_vals = sorted(set(r['epsilon_n'] for r in rows))

    _print(f"\n{'=' * 90}")
    _print("  P12: MULTI-SEED ROBUSTNESS SUMMARY (RMS = √MSE)")
    _print(f"{'=' * 90}")
    _print(f"{'ε_n':>6} {'Seeds':>6} {'Best RMS':>12} {'Median':>12} "
           f"{'Worst':>12} {'CV%':>8} "
           f"{'Exc':>4} {'Acc':>4} {'Mar':>4} {'Fail':>4}")
    _print("-" * 90)

    for eps in eps_vals:
        mses = [r['mse'] for r in rows
                if r['epsilon_n'] == eps and not math.isnan(r['mse'])]
        if not mses:
            _print(f"{eps:6.0f} {'—':>6}")
            continue
        rms = [math.sqrt(m) for m in mses]
        cv = np.std(rms) / np.mean(rms) * 100 if np.mean(rms) > 0 else 0
        cats = {c: 0 for c in QUALITY_ORDER}
        for m in mses:
            cats[classify(m)] += 1
        _print(f"{eps:6.0f} {len(mses):6d} {min(rms):12.2e} "
               f"{np.median(rms):12.2e} {max(rms):12.2e} {cv:8.1f} "
               f"{cats['Excellent']:4d} {cats['Acceptable']:4d} "
               f"{cats['Marginal']:4d} {cats['Failed']:4d}")

    # Quad current spread summary
    _print(f"\n  Quad current CV% (mean across all quads):")
    for eps in eps_vals:
        subset = [r for r in rows
                  if r['epsilon_n'] == eps and not math.isnan(r['mse'])]
        if len(subset) < 2:
            continue
        cvs = []
        for idx in QUAD_INDICES:
            qk = f'quad_{idx}'
            vals = [r[qk] for r in subset if not math.isnan(r.get(qk, float('nan')))]
            if vals and np.mean(vals) != 0:
                cvs.append(np.std(vals) / abs(np.mean(vals)) * 100)
        if cvs:
            _print(f"    ε_n={eps:2.0f}: mean CV = {np.mean(cvs):.1f}%, "
                   f"max CV = {max(cvs):.1f}% (q{QUAD_INDICES[np.argmax(cvs)]})")


def main():
    parser = argparse.ArgumentParser(
        description='P12: Multi-Seed Robustness Study')
    parser.add_argument('--seeds', type=int, default=DEFAULT_N_SEEDS,
                        help='Number of seeds (default: 20)')
    parser.add_argument('--particles', type=int, default=DEFAULT_PARTICLES,
                        help='Particles per run (default: 500)')
    parser.add_argument('--emittances', type=float, nargs='+',
                        default=DEFAULT_EMITTANCES,
                        help='Emittance values (default: 5 8 14)')
    parser.add_argument('--plots-only', action='store_true',
                        help='Regenerate plots from cached CSV')
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    seeds = [42 + 100 * i for i in range(args.seeds)]

    if not args.plots_only:
        run_seed_study(seeds, args.emittances, args.particles)

    rows = load_csv(CSV_PATH)
    if not rows:
        _print("No data to plot. Run without --plots-only first.")
        return

    _print(f"\nGenerating plots from {len(rows)} data points...")
    plot_mse_boxplot(rows, OUTDIR)
    plot_quad_spread(rows, OUTDIR)
    plot_quality_histogram(rows, OUTDIR)
    plot_twiss_scatter(rows, OUTDIR)
    print_summary(rows)

    _print(f"\n{'=' * 60}")
    _print("  P12 Complete")
    _print(f"{'=' * 60}")


if __name__ == '__main__':
    main()
