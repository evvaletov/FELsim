"""P13: Deterministic Beam Generation — Sobol vs Pseudo-Random

P12 showed that single-seed results are not representative: 50% failure at
epsilon_n=5, CV=200-400%. The variability originates from upstream stages 1-10
(random beam generation), not Stage 11 (which is unimodal per S8).

This study replaces the pseudo-random Gaussian beam with a Sobol quasi-random
sequence (inverse normal CDF). The Sobol distribution:
  - Is deterministic: identical output regardless of seed
  - Has better space-filling properties than pseudo-random
  - Eliminates the seed-dependent variability in stages 1-10

The study runs the same 20 seeds × 3 emittances as P12, but with method='sobol'.
If the hypothesis is correct, all seeds should produce identical (or near-identical)
results, confirming that upstream variability was purely from random sampling.

Author: Eremey Valetov
"""

import sys
import csv
import math
import argparse
import time
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

OUTDIR = Path(__file__).resolve().parent / 'results' / 'P13'

DEFAULT_EMITTANCES = [5, 8, 14]
DEFAULT_PARTICLES = 512  # Power of 2 for Sobol
DEFAULT_N_SEEDS = 20

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
    cols = ['seed', 'epsilon_n', 'beam_method', 'mse', 'alpha_x', 'alpha_y',
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


def run_comparison(seeds, emittances, n_particles):
    """Run each (seed, emittance) with both random and sobol beams."""
    OUTDIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUTDIR / 'sobol_vs_random.csv'
    header = csv_header()

    completed = set()
    if csv_path.exists():
        with open(csv_path) as f:
            for r in csv.DictReader(f):
                completed.add((int(float(r['seed'])), float(r['epsilon_n']),
                               r['beam_method']))
    else:
        with open(csv_path, 'w', newline='') as f:
            csv.writer(f).writerow(header)

    methods = ['sobol', 'random']
    total = len(seeds) * len(emittances) * len(methods)
    done = len(completed)

    _print(f"P13: Deterministic Beam — Sobol vs Pseudo-Random")
    _print(f"  {len(emittances)} emittances × {len(seeds)} seeds × {len(methods)} methods = {total} runs")
    _print(f"  {n_particles} particles, cold-start")
    _print(f"  Already completed: {done}")

    for eps_n in emittances:
        _print(f"\n── ε_n = {eps_n} ──")
        for method in methods:
            _print(f"  Method: {method}")
            for seed in seeds:
                if (seed, float(eps_n), method) in completed:
                    continue
                done += 1
                _print(f"    [{done}/{total}] seed={seed}", end='')

                kwargs = dict(BASELINE)
                kwargs['epsilon_n'] = eps_n
                kwargs['nb_particles'] = n_particles
                kwargs['seed'] = seed
                kwargs['beam_method'] = method

                try:
                    res = run_optimization(**kwargs)
                    rms = math.sqrt(res['mse'])
                    row = [seed, eps_n, method, res['mse'],
                           res['alpha_x'], res['alpha_y'],
                           res['beta_x'], res['beta_y'], res['disp_resid']]
                    row += [res['quad_currents'][idx] for idx in QUAD_INDICES]
                    row += [res['time_s'], res.get('nfev') or 0,
                            int(res.get('converged') or False)]
                    _print(f" RMS={rms:.2e} ({classify(res['mse'])})")
                except Exception as e:
                    _print(f" FAILED: {e}")
                    row = [seed, eps_n, method] + [float('nan')] * (len(header) - 3)

                with open(csv_path, 'a', newline='') as f:
                    csv.writer(f).writerow(row)


def print_summary(rows):
    """Print comparison table: Sobol vs Random at each emittance."""
    eps_vals = sorted(set(r['epsilon_n'] for r in rows))

    _print("\n" + "=" * 72)
    _print("P13 Summary: Sobol vs Pseudo-Random Beam Generation")
    _print("=" * 72)

    for eps in eps_vals:
        _print(f"\n── ε_n = {eps:.0f} ──")
        for method in ['sobol', 'random']:
            sub = [r for r in rows if r['epsilon_n'] == eps
                   and r.get('beam_method') == method]
            if not sub:
                continue
            mses = [r['mse'] for r in sub if not math.isnan(r['mse'])]
            if not mses:
                _print(f"  {method:7s}: no valid results")
                continue
            rmses = [math.sqrt(m) for m in mses]
            quals = [classify(m) for m in mses]
            n_exc = quals.count('Excellent')
            n_fail = quals.count('Failed')
            med = np.median(rmses)
            cv = np.std(rmses) / np.mean(rmses) * 100 if np.mean(rmses) > 0 else 0

            _print(f"  {method:7s}: {len(mses):2d} runs, "
                   f"median RMS={med:.2e}, CV={cv:.0f}%, "
                   f"Exc={n_exc}, Fail={n_fail}")

            # For Sobol: check if all results are identical
            if method == 'sobol' and len(mses) > 1:
                spread = max(rmses) - min(rmses)
                if spread < 1e-12:
                    _print(f"          → ALL IDENTICAL (spread={spread:.1e})")
                elif spread < 1e-6:
                    _print(f"          → Near-identical (spread={spread:.1e})")
                else:
                    _print(f"          → Variable (spread={spread:.1e})")


def plot_comparison(rows, outdir):
    """Side-by-side box plots: Sobol vs Random."""
    eps_vals = sorted(set(r['epsilon_n'] for r in rows))

    fig, axes = plt.subplots(1, len(eps_vals), figsize=(4 * len(eps_vals), 5),
                             sharey=True)
    if len(eps_vals) == 1:
        axes = [axes]

    colors = {'sobol': '#2ecc71', 'random': '#3498db'}

    for ax, eps in zip(axes, eps_vals):
        data, labels, box_colors = [], [], []
        for method in ['sobol', 'random']:
            sub = [r for r in rows if r['epsilon_n'] == eps
                   and r.get('beam_method') == method]
            mses = [r['mse'] for r in sub if not math.isnan(r['mse'])]
            rmses = [math.sqrt(m) for m in mses] if mses else [float('nan')]
            data.append(rmses)
            labels.append(method)
            box_colors.append(colors[method])

        bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, widths=0.5)
        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)

        # Overlay individual points
        jitter_rng = np.random.RandomState(0)
        for i, d in enumerate(data):
            x = jitter_rng.normal(i + 1, 0.04, len(d))
            ax.scatter(x, d, alpha=0.5, s=15, c='#2c3e50', zorder=3)

        ax.set_yscale('log')
        ax.set_title(f'$\\varepsilon_n = {eps:.0f}$')
        ax.set_ylabel('RMS' if ax == axes[0] else '')

        for name, thresh in MSE_THRESHOLDS.items():
            rms_thresh = math.sqrt(thresh)
            ax.axhline(rms_thresh, ls=':', color='gray', alpha=0.5)

    fig.suptitle('P13: Sobol vs Pseudo-Random Beam (20 seeds)')
    fig.tight_layout()
    fig.savefig(outdir / 'sobol_vs_random.pdf')
    fig.savefig(outdir / 'sobol_vs_random.png', dpi=150)
    plt.close(fig)
    _print(f"  Plot saved: {outdir / 'sobol_vs_random.pdf'}")


def main():
    parser = argparse.ArgumentParser(
        description='P13: Deterministic beam generation — Sobol vs pseudo-random')
    parser.add_argument('--seeds', type=int, default=DEFAULT_N_SEEDS,
                        help=f'Number of seeds (default {DEFAULT_N_SEEDS})')
    parser.add_argument('--particles', type=int, default=DEFAULT_PARTICLES,
                        help=f'Particles per run (default {DEFAULT_PARTICLES}, power of 2)')
    parser.add_argument('--emittances', type=float, nargs='+',
                        default=DEFAULT_EMITTANCES)
    parser.add_argument('--plots-only', action='store_true',
                        help='Regenerate plots from existing CSV')
    parser.add_argument('--sobol-only', action='store_true',
                        help='Run only Sobol beam (skip random comparison)')
    args = parser.parse_args()

    seeds = [42 + 100 * i for i in range(args.seeds)]

    if not args.plots_only:
        if args.sobol_only:
            # Quick validation: run 3 seeds with Sobol only
            _print("P13 quick validation: Sobol determinism check")
            OUTDIR.mkdir(parents=True, exist_ok=True)
            for eps_n in args.emittances:
                results = []
                for seed in seeds[:3]:
                    kwargs = dict(BASELINE)
                    kwargs['epsilon_n'] = eps_n
                    kwargs['nb_particles'] = args.particles
                    kwargs['seed'] = seed
                    kwargs['beam_method'] = 'sobol'
                    res = run_optimization(**kwargs)
                    rms = math.sqrt(res['mse'])
                    _print(f"  ε_n={eps_n}, seed={seed}: "
                           f"RMS={rms:.6e} ({classify(res['mse'])})")
                    results.append(res)

                # Check if results are identical
                mses = [r['mse'] for r in results]
                spread = max(mses) - min(mses)
                _print(f"  ε_n={eps_n}: MSE spread across seeds = {spread:.1e}"
                       f" {'(IDENTICAL)' if spread < 1e-15 else ''}")
            return

        run_comparison(seeds, args.emittances, args.particles)

    csv_path = OUTDIR / 'sobol_vs_random.csv'
    rows = load_csv(csv_path)
    if rows:
        print_summary(rows)
        plot_comparison(rows, OUTDIR)
    else:
        _print("No results found. Run without --plots-only first.")


if __name__ == '__main__':
    main()
