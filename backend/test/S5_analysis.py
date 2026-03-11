"""S5 Analysis: 2D Coupled Parameter Scan Results

Analyzes the 300-point S5 dataset (3 × 10×10 grids) and generates:
  1. Summary statistics table
  2. Feasibility boundary comparison with S4 1D scans
  3. Combined analysis figure for report

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
from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from UHM_beamline_opt_05ps_params import MSE_THRESHOLDS

SCAN_DIR = Path(__file__).resolve().parent / 'results' / 'params_05ps_2d'
S4_DIR = Path(__file__).resolve().parent / 'results' / 'params_05ps'
OUTDIR = Path(__file__).resolve().parent / 'results' / 'S5_analysis'

SCAN_CONFIGS = {
    's5a': {
        'param1': 'energy_std_percent', 'param2': 'h',
        'p1_label': r'$\sigma_E$ (%)', 'p2_label': r'$h$ ($\times 10^9$/s)',
        'p2_display_scale': 1e-9,
        'fixed_param': 'epsilon_n', 'fixed_value': 8,
        'title': r'S5a: $\sigma_E \times h$',
    },
    's5b': {
        'param1': 'energy_std_percent', 'param2': 'epsilon_n',
        'p1_label': r'$\sigma_E$ (%)', 'p2_label': r'$\varepsilon_n$ ($\pi$·mm·mrad)',
        'p2_display_scale': 1.0,
        'fixed_param': 'h', 'fixed_value': 5e9,
        'title': r'S5b: $\sigma_E \times \varepsilon_n$',
    },
    's5c': {
        'param1': 'h', 'param2': 'epsilon_n',
        'p1_label': r'$h$ ($\times 10^9$/s)', 'p2_label': r'$\varepsilon_n$ ($\pi$·mm·mrad)',
        'p1_display_scale': 1e-9, 'p2_display_scale': 1.0,
        'fixed_param': 'energy_std_percent', 'fixed_value': 0.5,
        'title': r'S5c: $h \times \varepsilon_n$',
    },
}


def read_csv(filepath):
    with open(filepath) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


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


def print_summary():
    """Print comprehensive summary statistics."""
    print("\n" + "=" * 72)
    print("  S5: 2D Coupled Parameter Scan — Analysis Summary")
    print("=" * 72)

    for name, cfg in SCAN_CONFIGS.items():
        csv_path = SCAN_DIR / f'scan_{name}.csv'
        if not csv_path.exists():
            print(f"\n  {name.upper()}: no data")
            continue

        rows = read_csv(csv_path)
        mses = [r['mse'] for r in rows if not math.isnan(r['mse'])]

        cats = {'Excellent': 0, 'Acceptable': 0, 'Marginal': 0, 'Failed': 0}
        for r in rows:
            cats[classify(r['mse'])] += 1

        print(f"\n  {cfg['title']} (fixed {cfg['fixed_param']}={cfg['fixed_value']})")
        print(f"  {'─' * 50}")
        print(f"  Points: {len(rows)}")
        print(f"  MSE range: [{min(mses):.2e}, {max(mses):.2e}]")
        print(f"  Median MSE: {np.median(mses):.2e}")
        print(f"  Quality: {cats['Excellent']} Excellent, {cats['Acceptable']} Acceptable, "
              f"{cats['Marginal']} Marginal, {cats['Failed']} Failed")
        print(f"  Success rate: {(cats['Excellent'] + cats['Acceptable']) / len(rows) * 100:.0f}%")

        # Identify failure regions
        if cats['Failed'] > 0:
            p2 = cfg['param2']
            scale = cfg.get('p2_display_scale', 1.0)
            failed_p2 = sorted(set(r[p2] * scale for r in rows
                                   if classify(r['mse']) == 'Failed'))
            print(f"  Failed at {p2}: {[f'{v:.3g}' for v in failed_p2]}")


def generate_combined_figure():
    """3-panel feasibility map with S4 1D overlay."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    cmap = ListedColormap(['#2ecc71', '#f1c40f', '#e67e22', '#e74c3c'])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    quality_map = {'Excellent': 0, 'Acceptable': 1, 'Marginal': 2, 'Failed': 3}

    for ax, (name, cfg) in zip(axes, SCAN_CONFIGS.items()):
        csv_path = SCAN_DIR / f'scan_{name}.csv'
        if not csv_path.exists():
            ax.set_title(f'{name}: no data')
            continue

        rows = read_csv(csv_path)
        p1_key = cfg['param1']
        p2_key = cfg['param2']
        p1_scale = cfg.get('p1_display_scale', 1.0)
        p2_scale = cfg.get('p2_display_scale', 1.0)

        p1_vals = sorted(set(r[p1_key] for r in rows))
        p2_vals = sorted(set(r[p2_key] for r in rows))
        n1, n2 = len(p1_vals), len(p2_vals)

        quality_grid = np.full((n2, n1), 3.0)
        p1_map = {v: i for i, v in enumerate(p1_vals)}
        p2_map = {v: i for i, v in enumerate(p2_vals)}

        for r in rows:
            i = p1_map.get(r[p1_key])
            j = p2_map.get(r[p2_key])
            if i is not None and j is not None:
                quality_grid[j, i] = quality_map[classify(r['mse'])]

        p1_d = np.array(p1_vals) * p1_scale
        p2_d = np.array(p2_vals) * p2_scale

        im = ax.pcolormesh(p1_d, p2_d, quality_grid, cmap=cmap, norm=norm,
                           shading='nearest')
        ax.set_xlabel(cfg['p1_label'], fontsize=11)
        ax.set_ylabel(cfg['p2_label'], fontsize=11)
        ax.set_title(cfg['title'], fontsize=12)

    cbar = fig.colorbar(im, ax=axes.tolist(), ticks=[0, 1, 2, 3],
                        fraction=0.02, pad=0.04)
    cbar.set_ticklabels(['Excellent\n(<1e-3)', 'Acceptable\n(<0.01)',
                         'Marginal\n(<0.1)', 'Failed\n(≥0.1)'])

    fig.suptitle('S5: 2D Parameter Scan Feasibility Maps', fontsize=14, y=1.02)
    plt.tight_layout()
    for ext in ['eps', 'png', 'pdf']:
        fig.savefig(OUTDIR / f'S5_feasibility_maps.{ext}', dpi=200, bbox_inches='tight')
    print(f"  Saved: S5_feasibility_maps.{{eps,png,pdf}}")
    plt.close(fig)


def generate_mse_landscape():
    """3-panel MSE contour plot (log scale)."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    for ax, (name, cfg) in zip(axes, SCAN_CONFIGS.items()):
        csv_path = SCAN_DIR / f'scan_{name}.csv'
        if not csv_path.exists():
            continue

        rows = read_csv(csv_path)
        p1_key = cfg['param1']
        p2_key = cfg['param2']
        p1_scale = cfg.get('p1_display_scale', 1.0)
        p2_scale = cfg.get('p2_display_scale', 1.0)

        p1_vals = sorted(set(r[p1_key] for r in rows))
        p2_vals = sorted(set(r[p2_key] for r in rows))
        n1, n2 = len(p1_vals), len(p2_vals)

        mse_grid = np.full((n2, n1), np.nan)
        p1_map = {v: i for i, v in enumerate(p1_vals)}
        p2_map = {v: i for i, v in enumerate(p2_vals)}

        for r in rows:
            i = p1_map.get(r[p1_key])
            j = p2_map.get(r[p2_key])
            if i is not None and j is not None:
                mse_grid[j, i] = r['mse']

        p1_d = np.array(p1_vals) * p1_scale
        p2_d = np.array(p2_vals) * p2_scale

        mse_plot = np.where(np.isnan(mse_grid), 1e2, np.clip(mse_grid, 1e-10, 1e5))
        cs = ax.contourf(p1_d, p2_d, mse_plot,
                         levels=np.logspace(-9, 5, 30),
                         norm=LogNorm(vmin=1e-9, vmax=1e5),
                         cmap='viridis')
        plt.colorbar(cs, ax=ax, label='MSE')

        # Threshold contours
        for thresh, color in [(1e-3, 'lime'), (0.01, 'yellow'), (0.1, 'orange')]:
            try:
                ax.contour(p1_d, p2_d, mse_plot, levels=[thresh],
                          colors=[color], linewidths=2, linestyles='--')
            except Exception:
                pass

        ax.set_xlabel(cfg['p1_label'], fontsize=11)
        ax.set_ylabel(cfg['p2_label'], fontsize=11)
        ax.set_title(cfg['title'], fontsize=12)

    fig.suptitle('S5: MSE Landscape (log scale)', fontsize=14, y=1.02)
    plt.tight_layout()
    for ext in ['eps', 'png', 'pdf']:
        fig.savefig(OUTDIR / f'S5_mse_landscape.{ext}', dpi=200, bbox_inches='tight')
    print(f"  Saved: S5_mse_landscape.{{eps,png,pdf}}")
    plt.close(fig)


def s4_comparison():
    """Compare S5 marginals with S4 1D scans."""
    print("\n── S4 vs S5 Comparison ──")
    print("  S4 sweeps one parameter at a time (others fixed at baseline).")
    print("  S5 sweeps two parameters simultaneously.")
    print("  If S4 and S5 marginals agree, coupling effects are negligible.\n")

    # S5a: σ_E × h, fixed ε_n=8. S4 energy sweep has h=5e9 (baseline).
    # Compare S5a at h=5e9 (closest grid point) vs S4 energy scan.
    s5a_path = SCAN_DIR / 'scan_s5a.csv'
    s4_energy_path = S4_DIR / 'scan_energy_spread.csv'

    if s5a_path.exists() and s4_energy_path.exists():
        s5a = read_csv(s5a_path)
        s4e = read_csv(s4_energy_path)

        # S5a row where h is closest to 5e9
        h_vals = sorted(set(r['h'] for r in s5a))
        h_baseline = min(h_vals, key=lambda x: abs(x - 5e9))

        s5a_at_baseline_h = sorted(
            [r for r in s5a if r['h'] == h_baseline],
            key=lambda r: r['energy_std_percent']
        )

        print(f"  σ_E marginal comparison (h≈{h_baseline/1e9:.1f}×10⁹ vs S4 h=5×10⁹):")
        print(f"  {'σ_E':>8} {'S4 MSE':>12} {'S5a MSE':>12} {'Match?':>8}")
        for s5r in s5a_at_baseline_h:
            se = s5r['energy_std_percent']
            # Find closest S4 point
            closest_s4 = min(s4e, key=lambda r: abs(r['param_value'] - se))
            if abs(closest_s4['param_value'] - se) < 0.05:
                s4_cat = classify(closest_s4['mse'])
                s5_cat = classify(s5r['mse'])
                match = 'Yes' if s4_cat == s5_cat else 'NO'
                print(f"  {se:8.2f} {closest_s4['mse']:12.2e} {s5r['mse']:12.2e} {match:>8}")

    # S5b/S5c: check emittance marginals
    s5b_path = SCAN_DIR / 'scan_s5b.csv'
    s4_emit_path = S4_DIR / 'scan_emittance.csv'

    if s5b_path.exists() and s4_emit_path.exists():
        s5b = read_csv(s5b_path)
        s4em = read_csv(s4_emit_path)

        se_vals = sorted(set(r['energy_std_percent'] for r in s5b))
        se_baseline = min(se_vals, key=lambda x: abs(x - 0.5))

        s5b_at_baseline = sorted(
            [r for r in s5b if r['energy_std_percent'] == se_baseline],
            key=lambda r: r['epsilon_n']
        )

        print(f"\n  ε_n marginal comparison (σ_E≈{se_baseline}% vs S4 σ_E=0.5%):")
        print(f"  {'ε_n':>8} {'S4 MSE':>12} {'S5b MSE':>12} {'Match?':>8}")
        for s5r in s5b_at_baseline:
            en = s5r['epsilon_n']
            closest_s4 = min(s4em, key=lambda r: abs(r['param_value'] - en))
            if abs(closest_s4['param_value'] - en) < 1.0:
                s4_cat = classify(closest_s4['mse'])
                s5_cat = classify(s5r['mse'])
                match = 'Yes' if s4_cat == s5_cat else 'NO'
                print(f"  {en:8.2f} {closest_s4['mse']:12.2e} {s5r['mse']:12.2e} {match:>8}")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print_summary()
    generate_combined_figure()
    generate_mse_landscape()
    s4_comparison()

    print("\n" + "=" * 60)
    print("  S5 Analysis Complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
