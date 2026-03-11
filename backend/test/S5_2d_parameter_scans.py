"""S5: 2D Coupled Parameter Scans

Maps feasibility surfaces for three parameter pairs at the 0.5 ps operating point:
  S5a: σ_E × h    (energy spread vs chirp)    at ε_n = 8
  S5b: σ_E × ε_n  (energy spread vs emittance) at h = 5e9
  S5c: h × ε_n    (chirp vs emittance)         at σ_E = 0.5%

Each scan runs a 10×10 grid (100 optimizations) with checkpoint/resume.

Author: Eremey Valetov
"""

import sys
import time
import argparse
import csv
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from UHM_beamline_opt_05ps_params import (
    run_optimization, BASELINE, QUAD_INDICES, MSE_THRESHOLDS,
)

# ── Constants ─────────────────────────────────────────────────────────────────
OUTDIR = Path(__file__).resolve().parent / 'results' / 'params_05ps_2d'

SCAN_CONFIGS = {
    's5a': {
        'param1': 'energy_std_percent',
        'param2': 'h',
        'p1_label': r'$\sigma_E$ (%)',
        'p2_label': r'$h$ (10$^9$/s)',
        'p1_range': np.linspace(0.1, 3.0, 10),
        'p2_range': np.linspace(0, 30e9, 10),
        'p2_display_scale': 1e-9,  # convert to 10^9/s for display
        'fixed': {'epsilon_n': 8},
        'title': r'S5a: $\sigma_E$ vs $h$ ($\varepsilon_n$ = 8)',
    },
    's5b': {
        'param1': 'energy_std_percent',
        'param2': 'epsilon_n',
        'p1_label': r'$\sigma_E$ (%)',
        'p2_label': r'$\varepsilon_n$ ($\pi$·mm·mrad)',
        'p1_range': np.linspace(0.1, 3.0, 10),
        'p2_range': np.linspace(3, 20, 10),
        'p2_display_scale': 1.0,
        'fixed': {'h': 5e9},
        'title': r'S5b: $\sigma_E$ vs $\varepsilon_n$ ($h$ = 5×10$^9$/s)',
    },
    's5c': {
        'param1': 'h',
        'param2': 'epsilon_n',
        'p1_label': r'$h$ (10$^9$/s)',
        'p2_label': r'$\varepsilon_n$ ($\pi$·mm·mrad)',
        'p1_range': np.linspace(0, 30e9, 10),
        'p2_range': np.linspace(3, 20, 10),
        'p1_display_scale': 1e-9,
        'p2_display_scale': 1.0,
        'fixed': {'energy_std_percent': 0.5},
        'title': r'S5c: $h$ vs $\varepsilon_n$ ($\sigma_E$ = 0.5%)',
    },
}


def _print(msg):
    print(msg, flush=True)


# ── CSV I/O ───────────────────────────────────────────────────────────────────

def csv_header(param1_name, param2_name):
    """Column names for a 2D scan CSV."""
    cols = [param1_name, param2_name, 'mse', 'alpha_x', 'alpha_y', 'beta_x', 'beta_y',
            'disp_resid', 'time_s']
    for idx in QUAD_INDICES:
        cols.append(f'quad_{idx}')
    return cols


def result_to_row(p1_val, p2_val, res):
    """Convert optimization result to a CSV row."""
    row = [p1_val, p2_val, res['mse'], res['alpha_x'], res['alpha_y'],
           res['beta_x'], res['beta_y'], res['disp_resid'], res['time_s']]
    for idx in QUAD_INDICES:
        row.append(res['quad_currents'].get(idx, 0.0))
    return row


def write_csv(filepath, header, rows):
    """Write CSV file."""
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def read_csv(filepath):
    """Read CSV file and return list of dicts."""
    with open(filepath) as f:
        reader = csv.DictReader(f)
        rows = []
        for r in reader:
            rows.append({k: float(v) for k, v in r.items()})
        return rows


def load_completed(filepath, param1_name, param2_name):
    """Load completed (p1, p2) pairs from existing CSV."""
    completed = set()
    if filepath.exists():
        rows = read_csv(filepath)
        for r in rows:
            p1 = round(r[param1_name], 10)
            p2 = round(r[param2_name], 10)
            completed.add((p1, p2))
    return completed


# ── Core scan ─────────────────────────────────────────────────────────────────

def run_2d_scan(scan_name, nb_particles=500, seed=42, grid_size=None):
    """Run a 2D parameter scan with checkpoint/resume."""
    cfg = SCAN_CONFIGS[scan_name]
    param1 = cfg['param1']
    param2 = cfg['param2']
    p1_range = cfg['p1_range']
    p2_range = cfg['p2_range']
    fixed = cfg['fixed']

    if grid_size is not None:
        p1_range = np.linspace(p1_range[0], p1_range[-1], grid_size)
        p2_range = np.linspace(p2_range[0], p2_range[-1], grid_size)

    total = len(p1_range) * len(p2_range)
    csv_path = OUTDIR / f'scan_{scan_name}.csv'

    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Load checkpoint
    completed = load_completed(csv_path, param1, param2)
    _print(f"\n{'='*72}")
    _print(f"  {scan_name.upper()}: {cfg['title']}")
    _print(f"{'='*72}")
    _print(f"  {param1}: {len(p1_range)} values [{p1_range[0]:.4g} .. {p1_range[-1]:.4g}]")
    _print(f"  {param2}: {len(p2_range)} values [{p2_range[0]:.4g} .. {p2_range[-1]:.4g}]")
    _print(f"  Fixed: {fixed}")
    _print(f"  Total: {total} points, {len(completed)} already completed")

    if len(completed) == total:
        _print("  All points already computed — skipping")
        return read_csv(csv_path)

    # If fresh start, write header
    if not csv_path.exists():
        header = csv_header(param1, param2)
        write_csv(csv_path, header, [])

    header = csv_header(param1, param2)
    done = len(completed)
    t_start = time.time()

    for i, p1 in enumerate(p1_range):
        for j, p2 in enumerate(p2_range):
            key = (round(p1, 10), round(p2, 10))
            if key in completed:
                continue

            done += 1
            elapsed = time.time() - t_start
            rate = elapsed / (done - len(completed)) if done > len(completed) else 0
            remaining = (total - done) * rate
            eta_min = remaining / 60

            _print(f"  [{done}/{total}] {param1}={p1:.4g}, {param2}={p2:.4g} "
                   f"(ETA: {eta_min:.0f} min)")

            # Build kwargs
            kwargs = dict(BASELINE)
            kwargs.update(fixed)
            kwargs[param1] = p1
            kwargs[param2] = p2

            try:
                res = run_optimization(
                    nb_particles=nb_particles, seed=seed, **kwargs)

                row = result_to_row(p1, p2, res)

                # Append to CSV
                with open(csv_path, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(row)

                mse = res['mse']
                quality = 'Failed'
                for label, thresh in sorted(MSE_THRESHOLDS.items(),
                                            key=lambda x: x[1]):
                    if mse < thresh:
                        quality = label
                        break
                _print(f"         RMS = {math.sqrt(mse):.2e} ({quality}), "
                       f"t = {res['time_s']:.1f}s")

            except Exception as e:
                _print(f"         FAILED: {e}")
                # Write NaN row
                row = [p1, p2] + [float('nan')] * (len(header) - 2)
                with open(csv_path, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(row)

    _print(f"\n  Scan complete. Total time: {(time.time()-t_start)/60:.1f} min")
    _print(f"  Saved: {csv_path}")

    return read_csv(csv_path)


# ── Plotting ──────────────────────────────────────────────────────────────────

def classify_mse(mse):
    """Classify MSE into quality categories. Returns 0–3 index."""
    if mse < MSE_THRESHOLDS['Excellent']:
        return 0  # Excellent
    elif mse < MSE_THRESHOLDS['Acceptable']:
        return 1  # Acceptable
    elif mse < MSE_THRESHOLDS['Marginal']:
        return 2  # Marginal
    else:
        return 3  # Failed


def plot_2d_scan(scan_name, rows=None):
    """Generate contour plots for a 2D scan."""
    cfg = SCAN_CONFIGS[scan_name]
    csv_path = OUTDIR / f'scan_{scan_name}.csv'

    if rows is None:
        if not csv_path.exists():
            _print(f"  No data for {scan_name} — skipping plots")
            return
        rows = read_csv(csv_path)

    param1 = cfg['param1']
    param2 = cfg['param2']
    p1_scale = cfg.get('p1_display_scale', 1.0)
    p2_scale = cfg.get('p2_display_scale', 1.0)

    # Extract unique grid values
    p1_vals = sorted(set(r[param1] for r in rows))
    p2_vals = sorted(set(r[param2] for r in rows))

    n1, n2 = len(p1_vals), len(p2_vals)
    if n1 < 2 or n2 < 2:
        _print(f"  Insufficient data for contour plot ({n1}×{n2})")
        return

    # Build 2D arrays
    mse_grid = np.full((n2, n1), np.nan)
    quality_grid = np.full((n2, n1), 3)  # default = Failed
    beta_x_dev = np.full((n2, n1), np.nan)
    beta_y_dev = np.full((n2, n1), np.nan)

    p1_map = {v: i for i, v in enumerate(p1_vals)}
    p2_map = {v: i for i, v in enumerate(p2_vals)}

    for r in rows:
        i = p1_map.get(r[param1])
        j = p2_map.get(r[param2])
        if i is not None and j is not None:
            mse = r['mse']
            mse_grid[j, i] = mse
            quality_grid[j, i] = classify_mse(mse) if not np.isnan(mse) else 3
            beta_x_dev[j, i] = abs(r['beta_x'] - 1.4) / 1.4 * 100 if not np.isnan(r['beta_x']) else np.nan
            beta_y_dev[j, i] = abs(r['beta_y'] - 0.2418) / 0.2418 * 100 if not np.isnan(r['beta_y']) else np.nan

    p1_display = np.array(p1_vals) * p1_scale
    p2_display = np.array(p2_vals) * p2_scale

    # ── Figure 1: RMS contour ──
    fig, ax = plt.subplots(figsize=(8, 6))

    mse_plot = np.where(np.isnan(mse_grid), 1e2, mse_grid)
    mse_plot = np.clip(mse_plot, 1e-10, 1e2)
    rms_plot = np.sqrt(mse_plot)

    cs = ax.contourf(p1_display, p2_display, rms_plot,
                     levels=np.logspace(-5, 1, 30),
                     norm=LogNorm(vmin=1e-5, vmax=1e1),
                     cmap='viridis')
    cbar = plt.colorbar(cs, ax=ax, label='RMS Twiss Mismatch')

    # Feasibility boundaries (sqrt of MSE thresholds)
    for thresh_val, (thresh_name, _) in zip(
            [math.sqrt(1e-3), math.sqrt(0.01), math.sqrt(0.1)],
            [('Excellent', 'lime'), ('Acceptable', 'yellow'), ('Marginal', 'orange')]):
        try:
            ct = ax.contour(p1_display, p2_display, rms_plot,
                           levels=[thresh_val], colors=[_[1]], linewidths=2)
            ax.clabel(ct, fmt={thresh_val: thresh_name}, fontsize=8)
        except Exception:
            pass

    ax.set_xlabel(cfg['p1_label'])
    ax.set_ylabel(cfg['p2_label'])
    ax.set_title(cfg['title'])

    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'{scan_name}_mse_landscape.{ext}', dpi=150)
    _print(f"  Saved: {scan_name}_mse_landscape.{{eps,png}}")
    plt.close(fig)

    # ── Figure 2: Twiss deviation contour ──
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for ax, data, label in [(ax1, beta_x_dev, r'$|\Delta\beta_x/\beta_{x,\mathrm{tgt}}|$ (%)'),
                             (ax2, beta_y_dev, r'$|\Delta\beta_y/\beta_{y,\mathrm{tgt}}|$ (%)')]:
        data_plot = np.where(np.isnan(data), 100, data)
        cs = ax.contourf(p1_display, p2_display, data_plot,
                        levels=np.linspace(0, 100, 20), cmap='RdYlGn_r')
        plt.colorbar(cs, ax=ax, label=label)
        ax.set_xlabel(cfg['p1_label'])
        ax.set_ylabel(cfg['p2_label'])
        ax.set_title(label)

    fig.suptitle(cfg['title'] + ' — Twiss deviation', fontsize=12)
    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'{scan_name}_twiss_deviation.{ext}', dpi=150)
    _print(f"  Saved: {scan_name}_twiss_deviation.{{eps,png}}")
    plt.close(fig)

    # ── Figure 3: Feasibility map (discrete quality categories) ──
    fig, ax = plt.subplots(figsize=(8, 6))

    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap = ListedColormap(['#2ecc71', '#f1c40f', '#e67e22', '#e74c3c'])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    im = ax.pcolormesh(p1_display, p2_display, quality_grid.astype(float),
                       cmap=cmap, norm=norm, shading='nearest')
    cbar = plt.colorbar(im, ax=ax, ticks=[0, 1, 2, 3])
    cbar.set_ticklabels([f'Excellent\n(<{math.sqrt(1e-3):.2e})',
                         f'Acceptable\n(<{math.sqrt(0.01):.2e})',
                         f'Marginal\n(<{math.sqrt(0.1):.2e})',
                         f'Failed\n(\u2265{math.sqrt(0.1):.2e})'])
    ax.set_xlabel(cfg['p1_label'])
    ax.set_ylabel(cfg['p2_label'])
    ax.set_title(cfg['title'] + ' — Feasibility')

    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'{scan_name}_feasibility.{ext}', dpi=150)
    _print(f"  Saved: {scan_name}_feasibility.{{eps,png}}")
    plt.close(fig)


def plot_combined_feasibility():
    """Combined 1×3 feasibility summary figure."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap = ListedColormap(['#2ecc71', '#f1c40f', '#e67e22', '#e74c3c'])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    for ax, scan_name in zip(axes, ['s5a', 's5b', 's5c']):
        cfg = SCAN_CONFIGS[scan_name]
        csv_path = OUTDIR / f'scan_{scan_name}.csv'

        if not csv_path.exists():
            ax.set_title(f'{scan_name}: No data')
            continue

        rows = read_csv(csv_path)
        param1 = cfg['param1']
        param2 = cfg['param2']
        p1_scale = cfg.get('p1_display_scale', 1.0)
        p2_scale = cfg.get('p2_display_scale', 1.0)

        p1_vals = sorted(set(r[param1] for r in rows))
        p2_vals = sorted(set(r[param2] for r in rows))
        n1, n2 = len(p1_vals), len(p2_vals)

        quality_grid = np.full((n2, n1), 3)
        p1_map = {v: i for i, v in enumerate(p1_vals)}
        p2_map = {v: i for i, v in enumerate(p2_vals)}

        for r in rows:
            i = p1_map.get(r[param1])
            j = p2_map.get(r[param2])
            if i is not None and j is not None:
                mse = r['mse']
                quality_grid[j, i] = classify_mse(mse) if not np.isnan(mse) else 3

        p1_display = np.array(p1_vals) * p1_scale
        p2_display = np.array(p2_vals) * p2_scale

        im = ax.pcolormesh(p1_display, p2_display, quality_grid.astype(float),
                          cmap=cmap, norm=norm, shading='nearest')
        ax.set_xlabel(cfg['p1_label'])
        ax.set_ylabel(cfg['p2_label'])
        ax.set_title(scan_name.upper(), fontsize=11)

        # Mark baseline point
        p1_base = BASELINE.get(param1, None)
        p2_base = BASELINE.get(param2, None)
        if p1_base is not None and p2_base is not None:
            ax.plot(p1_base * p1_scale, p2_base * p2_scale, 'w*',
                   markersize=15, markeredgecolor='k', markeredgewidth=1)

    # Shared colorbar
    cbar = fig.colorbar(im, ax=axes.tolist(), ticks=[0, 1, 2, 3],
                       fraction=0.02, pad=0.04)
    cbar.set_ticklabels(['Excellent', 'Acceptable', 'Marginal', 'Failed'])

    fig.suptitle('S5: 2D Parameter Scan Feasibility Summary', fontsize=13, y=1.02)
    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'S5_feasibility_summary.{ext}', dpi=150, bbox_inches='tight')
    _print(f"  Saved: S5_feasibility_summary.{{eps,png}}")
    plt.close(fig)


def generate_all_plots():
    """Generate all plots from existing CSV data."""
    _print("\n── Generating plots ──")
    for scan_name in ['s5a', 's5b', 's5c']:
        csv_path = OUTDIR / f'scan_{scan_name}.csv'
        if csv_path.exists():
            plot_2d_scan(scan_name)
        else:
            _print(f"  {scan_name}: no data")
    plot_combined_feasibility()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='S5: 2D coupled parameter scans')
    parser.add_argument('--s5a', action='store_true',
                        help='Run S5a (σ_E × h) scan')
    parser.add_argument('--s5b', action='store_true',
                        help='Run S5b (σ_E × ε_n) scan')
    parser.add_argument('--s5c', action='store_true',
                        help='Run S5c (h × ε_n) scan')
    parser.add_argument('--all', action='store_true',
                        help='Run all three scans')
    parser.add_argument('--plots-only', action='store_true',
                        help='Regenerate plots from existing data')
    parser.add_argument('--particles', type=int, default=500,
                        help='Particles per optimization (default: 500)')
    parser.add_argument('--grid', type=int, default=None,
                        help='Grid size per axis (default: 10)')
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)

    if args.plots_only:
        generate_all_plots()
        return

    if not any([args.s5a, args.s5b, args.s5c, args.all]):
        args.all = True

    _print("S5: 2D Coupled Parameter Scans")
    _print(f"Particles: {args.particles}, Grid: {args.grid or 10}")

    if args.s5a or args.all:
        rows = run_2d_scan('s5a', nb_particles=args.particles, grid_size=args.grid)
        plot_2d_scan('s5a', rows)

    if args.s5b or args.all:
        rows = run_2d_scan('s5b', nb_particles=args.particles, grid_size=args.grid)
        plot_2d_scan('s5b', rows)

    if args.s5c or args.all:
        rows = run_2d_scan('s5c', nb_particles=args.particles, grid_size=args.grid)
        plot_2d_scan('s5c', rows)

    plot_combined_feasibility()

    _print("\n" + "=" * 72)
    _print("  S5 Complete")
    _print("=" * 72)


if __name__ == "__main__":
    main()
