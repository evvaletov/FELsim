"""R2: Unified Comparison Table Across All Studies

Aggregates results from W4, S4, W8, W9, W10, W11, W12 into comparison tables
and summary plots. Generates LaTeX tables and EPS figures.

Data Sources:
  W4  — COSY FR0/FR1 cross-validation
  S4  — 1D parameter sensitivity scans
  W8  — RF-Track cross-validation & optimization
  W9  — COSY longitudinal study (R56, T566)
  W10 — Beam losses & transmission
  W11 — Throughput optimization
  W12 — Bunch compression feasibility

Author: Eremey Valetov
"""

import sys
import json
import csv
import math
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
sys.path.insert(0, str(Path(__file__).resolve().parent))

# ── Paths ─────────────────────────────────────────────────────────────────────
RESULTS = Path(__file__).resolve().parent / 'results'
OUTDIR = RESULTS / 'R2'

MSE_THRESHOLDS = {
    'Excellent': 1e-3,
    'Acceptable': 0.01,
    'Marginal': 0.1,
}


def _print(msg='', **kwargs):
    print(msg, flush=True, **kwargs)


def _sanitize(obj):
    """JSON serializer for NaN and numpy types."""
    if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
        return None
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def classify_mse(mse):
    if mse < MSE_THRESHOLDS['Excellent']:
        return 'Excellent'
    elif mse < MSE_THRESHOLDS['Acceptable']:
        return 'Acceptable'
    elif mse < MSE_THRESHOLDS['Marginal']:
        return 'Marginal'
    return 'Failed'


def load_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        _print(f"  WARNING: {path} not found")
        return None


def load_csv(path):
    try:
        with open(path) as f:
            return list(csv.DictReader(f))
    except FileNotFoundError:
        _print(f"  WARNING: {path} not found")
        return None


# ═══════════════════════════════════════════════════════════════════════════════
#  Table 1: Baseline Optimization Cross-Code Comparison
# ═══════════════════════════════════════════════════════════════════════════════

def table_1_baseline():
    """FELsim / COSY FR0 / COSY FR1 / RFT-opt comparison at ε_n=8."""
    _print("\n" + "=" * 72)
    _print("  Table 1: Baseline Optimization — Cross-Code Comparison (ε_n = 8)")
    _print("=" * 72)

    rows = []

    # FELsim baseline (from S4 scan at σ_E=0.5%)
    s4_energy = load_csv(RESULTS / 'params_05ps' / 'scan_energy_spread.csv')
    if s4_energy:
        for r in s4_energy:
            if abs(float(r['param_value']) - 0.5) < 0.01:
                rows.append({
                    'Code': 'FELsim',
                    'Model': 'Transfer matrix',
                    'MSE': float(r['mse']),
                    'β_x': float(r['beta_x']),
                    'α_x': float(r['alpha_x']),
                    'β_y': float(r['beta_y']),
                    'α_y': float(r['alpha_y']),
                })
                break

    # COSY FR0
    fr0 = load_json(RESULTS / 'cosy_s1_fr0.json')
    if fr0:
        t = fr0.get('twiss_undulator', {})
        rows.append({
            'Code': 'COSY FR0',
            'Model': 'DA map (no fringe)',
            'MSE': fr0['mse'],
            'β_x': t.get('beta_x', 0),
            'α_x': t.get('alpha_x', 0),
            'β_y': t.get('beta_y', 0),
            'α_y': t.get('alpha_y', 0),
        })

    # COSY FR1
    fr1 = load_json(RESULTS / 'cosy_s1_fr1_warm.json')
    if fr1:
        t = fr1.get('twiss_undulator', {})
        rows.append({
            'Code': 'COSY FR1',
            'Model': 'DA map (1st-order fringe)',
            'MSE': fr1['mse'],
            'β_x': t.get('beta_x', 0),
            'α_x': t.get('alpha_x', 0),
            'β_y': t.get('beta_y', 0),
            'α_y': t.get('alpha_y', 0),
        })

    # RFT-opt at ε_n=8
    w8 = load_csv(RESULTS / 'rftrack_opt' / 'comparison.csv')
    if w8:
        for r in w8:
            if abs(float(r['epsilon_n']) - 8) < 0.1 and r['method'] == 'RFT-opt':
                rows.append({
                    'Code': 'RF-Track opt',
                    'Model': 'Particle tracking',
                    'MSE': float(r['mse']),
                    'β_x': float(r['beta_x']),
                    'α_x': float(r['alpha_x']),
                    'β_y': float(r['beta_y']),
                    'α_y': float(r['alpha_y']),
                })
                break

    # Print table
    if rows:
        _print(f"\n  {'Code':<16s} {'Model':<24s} {'RMS':>10s} {'β_x':>8s} {'α_x':>8s} "
               f"{'β_y':>8s} {'α_y':>8s} {'Quality':<12s}")
        _print("  " + "-" * 94)
        for r in rows:
            q = classify_mse(r['MSE'])
            _print(f"  {r['Code']:<16s} {r['Model']:<24s} {math.sqrt(r['MSE']):10.2e} {r['β_x']:8.4f} "
                   f"{r['α_x']:8.4f} {r['β_y']:8.4f} {r['α_y']:8.4f} {q:<12s}")

        _print(f"\n  Targets: β_x = 1.400 m, α_x = 0.470, β_y = 0.2418 m, α_y ≈ 0")

    return rows


# ═══════════════════════════════════════════════════════════════════════════════
#  Table 2: Parameter Sensitivity Summary
# ═══════════════════════════════════════════════════════════════════════════════

def table_2_sensitivity():
    """Summarize S4 scan results by quality category."""
    _print("\n" + "=" * 72)
    _print("  Table 2: Parameter Sensitivity Summary (S4)")
    _print("=" * 72)

    scans = {
        'Energy spread (%)': RESULTS / 'params_05ps' / 'scan_energy_spread.csv',
        'Chirp (1/s)': RESULTS / 'params_05ps' / 'scan_chirp.csv',
        'Emittance (π·mm·mrad)': RESULTS / 'params_05ps' / 'scan_emittance.csv',
    }

    rows = []
    for name, path in scans.items():
        data = load_csv(path)
        if not data:
            continue

        counts = {'Excellent': 0, 'Acceptable': 0, 'Marginal': 0, 'Failed': 0}
        mse_values = []
        for r in data:
            mse = float(r['mse'])
            mse_values.append(mse)
            counts[classify_mse(mse)] += 1

        total = len(data)
        best_mse = min(mse_values)
        worst_mse = max(mse_values)
        param_range = f"{float(data[0]['param_value']):.2g} – {float(data[-1]['param_value']):.2g}"

        rows.append({
            'Scan': name,
            'Range': param_range,
            'N': total,
            **counts,
            'Best MSE': best_mse,
            'Worst MSE': worst_mse,
        })

        _print(f"\n  {name} ({param_range}, {total} points):")
        _print(f"    Excellent: {counts['Excellent']}, Acceptable: {counts['Acceptable']}, "
               f"Marginal: {counts['Marginal']}, Failed: {counts['Failed']}")
        _print(f"    Best MSE: {best_mse:.2e}, Worst MSE: {worst_mse:.2e}")

    return rows


# ═══════════════════════════════════════════════════════════════════════════════
#  Table 3: Bunch Length & Transmission Cross-Code
# ═══════════════════════════════════════════════════════════════════════════════

def table_3_bunch_transmission():
    """Cross-code comparison of σ_t and transmission."""
    _print("\n" + "=" * 72)
    _print("  Table 3: Bunch Length & Transmission — Cross-Code")
    _print("=" * 72)

    rows = []

    # W10 Part A: Transmission baseline
    w10a = load_json(RESULTS / 'W10' / 'part_a_results.json')
    if w10a:
        for entry in w10a:
            rows.append({
                'Source': 'W10',
                'Label': entry.get('label', '?'),
                'σ_t_in (ps)': entry.get('sigma_t_ps_initial', '?'),
                'σ_t_out (ps)': entry.get('sigma_z_ps', '?'),
                'Transmission': entry.get('transmission', '?'),
                'I_peak (A)': entry.get('I_peak_A', '?'),
            })

    # W11 comparison CSV
    w11 = load_csv(RESULTS / 'W11' / 'W11_comparison.csv')
    if w11:
        for r in w11:
            rows.append({
                'Source': 'W11',
                'Label': f"{r['scenario']}/{r['label']}/{r['code']}",
                'σ_t_in (ps)': '2.0' if '2ps' in r['scenario'] else '0.5',
                'σ_t_out (ps)': r['sigma_t_ps'],
                'Transmission': r['transmission'],
                'I_peak (A)': r['I_peak_A'],
            })

    if rows:
        _print(f"\n  {'Source':<6s} {'Label':<45s} {'σ_t in':>8s} {'σ_t out':>8s} "
               f"{'T':>6s} {'I_pk':>8s}")
        _print("  " + "-" * 86)
        for r in rows:
            _print(f"  {r['Source']:<6s} {r['Label']:<45s} {str(r['σ_t_in (ps)']):>8s} "
                   f"{str(r['σ_t_out (ps)']):>8s} {str(r['Transmission']):>6s} "
                   f"{str(r['I_peak (A)']):>8s}")

    return rows


# ═══════════════════════════════════════════════════════════════════════════════
#  Table 4: Compression Feasibility
# ═══════════════════════════════════════════════════════════════════════════════

def table_4_compression():
    """Compression results from W12."""
    _print("\n" + "=" * 72)
    _print("  Table 4: Compression Feasibility (W12)")
    _print("=" * 72)

    w12a = load_json(RESULTS / 'W12' / 'part_a_results.json')
    if not w12a:
        return []

    R56 = w12a.get('R56_m', 0)
    T566 = w12a.get('T566_m', 0)
    h_opt = w12a.get('h_opt', 0)
    floor = w12a.get('compression_floor_ps', 0)

    _print(f"\n  R56 = {R56*1e3:.2f} mm, T566 = {T566:.2e} m")
    _print(f"  Compression floor = {floor:.3f} ps")
    _print(f"  Optimal chirp = {h_opt:.3e} 1/s")

    sweep = w12a.get('sweep', [])
    key_points = []

    for pt in sweep:
        h = pt.get('h', 0)
        # Select key chirp values
        if abs(h) < 0.1 or abs(h - (-4.2e9)) < 1e8 or abs(h - (-8.3e9)) < 1e8 or \
           abs(h - (-1.1e10)) < 2e8:
            key_points.append(pt)

    if key_points:
        _print(f"\n  {'h (1/s)':>12s} {'φ (deg)':>8s} {'C_map':>6s} {'σ_z out':>8s} "
               f"{'σ_δ_eff':>8s} {'I_peak':>8s}")
        _print("  " + "-" * 58)
        for pt in key_points:
            _print(f"  {pt.get('h', 0):12.2e} {pt.get('phi_deg', 0):8.1f} "
                   f"{pt.get('C_map', 0):6.2f} {pt.get('sigma_z_out_ps', 0):8.3f} "
                   f"{pt.get('sigma_delta_eff_pct', 0):7.3f}% {pt.get('I_peak_out', 0):8.1f}")

    # W12 Part B: RF-Track validation
    w12b = load_json(RESULTS / 'W12' / 'part_b_results.json')
    if w12b:
        _print(f"\n  RF-Track validation (W12 Part B):")
        _print(f"  {'Label':<20s} {'map_C':>6s} {'rft_C':>6s} {'T':>6s} {'I_pk':>8s}")
        _print("  " + "-" * 50)
        for entry in w12b:
            _print(f"  {entry.get('label', '?'):<20s} {entry.get('map_C', 0):6.2f} "
                   f"{entry.get('rft_C', 0):6.3f} {entry.get('rft_transmission', 0):6.3f} "
                   f"{entry.get('rft_I_peak', 0):8.2f}")

    return {'w12a': w12a, 'w12b': w12b}


# ═══════════════════════════════════════════════════════════════════════════════
#  Table 5: Quad Current Comparison
# ═══════════════════════════════════════════════════════════════════════════════

def table_5_currents():
    """Compare quad currents across key configurations."""
    _print("\n" + "=" * 72)
    _print("  Table 5: Quad Current Comparison (A)")
    _print("=" * 72)

    configs = {}

    # FELsim S4 baseline (σ_E=0.5%, ε_n=8)
    s4 = load_csv(RESULTS / 'params_05ps' / 'scan_energy_spread.csv')
    if s4:
        for r in s4:
            if abs(float(r['param_value']) - 0.5) < 0.01:
                currents = {}
                for k, v in r.items():
                    if k.startswith('quad_'):
                        idx = int(k.split('_')[1])
                        currents[idx] = float(v)
                configs['FELsim'] = currents
                break

    # COSY FR0
    fr0 = load_json(RESULTS / 'cosy_s1_fr0.json')
    if fr0:
        configs['COSY FR0'] = {int(k): float(v) for k, v in fr0['currents'].items()}

    # COSY FR1
    fr1 = load_json(RESULTS / 'cosy_s1_fr1_warm.json')
    if fr1:
        configs['COSY FR1'] = {int(k): float(v) for k, v in fr1['currents'].items()}

    # W9 optimized
    w9 = load_json(RESULTS / 'W9' / 'part_a_longitudinal_map.json')
    if w9:
        configs['W9 COSY'] = {int(k): float(v) for k, v in w9['currents'].items()}

    if not configs:
        return {}

    # Get all quad indices
    all_indices = sorted(set().union(*[set(c.keys()) for c in configs.values()]))
    config_names = list(configs.keys())

    _print(f"\n  {'Quad':>6s}", end='')
    for name in config_names:
        _print(f"  {name:>12s}", end='')
    _print()
    _print("  " + "-" * (8 + 14 * len(config_names)))

    for idx in all_indices:
        _print(f"  {idx:>6d}", end='')
        for name in config_names:
            val = configs[name].get(idx, 0)
            _print(f"  {val:12.4f}", end='')
        _print()

    return configs


# ═══════════════════════════════════════════════════════════════════════════════
#  Summary Plots
# ═══════════════════════════════════════════════════════════════════════════════

def plot_mse_comparison():
    """3-panel RMS vs parameter plot combining all S4 scans."""
    scans = [
        ('scan_energy_spread.csv', r'$\sigma_E$ (%)', 'Energy Spread'),
        ('scan_chirp.csv', r'$h$ (1/s)', 'Chirp'),
        ('scan_emittance.csv', r'$\varepsilon_n$ ($\pi$·mm·mrad)', 'Emittance'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    for ax, (filename, xlabel, title) in zip(axes, scans):
        data = load_csv(RESULTS / 'params_05ps' / filename)
        if not data:
            ax.set_title(f'{title}: No data')
            continue

        params = [float(r['param_value']) for r in data]
        mses = [float(r['mse']) for r in data]

        ax.semilogy(params, mses, 'o-', color='C0', markersize=5)

        # Threshold lines
        colors = {'Excellent': 'green', 'Acceptable': 'orange', 'Marginal': 'red'}
        for label, thresh in MSE_THRESHOLDS.items():
            ax.axhline(thresh, ls='--', color=colors[label], lw=1, alpha=0.7, label=label)

        ax.set_xlabel(xlabel)
        ax.set_ylabel('RMS')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)

    fig.suptitle('S4: Parameter Sensitivity Summary', fontsize=13)
    plt.tight_layout()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'R2_mse_vs_parameter.{ext}', dpi=150)
    _print(f"  Saved: R2_mse_vs_parameter.{{eps,png}}")
    plt.close(fig)


def plot_cross_code_twiss():
    """Cross-code Twiss parameter bar chart."""
    codes = []
    beta_x_vals = []
    beta_y_vals = []
    alpha_x_vals = []

    # FELsim
    s4 = load_csv(RESULTS / 'params_05ps' / 'scan_energy_spread.csv')
    if s4:
        for r in s4:
            if abs(float(r['param_value']) - 0.5) < 0.01:
                codes.append('FELsim')
                beta_x_vals.append(float(r['beta_x']))
                beta_y_vals.append(float(r['beta_y']))
                alpha_x_vals.append(float(r['alpha_x']))
                break

    # COSY FR0/FR1
    for label, path in [('COSY FR0', 'cosy_s1_fr0.json'),
                         ('COSY FR1', 'cosy_s1_fr1_warm.json')]:
        data = load_json(RESULTS / path)
        if data:
            t = data.get('twiss_undulator', {})
            codes.append(label)
            beta_x_vals.append(t.get('beta_x', 0))
            beta_y_vals.append(t.get('beta_y', 0))
            alpha_x_vals.append(t.get('alpha_x', 0))

    # RFT-opt at ε_n=8
    w8 = load_csv(RESULTS / 'rftrack_opt' / 'comparison.csv')
    if w8:
        for r in w8:
            if abs(float(r['epsilon_n']) - 8) < 0.1 and r['method'] == 'RFT-opt':
                codes.append('RF-Track opt')
                beta_x_vals.append(float(r['beta_x']))
                beta_y_vals.append(float(r['beta_y']))
                alpha_x_vals.append(float(r['alpha_x']))
                break

    if not codes:
        return

    x = np.arange(len(codes))
    width = 0.25

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.bar(x - width/2, beta_x_vals, width, label=r'$\beta_x$', color='C0', alpha=0.8)
    ax1.bar(x + width/2, beta_y_vals, width, label=r'$\beta_y$', color='C1', alpha=0.8)
    ax1.axhline(1.4, color='C0', ls='--', lw=1, alpha=0.5)
    ax1.axhline(0.2418, color='C1', ls='--', lw=1, alpha=0.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(codes, rotation=15, ha='right')
    ax1.set_ylabel('Beta function (m)')
    ax1.set_title('Beta functions at undulator entrance')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')

    ax2.bar(x, alpha_x_vals, width * 2, color='C2', alpha=0.8)
    ax2.axhline(0.47, color='red', ls='--', lw=1, alpha=0.5, label='Target')
    ax2.set_xticks(x)
    ax2.set_xticklabels(codes, rotation=15, ha='right')
    ax2.set_ylabel(r'$\alpha_x$')
    ax2.set_title(r'$\alpha_x$ at undulator entrance')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')

    fig.suptitle('R2: Cross-Code Twiss Comparison (ε_n = 8)', fontsize=13)
    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'R2_cross_code_twiss.{ext}', dpi=150)
    _print(f"  Saved: R2_cross_code_twiss.{{eps,png}}")
    plt.close(fig)


def plot_compression():
    """Compression curve from W12."""
    w12a = load_json(RESULTS / 'W12' / 'part_a_results.json')
    if not w12a or 'sweep' not in w12a:
        return

    sweep = w12a['sweep']
    h_vals = [pt['h'] for pt in sweep]
    sigma_out = [pt['sigma_z_out_ps'] for pt in sweep]
    C_map = [pt['C_map'] for pt in sweep]
    I_peak = [pt['I_peak_out'] for pt in sweep]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5))

    h_display = np.array(h_vals) * 1e-9

    ax1.plot(h_display, sigma_out, 'o-', markersize=4, color='C0')
    ax1.axhline(0.5, color='red', ls='--', lw=1, label='0.5 ps target')
    ax1.axhline(2.0, color='gray', ls='--', lw=1, label='2 ps input')
    ax1.set_xlabel(r'$h$ (10$^9$/s)')
    ax1.set_ylabel(r'$\sigma_z$ out (ps)')
    ax1.set_title('Bunch length after transport')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.plot(h_display, C_map, 'o-', markersize=4, color='C1')
    ax2.axhline(4.0, color='red', ls='--', lw=1, label='C = 4')
    ax2.set_xlabel(r'$h$ (10$^9$/s)')
    ax2.set_ylabel('Compression ratio')
    ax2.set_title('Compression ratio vs chirp')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    ax3.plot(h_display, I_peak, 'o-', markersize=4, color='C2')
    ax3.set_xlabel(r'$h$ (10$^9$/s)')
    ax3.set_ylabel(r'$I_{\mathrm{peak}}$ (A)')
    ax3.set_title('Peak current vs chirp')
    ax3.grid(True, alpha=0.3)

    fig.suptitle('W12: Compression Feasibility', fontsize=13)
    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'R2_compression_curve.{ext}', dpi=150)
    _print(f"  Saved: R2_compression_curve.{{eps,png}}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
#  LaTeX Report
# ═══════════════════════════════════════════════════════════════════════════════

def generate_latex(t1_data, t2_data, t3_data, t4_data, t5_data):
    """Generate LaTeX tables."""
    OUTDIR.mkdir(parents=True, exist_ok=True)

    lines = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[margin=1in]{geometry}",
        r"\usepackage{booktabs,siunitx,longtable}",
        r"\title{R2: Unified Comparison Table — UH MkV FEL Beamline Studies}",
        r"\author{Eremey Valetov}",
        r"\date{\today}",
        r"\begin{document}",
        r"\maketitle",
        "",
    ]

    # Table 1
    if t1_data:
        lines.extend([
            r"\section{Baseline Optimization — Cross-Code (Table~1)}",
            r"\begin{table}[htbp]\centering",
            r"\caption{Cross-code comparison at $\varepsilon_n = 8\;\pi\cdot\text{mm}\cdot\text{mrad}$, "
            r"$\sigma_E = 0.5\%$, $h = 5\times10^9\;/\text{s}$.}",
            r"\begin{tabular}{llS[scientific-notation=true]S[round-precision=4]S[round-precision=4]"
            r"S[round-precision=4]S[round-precision=4]l}",
            r"\toprule",
            r"Code & Model & {RMS} & {$\beta_x$ (m)} & {$\alpha_x$} & {$\beta_y$ (m)} & {$\alpha_y$} & Quality \\",
            r"\midrule",
        ])
        for r in t1_data:
            q = classify_mse(r['MSE'])
            lines.append(
                f"  {r['Code']} & {r['Model']} & {math.sqrt(r['MSE']):.2e} & {r['β_x']:.4f} & "
                f"{r['α_x']:.4f} & {r['β_y']:.4f} & {r['α_y']:.4f} & {q} \\\\"
            )
        lines.extend([
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{table}",
            "",
        ])

    # Table 2
    if t2_data:
        lines.extend([
            r"\section{Parameter Sensitivity Summary (Table~2)}",
            r"\begin{table}[htbp]\centering",
            r"\caption{Quality distribution across S4 parameter scans.}",
            r"\begin{tabular}{lccccccS[scientific-notation=true]S[scientific-notation=true]}",
            r"\toprule",
            r"Scan & Range & $N$ & Excellent & Acceptable & Marginal & Failed & {Best MSE} & {Worst MSE} \\",
            r"\midrule",
        ])
        for r in t2_data:
            lines.append(
                f"  {r['Scan']} & {r['Range']} & {r['N']} & "
                f"{r['Excellent']} & {r['Acceptable']} & {r['Marginal']} & {r['Failed']} & "
                f"{r['Best MSE']:.2e} & {r['Worst MSE']:.2e} \\\\"
            )
        lines.extend([
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{table}",
            "",
        ])

    lines.extend([
        r"\section{Figures}",
        r"\begin{figure}[htbp]\centering",
        r"  \includegraphics[width=\textwidth]{R2/R2_mse_vs_parameter.eps}",
        r"  \caption{RMS vs.\ parameter for three S4 sensitivity scans.}",
        r"\end{figure}",
        "",
        r"\begin{figure}[htbp]\centering",
        r"  \includegraphics[width=0.85\textwidth]{R2/R2_cross_code_twiss.eps}",
        r"  \caption{Cross-code Twiss comparison at $\varepsilon_n = 8$.}",
        r"\end{figure}",
        "",
        r"\begin{figure}[htbp]\centering",
        r"  \includegraphics[width=\textwidth]{R2/R2_compression_curve.eps}",
        r"  \caption{Bunch compression feasibility (W12 Part A).}",
        r"\end{figure}",
        "",
        r"\end{document}",
    ])

    tex_path = Path(__file__).resolve().parent / 'R2_unified_comparison_report.tex'
    with open(tex_path, 'w') as f:
        f.write('\n'.join(lines))
    _print(f"\n  Saved: {tex_path}")


# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='R2: Unified comparison table across all studies')
    parser.add_argument('--plots-only', action='store_true',
                        help='Regenerate plots from existing data only')
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)

    _print("R2: Unified Comparison Table")

    if args.plots_only:
        plot_mse_comparison()
        plot_cross_code_twiss()
        plot_compression()
        return

    t1 = table_1_baseline()
    t2 = table_2_sensitivity()
    t3 = table_3_bunch_transmission()
    t4 = table_4_compression()
    t5 = table_5_currents()

    # Plots
    plot_mse_comparison()
    plot_cross_code_twiss()
    plot_compression()

    # LaTeX
    generate_latex(t1, t2, t3, t4, t5)

    # Save aggregated JSON
    summary = {
        'table_1_baseline': t1,
        'table_2_sensitivity': t2,
        'table_4_compression_R56_m': t4.get('w12a', {}).get('R56_m') if isinstance(t4, dict) else None,
    }
    with open(OUTDIR / 'R2_summary.json', 'w') as f:
        json.dump(summary, f, indent=2, default=_sanitize)
    _print(f"  Saved: {OUTDIR / 'R2_summary.json'}")

    _print("\n" + "=" * 72)
    _print("  R2 Complete")
    _print("=" * 72)


if __name__ == "__main__":
    main()
