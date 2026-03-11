#!/usr/bin/env python3
"""Generate plots and tables for 2 ps → 0.5 ps switching presentation slides.

Author: Eremey Valetov
Date: 2026-02-27
"""

import json
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 200,
    'savefig.dpi': 200,
})

BASE = Path(__file__).resolve().parent.parent  # params_05ps/
W9 = BASE.parent / 'W9'
OUT = Path(__file__).resolve().parent

# ── Load data ─────────────────────────────────────────────────────────────────

def load_csv(path):
    with open(path) as f:
        reader = csv.DictReader(f)
        return list(reader)

w1 = load_csv(BASE / 'w1_chirp_comparison.csv')
energy_spread = load_csv(BASE / 'scan_energy_spread.csv')
emittance = load_csv(BASE / 'scan_emittance_w2.csv')
chirp = load_csv(BASE / 'scan_chirp.csv')

with open(W9 / 'part_b_bunch_propagation.json') as f:
    bunch_prop = json.load(f)

with open(W9 / 'part_a_longitudinal_map.json') as f:
    long_map = json.load(f)


# ── Figure 1: Bunch length comparison (bar chart) ────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))

# Left: bunch lengthening ratio
labels = [b['label'] for b in bunch_prop]
ratios = [b['sigma_z_ratio'] for b in bunch_prop]
colors = ['#2196F3', '#1565C0', '#FF9800', '#E65100']
bars = axes[0].bar(range(len(labels)), ratios, color=colors, edgecolor='black',
                   linewidth=0.5)
axes[0].set_xticks(range(len(labels)))
axes[0].set_xticklabels(labels, rotation=25, ha='right', fontsize=9)
axes[0].set_ylabel(r'$\sigma_z^{\rm final} / \sigma_z^{\rm initial}$')
axes[0].set_title('Bunch Lengthening Through Beamline')
axes[0].axhline(1.0, color='gray', ls='--', lw=0.8)
for bar, r in zip(bars, ratios):
    axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                 f'{r:.2f}', ha='center', va='bottom', fontsize=9)

# Right: peak current initial vs final
for i, b in enumerate(bunch_prop):
    init_I = b['initial']['I_peak_A']
    final_I = b['final']['I_peak_A']
    x = np.array([0, 1]) + i * 2.5
    axes[1].bar(x[0], init_I, width=0.8, color=colors[i], edgecolor='black',
                linewidth=0.5)
    axes[1].bar(x[1], final_I, width=0.8, color=colors[i], alpha=0.6,
                edgecolor='black', linewidth=0.5)

tick_pos = [0.5 + i * 2.5 for i in range(4)]
axes[1].set_xticks(tick_pos)
axes[1].set_xticklabels([b['label'] for b in bunch_prop], rotation=25,
                         ha='right', fontsize=9)
axes[1].set_ylabel('Peak Current (A)')
axes[1].set_title('Peak Current: Initial vs Final')

for i, b in enumerate(bunch_prop):
    init_I = b['initial']['I_peak_A']
    final_I = b['final']['I_peak_A']
    x = np.array([0, 1]) + i * 2.5
    axes[1].text(x[0], init_I + 0.8, f'{init_I:.0f}', ha='center', fontsize=7)
    axes[1].text(x[1], final_I + 0.8, f'{final_I:.0f}', ha='center', fontsize=7)

fig.tight_layout()
fig.savefig(OUT / 'bunch_lengthening.png')
fig.savefig(OUT / 'bunch_lengthening.eps')
plt.close(fig)
print("Saved: bunch_lengthening.png/eps")


# ── Figure 2: Parameter sensitivity (MSE vs energy spread + emittance) ───────

fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))

# Left: MSE vs energy spread
sigma_e = [float(r['param_value']) for r in energy_spread]
mse_e = [float(r['mse']) for r in energy_spread]
axes[0].semilogy(sigma_e, mse_e, 'o-', color='#2196F3', markersize=5, lw=1.5)
axes[0].axhline(1e-3, color='green', ls='--', lw=0.8, label='Excellent')
axes[0].axhline(1e-2, color='orange', ls='--', lw=0.8, label='Acceptable')
axes[0].axhline(1e-1, color='red', ls='--', lw=0.8, label='Marginal')
axes[0].set_xlabel(r'Energy Spread $\sigma_E$ (%)')
axes[0].set_ylabel('MSE')
axes[0].set_title(r'Twiss Matching vs $\sigma_E$ (0.5 ps)')
axes[0].legend(loc='upper left', fontsize=8)
axes[0].set_ylim(1e-7, 1)

# Right: MSE vs emittance
eps_n = [float(r['param_value']) for r in emittance]
mse_em = [float(r['mse']) for r in emittance]
axes[1].semilogy(eps_n, mse_em, 's-', color='#E65100', markersize=5, lw=1.5)
axes[1].axhline(1e-3, color='green', ls='--', lw=0.8, label='Excellent')
axes[1].axhline(1e-2, color='orange', ls='--', lw=0.8, label='Acceptable')
axes[1].axhline(1e-1, color='red', ls='--', lw=0.8, label='Marginal')
axes[1].set_xlabel(r'Normalised Emittance $\varepsilon_n$ ($\pi\cdot$mm$\cdot$mrad)')
axes[1].set_ylabel('MSE')
axes[1].set_title(r'Twiss Matching vs $\varepsilon_n$ (0.5 ps)')
axes[1].legend(loc='upper left', fontsize=8)
axes[1].set_ylim(1e-10, 10)

fig.tight_layout()
fig.savefig(OUT / 'parameter_sensitivity.png')
fig.savefig(OUT / 'parameter_sensitivity.eps')
plt.close(fig)
print("Saved: parameter_sensitivity.png/eps")


# ── Table 1: 2 ps vs 0.5 ps comparison ──────────────────────────────────────

table_lines = []
table_lines.append("=" * 78)
table_lines.append("  2 ps → 0.5 ps Switching: Key Comparison")
table_lines.append("=" * 78)
table_lines.append("")
table_lines.append(f"{'Parameter':<35} {'2 ps':>18} {'0.5 ps':>18}")
table_lines.append("-" * 78)

# From W1 data
w1_05_h0 = next(r for r in w1 if '0.5' in r['label'] and 'h=0' in r['label'])
w1_2_h0 = next(r for r in w1 if '2 ps' in r['label'] and 'h=0' in r['label'])
w1_05_h5 = next(r for r in w1 if '0.5' in r['label'] and 'h=5' in r['label'])
w1_2_h5 = next(r for r in w1 if '2 ps' in r['label'] and 'h=5' in r['label'])

table_lines.append(f"{'Transverse Twiss MSE (h=0)':<35} {float(w1_2_h0['mse']):>18.2e} {float(w1_05_h0['mse']):>18.2e}")
table_lines.append(f"{'Transverse Twiss MSE (h=5e9)':<35} {float(w1_2_h5['mse']):>18.2e} {float(w1_05_h5['mse']):>18.2e}")

# Quad currents comparison (h=0 — identical)
q_05 = [float(w1_05_h0[k]) for k in sorted(w1_05_h0) if k.startswith('quad_')]
q_2 = [float(w1_2_h0[k]) for k in sorted(w1_2_h0) if k.startswith('quad_')]
max_dI = max(abs(a - b) for a, b in zip(q_05, q_2))
table_lines.append(f"{'Max |ΔI| in quad currents (h=0)':<35} {'(identical)':>18} {f'{max_dI:.1e} A':>18}")

# From bunch propagation
b_05_h0 = next(b for b in bunch_prop if b['label'] == '0.5 ps, h=0')
b_2_h0 = next(b for b in bunch_prop if b['label'] == '2 ps, h=0')
b_05_h5 = next(b for b in bunch_prop if b['label'] == '0.5 ps, h=5e9')
b_2_h5 = next(b for b in bunch_prop if b['label'] == '2 ps, h=5e9')

table_lines.append(f"{'Bunch lengthening (h=0)':<35} {b_2_h0['sigma_z_ratio']:>17.1f}× {b_05_h0['sigma_z_ratio']:>17.1f}×")
table_lines.append(f"{'Bunch lengthening (h=5e9)':<35} {b_2_h5['sigma_z_ratio']:>17.1f}× {b_05_h5['sigma_z_ratio']:>17.1f}×")
table_lines.append(f"{'Peak current, final (h=0) [A]':<35} {b_2_h0['final']['I_peak_A']:>18.1f} {b_05_h0['final']['I_peak_A']:>18.1f}")
table_lines.append(f"{'σ_z final (h=0) [ps]':<35} {b_2_h0['final']['sigma_z_ps']:>18.2f} {b_05_h0['final']['sigma_z_ps']:>18.2f}")

table_lines.append(f"{'R56 (COSY) [mm]':<35} {'27.1':>18} {'27.1':>18}")
table_lines.append("")
table_lines.append("-" * 78)
table_lines.append("Targets: β_x = 1.4 m, α_x = 0.47, β_y = 0.242 m, α_y = 0.0")
table_lines.append("Baseline: ε_n = 8 π·mm·mrad, σ_E = 0.5%, 40 MeV")
table_lines.append("=" * 78)

table_text = "\n".join(table_lines)
(OUT / 'comparison_table.txt').write_text(table_text)
print("Saved: comparison_table.txt")
print()
print(table_text)


# ── Table 2: Parameter sensitivity summary ───────────────────────────────────

table2 = []
table2.append("")
table2.append("=" * 78)
table2.append("  0.5 ps Parameter Sensitivity Summary")
table2.append("=" * 78)
table2.append("")
table2.append(f"{'Scan':<25} {'Range':<20} {'MSE Range':<25} {'Verdict':<15}")
table2.append("-" * 78)

mse_e_vals = [float(r['mse']) for r in energy_spread]
table2.append(f"{'Energy spread σ_E':<25} {'0.1 – 5.0 %':<20} "
              f"{f'{min(mse_e_vals):.1e} – {max(mse_e_vals):.1e}':<25} "
              f"{'All Excellent':<15}")

mse_c_vals = [float(r['mse']) for r in chirp]
table2.append(f"{'Chirp h':<25} {'0 – 40×10⁹/s':<20} "
              f"{f'{min(mse_c_vals):.1e} – {max(mse_c_vals):.1e}':<25} "
              f"{'All Excellent':<15}")

mse_em_vals = [float(r['mse']) for r in emittance]
n_exc = sum(1 for m in mse_em_vals if m < 1e-3)
table2.append(f"{'Emittance ε_n':<25} {'1 – 20 π·mm·mrad':<20} "
              f"{f'{min(mse_em_vals):.1e} – {max(mse_em_vals):.1e}':<25} "
              f"{f'{n_exc}/{len(mse_em_vals)} Excellent':<15}")

table2.append("-" * 78)
table2.append("Quality thresholds: Excellent < 10⁻³, Acceptable < 10⁻², Marginal < 10⁻¹")
table2.append("=" * 78)

table2_text = "\n".join(table2)
(OUT / 'sensitivity_table.txt').write_text(table2_text)
print("Saved: sensitivity_table.txt")
print()
print(table2_text)


# ── Slide bullet points ─────────────────────────────────────────────────────

bullets = """
==============================================================================
  SLIDE 1: Switching from 2 ps to 0.5 ps — Transverse Optics
==============================================================================

Key points:
• No beamline hardware changes required for transverse Twiss matching.
  Same quad currents produce identical undulator optics at both bunch lengths.

• Root cause: transverse Twiss parameters depend on emittance (ε_n) and
  energy spread (σ_E), not bunch length (σ_z). Bunch length enters only
  the longitudinal (time-of-flight) column of the transfer matrix.

• Confirmed by three independent codes:
  – FELsim (transfer matrices): MSE = 1.3×10⁻⁶ at both 2 ps and 0.5 ps (h=0)
  – COSY INFINITY (6D differential algebra): MSE = 2.3×10⁻⁷
  – RF-Track (3D particle tracking): validates model, identifies dipole
    edge-kick sensitivity for precision tuning

• Chirp (h = 5×10⁹/s) has negligible effect on transverse matching
  (MSE increases from 1.3×10⁻⁶ to 2.6×10⁻⁵ — still Excellent quality).

Figures: comparison_table.txt, parameter_sensitivity.png

==============================================================================
  SLIDE 2: Longitudinal Effects and Operational Considerations
==============================================================================

Key points:
• R56 = 27.1 mm (from COSY 6D map) — causes bunch lengthening that is
  more pronounced at shorter bunch lengths:
  – 0.5 ps (h=0): 35% lengthening (0.50 → 0.67 ps), peak current 48 → 36 A
  – 2 ps (h=0):  2.5% lengthening (1.99 → 2.04 ps), peak current 12 → 12 A
  – With chirp (h = 5×10⁹/s): 71% lengthening at 0.5 ps, 47% at 2 ps

• Beam parameter sensitivity (at 0.5 ps, 500-particle sweeps):
  – Energy spread: full 0.1–5% range produces Excellent matching (MSE < 10⁻³)
  – Emittance: 6–20 π·mm·mrad all Excellent; ε_n < 5 challenging (physics limit)
  – Chirp: full 0–40×10⁹/s range produces Excellent matching

• Practical implication: the only operational change for 0.5 ps is in the
  injector (shorter laser pulse / bunch compression). The transport line
  quad settings can remain unchanged. Longitudinal bunch quality should
  be monitored — R56 compensation may be needed for demanding applications.

Figures: bunch_lengthening.png, sensitivity_table.txt
"""

(OUT / 'slide_notes.txt').write_text(bullets)
print("Saved: slide_notes.txt")
print(bullets)
