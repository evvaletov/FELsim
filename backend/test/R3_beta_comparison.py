#!/usr/bin/env python3
"""Three-way beta function comparison: FELsim, COSY, RF-Track.

Applies the COSY FR3+MGE optimised currents to all three simulation backends
and compares beta_x(s), beta_y(s) along the beamline.

Author: Eremey Valetov
"""

import sys
import json
import math
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

from ebeam import beam
from beamline import lattice
from excelElements import ExcelElements

EXCEL_PATH = Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx'

# Beam parameters
ENERGY = 40       # MeV
EPSILON_N = 8     # pi.mm.mrad (normalized)
X_STD = 0.8       # mm
Y_STD = 0.8       # mm
FREQ = 2856e6     # Hz
BUNCH_SPREAD = 2  # ps
ENERGY_STD_PCT = 0.5
H_CHIRP = 5e9
NB_PARTICLES = 2000

# Undulator targets
K = 1.2
LAMBDA_U = 2.3e-2  # m

QUAD_INDICES = [1, 3, 10, 16, 18, 20, 27, 33, 35, 37, 39, 41, 43,
                50, 56, 58, 61, 63, 70, 76, 78, 80, 87, 93, 95, 97]

# Stage 5 mirror symmetry
MIRROR = {39: 37, 41: 35, 43: 33}


def compute_params():
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    norm = relat.gamma * relat.beta
    epsilon = EPSILON_N / norm
    beta_ym = relat.gamma * LAMBDA_U / (2 * np.pi * K)
    beta_0 = X_STD**2 / epsilon
    return {
        'beta_xm': 1.4, 'alpha_xm': 0.47,
        'beta_ym': beta_ym, 'alpha_ym': 0.0,
        'epsilon': epsilon, 'beta_0': beta_0,
        'gamma': relat.gamma, 'beta_rel': relat.beta, 'norm': norm,
    }


def generate_beam(params):
    np.random.seed(42)
    epsilon = params['epsilon']
    x_prime_std = epsilon / X_STD
    y_prime_std = epsilon / Y_STD
    tof_std = BUNCH_SPREAD * 1e-9 * FREQ
    energy_std = ENERGY_STD_PCT * 10

    ebeam_gen = beam()
    beam_dist = ebeam_gen.gen_6d_gaussian(
        0, [X_STD, x_prime_std, Y_STD, y_prime_std, tof_std, energy_std],
        NB_PARTICLES
    )
    tof_dist = beam_dist[:, 4] / FREQ
    beam_dist[:, 5] += H_CHIRP * tof_dist
    return beam_dist


def load_currents(json_path):
    with open(json_path) as f:
        data = json.load(f)
    currents = {int(k): v for k, v in data['currents'].items()}
    # Apply mirror symmetry for stage 5
    for dst, src in MIRROR.items():
        if src in currents and dst not in currents:
            currents[dst] = currents[src]
    return currents, data


def build_felsim_line(currents):
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    excel = ExcelElements(str(EXCEL_PATH))
    beamlineUH = excel.create_beamline()
    line = relat.changeBeamType("electron", ENERGY, beamlineUH)[:118]
    for idx, current in currents.items():
        if idx < len(line):
            line[idx].current = abs(current)
    return line


def propagate_felsim_particles(line, beam_dist):
    ebeam_calc = beam()
    particles = beam_dist.copy()
    s = 0.0
    results = []
    for idx, elem in enumerate(line):
        particles = np.array(elem.useMatrice(particles))
        s += getattr(elem, 'length', getattr(elem, 'L', 0))
        _, _, twiss_df = ebeam_calc.cal_twiss(particles, ddof=1)
        results.append({
            'index': idx, 's': s,
            'beta_x': twiss_df.loc['x', r'$\beta$ (m)'],
            'beta_y': twiss_df.loc['y', r'$\beta$ (m)'],
            'alpha_x': twiss_df.loc['x', r'$\alpha$'],
            'alpha_y': twiss_df.loc['y', r'$\alpha$'],
        })
    return results


def propagate_felsim_matrix(line, beta_0):
    bx, ax = beta_0, 0.0
    gx = (1 + ax**2) / bx
    by, ay = beta_0, 0.0
    gy = (1 + ay**2) / by
    s = 0.0
    results = []
    for idx, elem in enumerate(line):
        M = elem._compute_numeric_matrix()
        s += getattr(elem, 'length', getattr(elem, 'L', 0))
        m11, m12, m21, m22 = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
        bx, ax, gx = (m11**2*bx - 2*m11*m12*ax + m12**2*gx,
                       -m11*m21*bx + (m11*m22 + m12*m21)*ax - m12*m22*gx,
                       m21**2*bx - 2*m21*m22*ax + m22**2*gx)
        m11, m12, m21, m22 = M[2, 2], M[2, 3], M[3, 2], M[3, 3]
        by, ay, gy = (m11**2*by - 2*m11*m12*ay + m12**2*gy,
                       -m11*m21*by + (m11*m22 + m12*m21)*ay - m12*m22*gy,
                       m21**2*by - 2*m21*m22*ay + m22**2*gy)
        results.append({'index': idx, 's': s, 'beta_x': bx, 'beta_y': by,
                         'alpha_x': ax, 'alpha_y': ay})
    return results


def run_rftrack_evolution(currents, beam_dist):
    try:
        from rftrackAdapter import RFTrackAdapter
    except ImportError:
        print("  RF-Track adapter not available")
        return None

    try:
        sim = RFTrackAdapter(
            lattice_path=str(EXCEL_PATH),
            beam_energy=ENERGY,
            space_charge=False,
            aperture=0.5,
        )
    except ImportError:
        print("  RF-Track package not installed")
        return None

    sim.beamline = sim.beamline[:118]
    sim._build_lattice()

    for idx, current in currents.items():
        if idx < len(sim.beamline):
            sim._modify_element(idx, current=abs(current))
    sim._build_lattice()

    print("  Tracking particles through RF-Track...")
    evolution = sim.collect_evolution(beam_dist, checkpoint_elements='all')

    results = []
    for s_pos in evolution.s_positions:
        twiss = evolution.twiss.get(s_pos, {})
        tx = twiss.get('x', {})
        ty = twiss.get('y', {})
        if tx and ty:
            results.append({
                's': s_pos,
                'beta_x': tx.get('beta', np.nan),
                'beta_y': ty.get('beta', np.nan),
                'alpha_x': tx.get('alpha', np.nan),
                'alpha_y': ty.get('alpha', np.nan),
            })
    return results


def plot_comparison(felsim_part, felsim_mat, rftrack, targets, cosy_data,
                    output_path):
    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)

    s_fm = [r['s'] for r in felsim_mat]

    # Beta functions
    ax = axes[0]
    ax.plot(s_fm, [r['beta_x'] for r in felsim_mat],
            'b-', lw=1.5, label=r'FELsim $\beta_x$ (matrix)')
    ax.plot(s_fm, [r['beta_y'] for r in felsim_mat],
            'r-', lw=1.5, label=r'FELsim $\beta_y$ (matrix)')

    if felsim_part:
        s_fp = [r['s'] for r in felsim_part]
        ax.plot(s_fp, [r['beta_x'] for r in felsim_part],
                'b--', lw=1, alpha=0.5, label=r'FELsim $\beta_x$ (particles)')
        ax.plot(s_fp, [r['beta_y'] for r in felsim_part],
                'r--', lw=1, alpha=0.5, label=r'FELsim $\beta_y$ (particles)')

    if rftrack:
        s_rf = [r['s'] for r in rftrack]
        ax.plot(s_rf, [r['beta_x'] for r in rftrack],
                'b:', lw=1.5, alpha=0.8, label=r'RF-Track $\beta_x$')
        ax.plot(s_rf, [r['beta_y'] for r in rftrack],
                'r:', lw=1.5, alpha=0.8, label=r'RF-Track $\beta_y$')

    # COSY endpoint marker (from JSON metadata)
    if cosy_data and 'twiss_undulator' in cosy_data:
        tw = cosy_data['twiss_undulator']
        s_end = s_fm[-1]
        ax.plot(s_end, tw['beta_x'], 'bv', ms=10, zorder=5,
                label=f"COSY $\\beta_x$ = {tw['beta_x']:.2f} m")
        ax.plot(s_end, tw['beta_y'], 'r^', ms=10, zorder=5,
                label=f"COSY $\\beta_y$ = {tw['beta_y']:.4f} m")

    ax.axhline(targets['beta_xm'], color='b', ls=':', alpha=0.3,
               label=f"Target $\\beta_x$ = {targets['beta_xm']} m")
    ax.axhline(targets['beta_ym'], color='r', ls=':', alpha=0.3,
               label=f"Target $\\beta_y$ = {targets['beta_ym']:.4f} m")
    ax.set_ylabel(r'$\beta$ (m)', fontsize=12)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=7, ncol=2, loc='upper right')
    ax.grid(True, alpha=0.3)

    fr = cosy_data.get('fringe_field_order', '?') if cosy_data else '?'
    mge = cosy_data.get('mge', False) if cosy_data else False
    mse = cosy_data.get('mse', '?') if cosy_data else '?'
    rms_val = math.sqrt(mse) if isinstance(mse, (int, float)) else mse
    ax.set_title(
        f"Beta functions — COSY FR{fr}{'+MGE' if mge else ''} optimised currents "
        f"(RMS = {rms_val:.1f})\n"
        f"FELsim (matrix + particles) vs RF-Track vs COSY endpoint",
        fontsize=11
    )

    # Alpha functions
    ax = axes[1]
    ax.plot(s_fm, [r['alpha_x'] for r in felsim_mat],
            'b-', lw=1.5, label=r'FELsim $\alpha_x$ (matrix)')
    ax.plot(s_fm, [r['alpha_y'] for r in felsim_mat],
            'r-', lw=1.5, label=r'FELsim $\alpha_y$ (matrix)')

    if felsim_part:
        s_fp = [r['s'] for r in felsim_part]
        ax.plot(s_fp, [r['alpha_x'] for r in felsim_part],
                'b--', lw=1, alpha=0.5, label=r'FELsim $\alpha_x$ (particles)')
        ax.plot(s_fp, [r['alpha_y'] for r in felsim_part],
                'r--', lw=1, alpha=0.5, label=r'FELsim $\alpha_y$ (particles)')

    if rftrack:
        s_rf = [r['s'] for r in rftrack]
        ax.plot(s_rf, [r['alpha_x'] for r in rftrack],
                'b:', lw=1.5, alpha=0.8, label=r'RF-Track $\alpha_x$')
        ax.plot(s_rf, [r['alpha_y'] for r in rftrack],
                'r:', lw=1.5, alpha=0.8, label=r'RF-Track $\alpha_y$')

    if cosy_data and 'twiss_undulator' in cosy_data:
        tw = cosy_data['twiss_undulator']
        s_end = s_fm[-1]
        ax.plot(s_end, tw['alpha_x'], 'bv', ms=10, zorder=5,
                label=f"COSY $\\alpha_x$ = {tw['alpha_x']:.2f}")
        ax.plot(s_end, tw['alpha_y'], 'r^', ms=10, zorder=5,
                label=f"COSY $\\alpha_y$ = {tw['alpha_y']:.2f}")

    ax.axhline(targets['alpha_xm'], color='b', ls=':', alpha=0.3)
    ax.axhline(targets['alpha_ym'], color='r', ls=':', alpha=0.3)
    ax.set_ylabel(r'$\alpha$', fontsize=12)
    ax.set_xlabel('s (m)', fontsize=12)
    ax.legend(fontsize=7, ncol=2, loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: {output_path}")

    # Clipped version with y-axis cap
    base = Path(output_path)
    clipped_path = base.parent / (base.stem + '_clipped' + base.suffix)

    fig2, axes2 = plt.subplots(2, 1, figsize=(14, 9), sharex=True)

    ax = axes2[0]
    ax.plot(s_fm, [r['beta_x'] for r in felsim_mat],
            'b-', lw=1.5, label=r'FELsim $\beta_x$ (matrix)')
    ax.plot(s_fm, [r['beta_y'] for r in felsim_mat],
            'r-', lw=1.5, label=r'FELsim $\beta_y$ (matrix)')
    if felsim_part:
        s_fp = [r['s'] for r in felsim_part]
        ax.plot(s_fp, [r['beta_x'] for r in felsim_part],
                'b--', lw=1, alpha=0.5, label=r'FELsim $\beta_x$ (particles)')
        ax.plot(s_fp, [r['beta_y'] for r in felsim_part],
                'r--', lw=1, alpha=0.5, label=r'FELsim $\beta_y$ (particles)')
    if rftrack:
        s_rf = [r['s'] for r in rftrack]
        ax.plot(s_rf, [r['beta_x'] for r in rftrack],
                'b:', lw=1.5, alpha=0.8, label=r'RF-Track $\beta_x$')
        ax.plot(s_rf, [r['beta_y'] for r in rftrack],
                'r:', lw=1.5, alpha=0.8, label=r'RF-Track $\beta_y$')
    ax.axhline(targets['beta_xm'], color='b', ls=':', alpha=0.3)
    ax.axhline(targets['beta_ym'], color='r', ls=':', alpha=0.3)
    ax.set_ylabel(r'$\beta$ (m)', fontsize=12)
    ax.set_yscale('log')
    ax.set_ylim(1e-2, 1e4)
    ax.legend(fontsize=7, ncol=2, loc='upper left')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title(
        f"Beta functions (log scale) — COSY FR{fr}{'+MGE' if mge else ''} "
        f"optimised currents (RMS = {rms_val:.1f})\n"
        f"FELsim vs RF-Track — NOTE: solution is unstable ($\\beta_x$ diverges)",
        fontsize=11
    )

    ax = axes2[1]
    ax.plot(s_fm, [r['beta_x'] for r in felsim_mat],
            'b-', lw=1.5, label=r'FELsim $\beta_x$ (matrix)')
    ax.plot(s_fm, [r['beta_y'] for r in felsim_mat],
            'r-', lw=1.5, label=r'FELsim $\beta_y$ (matrix)')
    if felsim_part:
        ax.plot(s_fp, [r['beta_x'] for r in felsim_part],
                'b--', lw=1, alpha=0.5, label=r'FELsim $\beta_x$ (particles)')
        ax.plot(s_fp, [r['beta_y'] for r in felsim_part],
                'r--', lw=1, alpha=0.5, label=r'FELsim $\beta_y$ (particles)')
    if rftrack:
        ax.plot(s_rf, [r['beta_x'] for r in rftrack],
                'b:', lw=1.5, alpha=0.8, label=r'RF-Track $\beta_x$')
        ax.plot(s_rf, [r['beta_y'] for r in rftrack],
                'r:', lw=1.5, alpha=0.8, label=r'RF-Track $\beta_y$')
    ax.axhline(targets['beta_xm'], color='b', ls=':', alpha=0.3,
               label=f"Target $\\beta_x$ = {targets['beta_xm']} m")
    ax.axhline(targets['beta_ym'], color='r', ls=':', alpha=0.3,
               label=f"Target $\\beta_y$ = {targets['beta_ym']:.4f} m")
    ax.set_ylabel(r'$\beta$ (m)', fontsize=12)
    ax.set_ylim(0, 100)
    ax.set_xlabel('s (m)', fontsize=12)
    ax.legend(fontsize=7, ncol=2, loc='upper left')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(clipped_path, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {clipped_path}")


def print_endpoint_comparison(felsim_mat, felsim_part, rftrack, cosy_data, targets):
    print(f"\n{'='*80}")
    print("  Twiss at undulator entrance (s = end)")
    print(f"{'='*80}")
    header = f"  {'Parameter':<12} {'Target':>8} {'FELsim(M)':>10} {'FELsim(P)':>10}"
    if rftrack:
        header += f" {'RF-Track':>10}"
    if cosy_data and 'twiss_undulator' in cosy_data:
        header += f" {'COSY':>10}"
    print(header)
    print(f"  {'-'*70}")

    fm = felsim_mat[-1]
    fp = felsim_part[-1] if felsim_part else {}
    rf = rftrack[-1] if rftrack else {}
    tw = cosy_data.get('twiss_undulator', {}) if cosy_data else {}

    for name, key, target in [
        ('beta_x', 'beta_x', targets['beta_xm']),
        ('beta_y', 'beta_y', targets['beta_ym']),
        ('alpha_x', 'alpha_x', targets['alpha_xm']),
        ('alpha_y', 'alpha_y', targets['alpha_ym']),
    ]:
        row = f"  {name:<12} {target:8.4f} {fm.get(key, 0):10.4f} {fp.get(key, 0):10.4f}"
        if rftrack:
            row += f" {rf.get(key, 0):10.4f}"
        if tw:
            row += f" {tw.get(key, 0):10.4f}"
        print(row)
    print(f"{'='*80}")


def run_felsim_optimised(params, beam_dist):
    """Run fresh FELsim 11-stage optimization and propagate."""
    from beamOptimizer import beamOptimizer

    line = build_felsim_line({})  # unoptimised
    opti = beamOptimizer(line, beam_dist)
    axm, aym = params['alpha_xm'], params['alpha_ym']
    bxm, bym = params['beta_xm'], params['beta_ym']

    opti.calc("Nelder-Mead",
              {1: ["I", "current", lambda n: n], 3: ["I2", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}, "I2": {"bounds": (0, 10), "start": 1}},
              {8: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                   {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.0}],
               9: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                   {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]})
    opti.calc("Nelder-Mead",
              {10: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {15: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]})
    opti.calc("Nelder-Mead",
              {16: ["I", "current", lambda n: n], 18: ["I2", "current", lambda n: n],
               20: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 5},
               "I3": {"bounds": (0, 10), "start": 3}},
              {25: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
               26: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]})
    opti.calc("Nelder-Mead",
              {27: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {32: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]})
    opti.calc("Nelder-Mead",
              {37: ["I", "current", lambda n: n], 35: ["I2", "current", lambda n: n],
               33: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2},
               "I3": {"bounds": (0, 10), "start": 2}},
              {37: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["x", "envelope"], "goal": 2.0, "weight": 1},
                    {"measure": ["y", "envelope"], "goal": 2.0, "weight": 1}]})
    line[43].current = line[33].current
    line[41].current = line[35].current
    line[39].current = line[37].current
    opti.calc("Nelder-Mead",
              {50: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {55: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]})
    opti.calc("Nelder-Mead",
              {56: ["I", "current", lambda n: n], 58: ["I2", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2}},
              {59: [{"measure": ["x", "envelope"], "goal": 0.0, "weight": 1},
                    {"measure": ["y", "envelope"], "goal": 0.0, "weight": 1}]})
    opti.calc("Nelder-Mead",
              {61: ["I", "current", lambda n: n], 63: ["I2", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2}},
              {68: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
               69: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]})
    opti.calc("Nelder-Mead",
              {70: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {75: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]})
    opti.calc("Nelder-Mead",
              {76: ["I", "current", lambda n: n], 78: ["I2", "current", lambda n: n],
               80: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2},
               "I3": {"bounds": (0, 10), "start": 2}},
              {85: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
               86: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]})
    opti.calc("Nelder-Mead",
              {87: ["Ic", "current", lambda n: n], 93: ["I", "current", lambda n: n],
               95: ["I2", "current", lambda n: n], 97: ["I3", "current", lambda n: n]},
              {"Ic": {"bounds": (0, 10), "start": 4}, "I": {"bounds": (0, 10), "start": 2},
               "I2": {"bounds": (0, 10), "start": 2}, "I3": {"bounds": (0, 10), "start": 2}},
              {92: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 0.5}],
               117: [{"measure": ["x", "alpha"], "goal": axm, "weight": 1},
                     {"measure": ["y", "alpha"], "goal": aym, "weight": 1},
                     {"measure": ["x", "beta"], "goal": bxm, "weight": 1},
                     {"measure": ["y", "beta"], "goal": bym, "weight": 1}]})

    currents = {idx: line[idx].current for idx in QUAD_INDICES}
    return line, currents


def plot_three_panel(results_list, targets, output_path):
    """Three-panel plot: each code's beta(s) with its own best currents."""
    n = len(results_list)
    fig, axes = plt.subplots(n, 1, figsize=(14, 4*n), sharex=True)
    if n == 1:
        axes = [axes]

    colors = {'FELsim': ('tab:blue', 'tab:red'),
              'COSY': ('tab:blue', 'tab:red'),
              'RF-Track': ('tab:blue', 'tab:red')}

    for ax, res in zip(axes, results_list):
        data = res['data']
        s = [r['s'] for r in data]

        ax.plot(s, [r['beta_x'] for r in data],
                'b-', lw=1.5, label=r'$\beta_x$')
        ax.plot(s, [r['beta_y'] for r in data],
                'r-', lw=1.5, label=r'$\beta_y$')

        # COSY endpoint marker if available
        if res.get('cosy_endpoint'):
            tw = res['cosy_endpoint']
            ax.plot(s[-1], tw['beta_x'], 'bv', ms=10, zorder=5,
                    label=f"COSY endpoint $\\beta_x$ = {tw['beta_x']:.4f}")
            ax.plot(s[-1], tw['beta_y'], 'r^', ms=10, zorder=5,
                    label=f"COSY endpoint $\\beta_y$ = {tw['beta_y']:.4f}")

        ax.axhline(targets['beta_xm'], color='b', ls=':', alpha=0.3)
        ax.axhline(targets['beta_ym'], color='r', ls=':', alpha=0.3)
        ax.set_ylabel(r'$\beta$ (m)', fontsize=12)

        # Auto y-limit: reasonable range
        betas = [r['beta_x'] for r in data] + [r['beta_y'] for r in data]
        ymax = min(max(betas) * 1.3, 50)
        ax.set_ylim(0, max(ymax, 5))
        ax.legend(fontsize=8, ncol=2, loc='upper right')
        ax.grid(True, alpha=0.3)

        final = data[-1]
        rms_str = f"{math.sqrt(res.get('mse', 0)):.2e}"
        ax.set_title(
            f"{res['label']}  —  "
            f"$\\beta_x$={final['beta_x']:.4f} m, "
            f"$\\beta_y$={final['beta_y']:.4f} m, "
            f"$\\alpha_x$={final['alpha_x']:.4f}, "
            f"$\\alpha_y$={final['alpha_y']:.4f}  "
            f"(RMS = {rms_str})",
            fontsize=10
        )

    axes[-1].set_xlabel('s (m)', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Three-way beta function comparison: each code's own optimised currents")
    parser.add_argument('--skip-rftrack', action='store_true')
    parser.add_argument('--output', type=str, default=None)
    args = parser.parse_args()

    results_dir = Path(__file__).resolve().parent / 'results'
    params = compute_params()

    print("R3: Three-way beta function comparison")
    print(f"  beta_0 = {params['beta_0']:.4f} m")
    print(f"  Targets: bx={params['beta_xm']}, ax={params['alpha_xm']}, "
          f"by={params['beta_ym']:.4f}, ay={params['alpha_ym']}")

    beam_dist = generate_beam(params)
    all_results = []

    # ── 1. FELsim: its own optimised currents, its own propagation ────────
    print("\n[1/3] FELsim 11-stage Nelder-Mead optimization + propagation...")
    line, felsim_currents = run_felsim_optimised(params, beam_dist)
    felsim_part = propagate_felsim_particles(line, beam_dist)
    fp = felsim_part[-1]
    felsim_mse = (
        (fp['beta_x'] - params['beta_xm'])**2
        + (fp['beta_y'] - params['beta_ym'])**2
        + (fp['alpha_x'] - params['alpha_xm'])**2
        + (fp['alpha_y'] - params['alpha_ym'])**2
    ) / 4
    print(f"  bx={fp['beta_x']:.4f}, by={fp['beta_y']:.4f}, "
          f"ax={fp['alpha_x']:.4f}, ay={fp['alpha_y']:.4f}, RMS={math.sqrt(felsim_mse):.2e}")
    all_results.append({
        'label': 'FELsim (particle tracking, Nelder-Mead)',
        'data': felsim_part, 'mse': felsim_mse,
    })

    # ── 2. COSY: endpoint-only (no element-by-element available) ─────────
    cosy_json = results_dir / 'cosy_s1_fr3_postfix.json'
    if not cosy_json.exists():
        cosy_json = results_dir / 'cosy_s1_fr3_warm.json'
    cosy_endpoint = None
    if cosy_json.exists():
        _, cosy_data_fr3 = load_currents(str(cosy_json))
        cosy_endpoint = cosy_data_fr3.get('twiss_undulator', {})
        cosy_mse = cosy_data_fr3.get('mse', 0)
        print(f"\n[2/3] COSY FR3 endpoint ({cosy_json.name}): "
              f"bx={cosy_endpoint.get('beta_x',0):.4f}, "
              f"by={cosy_endpoint.get('beta_y',0):.4f}, "
              f"ax={cosy_endpoint.get('alpha_x',0):.4f}, "
              f"ay={cosy_endpoint.get('alpha_y',0):.4f}, RMS={math.sqrt(cosy_mse):.2e}")
        print("  (Element-by-element Twiss not available from COSY — endpoint only)")

    # ── 3. RF-Track: FELsim stages 1-10 + RF-Track stage 11, tracked ─────
    if not args.skip_rftrack:
        print("\n[3/3] RF-Track: FELsim base + RF-Track stage 11 optimization...")
        # Build RF-Track currents: FELsim stages 1-10, RF-Track stage 11
        rftrack_currents = dict(felsim_currents)
        # RF-Track optimised stage 11 for eps_n=8 (post-bug-fix, from comparison.csv)
        rftrack_currents[87] = 0.2563250425674134
        rftrack_currents[93] = 2.0890902626402377
        rftrack_currents[95] = 4.346161882778803
        rftrack_currents[97] = 7.039697750925562
        rftrack_data = run_rftrack_evolution(rftrack_currents, beam_dist)
        if rftrack_data:
            rf = rftrack_data[-1]
            rf_mse = (
                (rf['beta_x'] - params['beta_xm'])**2
                + (rf['beta_y'] - params['beta_ym'])**2
                + (rf['alpha_x'] - params['alpha_xm'])**2
                + (rf['alpha_y'] - params['alpha_ym'])**2
            ) / 4
            print(f"  bx={rf['beta_x']:.4f}, by={rf['beta_y']:.4f}, "
                  f"ax={rf['alpha_x']:.4f}, ay={rf['alpha_y']:.4f}, RMS={math.sqrt(rf_mse):.2e}")
            all_results.append({
                'label': 'RF-Track (particle tracking, hybrid NM)',
                'data': rftrack_data, 'mse': rf_mse,
            })

    # Summary
    print(f"\n{'='*72}")
    print(f"  {'Code':<40} {'beta_x':>8} {'beta_y':>8} {'alpha_x':>8} {'alpha_y':>8} {'RMS':>10}")
    print(f"  {'Target':<40} {params['beta_xm']:8.4f} {params['beta_ym']:8.4f} "
          f"{params['alpha_xm']:8.4f} {params['alpha_ym']:8.4f}")
    print(f"  {'-'*68}")
    for res in all_results:
        d = res['data'][-1]
        print(f"  {res['label']:<40} {d['beta_x']:8.4f} {d['beta_y']:8.4f} "
              f"{d['alpha_x']:8.4f} {d['alpha_y']:8.4f} {math.sqrt(res['mse']):10.2e}")
    if cosy_endpoint:
        print(f"  {'COSY FR3 (endpoint only)':<40} "
              f"{cosy_endpoint['beta_x']:8.4f} {cosy_endpoint['beta_y']:8.4f} "
              f"{cosy_endpoint['alpha_x']:8.4f} {cosy_endpoint['alpha_y']:8.4f} "
              f"{math.sqrt(cosy_mse):10.2e}")
    print(f"{'='*72}")

    # Single-figure overlay plot
    output = args.output or str(results_dir / 'R3' / 'beta_comparison.png')
    Path(output).parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)

    # Collect all s values for consistent x-axis
    s_end = 0

    # FELsim
    fp = all_results[0]['data']
    s_fp = [r['s'] for r in fp]
    s_end = s_fp[-1]
    axes[0].plot(s_fp, [r['beta_x'] for r in fp],
                 'b-', lw=1.5, label=r'FELsim $\beta_x$')
    axes[0].plot(s_fp, [r['beta_y'] for r in fp],
                 'r-', lw=1.5, label=r'FELsim $\beta_y$')
    axes[1].plot(s_fp, [r['alpha_x'] for r in fp],
                 'b-', lw=1.5, label=r'FELsim $\alpha_x$')
    axes[1].plot(s_fp, [r['alpha_y'] for r in fp],
                 'r-', lw=1.5, label=r'FELsim $\alpha_y$')

    # RF-Track (if available and sensible)
    if len(all_results) > 1 and all_results[1].get('mse', 999) < 1:
        rf = all_results[1]['data']
        s_rf = [r['s'] for r in rf]
        axes[0].plot(s_rf, [r['beta_x'] for r in rf],
                     'b:', lw=1.8, alpha=0.7, label=r'RF-Track $\beta_x$')
        axes[0].plot(s_rf, [r['beta_y'] for r in rf],
                     'r:', lw=1.8, alpha=0.7, label=r'RF-Track $\beta_y$')
        axes[1].plot(s_rf, [r['alpha_x'] for r in rf],
                     'b:', lw=1.8, alpha=0.7, label=r'RF-Track $\alpha_x$')
        axes[1].plot(s_rf, [r['alpha_y'] for r in rf],
                     'r:', lw=1.8, alpha=0.7, label=r'RF-Track $\alpha_y$')

    # COSY endpoint markers
    if cosy_endpoint:
        axes[0].plot(s_end, cosy_endpoint['beta_x'], 'gv', ms=12, zorder=5,
                     label=f"COSY FR3 $\\beta_x$ = {cosy_endpoint['beta_x']:.4f}")
        axes[0].plot(s_end, cosy_endpoint['beta_y'], 'g^', ms=12, zorder=5,
                     label=f"COSY FR3 $\\beta_y$ = {cosy_endpoint['beta_y']:.4f}")
        axes[1].plot(s_end, cosy_endpoint['alpha_x'], 'gv', ms=12, zorder=5,
                     label=f"COSY FR3 $\\alpha_x$ = {cosy_endpoint['alpha_x']:.4f}")
        axes[1].plot(s_end, cosy_endpoint['alpha_y'], 'g^', ms=12, zorder=5,
                     label=f"COSY FR3 $\\alpha_y$ = {cosy_endpoint['alpha_y']:.4f}")

    # Targets
    axes[0].axhline(params['beta_xm'], color='b', ls=':', alpha=0.3)
    axes[0].axhline(params['beta_ym'], color='r', ls=':', alpha=0.3)
    axes[0].set_ylabel(r'$\beta$ (m)', fontsize=12)
    axes[0].set_ylim(0, 50)
    axes[0].legend(fontsize=8, ncol=2, loc='upper right')
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title(
        r'Beta functions along beamline — each code with its own optimised currents'
        '\n(post-bug-fix; COSY = endpoint only)',
        fontsize=11)

    axes[1].axhline(params['alpha_xm'], color='b', ls=':', alpha=0.3)
    axes[1].axhline(params['alpha_ym'], color='r', ls=':', alpha=0.3)
    axes[1].set_ylabel(r'$\alpha$', fontsize=12)
    axes[1].set_ylim(-20, 20)
    axes[1].set_xlabel('s (m)', fontsize=12)
    axes[1].legend(fontsize=8, ncol=2, loc='upper right')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: {output}")
