#!/usr/bin/env python3
"""Generate seminar-quality figures for the UH FEL transport line.

Figures:
  1. Twiss parameter evolution beta(s), alpha(s) with element-type strip
  2. Transverse phase space at the undulator entrance
  3. Cross-code Twiss matching accuracy (FELsim vs COSY INFINITY)
  4. Multi-stage optimization convergence
  5. Beam envelope sigma(s) and dispersion D(s) with element strip

Author: Eremey Valetov
"""

import sys
import math
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

from ebeam import beam
from beamline import lattice, qpfLattice, qpdLattice, dipole, dipole_wedge
from excelElements import ExcelElements
from beamOptimizer import beamOptimizer

# ── Paths ──────────────────────────────────────────────────────────────────
EXCEL_PATH = backend_dir.parent / 'beam_excel' / 'Beamline_elements.xlsx'
RESULTS_DIR = Path(__file__).resolve().parent / 'results'
OUTPUT_DIR = RESULTS_DIR / 'seminar'

# ── Beam parameters ───────────────────────────────────────────────────────
ENERGY = 40        # MeV
EPSILON_N = 8      # pi.mm.mrad (normalised)
X_STD = 0.8        # mm
Y_STD = 0.8        # mm
FREQ = 2856e6      # Hz
BUNCH_SPREAD = 2   # ps
ENERGY_STD_PCT = 0.5
H_CHIRP = 5e9
NB_PARTICLES = 5000

# ── Undulator targets ─────────────────────────────────────────────────────
K_UND = 1.2
LAMBDA_U = 2.3e-2  # m
QUAD_INDICES = [1, 3, 10, 16, 18, 20, 27, 33, 35, 37, 39, 41, 43,
                50, 56, 58, 61, 63, 70, 76, 78, 80, 87, 93, 95, 97]
MIRROR = {39: 37, 41: 35, 43: 33}

# ── Stage labels for convergence plot ─────────────────────────────────────
STAGE_LABELS = [
    'Doublet',         'Chicane 1 corr.', 'Triplet A',
    'Chicane 2 corr.', 'Symm. triplet',   'Chicane 3 corr.',
    'Doublet B',       'Triplet B',       'Chicane 4 corr.',
    'Triplet C',       'Undulator match',
]

# ── Colour-blind-friendly palette (Tol bright) ───────────────────────────
C_BLUE = '#0077BB'
C_RED  = '#CC3311'
C_TEAL = '#009988'
C_ORANGE = '#EE7733'
C_GREY = '#BBBBBB'

# ── Seminar rcParams ─────────────────────────────────────────────────────
RCPARAMS = {
    'font.family': 'serif',
    'font.size': 13,
    'axes.labelsize': 15,
    'axes.titlesize': 15,
    'legend.fontsize': 10,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.grid': True,
    'grid.alpha': 0.3,
    'lines.linewidth': 2.5,
}


# ── Helpers ───────────────────────────────────────────────────────────────

def compute_params():
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    norm = relat.gamma * relat.beta
    epsilon = EPSILON_N / norm
    beta_ym = relat.gamma * LAMBDA_U / (2 * np.pi * K_UND)
    beta_0 = X_STD**2 / epsilon
    return {
        'beta_xm': 1.4, 'alpha_xm': 0.47,
        'beta_ym': beta_ym, 'alpha_ym': 0.0,
        'epsilon': epsilon, 'beta_0': beta_0,
        'gamma': relat.gamma, 'beta_rel': relat.beta, 'norm': norm,
    }


def generate_beam(params, n=NB_PARTICLES):
    np.random.seed(42)
    eps = params['epsilon']
    ebeam_gen = beam()
    dist = ebeam_gen.gen_6d_gaussian(
        0,
        [X_STD, eps / X_STD, Y_STD, eps / Y_STD,
         BUNCH_SPREAD * 1e-9 * FREQ, ENERGY_STD_PCT * 10],
        n,
    )
    dist[:, 5] += H_CHIRP * dist[:, 4] / FREQ
    return dist


def load_currents(json_path):
    with open(json_path) as f:
        data = json.load(f)
    currents = {int(k): v for k, v in data['currents'].items()}
    for dst, src in MIRROR.items():
        if src in currents and dst not in currents:
            currents[dst] = currents[src]
    return currents, data


def build_line():
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    excel = ExcelElements(str(EXCEL_PATH))
    line = relat.changeBeamType("electron", ENERGY, excel.create_beamline())[:118]
    return line


def _run_stage(opti, method, segmentVar, startPoint, objectives):
    """Run one optimizer stage and return its convergence data."""
    opti.calc(method, segmentVar, startPoint, objectives)
    return list(opti.plotMSE), list(opti.plotIterate)


def run_11_stage_optimization(line, beam_dist, params):
    """Run the 11-stage Nelder-Mead optimization. Returns (currents, convergence_data)."""
    axm, aym = params['alpha_xm'], params['alpha_ym']
    bxm, bym = params['beta_xm'], params['beta_ym']

    opti = beamOptimizer(line, beam_dist)
    all_mse = []
    stage_boundaries = [0]

    stages = [
        # (segmentVar, startPoint, objectives)
        ({1: ["I", "current", lambda n: n], 3: ["I2", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}, "I2": {"bounds": (0, 10), "start": 1}},
         {8: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
              {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.0}],
          9: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
              {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]}),
        ({10: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {15: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}),
        ({16: ["I", "current", lambda n: n], 18: ["I2", "current", lambda n: n],
          20: ["I3", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 5},
          "I3": {"bounds": (0, 10), "start": 3}},
         {25: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
          26: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]}),
        ({27: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {32: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}),
        ({37: ["I", "current", lambda n: n], 35: ["I2", "current", lambda n: n],
          33: ["I3", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2},
          "I3": {"bounds": (0, 10), "start": 2}},
         {37: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["x", "envelope"], "goal": 2.0, "weight": 1},
               {"measure": ["y", "envelope"], "goal": 2.0, "weight": 1}]}),
    ]

    for seg, sp, obj in stages:
        mse, itr = _run_stage(opti, "Nelder-Mead", seg, sp, obj)
        all_mse.extend(mse)
        stage_boundaries.append(len(all_mse))

    # Mirror symmetry
    line[43].current = line[33].current
    line[41].current = line[35].current
    line[39].current = line[37].current

    stages_2 = [
        ({50: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {55: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}),
        ({56: ["I", "current", lambda n: n], 58: ["I2", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2}},
         {59: [{"measure": ["x", "envelope"], "goal": 0.0, "weight": 1},
               {"measure": ["y", "envelope"], "goal": 0.0, "weight": 1}]}),
        ({61: ["I", "current", lambda n: n], 63: ["I2", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2}},
         {68: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
          69: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]}),
        ({70: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {75: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}),
        ({76: ["I", "current", lambda n: n], 78: ["I2", "current", lambda n: n],
          80: ["I3", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2},
          "I3": {"bounds": (0, 10), "start": 2}},
         {85: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
          86: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]}),
        ({87: ["Ic", "current", lambda n: n], 93: ["I", "current", lambda n: n],
          95: ["I2", "current", lambda n: n], 97: ["I3", "current", lambda n: n]},
         {"Ic": {"bounds": (0, 10), "start": 4}, "I": {"bounds": (0, 10), "start": 2},
          "I2": {"bounds": (0, 10), "start": 2}, "I3": {"bounds": (0, 10), "start": 2}},
         {92: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 0.5}],
          117: [{"measure": ["x", "alpha"], "goal": axm, "weight": 1},
                {"measure": ["y", "alpha"], "goal": aym, "weight": 1},
                {"measure": ["x", "beta"], "goal": bxm, "weight": 1},
                {"measure": ["y", "beta"], "goal": bym, "weight": 1}]}),
    ]

    for seg, sp, obj in stages_2:
        mse, itr = _run_stage(opti, "Nelder-Mead", seg, sp, obj)
        all_mse.extend(mse)
        stage_boundaries.append(len(all_mse))

    currents = {idx: line[idx].current for idx in QUAD_INDICES}
    convergence = {'mse': all_mse, 'boundaries': stage_boundaries}
    return currents, convergence


def propagate_particles_full(line, beam_dist):
    """Track particles; return (evolution, final_particles).

    Evolution includes sigma and dispersion for Figure 5.
    """
    ebeam_calc = beam()
    particles = beam_dist.copy()
    s = 0.0
    out = []
    for elem in line:
        particles = np.array(elem.useMatrice(particles))
        s += getattr(elem, 'length', getattr(elem, 'L', 0))
        _, _, tw = ebeam_calc.cal_twiss(particles, ddof=1)
        out.append({
            's': s,
            'beta_x': tw.loc['x', r'$\beta$ (m)'],
            'beta_y': tw.loc['y', r'$\beta$ (m)'],
            'alpha_x': tw.loc['x', r'$\alpha$'],
            'alpha_y': tw.loc['y', r'$\alpha$'],
            'sigma_x': np.std(particles[:, 0], ddof=1),
            'sigma_y': np.std(particles[:, 2], ddof=1),
            'disp_x': tw.loc['x', r'$D$ (m)'],
        })
    return out, particles


def element_spans(line):
    """Return list of (s_start, s_end, type_str) for non-drift elements."""
    spans = []
    s = 0.0
    for elem in line:
        length = getattr(elem, 'length', getattr(elem, 'L', 0))
        s_end = s + length
        if isinstance(elem, (qpfLattice, qpdLattice)):
            spans.append((s, s_end, 'quad'))
        elif isinstance(elem, dipole):
            spans.append((s, s_end, 'dipole'))
        elif isinstance(elem, dipole_wedge):
            spans.append((s, s_end, 'wedge'))
        s = s_end
    return spans


def draw_element_strip(ax, spans):
    """Draw a thin element-type strip on a dedicated axes."""
    colors = {'quad': C_BLUE, 'dipole': C_RED, 'wedge': C_ORANGE}
    ax.set_ylim(0, 1)
    for s0, s1, etype in spans:
        ax.axvspan(s0, s1, color=colors[etype], alpha=0.7)
    ax.set_yticks([])
    ax.set_ylabel('', fontsize=1)
    ax.legend(
        handles=[Patch(color=C_BLUE, alpha=0.7, label='Quadrupole'),
                 Patch(color=C_RED, alpha=0.7, label='Dipole'),
                 Patch(color=C_ORANGE, alpha=0.7, label='Wedge')],
        fontsize=8, ncol=3, loc='upper right', framealpha=0.8,
    )


def _save(fig, name):
    for fmt in ('pdf', 'png'):
        fig.savefig(OUTPUT_DIR / f'{name}.{fmt}')
    plt.close(fig)
    print(f"  -> {name}.pdf/png")


# ── Figure 1: Twiss evolution ────────────────────────────────────────────

def figure_twiss_evolution(params, twiss, line):
    print("Figure 1: Twiss evolution...")
    _, cosy_data = load_currents(str(RESULTS_DIR / 'cosy_s1_fr3_postfix.json'))
    tw_cosy = cosy_data['twiss_undulator']
    spans = element_spans(line)

    s = [r['s'] for r in twiss]
    s_end = s[-1]
    bx_data = [r['beta_x'] for r in twiss]
    by_data = [r['beta_y'] for r in twiss]
    ax_data = [r['alpha_x'] for r in twiss]
    ay_data = [r['alpha_y'] for r in twiss]
    beta_ylim = max(max(bx_data), max(by_data)) * 1.15
    alpha_ylim = max(max(abs(v) for v in ax_data),
                     max(abs(v) for v in ay_data)) * 1.15

    # Layout: 3 rows × 2 cols; element strip spans full width,
    # Twiss panels get a main (wide) + zoom (narrow) column
    fig = plt.figure(figsize=(17, 9))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 8, 8],
                          width_ratios=[3, 1],
                          hspace=0.03, wspace=0.30)

    # ── Element strip (spans both columns) ──
    ax_strip = fig.add_subplot(gs[0, 0])
    draw_element_strip(ax_strip, spans)
    # Hide the top-right cell
    ax_dummy = fig.add_subplot(gs[0, 1])
    ax_dummy.set_visible(False)

    # ── Helper to plot Twiss data on an axes ──
    def _plot_beta(ax, show_xlabel=False, show_legend=True, zoom=False):
        ax.plot(s, bx_data, color=C_BLUE, ls='-', label=r'$\beta_x$ (FELsim)')
        ax.plot(s, by_data, color=C_RED, ls='--', label=r'$\beta_y$ (FELsim)')
        ax.plot(s_end, tw_cosy['beta_x'], 'v', color=C_TEAL, ms=10, zorder=5,
                label=rf"$\beta_x$ COSY = {tw_cosy['beta_x']:.4f} m")
        ax.plot(s_end, tw_cosy['beta_y'], '^', color=C_TEAL, ms=10, zorder=5,
                label=rf"$\beta_y$ COSY = {tw_cosy['beta_y']:.4f} m")
        ax.axhline(params['beta_xm'], color=C_BLUE, ls=':', alpha=0.5, lw=1.2)
        ax.axhline(params['beta_ym'], color=C_RED, ls=':', alpha=0.5, lw=1.2)
        ax.set_ylabel(r'$\beta$ (m)')
        if show_xlabel:
            ax.set_xlabel('s (m)')
        if show_legend:
            ax.legend(fontsize=9, ncol=2, loc='upper left', framealpha=0.9)

    def _plot_alpha(ax, show_xlabel=True, show_legend=True, zoom=False):
        ax.plot(s, ax_data, color=C_BLUE, ls='-', label=r'$\alpha_x$ (FELsim)')
        ax.plot(s, ay_data, color=C_RED, ls='--', label=r'$\alpha_y$ (FELsim)')
        ax.plot(s_end, tw_cosy['alpha_x'], 'v', color=C_TEAL, ms=10, zorder=5,
                label=rf"$\alpha_x$ COSY = {tw_cosy['alpha_x']:.4f}")
        ax.plot(s_end, tw_cosy['alpha_y'], '^', color=C_TEAL, ms=10, zorder=5,
                label=rf"$\alpha_y$ COSY = {tw_cosy['alpha_y']:.5f}")
        ax.axhline(params['alpha_xm'], color=C_BLUE, ls=':', alpha=0.5, lw=1.2)
        ax.axhline(params['alpha_ym'], color=C_RED, ls=':', alpha=0.5, lw=1.2)
        ax.set_ylabel(r'$\alpha$')
        if show_xlabel:
            ax.set_xlabel('s (m)')
        if show_legend:
            ax.legend(fontsize=9, ncol=2, loc='upper left', framealpha=0.9)

    # ── Beta: main panel ──
    ax_bm = fig.add_subplot(gs[1, 0])
    _plot_beta(ax_bm, show_xlabel=False)
    ax_bm.set_ylim(0, beta_ylim)
    ax_bm.axvspan(10.8, s_end + 0.15, color=C_TEAL, alpha=0.04, zorder=0)

    # ── Beta: zoom panel ──
    ax_bz = fig.add_subplot(gs[1, 1])
    _plot_beta(ax_bz, show_legend=False)
    ax_bz.set_xlim(10.8, s_end + 0.15)
    ax_bz.set_ylim(-0.1, 2.8)
    ax_bz.set_title('Undulator matching', fontsize=11)
    ax_bz.tick_params(labelsize=9)
    ax_bz.set_xlabel('s (m)', fontsize=10)
    ax_bz.annotate(rf"target $\beta_x$ = {params['beta_xm']:.1f}",
                    xy=(11.9, params['beta_xm'] + 0.06),
                    fontsize=8, color=C_BLUE, alpha=0.8, va='bottom')
    ax_bz.annotate(rf"target $\beta_y$ = {params['beta_ym']:.4f}",
                    xy=(10.85, params['beta_ym'] + 0.08),
                    fontsize=8, color=C_RED, alpha=0.8, va='bottom')

    # ── Alpha: main panel ──
    ax_am = fig.add_subplot(gs[2, 0], sharex=ax_bm)
    _plot_alpha(ax_am)
    ax_am.set_ylim(-alpha_ylim, alpha_ylim)
    ax_am.axvspan(10.8, s_end + 0.15, color=C_TEAL, alpha=0.04, zorder=0)

    # ── Alpha: zoom panel ──
    ax_az = fig.add_subplot(gs[2, 1])
    _plot_alpha(ax_az, show_legend=False)
    ax_az.set_xlim(10.8, s_end + 0.15)
    ax_az.set_ylim(-1.5, 4.0)
    ax_az.set_title('Undulator matching', fontsize=11)
    ax_az.tick_params(labelsize=9)
    ax_az.annotate(rf"target $\alpha_x$ = {params['alpha_xm']:.2f}",
                    xy=(s_end + 0.1, params['alpha_xm']),
                    fontsize=8, color=C_BLUE, alpha=0.8, va='bottom')
    ax_az.annotate(rf"target $\alpha_y$ = {params['alpha_ym']:.1f}",
                    xy=(s_end + 0.1, params['alpha_ym'] + 0.15),
                    fontsize=8, color=C_RED, alpha=0.8, va='bottom')

    fig.suptitle(
        r'Twiss parameter evolution — UH FEL transport line'
        '\n'
        r'26 quadrupole currents optimised for undulator matching '
        f'(E = {ENERGY} MeV, '
        r'$\varepsilon_n$ = '
        f'{EPSILON_N} '
        r'$\pi{\cdot}$mm${\cdot}$mrad)',
        fontsize=14, y=0.995,
    )

    _save(fig, 'twiss_evolution')
    final = twiss[-1]
    print(f"  Final: bx={final['beta_x']:.4f}, by={final['beta_y']:.4f}, "
          f"ax={final['alpha_x']:.4f}, ay={final['alpha_y']:.4f}")


# ── Figure 2: Phase space at undulator ────────────────────────────────────

def figure_phase_space(final_particles):
    print("Figure 2: Phase space at undulator...")
    ebeam_calc = beam()
    _, _, tw = ebeam_calc.cal_twiss(final_particles, ddof=1)

    fig, axes = plt.subplots(1, 3, figsize=(14, 6.5),
                              gridspec_kw={'width_ratios': [1, 1, 0.05]})
    theta = np.linspace(0, 2 * np.pi, 300)

    for i, (plane, p1, p2, lbl1, lbl2) in enumerate([
        ('x', 0, 1, 'x (mm)', "x' (mrad)"),
        ('y', 2, 3, 'y (mm)', "y' (mrad)"),
    ]):
        ax = axes[i]
        beta = tw.loc[plane, r'$\beta$ (m)']
        alpha = tw.loc[plane, r'$\alpha$']
        eps = tw.loc[plane, r'$\epsilon$ ($\pi$.mm.mrad)']

        sc = ax.scatter(final_particles[:, p1], final_particles[:, p2],
                        c=final_particles[:, 5], cmap='coolwarm',
                        s=2.5, alpha=0.45, rasterized=True, zorder=1)

        u = np.sqrt(eps * beta) * np.cos(theta)
        v = np.sqrt(eps / beta) * (-alpha * np.cos(theta) + np.sin(theta))
        ax.plot(u, v, color='black', lw=1.5, alpha=0.85,
                label=rf'$1\sigma$ ellipse', zorder=3)
        ax.plot(2*u, 2*v, color='black', lw=1, ls='--', alpha=0.5,
                label=rf'$2\sigma$ ellipse', zorder=3)

        textstr = (rf'$\beta_{plane}$ = {beta:.4f} m'
                   '\n'
                   rf'$\alpha_{plane}$ = {alpha:.4f}'
                   '\n'
                   rf'$\varepsilon_{plane}$ = {eps:.3f} $\pi\cdot$mm$\cdot$mrad')
        ax.text(0.03, 0.97, textstr, transform=ax.transAxes,
                fontsize=9, va='top', ha='left',
                bbox=dict(boxstyle='round,pad=0.4', fc='white', alpha=0.85))

        ax.set_xlabel(lbl1)
        ax.set_ylabel(lbl2)
        title = 'Horizontal' if plane == 'x' else 'Vertical'
        ax.set_title(f'{title} phase space')
        ax.legend(fontsize=9, loc='lower right')

        # Tight symmetric limits from data
        xdata, ydata = final_particles[:, p1], final_particles[:, p2]
        xext = max(abs(xdata.min()), abs(xdata.max())) * 1.15
        yext = max(abs(ydata.min()), abs(ydata.max())) * 1.15
        ax.set_xlim(-xext, xext)
        ax.set_ylim(-yext, yext)

    # Colorbar in dedicated right axes
    fig.colorbar(sc, cax=axes[2], label=r'$\Delta K/K_0 \times 10^3$')

    # Note energy-position correlation in horizontal plane (from residual dispersion)
    corr_xE = np.corrcoef(final_particles[:, 0], final_particles[:, 5])[0, 1]
    if abs(corr_xE) > 0.05:
        axes[0].text(0.03, 0.03,
                     rf'$\rho(x, \delta)$ = {corr_xE:.2f} (residual $D_x$)',
                     transform=axes[0].transAxes, fontsize=8, va='bottom',
                     color='grey', alpha=0.7)

    fig.suptitle(
        'Transverse phase space at undulator entrance\n'
        f'FELsim particle tracking, N = {NB_PARTICLES}, '
        f'E = {ENERGY} MeV',
        fontsize=14,
    )
    plt.tight_layout()
    _save(fig, 'phase_space_undulator')


# ── Figure 3: Cross-code accuracy comparison ─────────────────────────────

def figure_accuracy(params, felsim_twiss_final):
    print("Figure 3: Cross-code accuracy...")
    _, cosy_data = load_currents(str(RESULTS_DIR / 'cosy_s1_fr3_postfix.json'))
    tc = cosy_data['twiss_undulator']

    labels = [r'$\beta_x$ (m)', r'$\beta_y$ (m)', r'$\alpha_x$', r'$\alpha_y$']
    tgt = np.array([params['beta_xm'], params['beta_ym'],
                     params['alpha_xm'], params['alpha_ym']])
    fel = np.array([felsim_twiss_final['beta_x'], felsim_twiss_final['beta_y'],
                     felsim_twiss_final['alpha_x'], felsim_twiss_final['alpha_y']])
    cosy = np.array([tc['beta_x'], tc['beta_y'], tc['alpha_x'], tc['alpha_y']])

    felsim_mse = np.mean((fel - tgt)**2)
    cosy_mse = cosy_data.get('mse', np.mean((cosy - tgt)**2))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5),
                                    gridspec_kw={'width_ratios': [3, 2]})

    x = np.arange(len(labels))
    w = 0.22
    ax1.bar(x - w, tgt, w, label='Target', color=C_GREY, edgecolor='black', lw=0.5)
    ax1.bar(x, fel, w,
            label=f'FELsim (RMS = {math.sqrt(felsim_mse):.1e})',
            color=C_BLUE, alpha=0.85, edgecolor='black', lw=0.5)
    ax1.bar(x + w, cosy, w,
            label=f'COSY (RMS = {math.sqrt(cosy_mse):.1e})',
            color=C_ORANGE, alpha=0.85, edgecolor='black', lw=0.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_ylabel('Value')
    ax1.legend(fontsize=10, loc='upper right')
    ax1.set_title('Twiss parameters at undulator', fontsize=13)

    for bar_set in [(x - w, tgt), (x, fel), (x + w, cosy)]:
        for xi, val in zip(*bar_set):
            if abs(val) > 0.005:
                ax1.annotate(f'{val:.4f}', xy=(xi, val),
                             xytext=(0, 4), textcoords='offset points',
                             ha='center', va='bottom', fontsize=9, rotation=40)

    dev_fel = np.maximum(np.abs(fel - tgt), 1e-12)
    dev_cosy = np.maximum(np.abs(cosy - tgt), 1e-12)
    w2 = 0.3
    ax2.bar(x - w2/2, dev_fel, w2, label='FELsim', color=C_BLUE,
            alpha=0.85, edgecolor='black', lw=0.5)
    ax2.bar(x + w2/2, dev_cosy, w2, label='COSY', color=C_ORANGE,
            alpha=0.85, edgecolor='black', lw=0.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels)
    ax2.set_ylabel('|Achieved - Target|')
    ax2.set_yscale('log')
    ax2.legend(fontsize=10)
    ax2.set_title('Absolute deviation from target', fontsize=13)

    fig.suptitle(
        'Cross-code Twiss matching accuracy — FELsim vs COSY INFINITY\n'
        f'E = {ENERGY} MeV, 26 quadrupole currents optimised',
        fontsize=14, y=1.01,
    )
    plt.tight_layout()
    _save(fig, 'accuracy_comparison')

    print(f"\n  {'Parameter':<10} {'Target':>10} {'FELsim':>10} {'COSY':>10}")
    print(f"  {'-'*42}")
    for l, t, f, c in zip(labels, tgt, fel, cosy):
        clean = l.replace('$', '').replace('\\', '')
        print(f"  {clean:<10} {t:10.6f} {f:10.6f} {c:10.6f}")
    print(f"\n  FELsim RMS: {math.sqrt(felsim_mse):.2e}   COSY RMS: {math.sqrt(cosy_mse):.2e}")


# ── Figure 4: Optimization convergence ───────────────────────────────────

def figure_convergence(convergence):
    print("Figure 4: Optimization convergence...")
    mse = convergence['mse']
    boundaries = convergence['boundaries']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 7.5),
                                    gridspec_kw={'height_ratios': [3, 1.2], 'hspace': 0.35})

    # ── Top: MSE vs iteration with stage boundaries ──
    iterations = np.arange(1, len(mse) + 1)
    ax1.semilogy(iterations, mse, color=C_BLUE, lw=0.8, alpha=0.6)

    cmap = plt.colormaps['tab10']
    for i in range(len(boundaries) - 1):
        i0, i1 = boundaries[i], boundaries[i + 1]
        if i1 > i0:
            ax1.axvspan(i0 + 0.5, i1 + 0.5, color=cmap(i % 10), alpha=0.08)

    ax1.set_ylabel('Stage objective')
    ax1.set_xlabel('Cumulative iteration')
    ax1.set_title('11-stage Nelder-Mead optimization convergence', fontsize=14)

    # Stage labels inside top of plot (use axes fraction for y)
    for i in range(len(boundaries) - 1):
        i0, i1 = boundaries[i], boundaries[i + 1]
        if i1 > i0:
            mid = (i0 + i1) / 2 + 0.5
            ax1.annotate(f'S{i+1}', xy=(mid, 0.97),
                         xycoords=('data', 'axes fraction'),
                         fontsize=7, ha='center', va='top', alpha=0.55)

    # ── Bottom: best MSE per stage (bar chart) ──
    best_per_stage = []
    for i in range(len(boundaries) - 1):
        i0, i1 = boundaries[i], boundaries[i + 1]
        if i1 > i0:
            best_per_stage.append(min(mse[i0:i1]))
        else:
            best_per_stage.append(np.nan)

    x_stages = np.arange(len(best_per_stage))
    colours = [cmap(i % 10) for i in x_stages]
    ax2.bar(x_stages, best_per_stage, color=colours, alpha=0.8,
            edgecolor='black', lw=0.5)
    ax2.set_yscale('log')
    ax2.set_xticks(x_stages)
    # Use stage names + numbers as x-tick labels
    tick_labels = [f'S{i+1}\n{STAGE_LABELS[i]}' if i < len(STAGE_LABELS) else f'S{i+1}'
                   for i in x_stages]
    ax2.set_xticklabels(tick_labels, fontsize=7.5)
    ax2.set_ylabel('Best objective')
    ax2.set_title('Best objective per stage', fontsize=13)

    _save(fig, 'convergence')


# ── Figure 5: Beam envelope and dispersion ───────────────────────────────

def figure_envelope_dispersion(twiss, line):
    print("Figure 5: Beam envelope and dispersion...")
    spans = element_spans(line)
    s = [r['s'] for r in twiss]

    fig, axes = plt.subplots(
        3, 1, figsize=(14, 8.5), sharex=True,
        gridspec_kw={'height_ratios': [1, 6, 6], 'hspace': 0.03},
    )

    draw_element_strip(axes[0], spans)

    # ── Beam envelope ──
    ax = axes[1]
    ax.plot(s, [r['sigma_x'] for r in twiss], color=C_BLUE, ls='-',
            label=r'$\sigma_x$ (horizontal)')
    ax.plot(s, [r['sigma_y'] for r in twiss], color=C_RED, ls='--',
            label=r'$\sigma_y$ (vertical)')
    ax.set_ylabel(r'RMS beam size $\sigma$ (mm)')
    ax.set_ylim(0, max(max(r['sigma_x'] for r in twiss),
                       max(r['sigma_y'] for r in twiss)) * 1.15)
    ax.legend(fontsize=10, loc='upper left', framealpha=0.9)

    # ── Dispersion ──
    ax = axes[2]
    disp_data = [r['disp_x'] for r in twiss]
    ax.plot(s, disp_data, color=C_BLUE, ls='-',
            label=r'$D_x$ (horizontal)')
    ax.axhline(0, color='black', ls='-', lw=0.5, alpha=0.3)
    ax.set_ylabel(r'Dispersion $D_x$ (m)')
    ax.set_xlabel('s (m)')
    ax.legend(fontsize=10, loc='upper left', framealpha=0.9)
    # Annotate residual dispersion at undulator
    dx_final = disp_data[-1]
    ax.annotate(rf'$D_x$ = {dx_final:.3f} m',
                xy=(s[-1], dx_final), xytext=(-120, 40),
                textcoords='offset points', fontsize=10,
                arrowprops=dict(arrowstyle='->', color=C_BLUE, lw=1.5),
                color=C_BLUE, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none', alpha=0.85))

    fig.suptitle(
        r'Beam envelope and dispersion — UH FEL transport line'
        '\n'
        f'Optimised 26-quadrupole transport '
        f'(E = {ENERGY} MeV)',
        fontsize=14, y=0.995,
    )

    _save(fig, 'envelope_dispersion')


# ── Main ──────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(RCPARAMS)
    params = compute_params()

    print("UH FEL Seminar Figures")
    print(f"  Targets: beta_x={params['beta_xm']}, beta_y={params['beta_ym']:.4f}, "
          f"alpha_x={params['alpha_xm']}, alpha_y={params['alpha_ym']}")
    print(f"  Output: {OUTPUT_DIR}\n")

    line = build_line()
    beam_dist = generate_beam(params)

    print("Running 11-stage Nelder-Mead optimization...")
    currents, convergence = run_11_stage_optimization(line, beam_dist, params)
    print(f"  Optimised {len(currents)} quadrupole currents")

    print("Propagating particles through optimised beamline...")
    twiss_evol, final_particles = propagate_particles_full(line, beam_dist)

    felsim_mse = np.mean([
        (twiss_evol[-1]['beta_x'] - params['beta_xm'])**2,
        (twiss_evol[-1]['beta_y'] - params['beta_ym'])**2,
        (twiss_evol[-1]['alpha_x'] - params['alpha_xm'])**2,
        (twiss_evol[-1]['alpha_y'] - params['alpha_ym'])**2,
    ])
    updated = {
        'fringe_field_order': 0,
        'mse': float(felsim_mse),
        'currents': {str(k): float(v) for k, v in currents.items()},
    }
    with open(RESULTS_DIR / 'felsim_nm_warm.json', 'w') as f:
        json.dump(updated, f, indent=2)
    print(f"  Updated felsim_nm_warm.json (RMS = {math.sqrt(felsim_mse):.2e})\n")

    figure_twiss_evolution(params, twiss_evol, line)
    figure_phase_space(final_particles)
    figure_accuracy(params, twiss_evol[-1])
    figure_convergence(convergence)
    figure_envelope_dispersion(twiss_evol, line)

    print(f"\nDone. All figures saved to {OUTPUT_DIR}")
