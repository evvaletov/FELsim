#!/usr/bin/env python3
"""Generate figures for the IPAC 2026 (TUPM005) optimisation talk.

Two parts:

  Part I  — multi-code modelling stack (context / credibility)
    P1  Multi-code capability matrix
    P2  RF-Track linac model vs elegant (0.06% at peak)
    P3  Space-charge capability (DA-FMM vs Xsuite frozen/PIC3D)
    P4  Cross-code transport-line Twiss match quality

  Part II — the optimisation result (core narrative)
    F1  Twiss(s) with the TUPM005 undulator target band
    F2  Objective-design ablation: A/B/C failure-rate bar
    F3  Per-seed undulator-RMS distribution (basin structure)
    F4  Per-stage MSE divergence, good vs trapped seed
    F5  Twiss ellipse at the undulator, good vs trapped basin
    F6  GD vs BO outlook (NM real, BO template pending S6)

All Part I/II data is read from disk; nothing is re-simulated except the
lightweight single-pass particle tracking used by F1 and F5.

Author: Eremey Valetov
"""

import sys
import csv
import json
import glob
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Ellipse

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))

# Reuse the seminar script's lattice/beam helpers and constants — single
# source of truth for the model, not a re-implementation.
from generate_seminar_figures import (
    build_line, compute_params, generate_beam, load_currents,
    propagate_particles_full, element_spans,
    QUAD_INDICES, STAGE_LABELS,
    C_BLUE, C_RED, C_TEAL, C_ORANGE, C_GREY,
)

RESULTS_DIR = THIS_DIR / 'results'
ABL_DIR = RESULTS_DIR / 'ablation_TUPM005'
OUTPUT_DIR = RESULTS_DIR / 'ipac'
STYLE_PATH = THIS_DIR / 'felsim_talk.mplstyle'

FAIL_RMS = 1e-2          # undulator-RMS threshold separating the two basins
CONFIGS = ('A', 'B', 'C')
CONFIG_DESC = {
    'A': 'A\noriginal\nobjective',
    'B': 'B\nper-measure\nrescaling',
    'C': 'C\ntypo +\nfinite envelope',
}
CONFIG_COLOR = {'A': C_BLUE, 'B': C_TEAL, 'C': C_RED}


# ── Data loaders ──────────────────────────────────────────────────────────

def load_ablation():
    """Return {config: [seed_dict, ...]} for the TUPM005-target ablation."""
    out = {}
    for cfg in CONFIGS:
        seeds = []
        for f in sorted(glob.glob(str(ABL_DIR / f'{cfg}_seed*.json'))):
            with open(f) as fh:
                seeds.append(json.load(fh))
        out[cfg] = seeds
    return out


def _read_csv(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    cols = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    return cols


def line_at_seed(seed_json_path):
    """Build the lattice with the converged currents of one ablation run."""
    currents, _ = load_currents(str(seed_json_path))
    line = build_line()
    for idx in QUAD_INDICES:
        if idx < len(line) and idx in currents:
            line[idx].current = currents[idx]
    return line


def _save(fig, name):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for fmt in ('png', 'pdf'):
        fig.savefig(OUTPUT_DIR / f'{name}.{fmt}')
    plt.close(fig)
    print(f"  -> {name}.png/pdf")


# ── Part II ───────────────────────────────────────────────────────────────

def f1_twiss(params):
    """F1: matched Twiss(s) with the TUPM005 undulator target band."""
    print("F1: Twiss with TUPM005 target band...")
    best = min(glob.glob(str(ABL_DIR / 'A_seed*.json')),
               key=lambda f: json.load(open(f))['undulator_rms'])
    line = line_at_seed(best)
    beam_dist = generate_beam(params, seed=json.load(open(best))['seed'])
    twiss, _ = propagate_particles_full(line, beam_dist)
    spans = element_spans(line)

    s = np.array([r['s'] for r in twiss])
    bx = np.array([r['beta_x'] for r in twiss])
    by = np.array([r['beta_y'] for r in twiss])

    fig, (axs, ax) = plt.subplots(
        2, 1, figsize=(12.8, 7.2), height_ratios=[1, 14], sharex=True)
    axs.axis('off')
    for s0, s1, et in spans:
        c = {'quad': C_BLUE, 'dipole': C_RED, 'wedge': C_ORANGE}[et]
        axs.axvspan(s0, s1, color=c, alpha=0.7)
    axs.set_xlim(s[0], s[-1])
    axs.legend(handles=[Patch(color=C_BLUE, alpha=.7, label='Quad'),
                        Patch(color=C_RED, alpha=.7, label='Dipole'),
                        Patch(color=C_ORANGE, alpha=.7, label='Wedge')],
               ncol=3, loc='lower center', fontsize=15,
               bbox_to_anchor=(0.5, 1.0), frameon=False)

    ax.plot(s, bx, color=C_BLUE, label=r'$\beta_x$')
    ax.plot(s, by, color=C_RED, label=r'$\beta_y$')

    bxm, bym = params['beta_xm'], params['beta_ym']
    band = 0.10
    ax.axhspan(bxm * (1 - band), bxm * (1 + band), xmin=0.96,
               color=C_BLUE, alpha=0.25)
    ax.axhspan(bym * (1 - band), bym * (1 + band), xmin=0.96,
               color=C_RED, alpha=0.25)
    ax.scatter([s[-1]], [bxm], color=C_BLUE, marker='*', s=320, zorder=5,
               edgecolor='k', linewidth=0.6)
    ax.scatter([s[-1]], [bym], color=C_RED, marker='*', s=320, zorder=5,
               edgecolor='k', linewidth=0.6)
    ax.annotate(f'undulator target\n'
                rf'$\beta_x={bxm:.3f}$ m, $\alpha_x={params["alpha_xm"]:.3f}$',
                xy=(s[-1], bxm), xytext=(s[-1] * 0.62, bxm + 2.4),
                fontsize=17, ha='left',
                arrowprops=dict(arrowstyle='->', color='k', lw=1.5))

    ax.set_xlim(s[0], s[-1])
    ax.set_xlabel('s (m)')
    ax.set_ylabel(r'$\beta$ (m)')
    ax.legend(loc='upper left')
    fig.suptitle('UH FEL transport line — matched to the undulator',
                 fontsize=22)
    _save(fig, 'F1_twiss_target')


def f2_failrate(abl):
    """F2: A/B/C failure-rate bar — the headline."""
    print("F2: failure-rate bar...")
    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    rates, medians, labels = [], [], []
    for cfg in CONFIGS:
        rms = np.array([s['undulator_rms'] for s in abl[cfg]])
        ok = rms[rms <= FAIL_RMS]
        rates.append(100.0 * np.mean(rms > FAIL_RMS))
        medians.append(np.median(ok) if ok.size else np.nan)
        labels.append(cfg)

    bars = ax.bar(labels, rates,
                  color=[CONFIG_COLOR[c] for c in CONFIGS],
                  width=0.6, edgecolor='k', linewidth=1.2)
    for b, r, m, cfg in zip(bars, rates, medians, CONFIGS):
        n = len(abl[cfg])
        ax.text(b.get_x() + b.get_width() / 2, r + 2.5,
                f'{r:.0f}%\n({int(round(r/100*n))}/{n} seeds)',
                ha='center', va='bottom', fontsize=19, fontweight='bold')
        ax.text(b.get_x() + b.get_width() / 2, 4,
                f'success RMS\n{m:.1e}', ha='center', va='bottom',
                fontsize=14, color='white')

    ax.set_ylim(0, 100)
    ax.set_ylabel('Nelder–Mead failure rate (%)')
    ax.set_xlabel('Objective configuration')
    ax.set_xticks(range(len(CONFIGS)))
    ax.set_xticklabels([CONFIG_DESC[c] for c in CONFIGS], fontsize=15)
    ax.set_title('Objective design dominates optimiser robustness\n'
                 '(20 seeds/config, TUPM005 undulator targets)',
                 fontsize=21)
    _save(fig, 'F2_failrate')


def f3_distribution(abl):
    """F3: per-seed undulator-RMS distribution (basin structure)."""
    print("F3: per-seed RMS distribution...")
    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    rng = np.random.default_rng(0)
    for i, cfg in enumerate(CONFIGS):
        rms = np.array([s['undulator_rms'] for s in abl[cfg]])
        x = i + (rng.random(rms.size) - 0.5) * 0.28
        good = rms <= FAIL_RMS
        ax.scatter(x[good], rms[good], color=CONFIG_COLOR[cfg], s=90,
                   edgecolor='k', linewidth=0.5, zorder=3)
        ax.scatter(x[~good], rms[~good], color=CONFIG_COLOR[cfg], s=90,
                   marker='X', edgecolor='k', linewidth=0.5, zorder=3)
        ax.boxplot(rms, positions=[i], widths=0.55, showfliers=False,
                   medianprops=dict(color='k', lw=2))

    ax.axhline(FAIL_RMS, color=C_GREY, ls='--', lw=2)
    ax.text(2.45, FAIL_RMS * 1.4, 'basin threshold', fontsize=15,
            color='dimgray', ha='right')
    ax.set_yscale('log')
    ax.set_xticks(range(len(CONFIGS)))
    ax.set_xticklabels([CONFIG_DESC[c] for c in CONFIGS], fontsize=15)
    ax.set_ylabel('Final undulator RMS')
    ax.set_title('Outcomes are bimodal — a basin problem, not variance',
                 fontsize=21)
    ax.legend(handles=[
        plt.Line2D([], [], marker='o', ls='', color=C_GREY, mec='k',
                   label='converged (good basin)'),
        plt.Line2D([], [], marker='X', ls='', color=C_GREY, mec='k',
                   label='trapped (bad basin)')],
        loc='center right', fontsize=16)
    _save(fig, 'F3_distribution')


def _seed_json(cfg, seed):
    return ABL_DIR / f'{cfg}_seed{seed}.json'


def _basin_seeds(abl, cfg='A'):
    """Return (good_seed_dict, bad_seed_dict) for same config."""
    seeds = sorted(abl[cfg], key=lambda s: s['undulator_rms'])
    return seeds[0], seeds[-1]


def f4_perstage(abl):
    """F4: per-stage MSE divergence, good vs trapped seed (same config)."""
    print("F4: per-stage MSE divergence...")
    good, bad = _basin_seeds(abl, 'A')
    stages = np.arange(1, 12)
    g = [st['final_mse'] for st in good['stage_traces']]
    b = [st['final_mse'] for st in bad['stage_traces']]

    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    ax.plot(stages, g, '-o', color=C_BLUE, ms=11,
            label=f"converged (seed {good['seed']}, "
                  f"RMS {good['undulator_rms']:.1e})")
    ax.plot(stages, b, '-s', color=C_RED, ms=11,
            label=f"trapped (seed {bad['seed']}, "
                  f"RMS {bad['undulator_rms']:.1e})")
    ax.axvspan(10.5, 11.5, color=C_ORANGE, alpha=0.18)
    ax.annotate('Stage 11\nundulator match\n'
                f'{b[-1] / g[-1]:.0f}× worse',
                xy=(11, b[-1]), xytext=(8.0, b[-1] * 0.6),
                fontsize=17, ha='center',
                arrowprops=dict(arrowstyle='->', color='k', lw=1.5))

    ax.set_yscale('log')
    ax.set_xticks(stages)
    ax.set_xticklabels([f'{i}. {STAGE_LABELS[i-1]}' for i in stages],
                       rotation=40, ha='right', fontsize=14)
    ax.set_xlabel('Sequential optimisation stage')
    ax.set_ylabel('Stage final MSE')
    ax.set_title('The sequential chain breaks at the final undulator match\n'
                 '(Config A, same lattice, two beam seeds)', fontsize=20)
    ax.legend(loc='lower left', fontsize=16)
    _save(fig, 'F4_perstage')


def _ellipse_xy(beta, alpha, emit, n=200):
    t = np.linspace(0, 2 * np.pi, n)
    # Parametric phase-space ellipse for given Twiss + geometric emittance.
    x = np.sqrt(emit * beta) * np.cos(t)
    xp = -np.sqrt(emit / beta) * (np.sin(t) + alpha * np.cos(t))
    return x, xp


def f5_ellipse(params, abl):
    """F5: vertical Twiss ellipse at the undulator — cost of trapping.

    The bad basin keeps beta_x near target but collapses beta_y, so the
    physical cost is a vertical-plane mismatch. Uses the exact converged
    Twiss recorded in the ablation JSON, at a common geometric emittance.
    """
    print("F5: undulator phase-space ellipse (y-plane)...")
    good, bad = _basin_seeds(abl, 'A')
    emit = (8e-6) / params['norm']        # eps_n = 8 mm.mrad -> geometric

    tgt, fg, fb = good['targets'], good['final_twiss'], bad['final_twiss']
    yt, ypt = _ellipse_xy(tgt['beta_y'], tgt['alpha_y'], emit)
    yg, ypg = _ellipse_xy(fg['beta_y'], fg['alpha_y'], emit)
    yb, ypb = _ellipse_xy(fb['beta_y'], fb['alpha_y'], emit)

    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    ax.plot(yt * 1e3, ypt * 1e3, color='k', lw=3, ls='--',
            label=f"design target ($\\beta_y$={tgt['beta_y']:.3f} m)")
    ax.plot(yg * 1e3, ypg * 1e3, color=C_BLUE, lw=3,
            label=f"converged ($\\beta_y$={fg['beta_y']:.3f} m, "
                  f"RMS {good['undulator_rms']:.0e})")
    ax.plot(yb * 1e3, ypb * 1e3, color=C_RED, lw=3,
            label=f"trapped ($\\beta_y$={fb['beta_y']:.3f} m, "
                  f"RMS {bad['undulator_rms']:.0e})")
    ax.set_xlabel('y (mm)')
    ax.set_ylabel("y' (mrad)")
    ax.set_title('A trapped seed collapses the vertical match\n'
                 '(phase-space ellipse at the undulator entrance)',
                 fontsize=20)
    ax.legend(loc='upper right', fontsize=15)
    ax.axhline(0, color=C_GREY, lw=1)
    ax.axvline(0, color=C_GREY, lw=1)
    ax.set_aspect('auto')
    _save(fig, 'F5_ellipse')


def f6_gd_vs_bo(abl, bo_results):
    """F6: GD vs BO. NM curves are real; BO is a gated template."""
    print("F6: GD vs BO outlook...")
    fig, ax = plt.subplots(figsize=(12.8, 7.2))

    width = 0.35
    xs = np.arange(len(CONFIGS))
    nm_med = []
    for cfg in CONFIGS:
        rms = np.array([s['undulator_rms'] for s in abl[cfg]])
        nm_med.append(np.median(rms))
    ax.bar(xs - width / 2, nm_med, width, color=[CONFIG_COLOR[c]
           for c in CONFIGS], edgecolor='k', label='Nelder–Mead (this work)')

    if bo_results and Path(bo_results).exists():
        with open(bo_results) as fh:
            bo = json.load(fh)
        bo_med = [bo[c]['median_rms'] for c in CONFIGS]
        ax.bar(xs + width / 2, bo_med, width, color=C_GREY,
               edgecolor='k', hatch='//', label='Bayesian optimisation')
    else:
        for x in xs:
            ax.bar(x + width / 2, max(nm_med), width, color='none',
                   edgecolor=C_GREY, hatch='..', linewidth=1.5)
        ax.text(np.mean(xs), max(nm_med) * 0.30,
                'BO — work in progress\n(git-bug 7f690aa S6)',
                fontsize=17, color='dimgray', ha='center', style='italic')

    ax.set_yscale('log')
    ax.set_xticks(xs)
    ax.set_xticklabels([CONFIG_DESC[c] for c in CONFIGS], fontsize=15)
    ax.set_ylabel('Median final undulator RMS')
    ax.set_title('Outlook: can Bayesian optimisation escape the bad basin?',
                 fontsize=21)
    ax.legend(loc='upper left', fontsize=16)
    _save(fig, 'F6_gd_vs_bo')


# ── Part I ────────────────────────────────────────────────────────────────

def p1_capability_matrix():
    """P1: multi-code capability matrix — the modelling stack at a glance."""
    print("P1: capability matrix...")
    blocks = ['Injector\n(α-magnet)', 'S-band linac', 'Diagnostic\nchicane',
              'Transport line', 'Space charge']
    codes = ['FELsim\n(1st order)', 'COSY\nINFINITY', 'RF-Track', 'elegant']
    # 0 = n/a, 1 = implemented, 2 = cross-validated
    grid = np.array([
        [1, 2, 1, 0],   # injector
        [1, 1, 2, 2],   # linac (RF-Track validated vs elegant, 0.06%)
        [1, 2, 1, 0],   # chicane
        [2, 2, 2, 0],   # transport line
        [0, 2, 2, 0],   # space charge (COSY SC + DA-FMM/Xsuite)
    ])
    cmap = {0: '#EEEEEE', 1: C_TEAL, 2: C_BLUE}
    txt = {0: '–', 1: 'model', 2: 'validated'}

    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    for i in range(len(blocks)):
        for j in range(len(codes)):
            v = grid[i, j]
            ax.add_patch(plt.Rectangle((j, len(blocks) - 1 - i), 1, 1,
                         facecolor=cmap[v], edgecolor='white', lw=3))
            ax.text(j + 0.5, len(blocks) - 1 - i + 0.5, txt[v],
                    ha='center', va='center', fontsize=16,
                    color='white' if v else 'gray',
                    fontweight='bold' if v == 2 else 'normal')
    ax.set_xlim(0, len(codes))
    ax.set_ylim(0, len(blocks))
    ax.set_xticks(np.arange(len(codes)) + 0.5)
    ax.set_xticklabels(codes, fontsize=16)
    ax.set_yticks(np.arange(len(blocks)) + 0.5)
    ax.set_yticklabels(blocks[::-1], fontsize=16)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title('FELsim multi-code modelling stack', fontsize=22)
    ax.legend(handles=[Patch(color=C_BLUE, label='cross-validated'),
                       Patch(color=C_TEAL, label='model in place'),
                       Patch(color='#EEEEEE', label='not applicable')],
              loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=15)
    _save(fig, 'P1_capability_matrix')


def p2_linac_vs_elegant():
    """P2: RF-Track linac model vs elegant phase scan."""
    print("P2: RF-Track vs elegant...")
    rft = _read_csv(THIS_DIR / 'rftrack_linac' / 'rftrack_linac_phase_scan.csv')
    ele = _read_csv(THIS_DIR / 'elegant_linac' / 'phase_scan_results.csv')

    # RF-Track phase convention is offset by +70 deg from elegant's RFCA.
    phi = (rft['phase_deg'] + 70.0) % 360.0
    order = np.argsort(phi)
    phi, k_rft = phi[order], rft['K_out_MeV'][order]
    # Restrict to the accelerating lobe through crest (elegant's raw scan
    # has non-physical wrap spikes in the far decelerating region).
    win = (phi >= 0) & (phi <= 180)
    ew = (ele['phase_deg'] >= 0) & (ele['phase_deg'] <= 180)

    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    ax.plot(phi[win], k_rft[win], '-', color=C_BLUE,
            label='RF-Track (TW_Structure, +70° aligned)')
    ax.plot(ele['phase_deg'][ew], ele['K_out_MeV'][ew], '--', color=C_RED,
            label='elegant (RFCA + TWLA)')

    pk_rft = float(np.nanmax(rft['K_out_MeV']))
    ie = int(np.argmax(ele['K_out_MeV']))
    pk_ele, pk_phase = ele['K_out_MeV'][ie], ele['phase_deg'][ie]
    rel = abs(pk_rft - pk_ele) / pk_ele * 100
    ax.scatter([pk_phase], [pk_rft], color='k', zorder=5, s=140)
    ax.annotate(f'peak agreement {rel:.2f}%\n'
                f'({pk_rft:.3f} vs {pk_ele:.3f} MeV)',
                xy=(pk_phase, pk_rft), xytext=(pk_phase + 18, pk_rft * 0.55),
                fontsize=18, arrowprops=dict(arrowstyle='->', lw=1.5))

    ax.set_xlabel('RF phase, elegant convention (deg) — accelerating lobe')
    ax.set_ylabel('Output kinetic energy (MeV)')
    ax.set_title('Linac model validated against an independent tracker\n'
                 '(accelerating lobe through crest; off-crest region '
                 'omitted — code phase conventions diverge there)',
                 fontsize=17)
    ax.legend(loc='lower center', fontsize=17)
    _save(fig, 'P2_linac_vs_elegant')


def p3_space_charge():
    """P3: space-charge engine comparison (emittance growth vs charge)."""
    print("P3: space-charge capability...")
    sc_dir = THIS_DIR / 'rftrack_linac' / 'sc_compare_output'
    engines = [('sweep_dafmm.csv', 'DA-FMM (COSY)', C_BLUE, '-o'),
               ('sweep_xsuite-frozen.csv', 'Xsuite frozen', C_TEAL, '-s'),
               ('sweep_xsuite-pic3d.csv', 'Xsuite PIC3D', C_RED, '-^')]

    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    for fname, label, color, style in engines:
        p = sc_dir / fname
        if not p.exists():
            continue
        d = _read_csv(p)
        ax.plot(d['Q_nC'], d['epsnx_growth'] * 100, style, color=color,
                ms=10, label=label)

    ax.set_xlabel('Bunch charge (nC)')
    ax.set_ylabel(r'Norm. emittance growth $\Delta\varepsilon_{nx}$ (%)')
    ax.set_title('Three space-charge models now in the FELsim stack',
                 fontsize=21)
    ax.legend(loc='upper left', fontsize=17)
    _save(fig, 'P3_space_charge')


def p4_cross_code():
    """P4: cross-code transport-line Twiss match quality."""
    print("P4: cross-code Twiss match...")
    with open(RESULTS_DIR / 'R2' / 'R2_summary.json') as fh:
        rows = json.load(fh)['table_1_baseline']

    codes = [r['Code'] for r in rows]
    mse = [r['MSE'] for r in rows]
    by = [r['β_y'] for r in rows]

    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    bars = ax.bar(codes, mse, color=[C_BLUE, C_TEAL, C_RED],
                  edgecolor='k', width=0.55)
    ax.set_yscale('log')
    for b, m, r in zip(bars, mse, rows):
        ax.text(b.get_x() + b.get_width() / 2, m * 1.4,
                f'MSE {m:.1e}\n' + r['Model'], ha='center', fontsize=14)
    ax.set_ylabel('Twiss-match MSE at the undulator')
    ax.set_title('Cross-code Twiss match at the undulator\n'
                 r'(COSY DA maps to $10^{-9}$; RF-Track fringe-limited '
                 r'in $\beta_y$; $\varepsilon_n=8$ baseline)', fontsize=18)
    ax.annotate(f'RF-Track $\\beta_y$ deficit\n({by[-1]:.3f} vs '
                f'{by[0]:.3f} m — missing fringe $\\varphi$)',
                xy=(2, mse[2]), xytext=(0.6, mse[2] * 0.12),
                fontsize=15, arrowprops=dict(arrowstyle='->', lw=1.4))
    ax.set_ylim(min(mse) * 0.05, max(mse) * 6)
    _save(fig, 'P4_cross_code')


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--beta-xm', type=float, default=1.267,
                        help='TUPM005 target beta_x in m (default 1.267)')
    parser.add_argument('--alpha-xm', type=float, default=0.560,
                        help='TUPM005 target alpha_x (default 0.560)')
    parser.add_argument('--bo-results', type=str, default=None,
                        help='JSON of BO results (S6); placeholder if absent')
    parser.add_argument('--only', nargs='*', default=None,
                        help='Subset of figure ids, e.g. F2 P1')
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.style.use(str(STYLE_PATH))
    params = compute_params(beta_xm=args.beta_xm, alpha_xm=args.alpha_xm)
    abl = load_ablation()

    print(f"IPAC TUPM005 figures — targets beta_x={args.beta_xm}, "
          f"alpha_x={args.alpha_xm}, beta_y={params['beta_ym']:.4f}")
    print(f"Output: {OUTPUT_DIR}\n")

    figures = {
        'P1': lambda: p1_capability_matrix(),
        'P2': lambda: p2_linac_vs_elegant(),
        'P3': lambda: p3_space_charge(),
        'P4': lambda: p4_cross_code(),
        'F1': lambda: f1_twiss(params),
        'F2': lambda: f2_failrate(abl),
        'F3': lambda: f3_distribution(abl),
        'F4': lambda: f4_perstage(abl),
        'F5': lambda: f5_ellipse(params, abl),
        'F6': lambda: f6_gd_vs_bo(abl, args.bo_results),
    }
    todo = args.only if args.only else list(figures)
    for key in todo:
        figures[key]()


if __name__ == '__main__':
    main()
