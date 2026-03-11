"""W10: Beam Losses & Bunch Compression Study

Quantifies particle losses through the UH MkV FEL transport line for both
0.5 ps and 2 ps operating modes. Demonstrates that bunch compression requires
a negative chirp (not just increased energy spread), and explores higher-charge
operation.

Part A: Transmission baseline (COSY + RF-Track, with/without apertures)
Part B: Compression via chirp — 6 scenarios including no-compression controls
Part C: Charge scan at both operating modes

Author: Eremey Valetov
"""

import sys
import json
import math
import argparse
import time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ebeam import beam
from beamline import lattice
from excelElements import ExcelElements
from beamOptimizer import beamOptimizer

from UHM_beamline_opt_05ps_params import (
    run_optimization, compute_twiss_targets, QUAD_INDICES,
    ENERGY, RF_FREQ, SEGMENTS,
)

# ── Constants ────────────────────────────────────────────────────────────────
C_LIGHT = 299792458.0
M_E_MEV = 0.51099895
GAMMA_REL = 1 + ENERGY / M_E_MEV
BETA_REL = np.sqrt(1 - 1 / GAMMA_REL**2)
P_C = GAMMA_REL * BETA_REL * M_E_MEV

# Apertures [m]
QUAD_HALF_APERTURE = 0.0135      # 27 mm bore / 2
DIPOLE_HALF_GAP = 0.00724        # 14.48 mm gap / 2
DIPOLE_HALF_WIDTH = 0.025        # 50 mm placeholder / 2
BEAM_PIPE_RADIUS = 0.0127        # 1" beam pipe

OUTDIR = Path(__file__).resolve().parent / 'results' / 'W10'
EXCEL_PATH = (Path(__file__).resolve().parent.parent.parent
              / 'beam_excel' / 'Beamline_elements.xlsx')

N_PARTICLES = 5000
SEED = 42


def _print(msg):
    print(msg, flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
#  Beam generation utilities
# ═══════════════════════════════════════════════════════════════════════════════

def create_felsim_beam(sigma_t_ps=2.0, sigma_E_pct=0.5, h_chirp=0.0,
                       epsilon_n=8, x_std=0.8, y_std=0.8,
                       N=N_PARTICLES, seed=SEED):
    """Generate 6D Gaussian beam in FELsim coordinates.

    FELsim: [x(mm), x'(mrad), y(mm), y'(mrad), ΔToF/T×10³, δW/W×10³]
    """
    _, _, _, _, relat = compute_twiss_targets()
    norm = relat.gamma * relat.beta
    epsilon = epsilon_n / norm

    np.random.seed(seed)
    tof_std = sigma_t_ps * 1e-9 * RF_FREQ
    energy_std = sigma_E_pct * 10  # FELsim units: δW/W × 10³

    beam_dist = beam().gen_6d_gaussian(
        0, [x_std, epsilon / x_std, y_std, epsilon / y_std, tof_std, energy_std],
        N)

    if h_chirp != 0:
        tof_dist = beam_dist[:, 4] / RF_FREQ  # convert to seconds
        beam_dist[:, 5] += h_chirp * tof_dist  # FELsim δW/W×10³ += h × Δt × 10³

    return beam_dist


def create_cosy_beam(sigma_t_ps=2.0, sigma_E_pct=0.5, h_chirp=0.0,
                     epsilon_n=8, x_std_mm=0.8, N=N_PARTICLES, seed=SEED):
    """Generate 6D Gaussian beam in COSY coordinates.

    COSY: [x(m), a(rad), y(m), b(rad), l(m), δ=ΔK/K₀]
    """
    _, _, _, _, relat = compute_twiss_targets()
    norm = relat.gamma * relat.beta
    epsilon_geom = epsilon_n / norm  # pi.mm.mrad

    sigma_x = x_std_mm * 1e-3  # m
    sigma_xp = epsilon_geom * 1e-6 / sigma_x  # rad
    sigma_z = sigma_t_ps * 1e-12 * BETA_REL * C_LIGHT  # m
    sigma_delta = sigma_E_pct / 100  # fractional ΔK/K₀

    rng = np.random.default_rng(seed)
    X = np.zeros((N, 6))
    X[:, 0] = rng.normal(0, sigma_x, N)
    X[:, 1] = rng.normal(0, sigma_xp, N)
    X[:, 2] = rng.normal(0, sigma_x, N)
    X[:, 3] = rng.normal(0, sigma_xp, N)
    X[:, 4] = rng.normal(0, sigma_z, N)
    X[:, 5] = rng.normal(0, sigma_delta, N)

    if h_chirp != 0:
        X[:, 5] += h_chirp * X[:, 4] / (BETA_REL * C_LIGHT)

    return X


def beam_sigma_z_ps(particles_felsim):
    """Return RMS bunch length in ps from FELsim-coordinate particles."""
    return np.std(particles_felsim[:, 4]) / RF_FREQ * 1e9


def peak_current(Q_C, sigma_t_ps, transmission=1.0):
    """Peak current = Q × T / (√(2π) × σ_t)."""
    sigma_t = sigma_t_ps * 1e-12
    if sigma_t <= 0:
        return 0.0
    return Q_C * transmission / (np.sqrt(2 * np.pi) * sigma_t)


# ═══════════════════════════════════════════════════════════════════════════════
#  FELsim optimization → quad currents
# ═══════════════════════════════════════════════════════════════════════════════

def get_optimized_currents():
    """Run FELsim 11-stage optimization and return quad currents dict."""
    _print("Running FELsim 11-stage optimization for quad currents...")
    t0 = time.perf_counter()
    res = run_optimization(
        bunch_spread=2.0, energy_std_percent=0.5, h=0,
        epsilon_n=8, nb_particles=500, seed=42,
    )
    elapsed = time.perf_counter() - t0
    _print(f"  RMS = {math.sqrt(res['mse']):.4e}, time = {elapsed:.1f} s")
    _print(f"  β_x={res['beta_x']:.4f}, α_x={res['alpha_x']:.4f}, "
           f"β_y={res['beta_y']:.4f}, α_y={res['alpha_y']:.4f}")
    return res['quad_currents']


def setup_felsim_line(currents):
    """Set up FELsim beamline with given quad currents."""
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    excel = ExcelElements(EXCEL_PATH)
    beamlineUH = excel.create_beamline()
    line = relat.changeBeamType("electron", ENERGY, beamlineUH)[:SEGMENTS]
    for idx, current in currents.items():
        if idx < len(line):
            line[idx].current = current
    return line


# ═══════════════════════════════════════════════════════════════════════════════
#  COSY tracking with apertures
# ═══════════════════════════════════════════════════════════════════════════════

def _apply_felsim_currents(beamline, currents):
    """Apply FELsim-indexed quad currents to a FELsim-indexed COSY beamline."""
    for idx, current in currents.items():
        if idx < len(beamline):
            beamline[idx]['current'] = current
    return beamline


def run_cosy_tracking(currents, beam_cosy, label=""):
    """Run COSY 3D particle tracking with aperture cuts.

    Uses parse_beamline_felsim_indexed() so FELsim quad indices map correctly
    to the COSY beamline (ExcelElements and BeamlineBuilder use different
    element indexing).

    Returns dict with transmission data and beam stats.
    """
    from cosyAdapter import COSYAdapter
    from cosyOptHelper import parse_beamline_felsim_indexed

    sim = COSYAdapter(
        excel_path=str(EXCEL_PATH),
        mode='particle_tracking',
        fringe_field_order=0,
        quad_aperture=QUAD_HALF_APERTURE * 2,
        dipole_aperture=DIPOLE_HALF_GAP * 2,
        config={'simulation': {'dimensions': 3, 'KE': ENERGY}},
    )
    sim.enable_aperture_cuts(dipole_half_width=DIPOLE_HALF_WIDTH)

    native = sim.get_native_simulator()
    native.use_enge_coeffs = False
    beamline = parse_beamline_felsim_indexed(str(EXCEL_PATH))
    _apply_felsim_currents(beamline, currents)
    native.beamline = beamline[:SEGMENTS]

    particles_felsim = native.transform_from_cosy_coordinates(beam_cosy)

    _print(f"  [{label}] Tracking {beam_cosy.shape[0]} particles through COSY...")
    evolution = sim.collect_evolution(particles_felsim, checkpoint_elements='all')

    # Extract final state
    n_initial = beam_cosy.shape[0]
    if evolution.s_positions:
        s_final = max(evolution.s_positions)
        ps_final = evolution.particles.get(s_final)
        n_final = ps_final.shape[0] if ps_final is not None else 0
    else:
        n_final = 0
        ps_final = None

    transmission = n_final / n_initial if n_initial > 0 else 0

    elem_transmission = []
    for ei in evolution.elements:
        elem_transmission.append({
            'index': ei.index,
            'type': ei.element_type,
            's_end': ei.s_end,
            'n_good': ei.parameters.get('n_good', n_initial),
            'transmission': ei.parameters.get('transmission', 1.0),
        })

    # Compute bunch length at undulator from final particles (FELsim col 4)
    sigma_z_ps = None
    twiss_final = {}
    if ps_final is not None and n_final > 1:
        tof_std = np.std(ps_final[:, 4])
        sigma_z_ps = tof_std / RF_FREQ * 1e9  # ps
        twiss_final = evolution.twiss.get(s_final, {})

    return {
        'label': label,
        'n_initial': n_initial,
        'n_final': n_final,
        'transmission': transmission,
        'sigma_z_ps': sigma_z_ps,
        'twiss': twiss_final,
        'element_transmission': elem_transmission,
        'evolution': evolution,
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  RF-Track tracking with apertures
# ═══════════════════════════════════════════════════════════════════════════════

def run_rftrack_tracking(currents, beam_felsim, label="",
                         space_charge=False, Q_C=60e-12):
    """Run RF-Track particle tracking with physical apertures."""
    try:
        from rftrackAdapter import RFTrackAdapter
    except ImportError:
        _print(f"  [{label}] RF-Track not available, skipping")
        return None

    sim = RFTrackAdapter(
        lattice_path=str(EXCEL_PATH),
        beam_energy=ENERGY,
        space_charge=space_charge,
    )
    sim.beamline = sim.beamline[:SEGMENTS]
    for idx, current in currents.items():
        if idx < len(sim.beamline):
            sim._modify_element(idx, current=current)

    sim.enable_physical_apertures()

    if space_charge:
        sim.set_space_charge(True, mesh=(32, 32, 64))

    sim._build_lattice()

    _print(f"  [{label}] Tracking {beam_felsim.shape[0]} particles through RF-Track"
           f"{' (SC)' if space_charge else ''}...")

    evolution = sim.collect_evolution(beam_felsim, checkpoint_elements='all')

    if evolution.s_positions:
        s_final = max(evolution.s_positions)
        ps_final = evolution.particles.get(s_final)
        n_final = ps_final.shape[0] if ps_final is not None else 0
    else:
        n_final = 0
        ps_final = None

    n_initial = beam_felsim.shape[0]
    transmission = n_final / n_initial if n_initial > 0 else 0

    elem_transmission = []
    for ei in evolution.elements:
        elem_transmission.append({
            'index': ei.index,
            'type': ei.element_type,
            's_end': ei.s_end,
            'n_good': ei.parameters.get('n_good', n_initial),
            'n_lost': ei.parameters.get('n_lost', 0),
            'transmission': ei.parameters.get('transmission', 1.0),
        })

    sigma_z_ps = None
    twiss_final = {}
    if ps_final is not None and n_final > 1:
        # σ_z from FELsim coords: col 4 is ΔToF/T×10³
        tof_std = np.std(ps_final[:, 4])
        sigma_z_ps = tof_std / RF_FREQ * 1e9  # ps
        twiss_final = evolution.twiss.get(s_final, {})

    return {
        'label': label,
        'n_initial': n_initial,
        'n_final': n_final,
        'transmission': transmission,
        'sigma_z_ps': sigma_z_ps,
        'twiss': twiss_final,
        'element_transmission': elem_transmission,
        'evolution': evolution,
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  Plotting utilities
# ═══════════════════════════════════════════════════════════════════════════════

def plot_transmission_profile(results_list, title, filename):
    """Plot element-by-element transmission for multiple scenarios."""
    fig, ax = plt.subplots(figsize=(12, 5))
    for res in results_list:
        if res is None:
            continue
        et = res['element_transmission']
        s_vals = [e['s_end'] for e in et]
        t_vals = [e['transmission'] * 100 for e in et]
        ax.plot(s_vals, t_vals, label=res['label'], linewidth=1.5)

    ax.set_xlabel('s (m)')
    ax.set_ylabel('Transmission (%)')
    ax.set_title(title)
    ax.set_ylim(bottom=0, top=105)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'{filename}.{ext}', dpi=150)
    plt.close(fig)
    _print(f"  Saved: {filename}.{{eps,png}}")


def plot_compression_summary(scenarios, filename):
    """Bar chart of compression ratio, transmission, and I_peak."""
    labels = [s['label'] for s in scenarios]
    compressions = [s.get('compression_ratio', 1.0) for s in scenarios]
    transmissions = [s['transmission'] * 100 for s in scenarios]
    I_peaks = [s.get('I_peak_A', 0) for s in scenarios]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

    x = np.arange(len(labels))
    colors = ['C0' if s.get('h_chirp', 0) != 0 else 'C1' for s in scenarios]

    ax1.bar(x, compressions, color=colors, alpha=0.7, edgecolor='k')
    ax1.axhline(1.0, color='r', ls='--', lw=1)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
    ax1.set_ylabel('Compression ratio (σ_z,in / σ_z,out)')
    ax1.set_title('Bunch Compression')
    ax1.grid(True, alpha=0.3, axis='y')

    ax2.bar(x, transmissions, color=colors, alpha=0.7, edgecolor='k')
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
    ax2.set_ylabel('Transmission (%)')
    ax2.set_title('Particle Transmission')
    ax2.set_ylim(0, 105)
    ax2.grid(True, alpha=0.3, axis='y')

    ax3.bar(x, I_peaks, color=colors, alpha=0.7, edgecolor='k')
    ax3.set_xticks(x)
    ax3.set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
    ax3.set_ylabel('I_peak (A)')
    ax3.set_title('Peak Current at Undulator')
    ax3.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'{filename}.{ext}', dpi=150)
    plt.close(fig)
    _print(f"  Saved: {filename}.{{eps,png}}")


def plot_charge_scan(results_2ps, results_05ps, filename):
    """Plot transmission and I_peak vs charge."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    for results, label, color in [(results_2ps, '2 ps (no chirp)', 'C0'),
                                   (results_05ps, '0.5 ps (C=4 chirp)', 'C1')]:
        charges = [r['Q_pC'] for r in results]
        trans = [r['transmission'] * 100 for r in results]
        I_peaks = [r.get('I_peak_A', 0) for r in results]
        ax1.plot(charges, trans, 'o-', label=label, color=color)
        ax2.plot(charges, I_peaks, 'o-', label=label, color=color)

    ax1.set_xlabel('Charge (pC)')
    ax1.set_ylabel('Transmission (%)')
    ax1.set_title('Transmission vs Charge')
    ax1.set_ylim(0, 105)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.set_xlabel('Charge (pC)')
    ax2.set_ylabel('I_peak (A)')
    ax2.set_title('Peak Current vs Charge')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'{filename}.{ext}', dpi=150)
    plt.close(fig)
    _print(f"  Saved: {filename}.{{eps,png}}")


# ═══════════════════════════════════════════════════════════════════════════════
#  Part A: Transmission Baseline
# ═══════════════════════════════════════════════════════════════════════════════

def part_a(currents):
    """Baseline transmission at 2 ps and 0.5 ps, COSY and RF-Track."""
    _print("\n" + "=" * 72)
    _print("  Part A: Transmission Baseline")
    _print("=" * 72)

    results = []

    scenarios = [
        ('A1: 2ps COSY',     2.0, 'cosy'),
        ('A2: 2ps RF-Track', 2.0, 'rftrack'),
        ('A3: 0.5ps COSY',   0.5, 'cosy'),
        ('A4: 0.5ps RF-Track', 0.5, 'rftrack'),
    ]

    for label, sigma_t_ps, code in scenarios:
        _print(f"\n── {label} ──")
        if code == 'cosy':
            beam_cosy = create_cosy_beam(sigma_t_ps=sigma_t_ps, sigma_E_pct=0.5,
                                         h_chirp=0, N=N_PARTICLES, seed=SEED)
            res = run_cosy_tracking(currents, beam_cosy, label=label)
        else:
            beam_felsim = create_felsim_beam(sigma_t_ps=sigma_t_ps, sigma_E_pct=0.5,
                                             h_chirp=0, N=N_PARTICLES, seed=SEED)
            res = run_rftrack_tracking(currents, beam_felsim, label=label)

        if res is not None:
            Q = 60e-12
            res['sigma_t_ps_initial'] = sigma_t_ps
            res['I_peak_A'] = peak_current(Q, res['sigma_z_ps'] or sigma_t_ps,
                                           res['transmission'])
            results.append(res)
            _print(f"  Transmission: {res['transmission']:.4f} "
                   f"({res['n_final']}/{res['n_initial']})")
            if res['sigma_z_ps'] is not None:
                _print(f"  σ_z at undulator: {res['sigma_z_ps']:.3f} ps")
            _print(f"  I_peak: {res['I_peak_A']:.1f} A")

    # Summary table
    _print("\n── Part A Summary ──")
    _print(f"{'Scenario':<22s}  {'T (%)':>8s}  {'N_final':>8s}  "
           f"{'σ_z (ps)':>10s}  {'I_peak (A)':>10s}")
    _print("-" * 64)
    for r in results:
        sz = f"{r['sigma_z_ps']:.3f}" if r['sigma_z_ps'] else "—"
        _print(f"{r['label']:<22s}  {r['transmission']*100:8.2f}  {r['n_final']:8d}  "
               f"{sz:>10s}  {r['I_peak_A']:10.1f}")

    # Identify limiting apertures
    for r in results:
        lost_elems = [(e['index'], e['type'], e['s_end'])
                      for e in r['element_transmission']
                      if e.get('n_lost', 0) > 0 or
                      (e.get('n_good', N_PARTICLES) < N_PARTICLES and
                       r['element_transmission'].index(e) > 0 and
                       e['n_good'] < r['element_transmission'][r['element_transmission'].index(e)-1].get('n_good', N_PARTICLES))]
        if lost_elems:
            _print(f"\n  Limiting apertures for {r['label']}:")
            for idx, etype, s_end in lost_elems[:5]:
                _print(f"    Element {idx} ({etype}) at s={s_end:.3f} m")

    # Transmission profile plot
    OUTDIR.mkdir(parents=True, exist_ok=True)
    plot_transmission_profile(results, 'Part A: Transmission Baseline',
                              'W10_partA_transmission')

    # Save data
    save_data = []
    for r in results:
        save_data.append({
            'label': r['label'],
            'transmission': r['transmission'],
            'n_initial': r['n_initial'],
            'n_final': r['n_final'],
            'sigma_z_ps': r['sigma_z_ps'],
            'I_peak_A': r['I_peak_A'],
            'sigma_t_ps_initial': r['sigma_t_ps_initial'],
        })
    with open(OUTDIR / 'part_a_results.json', 'w') as f:
        json.dump(save_data, f, indent=2)
    _print(f"\n  Saved: part_a_results.json")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  Part B: Compression via Chirp
# ═══════════════════════════════════════════════════════════════════════════════

def part_b(currents):
    """Bunch compression study: chirp vs energy spread increase."""
    _print("\n" + "=" * 72)
    _print("  Part B: Compression via Chirp")
    _print("=" * 72)

    Q = 60e-12
    sigma_t_initial = 2.0  # ps — beam always starts at 2 ps

    # Scenarios: (label, σ_E₀ [%], h [1/s], description)
    scenarios = [
        ('B1: baseline',     0.5, 0.0,       'No chirp, baseline'),
        ('B2: C=2 chirp',    0.5, -4.2e9,    'Moderate chirp, C≈2'),
        ('B3: C=4 chirp',    0.5, -8.3e9,    'Full chirp, C≈4'),
        ('B4: C=6 chirp',    0.5, -12.5e9,   'Strong chirp, C≈6'),
        ('B5: σ_E=2% only',  2.0, 0.0,       'Energy spread increase, no chirp'),
        ('B6: σ_E=3% only',  3.0, 0.0,       'Large energy spread, no chirp'),
    ]

    results = []

    # Scenarios that also get RF-Track validation
    rftrack_scenarios = {'B1', 'B3'}

    for label, sigma_E_pct, h_chirp, desc in scenarios:
        _print(f"\n── {label}: {desc} ──")
        _print(f"  σ_E₀={sigma_E_pct}%, h={h_chirp:.1e}/s")

        # Effective energy spread after chirp
        sigma_t_s = sigma_t_initial * 1e-12
        sigma_delta_0 = sigma_E_pct / 100
        if h_chirp != 0:
            sigma_delta_chirp = abs(h_chirp) * sigma_t_s
            sigma_delta_eff = np.sqrt(sigma_delta_0**2 + sigma_delta_chirp**2)
        else:
            sigma_delta_eff = sigma_delta_0
        _print(f"  σ_δ,eff = {sigma_delta_eff*100:.3f}%")

        # COSY tracking with apertures
        beam_cosy = create_cosy_beam(sigma_t_ps=sigma_t_initial,
                                     sigma_E_pct=sigma_E_pct,
                                     h_chirp=h_chirp,
                                     N=N_PARTICLES, seed=SEED)
        res = run_cosy_tracking(currents, beam_cosy, label=label)

        if res is None:
            _print(f"  COSY tracking failed for {label}")
            continue

        def _finish_result(res, sigma_t_initial, sigma_E_pct, h_chirp,
                           sigma_delta_eff, Q):
            if res['sigma_z_ps'] is not None and res['sigma_z_ps'] > 0:
                compression = sigma_t_initial / res['sigma_z_ps']
            else:
                compression = 0
            res['sigma_t_initial'] = sigma_t_initial
            res['sigma_E_pct'] = sigma_E_pct
            res['h_chirp'] = h_chirp
            res['sigma_delta_eff'] = sigma_delta_eff
            res['compression_ratio'] = compression
            res['I_peak_A'] = peak_current(Q, res['sigma_z_ps'] or sigma_t_initial,
                                           res['transmission'])
            return res

        res = _finish_result(res, sigma_t_initial, sigma_E_pct, h_chirp,
                             sigma_delta_eff, Q)
        results.append(res)

        compression = res['compression_ratio']
        _print(f"  Transmission: {res['transmission']:.4f}")
        if res['sigma_z_ps'] is not None:
            _print(f"  σ_z at undulator: {res['sigma_z_ps']:.3f} ps "
                   f"(C = {compression:.2f})")
        _print(f"  I_peak: {res['I_peak_A']:.1f} A")

        # RF-Track validation for selected scenarios
        scenario_tag = label.split(':')[0]
        if scenario_tag in rftrack_scenarios:
            rft_label = label.replace(':', ' RFT:')
            beam_felsim = create_felsim_beam(
                sigma_t_ps=sigma_t_initial, sigma_E_pct=sigma_E_pct,
                h_chirp=h_chirp, N=N_PARTICLES, seed=SEED)
            rft_res = run_rftrack_tracking(currents, beam_felsim, label=rft_label)
            if rft_res is not None:
                rft_res = _finish_result(rft_res, sigma_t_initial, sigma_E_pct,
                                         h_chirp, sigma_delta_eff, Q)
                results.append(rft_res)
                rft_c = rft_res['compression_ratio']
                _print(f"  [RF-Track] T={rft_res['transmission']:.4f}, "
                       f"σ_z={rft_res['sigma_z_ps']:.3f} ps (C={rft_c:.2f}), "
                       f"I_pk={rft_res['I_peak_A']:.1f} A")

    # Summary table
    _print("\n── Part B Summary ──")
    _print(f"{'Scenario':<20s}  {'σ_E₀':>5s}  {'h':>10s}  {'σ_δ,eff':>7s}  "
           f"{'T (%)':>6s}  {'σ_z,out':>8s}  {'C':>5s}  {'I_pk':>6s}")
    _print("-" * 82)
    for r in results:
        sz = f"{r['sigma_z_ps']:.3f}" if r['sigma_z_ps'] else "—"
        C_str = f"{r['compression_ratio']:.2f}" if r['compression_ratio'] > 0 else "—"
        _print(f"{r['label']:<20s}  {r['sigma_E_pct']:5.1f}  {r['h_chirp']:10.1e}  "
               f"{r['sigma_delta_eff']*100:6.3f}%  {r['transmission']*100:6.2f}  "
               f"{sz:>8s}  {C_str:>5s}  {r['I_peak_A']:6.1f}")

    # Key finding
    _print("\n── Key Finding ──")
    chirped = [r for r in results if r['h_chirp'] != 0]
    unchirped_high_E = [r for r in results if r['h_chirp'] == 0 and r['sigma_E_pct'] > 1]

    if chirped and unchirped_high_E:
        best_chirp = min(chirped, key=lambda r: r['sigma_z_ps'] or 999)
        worst_spread = max(unchirped_high_E, key=lambda r: r['sigma_z_ps'] or 0)
        _print(f"  Best chirped compression: {best_chirp['label']} → "
               f"σ_z = {best_chirp['sigma_z_ps']:.3f} ps (C={best_chirp['compression_ratio']:.2f})")
        if worst_spread['sigma_z_ps'] is not None:
            _print(f"  Energy spread only: {worst_spread['label']} → "
                   f"σ_z = {worst_spread['sigma_z_ps']:.3f} ps (ELONGATED)")
        _print("  → Chirp is necessary for compression; energy spread alone elongates the bunch.")

    # Plots
    OUTDIR.mkdir(parents=True, exist_ok=True)
    plot_compression_summary(results, 'W10_partB_compression')
    plot_transmission_profile(results, 'Part B: Compression Scenarios',
                              'W10_partB_transmission')

    # Save data
    save_data = []
    for r in results:
        save_data.append({
            'label': r['label'],
            'sigma_E_pct': r['sigma_E_pct'],
            'h_chirp': r['h_chirp'],
            'sigma_delta_eff': r['sigma_delta_eff'],
            'transmission': r['transmission'],
            'n_final': r['n_final'],
            'sigma_z_ps': r['sigma_z_ps'],
            'compression_ratio': r['compression_ratio'],
            'I_peak_A': r['I_peak_A'],
        })
    with open(OUTDIR / 'part_b_results.json', 'w') as f:
        json.dump(save_data, f, indent=2)
    _print(f"\n  Saved: part_b_results.json")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  Part C: Charge Scan
# ═══════════════════════════════════════════════════════════════════════════════

def part_c(currents):
    """Charge scan at both operating modes."""
    _print("\n" + "=" * 72)
    _print("  Part C: Charge Scan")
    _print("=" * 72)

    charges_pC = [20, 40, 60, 100, 150, 200, 300]
    results_2ps = []
    results_05ps = []

    for Q_pC in charges_pC:
        Q_C = Q_pC * 1e-12
        _print(f"\n── Q = {Q_pC} pC ──")

        # Mode 1: 2 ps, no chirp
        _print(f"  Mode 1: 2 ps, σ_E=0.5%, h=0")
        beam_cosy = create_cosy_beam(sigma_t_ps=2.0, sigma_E_pct=0.5,
                                     h_chirp=0, N=N_PARTICLES, seed=SEED)
        use_sc = Q_pC >= 100
        if use_sc:
            # RF-Track with space charge for high charge
            beam_felsim = create_felsim_beam(sigma_t_ps=2.0, sigma_E_pct=0.5,
                                             h_chirp=0, N=N_PARTICLES, seed=SEED)
            res = run_rftrack_tracking(
                currents, beam_felsim,
                label=f'{Q_pC}pC 2ps SC',
                space_charge=True, Q_C=Q_C)
        else:
            res = run_cosy_tracking(currents, beam_cosy,
                                    label=f'{Q_pC}pC 2ps')

        if res is not None:
            res['Q_pC'] = Q_pC
            res['Q_C'] = Q_C
            res['sigma_t_initial'] = 2.0
            res['I_peak_A'] = peak_current(Q_C, res['sigma_z_ps'] or 2.0,
                                           res['transmission'])
            results_2ps.append(res)
            _print(f"    T={res['transmission']:.4f}, I_pk={res['I_peak_A']:.1f} A")

        # Mode 2: 2 ps initial → 0.5 ps via C=4 chirp
        _print(f"  Mode 2: 2 ps initial, h=-8.3e9 (C=4 → 0.5 ps)")
        beam_cosy_chirp = create_cosy_beam(sigma_t_ps=2.0, sigma_E_pct=0.5,
                                           h_chirp=-8.3e9,
                                           N=N_PARTICLES, seed=SEED)
        if use_sc:
            beam_felsim_chirp = create_felsim_beam(sigma_t_ps=2.0,
                                                    sigma_E_pct=0.5,
                                                    h_chirp=-8.3e9,
                                                    N=N_PARTICLES, seed=SEED)
            res = run_rftrack_tracking(
                currents, beam_felsim_chirp,
                label=f'{Q_pC}pC 0.5ps SC',
                space_charge=True, Q_C=Q_C)
        else:
            res = run_cosy_tracking(currents, beam_cosy_chirp,
                                    label=f'{Q_pC}pC 0.5ps')

        if res is not None:
            res['Q_pC'] = Q_pC
            res['Q_C'] = Q_C
            res['sigma_t_initial'] = 2.0
            res['h_chirp'] = -8.3e9
            res['I_peak_A'] = peak_current(Q_C, res['sigma_z_ps'] or 0.5,
                                           res['transmission'])
            results_05ps.append(res)
            _print(f"    T={res['transmission']:.4f}, I_pk={res['I_peak_A']:.1f} A")

    # Summary table
    _print("\n── Part C Summary: 2 ps mode ──")
    _print(f"{'Q (pC)':>8s}  {'T (%)':>8s}  {'σ_z (ps)':>10s}  {'I_peak (A)':>10s}")
    _print("-" * 42)
    for r in results_2ps:
        sz = f"{r['sigma_z_ps']:.3f}" if r['sigma_z_ps'] else "—"
        _print(f"{r['Q_pC']:8d}  {r['transmission']*100:8.2f}  {sz:>10s}  "
               f"{r['I_peak_A']:10.1f}")

    _print("\n── Part C Summary: 0.5 ps mode (C=4 chirp) ──")
    _print(f"{'Q (pC)':>8s}  {'T (%)':>8s}  {'σ_z (ps)':>10s}  {'I_peak (A)':>10s}")
    _print("-" * 42)
    for r in results_05ps:
        sz = f"{r['sigma_z_ps']:.3f}" if r['sigma_z_ps'] else "—"
        _print(f"{r['Q_pC']:8d}  {r['transmission']*100:8.2f}  {sz:>10s}  "
               f"{r['I_peak_A']:10.1f}")

    # Plots
    OUTDIR.mkdir(parents=True, exist_ok=True)
    if results_2ps and results_05ps:
        plot_charge_scan(results_2ps, results_05ps, 'W10_partC_charge_scan')

    # Save data
    def _serialize(r):
        return {k: v for k, v in r.items()
                if k not in ('evolution', 'element_transmission', 'twiss')}

    save_data = {
        '2ps': [_serialize(r) for r in results_2ps],
        '0.5ps': [_serialize(r) for r in results_05ps],
    }
    with open(OUTDIR / 'part_c_results.json', 'w') as f:
        json.dump(save_data, f, indent=2, default=str)
    _print(f"\n  Saved: part_c_results.json")

    return results_2ps, results_05ps


# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="W10: Beam Losses & Bunch Compression Study")
    parser.add_argument('--part-a', action='store_true',
                        help='Part A: transmission baseline')
    parser.add_argument('--part-b', action='store_true',
                        help='Part B: compression via chirp')
    parser.add_argument('--part-c', action='store_true',
                        help='Part C: charge scan')
    parser.add_argument('--all', action='store_true',
                        help='Run all parts')
    parser.add_argument('--currents', type=str, default=None,
                        help='Path to pre-computed quad currents JSON (skip optimization)')
    args = parser.parse_args()

    if not any([args.part_a, args.part_b, args.part_c, args.all]):
        args.all = True

    _print("W10: Beam Losses & Bunch Compression Study")
    _print(f"E = {ENERGY} MeV, γ = {GAMMA_REL:.2f}, β = {BETA_REL:.6f}")
    _print(f"p₀c = {P_C:.3f} MeV, f_RF = {RF_FREQ / 1e6:.0f} MHz")
    _print(f"N = {N_PARTICLES} particles, seed = {SEED}")

    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Get optimized quad currents (shared across all parts)
    if args.currents:
        with open(args.currents) as f:
            data = json.load(f)
        currents = {int(k): float(v) for k, v in data.items()}
        _print(f"Loaded pre-computed currents from {args.currents}")
    else:
        currents = get_optimized_currents()
        # Save currents for reuse
        currents_path = OUTDIR / 'currents_felsim.json'
        with open(currents_path, 'w') as f:
            json.dump({str(k): float(v) for k, v in currents.items()}, f, indent=2)
        _print(f"Saved currents to {currents_path}")

    if args.part_a or args.all:
        part_a(currents)

    if args.part_b or args.all:
        part_b(currents)

    if args.part_c or args.all:
        part_c(currents)

    _print("\n" + "=" * 72)
    _print("  W10 Complete")
    _print("=" * 72)


if __name__ == "__main__":
    main()
