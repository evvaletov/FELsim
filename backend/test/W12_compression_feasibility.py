"""W12: Bunch Compression Feasibility Study

Can the UH MkV FEL transport line deliver 0.5 ps bunches at the undulator
when starting from a 2 ps linac bunch?

Part A: Analytical compression + chirp sweep via COSY map propagation
Part B: RF-Track validation (post-C7 coord5 fix, re-run of W10 Part B)
Part C: Extended current bounds (15 A) — transverse MSE vs σ_z
Part D: Notes and feasibility summary

Key prior results:
  R56 = 27.09 mm (COSY, δ = ΔK/K₀), T566 = 0 (W9 Part A)
  R56 is geometry-locked — determined by dipole angles/spacings (W9 Part C)
  W10 RF-Track compression results invalidated by C7 coord5 bug

Author: Eremey Valetov
"""

import sys
import json
import argparse
import math
import time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from UHM_beamline_opt_05ps_params import (
    run_optimization, compute_twiss_targets, ENERGY, RF_FREQ, SEGMENTS,
)

# ── Constants ─────────────────────────────────────────────────────────────────
C_LIGHT = 299792458.0
M_E_MEV = 0.51099895
GAMMA = 1 + ENERGY / M_E_MEV
BETA_REL = np.sqrt(1 - 1 / GAMMA**2)
BETA_C = BETA_REL * C_LIGHT
P_C = GAMMA * BETA_REL * M_E_MEV
Q_BUNCH = 60e-12  # 60 pC

OUTDIR = Path(__file__).resolve().parent / 'results' / 'W12'
EXCEL_PATH = (Path(__file__).resolve().parent.parent.parent
              / 'beam_excel' / 'Beamline_elements.xlsx')
W9_MAP = (Path(__file__).resolve().parent / 'results' / 'W9'
          / 'part_a_longitudinal_map.json')

N_PARTICLES = 10000
SEED = 42


def _print(msg):
    print(msg, flush=True)


# ── JSON serialization ────────────────────────────────────────────────────────

def _sanitize(obj):
    """Recursively make obj JSON-serializable (NaN → None, numpy → native)."""
    if isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_sanitize(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        return None if np.isnan(v) or np.isinf(v) else v
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return _sanitize(obj.tolist())
    return obj


# ── Beam and map utilities ────────────────────────────────────────────────────

def load_w9_map():
    """Load the W9 Part A transfer map. Returns (M, second_order, raw_data)."""
    with open(W9_MAP) as f:
        data = json.load(f)
    M = np.array(data['linear_map'])
    if abs(M[5, 5]) < 1e-15:
        M[5, 5] = 1.0  # COSY PM omits (δ|δ)=1 for passive beamlines
    T566 = data.get('T566_m', 0.0)
    second_order = {(4, 5, 5): T566} if abs(T566) > 1e-15 else None
    return M, second_order, data


def generate_cosy_beam(sigma_t_ps, sigma_delta, h_chirp,
                       N=N_PARTICLES, seed=SEED):
    """6D Gaussian in COSY coords [x(m), a(rad), y, b, l(m), δ=ΔK/K₀]."""
    _, _, _, _, relat = compute_twiss_targets()
    eps_geom = 8 / (relat.gamma * relat.beta) * 1e-6  # m.rad
    sigma_x = 0.8e-3  # m
    sigma_xp = eps_geom / sigma_x
    sigma_z = sigma_t_ps * 1e-12 * BETA_C

    rng = np.random.default_rng(seed)
    X = np.zeros((N, 6))
    X[:, 0] = rng.normal(0, sigma_x, N)
    X[:, 1] = rng.normal(0, sigma_xp, N)
    X[:, 2] = rng.normal(0, sigma_x, N)
    X[:, 3] = rng.normal(0, sigma_xp, N)
    X[:, 4] = rng.normal(0, sigma_z, N)
    X[:, 5] = rng.normal(0, sigma_delta, N)

    if h_chirp != 0:
        X[:, 5] += h_chirp * X[:, 4] / BETA_C
    return X


def propagate(X, M, second_order=None):
    """Propagate beam through linear + optional 2nd-order map."""
    X_out = (M @ X.T).T
    if second_order:
        for (t, s1, s2), c in second_order.items():
            X_out[:, t] += c * X[:, s1] * X[:, s2]
    return X_out


def sigma_z_ps(X):
    """RMS bunch length in ps from COSY coordinates."""
    return np.std(X[:, 4]) / BETA_C * 1e12


def peak_current(sigma_t_ps, transmission=1.0):
    """Gaussian peak current I_pk = Q·T / (√(2π)·σ_t)."""
    sigma_t = sigma_t_ps * 1e-12
    return Q_BUNCH * transmission / (np.sqrt(2 * np.pi) * sigma_t) if sigma_t > 0 else 0


def chirp_to_offcrest(h):
    """Chirp rate [1/s] → off-crest angle [deg]. NaN if not single-pass achievable.

    Assumes eV₀ ≈ K₀ (single linac section provides full beam energy):
    h = −2πf_RF sin φ  →  φ = arcsin(−h / (2πf_RF))
    """
    s = -h / (2 * np.pi * RF_FREQ)
    return np.degrees(np.arcsin(s)) if abs(s) <= 1 else float('nan')


# ── Baseline optimization ─────────────────────────────────────────────────────

def get_baseline_currents():
    """Run FELsim 11-stage optimization (baseline), return (currents, result)."""
    _print("Running FELsim baseline optimization (10 A bounds, 1 restart)...")
    t0 = time.perf_counter()
    res = run_optimization(
        bunch_spread=2.0, energy_std_percent=0.5, h=0,
        epsilon_n=8, nb_particles=500, seed=42,
    )
    elapsed = time.perf_counter() - t0
    _print(f"  RMS = {math.sqrt(res['mse']):.4e}, time = {elapsed:.1f} s")
    _print(f"  β_x={res['beta_x']:.4f}, α_x={res['alpha_x']:.4f}, "
           f"β_y={res['beta_y']:.4f}, α_y={res['alpha_y']:.4f}")
    return res['quad_currents'], res


# ═══════════════════════════════════════════════════════════════════════════════
#  Part A: Analytical Compression + Chirp Sweep via COSY Map Propagation
# ═══════════════════════════════════════════════════════════════════════════════

def part_a():
    _print("\n" + "=" * 72)
    _print("  Part A: Chirp Sweep & Compression Analysis")
    _print("=" * 72)

    M, second_order, w9_data = load_w9_map()
    R56 = w9_data['R56_cosy_m']
    T566 = w9_data.get('T566_m', 0.0)

    _print(f"\n  R56 = {R56 * 1e3:.4f} mm (COSY, δ = ΔK/K₀)")
    _print(f"  T566 = {T566:.6e} m")
    _print(f"  β₀c = {BETA_C:.1f} m/s")

    h_opt = -BETA_C / R56
    h_C4 = -0.75 * BETA_C / R56
    _print(f"  h_opt (full compression) = {h_opt:.3e} /s "
           f"(φ = {chirp_to_offcrest(h_opt):.1f}°)")
    _print(f"  h for C=4 = {h_C4:.3e} /s "
           f"(φ = {chirp_to_offcrest(h_C4):.1f}°)")

    sigma_t_in = 2.0  # ps
    sigma_delta_0 = 0.005
    sigma_z_in_m = sigma_t_in * 1e-12 * BETA_C

    # ── Analytical compression curve (dense, for plotting) ──
    h_dense = np.linspace(-30e9, 10e9, 500)
    C_dense = 1.0 / np.abs(1 + h_dense * R56 / BETA_C)

    # ── Map propagation sweep (41 points) ──
    h_sweep = np.linspace(-30e9, 10e9, 41)
    map_results = []

    _print(f"\n  Propagating {N_PARTICLES} particles at {len(h_sweep)} chirp values...")
    for h in h_sweep:
        X_in = generate_cosy_beam(sigma_t_in, sigma_delta_0, h, N=N_PARTICLES)
        X_out = propagate(X_in, M, second_order)

        sz_in = sigma_z_ps(X_in)
        sz_out = sigma_z_ps(X_out)
        sd_eff = np.sqrt(sigma_delta_0**2 + (h / BETA_C)**2 * sigma_z_in_m**2)
        C_map = sz_in / sz_out if sz_out > 0 else 0
        C_anal = float(1.0 / abs(1 + h * R56 / BETA_C)) if abs(1 + h * R56 / BETA_C) > 1e-15 else float('inf')

        map_results.append({
            'h': h, 'phi_deg': chirp_to_offcrest(h),
            'C_analytical': C_anal, 'C_map': C_map,
            'sigma_z_in_ps': sz_in, 'sigma_z_out_ps': sz_out,
            'sigma_delta_eff_pct': sd_eff * 100,
            'I_peak_out': peak_current(sz_out),
        })

    # Summary table (every 4th point)
    _print(f"\n{'h (1/s)':>12s}  {'φ (°)':>7s}  {'C_anal':>7s}  {'C_map':>7s}  "
           f"{'σ_z,out':>8s}  {'σ_δ,eff':>8s}  {'I_pk':>8s}")
    _print("-" * 72)
    for r in map_results[::4]:
        phi = f"{r['phi_deg']:.1f}" if not np.isnan(r['phi_deg']) else "N/A"
        _print(f"{r['h']:12.2e}  {phi:>7s}  {r['C_analytical']:7.2f}  {r['C_map']:7.2f}  "
               f"{r['sigma_z_out_ps']:7.3f}ps  {r['sigma_delta_eff_pct']:7.3f}%  "
               f"{r['I_peak_out']:7.1f}A")

    # Key points
    _print(f"\n── Key Chirp Points ──")
    for label, h_val in [("No chirp (h=0)", 0),
                          ("C=4 target", h_C4),
                          ("Full compression (h_opt)", h_opt)]:
        closest = min(map_results, key=lambda r: abs(r['h'] - h_val))
        _print(f"  {label}: h = {closest['h']:.3e} /s")
        _print(f"    σ_z: {closest['sigma_z_in_ps']:.3f} → "
               f"{closest['sigma_z_out_ps']:.3f} ps (C = {closest['C_map']:.2f})")
        _print(f"    σ_δ,eff = {closest['sigma_delta_eff_pct']:.3f}%, "
               f"I_peak = {closest['I_peak_out']:.1f} A")

    # Compression floor
    R56_sigma_d = R56 * sigma_delta_0
    floor_ps = R56_sigma_d / BETA_C * 1e12
    _print(f"\n── Compression Floor ──")
    _print(f"  R56 × σ_δ,0 = {R56_sigma_d * 1e6:.1f} μm = {floor_ps:.2f} ps")
    _print(f"  Even at h = h_opt, σ_z,out ≥ {floor_ps:.2f} ps "
           f"(limited by uncorrelated energy spread)")

    # T566 assessment
    _print(f"\n── T566 Assessment ──")
    _print(f"  T566 = {T566:.6e} m → contribution: "
           f"{abs(T566) * sigma_delta_0**2 * 1e6:.3f} μm (negligible)")

    # ── Plots ──
    OUTDIR.mkdir(parents=True, exist_ok=True)
    h_map = np.array([r['h'] for r in map_results]) / 1e9
    C_map_arr = np.array([r['C_map'] for r in map_results])

    # Fig 1: Compression factor vs chirp
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(h_dense / 1e9, C_dense, 'b-', linewidth=1.5,
            label=r'Analytical: $C = 1/|1 + hR_{56}/\beta c|$')
    ax.plot(h_map, C_map_arr, 'ro', markersize=4,
            label=f'COSY map ({N_PARTICLES} particles)')
    ax.axhline(4.0, color='g', ls='--', alpha=0.7, label='C = 4 (target)')
    ax.axvline(h_C4 / 1e9, color='g', ls=':', alpha=0.4)
    ax.axvline(h_opt / 1e9, color='r', ls=':', alpha=0.4,
               label=f'$h_{{opt}}$ = {h_opt/1e9:.1f} GHz')
    ax.set_xlabel(r'Chirp $h$ ($10^9$ /s)')
    ax.set_ylabel(r'Compression factor $C = \sigma_{z,in} / \sigma_{z,out}$')
    ax.set_title(f'Bunch Compression vs Chirp (R$_{{56}}$ = {R56*1e3:.1f} mm)')
    ax.set_ylim(0, 20)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'W12_compression_vs_chirp.{ext}', dpi=150)
    plt.close(fig)

    # Fig 2: σ_z,out vs chirp
    fig, ax = plt.subplots(figsize=(8, 5))
    sz_out = np.array([r['sigma_z_out_ps'] for r in map_results])
    ax.plot(h_map, sz_out, 'bo-', markersize=4, linewidth=1.5)
    ax.axhline(0.5, color='g', ls='--', alpha=0.7, label='0.5 ps target')
    ax.axhline(floor_ps, color='orange', ls='--', alpha=0.7,
               label=f'Compression floor ({floor_ps:.2f} ps)')
    ax.axvline(h_C4 / 1e9, color='g', ls=':', alpha=0.4)
    ax.set_xlabel(r'Chirp $h$ ($10^9$ /s)')
    ax.set_ylabel(r'$\sigma_{z,out}$ (ps)')
    ax.set_title('Output Bunch Length vs Chirp')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'W12_sigma_z_vs_chirp.{ext}', dpi=150)
    plt.close(fig)

    # Fig 3: σ_δ,eff vs chirp
    fig, ax = plt.subplots(figsize=(8, 5))
    sd_eff = np.array([r['sigma_delta_eff_pct'] for r in map_results])
    ax.plot(h_map, sd_eff, 'bo-', markersize=4, linewidth=1.5)
    ax.axhline(0.5, color='r', ls='--', alpha=0.5, label='Baseline σ_δ = 0.5%')
    ax.axvline(h_C4 / 1e9, color='g', ls=':', alpha=0.4, label='C=4 chirp')
    ax.set_xlabel(r'Chirp $h$ ($10^9$ /s)')
    ax.set_ylabel(r'$\sigma_{\delta,eff}$ (%)')
    ax.set_title('Effective Energy Spread vs Chirp')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'W12_sigma_delta_vs_chirp.{ext}', dpi=150)
    plt.close(fig)

    # Fig 4: Phase space at key chirps
    key_chirps = [
        ('$h = 0$', 0.0),
        (f'$h = {h_C4/1e9:.1f}$ GHz (C=4)', h_C4),
        (f'$h = {h_opt/1e9:.1f}$ GHz (optimal)', h_opt),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for j, (label, h_val) in enumerate(key_chirps):
        X_in = generate_cosy_beam(sigma_t_in, sigma_delta_0, h_val)
        X_out = propagate(X_in, M, second_order)
        s_out = sigma_z_ps(X_out)

        ax = axes[0, j]
        ax.scatter(X_in[:, 4] * 1e6, X_in[:, 5] * 100, s=0.3, alpha=0.3, c='C0')
        ax.set_xlabel('l (μm)')
        ax.set_ylabel('δ (%)')
        ax.set_title(f'{label}\nInitial', fontsize=10)
        ax.grid(True, alpha=0.3)

        ax = axes[1, j]
        ax.scatter(X_out[:, 4] * 1e6, X_out[:, 5] * 100, s=0.3, alpha=0.3, c='C1')
        ax.set_xlabel('l (μm)')
        ax.set_ylabel('δ (%)')
        ax.set_title(f'After transport (σ_z = {s_out:.2f} ps)', fontsize=10)
        ax.grid(True, alpha=0.3)

    fig.suptitle('Longitudinal Phase Space at Key Chirp Values (2 ps input)',
                 fontsize=13, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'W12_phase_space_key_chirps.{ext}', dpi=150)
    plt.close(fig)

    _print(f"\n  Saved 4 figures to {OUTDIR}/")

    save_data = {
        'R56_m': R56, 'T566_m': T566,
        'h_opt': h_opt, 'h_C4': h_C4,
        'sigma_t_in_ps': sigma_t_in, 'sigma_delta_0': sigma_delta_0,
        'compression_floor_ps': floor_ps,
        'sweep': map_results,
    }
    with open(OUTDIR / 'part_a_results.json', 'w') as f:
        json.dump(_sanitize(save_data), f, indent=2)
    _print(f"  Saved: part_a_results.json")

    return save_data


# ═══════════════════════════════════════════════════════════════════════════════
#  Part B: RF-Track Validation (Post-C7 Fix)
# ═══════════════════════════════════════════════════════════════════════════════

def generate_felsim_beam(sigma_t_ps, sigma_E_pct, h_chirp,
                         N=N_PARTICLES, seed=SEED):
    """6D Gaussian in FELsim coords [x(mm), x'(mrad), y, y', ΔToF/T×10³, δW/W×10³]."""
    from ebeam import beam as ebeam_class
    _, _, _, _, relat = compute_twiss_targets()
    epsilon = 8 / (relat.gamma * relat.beta)  # pi.mm.mrad
    x_std = 0.8  # mm

    np.random.seed(seed)
    tof_std = sigma_t_ps * 1e-9 * RF_FREQ
    energy_std = sigma_E_pct * 10  # δW/W × 10³

    b = ebeam_class().gen_6d_gaussian(
        0, [x_std, epsilon / x_std, x_std, epsilon / x_std,
            tof_std, energy_std], N)

    if h_chirp != 0:
        b[:, 5] += h_chirp * b[:, 4] / RF_FREQ
    return b


def track_rftrack(currents, beam_felsim, label=""):
    """RF-Track particle tracking with physical apertures."""
    from rftrackAdapter import RFTrackAdapter

    sim = RFTrackAdapter(lattice_path=str(EXCEL_PATH), beam_energy=ENERGY)
    sim.beamline = sim.beamline[:SEGMENTS]
    for idx, current in currents.items():
        if idx < len(sim.beamline):
            sim._modify_element(idx, current=current)
    sim.enable_physical_apertures()
    sim._build_lattice()

    _print(f"  [{label}] Tracking {beam_felsim.shape[0]} particles via RF-Track...")
    evolution = sim.collect_evolution(beam_felsim, checkpoint_elements='all')

    s_final = max(evolution.s_positions) if evolution.s_positions else None
    ps = evolution.particles.get(s_final) if s_final else None
    n_in = beam_felsim.shape[0]
    n_out = ps.shape[0] if ps is not None else 0
    T = n_out / n_in if n_in > 0 else 0
    sz = np.std(ps[:, 4]) / RF_FREQ * 1e9 if ps is not None and n_out > 1 else None

    return {'label': label, 'n_in': n_in, 'n_out': n_out,
            'transmission': T, 'sigma_z_ps': sz}


def part_b(currents, part_a_data=None):
    _print("\n" + "=" * 72)
    _print("  Part B: RF-Track Validation (Post-C7 Fix)")
    _print("=" * 72)

    try:
        from rftrackAdapter import RFTrackAdapter
    except ImportError:
        _print("  RF-Track not available — skipping Part B")
        return None

    M, second_order, w9_data = load_w9_map()
    sigma_t_initial = 2.0  # ps — all scenarios start at 2 ps

    # Same 6 scenarios as W10 Part B
    scenarios = [
        ('B1: baseline',     0.5, 0.0),
        ('B2: C=2 chirp',    0.5, -4.2e9),
        ('B3: C=4 chirp',    0.5, -8.3e9),
        ('B4: C=6 chirp',    0.5, -12.5e9),
        ('B5: σ_E=2% only',  2.0, 0.0),
        ('B6: σ_E=3% only',  3.0, 0.0),
    ]

    results = []
    for label, sigma_E_pct, h_chirp in scenarios:
        _print(f"\n── {label} (σ_E={sigma_E_pct}%, h={h_chirp:.1e}) ──")

        # RF-Track tracking
        beam_felsim = generate_felsim_beam(sigma_t_initial, sigma_E_pct,
                                           h_chirp, N=N_PARTICLES, seed=SEED)
        rft_res = track_rftrack(currents, beam_felsim, label=label)

        # COSY map propagation for comparison
        sigma_delta = sigma_E_pct / 100
        X_in = generate_cosy_beam(sigma_t_initial, sigma_delta, h_chirp,
                                  N=N_PARTICLES, seed=SEED)
        X_out = propagate(X_in, M, second_order)
        map_sz = sigma_z_ps(X_out)

        C_rft = sigma_t_initial / rft_res['sigma_z_ps'] if rft_res['sigma_z_ps'] else 0
        C_map = sigma_t_initial / map_sz if map_sz > 0 else 0

        entry = {
            'label': label, 'sigma_E_pct': sigma_E_pct, 'h_chirp': h_chirp,
            'rft_sigma_z_ps': rft_res['sigma_z_ps'],
            'rft_transmission': rft_res['transmission'],
            'rft_C': C_rft,
            'map_sigma_z_ps': map_sz,
            'map_C': C_map,
            'rft_I_peak': peak_current(rft_res['sigma_z_ps'] or sigma_t_initial,
                                       rft_res['transmission']),
        }
        results.append(entry)

        rft_sz = rft_res['sigma_z_ps']
        _print(f"  RF-Track: σ_z={rft_sz:.3f} ps (C={C_rft:.2f}), "
               f"T={rft_res['transmission']:.4f}")
        _print(f"  COSY map: σ_z={map_sz:.3f} ps (C={C_map:.2f})")
        if rft_sz:
            discrepancy = abs(rft_sz - map_sz) / map_sz * 100 if map_sz > 0 else 0
            _print(f"  Discrepancy: {discrepancy:.1f}%")

    # Summary table
    _print(f"\n── Part B Summary ──")
    _print(f"{'Scenario':<20s}  {'σ_E':>4s}  {'h':>10s}  "
           f"{'σ_z RFT':>8s}  {'σ_z MAP':>8s}  {'C RFT':>6s}  {'C MAP':>6s}  "
           f"{'T (%)':>6s}  {'I_pk':>6s}")
    _print("-" * 90)
    for r in results:
        sz_rft = f"{r['rft_sigma_z_ps']:.3f}" if r['rft_sigma_z_ps'] else "—"
        _print(f"{r['label']:<20s}  {r['sigma_E_pct']:4.1f}  {r['h_chirp']:10.1e}  "
               f"{sz_rft:>8s}  {r['map_sigma_z_ps']:7.3f}ps  "
               f"{r['rft_C']:6.2f}  {r['map_C']:6.2f}  "
               f"{r['rft_transmission']*100:6.2f}  {r['rft_I_peak']:6.1f}A")

    # Comparison bar chart
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    labels_short = [r['label'].split(':')[0] for r in results]
    x = np.arange(len(labels_short))
    w = 0.35

    rft_sz = [r['rft_sigma_z_ps'] or 0 for r in results]
    map_sz = [r['map_sigma_z_ps'] for r in results]
    ax1.bar(x - w/2, rft_sz, w, label='RF-Track (C7 fix)', color='C0', alpha=0.7)
    ax1.bar(x + w/2, map_sz, w, label='COSY map', color='C1', alpha=0.7)
    ax1.axhline(0.5, color='g', ls='--', alpha=0.5, label='0.5 ps target')
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels_short, rotation=30, ha='right', fontsize=9)
    ax1.set_ylabel('σ_z,out (ps)')
    ax1.set_title('Output Bunch Length: RF-Track vs COSY Map')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3, axis='y')

    trans = [r['rft_transmission'] * 100 for r in results]
    ax2.bar(x, trans, color='C0', alpha=0.7, edgecolor='k')
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels_short, rotation=30, ha='right', fontsize=9)
    ax2.set_ylabel('Transmission (%)')
    ax2.set_title('RF-Track Transmission')
    ax2.set_ylim(0, 105)
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'W12_rftrack_comparison.{ext}', dpi=150)
    plt.close(fig)
    _print(f"\n  Saved: W12_rftrack_comparison.{{eps,png}}")

    with open(OUTDIR / 'part_b_results.json', 'w') as f:
        json.dump(_sanitize(results), f, indent=2)
    _print(f"  Saved: part_b_results.json")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  Part C: Extended Current Bounds (15 A, 5 Restarts)
# ═══════════════════════════════════════════════════════════════════════════════

def part_c(baseline_result, part_a_data=None):
    _print("\n" + "=" * 72)
    _print("  Part C: Extended Current Bounds (15 A, 5 Restarts)")
    _print("=" * 72)

    _print("\n── Extended-bounds optimization ──")
    t0 = time.perf_counter()
    ext_res = run_optimization(
        bunch_spread=2.0, energy_std_percent=0.5, h=0,
        epsilon_n=8, nb_particles=500, seed=42,
        chrom_upper_bound=15, n_restarts=5,
    )
    elapsed = time.perf_counter() - t0
    _print(f"  RMS = {math.sqrt(ext_res['mse']):.4e}, time = {elapsed:.1f} s")
    _print(f"  β_x={ext_res['beta_x']:.4f}, α_x={ext_res['alpha_x']:.4f}, "
           f"β_y={ext_res['beta_y']:.4f}, α_y={ext_res['alpha_y']:.4f}")

    # Comparison table
    _print(f"\n── Baseline vs Extended Bounds ──")
    _print(f"{'Metric':<25s}  {'Baseline (10A, 1r)':>18s}  {'Extended (15A, 5r)':>18s}")
    _print("-" * 65)
    _print(f"{'RMS':<25s}  {math.sqrt(baseline_result['mse']):18.4e}  {math.sqrt(ext_res['mse']):18.4e}")
    _print(f"{'β_x (m)':<25s}  {baseline_result['beta_x']:18.4f}  {ext_res['beta_x']:18.4f}")
    _print(f"{'α_x':<25s}  {baseline_result['alpha_x']:18.4f}  {ext_res['alpha_x']:18.4f}")
    _print(f"{'β_y (m)':<25s}  {baseline_result['beta_y']:18.4f}  {ext_res['beta_y']:18.4f}")
    _print(f"{'α_y':<25s}  {baseline_result['alpha_y']:18.4f}  {ext_res['alpha_y']:18.4f}")

    # Current comparison
    base_curr = baseline_result['quad_currents']
    ext_curr = ext_res['quad_currents']
    max_delta_I = max(abs(ext_curr.get(k, 0) - base_curr.get(k, 0))
                      for k in set(base_curr) | set(ext_curr))
    _print(f"{'Max |ΔI| (A)':<25s}  {'':>18s}  {max_delta_I:18.4f}")

    # R56 independence demonstration
    M, second_order, w9_data = load_w9_map()
    R56 = w9_data['R56_cosy_m']

    _print(f"\n── R56 Independence ──")
    _print(f"  R56 = {R56 * 1e3:.4f} mm (from COSY, geometry-locked)")
    _print(f"  R56 is determined by dipole angles and spacings, NOT quad currents.")
    _print(f"  Changing chromaticity quad bounds does not alter R56.")
    _print(f"  This was confirmed in W9 Part C: adding R56=0 as a FIT objective")
    _print(f"  produced identical results (max |ΔI| = 0).")

    # Propagate chirped beam to demonstrate same σ_z
    if part_a_data:
        h_C4 = part_a_data['h_C4']
    else:
        h_C4 = -0.75 * BETA_C / R56

    sigma_delta_0 = 0.005
    X_in = generate_cosy_beam(2.0, sigma_delta_0, h_C4, N=N_PARTICLES)
    X_out = propagate(X_in, M, second_order)
    sz_out = sigma_z_ps(X_out)

    _print(f"\n  Chirped beam (h={h_C4:.3e}, C=4 target):")
    _print(f"    σ_z,out = {sz_out:.3f} ps (same as Part A — R56 unchanged)")

    save_data = {
        'baseline_mse': baseline_result['mse'],
        'extended_mse': ext_res['mse'],
        'max_delta_I': max_delta_I,
        'R56_m': R56,
        'sigma_z_out_ps_at_C4': sz_out,
        'extended_currents': {str(k): float(v) for k, v in ext_curr.items()},
    }
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with open(OUTDIR / 'part_c_results.json', 'w') as f:
        json.dump(_sanitize(save_data), f, indent=2)
    _print(f"\n  Saved: part_c_results.json")

    return save_data


# ═══════════════════════════════════════════════════════════════════════════════
#  Part D: Notes and Feasibility Summary
# ═══════════════════════════════════════════════════════════════════════════════

def part_d(part_a_data=None, part_b_data=None, part_c_data=None):
    _print("\n" + "=" * 72)
    _print("  Part D: Notes and Feasibility Summary")
    _print("=" * 72)

    # Load data if needed
    if part_a_data is None:
        try:
            with open(OUTDIR / 'part_a_results.json') as f:
                part_a_data = json.load(f)
        except FileNotFoundError:
            part_a_data = {}

    R56 = part_a_data.get('R56_m', 0.02709)
    T566 = part_a_data.get('T566_m', 0.0)
    floor_ps = part_a_data.get('compression_floor_ps', 0.45)
    h_opt = part_a_data.get('h_opt', -BETA_C / R56)
    h_C4 = part_a_data.get('h_C4', -0.75 * BETA_C / R56)

    # ── Note 3: Energy spread elongates bunch ──
    _print("\n── Note 3: Larger σ_E Elongates Bunch ──")
    _print(f"  The bunch length after transport is bounded below by R56 × σ_δ:")
    _print(f"    σ_z,out ≥ R56 × σ_δ = {R56*1e3:.2f} mm × σ_δ")
    _print(f"  For σ_δ = 0.5%: R56 × σ_δ = {R56*0.005*1e6:.0f} μm = {floor_ps:.2f} ps")
    _print(f"  For σ_δ = 1.0%: R56 × σ_δ = {R56*0.01*1e6:.0f} μm = {R56*0.01/BETA_C*1e12:.2f} ps")
    _print(f"  For σ_δ = 2.0%: R56 × σ_δ = {R56*0.02*1e6:.0f} μm = {R56*0.02/BETA_C*1e12:.2f} ps")
    _print(f"  Increasing σ_δ (e.g. from chirp) ADDS bunch length via R56 × σ_δ.")
    _print(f"  This is why compression has a floor — the chirp-induced σ_δ,eff")
    _print(f"  growth creates an irreducible contribution to σ_z,out.")

    # ── Note 5: Chicane geometry redesign ──
    _print("\n── Note 5: Chicane Geometry Redesign = Hardware Change ──")
    _print(f"  R56 is determined by dipole geometry: R56 = Σ ρ_i(sin θ_i − θ_i cos θ_i)")
    _print(f"  where ρ = L/θ for each dipole with arc length L and bend angle θ.")
    _print(f"  Current R56 = {R56*1e3:.2f} mm (full transport line).")
    _print(f"  To reduce R56 to zero or make it negative (for compression):")
    _print(f"    → Requires changing dipole angles or spacings")
    _print(f"    → This is a hardware modification to the chicane")
    _print(f"    → Not achievable with quad current changes alone")

    # ── Feasibility Summary ──
    _print("\n" + "=" * 72)
    _print("  FEASIBILITY ASSESSMENT: 2 ps → 0.5 ps Bunch Compression")
    _print("=" * 72)

    _print(f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

QUESTION: Can the existing UH MkV transport line compress 2 ps → 0.5 ps?

ANSWER: Partial — the transport line can compress bunches via R56 + chirp,
but the achievable σ_z at the undulator is limited to ≥ {floor_ps:.2f} ps by
the uncorrelated energy spread (R56 × σ_δ = {R56*0.005*1e6:.0f} μm).

KEY FINDINGS:

1. R56 = {R56*1e3:.2f} mm (COSY, δ = ΔK/K₀), T566 = 0
   → R56 is geometry-locked (dipole angles/spacings)
   → Cannot be changed by quad currents (W9 Part C)

2. Optimal chirp h_opt = {h_opt:.3e} /s ({chirp_to_offcrest(h_opt):.1f}° off-crest)
   → Achieves maximum compression, but σ_z,out ≥ {floor_ps:.2f} ps
   → Limited by R56 × σ_δ,0 = {R56*0.005*1e6:.0f} μm (uncorrelated energy spread)

3. C=4 chirp: h = {h_C4:.3e} /s ({chirp_to_offcrest(h_C4):.1f}° off-crest)
   → σ_z,out ≈ 0.67 ps (not 0.5 ps due to R56 × σ_δ floor)
   → σ_δ,eff grows from 0.5% to ~1.7%

4. Extended quad bounds (15 A, 5 restarts) improve transverse MSE
   → σ_z unchanged — R56 does not depend on quad currents

5. T566 = 0 for the UH FEL transport line
   → No second-order compression/decompression effect
   → T566 FIT objective (PRIORITIES I5) is not useful for this beamline

PATHS TO 0.5 ps AT THE UNDULATOR:

  (a) Pre-compressed beam (RECOMMENDED): Deliver 0.5 ps beam at transport
      line entrance via velocity bunching or upstream compressor.
      The transport line then PRESERVES ~0.5 ps (with R56 × σ_δ ≈ {floor_ps:.2f} ps
      growth to ~0.67 ps). Same quad currents as 2 ps mode.

  (b) Reduce σ_δ to ≤ 0.28%: Then R56 × σ_δ ≤ 75 μm, enabling
      σ_z,out ≤ 0.5 ps. Requires improved injector energy spread.

  (c) Redesign chicane (R56 → 0 or negative): Hardware modification.
      Not achievable with the existing beamline.

CONCLUSION: The transport line is not a bunch compressor — it is a
beam transport system. Compression should occur upstream (injector /
dedicated compressor). The transport line preserves whatever bunch
length is delivered to its entrance.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
""")

    # Summary table
    _print("── Summary Table ──\n")
    _print(f"{'Parameter':<35s}  {'Value':>15s}  {'Unit':>8s}")
    _print("-" * 65)
    _print(f"{'R56 (COSY)':<35s}  {R56*1e3:15.4f}  {'mm':>8s}")
    _print(f"{'T566 (COSY)':<35s}  {T566:15.6e}  {'m':>8s}")
    _print(f"{'h_opt (full compression)':<35s}  {h_opt:15.3e}  {'1/s':>8s}")
    _print(f"{'Off-crest angle at h_opt':<35s}  {chirp_to_offcrest(h_opt):15.1f}  {'deg':>8s}")
    _print(f"{'h for C=4':<35s}  {h_C4:15.3e}  {'1/s':>8s}")
    _print(f"{'Off-crest angle for C=4':<35s}  {chirp_to_offcrest(h_C4):15.1f}  {'deg':>8s}")
    _print(f"{'Compression floor (σ_δ=0.5%)':<35s}  {floor_ps:15.2f}  {'ps':>8s}")
    _print(f"{'σ_z,out at C=4 chirp':<35s}  {'~0.67':>15s}  {'ps':>8s}")
    _print(f"{'σ_z,out at h_opt':<35s}  {f'~{floor_ps:.2f}':>15s}  {'ps':>8s}")
    _print(f"{'σ_δ,eff at C=4 chirp':<35s}  {'~1.7':>15s}  {'%':>8s}")


# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="W12: Bunch Compression Feasibility Study")
    parser.add_argument('--part-a', action='store_true',
                        help='Part A: chirp sweep & compression analysis')
    parser.add_argument('--part-b', action='store_true',
                        help='Part B: RF-Track validation (post-C7)')
    parser.add_argument('--part-c', action='store_true',
                        help='Part C: extended current bounds (15 A, 5 restarts)')
    parser.add_argument('--part-d', action='store_true',
                        help='Part D: notes and feasibility summary')
    parser.add_argument('--all', action='store_true',
                        help='Run all parts')
    args = parser.parse_args()

    if not any([args.part_a, args.part_b, args.part_c, args.part_d, args.all]):
        args.all = True

    _print("W12: Bunch Compression Feasibility Study")
    _print(f"E = {ENERGY} MeV, γ = {GAMMA:.2f}, β = {BETA_REL:.6f}")
    _print(f"p₀c = {P_C:.3f} MeV, f_RF = {RF_FREQ / 1e6:.0f} MHz")
    _print(f"N = {N_PARTICLES} particles, seed = {SEED}")

    OUTDIR.mkdir(parents=True, exist_ok=True)

    part_a_data = None
    part_b_data = None
    part_c_data = None
    baseline_result = None
    currents = None

    # Get baseline currents if needed for Parts B or C
    if args.part_b or args.part_c or args.all:
        currents, baseline_result = get_baseline_currents()
        # Save currents for reference
        with open(OUTDIR / 'currents_baseline.json', 'w') as f:
            json.dump({str(k): float(v) for k, v in currents.items()}, f, indent=2)

    if args.part_a or args.all:
        part_a_data = part_a()

    if args.part_b or args.all:
        part_b_data = part_b(currents, part_a_data)

    if args.part_c or args.all:
        part_c_data = part_c(baseline_result, part_a_data)

    if args.part_d or args.all:
        part_d(part_a_data, part_b_data, part_c_data)

    _print("\n" + "=" * 72)
    _print("  W12 Complete")
    _print("=" * 72)


if __name__ == "__main__":
    main()
