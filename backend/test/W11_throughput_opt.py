"""W11: Throughput Optimization — Maximize Peak Current at Undulator

Extends the Twiss-only Stage 11 optimization to include transmission and bunch
length objectives. Two scenarios:
  - 2 ps → 2 ps (transport mode): preserve bunch length, maximize throughput
  - 2 ps → 0.5 ps (compression mode): compress via chirp, maximize peak current

Three-code comparison: FELsim (transfer matrices), RF-Track (particle tracking),
COSY INFINITY (DA maps with FIT).

Architecture:
  Stages 1-10 optimized by FELsim (fast, correct for quads/drifts).
  Stage 11 (4 quads: 87, 93, 95, 97) optimized with RF-Track particle tracking
  using a weighted scalar objective: w_t × MSE_Twiss + w_T × (1-T)² + w_σ × (σ_t/σ_target - 1)².

Author: Eremey Valetov
"""

import sys
import json
import time
import argparse
import csv
from pathlib import Path
import numpy as np
import scipy.optimize

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ebeam import beam
from simulatorBase import CoordinateSystem

sys.path.insert(0, str(Path(__file__).resolve().parent))
from UHM_beamline_opt_05ps_params import (
    run_optimization, compute_twiss_targets, QUAD_INDICES,
    ENERGY, RF_FREQ, SEGMENTS, MSE_THRESHOLDS, write_csv,
)
from UHM_rftrack_opt import (
    setup_rftrack_adapter, pretrack_prefix, _track_suffix, _compute_mse,
    S11_QUAD_INDICES, PREFIX_END, create_beam,
    rftrack_mse_cached, evaluate_felsim_in_rftrack,
)
from W10_beam_losses_compression import (
    create_felsim_beam, create_cosy_beam, peak_current,
    run_cosy_tracking, run_rftrack_tracking,
    beam_sigma_z_ps, C_LIGHT, BETA_REL, GAMMA_REL, P_C,
    M_E_MEV, QUAD_HALF_APERTURE, DIPOLE_HALF_GAP,
    DIPOLE_HALF_WIDTH, BEAM_PIPE_RADIUS,
)

try:
    import RF_Track as rft
    from rftrackAdapter import RFTrackAdapter
    _RFTRACK_AVAILABLE = True
except ImportError:
    _RFTRACK_AVAILABLE = False

# ── Constants ────────────────────────────────────────────────────────────────
OUTDIR = Path(__file__).resolve().parent / 'results' / 'W11'
EXCEL_PATH = (Path(__file__).resolve().parent.parent.parent
              / 'beam_excel' / 'Beamline_elements.xlsx')

N_PARTICLES = 500
SEED = 42

# Scenario definitions
SCENARIOS = {
    '2ps_transport': {
        'h_chirp': 0.0,
        'sigma_t_target_ps': 2.0,
        'weights': (1.0, 1.0, 0.3),  # (w_twiss, w_transmission, w_sigma)
        'description': '2 ps transport mode (preserve bunch length)',
    },
    '05ps_compress': {
        'h_chirp': -8.3e9,
        'sigma_t_target_ps': 0.5,
        'weights': (0.5, 0.8, 1.5),  # heavier weight on compression
        'description': '0.5 ps compression mode (chirp h = -8.3e9/s)',
    },
}


def _print(msg):
    print(msg, flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
#  RF-Track throughput objective
# ═══════════════════════════════════════════════════════════════════════════════

def throughput_objective(x, sim, beam_rft_cached, targets, ebeam_obj,
                         sigma_t_target_ps, weights, counter):
    """NM objective combining Twiss MSE, transmission, and bunch length.

    Parameters
    ----------
    x : array-like, shape (4,)
        Stage 11 quad currents [I_87, I_93, I_95, I_97].
    sim : RFTrackAdapter
        Adapter with beamline already configured for stages 1-10.
    beam_rft_cached : ndarray
        RF-Track coordinate beam state at element 87 (prefix exit).
    targets : dict
        Twiss targets (alpha_x, alpha_y, beta_x, beta_y).
    ebeam_obj : beam
        FELsim beam analysis object.
    sigma_t_target_ps : float
        Target RMS bunch length in ps.
    weights : tuple of 3 floats
        (w_twiss, w_transmission, w_sigma).
    counter : list of [int]
        Mutable evaluation counter.

    Returns
    -------
    float
        Weighted cost.
    """
    counter[0] += 1
    for idx, current in zip(S11_QUAD_INDICES, x):
        sim._modify_element(idx, current=current)

    # Track suffix 87:118 → final beam state
    ps_final = _track_suffix(sim, beam_rft_cached, PREFIX_END, SEGMENTS)
    if ps_final is None:
        return 1e6

    n_surviving = ps_final.shape[0]
    n_entering = beam_rft_cached.shape[0]
    T = n_surviving / n_entering

    # Twiss at undulator entrance
    _, _, twiss_f = ebeam_obj.cal_twiss(ps_final)
    bx = twiss_f.loc['x'][r"$\beta$ (m)"]
    ax = twiss_f.loc['x'][r"$\alpha$"]
    by = twiss_f.loc['y'][r"$\beta$ (m)"]
    ay = twiss_f.loc['y'][r"$\alpha$"]

    # Dispersion at element 92 (suffix 87:93)
    ps_92 = _track_suffix(sim, beam_rft_cached, PREFIX_END, 93)
    disp = 0.0
    if ps_92 is not None:
        _, _, tw92 = ebeam_obj.cal_twiss(ps_92)
        disp = tw92.loc['x'][r"$D$ (m)"]

    mse_twiss = _compute_mse(bx, ax, by, ay, disp, targets)

    # Bunch length from FELsim coordinates (col 4 = ΔToF/T × 10³)
    sigma_t_ps = np.std(ps_final[:, 4]) / RF_FREQ * 1e9

    w_t, w_T, w_s = weights
    cost = (w_t * mse_twiss
            + w_T * (1 - T)**2
            + w_s * (sigma_t_ps / sigma_t_target_ps - 1)**2)

    if counter[0] % 50 == 0:
        _print(f"    eval {counter[0]}: cost={cost:.4e} "
               f"(MSE={mse_twiss:.4e}, T={T:.3f}, σ_t={sigma_t_ps:.3f} ps)")

    return cost


# ═══════════════════════════════════════════════════════════════════════════════
#  Scenario runner
# ═══════════════════════════════════════════════════════════════════════════════

def run_scenario(scenario_name, felsim_currents, h_chirp, sigma_t_target_ps,
                 weights, n_restarts=5, nb_particles=500, seed=42,
                 chrom_upper_bound=15):
    """Run throughput optimization for a single scenario.

    Returns dict with results or None on failure.
    """
    _print(f"\n{'='*72}")
    _print(f"  Scenario: {scenario_name}")
    _print(f"  h = {h_chirp:.1e}/s, σ_t target = {sigma_t_target_ps} ps")
    _print(f"  weights = {weights}")
    _print(f"{'='*72}")

    beta_xm, alpha_xm, beta_ym, alpha_ym, _ = compute_twiss_targets()
    targets = {
        'alpha_x': alpha_xm, 'alpha_y': alpha_ym,
        'beta_x': beta_xm, 'beta_y': beta_ym,
    }

    qb = 10
    cb = chrom_upper_bound
    ebeam_obj = beam()

    # Create beam with appropriate chirp
    beam_dist = create_felsim_beam(
        sigma_t_ps=2.0, sigma_E_pct=0.5, h_chirp=h_chirp,
        epsilon_n=8, N=nb_particles, seed=seed)

    # Set up RF-Track with physical apertures
    sim = setup_rftrack_adapter(felsim_currents, aperture=0.5)
    sim.enable_physical_apertures()
    sim._build_lattice()

    beam_rft = sim.transform_coordinates(
        beam_dist, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK)

    # Pre-track prefix (elements 0:87)
    _print(f"  Pre-tracking prefix (elements 0:{PREFIX_END})...")
    t_prefix = time.perf_counter()
    beam_rft_cached = pretrack_prefix(sim, beam_rft)
    t_prefix = time.perf_counter() - t_prefix
    if beam_rft_cached is None:
        _print(f"  FAILED: prefix tracking lost too many particles")
        return None
    _print(f"  Prefix tracked in {t_prefix:.1f} s "
           f"({beam_rft_cached.shape[0]}/{beam_rft.shape[0]} surviving)")

    # Multi-start NM optimization
    felsim_s11 = [felsim_currents[i] for i in S11_QUAD_INDICES]
    rng = np.random.RandomState(seed + 999)
    starts = [felsim_s11]
    for _ in range(n_restarts - 1):
        starts.append([rng.uniform(0, cb), rng.uniform(0, qb),
                       rng.uniform(0, qb), rng.uniform(0, qb)])

    bounds = [(0, cb), (0, qb), (0, qb), (0, qb)]
    best_result = None
    best_x = None
    total_nfev = 0
    t0 = time.perf_counter()

    for i, x0 in enumerate(starts):
        counter = [0]
        result = scipy.optimize.minimize(
            throughput_objective, x0, method='Nelder-Mead',
            bounds=bounds,
            args=(sim, beam_rft_cached, targets, ebeam_obj,
                  sigma_t_target_ps, weights, counter),
        )
        total_nfev += counter[0]
        _print(f"    Restart {i+1}/{len(starts)}: cost = {result.fun:.4e} "
               f"(nfev={counter[0]})")

        if best_result is None or result.fun < best_result.fun:
            best_result = result
            best_x = result.x.copy()

    elapsed = time.perf_counter() - t0

    # Final evaluation with best solution
    for idx, current in zip(S11_QUAD_INDICES, best_x):
        sim._modify_element(idx, current=current)

    ps_final = _track_suffix(sim, beam_rft_cached, PREFIX_END, SEGMENTS)
    if ps_final is None:
        _print(f"  FAILED: final evaluation lost all particles")
        return None

    n_surviving = ps_final.shape[0]
    T = n_surviving / beam_rft_cached.shape[0]
    _, _, twiss_f = ebeam_obj.cal_twiss(ps_final)

    bx = twiss_f.loc['x'][r"$\beta$ (m)"]
    ax = twiss_f.loc['x'][r"$\alpha$"]
    by = twiss_f.loc['y'][r"$\beta$ (m)"]
    ay = twiss_f.loc['y'][r"$\alpha$"]

    ps_92 = _track_suffix(sim, beam_rft_cached, PREFIX_END, 93)
    disp = 0.0
    if ps_92 is not None:
        _, _, tw92 = ebeam_obj.cal_twiss(ps_92)
        disp = tw92.loc['x'][r"$D$ (m)"]

    mse_twiss = _compute_mse(bx, ax, by, ay, disp, targets)
    sigma_t_ps = np.std(ps_final[:, 4]) / RF_FREQ * 1e9
    Q = 60e-12
    I_peak = peak_current(Q, sigma_t_ps, T)

    opt_currents = dict(felsim_currents)
    for idx, current in zip(S11_QUAD_INDICES, best_x):
        opt_currents[idx] = float(current)

    result_data = {
        'scenario': scenario_name,
        'h_chirp': h_chirp,
        'sigma_t_target_ps': sigma_t_target_ps,
        'weights': list(weights),
        'mse_twiss': float(mse_twiss),
        'transmission': float(T),
        'sigma_t_ps': float(sigma_t_ps),
        'I_peak_A': float(I_peak),
        'beta_x': float(bx), 'alpha_x': float(ax),
        'beta_y': float(by), 'alpha_y': float(ay),
        'dispersion': float(disp),
        'cost': float(best_result.fun),
        'nfev': total_nfev,
        'time_s': elapsed,
        'quad_currents': {str(k): float(v) for k, v in opt_currents.items()},
        'n_particles_in': int(beam_rft_cached.shape[0]),
        'n_particles_out': int(n_surviving),
    }

    _print(f"\n  ── Result ──")
    _print(f"  Cost = {best_result.fun:.4e}")
    _print(f"  MSE_Twiss = {mse_twiss:.4e}")
    _print(f"  Transmission = {T:.4f} ({n_surviving}/{beam_rft_cached.shape[0]})")
    _print(f"  σ_t = {sigma_t_ps:.3f} ps (target: {sigma_t_target_ps})")
    _print(f"  I_peak = {I_peak:.1f} A")
    _print(f"  β_x={bx:.4f}, α_x={ax:.4f}, β_y={by:.4f}, α_y={ay:.4f}")
    _print(f"  Disp = {disp:.6f}")
    _print(f"  Time = {elapsed:.1f} s, nfev = {total_nfev}")

    return result_data


# ═══════════════════════════════════════════════════════════════════════════════
#  3-code comparison
# ═══════════════════════════════════════════════════════════════════════════════

def run_3code_comparison(scenario_name, felsim_currents, throughput_currents,
                         h_chirp, sigma_t_ps_initial=2.0):
    """Compare FELsim-only vs throughput-optimized currents across 3 codes.

    Tracks beams through: FELsim (transfer matrices), COSY (apertures),
    RF-Track (apertures).
    """
    _print(f"\n{'='*72}")
    _print(f"  3-Code Comparison: {scenario_name}")
    _print(f"{'='*72}")

    Q = 60e-12
    rows = []

    for label, currents in [('FELsim-only', felsim_currents),
                             ('Throughput-opt', throughput_currents)]:
        _print(f"\n── {label} ──")

        # FELsim transfer matrix evaluation
        from beamline import lattice
        from excelElements import ExcelElements
        relat = lattice(1, fringeType=None)
        relat.setE(E=ENERGY)
        excel = ExcelElements(EXCEL_PATH)
        beamlineUH = excel.create_beamline()
        line = relat.changeBeamType("electron", ENERGY, beamlineUH)[:SEGMENTS]
        for idx, current in currents.items():
            if idx < len(line):
                line[idx].current = current

        beam_dist = create_felsim_beam(
            sigma_t_ps=sigma_t_ps_initial, sigma_E_pct=0.5,
            h_chirp=h_chirp, epsilon_n=8, N=N_PARTICLES, seed=SEED)

        particles = beam_dist.copy()
        for elem in line:
            particles = np.array(elem.useMatrice(particles))
        eb = beam()
        _, _, twiss = eb.cal_twiss(particles)
        bx_f = twiss.loc['x'][r"$\beta$ (m)"]
        ax_f = twiss.loc['x'][r"$\alpha$"]
        by_f = twiss.loc['y'][r"$\beta$ (m)"]
        ay_f = twiss.loc['y'][r"$\alpha$"]
        sigma_t_felsim = np.std(particles[:, 4]) / RF_FREQ * 1e9

        _print(f"  FELsim: β_x={bx_f:.4f}, α_x={ax_f:.4f}, "
               f"β_y={by_f:.4f}, α_y={ay_f:.4f}, σ_t={sigma_t_felsim:.3f} ps")

        row_base = {'label': label, 'scenario': scenario_name}

        rows.append({**row_base, 'code': 'FELsim',
                     'beta_x': bx_f, 'alpha_x': ax_f,
                     'beta_y': by_f, 'alpha_y': ay_f,
                     'sigma_t_ps': sigma_t_felsim,
                     'transmission': 1.0, 'n_out': len(particles),
                     'I_peak_A': peak_current(Q, sigma_t_felsim, 1.0)})

        # COSY tracking
        beam_cosy = create_cosy_beam(
            sigma_t_ps=sigma_t_ps_initial, sigma_E_pct=0.5,
            h_chirp=h_chirp, N=N_PARTICLES, seed=SEED)
        cosy_res = run_cosy_tracking(currents, beam_cosy,
                                     label=f'{label} COSY')
        if cosy_res is not None:
            sz = cosy_res['sigma_z_ps'] or sigma_t_ps_initial
            _print(f"  COSY: T={cosy_res['transmission']:.4f}, σ_z={sz:.3f} ps")
            rows.append({**row_base, 'code': 'COSY',
                         'beta_x': cosy_res['twiss'].get('beta_x', np.nan),
                         'alpha_x': cosy_res['twiss'].get('alpha_x', np.nan),
                         'beta_y': cosy_res['twiss'].get('beta_y', np.nan),
                         'alpha_y': cosy_res['twiss'].get('alpha_y', np.nan),
                         'sigma_t_ps': sz,
                         'transmission': cosy_res['transmission'],
                         'n_out': cosy_res['n_final'],
                         'I_peak_A': peak_current(Q, sz, cosy_res['transmission'])})

        # RF-Track tracking
        beam_rft = create_felsim_beam(
            sigma_t_ps=sigma_t_ps_initial, sigma_E_pct=0.5,
            h_chirp=h_chirp, epsilon_n=8, N=N_PARTICLES, seed=SEED)
        rft_res = run_rftrack_tracking(currents, beam_rft,
                                       label=f'{label} RF-Track')
        if rft_res is not None:
            sz = rft_res['sigma_z_ps'] or sigma_t_ps_initial
            _print(f"  RF-Track: T={rft_res['transmission']:.4f}, σ_z={sz:.3f} ps")
            rows.append({**row_base, 'code': 'RF-Track',
                         'beta_x': rft_res['twiss'].get('beta_x', np.nan),
                         'alpha_x': rft_res['twiss'].get('alpha_x', np.nan),
                         'beta_y': rft_res['twiss'].get('beta_y', np.nan),
                         'alpha_y': rft_res['twiss'].get('alpha_y', np.nan),
                         'sigma_t_ps': sz,
                         'transmission': rft_res['transmission'],
                         'n_out': rft_res['n_final'],
                         'I_peak_A': peak_current(Q, sz, rft_res['transmission'])})

    return rows


# ═══════════════════════════════════════════════════════════════════════════════
#  Plotting
# ═══════════════════════════════════════════════════════════════════════════════

def plot_comparison(comparison_rows, scenario_name, filename):
    """Bar chart: FELsim-only vs Throughput-opt across 3 codes."""
    codes = ['FELsim', 'COSY', 'RF-Track']
    labels_set = list(dict.fromkeys(r['label'] for r in comparison_rows))
    colors = {'FELsim-only': '#4477AA', 'Throughput-opt': '#EE6677'}

    metrics = [
        ('transmission', 'Transmission', '%', lambda v: v * 100),
        ('sigma_t_ps', 'σ_t (ps)', 'ps', lambda v: v),
        ('I_peak_A', 'I_peak (A)', 'A', lambda v: v),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for ax, (key, title, ylabel, transform) in zip(axes, metrics):
        x = np.arange(len(codes))
        width = 0.35

        for j, label in enumerate(labels_set):
            vals = []
            for code in codes:
                matching = [r for r in comparison_rows
                            if r['label'] == label and r['code'] == code]
                if matching:
                    vals.append(transform(matching[0][key]))
                else:
                    vals.append(0)
            offset = (j - 0.5) * width
            ax.bar(x + offset, vals, width, label=label,
                   color=colors.get(label, f'C{j}'), alpha=0.8, edgecolor='k')

        ax.set_xticks(x)
        ax.set_xticklabels(codes)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')

    fig.suptitle(f'W11: {scenario_name}', fontsize=14)
    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'{filename}.{ext}', dpi=150)
    plt.close(fig)
    _print(f"  Saved: {filename}.{{eps,png}}")


def plot_summary(all_results, filename):
    """Summary comparing both scenarios."""
    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    scenarios = list(all_results.keys())
    x = np.arange(len(scenarios))

    metrics = [
        ('mse_twiss', 'MSE Twiss', True),
        ('transmission', 'Transmission', False),
        ('sigma_t_ps', 'σ_t (ps)', False),
        ('I_peak_A', 'I_peak (A)', False),
    ]

    for ax, (key, ylabel, use_log) in zip(axes, metrics):
        vals = [all_results[s][key] for s in scenarios]
        ax.bar(x, vals, color=['#4477AA', '#EE6677'], alpha=0.8, edgecolor='k')
        ax.set_xticks(x)
        ax.set_xticklabels([s.replace('_', '\n') for s in scenarios], fontsize=9)
        ax.set_ylabel(ylabel)
        if use_log and all(v > 0 for v in vals):
            ax.set_yscale('log')
        ax.grid(True, alpha=0.3, axis='y')

    fig.suptitle('W11: Throughput Optimization Summary', fontsize=14)
    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'{filename}.{ext}', dpi=150)
    plt.close(fig)
    _print(f"  Saved: {filename}.{{eps,png}}")


# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='W11: Throughput Optimization')
    parser.add_argument('--scenario', choices=['2ps', '05ps', 'both'],
                        default='both', help='Scenario to run (default: both)')
    parser.add_argument('--particles', type=int, default=500,
                        help='Particles per run (default: 500)')
    parser.add_argument('--n-restarts', type=int, default=5,
                        help='NM restarts for Stage 11 (default: 5)')
    parser.add_argument('--chrom-bound', type=float, default=15,
                        help='Chromaticity quad upper bound (default: 15 A)')
    parser.add_argument('--currents', type=str, default=None,
                        help='Path to pre-computed FELsim currents JSON')
    parser.add_argument('--skip-comparison', action='store_true',
                        help='Skip 3-code comparison (faster)')
    args = parser.parse_args()

    if not _RFTRACK_AVAILABLE:
        _print("ERROR: RF-Track not available")
        sys.exit(1)

    OUTDIR.mkdir(parents=True, exist_ok=True)

    _print("W11: Throughput Optimization")
    _print(f"  E = {ENERGY} MeV, particles = {args.particles}")

    # Get FELsim-optimized currents (baseline for Stages 1-10)
    if args.currents:
        with open(args.currents) as f:
            data = json.load(f)
        felsim_currents = {int(k): float(v) for k, v in data.items()}
        _print(f"  Loaded FELsim currents from {args.currents}")
    else:
        _print("  Running FELsim 11-stage optimization...")
        t0 = time.perf_counter()
        res = run_optimization(
            bunch_spread=2.0, energy_std_percent=0.5, h=0,
            epsilon_n=8, nb_particles=args.particles, seed=SEED,
            n_restarts=args.n_restarts, chrom_upper_bound=args.chrom_bound)
        _print(f"  FELsim MSE = {res['mse']:.4e} ({time.perf_counter()-t0:.1f} s)")
        felsim_currents = {int(k): float(v) for k, v in res['quad_currents'].items()}
        with open(OUTDIR / 'currents_felsim.json', 'w') as f:
            json.dump({str(k): v for k, v in felsim_currents.items()}, f, indent=2)

    # Also run Twiss-only RF-Track Stage 11 optimization for comparison
    _print("\n  Running Twiss-only RF-Track Stage 11 optimization...")
    beta_xm, alpha_xm, beta_ym, alpha_ym, _ = compute_twiss_targets()
    targets = {
        'alpha_x': alpha_xm, 'alpha_y': alpha_ym,
        'beta_x': beta_xm, 'beta_y': beta_ym,
    }

    # Determine which scenarios to run
    scenario_keys = []
    if args.scenario in ('2ps', 'both'):
        scenario_keys.append('2ps_transport')
    if args.scenario in ('05ps', 'both'):
        scenario_keys.append('05ps_compress')

    all_results = {}
    all_comparison_rows = []

    for skey in scenario_keys:
        scen = SCENARIOS[skey]

        # Run throughput optimization
        result = run_scenario(
            scenario_name=skey,
            felsim_currents=felsim_currents,
            h_chirp=scen['h_chirp'],
            sigma_t_target_ps=scen['sigma_t_target_ps'],
            weights=scen['weights'],
            n_restarts=args.n_restarts,
            nb_particles=args.particles,
            seed=SEED,
            chrom_upper_bound=args.chrom_bound,
        )

        if result is None:
            _print(f"  Scenario {skey} FAILED")
            continue

        all_results[skey] = result

        # Save scenario results
        with open(OUTDIR / f'scenario_{skey}.json', 'w') as f:
            json.dump(result, f, indent=2)
        _print(f"  Saved: scenario_{skey}.json")

        # 3-code comparison
        if not args.skip_comparison:
            throughput_currents = {int(k): float(v)
                                   for k, v in result['quad_currents'].items()}
            comparison_rows = run_3code_comparison(
                skey, felsim_currents, throughput_currents,
                h_chirp=scen['h_chirp'])
            all_comparison_rows.extend(comparison_rows)

            plot_comparison(comparison_rows, skey,
                            f'W11_{skey}_comparison')

    # Summary comparison table
    if all_results:
        _print(f"\n{'='*72}")
        _print("  W11 Summary")
        _print(f"{'='*72}")

        _print(f"\n{'Scenario':<20s}  {'MSE':>10s}  {'T':>6s}  {'σ_t':>8s}  "
               f"{'I_pk':>8s}  {'Cost':>10s}  {'nfev':>6s}")
        _print("-" * 72)
        for skey, r in all_results.items():
            _print(f"{skey:<20s}  {r['mse_twiss']:10.4e}  {r['transmission']:6.3f}  "
                   f"{r['sigma_t_ps']:8.3f}  {r['I_peak_A']:8.1f}  "
                   f"{r['cost']:10.4e}  {r['nfev']:6d}")

        # Save comparison CSV
        if all_comparison_rows:
            csv_path = OUTDIR / 'W11_comparison.csv'
            header = ['scenario', 'label', 'code', 'beta_x', 'alpha_x',
                      'beta_y', 'alpha_y', 'sigma_t_ps', 'transmission',
                      'n_out', 'I_peak_A']
            with open(csv_path, 'w', newline='') as f:
                w = csv.DictWriter(f, fieldnames=header, extrasaction='ignore')
                w.writeheader()
                w.writerows(all_comparison_rows)
            _print(f"\n  Saved: W11_comparison.csv")

        if len(all_results) > 1:
            plot_summary(all_results, 'W11_throughput_summary')

    _print(f"\n{'='*72}")
    _print("  W11 Complete")
    _print(f"{'='*72}")


if __name__ == "__main__":
    main()
