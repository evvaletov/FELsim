"""RF-Track Stage 11 Optimization: Hybrid FELsim→RF-Track undulator matching.

Runs FELsim stages 1-10 (linear transfer matrices), then optimizes Stage 11
using RF-Track particle tracking. Compares with pure FELsim optimization to
assess whether nonlinear effects change optimal quad currents.

Performance: elements 0:87 are static during Stage 11 optimization, so the
beam is pre-tracked through them once and cached. Each NM evaluation only
tracks elements 87:118 (31 elements) and 87:93 (6 elements).

Usage:
    python -u test/UHM_rftrack_opt.py --smoke          # quick smoke test
    python -u test/UHM_rftrack_opt.py                   # full comparison
    python -u test/UHM_rftrack_opt.py --emittance 8 14  # custom points

Author: Eremey Valetov
"""

import sys
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
    run_optimization, compute_twiss_targets, QUAD_INDICES, BASELINE,
    ENERGY, RF_FREQ, SEGMENTS, MSE_THRESHOLDS, write_csv,
)

try:
    import RF_Track as rft
    from rftrackAdapter import RFTrackAdapter
    _RFTRACK_AVAILABLE = True
except ImportError:
    _RFTRACK_AVAILABLE = False

# ── Constants ─────────────────────────────────────────────────────────────────

S11_QUAD_INDICES = [87, 93, 95, 97]
# Elements 0:87 are static during Stage 11 optimization
PREFIX_END = 87

EXCEL_PATH = (Path(__file__).resolve().parent.parent.parent
              / 'beam_excel' / 'Beamline_elements.xlsx')
OUTDIR = Path(__file__).resolve().parent / 'results' / 'rftrack_opt'


def _print(msg):
    """Print with immediate flush (for redirected output)."""
    print(msg, flush=True)


# ── Beam generation ──────────────────────────────────────────────────────────

def create_beam(epsilon_n=8, bunch_spread=0.5, energy_std_percent=0.5,
                h=5e9, x_std=0.8, y_std=0.8, nb_particles=500, seed=42):
    """Create FELsim-coordinate beam distribution (same RNG sequence as run_optimization)."""
    _, _, _, _, relat = compute_twiss_targets()
    norm = relat.gamma * relat.beta
    epsilon = epsilon_n / norm

    np.random.seed(seed)
    beam_dist = beam().gen_6d_gaussian(
        0, [x_std, epsilon / x_std, y_std, epsilon / y_std,
            bunch_spread * 1e-9 * RF_FREQ, energy_std_percent * 10],
        nb_particles)
    beam_dist[:, 5] += h * beam_dist[:, 4] / RF_FREQ
    return beam_dist


# ── RF-Track adapter setup ───────────────────────────────────────────────────

def setup_rftrack_adapter(currents, space_charge=False, aperture=0.5):
    """Create RF-Track adapter with quad currents applied (118 elements)."""
    sim = RFTrackAdapter(
        lattice_path=str(EXCEL_PATH),
        beam_energy=ENERGY,
        space_charge=space_charge,
        aperture=aperture,
    )
    sim.beamline = sim.beamline[:SEGMENTS]
    for idx, current in currents.items():
        if idx < len(sim.beamline):
            sim._modify_element(idx, current=current)
    sim._build_lattice()
    return sim


# ── Prefix caching ──────────────────────────────────────────────────────────

def pretrack_prefix(sim, beam_rft, n_prefix=PREFIX_END):
    """Track beam through static prefix (elements 0:87). Returns RF-Track-coord beam state.

    Elements 0:86 don't change during Stage 11 optimization, so this
    needs to be done only once per emittance point.
    """
    ps = sim.track_elements(beam_rft, 0, n_prefix)
    if ps.ndim != 2 or ps.shape[0] < 10:
        return None
    return ps


def _track_suffix(sim, beam_rft_cached, start_idx, end_idx):
    """Track through beamline[start_idx:end_idx] with analytical dipole corrections.

    Returns FELsim-coord particles or None on failure.
    """
    ps = sim.track_elements(beam_rft_cached, start_idx, end_idx)
    if ps.ndim != 2 or ps.shape[0] < 10:
        return None
    return sim.transform_coordinates(
        ps, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM
    )


# ── MSE computation ──────────────────────────────────────────────────────────

def _compute_mse(bx, ax, by, ay, disp, targets):
    """MSE formula matching beamOptimizer Stage 11 (5 objectives)."""
    return (0.5 * disp**2
            + (ax - targets['alpha_x'])**2
            + (ay - targets['alpha_y'])**2
            + (bx - targets['beta_x'])**2
            + (by - targets['beta_y'])**2) / 5


def rftrack_mse_cached(sim, beam_rft_cached, targets, ebeam_obj):
    """Compute Stage 11 MSE using cached prefix beam state.

    Two suffix trackings per call:
      1. Elements 87:118 → Twiss at undulator entrance
      2. Elements 87:93  → dispersion at element 92

    Returns (mse, details_dict) or (1e6, None) on failure.
    """
    # Suffix full (87:118) → final Twiss
    ps_final = _track_suffix(sim, beam_rft_cached, PREFIX_END, SEGMENTS)
    if ps_final is None:
        return 1e6, None

    _, _, twiss_f = ebeam_obj.cal_twiss(ps_final)
    bx = twiss_f.loc['x'][r"$\beta$ (m)"]
    ax = twiss_f.loc['x'][r"$\alpha$"]
    by = twiss_f.loc['y'][r"$\beta$ (m)"]
    ay = twiss_f.loc['y'][r"$\alpha$"]

    # Suffix dispersion (87:93) → dispersion at element 92
    ps_92 = _track_suffix(sim, beam_rft_cached, PREFIX_END, 93)
    if ps_92 is None:
        return 1e6, None

    _, _, twiss_92 = ebeam_obj.cal_twiss(ps_92)
    disp = twiss_92.loc['x'][r"$D$ (m)"]

    mse = _compute_mse(bx, ax, by, ay, disp, targets)

    return mse, {
        'mse': mse,
        'beta_x': bx, 'alpha_x': ax,
        'beta_y': by, 'alpha_y': ay,
        'disp': disp, 'ngood': ps_final.shape[0],
    }


def rftrack_mse_full(sim, beam_rft, targets, ebeam_obj):
    """Compute Stage 11 MSE with full lattice tracking (for validation)."""
    # Full lattice (0:SEGMENTS) with analytical dipole corrections
    ps = sim.track_elements(beam_rft, 0, SEGMENTS)
    if ps.ndim != 2 or ps.shape[0] < 10:
        return 1e6, None
    ngood = ps.shape[0]

    ps_felsim = sim.transform_coordinates(
        ps, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM
    )
    _, _, twiss_f = ebeam_obj.cal_twiss(ps_felsim)
    bx = twiss_f.loc['x'][r"$\beta$ (m)"]
    ax = twiss_f.loc['x'][r"$\alpha$"]
    by = twiss_f.loc['y'][r"$\beta$ (m)"]
    ay = twiss_f.loc['y'][r"$\alpha$"]

    # Sub-lattice (0:93) for dispersion
    ps2 = sim.track_elements(beam_rft, 0, 93)
    if ps2.ndim != 2 or ps2.shape[0] < 10:
        return 1e6, None

    ps_felsim2 = sim.transform_coordinates(
        ps2, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM
    )
    _, _, twiss_92 = ebeam_obj.cal_twiss(ps_felsim2)
    disp = twiss_92.loc['x'][r"$D$ (m)"]

    mse = _compute_mse(bx, ax, by, ay, disp, targets)

    return mse, {
        'mse': mse,
        'beta_x': bx, 'alpha_x': ax,
        'beta_y': by, 'alpha_y': ay,
        'disp': disp, 'ngood': ngood,
    }


# ── NM objective ─────────────────────────────────────────────────────────────

def _nm_objective(x, sim, beam_rft_cached, targets, ebeam_obj, counter):
    """NM objective: set Stage 11 quads → track suffix → MSE."""
    counter[0] += 1
    for idx, current in zip(S11_QUAD_INDICES, x):
        sim._modify_element(idx, current=current)
    # No _build_lattice() needed — we build suffix lattices on the fly
    mse, _ = rftrack_mse_cached(sim, beam_rft_cached, targets, ebeam_obj)
    return mse


# ── Stage 11 optimization ───────────────────────────────────────────────────

def run_rftrack_stage11(felsim_result, beam_dist, targets,
                        n_restarts=5, chrom_upper_bound=15,
                        space_charge=False, aperture=0.5):
    """Optimize Stage 11 with RF-Track, warm-started from FELsim solution."""
    qb = 10
    cb = chrom_upper_bound
    ebeam_obj = beam()

    all_currents = dict(felsim_result['quad_currents'])
    sim = setup_rftrack_adapter(all_currents, space_charge=space_charge,
                                aperture=aperture)
    beam_rft = sim.transform_coordinates(
        beam_dist, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK
    )

    # Pre-track static prefix (elements 0:87) — done once
    _print(f"    Pre-tracking prefix (elements 0:{PREFIX_END})...")
    t_prefix = time.perf_counter()
    beam_rft_cached = pretrack_prefix(sim, beam_rft)
    t_prefix = time.perf_counter() - t_prefix
    if beam_rft_cached is None:
        _print(f"    FAILED: prefix tracking lost too many particles")
        return None
    _print(f"    Prefix tracked in {t_prefix:.1f} s")

    # Multi-start: FELsim solution (warm) + random restarts
    felsim_s11 = [all_currents[i] for i in S11_QUAD_INDICES]
    rng = np.random.RandomState(42 + 999)
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
            _nm_objective, x0, method='Nelder-Mead',
            bounds=bounds,
            args=(sim, beam_rft_cached, targets, ebeam_obj, counter),
        )
        total_nfev += counter[0]
        _print(f"    Restart {i+1}/{len(starts)}: MSE = {result.fun:.4e} "
               f"(nfev={counter[0]})")

        if best_result is None or result.fun < best_result.fun:
            best_result = result
            best_x = result.x.copy()

    elapsed = time.perf_counter() - t0

    # Final evaluation with best solution
    for idx, current in zip(S11_QUAD_INDICES, best_x):
        sim._modify_element(idx, current=current)
    _, details = rftrack_mse_cached(sim, beam_rft_cached, targets, ebeam_obj)
    if details is None:
        return None

    opt_currents = dict(all_currents)
    for idx, current in zip(S11_QUAD_INDICES, best_x):
        opt_currents[idx] = float(current)
    details['quad_currents'] = opt_currents
    details['nfev'] = total_nfev
    details['time_s'] = elapsed
    details['converged'] = best_result.success

    return details


# ── Validation evaluation ────────────────────────────────────────────────────

def evaluate_felsim_in_rftrack(felsim_result, beam_dist, targets,
                                space_charge=False, aperture=0.5):
    """Evaluate FELsim-optimized currents in RF-Track (full lattice, no optimization)."""
    ebeam_obj = beam()
    sim = setup_rftrack_adapter(felsim_result['quad_currents'],
                                space_charge=space_charge, aperture=aperture)
    beam_rft = sim.transform_coordinates(
        beam_dist, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK
    )
    _, details = rftrack_mse_full(sim, beam_rft, targets, ebeam_obj)
    if details:
        details['quad_currents'] = dict(felsim_result['quad_currents'])
    return details


# ── Main comparison loop ─────────────────────────────────────────────────────

def run_comparison(emittance_points, nb_particles=500, seed=42,
                   n_restarts=5, chrom_upper_bound=15,
                   space_charge=False, aperture=0.5):
    """Run three-way comparison: FELsim / RF-Track validation / RF-Track optimized."""
    beta_xm, alpha_xm, beta_ym, alpha_ym, _ = compute_twiss_targets()
    targets = {
        'alpha_x': alpha_xm, 'alpha_y': alpha_ym,
        'beta_x': beta_xm, 'beta_y': beta_ym,
    }

    OUTDIR.mkdir(parents=True, exist_ok=True)

    header = ['epsilon_n', 'method', 'mse', 'alpha_x', 'alpha_y',
              'beta_x', 'beta_y', 'disp',
              'quad_87', 'quad_93', 'quad_95', 'quad_97',
              'time_s', 'nfev', 'ngood']
    rows = []

    for i, en in enumerate(emittance_points):
        _print(f"\n{'='*70}")
        _print(f"  [{i+1}/{len(emittance_points)}] ε_n = {en}")
        _print(f"{'='*70}")

        beam_dist = create_beam(epsilon_n=en, nb_particles=nb_particles, seed=seed)

        # ── FELsim (full 11 stages) ──────────────────────────────────────
        _print(f"  Running FELsim optimization (11 stages)...")
        t0 = time.perf_counter()
        felsim = run_optimization(
            epsilon_n=en, nb_particles=nb_particles, seed=seed,
            chrom_upper_bound=chrom_upper_bound, n_restarts=n_restarts)
        felsim_time = time.perf_counter() - t0

        fc = felsim['quad_currents']
        rows.append([
            en, 'FELsim', felsim['mse'],
            felsim['alpha_x'], felsim['alpha_y'],
            felsim['beta_x'], felsim['beta_y'], felsim['disp_resid'],
            fc[87], fc[93], fc[95], fc[97],
            felsim_time, felsim['nfev'], nb_particles,
        ])
        _print(f"    FELsim: MSE = {felsim['mse']:.4e} ({felsim_time:.1f} s)")

        # ── RF-Track validation (FELsim currents) ────────────────────────
        _print(f"  Evaluating FELsim currents in RF-Track...")
        t0 = time.perf_counter()
        rft_val = evaluate_felsim_in_rftrack(
            felsim, beam_dist, targets,
            space_charge=space_charge, aperture=aperture)
        val_time = time.perf_counter() - t0

        if rft_val:
            rows.append([
                en, 'RFT-val', rft_val['mse'],
                rft_val['alpha_x'], rft_val['alpha_y'],
                rft_val['beta_x'], rft_val['beta_y'], rft_val['disp'],
                fc[87], fc[93], fc[95], fc[97],
                val_time, 0, rft_val['ngood'],
            ])
            _print(f"    RFT-val: MSE = {rft_val['mse']:.4e} "
                   f"(ngood={rft_val['ngood']}, {val_time:.1f} s)")
        else:
            rows.append([en, 'RFT-val'] + [float('nan')] * (len(header) - 2))
            _print(f"    RFT-val: FAILED")

        # ── RF-Track Stage 11 optimization ───────────────────────────────
        _print(f"  Running RF-Track Stage 11 optimization "
               f"({n_restarts} restarts)...")
        rft_opt = run_rftrack_stage11(
            felsim, beam_dist, targets,
            n_restarts=n_restarts, chrom_upper_bound=chrom_upper_bound,
            space_charge=space_charge, aperture=aperture)

        if rft_opt:
            oc = rft_opt['quad_currents']
            rows.append([
                en, 'RFT-opt', rft_opt['mse'],
                rft_opt['alpha_x'], rft_opt['alpha_y'],
                rft_opt['beta_x'], rft_opt['beta_y'], rft_opt['disp'],
                oc[87], oc[93], oc[95], oc[97],
                rft_opt['time_s'], rft_opt['nfev'], rft_opt['ngood'],
            ])
            _print(f"    RFT-opt: MSE = {rft_opt['mse']:.4e} "
                   f"({rft_opt['time_s']:.1f} s, nfev={rft_opt['nfev']})")
        else:
            rows.append([en, 'RFT-opt'] + [float('nan')] * (len(header) - 2))
            _print(f"    RFT-opt: FAILED")

    # Save CSV
    csv_path = OUTDIR / 'comparison.csv'
    write_csv(csv_path, header, rows)
    _print(f"\nSaved {csv_path}")

    # Summary table
    _print_summary(rows, targets)

    # Plots
    _plot_comparison(rows, emittance_points)

    return rows


# ── Output ───────────────────────────────────────────────────────────────────

def _print_summary(rows, targets):
    """Print formatted comparison table."""
    _print(f"\n{'─'*120}")
    _print(f"{'ε_n':>5} {'Method':>10} {'MSE':>12} {'β_x':>8} {'β_y':>8} "
           f"{'α_x':>8} {'α_y':>10} {'Disp':>8} "
           f"{'I_87':>6} {'I_93':>6} {'I_95':>6} {'I_97':>6} {'Time':>7}")
    _print(f"{'':>5} {'Target':>10} {'':>12} "
           f"{targets['beta_x']:>8.4f} {targets['beta_y']:>8.4f} "
           f"{targets['alpha_x']:>8.4f} {targets['alpha_y']:>10.6f}")
    _print(f"{'─'*120}")
    for r in rows:
        def _fmt(v, w, f):
            try:
                if np.isnan(v):
                    return ' ' * (w - 3) + 'N/A'
            except (TypeError, ValueError):
                pass
            return f.format(v)

        _print(f"{r[0]:>5.0f} {r[1]:>10} "
               f"{_fmt(r[2], 12, '{:12.4e}')} "
               f"{_fmt(r[5], 8, '{:8.4f}')} {_fmt(r[6], 8, '{:8.4f}')} "
               f"{_fmt(r[3], 8, '{:8.4f}')} {_fmt(r[4], 10, '{:10.6f}')} "
               f"{_fmt(r[7], 8, '{:8.4f}')} "
               f"{_fmt(r[8], 6, '{:6.2f}')} {_fmt(r[9], 6, '{:6.2f}')} "
               f"{_fmt(r[10], 6, '{:6.2f}')} {_fmt(r[11], 6, '{:6.2f}')} "
               f"{_fmt(r[12], 7, '{:7.1f}')}")
    _print(f"{'─'*120}")


def _plot_comparison(rows, emittance_points):
    """Generate MSE comparison bar chart."""
    methods = ['FELsim', 'RFT-val', 'RFT-opt']
    colors = {'FELsim': '#4477AA', 'RFT-val': '#CCBB44', 'RFT-opt': '#EE6677'}
    labels = {
        'FELsim': 'FELsim (transfer matrix)',
        'RFT-val': 'RF-Track (FELsim currents)',
        'RFT-opt': 'RF-Track (optimised)',
    }

    x = np.arange(len(emittance_points))
    width = 0.25

    fig, ax = plt.subplots(figsize=(10, 5))
    for j, m in enumerate(methods):
        m_rows = [r for r in rows if r[1] == m]
        mse_vals = [r[2] for r in m_rows]
        offset = (j - 1) * width
        ax.bar(x + offset, mse_vals, width, label=labels[m], color=colors[m])

    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(e)) for e in emittance_points])
    ax.set_xlabel(r'Normalised emittance $\varepsilon_n$ ($\pi\cdot$mm$\cdot$mrad)')
    ax.set_ylabel('MSE')
    ax.set_title('Stage 11 MSE: FELsim vs RF-Track')
    for tl, thresh in MSE_THRESHOLDS.items():
        clr = {'Excellent': 'green', 'Acceptable': 'orange', 'Marginal': 'red'}
        ax.axhline(thresh, color=clr[tl], linestyle='--', alpha=0.6,
                    label=f'{tl} ({thresh})')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')
    fig.tight_layout()
    fig.savefig(OUTDIR / 'mse_comparison.eps', format='eps')
    fig.savefig(OUTDIR / 'mse_comparison.pdf')
    plt.close(fig)
    _print(f"  Saved {OUTDIR / 'mse_comparison.eps'}")

    _plot_twiss_comparison(rows, emittance_points)


def _plot_twiss_comparison(rows, emittance_points):
    """Plot Twiss parameter comparison (FELsim vs RF-Track optimised)."""
    beta_xm, alpha_xm, beta_ym, alpha_ym, _ = compute_twiss_targets()

    methods = ['FELsim', 'RFT-opt']
    colors = {'FELsim': '#4477AA', 'RFT-opt': '#EE6677'}

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    panels = [
        ('beta_x', r'$\beta_x$ (m)', beta_xm, 5),
        ('beta_y', r'$\beta_y$ (m)', beta_ym, 6),
        ('alpha_x', r'$\alpha_x$', alpha_xm, 3),
        ('alpha_y', r'$\alpha_y$', alpha_ym, 4),
    ]

    for ax, (key, ylabel, target, col_idx) in zip(axes.flat, panels):
        for m in methods:
            m_rows = [r for r in rows if r[1] == m]
            eps = [r[0] for r in m_rows]
            vals = [r[col_idx] for r in m_rows]
            ax.plot(eps, vals, 'o-', color=colors[m], markersize=5, label=m)
        ax.axhline(target, color='black', linestyle='--', alpha=0.5,
                    label=f'Target = {target:.4f}')
        ax.set_xlabel(r'$\varepsilon_n$ ($\pi\cdot$mm$\cdot$mrad)')
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.suptitle('Twiss at Undulator Entrance: FELsim vs RF-Track', fontsize=14)
    fig.tight_layout()
    fig.savefig(OUTDIR / 'twiss_comparison.eps', format='eps')
    fig.savefig(OUTDIR / 'twiss_comparison.pdf')
    plt.close(fig)
    _print(f"  Saved {OUTDIR / 'twiss_comparison.eps'}")


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='RF-Track Stage 11 Optimisation — Hybrid FELsim/RF-Track')
    parser.add_argument('--emittance', type=float, nargs='+', default=[5, 8, 14],
                        help='Emittance points (default: 5 8 14)')
    parser.add_argument('--particles', type=int, default=500,
                        help='Particles per run (default: 500)')
    parser.add_argument('--n-restarts', type=int, default=5,
                        help='NM restarts for Stage 11 (default: 5)')
    parser.add_argument('--chrom-bound', type=float, default=15,
                        help='Chromaticity quad upper bound (default: 15 A)')
    parser.add_argument('--space-charge', action='store_true',
                        help='Enable space charge in RF-Track')
    parser.add_argument('--aperture', type=float, default=0.5,
                        help='RF-Track aperture in metres (default: 0.5)')
    parser.add_argument('--smoke', action='store_true',
                        help='Quick smoke test: ε_n=8, 1 restart')
    args = parser.parse_args()

    if not _RFTRACK_AVAILABLE:
        _print("ERROR: RF-Track not available. Install with: pip install RF-Track")
        sys.exit(1)

    _print("RF-Track Stage 11 Optimisation")
    _print(f"  E = {ENERGY} MeV, particles = {args.particles}, "
           f"SC = {'ON' if args.space_charge else 'OFF'}, "
           f"aperture = {args.aperture} m")

    if args.smoke:
        _print("  Mode: smoke test (ε_n=8, 1 restart)")
        run_comparison(
            [8], nb_particles=args.particles, n_restarts=1,
            chrom_upper_bound=args.chrom_bound,
            space_charge=args.space_charge, aperture=args.aperture)
    else:
        _print(f"  Mode: full comparison (ε_n={args.emittance}, "
               f"{args.n_restarts} restarts)")
        run_comparison(
            args.emittance, nb_particles=args.particles,
            n_restarts=args.n_restarts,
            chrom_upper_bound=args.chrom_bound,
            space_charge=args.space_charge, aperture=args.aperture)


if __name__ == "__main__":
    main()
