"""C1C: Five-way comparison — MC-opt Stage 11 optimization via MultiCodeSimulator.

Extends C1B with MC-opt: Nelder-Mead optimization of Stage 11 quads using
MultiCodeSimulator (FELsim prefix + RF-Track suffix) as the forward model.

Methods:
  1. FELsim:  All elements in FELsim transfer matrices
  2. MC-val:  FELsim(0:87) → RF-Track(87:N), FELsim-optimised currents
  3. MC-opt:  FELsim(0:87) → RF-Track(87:N), NM-optimised Stage 11 currents
  4. RFT-val: Full RF-Track(0:N), FELsim-optimised currents
  5. RFT-opt: RF-Track prefix(0:87) cached + suffix(87:N) NM-optimised

MC-opt vs RFT-opt answers: do the optimal currents depend on which code
tracks the prefix? If they agree, the FELsim prefix is sufficient for
production optimization (faster, no prefix caching needed).

Usage:
    python -u test/C1C_multicode_optimization.py --smoke
    python -u test/C1C_multicode_optimization.py
    python -u test/C1C_multicode_optimization.py --emittance 8 14

Author: Eremey Valetov
"""

import sys
import time
import argparse
from pathlib import Path
import numpy as np
import scipy.optimize

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ebeam import beam
from multiCodeSimulator import MultiCodeSimulator, SimSection

from UHM_beamline_opt_05ps_params import (
    run_optimization, compute_twiss_targets,
    ENERGY, SEGMENTS, MSE_THRESHOLDS, write_csv,
)

try:
    import RF_Track
    _RFTRACK_AVAILABLE = True
except ImportError:
    _RFTRACK_AVAILABLE = False

from UHM_rftrack_opt import (
    S11_QUAD_INDICES, PREFIX_END, EXCEL_PATH,
    create_beam, evaluate_felsim_in_rftrack, run_rftrack_stage11,
    _compute_mse,
)
from C1B_hybrid_comparison import evaluate_multicode

JSON_PATH = Path(__file__).resolve().parent.parent.parent / 'var' / 'UH_FEL_beamline.json'
OUTDIR = Path(__file__).resolve().parent / 'results' / 'C1C'


def _print(msg):
    print(msg, flush=True)


# ── MC-opt: MultiCodeSimulator Stage 11 optimization ────────────────────

def _setup_multicode_pair(all_currents, space_charge=False, aperture=0.5):
    """Create mc_full and mc_disp sharing the same master beamline.

    mc_full: FELsim(0:87) + RF-Track(87:SEGMENTS) → undulator Twiss
    mc_disp: FELsim(0:93) → dispersion at element 92
    """
    rt_config = {}
    if space_charge:
        rt_config['space_charge'] = True
    if aperture != 0.5:
        rt_config['aperture'] = aperture

    mc_full = MultiCodeSimulator(
        sections=[
            SimSection("prefix", "felsim", (0, PREFIX_END)),
            SimSection("stage11", "rftrack", (PREFIX_END, SEGMENTS),
                       config=rt_config if rt_config else None),
        ],
        lattice_path=str(JSON_PATH),
        beam_energy=ENERGY,
    )

    # Apply quad currents
    for idx, current in all_currents.items():
        if idx < len(mc_full._master_beamline):
            mc_full._master_beamline[idx].current = current

    # Dispersion tracker: FELsim-only through elements 0:93
    # (element 87 is the chromaticity quad — a quad, not a dipole —
    # so FELsim handles it accurately)
    mc_disp = MultiCodeSimulator(
        sections=[
            SimSection("disp", "felsim", (0, 93)),
        ],
        lattice_path=str(JSON_PATH),
        beam_energy=ENERGY,
    )

    # Share master beamline so current mutations propagate to both
    mc_disp._master_beamline = mc_full._master_beamline

    return mc_full, mc_disp


def _mc_mse(mc_full, mc_disp, beam_dist, targets, ebeam_obj):
    """Compute MSE via MultiCodeSimulator (two trackings)."""
    result = mc_full.simulate(particles=beam_dist)
    if not result.success or result.final_particles is None:
        return 1e6, None
    ps_final = result.final_particles
    if ps_final.shape[0] < 10:
        return 1e6, None

    _, _, twiss_f = ebeam_obj.cal_twiss(ps_final)
    bx = twiss_f.loc['x'][r"$\beta$ (m)"]
    ax = twiss_f.loc['x'][r"$\alpha$"]
    by = twiss_f.loc['y'][r"$\beta$ (m)"]
    ay = twiss_f.loc['y'][r"$\alpha$"]

    result_93 = mc_disp.simulate(particles=beam_dist)
    if not result_93.success or result_93.final_particles is None:
        disp = 0.0
    elif result_93.final_particles.shape[0] < 10:
        disp = 0.0
    else:
        _, _, twiss_93 = ebeam_obj.cal_twiss(result_93.final_particles)
        disp = twiss_93.loc['x'][r"$D$ (m)"]

    mse = _compute_mse(bx, ax, by, ay, disp, targets)
    return mse, {
        'mse': mse, 'beta_x': bx, 'alpha_x': ax,
        'beta_y': by, 'alpha_y': ay, 'disp': disp,
        'ngood': ps_final.shape[0],
    }


def _mc_nm_objective(x, mc_full, mc_disp, beam_dist, targets, ebeam_obj, counter):
    """NM objective: mutate Stage 11 quads → mc_mse."""
    counter[0] += 1
    for idx, current in zip(S11_QUAD_INDICES, x):
        mc_full._master_beamline[idx].current = current
    # mc_disp shares _master_beamline, picks up changes automatically
    mse, _ = _mc_mse(mc_full, mc_disp, beam_dist, targets, ebeam_obj)
    return mse


def run_mc_stage11(felsim_result, beam_dist, targets,
                   n_restarts=5, chrom_upper_bound=15,
                   space_charge=False, aperture=0.5):
    """Optimize Stage 11 via MultiCodeSimulator, warm-started from FELsim."""
    qb = 10
    cb = chrom_upper_bound
    ebeam_obj = beam()

    all_currents = dict(felsim_result['quad_currents'])
    mc_full, mc_disp = _setup_multicode_pair(
        all_currents, space_charge=space_charge, aperture=aperture)

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
            _mc_nm_objective, x0, method='Nelder-Mead',
            bounds=bounds,
            args=(mc_full, mc_disp, beam_dist, targets, ebeam_obj, counter),
        )
        total_nfev += counter[0]
        _print(f"    Restart {i+1}/{len(starts)}: MSE = {result.fun:.4e} "
               f"(nfev={counter[0]})")

        if best_result is None or result.fun < best_result.fun:
            best_result = result
            best_x = result.x.copy()

    elapsed = time.perf_counter() - t0

    # Final evaluation
    for idx, current in zip(S11_QUAD_INDICES, best_x):
        mc_full._master_beamline[idx].current = current
    _, details = _mc_mse(mc_full, mc_disp, beam_dist, targets, ebeam_obj)
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


# ── Comparison loop ──────────────────────────────────────────────────────

def run_comparison(emittance_points, nb_particles=500, seed=42,
                   n_restarts=5, chrom_upper_bound=15,
                   space_charge=False, aperture=0.5):
    """Five-way comparison: FELsim / MC-val / MC-opt / RFT-val / RFT-opt."""
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

        # ── 1. FELsim ────────────────────────────────────────────────
        _print(f"  [1/5] FELsim optimization...")
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
        _print(f"    MSE = {felsim['mse']:.4e} ({felsim_time:.1f} s)")

        # ── 2. MC-val ────────────────────────────────────────────────
        _print(f"  [2/5] MultiCode validation (FELsim currents)...")
        t0 = time.perf_counter()
        mc_val = evaluate_multicode(felsim, beam_dist, targets,
                                    space_charge=space_charge, aperture=aperture)
        mc_val_time = time.perf_counter() - t0

        if mc_val:
            rows.append([
                en, 'MC-val', mc_val['mse'],
                mc_val['alpha_x'], mc_val['alpha_y'],
                mc_val['beta_x'], mc_val['beta_y'], mc_val['disp'],
                fc[87], fc[93], fc[95], fc[97],
                mc_val_time, 0, mc_val['ngood'],
            ])
            _print(f"    MSE = {mc_val['mse']:.4e} ({mc_val_time:.1f} s)")
        else:
            rows.append([en, 'MC-val'] + [float('nan')] * (len(header) - 2))
            _print(f"    FAILED")

        # ── 3. MC-opt ────────────────────────────────────────────────
        _print(f"  [3/5] MultiCode optimization ({n_restarts} restarts)...")
        mc_opt = run_mc_stage11(
            felsim, beam_dist, targets,
            n_restarts=n_restarts, chrom_upper_bound=chrom_upper_bound,
            space_charge=space_charge, aperture=aperture)

        if mc_opt:
            oc = mc_opt['quad_currents']
            rows.append([
                en, 'MC-opt', mc_opt['mse'],
                mc_opt['alpha_x'], mc_opt['alpha_y'],
                mc_opt['beta_x'], mc_opt['beta_y'], mc_opt['disp'],
                oc[87], oc[93], oc[95], oc[97],
                mc_opt['time_s'], mc_opt['nfev'], mc_opt['ngood'],
            ])
            _print(f"    MSE = {mc_opt['mse']:.4e} "
                   f"({mc_opt['time_s']:.1f} s, nfev={mc_opt['nfev']})")
        else:
            rows.append([en, 'MC-opt'] + [float('nan')] * (len(header) - 2))
            _print(f"    FAILED")

        # ── 4. RFT-val ──────────────────────────────────────────────
        _print(f"  [4/5] RF-Track validation (FELsim currents)...")
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
            _print(f"    MSE = {rft_val['mse']:.4e} ({val_time:.1f} s)")
        else:
            rows.append([en, 'RFT-val'] + [float('nan')] * (len(header) - 2))
            _print(f"    FAILED")

        # ── 5. RFT-opt ──────────────────────────────────────────────
        _print(f"  [5/5] RF-Track optimization ({n_restarts} restarts)...")
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
            _print(f"    MSE = {rft_opt['mse']:.4e} "
                   f"({rft_opt['time_s']:.1f} s, nfev={rft_opt['nfev']})")
        else:
            rows.append([en, 'RFT-opt'] + [float('nan')] * (len(header) - 2))
            _print(f"    FAILED")

    # Save
    csv_path = OUTDIR / 'comparison.csv'
    write_csv(csv_path, header, rows)
    _print(f"\nSaved {csv_path}")

    _print_summary(rows, targets)
    _plot_comparison(rows, emittance_points)

    return rows


# ── Output ───────────────────────────────────────────────────────────────

def _print_summary(rows, targets):
    _print(f"\n{'─'*130}")
    _print(f"{'ε_n':>5} {'Method':>10} {'MSE':>12} {'β_x':>8} {'β_y':>8} "
           f"{'α_x':>8} {'α_y':>10} {'Disp':>8} "
           f"{'I_87':>6} {'I_93':>6} {'I_95':>6} {'I_97':>6} "
           f"{'Time':>7} {'nfev':>6}")
    _print(f"{'':>5} {'Target':>10} {'':>12} "
           f"{targets['beta_x']:>8.4f} {targets['beta_y']:>8.4f} "
           f"{targets['alpha_x']:>8.4f} {targets['alpha_y']:>10.6f}")
    _print(f"{'─'*130}")

    def _fmt(v, w, f):
        try:
            if np.isnan(v):
                return ' ' * (w - 3) + 'N/A'
        except (TypeError, ValueError):
            pass
        return f.format(v)

    for r in rows:
        _print(f"{r[0]:>5.0f} {r[1]:>10} "
               f"{_fmt(r[2], 12, '{:12.4e}')} "
               f"{_fmt(r[5], 8, '{:8.4f}')} {_fmt(r[6], 8, '{:8.4f}')} "
               f"{_fmt(r[3], 8, '{:8.4f}')} {_fmt(r[4], 10, '{:10.6f}')} "
               f"{_fmt(r[7], 8, '{:8.4f}')} "
               f"{_fmt(r[8], 6, '{:6.2f}')} {_fmt(r[9], 6, '{:6.2f}')} "
               f"{_fmt(r[10], 6, '{:6.2f}')} {_fmt(r[11], 6, '{:6.2f}')} "
               f"{_fmt(r[12], 7, '{:7.1f}')} {_fmt(r[13], 6, '{:6.0f}')}")
    _print(f"{'─'*130}")


def _plot_comparison(rows, emittance_points):
    methods = ['FELsim', 'MC-val', 'MC-opt', 'RFT-val', 'RFT-opt']
    colors = {
        'FELsim': '#4477AA', 'MC-val': '#228833', 'MC-opt': '#66CCEE',
        'RFT-val': '#CCBB44', 'RFT-opt': '#EE6677',
    }
    labels = {
        'FELsim': 'FELsim (transfer matrix)',
        'MC-val': 'MC validation (FELsim currents)',
        'MC-opt': 'MC optimised (NM)',
        'RFT-val': 'RF-Track (FELsim currents)',
        'RFT-opt': 'RF-Track optimised (NM)',
    }

    x = np.arange(len(emittance_points))
    width = 0.16

    fig, ax = plt.subplots(figsize=(12, 5))
    for j, m in enumerate(methods):
        m_rows = [r for r in rows if r[1] == m]
        mse_vals = [r[2] for r in m_rows]
        offset = (j - 2) * width
        ax.bar(x + offset, mse_vals, width, label=labels[m], color=colors[m])

    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(e)) for e in emittance_points])
    ax.set_xlabel(r'Normalised emittance $\varepsilon_n$ ($\pi\cdot$mm$\cdot$mrad)')
    ax.set_ylabel('MSE')
    ax.set_title('Stage 11 MSE: Five-Way Comparison (MC-opt vs RFT-opt)')
    for tl, thresh in MSE_THRESHOLDS.items():
        clr = {'Excellent': 'green', 'Acceptable': 'orange', 'Marginal': 'red'}
        ax.axhline(thresh, color=clr[tl], linestyle='--', alpha=0.6,
                    label=f'{tl} ({thresh})')
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3, axis='y')
    fig.tight_layout()
    fig.savefig(OUTDIR / 'mse_comparison.pdf')
    plt.close(fig)
    _print(f"  Saved {OUTDIR / 'mse_comparison.pdf'}")


# ── CLI ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='C1C: Five-way comparison with MC-opt Stage 11 optimization')
    parser.add_argument('--emittance', type=float, nargs='+', default=[5, 8, 14])
    parser.add_argument('--particles', type=int, default=500)
    parser.add_argument('--n-restarts', type=int, default=5)
    parser.add_argument('--chrom-bound', type=float, default=15)
    parser.add_argument('--space-charge', action='store_true')
    parser.add_argument('--aperture', type=float, default=0.5)
    parser.add_argument('--smoke', action='store_true',
                        help='Quick smoke test: ε_n=8, 1 restart')
    args = parser.parse_args()

    if not _RFTRACK_AVAILABLE:
        _print("ERROR: RF-Track not available")
        sys.exit(1)

    _print("C1C: Five-Way Comparison (MC-opt vs RFT-opt)")
    _print(f"  E = {ENERGY} MeV, particles = {args.particles}, "
           f"SC = {'ON' if args.space_charge else 'OFF'}")

    if args.smoke:
        _print("  Mode: smoke test (ε_n=8, 1 restart)")
        run_comparison([8], nb_particles=args.particles, n_restarts=1,
                       chrom_upper_bound=args.chrom_bound,
                       space_charge=args.space_charge, aperture=args.aperture)
    else:
        _print(f"  Mode: full comparison (ε_n={args.emittance}, "
               f"{args.n_restarts} restarts)")
        run_comparison(args.emittance, nb_particles=args.particles,
                       n_restarts=args.n_restarts,
                       chrom_upper_bound=args.chrom_bound,
                       space_charge=args.space_charge, aperture=args.aperture)


if __name__ == "__main__":
    main()
