"""C1B: Four-way hybrid comparison — FELsim / MC-val / RFT-val / RFT-opt.

Extends C1 (UHM_rftrack_opt.py) with MultiCodeSimulator validation.
Compares four evaluation strategies for the UH FEL beamline at multiple
emittance points:

  1. FELsim:  All 137 elements in FELsim transfer matrices
  2. MC-val:  FELsim(0:87) → RF-Track(87:137), FELsim-optimised currents
  3. RFT-val: Full RF-Track(0:137), FELsim-optimised currents
  4. RFT-opt: RF-Track prefix(0:87) cached + suffix(87:137) NM-optimised

The MC-val method demonstrates production use of MultiCodeSimulator:
elements 0:87 tracked with FELsim transfer matrices (fast, exact for
linear optics), elements 87:137 tracked with RF-Track (analytical
sector-bend dipole model, 3D space charge capability).

Usage:
    python -u test/C1B_hybrid_comparison.py --smoke       # ε_n=8, 1 restart
    python -u test/C1B_hybrid_comparison.py               # ε_n={5,8,14}, 5 restarts
    python -u test/C1B_hybrid_comparison.py --emittance 8 # custom

Author: Eremey Valetov
"""

import sys
import time
import math
import argparse
import csv
from pathlib import Path
import numpy as np

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ebeam import beam
from simulatorBase import CoordinateSystem
from multiCodeSimulator import MultiCodeSimulator, SimSection

from UHM_beamline_opt_05ps_params import (
    run_optimization, compute_twiss_targets,
    ENERGY, RF_FREQ, SEGMENTS, MSE_THRESHOLDS, write_csv,
)

try:
    import RF_Track as rft
    from rftrackAdapter import RFTrackAdapter
    _RFTRACK_AVAILABLE = True
except ImportError:
    _RFTRACK_AVAILABLE = False

# Reuse UHM_rftrack_opt infrastructure
from UHM_rftrack_opt import (
    S11_QUAD_INDICES, PREFIX_END, EXCEL_PATH,
    create_beam, setup_rftrack_adapter,
    evaluate_felsim_in_rftrack, run_rftrack_stage11,
    _compute_mse,
)

JSON_PATH = Path(__file__).resolve().parent.parent.parent / 'var' / 'UH_FEL_beamline.json'
OUTDIR = Path(__file__).resolve().parent / 'results' / 'C1B'


def _print(msg):
    print(msg, flush=True)


# ── MultiCodeSimulator evaluation ───────────────────────────────────────

def evaluate_multicode(felsim_result, beam_dist, targets,
                       space_charge=False, aperture=0.5):
    """Evaluate FELsim-optimised currents via MultiCodeSimulator hybrid.

    FELsim(0:87) + RF-Track(87:137) — demonstrates production multi-code
    simulation without manual prefix caching.
    """
    ebeam_obj = beam()
    all_currents = dict(felsim_result['quad_currents'])

    rt_config = {}
    if space_charge:
        rt_config['space_charge'] = True
    if aperture != 0.5:
        rt_config['aperture'] = aperture

    # Create MultiCodeSimulator with FELsim prefix + RF-Track suffix
    mc = MultiCodeSimulator(
        sections=[
            SimSection("prefix", "felsim", (0, PREFIX_END)),
            SimSection("stage11", "rftrack", (PREFIX_END, SEGMENTS),
                       config=rt_config if rt_config else None),
        ],
        lattice_path=str(JSON_PATH),
        beam_energy=ENERGY,
    )

    # Apply quad currents to the master beamline
    for idx, current in all_currents.items():
        if idx < len(mc._master_beamline):
            mc._master_beamline[idx].current = current

    # Re-initialize simulators after current changes
    mc._init_simulators()

    result = mc.simulate(particles=beam_dist)
    if not result.success:
        return None

    ps_final = result.final_particles
    if ps_final is None or ps_final.shape[0] < 10:
        return None

    _, _, twiss_f = ebeam_obj.cal_twiss(ps_final)
    bx = twiss_f.loc['x'][r"$\beta$ (m)"]
    ax = twiss_f.loc['x'][r"$\alpha$"]
    by = twiss_f.loc['y'][r"$\beta$ (m)"]
    ay = twiss_f.loc['y'][r"$\alpha$"]

    # Dispersion: track through elements 0:93 via FELsim-only
    # (all 93 elements are in the FELsim prefix or early Stage 11)
    mc_disp = MultiCodeSimulator(
        sections=[
            SimSection("prefix", "felsim", (0, 93)),
        ],
        lattice_path=str(JSON_PATH),
        beam_energy=ENERGY,
    )
    for idx, current in all_currents.items():
        if idx < len(mc_disp._master_beamline):
            mc_disp._master_beamline[idx].current = current
    mc_disp._init_simulators()

    result_93 = mc_disp.simulate(particles=beam_dist)
    if result_93.success and result_93.final_particles is not None:
        _, _, twiss_93 = ebeam_obj.cal_twiss(result_93.final_particles)
        disp = twiss_93.loc['x'][r"$D$ (m)"]
    else:
        disp = 0.0

    mse = _compute_mse(bx, ax, by, ay, disp, targets)

    return {
        'mse': mse,
        'beta_x': bx, 'alpha_x': ax,
        'beta_y': by, 'alpha_y': ay,
        'disp': disp, 'ngood': ps_final.shape[0],
    }


# ── Main comparison loop ─────────────────────────────────────────────────

def run_comparison(emittance_points, nb_particles=500, seed=42,
                   n_restarts=5, chrom_upper_bound=15,
                   space_charge=False, aperture=0.5):
    """Run four-way comparison: FELsim / MC-val / RFT-val / RFT-opt."""
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

        # ── 1. FELsim (full 11 stages) ───────────────────────────────
        _print(f"  [1/4] FELsim optimization (11 stages)...")
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
        _print(f"    RMS = {math.sqrt(felsim['mse']):.4e} ({felsim_time:.1f} s)")

        # ── 2. MC-val (MultiCodeSimulator hybrid) ────────────────────
        _print(f"  [2/4] MultiCode validation (FELsim→RF-Track)...")
        t0 = time.perf_counter()
        mc_val = evaluate_multicode(felsim, beam_dist, targets)
        mc_time = time.perf_counter() - t0

        if mc_val:
            rows.append([
                en, 'MC-val', mc_val['mse'],
                mc_val['alpha_x'], mc_val['alpha_y'],
                mc_val['beta_x'], mc_val['beta_y'], mc_val['disp'],
                fc[87], fc[93], fc[95], fc[97],
                mc_time, 0, mc_val['ngood'],
            ])
            _print(f"    RMS = {math.sqrt(mc_val['mse']):.4e} "
                   f"(ngood={mc_val['ngood']}, {mc_time:.1f} s)")
        else:
            rows.append([en, 'MC-val'] + [float('nan')] * (len(header) - 2))
            _print(f"    FAILED")

        # ── 3. RFT-val (full RF-Track, FELsim currents) ─────────────
        _print(f"  [3/4] RF-Track validation (FELsim currents)...")
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
            _print(f"    RMS = {math.sqrt(rft_val['mse']):.4e} "
                   f"(ngood={rft_val['ngood']}, {val_time:.1f} s)")
        else:
            rows.append([en, 'RFT-val'] + [float('nan')] * (len(header) - 2))
            _print(f"    FAILED")

        # ── 4. RFT-opt (RF-Track Stage 11 optimization) ─────────────
        _print(f"  [4/4] RF-Track Stage 11 optimization "
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
            _print(f"    RMS = {math.sqrt(rft_opt['mse']):.4e} "
                   f"({rft_opt['time_s']:.1f} s, nfev={rft_opt['nfev']})")
        else:
            rows.append([en, 'RFT-opt'] + [float('nan')] * (len(header) - 2))
            _print(f"    FAILED")

    # Save results
    csv_path = OUTDIR / 'comparison.csv'
    write_csv(csv_path, header, rows)
    _print(f"\nSaved {csv_path}")

    _print_summary(rows, targets)
    _plot_comparison(rows, emittance_points)

    return rows


# ── Output ───────────────────────────────────────────────────────────────

def _print_summary(rows, targets):
    _print(f"\n{'─'*125}")
    _print(f"{'ε_n':>5} {'Method':>10} {'RMS':>12} {'β_x':>8} {'β_y':>8} "
           f"{'α_x':>8} {'α_y':>10} {'Disp':>8} "
           f"{'I_87':>6} {'I_93':>6} {'I_95':>6} {'I_97':>6} {'Time':>7}")
    _print(f"{'':>5} {'Target':>10} {'':>12} "
           f"{targets['beta_x']:>8.4f} {targets['beta_y']:>8.4f} "
           f"{targets['alpha_x']:>8.4f} {targets['alpha_y']:>10.6f}")
    _print(f"{'─'*125}")
    for r in rows:
        def _fmt(v, w, f):
            try:
                if np.isnan(v):
                    return ' ' * (w - 3) + 'N/A'
            except (TypeError, ValueError):
                pass
            return f.format(v)

        try:
            rms_val = float('nan') if np.isnan(r[2]) else math.sqrt(r[2])
        except (TypeError, ValueError):
            rms_val = math.sqrt(r[2])
        _print(f"{r[0]:>5.0f} {r[1]:>10} "
               f"{_fmt(rms_val, 12, '{:12.4e}')} "
               f"{_fmt(r[5], 8, '{:8.4f}')} {_fmt(r[6], 8, '{:8.4f}')} "
               f"{_fmt(r[3], 8, '{:8.4f}')} {_fmt(r[4], 10, '{:10.6f}')} "
               f"{_fmt(r[7], 8, '{:8.4f}')} "
               f"{_fmt(r[8], 6, '{:6.2f}')} {_fmt(r[9], 6, '{:6.2f}')} "
               f"{_fmt(r[10], 6, '{:6.2f}')} {_fmt(r[11], 6, '{:6.2f}')} "
               f"{_fmt(r[12], 7, '{:7.1f}')}")
    _print(f"{'─'*125}")


def _plot_comparison(rows, emittance_points):
    methods = ['FELsim', 'MC-val', 'RFT-val', 'RFT-opt']
    colors = {
        'FELsim': '#4477AA', 'MC-val': '#228833',
        'RFT-val': '#CCBB44', 'RFT-opt': '#EE6677',
    }
    labels = {
        'FELsim': 'FELsim (transfer matrix)',
        'MC-val': 'MultiCode (FELsim→RF-Track)',
        'RFT-val': 'RF-Track (FELsim currents)',
        'RFT-opt': 'RF-Track (optimised)',
    }

    x = np.arange(len(emittance_points))
    width = 0.2

    fig, ax = plt.subplots(figsize=(11, 5))
    for j, m in enumerate(methods):
        m_rows = [r for r in rows if r[1] == m]
        rms_vals = [math.sqrt(r[2]) for r in m_rows]
        offset = (j - 1.5) * width
        ax.bar(x + offset, rms_vals, width, label=labels[m], color=colors[m])

    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(e)) for e in emittance_points])
    ax.set_xlabel(r'Normalised emittance $\varepsilon_n$ ($\pi\cdot$mm$\cdot$mrad)')
    ax.set_ylabel('RMS Twiss Mismatch')
    ax.set_title('Stage 11 RMS Twiss Mismatch: Four-Way Comparison')
    for tl, thresh in MSE_THRESHOLDS.items():
        clr = {'Excellent': 'green', 'Acceptable': 'orange', 'Marginal': 'red'}
        rms_thresh = math.sqrt(thresh)
        ax.axhline(rms_thresh, color=clr[tl], linestyle='--', alpha=0.6,
                    label=f'{tl} ({rms_thresh:.2e})')
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3, axis='y')
    fig.tight_layout()
    fig.savefig(OUTDIR / 'mse_comparison.pdf')
    plt.close(fig)
    _print(f"  Saved {OUTDIR / 'mse_comparison.pdf'}")


# ── CLI ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='C1B: Four-way hybrid comparison (FELsim/MC/RFT-val/RFT-opt)')
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

    _print("C1B: Four-Way Hybrid Comparison")
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
