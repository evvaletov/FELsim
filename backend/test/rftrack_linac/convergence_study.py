"""
RF-Track convergence study for the SLAC 3-m TW S-band linac.

Three metaparameter sweeps through a single on-crest (autophased)
TW_Structure, tracking one 1-MeV electron plus two perturbations:

  1. nsteps   — number of ODE-integration steps along the structure
                (RF-Track default: ~872 for L=3.048 m, cell_len≈35 mm)
  2. epsabs   — absolute local error tolerance for the stepper
                ("to.l." in informal shorthand; RF-Track default: 1e-3)
  3. algorithm — GSL integrator: rk2 (default), rk4, rkf45, rkck,
                 rk8pd, rk2imp, rk4imp, bsimp, msadams, msbdf

Outputs:
  - convergence_results.csv  — full sweep table
  - convergence_study.pdf    — 2×2 grid (K_out, det(R_x) vs nsteps & epsabs)

Eremey Valetov, 2026-04-20
"""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import RF_Track as rft

MC2 = 0.510998950        # MeV
K_INJECT = 1.0           # MeV
P_INJECT = math.sqrt((K_INJECT + MC2)**2 - MC2**2)

LINAC_LENGTH = 3.048     # m
PEAK_GRADIENT = 13.3e6   # V/m
FREQ = 2856e6            # Hz
PHI_ADV = 2 * math.pi / 3
L_CELL_SYNC = 299792458.0 * PHI_ADV / (2 * math.pi * FREQ)
N_CELLS_AUTO = LINAC_LENGTH / L_CELL_SYNC

OUT_DIR = Path(__file__).resolve().parent / 'convergence_output'
OUT_DIR.mkdir(exist_ok=True)


@dataclass
class TrackResult:
    K_out: float
    det_Rx: float
    R11: float
    R12: float
    R21: float
    R22: float
    wall_s: float


def build_structure(nsteps=None, epsabs=None, algorithm=None):
    # n_first=0 — single Fourier coefficient of TM01 travelling wave.
    # (The standalone model and rftrackAdapter both use 0.)
    s = rft.TW_Structure(PEAK_GRADIENT, 0, FREQ, PHI_ADV, N_CELLS_AUTO)
    if algorithm is not None:
        s.set_odeint_algorithm(algorithm)
    if epsabs is not None:
        s.set_odeint_epsabs(float(epsabs))
    if nsteps is not None:
        s.set_nsteps(int(nsteps))
    s.set_phid(0.0)  # autophased — on-crest (phase-slippage-compensated)
    return s


def track_once(nsteps=None, epsabs=None, algorithm=None, dx=0.01) -> TrackResult:
    s = build_structure(nsteps=nsteps, epsabs=epsabs, algorithm=algorithm)
    lat = rft.Lattice()
    lat.append(s)

    ps = np.array([
        [0.0, 0.0, 0.0, 0.0, 0.0, P_INJECT],
        [dx,  0.0, 0.0, 0.0, 0.0, P_INJECT],
        [0.0, dx,  0.0, 0.0, 0.0, P_INJECT],
    ])
    bunch = rft.Bunch6d(MC2, 1.0, -1.0, ps)
    t0 = time.perf_counter()
    bout = lat.track(bunch)
    wall = time.perf_counter() - t0
    M = bout.get_phase_space('%x %xp %y %yp %t %Pc')
    if M.shape[0] < 3:
        return TrackResult(*(np.nan,)*6, wall_s=wall)
    K_out = math.sqrt(M[0, 5]**2 + MC2**2) - MC2
    R11 = (M[1, 0] - M[0, 0]) / dx
    R21 = (M[1, 1] - M[0, 1]) / dx
    R12 = (M[2, 0] - M[0, 0]) / dx
    R22 = (M[2, 1] - M[0, 1]) / dx
    det = R11 * R22 - R12 * R21
    return TrackResult(K_out, det, R11, R12, R21, R22, wall)


def sweep_nsteps(nsteps_values, algorithm='rk2', epsabs=1e-3):
    rows = []
    for n in nsteps_values:
        r = track_once(nsteps=n, epsabs=epsabs, algorithm=algorithm)
        rows.append({
            'param': 'nsteps', 'value': n,
            'K_out_MeV': r.K_out, 'det_Rx': r.det_Rx,
            'wall_s': r.wall_s,
        })
        print(f"  nsteps={n:6d}  K_out={r.K_out:.6f} MeV  "
              f"det(R_x)={r.det_Rx:.6f}  t={r.wall_s:.3f} s")
    return rows


def sweep_epsabs(epsabs_values, algorithm='rk2', nsteps=None):
    rows = []
    for eps in epsabs_values:
        r = track_once(nsteps=nsteps, epsabs=eps, algorithm=algorithm)
        rows.append({
            'param': 'epsabs', 'value': eps,
            'K_out_MeV': r.K_out, 'det_Rx': r.det_Rx,
            'wall_s': r.wall_s,
        })
        print(f"  epsabs={eps:.1e}  K_out={r.K_out:.6f} MeV  "
              f"det(R_x)={r.det_Rx:.6f}  t={r.wall_s:.3f} s")
    return rows


def sweep_algorithm(algorithms, epsabs=1e-3, nsteps=None):
    rows = []
    for algo in algorithms:
        r = track_once(nsteps=nsteps, epsabs=epsabs, algorithm=algo)
        rows.append({
            'param': 'algorithm', 'value': algo,
            'K_out_MeV': r.K_out, 'det_Rx': r.det_Rx,
            'wall_s': r.wall_s,
        })
        print(f"  {algo:10s}  K_out={r.K_out:.6f} MeV  "
              f"det(R_x)={r.det_Rx:.6f}  t={r.wall_s:.3f} s")
    return rows


def plot_results(n_rows, e_rows, a_rows, out_pdf):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    (ax_kn, ax_de), (ax_ke, ax_t) = axes

    # Reference: finest converged run (largest nsteps, smallest epsabs)
    n_arr = np.array([r['value'] for r in n_rows])
    n_K = np.array([r['K_out_MeV'] for r in n_rows])
    n_det = np.array([r['det_Rx'] for r in n_rows])
    n_t = np.array([r['wall_s'] for r in n_rows])
    K_ref = n_K[-1]
    det_ref = n_det[-1]

    e_arr = np.array([r['value'] for r in e_rows])
    e_K = np.array([r['K_out_MeV'] for r in e_rows])
    e_det = np.array([r['det_Rx'] for r in e_rows])
    e_t = np.array([r['wall_s'] for r in e_rows])

    # K_out vs nsteps (and convergence bound)
    ax_kn.semilogx(n_arr, np.abs(n_K - K_ref), 'b.-', lw=1.2)
    ax_kn.set_xlabel('nsteps (integration steps)')
    ax_kn.set_ylabel('|K_out − K_ref|  (MeV)')
    ax_kn.set_yscale('log')
    ax_kn.set_title(f'K_out convergence vs nsteps  (K_ref={K_ref:.6f} MeV)')
    ax_kn.grid(alpha=0.3)
    ax_kn.axhline(1e-4, color='gray', lw=0.5, ls=':')

    # det(R_x) vs nsteps
    ax_de.semilogx(n_arr, np.abs(n_det - det_ref), 'r.-', lw=1.2)
    ax_de.set_xlabel('nsteps')
    ax_de.set_ylabel('|det(R_x) − det_ref|')
    ax_de.set_yscale('log')
    ax_de.set_title(f'det(R_x) convergence vs nsteps  (det_ref={det_ref:.6f})')
    ax_de.grid(alpha=0.3)

    # K_out vs epsabs (note: epsabs only active under adaptive stepping)
    ax_ke.loglog(e_arr, np.abs(e_K - K_ref), 'g.-', lw=1.2,
                 label='|K_out − K_ref|')
    ax_ke.loglog(e_arr, np.abs(e_det - det_ref), 'm.-', lw=1.2,
                 label='|det(R_x) − det_ref|')
    ax_ke.set_xlabel('epsabs  (ODE local error tolerance)')
    ax_ke.set_ylabel('error magnitude')
    ax_ke.set_title('Convergence vs epsabs  (nsteps = default)')
    ax_ke.legend(fontsize=9)
    ax_ke.grid(alpha=0.3)
    ax_ke.invert_xaxis()

    # Wall time
    ax_t.semilogx(n_arr, n_t, 'b.-', lw=1.2, label='sweep nsteps')
    if len(e_t) > 0:
        ax_t.semilogx(e_arr, e_t, 'g.--', lw=1.0, label='sweep epsabs')
    ax_t.set_xlabel('parameter value (log scale)')
    ax_t.set_ylabel('wall time (s)')
    ax_t.set_title('Cost of the sweeps')
    ax_t.legend(fontsize=9)
    ax_t.grid(alpha=0.3)

    # Annotation block listing algorithm comparison
    if a_rows:
        algo_text = '  '.join(f"{r['value']}:{r['K_out_MeV']:.5f}"
                              for r in a_rows)
        fig.text(0.5, 0.01, f"Algorithm sweep (epsabs=1e-3, default nsteps):  {algo_text}",
                 ha='center', fontsize=7, family='monospace')

    fig.suptitle('RF-Track TW_Structure convergence — SLAC 3 m, 1 MeV → on-crest')
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(out_pdf, dpi=150)
    plt.close()
    print(f"  Saved {out_pdf}")


def write_csv(rows, out_path):
    with open(out_path, 'w') as f:
        f.write('param,value,K_out_MeV,det_Rx,wall_s\n')
        for r in rows:
            f.write(f"{r['param']},{r['value']},{r['K_out_MeV']:.10e},"
                    f"{r['det_Rx']:.10e},{r['wall_s']:.6e}\n")
    print(f"  Saved {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--skip-algo', action='store_true',
                        help='skip algorithm comparison sweep')
    parser.add_argument('--quick', action='store_true',
                        help='smaller sweeps for a quick check')
    args = parser.parse_args()

    if args.quick:
        nsteps_values = [50, 200, 872, 3000]
        epsabs_values = [1e-1, 1e-3, 1e-6, 1e-9]
    else:
        nsteps_values = [20, 50, 100, 200, 500, 872, 2000, 5000, 10000]
        epsabs_values = [1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10]

    algorithms = ['rk2', 'rk4', 'rkf45', 'rkck', 'rk8pd']

    # First, show defaults for reference
    s_ref = build_structure()
    print("RF-Track TW_Structure defaults:")
    print(f"  nsteps     = {s_ref.get_nsteps()}")
    print(f"  algorithm  = {s_ref.get_odeint_algorithm()}")
    print(f"  epsabs     = {s_ref.get_odeint_epsabs()}")
    print(f"  epsrel     = {s_ref.get_odeint_epsrel()}")
    print(f"  length     = {s_ref.get_length()} m, cell_len = {s_ref.get_cell_length():.6f} m")
    print()

    print("Sweep 1: nsteps (epsabs=1e-3, algo=rk2)")
    n_rows = sweep_nsteps(nsteps_values)

    print("\nSweep 2: epsabs (nsteps=default, algo=rk2)")
    e_rows = sweep_epsabs(epsabs_values)

    a_rows = []
    if not args.skip_algo:
        print("\nSweep 3: algorithm (epsabs=1e-3, nsteps=default)")
        a_rows = sweep_algorithm(algorithms)

    all_rows = n_rows + e_rows + a_rows
    write_csv(all_rows, OUT_DIR / 'convergence_results.csv')
    plot_results(n_rows, e_rows, a_rows, OUT_DIR / 'convergence_study.pdf')

    print("\nFinest converged reference:")
    print(f"  K_out  = {n_rows[-1]['K_out_MeV']:.6f} MeV  "
          f"(nsteps={n_rows[-1]['value']})")
    print(f"  det(R) = {n_rows[-1]['det_Rx']:.6f}")


if __name__ == '__main__':
    main()
