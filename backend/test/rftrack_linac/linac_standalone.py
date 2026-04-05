"""
Standalone RF-Track model of the SLAC 3-m S-band constant-gradient
travelling-wave linac.

Phase scan of a single 1 MeV electron, single Fourier coefficient
(peak gradient 13.3 MV/m). Benchmarks against the elegant TWLA
reference (29.3 MeV on-crest at 1 MeV injection).

Eremey Valetov, 2026-04-05
"""

import argparse
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import RF_Track as rft

# --- physical constants ---
MC2_MEV = 0.510998950  # electron rest mass [MeV]
CLIGHT = 299792458.0   # [m/s]

# --- SLAC 3-m S-band structure (constant-gradient TW) ---
FREQ_HZ = 2.856e9
PHASE_ADV = 2.0 * np.pi / 3.0       # 2π/3 mode
PEAK_GRADIENT = 13.3e6              # [V/m] at 12 MW (UH FEL operating point)
L_TARGET = 3.048                    # [m]
K_INJECT = 1.0                      # injection kinetic energy [MeV]

WORK_DIR = Path(__file__).parent


def build_structure(gradient_vpm=PEAK_GRADIENT, n_cells=None, phi_deg=0.0):
    """
    Construct an RF-Track TW_Structure for the SLAC linac.

    Uses a single Fourier coefficient (peak-field approximation) of the
    TM01 travelling wave. Cell length follows from synchronous phase
    velocity β_wave = 1:

        L_cell = c · Δφ / (2π · f)  = c / (3 · f)  for Δφ = 2π/3

    For f = 2856 MHz: L_cell ≈ 35.01 mm, so 87 cells → L ≈ 3.046 m,
    matching the nominal 3.048 m section length to better than 0.1%.

    Args:
        gradient_vpm: peak on-axis Ez [V/m]
        n_cells: number of cells (default: float that yields L=L_TARGET)
        phi_deg: RF phase offset [deg]; 0 = on-crest in RF-Track convention

    Returns:
        TW_Structure instance
    """
    l_cell = CLIGHT * PHASE_ADV / (2.0 * np.pi * FREQ_HZ)
    if n_cells is None:
        n_cells = L_TARGET / l_cell  # REAL-valued, fractional cells allowed

    # TW_Structure(a_0, n_first, freq, ph_advance, n_cells)
    # Scalar overload = single Fourier coefficient at harmonic n_first.
    # Positive n_cells = structure starts from middle of cell.
    tw = rft.TW_Structure(float(gradient_vpm), 0, FREQ_HZ, PHASE_ADV,
                          float(n_cells))
    tw.set_phid(phi_deg)  # phase in degrees
    return tw


def make_single_particle(K_mev=K_INJECT):
    """
    Create a single-electron Bunch6d at specified kinetic energy.

    RF-Track Bunch6d column order: [X XP Y YP T P] with units
    [mm mrad mm mrad mm/c MeV/c].
    """
    E_total = K_mev + MC2_MEV                      # total energy [MeV]
    P = np.sqrt(E_total**2 - MC2_MEV**2)            # momentum [MeV/c]
    # On-axis, on-time reference particle
    phase_space = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P]])
    # population = 1 electron per macro-particle, charge = -1 e
    bunch = rft.Bunch6d(MC2_MEV, 1.0, -1.0, phase_space)
    return bunch, P


def extract_output_energy(bunch_out):
    """Read the momentum of the (single) output particle and return K.

    Returns (NaN, NaN) if the particle was lost.
    """
    M = bunch_out.get_phase_space('%x %xp %y %yp %t %Pc')
    if M.shape[0] == 0:
        return float('nan'), float('nan')
    P_out = M[0, 5]  # MeV/c
    E_total_out = np.sqrt(P_out**2 + MC2_MEV**2)
    K_out = E_total_out - MC2_MEV
    return K_out, P_out


def track_one_phase(phi_deg, gradient_vpm=PEAK_GRADIENT, verbose=False):
    """Build a lattice with one TW_Structure at phi_deg and track a 1 MeV e-."""
    tw = build_structure(gradient_vpm=gradient_vpm, phi_deg=phi_deg)
    L = tw.get_length()

    lat = rft.Lattice()
    lat.append(tw)

    bunch_in, P_in = make_single_particle()
    bunch_out = lat.track(bunch_in)
    K_out, P_out = extract_output_energy(bunch_out)

    if verbose:
        print(f"  TW structure length: {L:.4f} m")
        print(f"  Input:  K={K_INJECT:.3f} MeV, P={P_in:.4f} MeV/c")
        print(f"  Output: K={K_out:.3f} MeV, P={P_out:.4f} MeV/c")
    return K_out, P_out, L


def phase_scan(phases_deg, gradient_vpm=PEAK_GRADIENT):
    """Run a full phase scan."""
    results = {'phase_deg': [], 'K_out': [], 'P_out': []}
    print(f"RF-Track TW phase scan, {len(phases_deg)} points, "
          f"E0={gradient_vpm/1e6:.2f} MV/m")
    for phi in phases_deg:
        K_out, P_out, L = track_one_phase(phi, gradient_vpm=gradient_vpm)
        results['phase_deg'].append(phi)
        results['K_out'].append(K_out)
        results['P_out'].append(P_out)
        status = 'LOST' if np.isnan(K_out) else f'K_out={K_out:9.4f} MeV'
        print(f"  φ={phi:7.2f}°  {status}")
    for k in results:
        results[k] = np.array(results[k])
    return results


def save_csv(results, path):
    header = 'phase_deg,K_out_MeV,P_out_MeVc'
    data = np.column_stack([
        results['phase_deg'], results['K_out'], results['P_out']
    ])
    np.savetxt(path, data, header=header, delimiter=',',
               fmt='%.8e', comments='')
    print(f"CSV saved: {path}")


def plot_scan(results, elegant_csv=None, out_path=None):
    """Plot RF-Track phase scan with elegant RFCA overlay.

    RF-Track's autophase() adjusts the structure phase so phid=0
    corresponds to on-crest (maximum energy gain). Elegant RFCA uses a
    fixed convention where PHASE=sin(phi), with peak at ~70° at 1 MeV
    injection (shifted from 90° by the low-β phase-slippage correction).
    We shift the RF-Track curve by +70° for visual alignment near the
    operating point.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    phi_rft = results['phase_deg']
    K_rft = results['K_out']
    # Map RFT phid -> elegant phase by + 70° (wrap to 0-360)
    phi_eleg_equiv = (phi_rft + 70.0) % 360.0
    ordr = np.argsort(phi_eleg_equiv)

    ax.plot(phi_eleg_equiv[ordr], K_rft[ordr], 'b.-', ms=5, lw=1.2,
            label='RF-Track TW_Structure (phid shifted by +70°)')

    if elegant_csv is not None and elegant_csv.exists():
        eleg = np.loadtxt(elegant_csv, delimiter=',', skiprows=1)
        # columns: phase_deg, K_out_MeV, dE_MeV, R11, R12, R21, R22,
        #          R56, R66, det_Rx
        ax.plot(eleg[:, 0], eleg[:, 1], 'r--', lw=1.5, alpha=0.7,
                label='elegant RFCA reference')

    ax.axhline(29.3, color='gray', lw=0.6, ls=':',
               label='elegant TWLA on-crest (29.3 MeV)')
    ax.axhline(40.0, color='green', lw=0.6, ls=':',
               label='40 MeV design target')
    ax.axvline(70, color='gray', lw=0.4, ls=':')
    ax.set_xlabel('Phase (deg), elegant convention')
    ax.set_ylabel('Output Kinetic Energy (MeV)')
    ax.set_xlim(0, 360)
    ax.set_title('SLAC 3-m TW Linac: RF-Track vs elegant Phase Scan\n'
                 f'f={FREQ_HZ/1e9:.3f} GHz, E0={PEAK_GRADIENT/1e6:.1f} MV/m, '
                 f'K_inj={K_INJECT} MeV, L={L_TARGET:.3f} m')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(alpha=0.3)
    if out_path:
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved: {out_path}")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--oncrest-only', action='store_true',
                    help='Only evaluate on-crest (single point)')
    ap.add_argument('--gradient-mv-per-m', type=float, default=13.3,
                    help='Peak gradient [MV/m] (default: 13.3 = 12 MW)')
    ap.add_argument('--nstep', type=int, default=73,
                    help='Number of phase-scan points')
    args = ap.parse_args()

    gradient_vpm = args.gradient_mv_per_m * 1e6

    if args.oncrest_only:
        print("=== On-crest sanity check ===")
        K_out, P_out, L = track_one_phase(0.0, gradient_vpm=gradient_vpm,
                                          verbose=True)
        print(f"  Expected (elegant TWLA at 1 MeV): ~29.3 MeV")
        print(f"  Delta from TWLA: {K_out - 29.3:+.3f} MeV")
        return

    phases = np.linspace(-180, 180, args.nstep)
    results = phase_scan(phases, gradient_vpm=gradient_vpm)

    # write outputs
    csv_path = WORK_DIR / 'rftrack_linac_phase_scan.csv'
    save_csv(results, csv_path)

    plot_path = WORK_DIR / 'rftrack_linac_phase_scan.pdf'
    elegant_csv = (WORK_DIR.parent / 'elegant_linac'
                   / 'phase_scan_results.csv')
    plot_scan(results, elegant_csv=elegant_csv, out_path=plot_path)

    # print peak (ignoring NaN/lost particles)
    Kvals = results['K_out']
    valid = ~np.isnan(Kvals)
    if valid.any():
        i_peak = int(np.nanargmax(Kvals))
        print(f"\nPeak: phi={results['phase_deg'][i_peak]:.2f}°  "
              f"K_out={Kvals[i_peak]:.3f} MeV")
        n_lost = int((~valid).sum())
        if n_lost:
            print(f"Lost particles: {n_lost}/{len(Kvals)} phase points")


if __name__ == '__main__':
    main()
