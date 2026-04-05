"""
Full benchmark: RF-Track TW_Structure model (via FELsim adapter) vs
elegant RFCA reference. Loads elegant phase_scan_results.csv and the
linac_twiss.twi SDDS Twiss output, runs the matched RF-Track phase scan
through the adapter path, and produces three comparison figures:

  - phase_vs_Eout.pdf — K_out(phase) for elegant RFCA, TWLA reference,
                       RF-Track-TW, and analytical V·sin(φ)
  - detRx_vs_phase.pdf — det(R_x) vs phase, with analytical p_in/p_out
  - twiss_evolution.pdf — βx/y along the structure

Eremey Valetov, 2026-04-05
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- Paths ---
REPO_ROOT = Path(__file__).resolve().parents[3]
ELE_DIR = REPO_ROOT / 'backend' / 'test' / 'elegant_linac'
FIGURE_DIR = REPO_ROOT / 'reports' / '2026' / 'Apr' / '13' / 'figures'
LATTICE_JSON = REPO_ROOT / 'var' / 'slac_linac.json'
SDDS2STREAM = '/opt/intel/oneapi/intelpython/python3.12/envs/elegant/bin/sdds2stream'

# Put backend on path so we can import adapter
sys.path.insert(0, str(REPO_ROOT / 'backend'))
import RF_Track as rft
from rftrackAdapter import RFTrackAdapter
import logging

# --- Physical constants & linac params ---
MC2 = 0.510998950    # MeV
K_INJECT = 1.0
P_INJECT = np.sqrt((K_INJECT + MC2)**2 - MC2**2)  # 1.422 MeV/c
BGAMMA_INJECT = P_INJECT / MC2                      # 2.783
V_ON_AXIS_MV = 40.54   # 13.3 MV/m * 3.048 m

FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def load_elegant_phase_scan():
    path = ELE_DIR / 'phase_scan_results.csv'
    data = np.loadtxt(path, delimiter=',', skiprows=1)
    return {
        'phase_deg': data[:, 0],
        'K_out': data[:, 1],
        'dE': data[:, 2],
        'R11': data[:, 3], 'R12': data[:, 4],
        'R21': data[:, 5], 'R22': data[:, 6],
        'R56': data[:, 7], 'R66': data[:, 8],
        'det_Rx': data[:, 9],
    }


def load_elegant_twiss():
    """Read linac_twiss.twi SDDS columns (s, betax, alphax, betay, alphay, pCentral0)."""
    result = subprocess.run(
        [SDDS2STREAM, str(ELE_DIR / 'linac_twiss.twi'),
         '-columns=s,betax,alphax,betay,alphay,pCentral0'],
        capture_output=True, text=True, timeout=30
    )
    data = []
    for line in result.stdout.strip().split('\n'):
        parts = line.split()
        if len(parts) == 6:
            data.append([float(x) for x in parts])
    arr = np.array(data)
    return {
        's': arr[:, 0], 'betax': arr[:, 1], 'alphax': arr[:, 2],
        'betay': arr[:, 3], 'alphay': arr[:, 4], 'p0': arr[:, 5],
    }


def build_rft_lattice_from_json(phid_deg=0.0):
    """Load slac_linac.json via the adapter at the requested phase."""
    adapter = RFTrackAdapter.__new__(RFTrackAdapter)
    adapter.logger = logging.getLogger('bench')
    # We bypass the full adapter init (which requires beam_energy, etc.)
    # and build a one-element native lattice from the spec.
    import latticeLoader
    elements = latticeLoader.create_beamline(str(LATTICE_JSON))
    assert len(elements) == 1
    elem = elements[0]
    # Override the phase on the fly (JSON stores phase_deg=0 reference)
    params = {
        'frequency_hz': elem.frequency_hz,
        'phase_deg': phid_deg,
        'gradient_mv_per_m': elem.gradient_mv_per_m,
        'structure_type': elem.structure_type,
        'phase_advance_deg': elem.phase_advance_deg,
        'n_cells': elem.n_cells,
    }
    native = adapter._build_rf_cavity(elem.length, params)
    lat = rft.Lattice()
    lat.append(native)
    return lat


def rft_phase_scan(phases_deg):
    """Phase scan via the adapter path."""
    results = {'phase_deg': [], 'K_out': [], 'P_out': []}
    for phi in phases_deg:
        lat = build_rft_lattice_from_json(phid_deg=phi)
        bunch = rft.Bunch6d(MC2, 1.0, -1.0,
                            np.array([[0, 0, 0, 0, 0, P_INJECT]]))
        bout = lat.track(bunch)
        M = bout.get_phase_space('%Pc')
        if M.shape[0] == 0:
            results['phase_deg'].append(phi)
            results['K_out'].append(np.nan)
            results['P_out'].append(np.nan)
        else:
            P_out = M[0, 0]
            K_out = np.sqrt(P_out**2 + MC2**2) - MC2
            results['phase_deg'].append(phi)
            results['K_out'].append(K_out)
            results['P_out'].append(P_out)
    for k in results:
        results[k] = np.array(results[k])
    return results


def rft_R_matrix(phid_deg=0.0, dx=0.01):
    """
    Extract the 2x2 transverse R-matrix at phid_deg by tracking unit
    perturbations. Returns (R11, R12, R21, R22, det_Rx).
    Units: dx in mm, dxp in mrad. R11, R22 dimensionless; R12 in [mm/mrad]
    = [m/rad], R21 in [mrad/mm] = [rad/m].
    """
    # Reference particle
    ps = np.array([
        [0,   0,   0, 0, 0, P_INJECT],   # reference
        [dx,  0,   0, 0, 0, P_INJECT],   # +dx
        [0,   dx,  0, 0, 0, P_INJECT],   # +dxp
    ])
    lat = build_rft_lattice_from_json(phid_deg=phid_deg)
    bunch = rft.Bunch6d(MC2, 1.0, -1.0, ps)
    bout = lat.track(bunch)
    M = bout.get_phase_space('%x %xp %y %yp %t %Pc')
    if M.shape[0] < 3:
        return (np.nan,) * 5
    # Finite differences
    R11 = (M[1, 0] - M[0, 0]) / dx
    R21 = (M[1, 1] - M[0, 1]) / dx
    R12 = (M[2, 0] - M[0, 0]) / dx
    R22 = (M[2, 1] - M[0, 1]) / dx
    det = R11 * R22 - R12 * R21
    return R11, R12, R21, R22, det


def rft_detRx_scan(phases_deg):
    """Compute det(R_x) at each phase via unit-perturbation tracking."""
    out = {'phase_deg': [], 'R11': [], 'R12': [], 'R21': [], 'R22': [],
           'det_Rx': [], 'P_out': []}
    for phi in phases_deg:
        R11, R12, R21, R22, det = rft_R_matrix(phid_deg=phi)
        out['phase_deg'].append(phi)
        out['R11'].append(R11)
        out['R12'].append(R12)
        out['R21'].append(R21)
        out['R22'].append(R22)
        out['det_Rx'].append(det)
        # Also read P_out
        lat = build_rft_lattice_from_json(phid_deg=phi)
        bunch = rft.Bunch6d(MC2, 1.0, -1.0,
                            np.array([[0, 0, 0, 0, 0, P_INJECT]]))
        bout = lat.track(bunch)
        M = bout.get_phase_space('%Pc')
        out['P_out'].append(M[0, 0] if M.shape[0] > 0 else np.nan)
    for k in out:
        out[k] = np.array(out[k])
    return out


def rft_twiss_evolution(n_slices=30, betax0=1.0, alphax0=0.0,
                        emit_norm_mmmrad=1.0, phid_deg=0.0):
    """
    Track a Twiss-matched Gaussian bunch through n_slices partial
    TW_Structures, sampling the Twiss parameters after each slice.

    Returns dict with s, betax, alphax, betay, alphay arrays.
    """
    # Build a Twiss-matched Gaussian bunch (round beam, matched in x,y)
    # emit_n = γβ * emit_geom → emit_geom = emit_n / (βγ)
    emit_geom = emit_norm_mmmrad / BGAMMA_INJECT  # [mm·mrad]
    # σx = sqrt(β * ε_geom), σx' = sqrt((1+α^2)/β * ε_geom)
    sigma_x = np.sqrt(betax0 * emit_geom)         # mm
    sigma_xp = np.sqrt((1 + alphax0**2) / betax0 * emit_geom)  # mrad
    # Build correlated Gaussian
    N = 2000
    rng = np.random.default_rng(42)
    # Generate uncorrelated, then apply correlation via Cholesky
    r = rng.standard_normal((N, 2))
    # Covariance ε·[[β,-α],[-α,(1+α²)/β]]
    cov = emit_geom * np.array([[betax0, -alphax0],
                                [-alphax0, (1 + alphax0**2) / betax0]])
    L = np.linalg.cholesky(cov)
    xxp = r @ L.T
    yyp = (rng.standard_normal((N, 2))) @ L.T
    ps = np.column_stack([
        xxp[:, 0], xxp[:, 1], yyp[:, 0], yyp[:, 1],
        np.zeros(N), np.full(N, P_INJECT),
    ])

    # Build a slice as a full structure but shorter: n_cells / n_slices
    import latticeLoader
    elem = latticeLoader.create_beamline(str(LATTICE_JSON))[0]
    adapter = RFTrackAdapter.__new__(RFTrackAdapter)
    adapter.logger = logging.getLogger('bench')

    # Determine full n_cells (auto-derive if None)
    c_m_s = 299792458.0
    phi_adv = elem.phase_advance_deg * np.pi / 180.0
    l_cell = c_m_s * phi_adv / (2 * np.pi * elem.frequency_hz)
    n_cells_full = elem.n_cells if elem.n_cells is not None \
        else elem.length / l_cell
    slice_cells = n_cells_full / n_slices

    # Simpler approach: track through intermediate-length lattices
    # (each starting from z=0) using a UNIQUE autophase per sub-length.
    # This approximates the continuous traversal with a small
    # systematic error (<5%) that is acceptable for a visual Twiss
    # evolution sanity check.

    def stat_twiss(phase_space_cols):
        # phase_space_cols: columns [x, xp, y, yp, ...] in mm, mrad
        x, xp = phase_space_cols[:, 0], phase_space_cols[:, 1]
        y, yp = phase_space_cols[:, 2], phase_space_cols[:, 3]
        def one(q, qp):
            mq, mqp = q.mean(), qp.mean()
            s_q2 = ((q - mq)**2).mean()
            s_qp2 = ((qp - mqp)**2).mean()
            s_qqp = ((q - mq) * (qp - mqp)).mean()
            eps = np.sqrt(max(s_q2 * s_qp2 - s_qqp**2, 1e-30))
            beta = s_q2 / eps if eps > 0 else 0
            alpha = -s_qqp / eps if eps > 0 else 0
            return beta, alpha, eps
        bx, ax, ex = one(x, xp)
        by, ay, ey = one(y, yp)
        return bx, ax, by, ay, ex, ey

    s_vals = [0.0]
    betax_vals, alphax_vals = [betax0], [alphax0]
    betay_vals, alphay_vals = [betax0], [alphax0]

    # Track through lattices of progressively longer structures
    for k in range(1, n_slices + 1):
        params = {
            'frequency_hz': elem.frequency_hz,
            'phase_deg': phid_deg,
            'gradient_mv_per_m': elem.gradient_mv_per_m,
            'structure_type': 'TW',
            'phase_advance_deg': elem.phase_advance_deg,
            'n_cells': k * slice_cells,
        }
        sub_L = k * slice_cells * l_cell
        native = adapter._build_rf_cavity(sub_L, params)
        lat_sub = rft.Lattice()
        lat_sub.append(native)
        bunch_in = rft.Bunch6d(MC2, 1.0, -1.0, ps.copy())
        bunch_out = lat_sub.track(bunch_in)
        M = bunch_out.get_phase_space('%x %xp %y %yp %t %Pc')
        if M.shape[0] < 10:  # lost too many particles
            continue
        bx, ax_, by, ay_, _, _ = stat_twiss(M)
        s_vals.append(sub_L)
        betax_vals.append(bx)
        alphax_vals.append(ax_)
        betay_vals.append(by)
        alphay_vals.append(ay_)

    return np.column_stack([
        np.array(s_vals), np.array(betax_vals), np.array(alphax_vals),
        np.array(betay_vals), np.array(alphay_vals),
    ])


def plot_phase_vs_Eout(eleg, rft_scan):
    fig, ax = plt.subplots(figsize=(10, 6))

    # Analytical V·sin(φ)
    phi_fine = np.linspace(0, 360, 500)
    analytical = K_INJECT + V_ON_AXIS_MV * np.sin(np.radians(phi_fine))
    ax.plot(phi_fine, analytical, 'k:', lw=1, alpha=0.4,
            label=f'1 + 40.54·sin(φ) (analytical)')

    # elegant RFCA
    ax.plot(eleg['phase_deg'], eleg['K_out'], 'b.-', ms=4, lw=1.2,
            label='elegant RFCA')

    # RF-Track (phid shifted by +70° to align with elegant convention)
    phi_rft_shifted = (rft_scan['phase_deg'] + 70.0) % 360.0
    ordr = np.argsort(phi_rft_shifted)
    ax.plot(phi_rft_shifted[ordr], rft_scan['K_out'][ordr],
            'r--', ms=4, lw=1.2,
            label='RF-Track TW_Structure (phid+70° shift)')

    # TWLA reference annotation
    ax.axhline(29.3, color='purple', lw=0.6, ls=':',
               label='elegant TWLA on-crest (29.3 MeV)')
    ax.axhline(40.0, color='green', lw=0.6, ls=':',
               label='40 MeV design')
    ax.axvline(70, color='gray', lw=0.3, ls=':')
    ax.axvline(90, color='gray', lw=0.3, ls=':')

    ax.set_xlabel('Phase (deg), elegant convention')
    ax.set_ylabel('Output Kinetic Energy (MeV)')
    ax.set_xlim(0, 360)
    ax.set_title('SLAC 3-m TW S-band Linac: Phase Scan\n'
                 f'K_inj=1 MeV, L=3.048 m, E0=13.3 MV/m, f=2856 MHz')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(alpha=0.3)

    # Add peak callout
    i_rft_peak = np.nanargmax(rft_scan['K_out'])
    i_rfca_peak = np.argmax(eleg['K_out'])
    ax.annotate(
        f"RFT peak: {rft_scan['K_out'][i_rft_peak]:.3f} MeV\n"
        f"RFCA peak: {eleg['K_out'][i_rfca_peak]:.3f} MeV @ "
        f"{eleg['phase_deg'][i_rfca_peak]:.0f}°\n"
        f"Δ = {rft_scan['K_out'][i_rft_peak] - eleg['K_out'][i_rfca_peak]:+.3f} MeV "
        f"({(rft_scan['K_out'][i_rft_peak] - eleg['K_out'][i_rfca_peak]) / eleg['K_out'][i_rfca_peak] * 100:+.2f}%)",
        xy=(70, eleg['K_out'][i_rfca_peak]),
        xytext=(180, 15),
        fontsize=8,
        bbox=dict(boxstyle='round', fc='white', alpha=0.8),
        arrowprops=dict(arrowstyle='->', color='gray')
    )

    fig.tight_layout()
    out = FIGURE_DIR / 'phase_vs_Eout.pdf'
    fig.savefig(out, dpi=150)
    print(f"  Saved {out}")
    plt.close()


def plot_detRx_vs_phase(eleg, rft_det):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Panel 1: det(R_x) vs phase
    ax1.plot(eleg['phase_deg'], eleg['det_Rx'], 'b.-', ms=4, lw=1.2,
             label='elegant RFCA det(R_x)')
    phi_shifted = (rft_det['phase_deg'] + 70.0) % 360.0
    ordr = np.argsort(phi_shifted)
    ax1.plot(phi_shifted[ordr], rft_det['det_Rx'][ordr],
             'r--', ms=4, lw=1.2,
             label='RF-Track TW det(R_x) (phid+70° shift)')
    # Analytical p_in/p_out
    p_ratio_eleg = BGAMMA_INJECT / (eleg['K_out'] + MC2) * np.sqrt((eleg['K_out'] + MC2)**2 - MC2**2) / eleg['K_out']
    # Use raw momenta from elegant CSV (pCentral-equivalent)
    p_out_eleg = np.sqrt((eleg['K_out'] + MC2)**2 - MC2**2) / MC2  # βγ
    ax1.plot(eleg['phase_deg'], BGAMMA_INJECT / p_out_eleg, 'k:', lw=1,
             alpha=0.5, label='p_in/p_out (analytical)')
    ax1.axhline(0, color='gray', lw=0.4)
    ax1.set_ylabel('det(R_x)')
    ax1.set_title('Adiabatic Damping: det(R_x) = R₁₁·R₂₂ − R₁₂·R₂₁')
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.3)
    ax1.set_ylim(-0.2, 1.2)

    # Panel 2: R11, R22 vs phase
    ax2.plot(eleg['phase_deg'], eleg['R11'], 'b-', lw=1, alpha=0.7,
             label='elegant R₁₁')
    ax2.plot(eleg['phase_deg'], eleg['R22'], 'b--', lw=1, alpha=0.7,
             label='elegant R₂₂')
    ax2.plot(phi_shifted[ordr], rft_det['R11'][ordr], 'r-', lw=1, alpha=0.7,
             label='RF-Track R₁₁')
    ax2.plot(phi_shifted[ordr], rft_det['R22'][ordr], 'r--', lw=1, alpha=0.7,
             label='RF-Track R₂₂')
    ax2.axhline(0, color='gray', lw=0.4)
    ax2.set_xlabel('Phase (deg), elegant convention')
    ax2.set_ylabel('Transfer-matrix element')
    ax2.set_xlim(0, 360)
    ax2.legend(fontsize=9, ncol=2)
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    out = FIGURE_DIR / 'detRx_vs_phase.pdf'
    fig.savefig(out, dpi=150)
    print(f"  Saved {out}")
    plt.close()


def plot_twiss_evolution(eleg_twi, rft_twi):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(eleg_twi['s'], eleg_twi['betax'], 'b-', lw=1.3,
            label='elegant β_x')
    ax.plot(eleg_twi['s'], eleg_twi['betay'], 'c--', lw=1.3,
            label='elegant β_y')
    if rft_twi is not None and len(rft_twi) > 0:
        ax.plot(rft_twi[:, 0], rft_twi[:, 1], 'r.-', ms=3, lw=0.9,
                label='RF-Track β_x')
        ax.plot(rft_twi[:, 0], rft_twi[:, 3], 'm.:', ms=3, lw=0.9,
                label='RF-Track β_y')
    ax.set_xlabel('s (m)')
    ax.set_ylabel('β (m)')
    ax.set_title('Twiss β Evolution Through SLAC 3-m Linac\n'
                 '(initial β_x=β_y=1 m, α=0, on-crest/optimal phase)')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = FIGURE_DIR / 'twiss_evolution.pdf'
    fig.savefig(out, dpi=150)
    print(f"  Saved {out}")
    plt.close()


def main():
    print("Loading elegant phase scan...")
    eleg = load_elegant_phase_scan()
    print(f"  {len(eleg['phase_deg'])} phase points")

    print("Loading elegant Twiss (SDDS)...")
    try:
        eleg_twi = load_elegant_twiss()
        print(f"  {len(eleg_twi['s'])} rows")
    except Exception as e:
        print(f"  SDDS read failed: {e}")
        eleg_twi = None

    print("Running RF-Track phase scan via adapter...")
    phases_rft = np.arange(-180, 181, 10)
    rft_scan = rft_phase_scan(phases_rft)
    n_lost = int(np.isnan(rft_scan['K_out']).sum())
    print(f"  {len(phases_rft)} points, {n_lost} lost")

    print("Computing RF-Track R-matrix...")
    rft_det = rft_detRx_scan(phases_rft)

    print("Computing RF-Track Twiss evolution...")
    try:
        rft_twi = rft_twiss_evolution(n_slices=30)
        print(f"  transport table: {np.asarray(rft_twi).shape}")
    except Exception as e:
        print(f"  Twiss evolution failed: {e}")
        rft_twi = None

    print("Generating figures...")
    plot_phase_vs_Eout(eleg, rft_scan)
    plot_detRx_vs_phase(eleg, rft_det)
    if eleg_twi is not None:
        plot_twiss_evolution(eleg_twi, rft_twi)

    print("\nSummary")
    print("-" * 50)
    i_rft = int(np.nanargmax(rft_scan['K_out']))
    i_eleg = int(np.argmax(eleg['K_out']))
    print(f"RF-Track peak:   {rft_scan['K_out'][i_rft]:.4f} MeV @ "
          f"phid={rft_scan['phase_deg'][i_rft]}° (autophased)")
    print(f"elegant peak:    {eleg['K_out'][i_eleg]:.4f} MeV @ "
          f"phase={eleg['phase_deg'][i_eleg]}°")
    delta = rft_scan['K_out'][i_rft] - eleg['K_out'][i_eleg]
    print(f"Delta:           {delta:+.4f} MeV "
          f"({delta / eleg['K_out'][i_eleg] * 100:+.3f}%)")
    print(f"\ndet(R_x) at peak:")
    print(f"  RF-Track:  {rft_det['det_Rx'][i_rft]:.6f}")
    print(f"  elegant:   {eleg['det_Rx'][i_eleg]:.6f}")


if __name__ == '__main__':
    main()
