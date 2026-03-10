#!/usr/bin/env python3
"""Phase 1 diagnostics: RF-Track vs FELsim/COSY model alignment.

D5: Asymmetric edge kick (Tx vs Ty) for all DPW elements
D2: Twiss dispersion subtraction comparison
D7: Total lattice path length comparison

Author: Eremey Valetov
"""

import sys
import json
from pathlib import Path
import numpy as np

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

from ebeam import beam
from beamline import lattice
from excelElements import ExcelElements
from simulatorBase import CoordinateSystem

EXCEL_PATH = backend_dir.parent / 'beam_excel' / 'Beamline_elements.xlsx'
ENERGY = 40
EPSILON_N = 8
X_STD = 0.8
Y_STD = 0.8
FREQ = 2856e6
BUNCH_SPREAD = 2
ENERGY_STD_PCT = 0.5
H_CHIRP = 5e9
NB_PARTICLES = 2000
SEGMENTS = 118

RESULTS_DIR = Path(__file__).resolve().parent / 'results' / 'D1'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def build_felsim_line():
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    excel = ExcelElements(str(EXCEL_PATH))
    beamlineUH = excel.create_beamline()
    line = relat.changeBeamType("electron", ENERGY, beamlineUH)[:SEGMENTS]
    return line, relat


def generate_beam(relat):
    np.random.seed(42)
    norm = relat.gamma * relat.beta
    epsilon = EPSILON_N / norm
    beam_dist = beam().gen_6d_gaussian(
        0, [X_STD, epsilon / X_STD, Y_STD, epsilon / Y_STD,
            BUNCH_SPREAD * 1e-9 * FREQ, ENERGY_STD_PCT * 10],
        NB_PARTICLES)
    beam_dist[:, 5] += H_CHIRP * beam_dist[:, 4] / FREQ
    return beam_dist


# ═══════════════════════════════════════════════════════════════════════════════
# D5: Asymmetric edge kick analysis — Tx vs Ty for each DPW
# ═══════════════════════════════════════════════════════════════════════════════

def diagnostic_D5(line):
    print("=" * 72)
    print("  D5: DPW asymmetric edge kick (Tx vs Ty)")
    print("=" * 72)
    print(f"  {'idx':>4} {'name':<10} {'η (deg)':>8} {'R (m)':>8} {'g (mm)':>8} "
          f"{'le (mm)':>8} {'K':>8} {'φ (mrad)':>8} "
          f"{'M21':>10} {'M43_fel':>10} {'M43_rft':>10} {'Δ(%)':>8}")
    print(f"  {'-' * 70}")

    dpw_data = []
    for idx, elem in enumerate(line):
        cls_name = type(elem).__name__
        if cls_name != 'dipole_wedge':
            continue

        eta_deg = elem.angle
        eta = np.radians(eta_deg)
        R = elem.dipole_length / (abs(elem.dipole_angle) * np.pi / 180)
        g = elem.pole_gap  # mm
        le = elem.length    # mm (wedge length)

        # Triangle model fringe correction
        K = le / (6.0 * g)
        h = 1.0 / R
        phi = K * g * h * (1 + np.sin(eta)**2) / np.cos(eta)

        Tx = np.tan(eta)
        Ty = np.tan(eta - phi)

        # FELsim matrix elements
        M21_felsim = Tx / R       # horizontal kick
        M43_felsim = -Ty / R      # vertical kick (with fringe correction)

        # RF-Track: thin-lens quad with K1L = -K0 * tan(eta)
        K0 = abs(np.radians(elem.dipole_angle) / elem.dipole_length)
        K1L = -K0 * np.tan(eta)
        # Thin-lens quad: Δx' = -K1L * x → effective M21 = -K1L = K0*tan(η)
        # Thin-lens quad: Δy' = +K1L * y → effective M43 = +K1L = -K0*tan(η)
        M43_rftrack = K1L  # = -K0 * tan(eta) = -Tx/R

        delta_pct = (M43_rftrack - M43_felsim) / abs(M43_felsim) * 100 if M43_felsim != 0 else 0

        name = getattr(elem, 'name', '') or ''
        print(f"  {idx:4d} {name:<10} {eta_deg:8.3f} {R:8.4f} {g:8.3f} "
              f"{le:8.3f} {K:8.4f} {phi*1e3:8.3f} "
              f"{M21_felsim:10.6f} {M43_felsim:10.6f} {M43_rftrack:10.6f} {delta_pct:8.2f}")

        dpw_data.append({
            'index': idx, 'name': name, 'eta_deg': eta_deg,
            'R_m': R, 'pole_gap_mm': g, 'wedge_length_mm': le,
            'K_triangle': K, 'phi_rad': phi,
            'Tx': Tx, 'Ty': Ty,
            'M21_felsim': M21_felsim, 'M43_felsim': M43_felsim,
            'M43_rftrack': M43_rftrack,
            'delta_pct': delta_pct,
        })

    # Cumulative effect: multiply all DPW matrices through
    cumul_felsim = np.eye(6)
    cumul_rftrack = np.eye(6)
    for d in dpw_data:
        M_f = np.eye(6)
        M_f[1, 0] = d['M21_felsim']
        M_f[3, 2] = d['M43_felsim']
        cumul_felsim = M_f @ cumul_felsim

        M_r = np.eye(6)
        M_r[1, 0] = -d['M43_rftrack']  # RF-Track M21 = -K1L = K0*tan(η) = Tx/R
        M_r[3, 2] = d['M43_rftrack']
        cumul_rftrack = M_r @ cumul_rftrack

    print(f"\n  Cumulative edge-kick matrices (DPW only):")
    print(f"  {'':>14} {'M[1,0]':>10} {'M[3,2]':>10}")
    print(f"  {'FELsim':<14} {cumul_felsim[1,0]:10.6f} {cumul_felsim[3,2]:10.6f}")
    print(f"  {'RF-Track':<14} {cumul_rftrack[1,0]:10.6f} {cumul_rftrack[3,2]:10.6f}")
    diff_10 = (cumul_rftrack[1,0] - cumul_felsim[1,0]) / abs(cumul_felsim[1,0]) * 100
    diff_32 = (cumul_rftrack[3,2] - cumul_felsim[3,2]) / abs(cumul_felsim[3,2]) * 100
    print(f"  {'Δ (%)':<14} {diff_10:10.2f} {diff_32:10.2f}")

    return dpw_data


# ═══════════════════════════════════════════════════════════════════════════════
# D2: Twiss dispersion subtraction comparison
# ═══════════════════════════════════════════════════════════════════════════════

def diagnostic_D2(line, beam_dist, currents):
    print(f"\n{'=' * 72}")
    print("  D2: Twiss computation — dispersion-corrected vs raw")
    print("=" * 72)

    try:
        from rftrackAdapter import RFTrackAdapter
    except ImportError:
        print("  RF-Track adapter not available — skipping D2")
        return

    # Apply currents
    for idx, current in currents.items():
        if idx < len(line):
            line[idx].current = abs(current)

    # Propagate with FELsim particles and sample at chicane midpoint and end
    ebeam_calc = beam()
    sim = RFTrackAdapter.__new__(RFTrackAdapter)
    sim.logger = __import__('logging').getLogger('diag')

    checkpoints = {
        'After chicane 3 (idx 55)': 55,    # inside chicane, high dispersion
        'Before stage 11 (idx 87)': 87,     # start of stage 11
        'Undulator entrance (idx 117)': 117  # final
    }

    for label, stop_idx in checkpoints.items():
        particles = beam_dist.copy()
        for i, elem in enumerate(line):
            particles = np.array(elem.useMatrice(particles))
            if i == stop_idx:
                break

        _, _, twiss_corrected = ebeam_calc.cal_twiss(particles, ddof=1)
        twiss_raw = sim._calculate_twiss(particles)

        disp_x = twiss_corrected.loc['x', r'$D$ (m)']
        print(f"\n  {label} (D_x = {disp_x:.4f} m):")
        print(f"  {'Param':<10} {'Corrected':>12} {'Raw':>12} {'Δ (%)':>10}")
        print(f"  {'-' * 48}")

        for param, plane, raw_key in [
            ('beta_x', 'x', 'beta'), ('beta_y', 'y', 'beta'),
            ('alpha_x', 'x', 'alpha'), ('alpha_y', 'y', 'alpha'),
        ]:
            twiss_col = {
                'beta_x': r'$\beta$ (m)', 'beta_y': r'$\beta$ (m)',
                'alpha_x': r'$\alpha$', 'alpha_y': r'$\alpha$',
            }[param]
            corr = twiss_corrected.loc[plane, twiss_col]
            raw = twiss_raw[plane][raw_key]
            delta_pct = (raw - corr) / abs(corr) * 100 if abs(corr) > 1e-6 else 0
            print(f"  {param:<10} {corr:12.4f} {raw:12.4f} {delta_pct:10.2f}")

    print(f"\n  NOTE: collect_evolution() uses _calculate_twiss (raw), not")
    print(f"  ebeam.cal_twiss (corrected). This affects beta(s) plots in")
    print(f"  chicane regions but NOT the optimization MSE.")


# ═══════════════════════════════════════════════════════════════════════════════
# D7: Total lattice path length comparison
# ═══════════════════════════════════════════════════════════════════════════════

def diagnostic_D7(line):
    print(f"\n{'=' * 72}")
    print("  D7: Total lattice path length comparison")
    print("=" * 72)

    try:
        from rftrackAdapter import RFTrackAdapter
    except ImportError:
        print("  RF-Track adapter not available — skipping D7")
        return

    # FELsim total length (beamline element lengths)
    felsim_len = sum(getattr(e, 'length', getattr(e, 'L', 0)) for e in line)

    # RF-Track: beamline element lengths vs internal lattice length
    sim = RFTrackAdapter(
        lattice_path=str(EXCEL_PATH),
        beam_energy=ENERGY,
        space_charge=False,
        aperture=0.5,
    )
    sim.beamline = sim.beamline[:SEGMENTS]
    sim._build_lattice()
    rftrack_beamline_len = sum(e.length for e in sim.beamline)
    rftrack_lattice_len = sim._lattice.get_length()

    print(f"  FELsim element sum:     {felsim_len:.6f} m")
    print(f"  RF-Track beamline sum:  {rftrack_beamline_len:.6f} m")
    print(f"  RF-Track _lattice len:  {rftrack_lattice_len:.6f} m")
    print(f"  Δ (beamline):           {rftrack_beamline_len - felsim_len:.6e} m")
    print(f"  Δ (_lattice):           {rftrack_lattice_len - felsim_len:.6e} m")

    # Count DPW elements and their length contribution
    n_dpw = sum(1 for e in line if type(e).__name__ == 'dipole_wedge')
    dpw_felsim_total = sum(e.length for e in line if type(e).__name__ == 'dipole_wedge')
    dpw_rft_total = n_dpw * sim.DPW_THIN_LENS_LENGTH
    print(f"\n  DPW thin-lens issue:")
    print(f"    DPW count:              {n_dpw}")
    print(f"    FELsim DPW length each: {line[6].length * 1e3:.3f} mm")
    print(f"    RF-Track DPW length:    {sim.DPW_THIN_LENS_LENGTH:.1e} m")
    print(f"    FELsim DPW total:       {dpw_felsim_total:.6f} m")
    print(f"    RF-Track DPW total:     {dpw_rft_total:.2e} m")
    print(f"    Missing drift:          {dpw_felsim_total - dpw_rft_total:.6f} m")
    print(f"    This explains {(dpw_felsim_total - dpw_rft_total)/(felsim_len - rftrack_lattice_len)*100:.0f}% "
          f"of the _lattice length discrepancy")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Phase 1 Diagnostics: RF-Track vs FELsim model alignment")
    print("=" * 72)

    line, relat = build_felsim_line()
    beam_dist = generate_beam(relat)

    # Load best FELsim currents (for D2)
    currents_path = Path(__file__).resolve().parent / 'results' / 'cosy_s1_fr3_postfix.json'
    if currents_path.exists():
        with open(currents_path) as f:
            data = json.load(f)
        currents = {int(k): v for k, v in data['currents'].items()}
    else:
        currents = {}

    # D5: asymmetric edge kicks
    dpw_data = diagnostic_D5(line)

    # D7: path lengths
    diagnostic_D7(line)

    # D2: Twiss methods (needs optimised currents for meaningful dispersion)
    if currents:
        diagnostic_D2(line, beam_dist, currents)
    else:
        print("\n  Skipping D2: no optimised currents found")

    # Save D5 data
    output_path = RESULTS_DIR / 'dpw_edge_kicks.json'
    with open(output_path, 'w') as f:
        json.dump(dpw_data, f, indent=2)
    print(f"\nD5 data saved: {output_path}")
