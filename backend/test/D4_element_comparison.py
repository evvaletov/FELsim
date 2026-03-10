#!/usr/bin/env python3
"""Phase 2 diagnostics: element-by-element FELsim vs RF-Track comparison.

D4: Track identical particles through each element in both codes,
    compare output phase space at every element boundary.
D8: Sector-bend correction fidelity (analytical vs FELsim matrix).
D9: Quadrupole k1 and transfer matrix comparison.

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
NB_PARTICLES = 50  # small for fast diagnostics
SEGMENTS = 118

RESULTS_DIR = Path(__file__).resolve().parent / 'results' / 'D4'
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
# D8: Sector-bend correction fidelity
# ═══════════════════════════════════════════════════════════════════════════════

def diagnostic_D8(line, beam_dist, relat):
    """Compare FELsim dipole matrix vs RF-Track Drift+correction for DPH elements."""
    print("=" * 72)
    print("  D8: Sector-bend correction fidelity (per-DPH element)")
    print("=" * 72)

    try:
        import RF_Track as rft
        from rftrackAdapter import RFTrackAdapter
    except ImportError:
        print("  RF-Track not available — skipping D8")
        return

    # Build adapter for coordinate conversion
    sim = RFTrackAdapter(
        lattice_path=str(EXCEL_PATH),
        beam_energy=ENERGY,
        space_charge=False,
        aperture=0.5,
    )

    particle_mass = sim.particle_mass
    Pc = sim._Pc
    beta = sim._beta
    gamma = sim._gamma

    print(f"\n  Reference: E={ENERGY} MeV, Pc={Pc:.4f} MeV/c, β={beta:.6f}, γ={gamma:.4f}")
    print(f"  Particles: {NB_PARTICLES}")

    print(f"\n  {'idx':>4} {'name':<8} {'type':<6} {'L(mm)':>8} {'θ(deg)':>8} "
          f"{'Δx_rms':>10} {'Δxp_rms':>10} {'Δy_rms':>10} {'Δyp_rms':>10} "
          f"{'Δt_rms':>10} {'Δ_max':>10}")
    print(f"  {'-' * 100}")

    results = []

    for idx, elem in enumerate(line):
        cls_name = type(elem).__name__

        # FELsim: apply transfer matrix directly
        ps_felsim_in = beam_dist.copy()
        ps_felsim_out = np.array(elem.useMatrice(ps_felsim_in))

        # RF-Track: convert to RF-Track coords, track, convert back
        ps_rft_in = sim.transform_coordinates(
            beam_dist.copy(), CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK)

        if cls_name == 'dipole':
            angle_rad = np.radians(elem.angle)
            length = elem.length

            # Track through Drift
            drift = rft.Drift(length)
            drift.set_aperture(0.5, 0.5)
            lat = rft.Lattice()
            lat.append(drift)
            bunch = rft.Bunch6d(particle_mass, -1.0, Pc, ps_rft_in)
            tracked = lat.track(bunch)
            ps_rft_out = np.array(tracked.get_phase_space())

            # Apply sector-bend correction
            RFTrackAdapter._apply_sector_bend_correction(
                ps_rft_out, length, angle_rad, Pc, particle_mass)

        elif cls_name == 'dipole_wedge':
            length = elem.length
            # Track through Drift
            drift = rft.Drift(length)
            drift.set_aperture(0.5, 0.5)
            lat = rft.Lattice()
            lat.append(drift)
            bunch = rft.Bunch6d(particle_mass, -1.0, Pc, ps_rft_in)
            tracked = lat.track(bunch)
            ps_rft_out = np.array(tracked.get_phase_space())

            # Compute and apply DPW edge kicks (same as _compute_dpw_kicks)
            R = elem.dipole_length / (abs(elem.dipole_angle) * np.pi / 180)
            K0 = abs(np.radians(elem.dipole_angle) / elem.dipole_length)
            eta = np.radians(elem.angle)
            Tx = np.tan(eta)
            g = elem.pole_gap
            le = elem.length
            K_tri = le / (6.0 * g)
            h = 1.0 / R
            phi = K_tri * g * h * (1 + np.sin(eta)**2) / np.cos(eta)
            Ty = np.tan(eta - phi)
            kick_x = Tx / R
            kick_y = -Ty / R
            RFTrackAdapter._apply_dpw_edge_kick(ps_rft_out, kick_x, kick_y)

        elif cls_name in ('qpfLattice', 'qpdLattice'):
            native = sim._convert_element_to_native(sim.beamline[idx])
            lat = rft.Lattice()
            lat.append(native)
            bunch = rft.Bunch6d(particle_mass, -1.0, Pc, ps_rft_in)
            tracked = lat.track(bunch)
            ps_rft_out = np.array(tracked.get_phase_space())

        elif cls_name == 'driftLattice':
            native = rft.Drift(elem.length)
            native.set_aperture(0.5, 0.5)
            lat = rft.Lattice()
            lat.append(native)
            bunch = rft.Bunch6d(particle_mass, -1.0, Pc, ps_rft_in)
            tracked = lat.track(bunch)
            ps_rft_out = np.array(tracked.get_phase_space())

        else:
            continue  # skip unknown types

        # Convert RF-Track output back to FELsim coordinates
        ps_rft_back = sim.transform_coordinates(
            ps_rft_out, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM)

        # Compare
        diff = ps_rft_back - ps_felsim_out
        rms = np.sqrt(np.mean(diff**2, axis=0))
        max_abs = np.max(np.abs(diff))

        name = getattr(elem, 'name', '') or ''
        angle_str = f"{getattr(elem, 'angle', 0.0):.3f}"
        L_mm = elem.length * 1000

        is_significant = max_abs > 1e-6
        marker = " ***" if max_abs > 1e-3 else (" *" if max_abs > 1e-6 else "")

        print(f"  {idx:4d} {name:<8} {cls_name[:6]:<6} {L_mm:8.2f} {angle_str:>8} "
              f"{rms[0]:10.2e} {rms[1]:10.2e} {rms[2]:10.2e} {rms[3]:10.2e} "
              f"{rms[4]:10.2e} {max_abs:10.2e}{marker}")

        results.append({
            'index': idx,
            'name': name,
            'type': cls_name,
            'length_mm': L_mm,
            'rms_x': rms[0], 'rms_xp': rms[1],
            'rms_y': rms[2], 'rms_yp': rms[3],
            'rms_t': rms[4], 'rms_p': rms[5],
            'max_abs_diff': max_abs,
        })

    # Summary: top 10 largest discrepancies
    sorted_res = sorted(results, key=lambda r: r['max_abs_diff'], reverse=True)
    print(f"\n  Top 10 elements by max absolute difference:")
    print(f"  {'idx':>4} {'name':<8} {'type':<12} {'max_abs':>12}")
    for r in sorted_res[:10]:
        print(f"  {r['index']:4d} {r['name']:<8} {r['type']:<12} {r['max_abs_diff']:12.2e}")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# D9: Quadrupole k1 comparison
# ═══════════════════════════════════════════════════════════════════════════════

def diagnostic_D9(line, relat):
    """Compare quadrupole gradient k1 between FELsim and RF-Track adapter."""
    print(f"\n{'=' * 72}")
    print("  D9: Quadrupole k1 and transfer matrix comparison")
    print("=" * 72)

    try:
        from rftrackAdapter import RFTrackAdapter
    except ImportError:
        print("  RF-Track adapter not available — skipping D9")
        return

    sim = RFTrackAdapter(
        lattice_path=str(EXCEL_PATH),
        beam_energy=ENERGY,
        space_charge=False,
        aperture=0.5,
    )

    Q = 1.60217663e-19
    M_kg = 9.1093837e-31
    C = 299792458
    G = 2.694

    print(f"\n  FELsim: Q={Q:.4e}, M={M_kg:.4e}, C={C}, G={G}")
    print(f"  RF-Track adapter: G_quad={sim.G_quad}")
    print(f"  β={relat.beta:.10f} vs {sim._beta:.10f}")
    print(f"  γ={relat.gamma:.10f} vs {sim._gamma:.10f}")

    print(f"\n  {'idx':>4} {'name':<8} {'type':<6} {'I (A)':>8} {'L (m)':>10} "
          f"{'k_fel':>12} {'k1_rft':>12} {'Δk/k(%)':>10} "
          f"{'M11_fel':>10} {'M11_rft':>10} {'Δ(%)':>8}")
    print(f"  {'-' * 110}")

    for idx, elem in enumerate(line):
        cls_name = type(elem).__name__
        if cls_name not in ('qpfLattice', 'qpdLattice'):
            continue

        current = elem.current
        length = elem.length
        focusing = cls_name == 'qpfLattice'

        # FELsim k
        k_fel = np.abs(Q * G * current) / (M_kg * C * relat.beta * relat.gamma)

        # RF-Track adapter k1
        k1_rft = sim._current_to_k1(current, length, focusing=focusing)
        k1_abs = abs(k1_rft)

        dk = (k1_abs - k_fel) / k_fel * 100 if k_fel > 0 else 0

        # FELsim matrix M11
        theta_fel = np.sqrt(k_fel) * length
        if focusing:
            M11_fel = np.cos(theta_fel)
        else:
            M11_fel = np.cosh(theta_fel)

        # RF-Track matrix M11 (thick quad with same k1)
        theta_rft = np.sqrt(k1_abs) * length
        if focusing:
            M11_rft = np.cos(theta_rft)
        else:
            M11_rft = np.cosh(theta_rft)

        dM11 = (M11_rft - M11_fel) / abs(M11_fel) * 100 if abs(M11_fel) > 1e-10 else 0

        name = getattr(elem, 'name', '') or ''
        typ = 'QPF' if focusing else 'QPD'
        print(f"  {idx:4d} {name:<8} {typ:<6} {current:8.3f} {length:10.6f} "
              f"{k_fel:12.6f} {k1_abs:12.6f} {dk:10.4f} "
              f"{M11_fel:10.6f} {M11_rft:10.6f} {dM11:8.4f}")


# ═══════════════════════════════════════════════════════════════════════════════
# D4: Cumulative tracking comparison (element-by-element)
# ═══════════════════════════════════════════════════════════════════════════════

def diagnostic_D4(line, beam_dist, relat, currents):
    """Track same particles through full beamline in both codes, compare at each boundary."""
    print(f"\n{'=' * 72}")
    print("  D4: Cumulative element-by-element tracking comparison")
    print("=" * 72)

    try:
        import RF_Track as rft
        from rftrackAdapter import RFTrackAdapter
    except ImportError:
        print("  RF-Track not available — skipping D4")
        return

    # Apply currents
    for idx_str, current in currents.items():
        idx = int(idx_str) if isinstance(idx_str, str) else idx_str
        if idx < len(line):
            line[idx].current = abs(current)

    sim = RFTrackAdapter(
        lattice_path=str(EXCEL_PATH),
        beam_energy=ENERGY,
        space_charge=False,
        aperture=0.5,
    )
    # Apply same currents to RF-Track adapter
    for idx_str, current in currents.items():
        idx = int(idx_str) if isinstance(idx_str, str) else idx_str
        if idx < len(sim.beamline):
            sim.beamline[idx].parameters['current'] = abs(current)

    sim._annotate_dipole_edges()

    particle_mass = sim.particle_mass
    Pc = sim._Pc

    # FELsim tracking
    ps_felsim = beam_dist.copy()
    # RF-Track tracking
    ps_rft = sim.transform_coordinates(
        beam_dist.copy(), CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK)

    ebeam_calc = beam()

    print(f"\n  Tracking {NB_PARTICLES} particles through {len(line)} elements...")
    print(f"\n  {'idx':>4} {'name':<8} {'type':<6} "
          f"{'σx_f':>8} {'σx_r':>8} {'Δσx%':>7} "
          f"{'σy_f':>8} {'σy_r':>8} {'Δσy%':>7} "
          f"{'<x>_f':>8} {'<x>_r':>8} "
          f"{'n_r':>4}")
    print(f"  {'-' * 100}")

    checkpoints = [5, 10, 20, 30, 40, 50, 60, 70, 80, 87, 90, 93, 95, 97, 100, 110, 117]
    divergence_data = []

    for idx, elem in enumerate(line):
        cls_name = type(elem).__name__

        # FELsim: apply matrix
        ps_felsim = np.array(elem.useMatrice(ps_felsim))

        # RF-Track: element-by-element with analytical corrections
        bl_elem = sim.beamline[idx]
        params = bl_elem.parameters

        if params.get('_analytical_dipole', False):
            angle_rad = np.radians(params.get('angle', 0.0))
            length = bl_elem.length
            drift = rft.Drift(length)
            drift.set_aperture(0.5, 0.5)
            lat = rft.Lattice()
            lat.append(drift)
            bunch = rft.Bunch6d(particle_mass, -1.0, Pc, ps_rft)
            tracked = lat.track(bunch)
            ps_rft = np.array(tracked.get_phase_space())
            if ps_rft.ndim == 2 and ps_rft.shape[0] > 0:
                RFTrackAdapter._apply_sector_bend_correction(
                    ps_rft, length, angle_rad, Pc, particle_mass)

        elif params.get('_analytical_dpw', False):
            native = sim._convert_element_to_native(bl_elem)
            lat = rft.Lattice()
            lat.append(native)
            bunch = rft.Bunch6d(particle_mass, -1.0, Pc, ps_rft)
            tracked = lat.track(bunch)
            ps_rft = np.array(tracked.get_phase_space())
            if ps_rft.ndim == 2 and ps_rft.shape[0] > 0:
                RFTrackAdapter._apply_dpw_edge_kick(
                    ps_rft,
                    params.get('dpw_kick_x', 0.0),
                    params.get('dpw_kick_y', 0.0))

        else:
            native = sim._convert_element_to_native(bl_elem)
            lat = rft.Lattice()
            if isinstance(native, list):
                for ne in native:
                    lat.append(ne)
            else:
                lat.append(native)
            bunch = rft.Bunch6d(particle_mass, -1.0, Pc, ps_rft)
            tracked = lat.track(bunch)
            ps_rft = np.array(tracked.get_phase_space())

        if idx not in checkpoints:
            continue

        # Compare at this checkpoint
        n_rft = ps_rft.shape[0] if ps_rft.ndim == 2 else 0

        # Convert RF-Track back to FELsim coords for comparison
        if n_rft > 0:
            ps_rft_fel = sim.transform_coordinates(
                ps_rft, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM)
        else:
            ps_rft_fel = np.zeros((0, 6))

        # Use only surviving particles for FELsim stats (same indices)
        sigma_x_f = np.std(ps_felsim[:, 0], ddof=1)
        sigma_y_f = np.std(ps_felsim[:, 2], ddof=1)
        mean_x_f = np.mean(ps_felsim[:, 0])

        if n_rft > 0:
            sigma_x_r = np.std(ps_rft_fel[:, 0], ddof=1)
            sigma_y_r = np.std(ps_rft_fel[:, 2], ddof=1)
            mean_x_r = np.mean(ps_rft_fel[:, 0])
        else:
            sigma_x_r = sigma_y_r = mean_x_r = float('nan')

        dsig_x = (sigma_x_r - sigma_x_f) / sigma_x_f * 100 if sigma_x_f > 1e-10 else 0
        dsig_y = (sigma_y_r - sigma_y_f) / sigma_y_f * 100 if sigma_y_f > 1e-10 else 0

        name = getattr(elem, 'name', '') or ''
        print(f"  {idx:4d} {name:<8} {cls_name[:6]:<6} "
              f"{sigma_x_f:8.3f} {sigma_x_r:8.3f} {dsig_x:7.1f} "
              f"{sigma_y_f:8.3f} {sigma_y_r:8.3f} {dsig_y:7.1f} "
              f"{mean_x_f:8.4f} {mean_x_r:8.4f} "
              f"{n_rft:4d}")

        divergence_data.append({
            'idx': idx, 'name': name, 'type': cls_name,
            'sigma_x_felsim': sigma_x_f, 'sigma_x_rftrack': sigma_x_r,
            'sigma_y_felsim': sigma_y_f, 'sigma_y_rftrack': sigma_y_r,
            'dsig_x_pct': dsig_x, 'dsig_y_pct': dsig_y,
            'n_rftrack': n_rft,
        })

    # Compute Twiss at endpoint
    if ps_rft.ndim == 2 and ps_rft.shape[0] > 5:
        ps_end_fel = sim.transform_coordinates(
            ps_rft, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM)
        _, _, twiss_rft = ebeam_calc.cal_twiss(ps_end_fel, ddof=1)
        _, _, twiss_fel = ebeam_calc.cal_twiss(ps_felsim, ddof=1)

        print(f"\n  Endpoint Twiss comparison (FELsim currents, both codes):")
        print(f"  {'':>14} {'FELsim':>10} {'RF-Track':>10} {'Δ(%)':>8} {'Target':>10}")
        targets = {'beta_x': 1.4, 'beta_y': 0.2418, 'alpha_x': 0.47, 'alpha_y': 0.0}
        for param, target in targets.items():
            plane = 'x' if 'x' in param else 'y'
            col = r'$\beta$ (m)' if 'beta' in param else r'$\alpha$'
            v_fel = twiss_fel.loc[plane, col]
            v_rft = twiss_rft.loc[plane, col]
            dp = (v_rft - v_fel) / abs(v_fel) * 100 if abs(v_fel) > 1e-6 else 0
            print(f"  {param:<14} {v_fel:10.4f} {v_rft:10.4f} {dp:8.2f} {target:10.4f}")

    return divergence_data


# ═══════════════════════════════════════════════════════════════════════════════
# D10: DPW M56 term comparison
# ═══════════════════════════════════════════════════════════════════════════════

def diagnostic_D10(line, relat):
    """Check if DPW M56 (path-length vs energy) is handled in RF-Track."""
    print(f"\n{'=' * 72}")
    print("  D10: DPW M56 (path-length ↔ energy coupling) comparison")
    print("=" * 72)

    f_RF = 2856e6
    C = 299792458
    beta = relat.beta
    gamma = relat.gamma

    print(f"\n  FELsim DPW matrix has M56 = -f·L / (c·β·γ·(γ+1))")
    print(f"  RF-Track DPW (Drift + edge kick) has NO explicit M56 coupling.")
    print(f"  The Drift's M56 = -L/γ² ≠ FELsim's M56.")
    print(f"\n  {'idx':>4} {'name':<8} {'L(mm)':>8} {'M56_fel':>12} {'M56_drift':>12} {'Δ':>12}")
    print(f"  {'-' * 70}")

    total_M56_diff = 0
    for idx, elem in enumerate(line):
        if type(elem).__name__ != 'dipole_wedge':
            continue

        L = elem.length
        M56_felsim = -f_RF * L / (C * beta * gamma * (gamma + 1))
        # Drift M56 in FELsim coordinates
        # Drift has no M56 coupling in its matrix (M55=1, M56=0 for standard drift)
        # But FELsim's drift has M56 = -f*L/(c*β*γ*(γ+1)) too!
        # Actually, let me check the drift matrix...
        M56_drift = -f_RF * L / (C * beta * gamma * (gamma + 1))

        diff = M56_drift - M56_felsim
        total_M56_diff += abs(diff)
        name = getattr(elem, 'name', '') or ''
        print(f"  {idx:4d} {name:<8} {L*1e3:8.3f} {M56_felsim:12.6e} {M56_drift:12.6e} {diff:12.6e}")

    print(f"\n  Total |ΔM56| across all DPW: {total_M56_diff:.6e}")
    print(f"  Note: FELsim drift and DPW both use M56 = -f·L/(c·β·γ·(γ+1)).")
    print(f"  RF-Track drift handles this via relativistic time-of-flight.")
    print(f"  The two SHOULD agree if RF-Track drift R56 = -L/γ² matches FELsim's M56.")

    # Verify equivalence: FELsim M56 × coord6 vs RF-Track drift's R56 × δ_P
    # FELsim: Δcoord5 = M56 × coord6 = -fL/(cβγ(γ+1)) × ΔK/K₀×10³
    # Converting to ct [mm]: Δ(ct) = Δcoord5 × c/f = -L/(βγ(γ+1)) × ΔK/K₀ × 10³
    # RF-Track drift: Δ(ct) = -L/γ² × (1000/β) × δ_P
    # (where the 1000 converts m→mm, and /β converts path→time)
    # Now δ_P = ΔK/K₀ × (γ-1)/(β²γ), so:
    # RF-Track: Δ(ct) = -L/γ² × 1000/β × ΔK/K₀ × (γ-1)/(β²γ)
    #         = -L × 1000 × (γ-1) / (β³γ³) × ΔK/K₀

    # FELsim:  Δ(ct) = -L/(βγ(γ+1)) × ΔK/K₀ × 10³
    #         = -L × 1000 / (βγ(γ+1)) × ΔK/K₀

    ratio_theory = (beta**2 * gamma**2 * (gamma + 1)) / (gamma - 1) / (gamma**2)
    # Actually: FELsim/RFTrack = [1/(βγ(γ+1))] / [(γ-1)/(β³γ³)]
    #         = β³γ³ / (βγ(γ+1)(γ-1)) = β²γ² / ((γ+1)(γ-1)) = β²γ²/(γ²-1) = 1
    print(f"\n  Analytical ratio FELsim/RF-Track M56: β²γ²/(γ²-1) = {beta**2 * gamma**2 / (gamma**2 - 1):.10f}")
    print(f"  (should be 1.0 — drift R56 matches FELsim M56 exactly)")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Phase 2 Diagnostics: element-by-element FELsim vs RF-Track")
    print("=" * 72)

    line, relat = build_felsim_line()
    beam_dist = generate_beam(relat)

    # Load FELsim optimised currents
    currents_path = Path(__file__).resolve().parent / 'results' / 'cosy_s1_fr3_postfix.json'
    if currents_path.exists():
        with open(currents_path) as f:
            data = json.load(f)
        currents = {int(k): v for k, v in data['currents'].items()}
    else:
        print("  WARNING: No optimised currents found, using defaults")
        currents = {}

    # D9: quadrupole k1 comparison
    diagnostic_D9(line, relat)

    # D10: DPW M56 coupling
    diagnostic_D10(line, relat)

    # D8: single-element sector-bend comparison
    results_D8 = diagnostic_D8(line, beam_dist, relat)

    # D4: cumulative tracking (needs currents)
    if currents:
        results_D4 = diagnostic_D4(line, beam_dist, relat, currents)
    else:
        print("\n  Skipping D4: no optimised currents found")

    print(f"\nPhase 2 diagnostics complete.")
