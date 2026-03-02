"""C6 QA: Comprehensive post-fix verification for coordinate convention (ΔK/K₀),
fringeField M56 sign, R56 β₀ factor, and sector-bend corrections.

Tier 1 (T1.1–T1.13): Unit-level physics checks (always runs; inline fallback
    for coordinate transforms if RF-Track unavailable)
Tier 2 (T2.1–T2.4): Integration tests (RF-Track tracking required)

Author: Eremey Valetov
"""

import sys
import os
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from physicalConstants import PhysicalConstants
from beamline import driftLattice, qpfLattice, qpdLattice, dipole, dipole_wedge, beamline

# RF-Track gated imports
try:
    from rftrackAdapter import RFTrackAdapter, _RFTRACK_AVAILABLE, CoordinateSystem
    HAS_RFTRACK = _RFTRACK_AVAILABLE
except ImportError:
    HAS_RFTRACK = False

SEED = 42
PC = PhysicalConstants


def _relativistic(KE, E0):
    gamma = 1 + KE / E0
    beta = np.sqrt(1 - 1 / gamma**2)
    return gamma, beta


# =========================================================================
# Inline coordinate transforms (fallback when RF-Track unavailable)
# =========================================================================

def _felsim_to_rftrack(p, KE, m):
    r = p.copy()
    K = KE * (1.0 + p[:, 5] * 1e-3)
    E = K + m
    r[:, 5] = np.sqrt(E**2 - m**2)
    return r


def _rftrack_to_felsim(p, KE, m):
    r = p.copy()
    K = np.sqrt(p[:, 5]**2 + m**2) - m
    r[:, 5] = (K / KE - 1.0) * 1e3
    return r


def _cosy_to_rftrack(p, KE, m):
    gamma, beta = _relativistic(KE, m)
    r = np.zeros_like(p)
    r[:, 0:4] = p[:, 0:4] * 1e3
    r[:, 4] = p[:, 4] / (beta * PC.C) * 1e3
    E0 = KE + m
    E = E0 + KE * p[:, 5]
    r[:, 5] = np.sqrt(E**2 - m**2)
    return r


def _rftrack_to_cosy(p, KE, m):
    gamma, beta = _relativistic(KE, m)
    r = np.zeros_like(p)
    r[:, 0:4] = p[:, 0:4] * 1e-3
    r[:, 4] = p[:, 4] * (beta * PC.C) * 1e-3
    K = np.sqrt(p[:, 5]**2 + m**2) - m
    r[:, 5] = K / KE - 1.0
    return r


# =========================================================================
# Tier 1 — Unit-level physics
# =========================================================================

def test_T1_1_m56_sign():
    """All 6 element types produce M56 < 0."""
    E = 45.0
    elements = [
        ("driftLattice",  driftLattice(0.5)),
        ("qpfLattice",    qpfLattice(2.0, 0.1)),
        ("qpdLattice",    qpdLattice(2.0, 0.1)),
        ("dipole",        dipole(0.2, 15.0)),
        ("dipole_wedge",  dipole_wedge(0.01, 5.0, 0.2, 15.0)),
        ("fringeField",   beamline.fringeField(0.01, 0.5)),
    ]
    for _, elem in elements:
        elem.setE(E)

    for name, elem in elements:
        M = elem._compute_numeric_matrix()
        m56 = M[4, 5]
        if m56 >= 0:
            return False, f"{name}: M56 = {m56:.6e} >= 0"
    return True, "All 6 element types: M56 < 0"


def test_T1_2_fringe_symbolic_numeric():
    """fringeField symbolic and numeric matrices agree to < 1e-14."""
    ff = beamline.fringeField(0.015, 0.3)
    ff.setE(45.0)
    M_num = ff._compute_numeric_matrix()
    M_sym = np.array(ff._compute_symbolic_matrix().tolist(), dtype=float)
    err = np.max(np.abs(M_num - M_sym))
    return err < 1e-14, f"max |numeric - symbolic| = {err:.2e}"


def test_T1_3_m56_quantitative():
    """M56 matches -(Lf)/(Cbgamma(gamma+1)) to < 1e-10 relative."""
    E = 45.0

    max_rel = 0
    for L in [0.05, 0.1, 0.5, 1.0]:
        d = driftLattice(L)
        d.setE(E)
        # Use the element's own constants (beamline.py E0 differs slightly from CODATA)
        expected = -(L * d.f) / (d.C * d.beta * d.gamma * (d.gamma + 1))
        M = d._compute_numeric_matrix()
        m56 = M[4, 5]
        rel_err = abs((m56 - expected) / expected)
        max_rel = max(max_rel, rel_err)

    return max_rel < 1e-10, f"max relative error = {max_rel:.2e}"


# --- T1.4–T1.8: Coordinate transforms ---

def test_T1_4_felsim_rftrack_roundtrip():
    """FELsim<->RF-Track round-trip error < 1e-11 at multiple energies."""
    np.random.seed(SEED)
    m = PC.E0_electron
    particles = np.random.normal(0, [1.0, 0.1, 1.0, 0.1, 0.5, 0.3], (100, 6))

    max_err = 0
    for KE in [10.0, 45.0, 100.0, 500.0]:
        if HAS_RFTRACK:
            sim = RFTrackAdapter(beam_energy=KE)
            rft = sim.transform_coordinates(particles, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK)
            back = sim.transform_coordinates(rft, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM)
        else:
            rft = _felsim_to_rftrack(particles, KE, m)
            back = _rftrack_to_felsim(rft, KE, m)
        err = np.max(np.abs(particles - back))
        max_err = max(max_err, err)

    return max_err < 1e-11, f"max round-trip error = {max_err:.2e}"


def test_T1_5_cosy_rftrack_roundtrip():
    """COSY<->RF-Track delta round-trip error < 1e-12."""
    np.random.seed(SEED)
    m = PC.E0_electron
    KE = 45.0
    particles_cosy = np.random.normal(0, [1e-3, 1e-4, 1e-3, 1e-4, 1e-4, 1e-3], (100, 6))

    if HAS_RFTRACK:
        sim = RFTrackAdapter(beam_energy=KE)
        rft = sim.transform_coordinates(particles_cosy, CoordinateSystem.COSY, CoordinateSystem.RFTRACK)
        back = sim.transform_coordinates(rft, CoordinateSystem.RFTRACK, CoordinateSystem.COSY)
    else:
        rft = _cosy_to_rftrack(particles_cosy, KE, m)
        back = _rftrack_to_cosy(rft, KE, m)

    err = np.max(np.abs(particles_cosy - back))
    return err < 1e-12, f"max round-trip error = {err:.2e}"


def test_T1_6_four_hop_chain():
    """FELSIM->RFTRACK->COSY->RFTRACK->FELSIM error < 1e-10."""
    np.random.seed(SEED)
    m = PC.E0_electron
    KE = 45.0
    p0 = np.random.normal(0, [1.0, 0.1, 1.0, 0.1, 0.5, 0.3], (100, 6))

    if HAS_RFTRACK:
        sim = RFTrackAdapter(beam_energy=KE)
        p1 = sim.transform_coordinates(p0, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK)
        p2 = sim.transform_coordinates(p1, CoordinateSystem.RFTRACK, CoordinateSystem.COSY)
        p3 = sim.transform_coordinates(p2, CoordinateSystem.COSY, CoordinateSystem.RFTRACK)
        p4 = sim.transform_coordinates(p3, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM)
    else:
        p1 = _felsim_to_rftrack(p0, KE, m)
        p2 = _rftrack_to_cosy(p1, KE, m)
        p3 = _cosy_to_rftrack(p2, KE, m)
        p4 = _rftrack_to_felsim(p3, KE, m)

    err = np.max(np.abs(p0 - p4))
    return err < 1e-10, f"4-hop chain error = {err:.2e}"


def test_T1_7_cosy_delta_felsim_coord6():
    """COSY delta = FELsim coord6/1000 (DK/K0 consistency)."""
    m = PC.E0_electron
    KE = 45.0

    coord6_values = np.array([1.0, 2.5, -0.5, 5.0, -3.0])
    delta_expected = coord6_values * 1e-3

    N = len(coord6_values)
    p_felsim = np.zeros((N, 6))
    p_felsim[:, 5] = coord6_values

    p_rft = _felsim_to_rftrack(p_felsim, KE, m)
    p_cosy = _rftrack_to_cosy(p_rft, KE, m)

    err = np.max(np.abs(p_cosy[:, 5] - delta_expected))
    return err < 1e-14, f"max |delta_COSY - coord6/1000| = {err:.2e}"


def test_T1_8_old_vs_new_convention():
    """Document error of confusing DK/K0 with DE/E0 (total energy)."""
    # If coord6 = DK/K0 but interpreted as DE/E0 (E = total energy):
    # error factor = m/K0 (rest mass / kinetic energy)
    #
    # 40 MeV electron: 0.511/40 = 1.28%
    # 250 MeV proton: 938.27/250 = 375%

    err_e = PC.E0_electron / 40.0 * 100
    err_p = PC.E0_proton / 250.0 * 100

    details = (f"40 MeV e-: {err_e:.2f}% error (expect ~1.28%), "
               f"250 MeV p: {err_p:.1f}% error (expect ~375%)")

    passed = (abs(err_e - 1.28) < 0.05 and abs(err_p - 375) < 2)
    return passed, details


# --- T1.9–T1.13: Analytical sector-bend physics ---

def test_T1_9_r56_sector_bend():
    """Dipole R56 matches -(rho*sin(theta) - L) formula at 5 angles."""
    E = 45.0
    L = 0.2

    max_rel = 0
    for angle_deg in [5.0, 10.0, 15.0, 20.0, 45.0]:
        d = dipole(L, angle_deg)
        d.setE(E)
        M = d._compute_numeric_matrix()
        m56 = M[4, 5]

        # Use element's own constants for consistency
        theta = angle_deg * np.pi / 180
        rho = L / theta
        S = np.sin(theta)
        expected = -d.f * (L - rho * S) / (d.C * d.beta * d.gamma * (d.gamma + 1))

        if abs(expected) > 0:
            rel_err = abs((m56 - expected) / expected)
            max_rel = max(max_rel, rel_err)

    return max_rel < 1e-12, f"max relative error = {max_rel:.2e}"


def test_T1_10_r56_beta_factor():
    """1/beta0 - 1: negligible for electrons, ~63% for 250 MeV proton."""
    _, beta_e = _relativistic(45.0, PC.E0_electron)
    factor_e = (1 / beta_e - 1) * 100

    _, beta_p = _relativistic(250.0, PC.E0_proton)
    factor_p = (1 / beta_p - 1) * 100

    details = (f"e-(45 MeV): 1/beta0-1 = {factor_e:.4f}%, "
               f"p(250 MeV): 1/beta0-1 = {factor_p:.1f}%")

    passed = (factor_e < 0.01 and factor_p > 50)
    return passed, details


def test_T1_11_dispersion_correction():
    """Dispersion terms rho*(1-cos theta)*delta and sin(theta)*delta agree with 2sin^2(theta/2) form."""
    L = 0.2
    theta = 15.0 * np.pi / 180
    rho = L / theta
    delta = 0.01

    dx_1mc = rho * (1 - np.cos(theta)) * delta
    dx_2s2 = rho * 2 * np.sin(theta / 2)**2 * delta
    dxp = np.sin(theta) * delta

    err = abs(dx_1mc - dx_2s2)
    details = (f"dx = {dx_1mc * 1000:.6f} mm, "
               f"dxp = {dxp * 1000:.6f} mrad, "
               f"|1-cos vs 2sin^2| = {err:.2e}")
    return err < 1e-16, details


def test_T1_12_trig_stability():
    """2sin^2(theta/2) more stable than 1-cos(theta) for small angles."""
    theta = 1e-8
    ref = theta**2 / 2  # exact to machine precision for theta << 1

    val_1mc = 1 - np.cos(theta)
    val_2s2 = 2 * np.sin(theta / 2)**2

    err_1mc = abs(val_1mc - ref) / ref
    err_2s2 = abs(val_2s2 - ref) / ref

    details = (f"theta={theta}: "
               f"1-cos rel err = {err_1mc:.2e}, "
               f"2sin^2 rel err = {err_2s2:.2e}")
    passed = err_2s2 <= err_1mc
    return passed, details


def test_T1_13_correction_det():
    """Correction matrix M_corr = M_sector * M_drift^-1 has det = 1 at 8 angles."""
    L = 0.2
    max_det_err = 0

    for angle_deg in [5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0, 90.0]:
        theta = angle_deg * np.pi / 180
        rho = L / theta

        C = np.cos(theta)
        S = np.sin(theta)

        # M_corr = M_sector @ inv(M_drift) for 2x2 horizontal submatrix
        R11 = C
        R12 = rho * S - L * C
        R21 = -S / rho
        R22 = L * S / rho + C

        det = R11 * R22 - R12 * R21
        max_det_err = max(max_det_err, abs(det - 1.0))

    return max_det_err < 1e-14, f"max |det - 1| = {max_det_err:.2e}"


# =========================================================================
# Tier 2 — Integration (RF-Track tracking required)
# =========================================================================

def test_T2_1_drift_velocity_dispersion():
    """RF-Track Drift includes 1/gamma^2 velocity dispersion."""
    if not HAS_RFTRACK:
        return None, "RF-Track not available"

    import RF_Track as rft

    KE = 45.0
    m = PC.E0_electron
    gamma, beta = _relativistic(KE, m)
    Pc = PC.momentum(KE, m)
    L = 1.0

    N = 11
    delta_values = np.linspace(-0.005, 0.005, N)
    ps = np.zeros((N, 6))
    ps[:, 5] = Pc * (1 + delta_values)

    lat = rft.Lattice()
    d = rft.Drift(L)
    d.set_aperture(0.05, 0.05)
    lat.append(d)

    bunch = rft.Bunch6d(m, -1.0, Pc, ps)
    bunch_out = lat.track(bunch)
    ps_out = np.array(bunch_out.get_phase_space())

    delta = (ps[:, 5] - Pc) / Pc
    dt = ps_out[:, 4] - ps[:, 4]  # mm/c

    R56_tracked = np.polyfit(delta, dt, 1)[0]
    # Expected: -L/(beta^2 * gamma^2 * beta) in mm/c per unit delta
    R56_expected = -L / (beta**2 * gamma**2 * beta) * 1000

    rel_err = abs(R56_tracked - R56_expected) / abs(R56_expected)
    details = (f"R56_tracked = {R56_tracked:.4f} mm/c, "
               f"R56_expected = {R56_expected:.4f} mm/c, "
               f"rel err = {rel_err:.2e}")
    return rel_err < 0.01, details


def test_T2_2_dipole_r56():
    """Dipole R56 (drift + correction) within 5% of analytical."""
    if not HAS_RFTRACK:
        return None, "RF-Track not available"

    import RF_Track as rft

    KE = 45.0
    m = PC.E0_electron
    gamma, beta = _relativistic(KE, m)
    Pc = PC.momentum(KE, m)
    L = 0.2
    angle_deg = 15.0
    theta = angle_deg * np.pi / 180
    rho = L / theta

    N = 201
    np.random.seed(SEED)
    delta_values = np.linspace(-0.003, 0.003, N)
    ps = np.zeros((N, 6))
    ps[:, 5] = Pc * (1 + delta_values)

    # Track through drift
    lat = rft.Lattice()
    d = rft.Drift(L)
    d.set_aperture(0.05, 0.05)
    lat.append(d)

    bunch = rft.Bunch6d(m, -1.0, Pc, ps)
    bunch_out = lat.track(bunch)
    ps_out = np.array(bunch_out.get_phase_space())

    # Apply sector-bend correction
    ps_corrected = RFTrackAdapter._apply_sector_bend_correction(
        ps_out.copy(), L, theta, Pc, m
    )

    # Extract R56 from corrected particles (mm/c per Dp/p)
    delta = (ps[:, 5] - Pc) / Pc
    dt = ps_corrected[:, 4] - ps[:, 4]
    R56_tracked = np.polyfit(delta, dt, 1)[0]

    # Expected total R56 in mm/c per Dp/p:
    # velocity: -L/(beta^2 * gamma^2 * beta) * 1000
    # geometric (correction): -(rho*sin(theta) - L) / beta * 1000
    R56_vel = -L / (beta**2 * gamma**2 * beta) * 1000
    R56_geo = -(rho * np.sin(theta) - L) / beta * 1000
    R56_expected = R56_vel + R56_geo

    rel_err = abs(R56_tracked - R56_expected) / abs(R56_expected)
    details = (f"R56_tracked = {R56_tracked:.2f} mm/c, "
               f"R56_expected = {R56_expected:.2f} mm/c, "
               f"rel err = {rel_err:.2%}")
    return rel_err < 0.05, details


def test_T2_3_full_beamline_r56():
    """Full beamline R56: compound RF-Track tracking within 10% of analytical."""
    if not HAS_RFTRACK:
        return None, "RF-Track not available"

    import RF_Track as rft

    KE = 45.0
    m = PC.E0_electron
    gamma, beta = _relativistic(KE, m)
    Pc = PC.momentum(KE, m)

    L_drift = 0.5
    L_dipole = 0.2
    angle_deg = 15.0
    theta = angle_deg * np.pi / 180
    rho = L_dipole / theta

    N = 201
    delta_values = np.linspace(-0.003, 0.003, N)
    ps = np.zeros((N, 6))
    ps[:, 5] = Pc * (1 + delta_values)

    def _make_drift_lattice(length):
        d = rft.Drift(length)
        d.set_aperture(0.05, 0.05)
        lat = rft.Lattice()
        lat.append(d)
        return lat

    # Track compound: drift1 → dipole(as drift+correction) → drift2
    bunch = rft.Bunch6d(m, -1.0, Pc, ps)
    ps1 = np.array(_make_drift_lattice(L_drift).track(bunch).get_phase_space())

    bunch = rft.Bunch6d(m, -1.0, Pc, ps1)
    ps2 = np.array(_make_drift_lattice(L_dipole).track(bunch).get_phase_space())
    ps2 = RFTrackAdapter._apply_sector_bend_correction(ps2, L_dipole, theta, Pc, m)

    bunch = rft.Bunch6d(m, -1.0, Pc, ps2)
    ps3 = np.array(_make_drift_lattice(L_drift).track(bunch).get_phase_space())

    # Extract R56 from tracked particles (mm/c per Δp/p)
    dt = ps3[:, 4] - ps[:, 4]
    R56_tracked = np.polyfit(delta_values, dt, 1)[0]

    # Analytical: velocity dispersion (all elements) + geometric (dipole correction)
    L_total = 2 * L_drift + L_dipole
    R56_vel = -L_total / (beta**2 * gamma**2 * beta) * 1000
    R56_geo = -(rho * np.sin(theta) - L_dipole) / beta * 1000
    R56_expected = R56_vel + R56_geo

    rel_err = abs(R56_tracked - R56_expected) / abs(R56_expected)
    details = (f"R56_tracked = {R56_tracked:.4f} mm/c, "
               f"R56_expected = {R56_expected:.4f} mm/c, "
               f"rel err = {rel_err:.2%}")
    return rel_err < 0.10, details


def test_T2_4_w8_smoke():
    """W8 smoke: 1-restart optimization, MSE < 0.1."""
    if not HAS_RFTRACK:
        return None, "RF-Track not available"

    try:
        sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
        from UHM_beamline_opt_05ps_params import run_optimization
    except ImportError as e:
        return None, f"Cannot import optimization params: {e}"

    result = run_optimization(epsilon_n=8, n_restarts=1, seed=SEED)
    mse = result['mse']
    details = f"MSE = {mse:.4f} (threshold: 0.1), nfev = {result['nfev']}"
    return mse < 0.1, details


# =========================================================================
# Test runner
# =========================================================================

ALL_TESTS = [
    ("T1.1",  "M56 sign consistency",         test_T1_1_m56_sign),
    ("T1.2",  "Symbolic vs numeric fringe",    test_T1_2_fringe_symbolic_numeric),
    ("T1.3",  "M56 quantitative values",       test_T1_3_m56_quantitative),
    ("T1.4",  "FELSIM<->RFTRACK round-trip",   test_T1_4_felsim_rftrack_roundtrip),
    ("T1.5",  "COSY<->RFTRACK round-trip",     test_T1_5_cosy_rftrack_roundtrip),
    ("T1.6",  "Full 4-hop chain",              test_T1_6_four_hop_chain),
    ("T1.7",  "COSY delta = coord6/1000",      test_T1_7_cosy_delta_felsim_coord6),
    ("T1.8",  "Old vs new convention size",     test_T1_8_old_vs_new_convention),
    ("T1.9",  "R56 sector-bend (5 angles)",    test_T1_9_r56_sector_bend),
    ("T1.10", "R56 beta0 factor",              test_T1_10_r56_beta_factor),
    ("T1.11", "Dispersion correction terms",   test_T1_11_dispersion_correction),
    ("T1.12", "2sin^2(theta/2) stability",     test_T1_12_trig_stability),
    ("T1.13", "Correction matrix det = 1",     test_T1_13_correction_det),
    ("T2.1",  "Drift velocity dispersion",     test_T2_1_drift_velocity_dispersion),
    ("T2.2",  "Dipole R56 (drift+correction)", test_T2_2_dipole_r56),
    ("T2.3",  "Full beamline R56",             test_T2_3_full_beamline_r56),
    ("T2.4",  "W8 smoke (1 restart, en=8)",    test_T2_4_w8_smoke),
]


def run_all():
    print("=" * 72)
    print("C6 QA: Comprehensive Post-Fix Verification")
    print("=" * 72)
    print(f"RF-Track available: {HAS_RFTRACK}")
    print()

    results = []
    t0 = time.time()

    for test_id, description, func in ALL_TESTS:
        try:
            passed, details = func()
        except Exception as e:
            passed, details = False, f"EXCEPTION: {e}"
            import traceback
            traceback.print_exc()

        if passed is None:
            status = "SKIP"
        elif passed:
            status = "PASS"
        else:
            status = "FAIL"

        results.append((test_id, description, status, details))
        print(f"  [{status}] {test_id}: {description}")
        print(f"         {details}")

    elapsed = time.time() - t0

    # Summary
    n_pass = sum(1 for _, _, s, _ in results if s == "PASS")
    n_fail = sum(1 for _, _, s, _ in results if s == "FAIL")
    n_skip = sum(1 for _, _, s, _ in results if s == "SKIP")
    n_total = len(results)

    print()
    print("=" * 72)
    print(f"Results: {n_pass} passed, {n_fail} failed, {n_skip} skipped "
          f"(of {n_total}) in {elapsed:.1f}s")
    print("=" * 72)

    if n_fail > 0:
        print("\nFailed tests:")
        for tid, desc, status, details in results:
            if status == "FAIL":
                print(f"  {tid}: {desc}")
                print(f"    {details}")

    return n_fail == 0


if __name__ == "__main__":
    success = run_all()
    sys.exit(0 if success else 1)
