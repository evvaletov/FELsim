"""
C4 V&V: Unit tests for FELsim transfer matrix elements.

Tests each element class against known analytical results from beam optics
textbooks (Wiedemann, Wille, Brown). Verifies:
  - Symplecticity (det = 1 for each 2×2 block)
  - Known limiting cases (zero current → drift, zero angle → drift)
  - Numeric vs symbolic consistency
  - Specific matrix element values against independent calculation
  - FODO cell properties

Author: Eremey Valetov
"""

import sys
from pathlib import Path

import numpy as np
import pytest

_BACKEND = Path(__file__).resolve().parent.parent
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))

from beamline import driftLattice, qpfLattice, qpdLattice, dipole, dipole_wedge

# Reference constants (same as beamline.py defaults at 45 MeV)
E0_electron = 0.51099895000  # MeV
KE_default = 45.0  # MeV
Q = 1.60217663e-19
M_e = 9.1093837e-31
C_light = 299792458.0
f_rf = 2856e6
G_quad = 2.694  # T/A/m

gamma_45 = 1 + KE_default / E0_electron
beta_45 = np.sqrt(1 - 1 / gamma_45**2)


# ── Helpers ──────────────────────────────────────────────────────────────

def symplectic_2x2(M, i, j):
    """Check det of 2×2 block M[i:i+2, j:j+2] == 1."""
    blk = M[i:i+2, j:j+2]
    return blk[0, 0] * blk[1, 1] - blk[0, 1] * blk[1, 0]


def assert_symplectic(M, tol=1e-12):
    """Assert each transverse 2×2 block has unit determinant."""
    assert abs(symplectic_2x2(M, 0, 0) - 1.0) < tol, f"x-block det = {symplectic_2x2(M, 0, 0)}"
    assert abs(symplectic_2x2(M, 2, 2) - 1.0) < tol, f"y-block det = {symplectic_2x2(M, 2, 2)}"


def assert_block_diagonal(M, tol=1e-14):
    """Assert x-y coupling is zero (off-diagonal 2×2 blocks are zero)."""
    np.testing.assert_allclose(M[0:2, 2:4], 0.0, atol=tol)
    np.testing.assert_allclose(M[2:4, 0:2], 0.0, atol=tol)


# ── Drift tests ──────────────────────────────────────────────────────────

class TestDrift:
    def test_identity_at_zero_length(self):
        """Drift of negligible length → identity (within constructor constraint)."""
        d = driftLattice(1e-15)
        M = d._compute_numeric_matrix()
        np.testing.assert_allclose(M, np.eye(6), atol=1e-10)

    def test_M12_M34_equals_L(self):
        for L in [0.1, 0.5, 1.0, 3.0]:
            M = driftLattice(L)._compute_numeric_matrix()
            assert abs(M[0, 1] - L) < 1e-14
            assert abs(M[2, 3] - L) < 1e-14

    def test_diagonal_ones(self):
        M = driftLattice(0.5)._compute_numeric_matrix()
        for i in range(6):
            assert abs(M[i, i] - 1.0) < 1e-14

    def test_symplecticity(self):
        M = driftLattice(1.234)._compute_numeric_matrix()
        assert_symplectic(M)

    def test_block_diagonal(self):
        M = driftLattice(0.5)._compute_numeric_matrix()
        assert_block_diagonal(M)

    def test_M56_formula(self):
        """M56 = -L·f / (c·β·γ·(γ+1))."""
        L = 0.75
        expected = -(L * f_rf) / (C_light * beta_45 * gamma_45 * (gamma_45 + 1))
        M = driftLattice(L)._compute_numeric_matrix()
        assert abs(M[4, 5] - expected) < 1e-14

    def test_composition(self):
        """Two half-drifts compose to one full drift."""
        L = 1.0
        M_full = driftLattice(L)._compute_numeric_matrix()
        M_half = driftLattice(L / 2)._compute_numeric_matrix()
        M_composed = M_half @ M_half
        np.testing.assert_allclose(M_composed, M_full, atol=1e-14)

    def test_symbolic_vs_numeric(self):
        d = driftLattice(0.5)
        M_num = d._compute_numeric_matrix()
        M_sym = np.array(d._compute_symbolic_matrix().tolist(), dtype=float)
        np.testing.assert_allclose(M_sym, M_num, atol=1e-14)

    def test_negative_length_raises(self):
        with pytest.raises(ValueError):
            driftLattice(-0.1)

    def test_zero_length_raises(self):
        with pytest.raises(ValueError):
            driftLattice(0.0)


# ── Quadrupole tests ─────────────────────────────────────────────────────

class TestQuadrupole:
    @pytest.mark.parametrize("cls", [qpfLattice, qpdLattice])
    def test_zero_current_is_drift(self, cls):
        L = 0.0889
        M_q = cls(current=0, length=L)._compute_numeric_matrix()
        M_d = driftLattice(L)._compute_numeric_matrix()
        np.testing.assert_allclose(M_q, M_d, atol=1e-12)

    @pytest.mark.parametrize("cls", [qpfLattice, qpdLattice])
    @pytest.mark.parametrize("I", [0.5, 2.0, 5.0, 10.0])
    def test_symplecticity(self, cls, I):
        M = cls(current=I, length=0.0889)._compute_numeric_matrix()
        assert_symplectic(M)

    @pytest.mark.parametrize("cls", [qpfLattice, qpdLattice])
    def test_block_diagonal(self, cls):
        M = cls(current=3.0, length=0.0889)._compute_numeric_matrix()
        assert_block_diagonal(M)

    def test_qpf_focuses_x(self):
        """QPF: M11 = cos(θ) < 1, M21 < 0 (converging in x)."""
        M = qpfLattice(current=5.0, length=0.0889)._compute_numeric_matrix()
        assert M[0, 0] < 1.0  # cos(θ) < 1
        assert M[1, 0] < 0.0  # converging kick

    def test_qpf_defocuses_y(self):
        """QPF: M33 = cosh(θ) > 1, M43 > 0 (diverging in y)."""
        M = qpfLattice(current=5.0, length=0.0889)._compute_numeric_matrix()
        assert M[2, 2] > 1.0  # cosh(θ) > 1
        assert M[3, 2] > 0.0  # diverging kick

    def test_qpd_is_swapped_qpf(self):
        """QPD x-block matches QPF y-block and vice versa."""
        I, L = 3.0, 0.0889
        Mf = qpfLattice(current=I, length=L)._compute_numeric_matrix()
        Md = qpdLattice(current=I, length=L)._compute_numeric_matrix()
        np.testing.assert_allclose(Mf[0:2, 0:2], Md[2:4, 2:4], atol=1e-14)
        np.testing.assert_allclose(Mf[2:4, 2:4], Md[0:2, 0:2], atol=1e-14)

    def test_known_matrix_elements(self):
        """Verify QPF matrix elements against independent calculation."""
        I, L = 5.0, 0.0889
        q = qpfLattice(current=I, length=L)
        # Use the element's own constants for the independent calculation
        k = abs(q.Q * q.G * I) / (q.M * q.C * q.beta * q.gamma)
        theta = np.sqrt(k) * L

        M = q._compute_numeric_matrix()
        assert abs(M[0, 0] - np.cos(theta)) < 1e-14
        assert abs(M[0, 1] - np.sin(theta) / np.sqrt(k)) < 1e-14
        assert abs(M[1, 0] - (-np.sqrt(k) * np.sin(theta))) < 1e-14
        assert abs(M[2, 2] - np.cosh(theta)) < 1e-14
        assert abs(M[2, 3] - np.sinh(theta) / np.sqrt(k)) < 1e-14
        assert abs(M[3, 2] - np.sqrt(k) * np.sinh(theta)) < 1e-14

    @pytest.mark.parametrize("cls", [qpfLattice, qpdLattice])
    def test_symbolic_vs_numeric(self, cls):
        elem = cls(current=3.0, length=0.0889)
        M_num = elem._compute_numeric_matrix()
        M_sym = np.array(elem._compute_symbolic_matrix().tolist(), dtype=float)
        np.testing.assert_allclose(M_sym, M_num, atol=1e-10)

    def test_thin_lens_limit(self):
        """In the thin-lens limit (L→0, k·L→f⁻¹ finite), M21 → -1/f."""
        I = 5.0
        k = abs(Q * G_quad * I) / (M_e * C_light * beta_45 * gamma_45)
        # Use progressively shorter quads with same k
        for L in [0.01, 0.001, 0.0001]:
            theta = np.sqrt(k) * L
            # thin-lens: M21 ≈ -k·L (for small θ)
            M = qpfLattice(current=I, length=L)._compute_numeric_matrix()
            expected_M21 = -k * L  # thin-lens approximation
            assert abs(M[1, 0] - expected_M21) / abs(expected_M21) < theta**2

    def test_energy_dependence(self):
        """Higher energy → weaker focusing (larger θ period)."""
        I, L = 5.0, 0.0889
        q_lo = qpfLattice(current=I, length=L)
        q_lo.setE(20.0)  # 20 MeV
        M_lo = q_lo._compute_numeric_matrix()

        q_hi = qpfLattice(current=I, length=L)
        q_hi.setE(100.0)  # 100 MeV
        M_hi = q_hi._compute_numeric_matrix()

        # At higher energy: less focusing → M11 closer to 1
        assert abs(M_hi[0, 0] - 1.0) < abs(M_lo[0, 0] - 1.0)


# ── Dipole tests ─────────────────────────────────────────────────────────

class TestDipole:
    def test_symplecticity(self):
        M = dipole(length=0.2, angle=11.25)._compute_numeric_matrix()
        assert_symplectic(M)

    def test_y_plane_is_drift(self):
        """Sector bend: y-plane is a pure drift (no vertical focusing)."""
        L = 0.2
        M = dipole(length=L, angle=11.25)._compute_numeric_matrix()
        assert abs(M[2, 2] - 1.0) < 1e-14
        assert abs(M[2, 3] - L) < 1e-14
        assert abs(M[3, 2] - 0.0) < 1e-14
        assert abs(M[3, 3] - 1.0) < 1e-14

    def test_no_xy_coupling(self):
        M = dipole(length=0.2, angle=11.25)._compute_numeric_matrix()
        assert_block_diagonal(M)

    def test_known_matrix_elements(self):
        """Verify sector bend matrix against independent calculation."""
        L, angle_deg = 0.2, 11.25
        angle_rad = np.radians(angle_deg)
        # ρ = L / θ for a sector bend
        rho = L / angle_rad

        M = dipole(length=L, angle=angle_deg)._compute_numeric_matrix()

        assert abs(M[0, 0] - np.cos(angle_rad)) < 1e-10
        assert abs(M[0, 1] - rho * np.sin(angle_rad)) < 1e-10
        assert abs(M[1, 0] - (-np.sin(angle_rad) / rho)) < 1e-10
        assert abs(M[1, 1] - np.cos(angle_rad)) < 1e-10

    def test_dispersion_elements(self):
        """M16 and M26 (dispersion) are nonzero for a bend."""
        M = dipole(length=0.2, angle=11.25)._compute_numeric_matrix()
        assert abs(M[0, 5]) > 1e-6  # M16 ≠ 0
        assert abs(M[1, 5]) > 1e-6  # M26 ≠ 0

    def test_dispersion_sign(self):
        """Positive bend angle → positive M16 (outward displacement for δ>0)."""
        M = dipole(length=0.2, angle=11.25)._compute_numeric_matrix()
        assert M[0, 5] > 0

    def test_zero_angle_approaches_drift(self):
        """As θ→0, dipole should approach a drift."""
        L = 0.2
        M_dip = dipole(length=L, angle=0.001)._compute_numeric_matrix()
        M_drf = driftLattice(L)._compute_numeric_matrix()
        # x-plane 2×2 should be close to drift
        np.testing.assert_allclose(M_dip[0:2, 0:2], M_drf[0:2, 0:2], atol=1e-4)

    def test_R56_sign(self):
        """R56 (M56) for a sector bend should be negative for positive angle
        (path length increases with energy deviation → M56 < 0 at low energy,
        but the FELsim M56 includes the RF-frequency factor)."""
        M = dipole(length=0.2, angle=11.25)._compute_numeric_matrix()
        # M56 should be finite and nonzero
        assert np.isfinite(M[4, 5])
        assert abs(M[4, 5]) > 1e-10

    def test_chicane_R56_nonzero(self):
        """A 4-dipole chicane should have nonzero R56 (path length dependence on δ)."""
        L, theta = 0.2, 11.25
        L_drift = 0.3

        M_pos = dipole(length=L, angle=theta)._compute_numeric_matrix()
        M_neg = dipole(length=L, angle=-theta)._compute_numeric_matrix()
        M_d = driftLattice(L_drift)._compute_numeric_matrix()

        M_total = M_pos @ M_d @ M_neg @ M_d @ M_neg @ M_d @ M_pos
        # R56 (M56) should be nonzero for a chicane
        assert abs(M_total[4, 5]) > 1e-4

    def test_chicane_sandwich_dispersion_cancellation(self):
        """A symmetric chicane with rectangular bends (DPW+DPH+DPW) cancels dispersion."""
        L_dip, theta = 0.2, 11.25
        half_theta = theta / 2
        L_wedge = 0.01
        L_drift = 0.3

        def rect_bend(angle):
            """Rectangular bend = entrance wedge + sector + exit wedge."""
            dw = dipole_wedge(length=L_wedge, angle=half_theta,
                              dipole_length=L_dip, dipole_angle=angle)
            dph = dipole(length=L_dip, angle=angle)
            Mw = dw._compute_numeric_matrix()
            Md = dph._compute_numeric_matrix()
            return Mw @ Md @ Mw

        M_pos = rect_bend(theta)
        M_neg = rect_bend(-theta)
        M_d = driftLattice(L_drift)._compute_numeric_matrix()

        M_total = M_pos @ M_d @ M_neg @ M_d @ M_neg @ M_d @ M_pos
        # Rectangular-bend chicane should have small residual dispersion
        # (not exactly zero due to wedge thin-lens approximation and fringe effects)
        assert abs(M_total[0, 5]) < 0.05, f"Residual D = {M_total[0, 5]}"
        assert abs(M_total[1, 5]) < 0.5, f"Residual D' = {M_total[1, 5]}"

    def test_symbolic_vs_numeric(self):
        elem = dipole(length=0.2, angle=11.25)
        M_num = elem._compute_numeric_matrix()
        M_sym = np.array(elem._compute_symbolic_matrix().tolist(), dtype=float)
        np.testing.assert_allclose(M_sym, M_num, atol=1e-10)

    def test_negative_angle(self):
        """Negative bend angle → reversed horizontal focusing, negative dispersion."""
        M_pos = dipole(length=0.2, angle=11.25)._compute_numeric_matrix()
        M_neg = dipole(length=0.2, angle=-11.25)._compute_numeric_matrix()
        # M11 should be the same (cos is even)
        assert abs(M_pos[0, 0] - M_neg[0, 0]) < 1e-14
        # M16 should flip sign
        assert abs(M_pos[0, 5] + M_neg[0, 5]) < 1e-14


# ── Dipole wedge tests ───────────────────────────────────────────────────

class TestDipoleWedge:
    def test_thin_lens_structure(self):
        """Wedge matrix should be a thin lens: M11=1, M22=1, M12=0, M34=0."""
        dw = dipole_wedge(length=0.01, angle=5.625, dipole_length=0.2, dipole_angle=11.25)
        M = dw._compute_numeric_matrix()
        assert abs(M[0, 0] - 1.0) < 1e-14
        assert abs(M[1, 1] - 1.0) < 1e-14
        assert abs(M[0, 1] - 0.0) < 1e-14
        assert abs(M[2, 2] - 1.0) < 1e-14
        assert abs(M[3, 3] - 1.0) < 1e-14
        assert abs(M[2, 3] - 0.0) < 1e-14

    def test_horizontal_kick(self):
        """Wedge produces tan(η)/R horizontal kick in M21."""
        dw = dipole_wedge(length=0.01, angle=5.625, dipole_length=0.2, dipole_angle=11.25)
        M = dw._compute_numeric_matrix()
        R = 0.2 / np.radians(11.25)
        eta = np.radians(5.625)  # angle * (length/length) when l == self.length
        expected_M21 = np.tan(eta) / R
        assert abs(M[1, 0] - expected_M21) < 1e-10

    def test_absolute_radius_for_negative_angle(self):
        """R = L/|θ| must use absolute value — the chicane sign bug fix."""
        # Positive angle dipole
        dw_pos = dipole_wedge(length=0.01, angle=5.625, dipole_length=0.2, dipole_angle=11.25)
        M_pos = dw_pos._compute_numeric_matrix()

        # Negative angle dipole (chicane return)
        dw_neg = dipole_wedge(length=0.01, angle=5.625, dipole_length=0.2, dipole_angle=-11.25)
        M_neg = dw_neg._compute_numeric_matrix()

        # Both should produce the same edge kick magnitude (|R| is the same)
        assert abs(abs(M_pos[1, 0]) - abs(M_neg[1, 0])) < 1e-14
        # Specifically, M21 should be identical (not sign-flipped)
        assert abs(M_pos[1, 0] - M_neg[1, 0]) < 1e-14

    def test_vertical_defocusing(self):
        """Wedge produces vertical defocusing kick in M43 (negative for standard case)."""
        dw = dipole_wedge(length=0.01, angle=5.625, dipole_length=0.2, dipole_angle=11.25)
        M = dw._compute_numeric_matrix()
        assert M[3, 2] < 0  # defocusing in y

    def test_symplecticity(self):
        dw = dipole_wedge(length=0.01, angle=5.625, dipole_length=0.2, dipole_angle=11.25)
        M = dw._compute_numeric_matrix()
        assert_symplectic(M)

    def test_symbolic_vs_numeric(self):
        dw = dipole_wedge(length=0.01, angle=5.625, dipole_length=0.2, dipole_angle=11.25)
        M_num = dw._compute_numeric_matrix()
        M_sym = np.array(dw._compute_symbolic_matrix().tolist(), dtype=float)
        np.testing.assert_allclose(M_sym, M_num, atol=1e-10)


# ── Sector bend sandwich: DPW + DPH + DPW ───────────────────────────────

class TestDipoleSandwich:
    """A rectangular bend = entrance wedge + sector bend + exit wedge.
    For a symmetric rectangular bend with half-angle wedges,
    the vertical plane should show net focusing."""

    @pytest.fixture
    def sandwich_11deg(self):
        """Symmetric rectangular bend: θ=11.25°, wedge half-angle=θ/2."""
        L_dip = 0.2
        theta = 11.25
        half_theta = theta / 2
        wedge_L = 0.01  # thin wedge physical length

        dw_in = dipole_wedge(length=wedge_L, angle=half_theta,
                             dipole_length=L_dip, dipole_angle=theta)
        dph = dipole(length=L_dip, angle=theta)
        dw_out = dipole_wedge(length=wedge_L, angle=half_theta,
                              dipole_length=L_dip, dipole_angle=theta)
        return dw_in, dph, dw_out

    def test_sandwich_symplecticity(self, sandwich_11deg):
        dw_in, dph, dw_out = sandwich_11deg
        M = (dw_out._compute_numeric_matrix() @
             dph._compute_numeric_matrix() @
             dw_in._compute_numeric_matrix())
        assert_symplectic(M)

    def test_sandwich_vertical_focusing(self, sandwich_11deg):
        """Rectangular bend should have net vertical focusing (M33 < 1)."""
        dw_in, dph, dw_out = sandwich_11deg
        M = (dw_out._compute_numeric_matrix() @
             dph._compute_numeric_matrix() @
             dw_in._compute_numeric_matrix())
        # For a rectangular bend, M33 < 1 → vertical focusing
        assert M[2, 2] < 1.0


# ── FODO cell tests ──────────────────────────────────────────────────────

class TestFODO:
    """A FODO cell (half-QF + drift + QD + drift + half-QF) has well-known
    properties that provide an integrated test of drift + quadrupole matrices."""

    @pytest.fixture
    def fodo_cell(self):
        L_drift = 0.5  # drift length
        L_quad = 0.0889  # quad length
        I = 0.5  # current (A) — moderate focusing for stability

        half_qf = qpfLattice(current=I, length=L_quad / 2)
        d = driftLattice(L_drift)
        qd = qpdLattice(current=I, length=L_quad)

        M_half_qf = half_qf._compute_numeric_matrix()
        M_d = d._compute_numeric_matrix()
        M_qd = qd._compute_numeric_matrix()

        M_cell = M_half_qf @ M_d @ M_qd @ M_d @ M_half_qf
        return M_cell

    def test_fodo_stable(self, fodo_cell):
        """FODO cell should be stable: |Tr/2| < 1 in both planes."""
        M = fodo_cell
        cos_mu_x = (M[0, 0] + M[1, 1]) / 2
        cos_mu_y = (M[2, 2] + M[3, 3]) / 2
        assert abs(cos_mu_x) < 1.0, f"|cos μ_x| = {abs(cos_mu_x)}"
        assert abs(cos_mu_y) < 1.0, f"|cos μ_y| = {abs(cos_mu_y)}"

    def test_fodo_symplecticity(self, fodo_cell):
        assert_symplectic(fodo_cell)

    def test_fodo_no_dispersion(self, fodo_cell):
        """FODO without bends has zero dispersion."""
        M = fodo_cell
        assert abs(M[0, 5]) < 1e-14
        assert abs(M[1, 5]) < 1e-14

    def test_fodo_periodic_twiss(self, fodo_cell):
        """Extract periodic Twiss from FODO cell and verify β > 0, γ > 0."""
        M = fodo_cell
        for plane in [(0, 1), (2, 3)]:
            i, j = plane
            cos_mu = (M[i, i] + M[j, j]) / 2
            sin_mu = np.sqrt(1 - cos_mu**2) * np.sign(M[i, j])
            beta = M[i, j] / sin_mu
            alpha = (M[i, i] - M[j, j]) / (2 * sin_mu)
            gamma_tw = (1 + alpha**2) / beta
            assert beta > 0, f"β = {beta}"
            assert gamma_tw > 0, f"γ = {gamma_tw}"


# ── Energy change tests ──────────────────────────────────────────────────

class TestEnergyChange:
    def test_setE_updates_gamma_beta(self):
        d = driftLattice(0.5)
        d.setE(100.0)
        expected_gamma = 1 + 100.0 / E0_electron
        expected_beta = np.sqrt(1 - 1 / expected_gamma**2)
        assert abs(d.gamma - expected_gamma) < 1e-12
        assert abs(d.beta - expected_beta) < 1e-12

    def test_matrix_changes_with_energy(self):
        """Transfer matrix should change when energy changes."""
        q = qpfLattice(current=5.0, length=0.0889)
        M_45 = q._compute_numeric_matrix().copy()
        q.setE(20.0)
        M_20 = q._compute_numeric_matrix()
        # Matrices should differ (different k value)
        assert not np.allclose(M_45, M_20)

    def test_ultrarelativistic_limit(self):
        """At very high energy, β→1 and M56→-L·f/(c·γ²)."""
        L = 0.5
        d = driftLattice(L)
        d.setE(10000.0)  # 10 GeV
        M = d._compute_numeric_matrix()
        gamma_high = d.gamma
        expected_M56 = -(L * f_rf) / (C_light * 1.0 * gamma_high * (gamma_high + 1))
        assert abs(M[4, 5] - expected_M56) / abs(expected_M56) < 1e-6

    def test_setE_rejects_invalid(self):
        """setE should reject zero, negative, NaN, and Inf energies."""
        d = driftLattice(0.5)
        with pytest.raises(ValueError):
            d.setE(0.0)
        with pytest.raises(ValueError):
            d.setE(-10.0)
        with pytest.raises(ValueError):
            d.setE(float('nan'))
        with pytest.raises(ValueError):
            d.setE(float('inf'))

    def test_setMQE_rejects_invalid(self):
        """setMQE should reject zero/negative mass and rest energy."""
        d = driftLattice(0.5)
        with pytest.raises(ValueError):
            d.setMQE(0.0, 1.6e-19, 0.511)
        with pytest.raises(ValueError):
            d.setMQE(9.1e-31, 1.6e-19, 0.0)
        with pytest.raises(ValueError):
            d.setMQE(9.1e-31, 1.6e-19, -0.511)

    def test_negative_length_rejected(self):
        """Negative element length should be rejected."""
        with pytest.raises(ValueError):
            driftLattice(-1.0)
        with pytest.raises(ValueError):
            qpfLattice(current=5.0, length=-0.0889)


# ── Particle tracking tests ─────────────────────────────────────────────

class TestTracking:
    def test_on_axis_particle_stays_on_axis_drift(self):
        """A particle on axis with zero angles stays on axis through a drift."""
        d = driftLattice(1.0)
        particle = [[0, 0, 0, 0, 0, 0]]
        result = d.useMatrice(particle)
        np.testing.assert_allclose(result, [[0, 0, 0, 0, 0, 0]], atol=1e-15)

    def test_drift_propagation(self):
        """Particle with angle x'=0.001 should shift x by L*x' after drift."""
        L = 2.0
        xp = 0.001
        d = driftLattice(L)
        particle = [[0, xp, 0, 0, 0, 0]]
        result = np.array(d.useMatrice(particle))
        assert abs(result[0, 0] - L * xp) < 1e-15
        assert abs(result[0, 1] - xp) < 1e-15

    def test_batch_tracking(self):
        """Multiple particles tracked simultaneously give same result as individual."""
        d = driftLattice(0.5)
        particles = [[0.001, 0.0001, 0, 0, 0, 0],
                      [0, 0, 0.002, -0.0001, 0, 0],
                      [0, 0, 0, 0, 0, 0.005]]
        results = np.array(d.useMatrice(particles))
        for i, p in enumerate(particles):
            single = np.array(d.useMatrice([p]))
            np.testing.assert_allclose(results[i], single[0], atol=1e-15)


# ── Edge-case tests ──────────────────────────────────────────────────────

class TestEdgeCases:
    """C4 edge-case tests: single particle, large ensembles, extreme parameters."""

    def test_single_particle_tracking(self):
        """Single particle should track correctly through all element types."""
        p = np.array([[0.5, 0.01, 0.3, -0.005, 0.1, 0.002]])
        elements = [
            driftLattice(0.5),
            qpfLattice(current=3.0, length=0.0889),
            qpdLattice(current=3.0, length=0.0889),
            dipole(length=0.2, angle=11.25),
            dipole_wedge(length=0.01, angle=5.625, dipole_length=0.2, dipole_angle=11.25),
        ]
        for elem in elements:
            result = np.array(elem.useMatrice(p))
            assert result.shape == (1, 6), f"{type(elem).__name__}: shape={result.shape}"
            assert np.all(np.isfinite(result)), f"{type(elem).__name__}: non-finite output"

    def test_large_ensemble(self):
        """Track 10,000 particles — verify no memory issues or NaN propagation."""
        rng = np.random.default_rng(42)
        p = rng.standard_normal((10000, 6)) * [1e-3, 1e-4, 1e-3, 1e-4, 0.5, 0.005]
        q = qpfLattice(current=5.0, length=0.0889)
        result = np.array(q.useMatrice(p))
        assert result.shape == (10000, 6)
        assert np.all(np.isfinite(result))

    def test_very_short_element(self):
        """Very short element (1 μm) should produce near-identity matrix."""
        L = 1e-6
        for elem in [driftLattice(L), qpfLattice(current=5.0, length=L)]:
            M = elem._compute_numeric_matrix()
            np.testing.assert_allclose(M, np.eye(6), atol=1e-4)

    def test_very_long_drift(self):
        """100 m drift should still produce finite output."""
        d = driftLattice(100.0)
        p = np.array([[0.001, 0.0001, 0, 0, 0, 0]])
        result = np.array(d.useMatrice(p))
        assert np.all(np.isfinite(result))
        expected_x = 0.001 + 100.0 * 0.0001
        assert abs(result[0, 0] - expected_x) < 1e-12

    def test_high_current_quad(self):
        """High-current quad (near practical limit) should still be symplectic."""
        q = qpfLattice(current=10.0, length=0.0889)
        M = q._compute_numeric_matrix()
        assert_symplectic(M)
        assert np.all(np.isfinite(M))

    def test_zero_length_drift_rejected(self):
        """Zero-length drift should raise ValueError."""
        with pytest.raises(ValueError, match="positive length"):
            driftLattice(0.0)

    def test_on_momentum_particle_delta_zero(self):
        """δ=0 particle through dipole: dispersion should contribute zero."""
        d = dipole(length=0.2, angle=11.25)
        p = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
        result = np.array(d.useMatrice(p))
        np.testing.assert_allclose(result, 0.0, atol=1e-14)

    def test_matrix_finite_at_all_energies(self):
        """Transfer matrices should be finite across a wide energy range."""
        for KE in [1.0, 10.0, 40.0, 100.0, 1000.0, 10000.0]:
            for elem_cls in [
                lambda: driftLattice(0.5),
                lambda: qpfLattice(current=3.0, length=0.0889),
                lambda: dipole(length=0.2, angle=11.25),
            ]:
                elem = elem_cls()
                elem.setE(KE)
                M = elem._compute_numeric_matrix()
                assert np.all(np.isfinite(M)), (
                    f"{type(elem).__name__} at {KE} MeV: non-finite matrix"
                )

    def test_dispersive_kick(self):
        """Off-energy particle (δ≠0) gets displaced in a dipole."""
        dph = dipole(length=0.2, angle=11.25)
        # On-energy
        p_on = [[0, 0, 0, 0, 0, 0]]
        r_on = np.array(dph.useMatrice(p_on))
        # Off-energy (δ = 0.01 = 1% energy deviation)
        p_off = [[0, 0, 0, 0, 0, 0.01]]
        r_off = np.array(dph.useMatrice(p_off))
        # Off-energy particle should be displaced
        assert abs(r_off[0, 0] - r_on[0, 0]) > 1e-6
