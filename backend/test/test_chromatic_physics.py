"""Chromatic transfer matrix and aperture physics tests.

Author: Eremey Valetov
"""

import numpy as np
import pytest

from beamline import qpfLattice, qpdLattice, dipole, dipole_wedge, driftLattice, lattice
from ebeam import beam

ENERGY_MEV = 40.0
QUAD_LENGTH = 0.0889
QUAD_CURRENT = 5.0
DIPOLE_LENGTH = 0.3
DIPOLE_ANGLE = 15.0
WEDGE_ANGLE = 7.5
POLE_GAP = 0.014478
N_PARTICLES = 200


@pytest.fixture
def rng():
    return np.random.default_rng(42)


@pytest.fixture
def particles(rng):
    return rng.standard_normal((N_PARTICLES, 6)) * [0.8, 0.1, 0.8, 0.1, 0.5, 5.0]


@pytest.fixture
def mono_particles(rng):
    """On-momentum particles (delta=0)."""
    p = rng.standard_normal((N_PARTICLES, 6)) * [0.8, 0.1, 0.8, 0.1, 0.5, 0.0]
    return p


def _make_quad(cls, current=QUAD_CURRENT, chromatic=False):
    q = cls(current=current, length=QUAD_LENGTH)
    q.setE(ENERGY_MEV)
    q.chromatic = chromatic
    return q


def _make_dipole(angle=DIPOLE_ANGLE, length=DIPOLE_LENGTH, chromatic=False):
    d = dipole(length=length, angle=angle)
    d.setE(ENERGY_MEV)
    d.chromatic = chromatic
    return d


def _make_wedge(chromatic=False):
    w = dipole_wedge(length=0.01, angle=WEDGE_ANGLE,
                     dipole_length=DIPOLE_LENGTH, dipole_angle=DIPOLE_ANGLE,
                     pole_gap=POLE_GAP)
    w.setE(ENERGY_MEV)
    w.chromatic = chromatic
    return w


def _extract_4x4(element, delta):
    """Numerically extract the transverse 4x4 transfer matrix at a given delta."""
    M = np.zeros((4, 4))
    r0 = np.zeros((1, 6))
    r0[0, 5] = delta
    base = np.array(element.useMatrice(r0))
    for i in range(4):
        e = np.zeros((1, 6))
        e[0, i] = 1e-6
        e[0, 5] = delta
        out = np.array(element.useMatrice(e))
        M[:, i] = (out[0, :4] - base[0, :4]) / 1e-6
    return M


# ---------------------------------------------------------------------------
# Quadrupole tests
# ---------------------------------------------------------------------------

class TestChromaticQuad:

    @pytest.mark.parametrize("cls", [qpfLattice, qpdLattice])
    def test_on_momentum_consistency(self, cls, mono_particles):
        q_std = _make_quad(cls, chromatic=False)
        q_chr = _make_quad(cls, chromatic=True)
        r_std = np.array(q_std.useMatrice(mono_particles))
        r_chr = np.array(q_chr.useMatrice(mono_particles))
        np.testing.assert_allclose(r_chr, r_std, atol=1e-12)

    @pytest.mark.parametrize("cls", [qpfLattice, qpdLattice])
    def test_chromaticity_effect(self, cls, particles):
        q_std = _make_quad(cls, chromatic=False)
        q_chr = _make_quad(cls, chromatic=True)
        r_std = np.array(q_std.useMatrice(particles))
        r_chr = np.array(q_chr.useMatrice(particles))
        assert np.max(np.abs(r_chr - r_std)) > 1e-6

    @pytest.mark.parametrize("cls", [qpfLattice, qpdLattice])
    @pytest.mark.parametrize("delta", [0.0, 5.0, -5.0, 10.0])
    def test_symplecticity(self, cls, delta):
        q = _make_quad(cls, chromatic=True)
        M = _extract_4x4(q, delta)
        det_x = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
        det_y = M[2, 2] * M[3, 3] - M[2, 3] * M[3, 2]
        np.testing.assert_allclose(det_x, 1.0, atol=1e-8)
        np.testing.assert_allclose(det_y, 1.0, atol=1e-8)

    @pytest.mark.parametrize("cls", [qpfLattice, qpdLattice])
    def test_zero_current_fallback(self, cls):
        p = np.array([[1.0, 0.1, 2.0, 0.2, 0.5, 5.0]])
        q_std = _make_quad(cls, current=0.0, chromatic=False)
        q_chr = _make_quad(cls, current=0.0, chromatic=True)
        np.testing.assert_allclose(q_chr.useMatrice(p), q_std.useMatrice(p), atol=1e-14)

    def test_rftrack_formula(self):
        """k_felsim = k0*bg0/bg must equal k_rft = k0*P0/P."""
        q = _make_quad(qpfLattice)
        E0 = q.E0
        K0 = q.E
        P0 = q.beta * q.gamma * E0  # MeV/c

        delta_pct = 0.5
        delta_coord6 = delta_pct * 10
        K_p = K0 * (1 + delta_coord6 * 1e-3)
        E_total_p = K_p + E0
        P_p = np.sqrt(E_total_p**2 - E0**2)

        k0 = np.abs(q.Q * q.G * q.current / (q.M * q.C * q.beta * q.gamma))
        k_rft = k0 * P0 / P_p

        gamma_p = E_total_p / E0
        bg_p = np.sqrt(gamma_p**2 - 1)
        k_felsim = k0 * (q.beta * q.gamma) / bg_p

        np.testing.assert_allclose(k_felsim, k_rft, rtol=1e-12)


# ---------------------------------------------------------------------------
# Dipole sector-bend tests
# ---------------------------------------------------------------------------

class TestChromaticDipole:

    def test_on_momentum_consistency(self, mono_particles):
        d_std = _make_dipole(chromatic=False)
        d_chr = _make_dipole(chromatic=True)
        np.testing.assert_allclose(
            d_chr.useMatrice(mono_particles),
            d_std.useMatrice(mono_particles),
            atol=1e-10,
        )

    def test_chromaticity_effect(self, particles):
        d_std = _make_dipole(chromatic=False)
        d_chr = _make_dipole(chromatic=True)
        diff = np.max(np.abs(
            np.array(d_chr.useMatrice(particles)) -
            np.array(d_std.useMatrice(particles))
        ))
        assert diff > 1e-6

    @pytest.mark.parametrize("delta", [0.0, 5.0, -5.0, 10.0])
    def test_symplecticity(self, delta):
        d = _make_dipole(chromatic=True)
        M = _extract_4x4(d, delta)
        det_x = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
        det_y = M[2, 2] * M[3, 3] - M[2, 3] * M[3, 2]
        np.testing.assert_allclose(det_x, 1.0, atol=1e-8)
        np.testing.assert_allclose(det_y, 1.0, atol=1e-8)

    def test_zero_angle_fallback(self):
        p = np.array([[1.0, 0.1, 2.0, 0.2, 0.5, 5.0]])
        d_std = _make_dipole(angle=0.0, chromatic=False)
        d_chr = _make_dipole(angle=0.0, chromatic=True)
        np.testing.assert_allclose(
            d_chr.useMatrice(p), d_std.useMatrice(p), atol=1e-14
        )


# ---------------------------------------------------------------------------
# Dipole wedge tests
# ---------------------------------------------------------------------------

class TestChromaticWedge:

    def test_on_momentum_consistency(self, mono_particles):
        w_std = _make_wedge(chromatic=False)
        w_chr = _make_wedge(chromatic=True)
        np.testing.assert_allclose(
            w_chr.useMatrice(mono_particles),
            w_std.useMatrice(mono_particles),
            atol=1e-10,
        )

    def test_chromaticity_effect(self, particles):
        w_std = _make_wedge(chromatic=False)
        w_chr = _make_wedge(chromatic=True)
        diff = np.max(np.abs(
            np.array(w_chr.useMatrice(particles)) -
            np.array(w_std.useMatrice(particles))
        ))
        assert diff > 1e-6


# ---------------------------------------------------------------------------
# Aperture tests
# ---------------------------------------------------------------------------

class TestAperture:

    def test_clips_particles(self):
        elem = driftLattice(length=0.1)
        elem.aperture_x = 5.0  # mm
        elem.aperture_y = 5.0
        p = np.array([
            [1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [10.0, 0.0, 1.0, 0.0, 0.0, 0.0],  # outside x
            [1.0, 0.0, 10.0, 0.0, 0.0, 0.0],  # outside y
        ])
        surviving = elem.apply_aperture(p)
        assert surviving.shape[0] == 1
        np.testing.assert_allclose(surviving[0], p[0])

    def test_no_clip_within_bounds(self):
        elem = driftLattice(length=0.1)
        elem.aperture_x = 5.0
        elem.aperture_y = 5.0
        p = np.array([
            [1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [-4.9, 0.0, 4.9, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ])
        surviving = elem.apply_aperture(p)
        assert surviving.shape[0] == 3

    def test_none_passes_all(self):
        elem = driftLattice(length=0.1)
        assert elem.aperture_x is None and elem.aperture_y is None
        p = np.array([[100.0, 0.0, 200.0, 0.0, 0.0, 0.0]])
        surviving = elem.apply_aperture(p)
        assert surviving.shape[0] == 1

    @pytest.mark.parametrize("cls", [qpfLattice, qpdLattice])
    def test_quad_default_aperture(self, cls):
        q = cls(current=1.0)
        assert q.aperture_x == 13.5
        assert q.aperture_y == 13.5


# ---------------------------------------------------------------------------
# Integration: FODO cell with chromatic on/off
# ---------------------------------------------------------------------------

class TestChromaticBeamlineTracking:

    def test_fodo_chromatic_differs(self, particles):
        """Track through a short FODO cell; chromatic and standard should differ."""
        qf = qpfLattice(current=5.0, length=QUAD_LENGTH)
        d1 = driftLattice(length=0.2)
        qd = qpdLattice(current=5.0, length=QUAD_LENGTH)
        d2 = driftLattice(length=0.2)
        elements = [qf, d1, qd, d2]
        for el in elements:
            el.setE(ENERGY_MEV)

        def _track(p, chromatic):
            out = p.copy()
            for el in elements:
                el.chromatic = chromatic
                out = np.array(el.useMatrice(out))
            return out

        r_std = _track(particles, chromatic=False)
        r_chr = _track(particles, chromatic=True)
        assert np.max(np.abs(r_chr - r_std)) > 1e-6
