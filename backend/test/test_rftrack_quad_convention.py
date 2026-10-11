#!/usr/bin/env python3
"""Pin the RF-Track Quadrupole strength convention (docs/rftrack-bunch6d.md).

RF-Track's integrated strength is S = (P/q) k1 L [MV/c/m] (reference manual
4.2.2): a fixed S gives a particle of momentum P and charge q the focusing
k1 = S q / (P L), so it scales as 1/P, and set_K1(P/q, k1) or
set_K1L(P/q, k1 L) set it for a reference momentum. The adapter calls
set_K1L(Pc/q, k1 L) with FELsim's k1, so each particle sees k1 Pc/P.

These tests first pinned S = k1 L. Their bunch was Bunch6d(m, -1, p0, X),
whose single-particle charge is p0 (the constructor is (mass, population,
charge, X)), so P/q was 1 at every p0; see
docs/rftrack-bunch6d.md. The bunches here are
built per particle, [X XP Y YP T P MASS Q N ID].
"""
import math

import numpy as np
import pytest

rft = pytest.importorskip("RF_Track")

from physicalConstants import PhysicalConstants
from simulatorBase import BeamlineElement, CoordinateSystem

MASS_MEV = 0.51099895
CHARGE_E = -1.0
TOL = 1e-9


def r21(length, p0_mev, strength=None, k1=None, p_ref=None):
    """dxp/dx [1/m] through one quadrupole, set by set_strength(strength) or
    by set_K1(p_ref/q, k1), for an electron of momentum p0_mev."""
    q = rft.Quadrupole()
    q.set_length(length)
    if strength is not None:
        q.set_strength(strength)
    else:
        q.set_K1(p_ref / CHARGE_E, k1)
    lat = rft.Lattice()
    lat.append(q)
    row = [[1.0, 0.0, 0.0, 0.0, 0.0, p0_mev, MASS_MEV, CHARGE_E, 0.0, 0.0]]   # x = 1 mm
    out = np.array(lat.track(rft.Bunch6d(np.array(row))).get_phase_space())
    return out[0, 1] / row[0][0]


def analytic(k1, length):
    """Thick-lens M21. Focusing: -sqrt(k) sin(sqrt(k) L). Defocusing: the
    hyperbolic branch is POSITIVE, +sqrt(|k|) sinh(sqrt(|k|) L)."""
    s = math.sqrt(abs(k1))
    return -s * math.sin(s * length) if k1 > 0 else +s * math.sinh(s * length)


@pytest.mark.parametrize("p0", [10.0, 45.0, 100.0, 400.0])
def test_set_strength_is_momentum_dependent(p0):
    """A fixed S is a field: a particle of momentum p0 sees k1 = S q/(p0 L)."""
    strength = (45.0 / CHARGE_E) * 2.0 * 0.1       # k1 = 2 at 45 MeV/c
    assert r21(0.1, p0, strength=strength) == pytest.approx(analytic(2.0 * 45.0 / p0, 0.1), abs=TOL)


@pytest.mark.parametrize("p0", [10.0, 45.0, 100.0, 400.0])
def test_set_K1_takes_the_rigidity(p0):
    """set_K1(P/q, k1) gives k1 to a particle of that momentum, at any momentum."""
    assert r21(0.1, p0, k1=2.0, p_ref=p0) == pytest.approx(analytic(2.0, 0.1), abs=TOL)


@pytest.mark.parametrize("length", [0.2, 0.1, 0.05, 0.01, 0.001])
def test_set_strength_is_integrated(length):
    """Fixed argument, varying length: the argument is (P/q) k1 L, so k1 = arg q/(P L)."""
    assert r21(length, 45.0, strength=(45.0 / CHARGE_E) * 0.2) == pytest.approx(
        analytic(0.2 / length, length), abs=TOL)


@pytest.mark.parametrize("strength", [0.05, 0.1, 0.2, 0.4, 0.8, -0.2])
def test_set_strength_matches_analytic_thick_lens(strength):
    assert r21(0.1, 45.0, strength=(45.0 / CHARGE_E) * strength) == pytest.approx(
        analytic(strength / 0.1, 0.1), abs=TOL)


def transverse_matrix(adapter, P):
    """4x4 matrix [mm, mrad] of the adapter's line for a particle of momentum
    P [MeV/c], from unit offsets (the elements are matrix-based, so exact)."""
    rows = np.zeros((5, 6))
    rows[:, 5] = P
    rows[1:, :4] = np.eye(4)
    p = adapter.transform_coordinates(rows, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM)
    out = adapter.simulate(p).final_particles
    return (out[1:, :4] - out[0, :4]).T


def thick_quad(k, length):
    from scipy.linalg import expm
    a = np.zeros((4, 4))
    a[0, 1] = a[2, 3] = 1.0
    a[1, 0], a[3, 2] = -k, k
    return expm(a * length)


def test_adapter_convention_matches():
    """What rftrackAdapter builds for a QUAD_F, end to end: FELsim's
    k1 = G I / (B rho) at the reference momentum (F3QC-sized, 2 A)."""
    from rftrackAdapter import RFTrackAdapter
    length, current = 0.0889, 2.0
    ad = RFTrackAdapter(beam_energy=45.0, aperture=None)
    ad.set_beamline([BeamlineElement('QUAD_F', length, current=current)])
    k1 = PhysicalConstants.G_quad_default * current * PhysicalConstants.C * 1e-6 / ad._Pc
    assert k1 == pytest.approx(35.5, abs=0.05)
    np.testing.assert_allclose(transverse_matrix(ad, ad._Pc), thick_quad(k1, length),
                               rtol=0, atol=TOL)


@pytest.mark.parametrize("delta", [-0.02, 0.0, 0.03])
@pytest.mark.parametrize("charge", [-1.0, 1.0])
@pytest.mark.parametrize("kind", ['QUAD_F', 'QUAD_D'])
def test_adapter_quadrupole_is_chromatic_and_normalised(kind, charge, delta):
    """k = k1/(1 + delta) for momentum deviation delta, and QUAD_F focuses in
    x whatever the sign of particle_charge (normalised strengths, as FELsim
    and the Xsuite adapter)."""
    from rftrackAdapter import RFTrackAdapter
    length, current = 0.0889, 2.0
    ad = RFTrackAdapter(beam_energy=45.0, aperture=None, particle_charge=charge)
    ad.set_beamline([BeamlineElement(kind, length, current=current)])
    k1 = PhysicalConstants.G_quad_default * current * PhysicalConstants.C * 1e-6 / ad._Pc
    k1 = k1 if kind == 'QUAD_F' else -k1
    np.testing.assert_allclose(transverse_matrix(ad, ad._Pc * (1 + delta)),
                               thick_quad(k1 / (1 + delta), length), rtol=0, atol=TOL)


def test_set_particle_type_rebuilds_the_strengths():
    """The strengths are set with P/q, so a new charge after set_beamline
    must reach them: QUAD_F still focuses."""
    from rftrackAdapter import RFTrackAdapter
    length, current = 0.0889, 2.0
    ad = RFTrackAdapter(beam_energy=45.0, aperture=None)
    ad.set_beamline([BeamlineElement('QUAD_F', length, current=current)])
    ad.set_particle_type(ad.particle_mass, 1.0)
    k1 = PhysicalConstants.G_quad_default * current * PhysicalConstants.C * 1e-6 / ad._Pc
    np.testing.assert_allclose(transverse_matrix(ad, ad._Pc), thick_quad(k1, length),
                               rtol=0, atol=TOL)
