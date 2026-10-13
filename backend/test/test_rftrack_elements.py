"""RF-Track adapter elements whose strength acts through the particle's charge
or momentum, against their analytic transfer (docs/rftrack-bunch6d.md): the
sextupole's thin-lens kick and the hard-edge solenoid.

Author: Eremey Valetov
"""
import numpy as np
import pytest

pytest.importorskip('RF_Track')

from physicalConstants import PhysicalConstants
from rftrackAdapter import RFTrackAdapter
from simulatorBase import BeamlineElement, CoordinateSystem

C_T = PhysicalConstants.C * 1e-6     # MeV/c per T m


def track(adapter, rows):
    """Final FELsim particles for RF-Track rows [x xp y yp t P]."""
    p = adapter.transform_coordinates(np.asarray(rows, float), CoordinateSystem.RFTRACK,
                                      CoordinateSystem.FELSIM)
    return adapter.simulate(p).final_particles


@pytest.mark.parametrize('charge', [-1.0, 1.0])
def test_sextupole_is_the_normalised_thin_lens(charge):
    """'strength' is k2 L [m^-2]: dx' = -(k2 L / 2)(x^2 - y^2)/(1 + delta) with
    x, y in m and dx' in rad, for either charge; the residual is the O(L)
    thick-element term."""
    k2l, length = 50.0, 1e-4
    ad = RFTrackAdapter(beam_energy=45.0, aperture=None, particle_charge=charge)
    ad.set_beamline([BeamlineElement('SEXTUPOLE', length, strength=k2l)])
    P0 = ad._Pc
    rows = [[1.0, 0, 0.0, 0, 0, P0], [2.0, 0, 1.0, 0, 0, P0], [1.0, 0, 0.0, 0, 0, 1.02 * P0]]
    out = track(ad, rows)
    x, y, P = (np.array([r[i] for r in rows]) for i in (0, 2, 5))
    thin = -0.5 * k2l * (x**2 - y**2) * 1e-3 / (P / P0)
    np.testing.assert_allclose(out[:, 1], thin, rtol=1e-5)


def solenoid_matrix(length, k):
    C, S = np.cos(k * length), np.sin(k * length)
    return np.array([[C * C, S * C / k, S * C, S * S / k],
                     [-k * S * C, C * C, -k * S * S, S * C],
                     [-S * C, -S * S / k, C * C, S * C / k],
                     [k * S * S, -S * C, -k * S * C, C * C]])


@pytest.mark.parametrize('aperture', [None, 0.05])
@pytest.mark.parametrize('momentum', [1.0, 1.03])
@pytest.mark.parametrize('charge', [-1.0, 1.0])
def test_solenoid_is_the_hard_edge_matrix(charge, momentum, aperture):
    """k = q c B0 / (2 P), signed with the charge (the rotation reverses for
    a positive beam); the 4x4 matrix in mm and mrad from central differences."""
    length, field = 0.2, 0.1
    ad = RFTrackAdapter(beam_energy=45.0, aperture=aperture, particle_charge=charge)
    ad.set_beamline([BeamlineElement('SOLENOID', length, field=field)])
    P = ad._Pc * momentum
    eps = 1e-3
    rows = np.zeros((8, 6))
    rows[:, 5] = P
    for i in range(4):
        rows[2 * i, i], rows[2 * i + 1, i] = eps, -eps
    out = track(ad, rows)
    R = np.column_stack([(out[2 * i, :4] - out[2 * i + 1, :4]) / (2 * eps) for i in range(4)])
    k = charge * C_T * field / (2 * P)
    np.testing.assert_allclose(R, solenoid_matrix(length, k), rtol=0, atol=1e-6)


def test_solenoid_without_field_is_a_drift():
    ad = RFTrackAdapter(beam_energy=45.0, aperture=None)
    ad.set_beamline([BeamlineElement('SOLENOID', 0.2)])
    out = track(ad, [[1.0, 2.0, -1.0, 0.5, 0.0, ad._Pc]])
    np.testing.assert_allclose(out[0, :4], [1.4, 2.0, -0.9, 0.5], rtol=1e-12)


@pytest.mark.parametrize('charge', [0, 0.0, np.float64(0.0)])
def test_neutral_particles_are_rejected(charge):
    """Quadrupoles (strength set with P/q) and space charge all
    act through the charge."""
    with pytest.raises(ValueError, match='particle_charge'):
        RFTrackAdapter(beam_energy=45.0, particle_charge=charge)
    ad = RFTrackAdapter(beam_energy=45.0)
    ad.set_beamline([BeamlineElement('QUAD_F', 0.1, current=1.0)])
    with pytest.raises(ValueError, match='particle_charge'):
        ad.set_particle_type(ad.particle_mass, charge)
    assert ad.particle_charge == -1.0
