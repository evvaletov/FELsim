"""The RF-Track adapter's bunches: one constructor, the
charge, mass, weight and id of every macro-particle through element-by-element
tracking, the reference particle behind the longitudinal output, and the
lost-particle table.

Author: Eremey Valetov
"""
import ast
import copy
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip('RF_Track')

import rftrackAdapter
from physicalConstants import PhysicalConstants
from rftrackAdapter import RFTrackAdapter
from simulatorBase import BeamlineElement, CoordinateSystem

JSON_PATH = Path(__file__).resolve().parents[2] / 'var' / 'UH_FEL_beamline.json'
LINAC_PATH = JSON_PATH.with_name('slac_linac.json')
E = PhysicalConstants.Q


def beam(n, sigmas, seed=5):
    rng = np.random.default_rng(seed)
    return np.column_stack([rng.normal(0.0, s, n) for s in sigmas])


def lossy_dipole_line():
    """A DPW-DPH-DPW triplet between two drifts; the wedges take their bend
    from the DPH."""
    wedge = dict(angle=5.0, pole_gap=0.014478)
    return [BeamlineElement('DRIFT', 0.5),
            BeamlineElement('DIPOLE_WEDGE', 0.01, **wedge),
            BeamlineElement('DIPOLE', 0.2, angle=10.0),
            BeamlineElement('DIPOLE_WEDGE', 0.01, **wedge),
            BeamlineElement('DRIFT', 0.5)]


def lossy_adapter(**kw):
    ad = RFTrackAdapter(beam_energy=45.0, aperture=0.01, **kw)
    ad.set_beamline(lossy_dipole_line())
    return ad


LOSSY_BEAM = (6.0, 6.0, 6.0, 6.0, 5.0, 5.0)


# ── one constructor ─────────────────────────────────────────────────────────

def test_one_bunch_constructor():
    """Every RF-Track bunch the adapter builds goes through _make_bunch."""
    tree = ast.parse(Path(rftrackAdapter.__file__).read_text())
    owners = []
    for fn in ast.walk(tree):
        if isinstance(fn, ast.FunctionDef):
            for node in ast.walk(fn):
                if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                        and node.func.attr in ('Bunch6d', 'Bunch6dT')):
                    owners.append(fn.name)
    assert owners == ['_make_bunch']


@pytest.mark.parametrize('charge', [-1.0, 1.0])
def test_bunch_carries_charge_mass_and_weight(charge):
    ad = RFTrackAdapter(beam_energy=45.0, aperture=None, particle_charge=charge,
                        bunch_charge_nc=0.3)
    ad.set_beamline([BeamlineElement('DRIFT', 0.5)])
    ad.simulate(beam(1000, (1, 1, 1, 1, 1, 1)))
    b = ad.get_bunch()
    assert b is ad._bunch and b.size() == 1000
    assert np.all(b.get_phase_space('%Q') == charge)
    assert np.all(b.get_phase_space('%m') == ad.particle_mass)
    assert b.get_phase_space('%N').sum() * E == pytest.approx(0.3e-9, rel=1e-12)
    np.testing.assert_array_equal(b.get_phase_space('%id').ravel(), np.arange(1000))


def test_survivors_keep_their_id_and_weight_through_collect_evolution(monkeypatch):
    """collect_evolution tracks the line element by element, one bunch
    rebuilt per element: ids and weights survive losses in several of them,
    and the survivors are simulate's."""
    ad = lossy_adapter(bunch_charge_nc=0.5)
    X = beam(3000, LOSSY_BEAM)
    tracked = []
    track = ad._track_rows

    def spy(lat, rows, *args):
        out = track(lat, rows, *args)
        tracked.append(out)
        return out

    monkeypatch.setattr(ad, '_track_rows', spy)
    ad.collect_evolution(X, checkpoint_elements=[])
    monkeypatch.undo()
    result = ad.simulate(X)
    n_good = result.metadata['num_good']
    assert 0 < n_good < 3000 and len(tracked) == 5
    assert len({t.shape[0] for t in tracked}) > 2           # losses in several elements
    final = tracked[-1]
    assert final[0, 9] == RFTrackAdapter.REFERENCE_ID
    assert final[1:, 8].sum() * E == pytest.approx(0.5e-9 * n_good / 3000, rel=1e-12)
    ids = result.metadata['survivor_ids']
    assert sorted(ids) == sorted(final[1:, 9])
    alone = ad.simulate(X[ids[:20]]).final_particles        # no space charge: independent
    np.testing.assert_array_equal(alone, result.final_particles[:20])


def test_collect_evolution_survives_total_loss():
    ad = RFTrackAdapter(lattice_path=str(JSON_PATH), beam_energy=45.0)
    ad.enable_physical_apertures()
    X = beam(200, (5.0, 1.0, 5.0, 1.0, 5.7, 5.0))
    assert ad.simulate(X).metadata['num_good'] == 0
    evolution = ad.collect_evolution(X)
    n_good = np.array([e.parameters['n_good'] for e in evolution.elements])
    assert len(n_good) == len(ad.beamline) and n_good[-1] == 0
    assert np.all(np.diff(n_good) <= 0)


# ── the reference particle and coord5 ───────────────────────────────────────

def drift_quad_line():
    return [BeamlineElement('DRIFT', 0.7), BeamlineElement('QUAD_F', 0.0889, current=2.0),
            BeamlineElement('DRIFT', 0.3), BeamlineElement('QUAD_D', 0.0889, current=1.5),
            BeamlineElement('DRIFT', 1.1)]


def test_reference_stays_at_zero_and_coord5_is_kept():
    """coord5 is the arrival time relative to the reference particle: the
    reference stays at 0, and an on-axis on-momentum particle launched at
    coord5 = c arrives at c."""
    ad = RFTrackAdapter(beam_energy=45.0, aperture=None)
    ad.set_beamline(drift_quad_line())
    p = np.zeros((3, 6))
    p[1, 4], p[2, 4] = 3.0, -7.5
    out = ad.simulate(p).final_particles
    np.testing.assert_allclose(out[:, 4], p[:, 4], rtol=0, atol=1e-9)
    assert out[0, 4] == 0.0
    evolution = ad.collect_evolution(p)
    for s in evolution.s_positions:
        np.testing.assert_allclose(evolution.particles[s][:, 4], p[:, 4], rtol=0, atol=1e-9)


def test_off_momentum_particle_drifts_relative_to_the_reference():
    ad = RFTrackAdapter(beam_energy=45.0, aperture=None)
    ad.set_beamline([BeamlineElement('DRIFT', 2.0)])
    p = np.zeros((2, 6))
    p[1, 5] = -20.0                                      # dK/K = -2 %
    P = ad.transform_coordinates(p, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK)[:, 5]
    beta = P / np.hypot(P, ad.particle_mass)
    expected = 2.0e3 * (1 / beta[1] - 1 / beta[0]) / ad._t_scale
    out = ad.simulate(p).final_particles
    assert out[1, 4] == pytest.approx(expected, rel=1e-9)


def test_track_elements_chains():
    """track_elements(a, b) then (b, c) equals track_elements(a, c), cuts
    between the wedges and the dipole of a triplet (99, 100, 101) included. t
    is relative to the reference, so its rounding is that of the arrival
    time: it is compared on the scale of the reference's time of flight."""
    ad = RFTrackAdapter(lattice_path=str(JSON_PATH), beam_energy=45.0)
    X = beam(2000, (0.3, 0.03, 0.3, 0.03, 5.7, 5.0))
    x0 = ad.transform_coordinates(X, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK)
    n = len(ad.beamline)
    full = ad.track_elements(x0, 0, n)
    tof = sum(e.length for e in ad.beamline) * 1e3 / ad._beta       # mm/c
    for b in (23, 46, 99, 100, 101):
        part = ad.track_elements(ad.track_elements(x0, 0, b), b, n)
        assert part.shape == full.shape and 0 < full.shape[0] < 2000
        scale = np.abs(full).max(axis=0)
        scale[4] = tof
        assert np.all(np.abs(part - full).max(axis=0) <= 1e-12 * scale)


def assert_chains(ad, x0, cuts, rtol):
    """track_elements over consecutive ranges equals one call over their
    union: same survivors, coordinates within rtol of each column's scale."""
    full = ad.track_elements(x0, cuts[0], cuts[-1])
    part = x0
    for a, b in zip(cuts[:-1], cuts[1:]):
        part = ad.track_elements(part, a, b)
    assert part.shape == full.shape
    assert np.all(np.abs(part - full).max(axis=0) <= rtol * np.abs(full).max(axis=0))


def two_cavity_line():
    rfc = RFTrackAdapter(lattice_path=str(LINAC_PATH), beam_energy=1.0).beamline[0]
    return [BeamlineElement('DRIFT', 0.3), copy.deepcopy(rfc), BeamlineElement('DRIFT', 1.0),
            copy.deepcopy(rfc), BeamlineElement('DRIFT', 0.5)]


@pytest.mark.parametrize('cut', [1, 2, 3, 4])
def test_track_elements_chains_through_rf_cavities(cut):
    """A range that starts downstream of a cavity phases its own cavities on
    the accelerated reference and measures t from it."""
    ad = RFTrackAdapter(beam_energy=1.0, aperture=None)
    ad.set_beamline(two_cavity_line())
    X = beam(50, (0.5, 0.05, 0.5, 0.05, 2.0, 5.0))
    x0 = ad.transform_coordinates(X, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK)
    assert_chains(ad, x0, [0, cut, 5], 1e-12)


# ── lost particles ──────────────────────────────────────────────────────────

def test_lost_particles_on_a_dipole_line():
    """One row per lost particle, ids distinct and disjoint from the
    survivors', and S (the start in the adapter's line of the element that
    lost it: a drift, a wedge or the dipole) sorts the losses into the
    elements as collect_evolution counts them."""
    ad = lossy_adapter()
    X = beam(3000, LOSSY_BEAM)
    result = ad.simulate(X)
    lost = result.metadata['lost_particles']
    assert lost.shape == (result.metadata['num_lost'], 11)
    ids = lost[:, 10].astype(int)
    survivors = result.metadata['survivor_ids']
    assert len(set(ids)) == len(ids) and set(ids).isdisjoint(survivors)
    assert np.all((ids >= 0) & (ids < 3000)) and len(ids) + len(survivors) == 3000

    evolution = ad.collect_evolution(X)
    starts = np.array([e.s_start for e in evolution.elements])
    per_element = np.array([np.sum(np.isclose(lost[:, 6], s, rtol=0, atol=1e-9)) for s in starts])
    assert per_element.sum() == len(ids)
    np.testing.assert_array_equal(per_element, [e.parameters['n_lost'] for e in evolution.elements])
    assert np.count_nonzero(per_element) >= 2


def quad_solenoid_line(with_dipole):
    line = [BeamlineElement('DRIFT', 0.3), BeamlineElement('QUAD_D', 0.0889, current=10.0),
            BeamlineElement('DRIFT', 0.2), BeamlineElement('QUAD_F', 0.0889, current=10.0),
            BeamlineElement('QUAD_D', 0.0889, current=10.0),
            BeamlineElement('SOLENOID', 0.2, field=0.5), BeamlineElement('DRIFT', 0.4)]
    if with_dipole:
        line[3:3] = lossy_dipole_line()[1:4]
    return line


@pytest.mark.parametrize('with_dipole', [False, True])
def test_lost_particles_in_quadrupoles_and_solenoids(with_dipole):
    """A quadrupole or solenoid reports a particle that leaves its aperture
    inside it at its exit, where the next element starts; the table books it
    to the element that lost it, as collect_evolution does."""
    ad = RFTrackAdapter(beam_energy=45.0, aperture=0.01)
    ad.set_beamline(quad_solenoid_line(with_dipole))
    X = beam(3000, (3.0, 1.0, 3.0, 1.0, 5.0, 5.0))
    result = ad.simulate(X)
    lost = result.metadata['lost_particles']
    assert lost.shape == (result.metadata['num_lost'], 11)
    evolution = ad.collect_evolution(X)
    starts = [e.s_start for e in evolution.elements]
    per_element = [np.sum(np.isclose(lost[:, 6], s, rtol=0, atol=1e-9)) for s in starts]
    n_lost = [e.parameters['n_lost'] for e in evolution.elements]
    assert per_element == n_lost
    thick = [i for i, e in enumerate(ad.beamline) if e.element_type in ('QUAD_F', 'QUAD_D', 'SOLENOID')]
    assert all(n_lost[i] > 0 and n_lost[i + 1] > 0 for i in thick)


def test_no_loss_no_table():
    ad = RFTrackAdapter(beam_energy=45.0, aperture=None)
    ad.set_beamline(lossy_dipole_line())
    metadata = ad.simulate(beam(100, (1, 1, 1, 1, 1, 1))).metadata
    assert metadata['lost_particles'] is None and metadata['lost_elements'] is None
