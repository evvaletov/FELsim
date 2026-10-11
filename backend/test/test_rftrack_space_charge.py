"""RF-Track adapter space charge (docs/rftrack-bunch6d.md): the known answer for a
round Gaussian bunch in a drift, the bunch-charge input, and the lifetime of
RF-Track's global space-charge engine.

A cold 45 MeV bunch, 1 mm rms round and 0.6 mm rms long, drifts 1 m. Its
rms-equivalent linear space-charge kick gives the slope <x x'>/<x x> =
K L / (4 sigma^2) with the perveance K = 2 I / (I_A (beta gamma)^3), for the
peak current I of the Gaussian line density reduced by 1/sqrt(2), the rms
equivalent of a round Gaussian beam: 5.853e-3 mrad/mm at 1 nC. RF-Track's
PIC_FreeSpace on 32^3 gives 0.987 of it standalone.

Author: Eremey Valetov
"""
import json
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad

pytest.importorskip('RF_Track')

from beamline import driftLattice
from multiCodeSimulator import MultiCodeSimulator, SimSection
from physicalConstants import PhysicalConstants
from rftrackAdapter import RFTrackAdapter
from simulatorBase import BeamlineElement, CoordinateSystem
from simulatorFactory import SimulatorFactory

BACKEND = Path(__file__).resolve().parents[1]
JSON_PATH = BACKEND.parent / 'var' / 'UH_FEL_beamline.json'

ENERGY = 45.0
N = 20000
SIGMA = 1.0         # mm
SIGMA_CT = 0.6      # mm
MESH = (32, 32, 32)
WINDOW = 0.03


def analytic_slope(charge_nc, length=1.0):
    """<x x'>/<x x> [mrad/mm] after a drift of `length` m."""
    me = PhysicalConstants.E0_electron
    gamma = 1.0 + ENERGY / me
    bg = np.sqrt(gamma**2 - 1.0)
    lam = charge_nc * 1e-9 / (np.sqrt(2 * np.pi) * SIGMA_CT * 1e-3) / np.sqrt(2)
    current = lam * bg / gamma * PhysicalConstants.C
    perveance = 2 * current / (17045.0 * bg**3)
    return perveance * length / (4 * (SIGMA * 1e-3)**2)


def slope(p):
    return np.dot(p[:, 0], p[:, 1]) / np.dot(p[:, 0], p[:, 0])


def in_window(value, charge_nc=1.0):
    return abs(value / analytic_slope(charge_nc) - 1.0) <= WINDOW


def drift_line(n=20, length=0.05):
    return [BeamlineElement('DRIFT', length) for _ in range(n)]


def triplet_line():
    """A DPW-DPH-DPW triplet in 1 m of drift, the wedges bent as the DPH:
    RF-Track SBends, and a TransferLine over each wedge's 1 cm."""
    wedge = dict(angle=0.25, pole_gap=0.014478)
    return (drift_line(10)
            + [BeamlineElement('DIPOLE_WEDGE', 0.01, **wedge),
               BeamlineElement('DIPOLE', 0.05, angle=0.5),
               BeamlineElement('DIPOLE_WEDGE', 0.01, **wedge)]
            + drift_line(9))


@pytest.fixture(scope='module')
def bunch():
    """FELsim coordinates of the cold round bunch."""
    rng = np.random.default_rng(20260924)
    p = np.zeros((N, 6))
    p[:, 0] = rng.normal(0.0, SIGMA, N)
    p[:, 2] = rng.normal(0.0, SIGMA, N)
    t_scale = PhysicalConstants.C / PhysicalConstants.f_RF_default
    p[:, 4] = rng.normal(0.0, SIGMA_CT, N) / t_scale
    return p


def sc_adapter(line, charge_nc=1.0, **kw):
    ad = RFTrackAdapter(beam_energy=ENERGY, aperture=None, space_charge=True,
                        bunch_charge_nc=charge_nc, sc_mesh=MESH, **kw)
    ad.set_beamline(line)
    return ad


def last_sample(evolution):
    return evolution.particles[max(evolution.s_positions)]


def test_analytic_slope():
    assert analytic_slope(1.0) == pytest.approx(5.8533e-3, rel=1e-4)


def test_known_answer_simulate(bunch):
    assert in_window(slope(sc_adapter(drift_line()).simulate(bunch).final_particles))


def test_known_answer_collect_evolution(bunch):
    assert in_window(slope(last_sample(sc_adapter(drift_line()).collect_evolution(bunch))))


def test_known_answer_dipole_line(bunch):
    """The kick on a line with a dipole and its wedges: SC on minus SC off
    removes the dipole and edge optics."""
    on = sc_adapter(triplet_line())
    off = RFTrackAdapter(beam_energy=ENERGY, aperture=None)
    off.set_beamline(triplet_line())
    assert in_window(slope(on.simulate(bunch).final_particles)
                     - slope(off.simulate(bunch).final_particles))
    assert in_window(slope(last_sample(on.collect_evolution(bunch)))
                     - slope(last_sample(off.collect_evolution(bunch))))


@pytest.mark.parametrize('column, shift_mm', [(4, 10.0), (4, 50.0), (0, 20.0)])
def test_known_answer_wherever_the_bunch_sits(bunch, column, shift_mm):
    """The reference particle stays at x = 0, t = 0 and is tracked outside
    the space-charge bunch, so a bunch away from it (here 17 and 83 sigma_ct
    late, or 20 sigma_x across) gets the same kick. Inside the bunch it
    would stretch the PIC mesh to reach it."""
    def centred_slope(p):
        x, xp = p[:, 0] - p[:, 0].mean(), p[:, 1] - p[:, 1].mean()
        return np.dot(x, xp) / np.dot(x, x)

    moved = bunch.copy()
    moved[:, column] += shift_mm / (PhysicalConstants.C / PhysicalConstants.f_RF_default
                                    if column == 4 else 1.0)
    here = centred_slope(sc_adapter(drift_line()).simulate(bunch).final_particles)
    there = centred_slope(sc_adapter(drift_line()).simulate(moved).final_particles)
    assert there == pytest.approx(here, rel=1e-3)


def test_track_elements_chains_with_losses(bunch):
    """Chained ranges keep each macro-particle's share of the bunch charge:
    given n_initial, the survivors of a lossy first range carry only their
    own charge into the next, as on one call over both. The second range
    also keeps the whole line's time origin, which the PIC kicks depend on."""
    ad = RFTrackAdapter(beam_energy=ENERGY, aperture=0.0015, space_charge=True,
                        bunch_charge_nc=1.0, sc_mesh=MESH)
    ad.set_beamline(drift_line())
    x0 = ad.transform_coordinates(bunch, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK)
    full = ad.track_elements(x0, 0, 20)
    mid = ad.track_elements(x0, 0, 10)
    assert len(mid) < 0.8 * len(bunch)
    part = ad.track_elements(mid, 10, 20, n_initial=len(bunch))
    assert part.shape == full.shape
    assert np.all(np.abs(part - full).max(axis=0) <= 1e-12 * np.abs(full).max(axis=0))
    with pytest.raises(ValueError, match='n_initial'):
        ad.track_elements(mid, 10, 20, n_initial=len(mid) - 1)


def test_track_elements_without_n_initial_warns_once(bunch, caplog):
    """With space charge on, a range tracked without n_initial takes its
    particles as the whole bunch, and says so once per adapter, unless it
    starts at the first element, where they are the whole bunch."""
    with caplog.at_level('WARNING', logger='rftrackAdapter'):
        off = RFTrackAdapter(beam_energy=ENERGY, aperture=None, bunch_charge_nc=1.0)
        off.set_beamline(drift_line(2))
        x0 = off.transform_coordinates(bunch[:200], CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK)
        off.track_elements(x0, 0, 2)
        ad = sc_adapter(drift_line(2))
        ad.track_elements(x0, 0, 1, n_initial=len(x0))
        ad.track_elements(x0, 0, 1)
        assert 'n_initial' not in caplog.text
        ad.track_elements(x0, 1, 2)
        ad.track_elements(x0, 1, 2)
    assert caplog.text.count('taken as the whole bunch') == 1


def test_linear_in_charge(bunch):
    s = {q: slope(sc_adapter(drift_line(), q).simulate(bunch).final_particles)
         for q in (0.5, 1.0, 2.0)}
    assert s[0.5] / s[1.0] == pytest.approx(0.5, abs=0.005)
    assert s[2.0] / s[1.0] == pytest.approx(2.0, abs=0.02)


def test_zero_charge_bunch_feels_nothing(bunch):
    """With the engine installed and SC steps on every element, a bunch that
    carries no charge is not kicked at all."""
    ad = sc_adapter(drift_line())
    assert all(ad._lattice[i].get_sc_nsteps() > 0 for i in range(ad._lattice.size()))
    rows = ad._bunch_rows(ad.transform_coordinates(bunch, CoordinateSystem.FELSIM,
                                                   CoordinateSystem.RFTRACK))
    rows[:, 8] = 0.0
    ad._install_space_charge()
    out = np.asarray(ad._lattice.track(ad._make_bunch(rows)).get_phase_space('%xp %yp'))
    assert np.all(out == 0.0)


def bends(total_deg, n=20, length=0.05):
    return [BeamlineElement('DIPOLE', length, angle=total_deg / n) for _ in range(n)]


def sc_kick(line, p, sc_nsteps=1):
    """The x and y slopes with space charge on minus off (1 nC)."""
    slopes = []
    for on in (True, False):
        ad = RFTrackAdapter(beam_energy=ENERGY, aperture=None, space_charge=on,
                            bunch_charge_nc=1.0 if on else 0.0, sc_mesh=MESH)
        ad.sc_nsteps = sc_nsteps
        ad.set_beamline(line)
        q = ad.simulate(p).final_particles
        slopes.append([slope(q), slope(q[:, 2:])])
    return np.subtract(*slopes)


def test_known_answer_in_a_weak_bend(bunch):
    """Space charge inside SBends: 1.5 deg over 1 m kicks as the drift does,
    within 0.5 %."""
    bend, drift = sc_kick(bends(1.5), bunch), sc_kick(drift_line(), bunch)
    np.testing.assert_allclose(bend, drift, rtol=5e-3)
    assert in_window(bend[0]) and in_window(bend[1])


def test_one_sbend_takes_its_space_charge_steps(bunch):
    """One 1 m SBend with 20 SC steps equals 20 slices of it (30 deg)."""
    np.testing.assert_allclose(sc_kick([BeamlineElement('DIPOLE', 1.0, angle=30.0)], bunch, 20),
                               sc_kick(bends(30.0), bunch), rtol=1e-3)


def sc_kick_per_offset(line, p, sc_nsteps):
    """<x0 dx'>/<x0 x0> and <y0 dy'>/<y0 y0>, dx' and dy' the exit slopes of
    each particle with space charge (1 nC) minus without."""
    out = []
    for on in (True, False):
        ad = RFTrackAdapter(beam_energy=ENERGY, aperture=None, space_charge=on,
                            bunch_charge_nc=1.0 if on else 0.0, sc_mesh=MESH)
        ad.sc_nsteps = sc_nsteps
        ad.set_beamline(line)
        result = ad.simulate(p)
        q = np.full_like(p, np.nan)
        q[result.metadata['survivor_ids']] = result.final_particles
        out.append(q)
    d = out[0] - out[1]
    return np.array([np.dot(p[:, i], d[:, i + 1]) / np.dot(p[:, i], p[:, i]) for i in (0, 2)])


def test_strong_bend_space_charge_follows_the_beam_size():
    """One 1 m, 30 deg SBend with 20 space-charge steps. Its weak focusing
    shrinks a cold round beam in x as c = cos(h s), so the rms-equivalent
    linear force grows as 2/(c (c + 1)) in x and 2/(c + 1) in y. Carried to
    the exit, the kick per initial offset is, in linear theory, 0.978 of
    the drift's in x and 1.024 in y; RF-Track gives 0.968 and 1.023. Stage
    1 (a drift, then the bend's linear map) gave 1.006 and 1.000. A 6 mm
    long bunch keeps R51's x-z tilt from thinning the beam."""
    th = np.radians(30.0)
    fx = quad(lambda u: 2 * np.cos(th * (1 - u)) / (1 + np.cos(th * u)), 0, 1)[0]
    fy = quad(lambda u: 2 / (1 + np.cos(th * u)), 0, 1)[0]
    rng = np.random.default_rng(20260924)
    p = np.zeros((N, 6))
    p[:, 0] = rng.normal(0.0, SIGMA, N)
    p[:, 2] = rng.normal(0.0, SIGMA, N)
    p[:, 4] = rng.normal(0.0, 6.0, N) / (PhysicalConstants.C / PhysicalConstants.f_RF_default)
    ratio = (sc_kick_per_offset([BeamlineElement('DIPOLE', 1.0, angle=30.0)], p, 20)
             / sc_kick_per_offset([BeamlineElement('DRIFT', 1.0)], p, 20))
    assert ratio[0] == pytest.approx(fx, rel=0.02)
    assert ratio[1] == pytest.approx(fy, rel=5e-3)


def test_space_charge_reaches_every_lattice(bunch, monkeypatch):
    """Every element of every lattice tracked takes SC steps, except the
    1e-12 m SBends that carry the wedges' edges."""
    ad = sc_adapter(triplet_line())
    seen = []
    track = ad._track_rows

    def spy(lat, rows, *args):
        out = track(lat, rows, *args)
        seen.extend((lat[i].get_length(), lat[i].get_sc_nsteps()) for i in range(lat.size()))
        return out

    monkeypatch.setattr(ad, '_track_rows', spy)
    ad.simulate(bunch[:200])
    ad.collect_evolution(bunch[:200])
    ad.track_elements(ad.transform_coordinates(bunch[:200], CoordinateSystem.FELSIM,
                                               CoordinateSystem.RFTRACK), 0, len(ad.beamline))
    n_native = ad._lattice.size()
    assert len(seen) == 3 * n_native and n_native == len(ad.beamline) + 2
    thin = [n for length, n in seen if length == RFTrackAdapter.THIN_WEDGE_LENGTH]
    assert len(thin) == 6 and not any(thin)
    assert all(n > 0 for length, n in seen if length > RFTrackAdapter.THIN_WEDGE_LENGTH)


def steps_tracked(ad, p, monkeypatch):
    """(length, SC steps) of every native element of the lattice simulate tracks."""
    seen = []
    track = ad._track_rows

    def spy(lat, rows, *args):
        out = track(lat, rows, *args)
        seen.extend((lat[i].get_length(), lat[i].get_sc_nsteps()) for i in range(lat.size()))
        return out

    monkeypatch.setattr(ad, '_track_rows', spy)
    ad.simulate(p)
    return seen


def step_line():
    return [BeamlineElement('DRIFT', 0.9), BeamlineElement('DRIFT', 0.04),
            BeamlineElement('QUAD_F', 0.0889, current=1.0)]


def test_space_charge_steps_follow_the_step_length(bunch, monkeypatch):
    """ceil(L / sc_ds_max) steps per element, at least one; none on the
    wedges' thin SBends."""
    ad = sc_adapter(step_line() + triplet_line(), sc_ds_max=0.1)
    seen = steps_tracked(ad, bunch[:200], monkeypatch)
    assert [n for _, n in seen[:3]] == [9, 1, 1]
    assert all(n == (0 if length == RFTrackAdapter.THIN_WEDGE_LENGTH else
                     max(1, int(np.ceil(length / 0.1 - 1e-9))))
               for length, n in seen)


def test_default_step_length_splits_a_long_drift(bunch, monkeypatch):
    ad = sc_adapter(step_line())
    assert ad.sc_ds_max > 0
    seen = steps_tracked(ad, bunch[:200], monkeypatch)
    assert seen[0][1] == int(np.ceil(0.9 / ad.sc_ds_max - 1e-9)) > 1


def test_explicit_steps_per_element_override_the_step_length(bunch, monkeypatch):
    ad = sc_adapter(step_line() + triplet_line(), sc_ds_max=0.1)
    ad.sc_nsteps = 3
    seen = steps_tracked(ad, bunch[:200], monkeypatch)
    assert all(n == (0 if length == RFTrackAdapter.THIN_WEDGE_LENGTH else 3)
               for length, n in seen)


def test_no_step_length_gives_one_step_per_element(bunch, monkeypatch):
    ad = sc_adapter(step_line(), sc_ds_max=None)
    assert [n for _, n in steps_tracked(ad, bunch[:200], monkeypatch)] == [1, 1, 1]


@pytest.mark.parametrize('ds', [0.0, -0.1, np.nan, np.inf, 'x', True, np.True_])
def test_step_length_must_be_positive(ds):
    with pytest.raises(ValueError, match='sc_ds_max'):
        RFTrackAdapter(beam_energy=ENERGY, sc_ds_max=ds)


@pytest.mark.parametrize('n', [0.9, 2.5, 3.0, -1, True, False, np.True_, '3'])
def test_steps_per_element_must_be_a_non_negative_integer(n):
    ad = RFTrackAdapter(beam_energy=ENERGY)
    with pytest.raises(ValueError, match='sc_nsteps'):
        ad.sc_nsteps = n
    assert ad.sc_nsteps is None


@pytest.mark.parametrize('n', [0, 4, np.int64(4), None])
def test_steps_per_element_accepts_integers_and_none(n):
    ad = RFTrackAdapter(beam_energy=ENERGY)
    ad.sc_nsteps = n
    assert ad.sc_nsteps == n
    assert ad._sc_steps(0.9) == (int(np.ceil(0.9 / ad.sc_ds_max - 1e-9)) if n is None else n)


def test_step_length_through_factory_and_multicode(bunch, caplog):
    with caplog.at_level('WARNING'):
        fac = SimulatorFactory.create('rftrack', beam_energy=ENERGY, sc_ds_max=0.07)
        mc = MultiCodeSimulator(sections=[SimSection('s', 'rftrack', (0, 4),
                                                     {'space_charge': True, 'bunch_charge_nc': 1.0,
                                                      'sc_ds_max': 0.07})],
                                lattice_path=str(JSON_PATH), beam_energy=ENERGY)
    assert fac.sc_ds_max == 0.07
    assert [sim.sc_ds_max for sim in mc._simulators.values()] == [0.07]
    assert 'sc_ds_max' not in caplog.text


# ── bunch-charge input ─────────────────────────────────────────────────────

def test_zero_charge_space_charge_raises():
    with pytest.raises(ValueError, match='bunch_charge_nc'):
        RFTrackAdapter(beam_energy=ENERGY, space_charge=True)
    ad = RFTrackAdapter(beam_energy=ENERGY)
    ad.set_beamline(drift_line(2))
    with pytest.raises(ValueError, match='bunch_charge_nc'):
        ad.set_space_charge(True)
    assert not ad.space_charge_enabled
    with pytest.raises(ValueError, match='bunch_charge_nc'):
        ad.set_space_charge(True, bunch_charge_nc=0.0)
    with pytest.raises(ValueError, match='bunch_charge_nc'):
        RFTrackAdapter(beam_energy=ENERGY, bunch_charge_nc=-1.0)
    with pytest.raises(ValueError, match='bunch_charge_nc'):
        SimulatorFactory.create('rftrack', beam_energy=ENERGY, space_charge=True)
    ad = RFTrackAdapter(beam_energy=ENERGY, bunch_charge_nc=1.0, space_charge=True)
    ad.set_beamline(drift_line(2))
    ad.bunch_charge_C = 0.0
    with pytest.raises(ValueError, match='bunch_charge_nc'):
        ad.simulate(np.zeros((3, 6)))


@pytest.mark.parametrize('charge', [-1.0, np.inf, np.nan, None, 'one'])
def test_bunch_charge_must_be_a_finite_magnitude(charge):
    """None too in the constructor; set_space_charge takes None as unchanged."""
    with pytest.raises(ValueError, match='bunch_charge_nc'):
        RFTrackAdapter(beam_energy=ENERGY, bunch_charge_nc=charge)
    ad = RFTrackAdapter(beam_energy=ENERGY, bunch_charge_nc=1.0)
    if charge is not None:
        with pytest.raises(ValueError, match='bunch_charge_nc'):
            ad.set_space_charge(False, bunch_charge_nc=charge)
    assert ad.bunch_charge_C == 1e-9


def test_charge_through_set_space_charge_and_factory(bunch):
    ad = RFTrackAdapter(beam_energy=ENERGY, aperture=None)
    ad.set_beamline(drift_line())
    ad.set_space_charge(True, mesh=MESH, bunch_charge_nc=1.0)
    result = ad.simulate(bunch)
    assert result.metadata['bunch_charge_nc'] == 1.0
    assert in_window(slope(result.final_particles))
    assert ad.get_capabilities()['bunch_charge_nc'] == 1.0

    fac = SimulatorFactory.create('rftrack', beam_energy=ENERGY, aperture=None,
                                  bunch_charge_nc=1.0, space_charge=True, sc_mesh=MESH)
    fac.set_beamline(drift_line())
    assert in_window(slope(fac.simulate(bunch).final_particles))


def test_charge_through_multicode_sections(bunch):
    """The charge is a creation key: sections with different charges get
    their own adapter and track with their own charge, and a section with
    space charge but no charge raises instead of inheriting one."""
    def run(configs):
        mc = MultiCodeSimulator(sections=[SimSection(f's{i}', 'rftrack', rng, cfg)
                                          for i, (rng, cfg) in enumerate(configs)],
                                lattice_path=str(JSON_PATH), beam_energy=ENERGY)
        for sim in mc._simulators.values():
            sim.default_aperture = None
        return mc

    mc = run([((0, 4), {'space_charge': True, 'bunch_charge_nc': 1.0, 'sc_mesh': MESH}),
              ((4, 6), {'space_charge': True, 'bunch_charge_nc': 2.0, 'sc_mesh': MESH})])
    assert len(mc._simulators) == 2
    charges = []
    for sim in mc._simulators.values():
        def spy(particles=None, mode=None, sim=sim, real=sim.simulate):
            charges.append(sim.bunch_charge_C * 1e9)
            return real(particles=particles, mode=mode)
        sim.simulate = spy
    mc.simulate(bunch[:500])
    assert charges == [1.0, 2.0]

    mc = run([((0, 4), {'space_charge': True, 'bunch_charge_nc': 3.0, 'sc_mesh': MESH}),
              ((4, 6), {'space_charge': True})])
    with pytest.raises(ValueError, match='bunch_charge_nc'):
        mc.simulate(bunch[:10])


def test_negative_charge_in_a_multicode_section_raises():
    """A section's charge reaches its adapter's constructor, which rejects a
    negative one whether or not the section turns space charge on."""
    for config in ({'bunch_charge_nc': -2.0}, {'space_charge': True, 'bunch_charge_nc': -2.0}):
        with pytest.raises(ValueError, match='bunch_charge_nc'):
            MultiCodeSimulator(sections=[SimSection('s', 'rftrack', (0, 4), config)],
                               lattice_path=str(JSON_PATH), beam_energy=ENERGY)


def test_multicode_section_gives_the_known_answer(bunch):
    """A drift-only line through MultiCodeSimulator: the section's charge
    reaches the adapter's bunch."""
    ad = sc_adapter(drift_line())
    ref = slope(ad.simulate(bunch).final_particles)
    mc = MultiCodeSimulator(sections=[SimSection('rft', 'rftrack', (0, 20),
                                                 {'space_charge': True, 'bunch_charge_nc': 1.0,
                                                  'sc_mesh': MESH, 'aperture': None})],
                            beam_energy=ENERGY)
    mc._master_beamline = [driftLattice(0.05) for _ in range(20)]
    mc._init_simulators()
    assert slope(mc.simulate(bunch).final_particles) == pytest.approx(ref, rel=1e-12)


def test_multicode_sections_keep_the_charge_of_the_survivors(bunch):
    """A lossy first section and a space-charge second one: the survivors
    keep the charge they carry on one continuous track of the line, as
    n_initial does for track_elements. Spread over the survivors, the bunch
    charge would kick them harder by the fraction lost."""
    ad = RFTrackAdapter(beam_energy=ENERGY, aperture=0.0015, space_charge=True,
                        bunch_charge_nc=1.0, sc_mesh=MESH)
    ad.set_beamline(drift_line())
    ref = ad.simulate(bunch).final_particles
    cfg = {'space_charge': True, 'bunch_charge_nc': 1.0, 'sc_mesh': MESH, 'aperture': 0.0015}
    mc = MultiCodeSimulator(sections=[SimSection('lossy', 'rftrack', (0, 10), cfg),
                                      SimSection('sc', 'rftrack', (10, 20), cfg)],
                            beam_energy=ENERGY)
    mc._master_beamline = [driftLattice(0.05) for _ in range(20)]
    mc._init_simulators()
    result = mc.simulate(bunch)
    assert result.metadata['sections'][1]['num_particles_in'] < 0.8 * len(bunch)
    assert result.final_particles.shape == ref.shape
    assert slope(result.final_particles) == pytest.approx(slope(ref), rel=1e-3)


# ── the global engine (each case in its own process) ────────────────────────

_PRELUDE = textwrap.dedent(f'''
    import gc, json, logging, sys
    import numpy as np
    sys.dont_write_bytecode = True
    sys.path.insert(0, {str(BACKEND)!r})
    logging.disable(logging.WARNING)
    import RF_Track as rft
    from rftrackAdapter import RFTrackAdapter
    from simulatorBase import BeamlineElement

    rng = np.random.default_rng(3)
    beam = np.zeros((2000, 6))
    beam[:, 0] = rng.normal(0, 1, 2000)
    beam[:, 2] = rng.normal(0, 1, 2000)
    beam[:, 4] = rng.normal(0, 5.7, 2000)

    def adapter(space_charge=True, charge=1.0, mesh=16):
        ad = RFTrackAdapter(beam_energy=45.0, aperture=None, space_charge=space_charge,
                            bunch_charge_nc=charge, sc_mesh=(mesh, mesh, mesh))
        ad.set_beamline([BeamlineElement('DRIFT', 0.1) for _ in range(10)])
        return ad

    def run(ad):
        return ad.simulate(beam).final_particles.tolist()

    def junk():
        gc.collect()
        return [rft.SpaceCharge_P2P() for _ in range(50)]

    def result(value):
        print('RESULT ' + json.dumps(value), flush=True)
''')


def _in_process(body):
    proc = subprocess.run([sys.executable, '-c', _PRELUDE + textwrap.dedent(body)],
                          capture_output=True, text=True, timeout=600, cwd=str(BACKEND))
    assert proc.returncode == 0, (proc.returncode, proc.stderr[-2000:])
    line, = [l for l in proc.stdout.splitlines() if l.startswith('RESULT ')]
    return json.loads(line[len('RESULT '):])


def _same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.shape == b.shape
    np.testing.assert_allclose(a, b, rtol=0, atol=1e-12 * np.abs(b).max())


def test_engine_after_space_charge_is_switched_off():
    out = _in_process('''
        ad = adapter()
        run(ad)
        ad.set_space_charge(False)
        j = junk()
        result([run(ad), run(adapter(space_charge=False, charge=0.0))])
    ''')
    _same(out[0], out[1])


def test_engine_after_its_adapter_is_collected():
    alone = _in_process('result(run(adapter(charge=2.0, mesh=32)))')
    out = _in_process('''
        a = adapter()
        run(a)
        del a
        j = junk()
        result(run(adapter(charge=2.0, mesh=32)))
    ''')
    _same(out, alone)


def test_adapters_alternate():
    """Each adapter tracks with its own engine. The mesh shows in the result
    (RF-Track's default engine would give 16^3 and 32^3 the same answer), so
    an engine left installed by the other adapter would too."""
    alone_a = _in_process('result(run(adapter(charge=1.0, mesh=16)))')
    alone_b = _in_process('result(run(adapter(charge=2.0, mesh=32)))')
    out = _in_process('''
        a, b = adapter(charge=1.0, mesh=16), adapter(charge=2.0, mesh=32)
        result([run(a), run(b), run(a), run(b), run(adapter(charge=1.0, mesh=32))])
    ''')
    for got, want in zip(out, (alone_a, alone_b, alone_a, alone_b)):
        _same(got, want)
    assert abs(slope(np.asarray(out[4])) / slope(np.asarray(alone_a)) - 1.0) > 0.01


def test_lattice_outlives_its_adapter():
    """rft.cvar.SC_engine does not own the engine, and a freed one left
    installed crashes the next track with space-charge steps. The class keeps
    the last installed engine, so an adapter's lattice can be tracked after
    the adapter has gone."""
    out = _in_process('''
        a = adapter()
        run(a)
        lat = a.get_lattice()
        del a
        j = junk()
        n = len(beam)
        rows = np.column_stack([beam[:, :5], np.full(n, 45.5), np.full(n, rft.electronmass),
                                np.full(n, -1.0), np.full(n, 1e-9 / 1.602176634e-19 / n),
                                np.arange(n)])
        result(lat.track(rft.Bunch6d(rows)).size())
    ''')
    assert out == 2000


def test_plain_adapter_after_a_collected_space_charge_adapter():
    out = _in_process('''
        a = adapter()
        run(a)
        del a
        j = junk()
        result(len(run(adapter(space_charge=False, charge=0.0))))
    ''')
    assert out == 2000
