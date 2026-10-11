"""The RF-Track adapter's dipoles (docs/rftrack-bunch6d.md): a DPH is RF-Track's
SBend, a DPW a thin SBend that carries its edge followed by a TransferLine
over the wedge's length, all in one Lattice.

Linear maps are central differences in FELsim units (mm, mrad, 1e-3 ToF/T_RF,
1e-3 dK/K) with steps 1e-3 in x, x', y, y' and 1e-2 in coord5 and coord6; a
smaller step loses R56 to the round-off of the absolute arrival time.

Author: Eremey Valetov
"""
import ast
import math
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip('RF_Track')

import RF_Track
import latticeLoader
import rftrackAdapter
from beamline import dipole, dipole_wedge
from multiCodeSimulator import MultiCodeSimulator, SimSection
from rftrackAdapter import RFTrackAdapter
from simulatorBase import BeamlineElement, CoordinateSystem

REPO = Path(__file__).resolve().parents[2]
JSON_PATH = REPO / 'var' / 'UH_FEL_beamline.json'
XLSX3_PATH = REPO / 'beam_excel' / 'Beamline_elements_3.xlsx'
XLSX_PATH = REPO / 'beam_excel' / 'Beamline_elements.xlsx'
FEL, RFT = CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK
STEP = np.array([1e-3, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2])

DPH_CASES = [(0.0889, 1.5), (0.04064, -4.0), (0.04064, 5.0), (0.037389, 11.25),
             (0.037389, -11.25), (0.2, 22.5)]
SPECIES = {'electron': (None, -1.0), 'positron': (None, 1.0),
           'proton': (RF_Track.protonmass, 1.0)}
# (dipole_length, dipole_angle) of the UH wedges: R = 3.39573, 0.582125,
# 0.4657, 0.190421 and 0.203718 m
WEDGE_BENDS = [(0.0889, 1.5), (0.04064, -4.0), (0.04064, 5.0), (0.037389, 11.25), (0.04, -11.25)]
ETAS = [0.0, 0.75, 1.5, 2.018, 2.536, 11.25, -11.25]


def fd_map(track, step=STEP):
    """6x6 map of track: (n, 6) -> (n, 6), FELsim coordinates in and out."""
    X = np.zeros((12, 6))
    X[0::2][np.diag_indices(6)] = step
    X[1::2][np.diag_indices(6)] = -step
    Y = track(X)
    return (Y[0::2] - Y[1::2]).T / (2 * step)


def line_map(ad, a=0, b=None, step=STEP):
    """Map of the adapter's beamline[a:b], tracked by track_elements."""
    b = len(ad.beamline) if b is None else b

    def track(X):
        out = ad.track_elements(ad.transform_coordinates(X, FEL, RFT), a, b)
        return ad.transform_coordinates(out, RFT, FEL)
    return fd_map(track, step)


def felsim_map(elements):
    M = np.eye(6)
    for e in elements:
        M = e._compute_numeric_matrix() @ M
    return M


def on_beam(elem, E, species='electron'):
    if species == 'proton':
        elem.changeBeamType('proton', E)
    else:
        elem.setE(E)
    return elem


def adapter(E=45.0, species='electron', **kw):
    mass, charge = SPECIES[species]
    kw.setdefault('aperture', None)
    return RFTrackAdapter(beam_energy=E, particle_mass=mass, particle_charge=charge, **kw)


def reference_exit(ad):
    return ad._track_rows(ad._lattice, ad._bunch_rows(np.empty((0, 6))))[0]


# ── D1: DPH ─────────────────────────────────────────────────────────────────

@pytest.mark.parametrize('species, E', [('electron', 45.0), ('electron', 40.0),
                                        ('positron', 45.0), ('positron', 40.0),
                                        ('proton', 45.0)])
@pytest.mark.parametrize('L, deg', DPH_CASES)
def test_dph_is_felsims_dipole(species, E, L, deg):
    """Every entry of the 6x6 map within 1e-8 max(1, |R_ij|) of FELsim's
    dipole, R56 within 1e-6; the reference particle stays on axis."""
    ad = adapter(E, species)
    ad.set_beamline([BeamlineElement('DIPOLE', L, angle=deg)])
    M = line_map(ad)
    Mf = on_beam(dipole(length=L, angle=deg), E, species)._compute_numeric_matrix()
    bound = 1e-8 * np.maximum(1.0, np.abs(Mf))
    bound[4, 5] = np.inf
    assert np.all(np.abs(M - Mf) <= bound)
    assert M[4, 5] == pytest.approx(Mf[4, 5], rel=1e-6)
    assert np.abs(reference_exit(ad)[:4]).max() < 1e-12


# ── D2: DPW ─────────────────────────────────────────────────────────────────

@pytest.mark.parametrize('eta', ETAS)
def test_dpw_is_felsims_wedge(eta):
    """M21 = tan(eta)/R and M43 = -tan(eta - phi)/R within 1e-9 (relative),
    M56 within 1e-6, the rest the identity within 1e-10, for the UH wedge
    bends, both charge signs and a wedge without a pole gap."""
    for E in (40.0, 45.0):
        for species in ('electron', 'positron'):
            for (dl, da), gap in zip(WEDGE_BENDS, (0.014478, 0.014478, 0.0, 0.0127, 0.0127)):
                ad = adapter(E, species)
                ad.set_beamline([BeamlineElement('DIPOLE_WEDGE', 0.01, angle=eta, dipole_length=dl,
                                                 dipole_angle=da, pole_gap=gap)])
                M = line_map(ad)
                Mf = on_beam(dipole_wedge(0.01, angle=eta, dipole_length=dl, dipole_angle=da,
                                          pole_gap=gap), E)._compute_numeric_matrix()
                for i, j in ((1, 0), (3, 2), (4, 5)):
                    tol = (1e-6 if (i, j) == (4, 5) else 1e-9) * abs(Mf[i, j]) or 1e-10
                    assert abs(M[i, j] - Mf[i, j]) <= tol, (dl, da, gap, i, j)
                rest = np.abs(M - Mf)
                rest[1, 0] = rest[3, 2] = rest[4, 5] = 0.0
                assert rest.max() <= 1e-10


def test_dpw_takes_its_own_bend():
    """A DPW's R is its own dipole_length/dipole_angle, as in FELsim, even
    beside a DPH of another curvature."""
    ad = adapter()
    ad.set_beamline([BeamlineElement('DIPOLE_WEDGE', 0.01, angle=2.018, dipole_length=0.0889,
                                     dipole_angle=1.5, pole_gap=0.014478),
                     BeamlineElement('DIPOLE', 0.037389, angle=11.25)])
    Mf = felsim_map([on_beam(dipole_wedge(0.01, angle=2.018, dipole_length=0.0889,
                                          dipole_angle=1.5, pole_gap=0.014478), 45.0),
                     on_beam(dipole(length=0.037389, angle=11.25), 45.0)])
    M = line_map(ad)
    assert np.all(np.abs(M - Mf) <= 1e-8 * np.maximum(1.0, np.abs(Mf)))


def test_dpw_without_bend_takes_the_dph_beside_it(caplog):
    """A generic DPW without dipole_length and dipole_angle takes them from
    the DPH beside it, the one before it first (a wedge after a DPH is its
    exit face), with one warning per wedge. With no DPH beside it, with only one of the two, or with a
    bend over no length, it raises."""
    wedge = dict(angle=5.0, pole_gap=0.014478)
    line = [BeamlineElement('DIPOLE_WEDGE', 0.01, **wedge),
            BeamlineElement('DIPOLE', 0.2, angle=10.0),
            BeamlineElement('DIPOLE_WEDGE', 0.01, **wedge),
            BeamlineElement('DIPOLE', 0.04064, angle=-4.0),
            BeamlineElement('DIPOLE_WEDGE', 0.01, **wedge)]
    bends = [(0.2, 10.0), None, (0.2, 10.0), None, (0.04064, -4.0)]
    Mf = felsim_map([on_beam(dipole_wedge(0.01, dipole_length=b[0], dipole_angle=b[1], **wedge)
                             if b else dipole(length=e.length, angle=e.parameters['angle']), 45.0)
                     for e, b in zip(line, bends)])
    ad = adapter()
    with caplog.at_level('WARNING', logger='rftrackAdapter'):
        ad.set_beamline(line)
        M = line_map(ad)
    assert np.all(np.abs(M - Mf) <= 1e-8 * np.maximum(1.0, np.abs(Mf)))
    assert caplog.text.count('has no dipole_length and dipole_angle') == 3

    with pytest.raises(ValueError, match='DIPOLE_WEDGE'):
        adapter().set_beamline([BeamlineElement('DRIFT', 0.1),
                                BeamlineElement('DIPOLE_WEDGE', 0.01, **wedge)])
    with pytest.raises(ValueError, match='only one'):
        adapter().set_beamline([BeamlineElement('DIPOLE_WEDGE', 0.01, dipole_length=0.2, **wedge),
                                line[1]])
    with pytest.raises(ValueError, match='positive'):
        adapter().set_beamline([BeamlineElement('DIPOLE_WEDGE', 0.01, dipole_length=0.0,
                                                dipole_angle=10.0, **wedge)])


# ── D3: triplets and chicanes ───────────────────────────────────────────────

def triplets(line):
    names = [type(e).__name__ for e in line]
    return [i for i in range(1, len(line) - 1)
            if names[i] == 'dipole' and names[i - 1] == names[i + 1] == 'dipole_wedge']


@pytest.mark.parametrize('path', [JSON_PATH, XLSX3_PATH], ids=['json', 'xlsx3'])
@pytest.mark.parametrize('E', [40.0, 45.0])
def test_triplets_and_chicanes(path, E):
    """Every DPW-DPH-DPW triplet, and the JSON line's chicane sections: each
    entry within 1e-9 absolute or 1e-8 relative of FELsim's matrix product,
    R56 within 1e-6. R51 to R54 come from a 1e-2 step: at 1e-3 one ulp of the
    arrival time over the 1.66 m FC1 section is already 1.1e-9 in R51."""
    line = latticeLoader.create_beamline(str(path))
    for e in line:
        e.setE(E)
    ad = adapter(E, lattice_path=str(path))
    ranges = [(i - 1, i + 2) for i in triplets(line)]
    assert len(ranges) == 18
    if path == JSON_PATH:
        ranges += [(6, 15), (23, 32), (46, 55), (66, 75), (83, 92), (99, 116), (121, 137)]
    for a, b in ranges:
        M, Mf = line_map(ad, a, b), felsim_map(line[a:b])
        M[4, :4] = line_map(ad, a, b, np.full(6, 1e-2))[4, :4]
        bound = np.maximum(1e-9, 1e-8 * np.abs(Mf))
        bound[4, 5] = np.inf
        assert np.all(np.abs(M - Mf) <= bound), (a, b)
        assert M[4, 5] == pytest.approx(Mf[4, 5], rel=1e-6), (a, b)


# ── D6 (i), (ii): off momentum ──────────────────────────────────────────────

def exact_sector(L, theta, x0, xp0, p_ratio):
    """Exact hard-edge sector bend in the horizontal plane: x [m] and x' [rad]
    at the exit face for a particle entering at x0, x'0 with momentum
    p_ratio P0 (theta > 0; +x away from the centre of curvature)."""
    from scipy.optimize import brentq
    rho0 = L / theta
    rho = rho0 * p_ratio
    a = math.atan(xp0)
    cx, cz = x0 - rho * math.cos(a), rho * math.sin(a)

    def pos(u):
        return cx + rho * math.cos(a - u), cz - rho * math.sin(a - u)

    def face(u):
        x, z = pos(u)
        return (x + rho0) * math.sin(theta) - z * math.cos(theta)

    u = brentq(face, -0.5 * theta, 3 * theta + 0.5)
    x, z = pos(u)
    return (x + rho0) * math.cos(theta) + z * math.sin(theta) - rho0, math.tan(a - u + theta)


@pytest.mark.parametrize('L, deg', [(0.037389, 11.25), (0.04064, -4.0), (0.0889, 1.5)])
def test_dph_off_momentum_against_the_exact_bend(L, deg):
    """Over |dK/K| <= 2 % and |x|, |x'| <= 5 mm, mrad, the SBend is the exact
    hard-edge bend within 1e-2 mm and mrad (measured 7.1e-3 at 11.25 deg;
    FELsim's chromatic dipole is off by 0.1 there, its K = h^2 for h0 h)."""
    ad = adapter()
    ad.set_beamline([BeamlineElement('DIPOLE', L, angle=deg)])
    X = np.array([[x0, xp0, 0.0, 0.0, 0.0, dk] for dk in (-20.0, -10.0, 0.0, 10.0, 20.0)
                  for x0, xp0 in ((0.0, 0.0), (5.0, 0.0), (0.0, 5.0), (-5.0, 3.0))])
    out = ad.simulate(X).final_particles
    P = ad.transform_coordinates(X, FEL, RFT)[:, 5] / ad._Pc
    s = math.copysign(1.0, deg)
    exact = np.array([exact_sector(L, abs(math.radians(deg)), s * x0 * 1e-3, s * xp0 * 1e-3, p)
                      for (x0, xp0), p in zip(X[:, :2], P)]) * s * 1e3
    assert np.abs(out[:, :2] - exact).max() <= 1e-2


def test_dpw_off_momentum_is_felsims_chromatic_wedge():
    """At dK/K = +-2 % the edge kicks equal FELsim's chromatic DPW (R ~ P)
    within 1e-5 relative (measured 1.4e-6)."""
    for eta, (dl, da), gap in ((11.25, (0.037389, 11.25), 0.0127), (2.018, (0.04064, -4.0), 0.014478)):
        ad = adapter()
        ad.set_beamline([BeamlineElement('DIPOLE_WEDGE', 0.01, angle=eta, dipole_length=dl,
                                         dipole_angle=da, pole_gap=gap)])
        w = on_beam(dipole_wedge(0.01, angle=eta, dipole_length=dl, dipole_angle=da, pole_gap=gap), 45.0)
        w.chromatic = True
        for dk in (-20.0, 20.0):
            X = np.zeros((3, 6))
            X[:, 5] = dk
            X[1, 0] = X[2, 2] = 1.0              # 1 mm in x, in y
            out, ref = ad.simulate(X).final_particles, w.useMatrice(X)
            kicks = [out[1, 1] - out[0, 1], out[2, 3] - out[0, 3]]
            assert kicks == pytest.approx([ref[1, 1] - ref[0, 1], ref[2, 3] - ref[0, 3]], rel=1e-5)


# ── D7: length and s ────────────────────────────────────────────────────────

def test_length_and_s_are_felsims():
    ad = adapter(lattice_path=str(JSON_PATH))
    lengths = [e.length for e in latticeLoader.create_beamline(str(JSON_PATH))]
    assert ad._lattice.get_length() == pytest.approx(sum(lengths), abs=1e-9)
    evolution = ad.collect_evolution(np.zeros((2, 6)))
    np.testing.assert_allclose([e.s_end for e in evolution.elements], np.cumsum(lengths),
                               rtol=0, atol=1e-12)


@pytest.mark.parametrize('physical', [False, True])
def test_loss_inside_a_dph_is_booked_to_it(physical):
    """A particle that reaches the pole gap inside the DPH is lost there, in
    the table and in collect_evolution alike."""
    ad = RFTrackAdapter(beam_energy=45.0, aperture=0.00724)
    if physical:
        ad.enable_physical_apertures()
    ad.set_beamline([BeamlineElement('DIPOLE_WEDGE', 0.01, angle=1.5, dipole_length=0.0889,
                                     dipole_angle=1.5, pole_gap=0.014478),
                     BeamlineElement('DIPOLE', 0.0889, angle=1.5, pole_gap=0.014478),
                     BeamlineElement('DIPOLE_WEDGE', 0.01, angle=1.5, dipole_length=0.0889,
                                     dipole_angle=1.5, pole_gap=0.014478),
                     BeamlineElement('DRIFT', 0.2)])
    p = np.zeros((2, 6))
    p[0, 2], p[0, 3] = 7.0, 5.0          # y = 7.44 mm at the DPH exit
    result = ad.simulate(p)
    assert result.metadata['survivor_ids'].tolist() == [1]
    assert result.metadata['lost_elements'].tolist() == [1]
    assert result.metadata['lost_particles'][0, 6] == pytest.approx(0.01, abs=1e-12)
    n_lost = [e.parameters['n_lost'] for e in ad.collect_evolution(p).elements]
    assert n_lost == [0, 1, 0, 0]


def native_apertures(ad, i):
    """(shape, ax, ay) [m] of the RF-Track elements of beamline[i]."""
    return [(e.get_aperture_shape(), e.get_aperture_x(), e.get_aperture_y())
            for e in (ad._lattice[int(k)] for k in np.flatnonzero(ad._lattice_index == i))]


@pytest.mark.parametrize('path', [JSON_PATH, XLSX_PATH], ids=['json', 'xlsx'])
def test_dipole_apertures_are_felsims(path):
    """With physical apertures every dipole is cut where FELsim's
    apply_aperture cuts it: in y only, at its own pole_gap/2 (12.7 mm gaps
    in the FC1 and FC2 magnets of the xlsx line, no cut in the JSON line's
    DPHs, which carry no pole gap), a rectangle without bound in x."""
    line = latticeLoader.create_beamline(str(path))
    ad = adapter(lattice_path=str(path))
    ad.enable_physical_apertures()
    cuts = set()
    for i, e in enumerate(line):
        if type(e).__name__ not in ('dipole', 'dipole_wedge'):
            continue
        ax, ay = (math.inf if a is None else a * 1e-3 for a in (e.aperture_x, e.aperture_y))
        got = native_apertures(ad, i)
        if ax == ay == math.inf:
            assert got and all(g[0] == 'none' for g in got), (i, got)
        else:
            assert got and all(g[0] == 'rectangular' and g[1:] == pytest.approx((ax, ay), rel=1e-12)
                               for g in got), (i, got, ax, ay)
        cuts.add((type(e).__name__, round(ay * 1e3, 3)))
    assert ('dipole', 6.35) in cuts if path == XLSX_PATH else ('dipole', math.inf) in cuts


def test_dipole_apertures_cut_in_y_only():
    """With physical apertures a particle 40 mm off axis in x crosses the
    dipoles, which a drift's beam pipe would cut, and one 6.5 mm off in y is
    lost in the entrance wedge of a 12.7 mm gap. Without, dipoles have the
    default aperture, the circle of every other element."""
    wedge = BeamlineElement('DIPOLE_WEDGE', 0.01, angle=1.5, dipole_length=0.0889,
                            dipole_angle=1.5, pole_gap=0.0127)
    line = [wedge, BeamlineElement('DIPOLE', 0.0889, angle=1.5, pole_gap=0.0127), wedge]
    ad = RFTrackAdapter(beam_energy=45.0)
    ad.set_beamline(line + [BeamlineElement('DRIFT', 0.1)])
    assert all(native_apertures(ad, i) and set(native_apertures(ad, i)) == {('circular', 0.05, 0.05)}
               for i in range(4))
    ad.enable_physical_apertures()
    assert native_apertures(ad, 3) == [('circular', 0.0127, 0.0127)]
    ad.set_beamline(line)
    p = np.zeros((3, 6))
    p[0, 0], p[1, 2], p[2, 2] = 40.0, 6.2, 6.5
    metadata = ad.simulate(p).metadata
    assert sorted(metadata['survivor_ids']) == [0, 1]
    assert metadata['lost_elements'].tolist() == [0]


def test_physical_dipole_apertures_can_be_set():
    """enable_physical_apertures' dipole_half_gap and dipole_half_width
    apply to every dipole, a DPH without a pole gap included."""
    ad = RFTrackAdapter(beam_energy=45.0)
    ad.enable_physical_apertures(dipole_half_gap=0.005, dipole_half_width=0.03)
    ad.set_beamline([BeamlineElement('DIPOLE_WEDGE', 0.01, angle=1.5, dipole_length=0.0889,
                                     dipole_angle=1.5, pole_gap=0.0127),
                     BeamlineElement('DIPOLE', 0.0889, angle=1.5)])
    assert set(native_apertures(ad, 0) + native_apertures(ad, 1)) == {('rectangular', 0.03, 0.005)}


# ── D11, D12: compatibility, one Lattice ────────────────────────────────────

def test_dipole_slices_is_deprecated_and_inert():
    """dipole_slices warns, in a script too (a FutureWarning, which Python
    shows by default), has no effect, and is quiet when not given."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        RFTrackAdapter(beam_energy=45.0, aperture=None)
        RFTrackAdapter(beam_energy=45.0, aperture=None, dipole_slices=None)
    assert not [w for w in caught if 'dipole_slices' in str(w.message)]
    script = subprocess.run([sys.executable, '-c', 'from rftrackAdapter import RFTrackAdapter; '
                             'RFTrackAdapter(beam_energy=45.0, dipole_slices=4)'],
                            cwd=REPO / 'backend', capture_output=True, text=True, timeout=120)
    assert script.returncode == 0 and 'dipole_slices is deprecated' in script.stderr

    with pytest.warns(FutureWarning, match='dipole_slices'):
        ad = RFTrackAdapter(beam_energy=45.0, aperture=None, dipole_slices=0)
    plain = RFTrackAdapter(beam_energy=45.0, aperture=None)
    line = [BeamlineElement('DIPOLE', 0.2, angle=10.0), BeamlineElement('DRIFT', 0.3)]
    for a in (ad, plain):
        a.set_beamline(line)
    X = np.array([[1.0, 0.1, 1.0, 0.1, 0.5, 5.0]])
    np.testing.assert_array_equal(ad.simulate(X).final_particles, plain.simulate(X).final_particles)

    mc = MultiCodeSimulator(sections=[SimSection('s', 'rftrack', (0, 4), {'dipole_slices': 10})],
                            lattice_path=str(JSON_PATH), beam_energy=45.0)
    with pytest.warns(FutureWarning, match='dipole_slices'):
        mc.simulate(np.zeros((2, 6)))
    sim, = mc._simulators.values()
    assert sim.dipole_slices == 10


def test_one_lattice():
    """No trace of the analytical dipole path, and simulate tracks
    self._lattice only."""
    gone = {'DPW_THIN_LENS_LENGTH', 'SECTOR_BEND_MIN_ANGLE', '_annotate_dipole_edges',
            '_compute_dpw_kicks', '_build_sliced_dipole', '_apply_sector_bend_correction',
            '_apply_dpw_edge_kick', '_track_segmented', '_track_range', '_track_analytical_dipole',
            '_track_analytical_dpw', '_annotate_analytical_dipoles', '_analytical_dipole',
            '_analytical_dpw', 'dpw_kick_x', 'dpw_kick_y', 'dipole_K0'}
    names = set()
    for n in ast.walk(ast.parse(Path(rftrackAdapter.__file__).read_text())):
        if isinstance(n, ast.Name):
            names.add(n.id)
        elif isinstance(n, ast.Attribute):
            names.add(n.attr)
        elif isinstance(n, ast.FunctionDef):
            names.add(n.name)
        elif isinstance(n, ast.Constant) and isinstance(n.value, str):
            names.add(n.value)
    assert not gone & names

    ad = adapter(lattice_path=str(JSON_PATH))
    built, tracked = [], []
    new_lattice, track_rows = ad._new_lattice, ad._track_rows
    ad._new_lattice = lambda *a: built.append(a) or new_lattice(*a)
    ad._track_rows = lambda lat, *a: tracked.append(lat) or track_rows(lat, *a)
    ad.simulate(np.zeros((3, 6)))
    assert built == [] and tracked == [ad._lattice]


@pytest.mark.parametrize('setter', ['set_beam_energy', 'set_particle_type'])
def test_a_new_reference_momentum_rebuilds_the_magnets(setter):
    """Dipoles, wedges and sextupoles are set for the reference momentum:
    after set_beam_energy (45 -> 40 MeV) or set_particle_type (electron ->
    positron) a built adapter tracks as a fresh one."""
    line = [BeamlineElement('DIPOLE_WEDGE', 0.01, angle=2.018, dipole_length=0.037389,
                            dipole_angle=11.25, pole_gap=0.0127),
            BeamlineElement('DIPOLE', 0.037389, angle=11.25),
            BeamlineElement('DIPOLE_WEDGE', 0.01, angle=2.018, dipole_length=0.037389,
                            dipole_angle=11.25, pole_gap=0.0127),
            BeamlineElement('DRIFT', 0.3), BeamlineElement('SEXTUPOLE', 0.1, strength=50.0),
            BeamlineElement('DRIFT', 0.2)]
    ad = adapter(45.0)
    ad.set_beamline(line)
    if setter == 'set_beam_energy':
        ad.set_beam_energy(40.0)
        fresh = adapter(40.0)
    else:
        ad.set_particle_type(ad.particle_mass, 1.0)
        fresh = adapter(45.0, 'positron')
    fresh.set_beamline(line)
    X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.5, -1.0, 0.3, 0.2, 5.0]])
    np.testing.assert_allclose(ad.simulate(X).final_particles, fresh.simulate(X).final_particles,
                               rtol=0, atol=1e-12)


def test_magnets_after_an_rf_cavity_warn_once(caplog):
    """Dipoles and sextupoles are set for the initial momentum; downstream of
    a cavity the adapter says so, once."""
    rfc = adapter(1.0, lattice_path=str(JSON_PATH.with_name('slac_linac.json'))).beamline[0]
    ad = adapter(1.0)
    with caplog.at_level('WARNING', logger='rftrackAdapter'):
        ad.set_beamline([BeamlineElement('DIPOLE', 0.1, angle=1.0), BeamlineElement('DRIFT', 0.1)])
        assert 'RF cavity' not in caplog.text
        for _ in range(2):
            ad.set_beamline([rfc, BeamlineElement('DRIFT', 0.1),
                             BeamlineElement('DIPOLE', 0.1, angle=1.0)])
    assert caplog.text.count('after an RF cavity') == 1
    assert 'DIPOLE' in caplog.text
