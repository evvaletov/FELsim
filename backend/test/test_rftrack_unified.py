"""
Test unified plotting with RF-Track adapter.

Author: Eremey Valetov
"""
import numpy as np
import pytest
from pathlib import Path

# Check RF-Track availability
try:
    from rftrackAdapter import RFTrackAdapter, _RFTRACK_AVAILABLE
    if not _RFTRACK_AVAILABLE:
        pytest.skip("RF-Track not installed", allow_module_level=True)
except ImportError as e:
    pytest.skip(f"Import error: {e}", allow_module_level=True)

from beamEvolution import BeamEvolution
from evolutionPlotter import EvolutionPlotter
from simulatorBase import CoordinateSystem, SimulationMode, BeamlineElement
from simulatorFactory import SimulatorFactory

try:
    from ebeam import beam
    FELSIM_AVAILABLE = True
except ImportError:
    FELSIM_AVAILABLE = False

EXCEL_PATH = Path("../../beam_excel/Beamline_elements_3.xlsx")
BEAM_ENERGY = 45.0
MAX_ELEMENTS = 30


def _create_test_particles(n=1000):
    if FELSIM_AVAILABLE:
        ebeam = beam()
        std_dev = [1.0, 0.1, 1.0, 0.1, 1.0, 0.5]
        return ebeam.gen_6d_gaussian(0, std_dev, n)
    std_dev = np.array([1.0, 0.1, 1.0, 0.1, 1.0, 0.5])
    return np.random.randn(n, 6) * std_dev


def _create_manual_beamline(sim):
    elements = []
    for _ in range(5):
        elements.extend([
            BeamlineElement('DRIFT', 0.2),
            BeamlineElement('QUAD_F', 0.08, current=2.0),
            BeamlineElement('DRIFT', 0.4),
            BeamlineElement('QUAD_D', 0.08, current=2.0),
            BeamlineElement('DRIFT', 0.2),
        ])
    sim.set_beamline(elements)


@pytest.fixture(scope="module")
def rftrack_sim():
    if EXCEL_PATH.exists():
        sim = RFTrackAdapter(
            excel_path=str(EXCEL_PATH),
            beam_energy=BEAM_ENERGY
        )
        if len(sim.beamline) > MAX_ELEMENTS:
            sim.beamline = sim.beamline[:MAX_ELEMENTS]
            sim._build_lattice()
    else:
        sim = RFTrackAdapter(beam_energy=BEAM_ENERGY)
        _create_manual_beamline(sim)
    return sim


@pytest.fixture(scope="module")
def rftrack_particles():
    return _create_test_particles(500)


@pytest.fixture(scope="module")
def rftrack_evolution(rftrack_sim, rftrack_particles):
    n_elements = len(rftrack_sim.beamline)
    checkpoint_list = list(range(min(n_elements, 20)))
    return rftrack_sim.collect_evolution(
        rftrack_particles.copy(),
        checkpoint_elements=checkpoint_list
    )


class TestRFTrackAdapter:
    def test_creation(self, rftrack_sim):
        assert rftrack_sim.name == "RF-Track"
        assert rftrack_sim.native_coordinates == CoordinateSystem.RFTRACK

    def test_particle_tracking(self, rftrack_sim, rftrack_particles):
        result = rftrack_sim.simulate(rftrack_particles.copy())
        assert result.success
        assert result.final_particles is not None
        assert result.metadata['num_good'] > 0


class TestRFTrackEvolution:
    def test_data_consistency(self, rftrack_evolution):
        ev = rftrack_evolution
        assert len(ev.s_positions) == len(ev.particles)
        assert len(ev.s_positions) == len(ev.twiss)
        assert ev.total_length > 0

    def test_twiss_fields(self, rftrack_evolution):
        ev = rftrack_evolution
        assert len(ev.s_positions) > 0
        s0 = ev.s_positions[0]
        assert s0 in ev.twiss
        twiss_x = ev.twiss[s0].get('x', {})
        for field in ['beta', 'alpha', 'gamma', 'emittance']:
            assert field in twiss_x, f"Missing Twiss field: {field}"

    def test_s_positions_monotonic(self, rftrack_evolution):
        s_sorted = sorted(rftrack_evolution.s_positions)
        if len(s_sorted) > 1:
            spacings = np.diff(s_sorted)
            assert np.all(spacings > 0)

    def test_twiss_dataframe(self, rftrack_evolution):
        df = rftrack_evolution.get_twiss_evolution()
        assert df.shape[0] > 0
        assert 's' in df.columns
        assert 'beta_x' in df.columns


@pytest.mark.visual
class TestRFTrackPlotting:
    def test_plotter(self, rftrack_evolution):
        plotter = EvolutionPlotter(axis_mode='local')
        plotter.plot(
            rftrack_evolution,
            show_phase_space=True,
            show_envelope=True,
            show_schematic=True,
            interactive=True,
            scatter=False,
            envelope_ylim=(0.0, 1.4)
        )
