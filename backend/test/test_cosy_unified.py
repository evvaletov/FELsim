"""
Test unified plotting with COSY adapter.

Author: Eremey Valetov
"""
import numpy as np
import pytest
from pathlib import Path

# Check COSY availability
try:
    from cosyAdapter import COSYAdapter, _COSY_AVAILABLE
    if not _COSY_AVAILABLE:
        pytest.skip("COSY components not available", allow_module_level=True)
except ImportError as e:
    pytest.skip(f"Import error: {e}", allow_module_level=True)

try:
    from beamEvolution import BeamEvolution
    from evolutionPlotter import EvolutionPlotter
    from simulatorBase import CoordinateSystem, SimulationMode
except ImportError as e:
    pytest.skip(f"Import error: {e}", allow_module_level=True)

try:
    from felsimAdapter import FELsimAdapter
    from ebeam import beam
    FELSIM_AVAILABLE = True
except ImportError:
    FELSIM_AVAILABLE = False

EXCEL_PATH = Path(__file__).resolve().parent.parent.parent / "beam_excel" / "Beamline_elements_3.xlsx"
BEAM_ENERGY = 45.0


def _create_test_particles(n=1000):
    ebeam = beam()
    std_dev = [1.0, 0.1, 1.0, 0.1, 1.0, 0.5]
    return ebeam.gen_6d_gaussian(0, std_dev, n)


@pytest.fixture(scope="module")
def cosy_sim():
    if not EXCEL_PATH.exists():
        pytest.skip(f"Excel file not found: {EXCEL_PATH}")
    return COSYAdapter(
        excel_path=str(EXCEL_PATH),
        mode='particle_tracking',
        debug=True
    )


@pytest.fixture(scope="module")
def particles_felsim():
    return _create_test_particles(500)


@pytest.fixture(scope="module")
def cosy_evolution(cosy_sim, particles_felsim):
    return cosy_sim.collect_evolution(
        particles_felsim.copy(),
        checkpoint_elements=list(range(1, 21))
    )


class TestCOSYAdapter:
    def test_creation(self, cosy_sim):
        assert cosy_sim.name is not None
        assert cosy_sim.native_coordinates == CoordinateSystem.COSY

    def test_beamline_loaded(self, cosy_sim):
        beamline = cosy_sim.get_beamline()
        assert len(beamline) > 0

    def test_coordinate_roundtrip(self, cosy_sim, particles_felsim):
        native_sim = cosy_sim.get_native_simulator()
        particles_cosy = native_sim.transform_to_cosy_coordinates(
            particles_felsim, energy=BEAM_ENERGY
        )
        particles_back = native_sim.transform_from_cosy_coordinates(
            particles_cosy, energy=BEAM_ENERGY
        )
        assert np.max(np.abs(particles_felsim - particles_back)) < 1e-10


class TestCOSYEvolution:
    def test_data_consistency(self, cosy_evolution):
        ev = cosy_evolution
        assert len(ev.s_positions) == len(ev.particles)
        assert len(ev.s_positions) == len(ev.twiss)
        assert ev.total_length > 0

    def test_twiss_fields(self, cosy_evolution):
        ev = cosy_evolution
        assert len(ev.s_positions) > 0
        s0 = ev.s_positions[0]
        assert s0 in ev.twiss
        twiss_x = ev.twiss[s0].get('x', {})
        for field in ['beta', 'alpha', 'gamma', 'emittance']:
            assert field in twiss_x, f"Missing Twiss field: {field}"

    def test_s_positions_monotonic(self, cosy_evolution):
        s_sorted = sorted(cosy_evolution.s_positions)
        if len(s_sorted) > 1:
            spacings = np.diff(s_sorted)
            assert np.all(spacings > 0)

    def test_twiss_dataframe(self, cosy_evolution):
        df = cosy_evolution.get_twiss_evolution()
        assert df.shape[0] > 0
        assert 's' in df.columns
        assert 'beta_x' in df.columns
        assert df['s'].min() >= 0

    def test_selective_checkpoints(self, cosy_sim, particles_felsim):
        beamline = cosy_sim.get_beamline()
        if len(beamline) < 5:
            pytest.skip("Beamline too short")

        checkpoint_list = list(range(1, 21))
        evolution = cosy_sim.collect_evolution(
            particles_felsim.copy(),
            checkpoint_elements=checkpoint_list
        )
        expected = len(checkpoint_list) + 1  # +1 for s=0
        assert len(evolution.s_positions) == expected


@pytest.mark.visual
class TestCOSYPlotting:
    def test_plotter(self, cosy_evolution):
        plotter = EvolutionPlotter(axis_mode='local')
        plotter.plot(
            cosy_evolution,
            show_phase_space=True,
            show_envelope=True,
            show_schematic=True,
            interactive=True,
            scatter=False
        )
