"""
Test unified plotting with FELsim adapter using the real beamline from Excel.
Limited to first N elements to avoid numerical instabilities.

Author: Eremey Valetov
"""
import numpy as np
import pytest
from pathlib import Path

from felsimAdapter import FELsimAdapter
from evolutionPlotter import EvolutionPlotter
from ebeam import beam

EXCEL_PATH = Path(__file__).resolve().parent.parent.parent / "beam_excel" / "Beamline_elements.xlsx"
MAX_ELEMENTS = 50
BEAM_ENERGY_MEV = 45.0


def _load_beamline(max_elements=MAX_ELEMENTS):
    if not EXCEL_PATH.exists():
        pytest.skip(f"Excel file not found: {EXCEL_PATH}")
    temp_sim = FELsimAdapter(excel_path=str(EXCEL_PATH))
    temp_sim.set_beam_energy(BEAM_ENERGY_MEV)
    limited = temp_sim.get_native_beamline()[:max_elements]
    for el in limited:
        el.setE(BEAM_ENERGY_MEV)
    return limited


def _create_test_particles(n=1000):
    ebeam = beam()
    std_dev = np.array([1.0, 0.1, 1.0, 0.1, 1.0, 0.5])
    return ebeam.gen_6d_gaussian(0, std_dev, n)


@pytest.fixture(scope="module")
def felsim_beamline():
    return _load_beamline()


@pytest.fixture(scope="module")
def felsim_sim(felsim_beamline):
    sim = FELsimAdapter()
    sim.set_beam_energy(BEAM_ENERGY_MEV)
    sim._native_beamline = felsim_beamline
    return sim


@pytest.fixture(scope="module")
def felsim_particles():
    return _create_test_particles(500)


@pytest.fixture(scope="module")
def felsim_evolution(felsim_sim, felsim_particles):
    return felsim_sim.collect_evolution(felsim_particles.copy(), interval=0.05)


class TestFELsimEvolution:
    def test_evolution_collection(self, felsim_evolution):
        ev = felsim_evolution
        assert len(ev.s_positions) > 0
        assert len(ev.s_positions) == len(ev.particles)
        assert len(ev.s_positions) == len(ev.twiss)

    def test_twiss_fields(self, felsim_evolution):
        ev = felsim_evolution
        s0 = ev.s_positions[0]
        twiss_x = ev.twiss[s0]['x']
        for field in ['beta', 'alpha', 'gamma', 'emittance', 'dispersion']:
            assert field in twiss_x, f"Missing Twiss field: {field}"

    def test_twiss_dataframe(self, felsim_evolution):
        df = felsim_evolution.get_twiss_evolution()
        assert df.shape[0] > 0
        assert 's' in df.columns
        assert 'beta_x' in df.columns
        assert 'dispersion_x' in df.columns
        assert df['beta_x'].min() >= 0


@pytest.mark.visual
class TestFELsimPlotting:
    def test_plotter(self, felsim_evolution):
        plotter = EvolutionPlotter(axis_mode='local')
        plotter.plot(
            felsim_evolution,
            show_phase_space=True,
            show_envelope=True,
            show_schematic=True,
            interactive=True,
            scatter=False
        )

    def test_axis_modes(self, felsim_sim):
        particles = _create_test_particles(1000)
        ev = felsim_sim.collect_evolution(particles, interval=0.03)
        EvolutionPlotter(axis_mode='local').plot(ev)
        EvolutionPlotter(axis_mode='global').plot(ev)
