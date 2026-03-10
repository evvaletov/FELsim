"""I7: Multi-code simulation framework tests.

Validates that MultiCodeSimulator correctly chains simulator sections
with coordinate transforms at handoff points.

Author: Eremey Valetov
"""

import sys
from pathlib import Path

import numpy as np
import pytest

_BACKEND = Path(__file__).resolve().parent.parent
_PROJECT_ROOT = _BACKEND.parent
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))

JSON_PATH = _PROJECT_ROOT / "var" / "UH_FEL_beamline.json"

from multiCodeSimulator import MultiCodeSimulator, SimSection
from simulatorBase import CoordinateSystem, SimulationResult
from simulatorFactory import SimulatorFactory, CoordinateTransformer


# ── Helpers ──────────────────────────────────────────────────────────────

def load_beamline(path):
    import latticeLoader
    return latticeLoader.create_beamline(str(path))


def propagate_single_pass(beamline, particles, energy_mev=40.0):
    """Single-pass FELsim tracking through full beamline."""
    state = particles.copy()
    for elem in beamline:
        elem.setE(energy_mev)
        state = np.array(elem.useMatrice(state))
    return state


# ── Fixtures ─────────────────────────────────────────────────────────────

@pytest.fixture
def beamline():
    if not JSON_PATH.exists():
        pytest.skip("UH_FEL_beamline.json not found")
    return load_beamline(JSON_PATH)


@pytest.fixture
def particles():
    np.random.seed(42)
    return np.random.randn(100, 6) * [1e-3, 1e-4, 1e-3, 1e-4, 1e-3, 0.005]


# ── Unit tests ───────────────────────────────────────────────────────────

class TestSimSection:
    """SimSection dataclass basics."""

    def test_creation(self):
        s = SimSection("prefix", "felsim", (0, 60))
        assert s.name == "prefix"
        assert s.simulator_key == "felsim"
        assert s.element_range == (0, 60)
        assert s.config == {}

    def test_with_config(self):
        s = SimSection("suffix", "cosy", (60, 137), config={"mode": "particle_tracking"})
        assert s.config["mode"] == "particle_tracking"


class TestMultiCodeInit:
    """MultiCodeSimulator initialisation."""

    def test_no_sections_raises(self):
        mc = MultiCodeSimulator(sections=[])
        with pytest.raises(ValueError, match="No sections"):
            mc.simulate(np.zeros((10, 6)))

    def test_no_particles_raises(self):
        mc = MultiCodeSimulator(sections=[
            SimSection("a", "felsim", (0, 10))
        ])
        with pytest.raises(ValueError, match="particles required"):
            mc.simulate()

    def test_from_config(self):
        config = {
            "beam_energy_mev": 40.0,
            "sections": [
                {"name": "first", "simulator": "felsim", "elements": [0, 60]},
                {"name": "second", "simulator": "felsim", "elements": [60, 137]},
            ]
        }
        mc = MultiCodeSimulator.from_config(config)
        assert len(mc.sections) == 2
        assert mc.beam_energy == 40.0


class TestCoordinateTransformerRoundtrips:
    """Round-trip tests for all 6 pairwise coordinate transforms."""

    @pytest.fixture
    def felsim_particles(self):
        np.random.seed(99)
        return np.random.randn(50, 6) * [1.0, 0.1, 1.0, 0.1, 5.0, 1.0]

    @pytest.mark.parametrize("from_sys,to_sys", [
        (CoordinateSystem.FELSIM, CoordinateSystem.COSY),
        (CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK),
        (CoordinateSystem.COSY, CoordinateSystem.RFTRACK),
    ])
    def test_roundtrip(self, felsim_particles, from_sys, to_sys):
        energy = 45.0

        if from_sys == CoordinateSystem.FELSIM:
            original = felsim_particles
        elif from_sys == CoordinateSystem.COSY:
            original = CoordinateTransformer.transform(
                felsim_particles, CoordinateSystem.FELSIM, CoordinateSystem.COSY, energy
            )
        else:
            original = CoordinateTransformer.transform(
                felsim_particles, CoordinateSystem.FELSIM, from_sys, energy
            )

        intermediate = CoordinateTransformer.transform(original, from_sys, to_sys, energy)
        recovered = CoordinateTransformer.transform(intermediate, to_sys, from_sys, energy)

        np.testing.assert_allclose(recovered, original, atol=1e-10, rtol=1e-10)

    def test_identity_transform(self, felsim_particles):
        result = CoordinateTransformer.transform(
            felsim_particles, CoordinateSystem.FELSIM, CoordinateSystem.FELSIM, 45.0
        )
        np.testing.assert_array_equal(result, felsim_particles)


class TestFELsimSplitEquivalence:
    """Core MVP test: splitting FELsim beamline into 2 sections
    via MultiCodeSimulator must produce identical results to single-pass."""

    SPLIT_POINT = 60
    ENERGY = 40.0

    def test_two_section_matches_single_pass(self, beamline, particles):
        """MultiCode(felsim+felsim) == single FELsim pass."""
        n_elem = len(beamline)

        # Single-pass reference
        ref = propagate_single_pass(beamline, particles, self.ENERGY)

        # Multi-code: two FELsim sections
        mc = MultiCodeSimulator(
            sections=[
                SimSection("prefix", "felsim", (0, self.SPLIT_POINT)),
                SimSection("suffix", "felsim", (self.SPLIT_POINT, n_elem)),
            ],
            lattice_path=str(JSON_PATH),
            beam_energy=self.ENERGY,
        )

        result = mc.simulate(particles=particles)
        assert result.success
        assert result.final_particles is not None

        np.testing.assert_allclose(
            result.final_particles, ref, atol=1e-12,
            err_msg="MultiCode 2-section FELsim != single-pass FELsim"
        )

    def test_three_section_matches_single_pass(self, beamline, particles):
        """MultiCode(felsim+felsim+felsim) == single FELsim pass."""
        n_elem = len(beamline)
        sp1, sp2 = 40, 90

        ref = propagate_single_pass(beamline, particles, self.ENERGY)

        mc = MultiCodeSimulator(
            sections=[
                SimSection("sec1", "felsim", (0, sp1)),
                SimSection("sec2", "felsim", (sp1, sp2)),
                SimSection("sec3", "felsim", (sp2, n_elem)),
            ],
            lattice_path=str(JSON_PATH),
            beam_energy=self.ENERGY,
        )

        result = mc.simulate(particles=particles)
        assert result.success
        np.testing.assert_allclose(
            result.final_particles, ref, atol=1e-12,
            err_msg="MultiCode 3-section FELsim != single-pass FELsim"
        )

    def test_metadata_sections(self, beamline, particles):
        """Result metadata should record section info."""
        n_elem = len(beamline)

        mc = MultiCodeSimulator(
            sections=[
                SimSection("prefix", "felsim", (0, self.SPLIT_POINT)),
                SimSection("suffix", "felsim", (self.SPLIT_POINT, n_elem)),
            ],
            lattice_path=str(JSON_PATH),
            beam_energy=self.ENERGY,
        )

        result = mc.simulate(particles=particles)
        assert result.metadata['num_sections'] == 2
        assert len(result.metadata['sections']) == 2
        assert result.metadata['sections'][0]['name'] == 'prefix'
        assert result.metadata['sections'][1]['name'] == 'suffix'

    def test_particle_count_preserved(self, beamline, particles):
        """All particles should survive through the pipeline."""
        n_elem = len(beamline)

        mc = MultiCodeSimulator(
            sections=[
                SimSection("prefix", "felsim", (0, self.SPLIT_POINT)),
                SimSection("suffix", "felsim", (self.SPLIT_POINT, n_elem)),
            ],
            lattice_path=str(JSON_PATH),
            beam_energy=self.ENERGY,
        )

        result = mc.simulate(particles=particles)
        assert result.final_particles.shape == particles.shape


class TestFactoryRegistration:
    """MultiCodeSimulator accessible via SimulatorFactory."""

    def test_multicode_in_available(self):
        available = SimulatorFactory.get_available_simulators()
        assert 'multicode' in available

    def test_create_multicode(self):
        mc = SimulatorFactory.create('multicode', sections=[])
        assert isinstance(mc, MultiCodeSimulator)
