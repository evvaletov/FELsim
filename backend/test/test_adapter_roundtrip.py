"""
C4 V&V: Adapter round-trip and regression tests.

Tests that all three lattice formats (Excel, JSON, YAML) produce identical
simulation results, and that key results remain stable across code changes.

Author: Eremey Valetov
"""

import sys
import os
from pathlib import Path

import numpy as np
import pytest

_BACKEND = Path(__file__).resolve().parent.parent
_PROJECT_ROOT = _BACKEND.parent
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))

EXCEL_PATH = _PROJECT_ROOT / "beam_excel" / "Beamline_elements.xlsx"
JSON_PATH = _PROJECT_ROOT / "var" / "UH_FEL_beamline.json"
YAML_PATH = _PROJECT_ROOT / "var" / "UH_FEL_beamline.yaml"

from beamline import driftLattice, qpfLattice, qpdLattice, dipole, dipole_wedge


# ── Helpers ──────────────────────────────────────────────────────────────

def load_beamline(path):
    """Load beamline via unified latticeLoader."""
    import latticeLoader
    return latticeLoader.create_beamline(str(path))


def cumulative_matrix(beamline):
    """Compute the full beamline transfer matrix M = M_N @ ... @ M_2 @ M_1."""
    M = np.eye(6)
    for elem in beamline:
        M = elem._compute_numeric_matrix() @ M
    return M


def propagate_particles(beamline, particles):
    """Track particles through beamline, return final state."""
    state = particles.copy()
    for elem in beamline:
        state = np.array(elem.useMatrice(state))
    return state


def set_energy(beamline, KE):
    """Set kinetic energy for all elements."""
    for elem in beamline:
        elem.setE(KE)


# ── Format availability ──────────────────────────────────────────────────

try:
    import tracked_dict  # noqa: F401
    _HAS_TRACKED_DICT = True
except ImportError:
    _HAS_TRACKED_DICT = False

_YAML_MARK = pytest.mark.skipif(
    not _HAS_TRACKED_DICT,
    reason="tracked_dict requires Python >=3.10",
)

ALL_PATHS = [
    pytest.param(EXCEL_PATH, id="Excel"),
    pytest.param(JSON_PATH, id="JSON"),
    pytest.param(YAML_PATH, id="YAML", marks=_YAML_MARK),
]


def available_paths():
    return [p for p in ALL_PATHS if p.values[0].exists()]


# ── Basic loading tests ──────────────────────────────────────────────────

class TestLoading:
    @pytest.mark.parametrize("path", ALL_PATHS)
    def test_loads_nonempty(self, path):
        if not path.exists():
            pytest.skip(f"File not found: {path}")
        bl = load_beamline(path)
        assert len(bl) > 100  # UH FEL has 118+ elements

    @pytest.mark.parametrize("path", ALL_PATHS)
    def test_all_elements_have_positive_length(self, path):
        if not path.exists():
            pytest.skip(f"File not found: {path}")
        bl = load_beamline(path)
        for i, elem in enumerate(bl):
            assert elem.length > 0, f"Element {i} ({type(elem).__name__}) has length {elem.length}"

    @pytest.mark.parametrize("path", ALL_PATHS)
    def test_element_names_populated(self, path):
        if not path.exists():
            pytest.skip(f"File not found: {path}")
        bl = load_beamline(path)
        named = sum(1 for e in bl if e.name is not None)
        # At least the quads and dipoles should have names
        assert named > 20, f"Only {named} elements have names"


# ── Cross-format equivalence ─────────────────────────────────────────────

class TestCrossFormatEquivalence:
    """Verify that all three formats produce identical beamline physics."""

    @pytest.fixture
    def beamlines(self):
        """Load beamlines from all available formats."""
        bls = {}
        _paths = [EXCEL_PATH, JSON_PATH]
        if _HAS_TRACKED_DICT:
            _paths.append(YAML_PATH)
        for path in _paths:
            if path.exists():
                bls[path.stem] = load_beamline(path)
        if len(bls) < 2:
            pytest.skip("Need at least 2 formats for cross-comparison")
        return bls

    def test_total_length(self, beamlines):
        """Total beamline length should be identical across formats."""
        lengths = {name: sum(e.length for e in bl) for name, bl in beamlines.items()}
        ref_name, ref_len = next(iter(lengths.items()))
        for name, length in lengths.items():
            assert abs(length - ref_len) < 1e-6, (
                f"Total length mismatch: {ref_name}={ref_len:.6f} vs {name}={length:.6f}"
            )

    def test_cumulative_matrix(self, beamlines):
        """Full beamline transfer matrix should be identical across formats."""
        matrices = {name: cumulative_matrix(bl) for name, bl in beamlines.items()}
        ref_name, ref_M = next(iter(matrices.items()))
        for name, M in matrices.items():
            diff = np.max(np.abs(M - ref_M))
            assert diff < 1e-10, (
                f"Cumulative matrix diff ({ref_name} vs {name}): {diff:.2e}"
            )

    def test_multi_particle_propagation(self, beamlines):
        """Multiple particles should produce identical output across formats."""
        np.random.seed(42)
        particles = np.random.randn(10, 6) * [1e-3, 1e-4, 1e-3, 1e-4, 1e-3, 0.005]

        results = {}
        for name, bl in beamlines.items():
            results[name] = propagate_particles(bl, particles)

        ref_name, ref_out = next(iter(results.items()))
        for name, out in results.items():
            diff = np.max(np.abs(out - ref_out))
            assert diff < 1e-10, (
                f"Propagation diff ({ref_name} vs {name}): {diff:.2e}"
            )

    def test_element_type_counts(self, beamlines):
        """Element type counts should be close across formats.

        Excel may have one extra drift from drift splitting at spectrometer
        dipole boundaries, so we allow ±1 drift difference.
        """
        from collections import Counter
        counts = {name: Counter(type(e).__name__ for e in bl)
                  for name, bl in beamlines.items()}

        ref_name, ref_counts = next(iter(counts.items()))
        for name, c in counts.items():
            all_types = set(ref_counts) | set(c)
            for t in all_types:
                diff = abs(ref_counts.get(t, 0) - c.get(t, 0))
                if t == "driftLattice":
                    assert diff <= 1, (
                        f"Drift count mismatch: {ref_name}={ref_counts[t]} vs {name}={c[t]}"
                    )
                else:
                    assert diff == 0, (
                        f"{t} count mismatch: {ref_name}={ref_counts.get(t, 0)} vs {name}={c.get(t, 0)}"
                    )


# ── Energy dependence ────────────────────────────────────────────────────

class TestEnergyDependence:
    @pytest.fixture
    def bl(self):
        candidates = [JSON_PATH, EXCEL_PATH]
        if _HAS_TRACKED_DICT:
            candidates.insert(0, YAML_PATH)
        path = next((p for p in candidates if p.exists()), None)
        if path is None:
            pytest.skip("No lattice file found")
        return load_beamline(path)

    def test_matrix_changes_with_energy(self, bl):
        """Transfer matrices should change when energy is set."""
        set_energy(bl, 45.0)
        M_45 = cumulative_matrix(bl)
        set_energy(bl, 30.0)
        M_30 = cumulative_matrix(bl)
        assert not np.allclose(M_45, M_30), "Matrices unchanged with energy change"
        # Reset to default
        set_energy(bl, 45.0)

    def test_cross_format_at_different_energies(self):
        """Formats should agree at non-default energies too."""
        bls = {}
        _paths = [EXCEL_PATH, JSON_PATH]
        if _HAS_TRACKED_DICT:
            _paths.append(YAML_PATH)
        for path in _paths:
            if path.exists():
                bls[path.stem] = load_beamline(path)
        if len(bls) < 2:
            pytest.skip("Need at least 2 formats")

        for KE in [20.0, 40.0, 100.0]:
            matrices = {}
            for name, bl in bls.items():
                set_energy(bl, KE)
                matrices[name] = cumulative_matrix(bl)
                set_energy(bl, 45.0)  # reset

            ref_name, ref_M = next(iter(matrices.items()))
            for name, M in matrices.items():
                diff = np.max(np.abs(M - ref_M))
                assert diff < 1e-6, (
                    f"KE={KE}: matrix diff ({ref_name} vs {name}): {diff:.2e}"
                )


# ── Regression tests with frozen reference values ────────────────────────

class TestRegression:
    """Frozen reference values from validated studies.
    These catch silent changes that alter physics results."""

    @pytest.fixture
    def bl_40MeV(self):
        """UH FEL beamline at 40 MeV (the operational energy)."""
        candidates = [JSON_PATH, EXCEL_PATH]
        if _HAS_TRACKED_DICT:
            candidates.insert(0, YAML_PATH)
        path = next((p for p in candidates if p.exists()), None)
        if path is None:
            pytest.skip("No lattice file found")
        bl = load_beamline(path)
        set_energy(bl, 40.0)
        return bl

    def test_element_count(self, bl_40MeV):
        """UH FEL beamline should have exactly 118 elements (or close, with drift splitting)."""
        # The exact count depends on format and drift splitting but should be ≥ 118
        assert len(bl_40MeV) >= 118

    def test_total_length(self, bl_40MeV):
        """Total beamline length should be ~14.76 m (full line including spectrometer)."""
        total_L = sum(e.length for e in bl_40MeV)
        assert 14.0 < total_L < 16.0, f"Total length = {total_L:.4f} m"

    def test_quad_count(self, bl_40MeV):
        """UH FEL has 23 quadrupoles in the transport line."""
        n_quads = sum(1 for e in bl_40MeV if isinstance(e, (qpfLattice, qpdLattice)))
        # Some formats may have different quad counts due to splitting, but ~23
        assert 20 <= n_quads <= 30, f"Quad count = {n_quads}"

    def test_dipole_count(self, bl_40MeV):
        """UH FEL has chicane dipoles + spectrometer dipoles."""
        n_dph = sum(1 for e in bl_40MeV if isinstance(e, dipole))
        n_dpw = sum(1 for e in bl_40MeV if isinstance(e, dipole_wedge))
        assert n_dph >= 4, f"DPH count = {n_dph}"  # at least 4 chicane
        assert n_dpw >= 8, f"DPW count = {n_dpw}"  # at least 8 wedges

    def test_cumulative_matrix_determinant(self, bl_40MeV):
        """Full beamline transfer matrix should have det ≈ 1."""
        M = cumulative_matrix(bl_40MeV)
        det = np.linalg.det(M)
        assert abs(det - 1.0) < 1e-6, f"det(M) = {det}"

    def test_on_axis_particle(self, bl_40MeV):
        """On-axis particle should remain near axis (stable beamline)."""
        particle = np.array([[0, 0, 0, 0, 0, 0]])
        result = propagate_particles(bl_40MeV, particle)
        # On-axis particle stays exactly on axis for linear optics
        np.testing.assert_allclose(result, 0.0, atol=1e-14)

    def test_off_energy_dispersion(self, bl_40MeV):
        """Off-energy particle (δ=1%) should be displaced at exit due to chicane dispersion."""
        on_axis = np.array([[0, 0, 0, 0, 0, 0]])
        off_energy = np.array([[0, 0, 0, 0, 0, 0.01]])
        r_on = propagate_particles(bl_40MeV, on_axis)
        r_off = propagate_particles(bl_40MeV, off_energy)
        # Should have some dispersion from the chicane
        x_disp = abs(r_off[0, 0] - r_on[0, 0])
        assert x_disp > 1e-6, f"No dispersion detected: Δx = {x_disp}"

    def test_gradient_constant(self, bl_40MeV):
        """All quads should use the standard gradient G = 2.694 T/A/m."""
        for elem in bl_40MeV:
            if isinstance(elem, (qpfLattice, qpdLattice)):
                assert abs(elem.G - 2.694) < 1e-10, (
                    f"Quad {elem.name}: G = {elem.G}"
                )

    def test_rf_frequency(self, bl_40MeV):
        """All elements should use f = 2856 MHz."""
        for elem in bl_40MeV:
            assert abs(elem.f - 2856e6) < 1, f"Element {elem.name}: f = {elem.f}"


# ── Adapter-level tests ──────────────────────────────────────────────────

class TestFELsimAdapter:
    """Test the FELsimAdapter interface for consistency across formats."""

    @pytest.fixture(params=[
        pytest.param("excel", id="Excel"),
        pytest.param("json", id="JSON"),
        pytest.param("yaml", id="YAML", marks=_YAML_MARK),
    ])
    def adapter(self, request):
        from felsimAdapter import FELsimAdapter
        paths = {"excel": EXCEL_PATH, "json": JSON_PATH, "yaml": YAML_PATH}
        path = paths[request.param]
        if not path.exists():
            pytest.skip(f"File not found: {path}")
        if request.param == "excel":
            return FELsimAdapter(excel_path=str(path))
        else:
            return FELsimAdapter(lattice_path=str(path))

    def test_simulate_returns_success(self, adapter):
        particles = adapter.generate_particles(num_particles=50)
        result = adapter.simulate(particles)
        assert result.success

    def test_simulate_preserves_particle_count(self, adapter):
        particles = adapter.generate_particles(num_particles=50)
        result = adapter.simulate(particles)
        assert result.final_particles.shape == particles.shape

    def test_simulate_particles_finite(self, adapter):
        particles = adapter.generate_particles(num_particles=50)
        result = adapter.simulate(particles)
        assert np.all(np.isfinite(result.final_particles))

    def test_generate_particles_shape(self, adapter):
        for n in [10, 100, 500]:
            p = adapter.generate_particles(num_particles=n)
            assert p.shape == (n, 6)

    def test_simulation_results_consistent(self):
        """All formats should give the same simulation result."""
        from felsimAdapter import FELsimAdapter

        adapters = {}
        _fmt_paths = [("Excel", EXCEL_PATH), ("JSON", JSON_PATH)]
        if _HAS_TRACKED_DICT:
            _fmt_paths.append(("YAML", YAML_PATH))
        for name, path in _fmt_paths:
            if path.exists():
                if name == "Excel":
                    adapters[name] = FELsimAdapter(excel_path=str(path))
                else:
                    adapters[name] = FELsimAdapter(lattice_path=str(path))

        if len(adapters) < 2:
            pytest.skip("Need at least 2 formats")

        # Use same seed for reproducible particles
        np.random.seed(42)
        particles = np.random.randn(20, 6) * [1e-3, 1e-4, 1e-3, 1e-4, 1e-3, 0.005]

        results = {}
        for name, adapter in adapters.items():
            results[name] = adapter.simulate(particles)

        ref_name, ref_result = next(iter(results.items()))
        for name, result in results.items():
            diff = np.max(np.abs(result.final_particles - ref_result.final_particles))
            assert diff < 1e-10, (
                f"Simulation diff ({ref_name} vs {name}): {diff:.2e}"
            )
