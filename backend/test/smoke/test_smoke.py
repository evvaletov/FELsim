"""
FELsim smoke tests — fast verification of core code paths.

Run with:
    python run_smoke.py          # dynamic CLI runner
    python -m pytest test_smoke.py -v  # pytest mode

Author: Eremey Valetov
"""

import sys
import os
from pathlib import Path

# Prevent RF-Track C library from loading — it prints a startup banner on import
# and a license notice at interpreter exit. Smoke tests don't exercise RF-Track.
if 'RF_Track' not in sys.modules:
    sys.modules['RF_Track'] = None

import numpy as np
import pytest

# Ensure backend/ is importable
_BACKEND = Path(__file__).resolve().parent.parent.parent
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))

_PROJECT_ROOT = _BACKEND.parent
EXCEL_PATH = _PROJECT_ROOT / "beam_excel" / "Beamline_elements.xlsx"
JSON_PATH = _PROJECT_ROOT / "var" / "UH_FEL_beamline.json"
YAML_PATH = _PROJECT_ROOT / "var" / "UH_FEL_beamline.yaml"


# ---------------------------------------------------------------------------
# 1. Imports
# ---------------------------------------------------------------------------

def test_imports():
    """All core modules import without errors."""
    import beamline
    import ebeam
    import physicalConstants
    import excelElements
    import latticeLoaderBase
    import latticeLoader
    import jsonLatticeLoader
    import yamlLatticeLoader
    import excelToJson
    import excelToYaml
    import simulatorBase
    import felsimAdapter
    import simulatorFactory


# ---------------------------------------------------------------------------
# 2. Physical constants & relativistic calculations
# ---------------------------------------------------------------------------

def test_physical_constants():
    """PhysicalConstants values match CODATA 2018 and relativistic formulas are correct."""
    from physicalConstants import PhysicalConstants

    assert PhysicalConstants.C == 299_792_458
    assert abs(PhysicalConstants.E0_electron - 0.51099895) < 1e-6
    assert abs(PhysicalConstants.Q - 1.602176634e-19) < 1e-28

    gamma, beta = PhysicalConstants.relativistic_parameters(45.0, PhysicalConstants.E0_electron)
    assert gamma > 1.0
    assert 0 < beta < 1.0
    assert abs(gamma - (1 + 45.0 / PhysicalConstants.E0_electron)) < 1e-10

    # Rest energy from mass should match tabulated value
    computed = PhysicalConstants.compute_rest_energy(PhysicalConstants.M_e)
    assert abs(computed - PhysicalConstants.E0_electron) < 1e-5


# ---------------------------------------------------------------------------
# 3. Beam generation
# ---------------------------------------------------------------------------

def test_beam_generation():
    """ebeam generates a 6D Gaussian distribution of the correct shape."""
    from ebeam import beam

    eb = beam()
    std = np.array([1.0, 0.1, 1.0, 0.1, 1.0, 0.1])
    particles = eb.gen_6d_gaussian(0, std, num_particles=200)

    assert particles.shape == (200, 6)
    assert np.all(np.isfinite(particles))


# ---------------------------------------------------------------------------
# 4. Drift transfer matrix
# ---------------------------------------------------------------------------

def test_drift_matrix():
    """Drift transfer matrix is correct: unit diagonal, M12=M34=L, det=1."""
    from beamline import driftLattice

    L = 0.5
    d = driftLattice(L)
    M = d._compute_numeric_matrix()

    assert M.shape == (6, 6)
    assert abs(M[0, 1] - L) < 1e-14
    assert abs(M[2, 3] - L) < 1e-14
    assert abs(M[0, 0] - 1.0) < 1e-14
    assert abs(M[1, 1] - 1.0) < 1e-14
    # Determinant of 6×6 transfer matrix should be 1
    assert abs(np.linalg.det(M) - 1.0) < 1e-10


# ---------------------------------------------------------------------------
# 5. Quadrupole transfer matrices
# ---------------------------------------------------------------------------

def test_quadrupole_matrix():
    """QPF/QPD matrices are symplectic, and zero-current quad reduces to drift."""
    from beamline import qpfLattice, qpdLattice, driftLattice

    L = 0.0889

    # Zero-current quad should match drift
    qpf_zero = qpfLattice(current=0, length=L)
    M_qpf0 = qpf_zero._compute_numeric_matrix()
    M_drift = driftLattice(L)._compute_numeric_matrix()
    np.testing.assert_allclose(M_qpf0, M_drift, atol=1e-12)

    qpd_zero = qpdLattice(current=0, length=L)
    M_qpd0 = qpd_zero._compute_numeric_matrix()
    np.testing.assert_allclose(M_qpd0, M_drift, atol=1e-12)

    # Non-zero current: det should be 1
    for cls, I in [(qpfLattice, 2.0), (qpdLattice, 2.0)]:
        M = cls(current=I, length=L)._compute_numeric_matrix()
        assert abs(np.linalg.det(M) - 1.0) < 1e-10
        # Transverse 2×2 blocks must be symplectic: M11*M22 - M12*M21 = 1
        assert abs(M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0] - 1.0) < 1e-10
        assert abs(M[2, 2] * M[3, 3] - M[2, 3] * M[3, 2] - 1.0) < 1e-10


# ---------------------------------------------------------------------------
# 6. Dipole transfer matrices
# ---------------------------------------------------------------------------

def test_dipole_matrix():
    """Dipole and dipole_wedge matrices have det ≈ 1 and finite entries."""
    from beamline import dipole, dipole_wedge

    d = dipole(length=0.2, angle=15.0)
    M = d._compute_numeric_matrix()
    assert M.shape == (6, 6)
    assert np.all(np.isfinite(M))
    assert abs(np.linalg.det(M) - 1.0) < 1e-8

    dw = dipole_wedge(length=0.01, angle=7.5, dipole_length=0.2, dipole_angle=15.0)
    Mw = dw._compute_numeric_matrix()
    assert Mw.shape == (6, 6)
    assert np.all(np.isfinite(Mw))
    assert abs(np.linalg.det(Mw) - 1.0) < 1e-8


# ---------------------------------------------------------------------------
# 7. Excel loading
# ---------------------------------------------------------------------------

def test_excel_loading():
    """ExcelElements loads beamline from the standard Excel file."""
    if not EXCEL_PATH.exists():
        pytest.skip(f"Excel file not found: {EXCEL_PATH}")

    from excelElements import ExcelElements

    ee = ExcelElements(str(EXCEL_PATH))
    bl = ee.create_beamline()
    assert len(bl) > 0
    # Every element should have a positive length
    for elem in bl:
        assert elem.length > 0


# ---------------------------------------------------------------------------
# 8. JSON loading
# ---------------------------------------------------------------------------

def test_json_loading():
    """JsonLatticeLoader loads beamline from JSON."""
    if not JSON_PATH.exists():
        pytest.skip(f"JSON file not found: {JSON_PATH}")

    from jsonLatticeLoader import JsonLatticeLoader

    loader = JsonLatticeLoader(str(JSON_PATH))
    bl = loader.create_beamline()
    assert len(bl) > 0
    for elem in bl:
        assert elem.length > 0


# ---------------------------------------------------------------------------
# 9. YAML loading
# ---------------------------------------------------------------------------

def test_yaml_loading():
    """YamlLatticeLoader loads beamline from YAML."""
    if not YAML_PATH.exists():
        pytest.skip(f"YAML file not found: {YAML_PATH}")

    from yamlLatticeLoader import YamlLatticeLoader

    loader = YamlLatticeLoader(str(YAML_PATH))
    bl = loader.create_beamline()
    assert len(bl) > 0
    for elem in bl:
        assert elem.length > 0


# ---------------------------------------------------------------------------
# 10. Excel → JSON round-trip
# ---------------------------------------------------------------------------

def test_excel_to_json():
    """Excel → JSON round-trip: total length matches and conversion succeeds."""
    if not EXCEL_PATH.exists():
        pytest.skip(f"Excel file not found: {EXCEL_PATH}")

    from excelElements import ExcelElements
    from excelToJson import convert
    from jsonLatticeLoader import JsonLatticeLoader
    import tempfile

    ee = ExcelElements(str(EXCEL_PATH))
    bl_excel = ee.create_beamline()

    with tempfile.NamedTemporaryFile(suffix=".json", mode="w", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        convert(EXCEL_PATH, tmp_path)
        loader = JsonLatticeLoader(tmp_path)
        bl_json = loader.create_beamline()
    finally:
        os.unlink(tmp_path)

    assert len(bl_json) > 0
    # Total beamline length must match (drift splitting can change element count)
    len_excel = sum(e.length for e in bl_excel)
    len_json = sum(e.length for e in bl_json)
    assert abs(len_excel - len_json) < 1e-6, (
        f"Total length mismatch: Excel={len_excel:.6f}, JSON={len_json:.6f}"
    )


# ---------------------------------------------------------------------------
# 11. Beamline propagation
# ---------------------------------------------------------------------------

def test_beamline_propagation():
    """Particles through full beamline produce finite output."""
    if not EXCEL_PATH.exists():
        pytest.skip(f"Excel file not found: {EXCEL_PATH}")

    from excelElements import ExcelElements
    from ebeam import beam

    ee = ExcelElements(str(EXCEL_PATH))
    bl = ee.create_beamline()

    eb = beam()
    std = np.array([1.0, 0.1, 1.0, 0.1, 1.0, 0.1])
    particles = eb.gen_6d_gaussian(0, std, num_particles=50)

    current = particles.copy()
    for seg in bl:
        current = np.array(seg.useMatrice(current))

    assert current.shape == particles.shape
    assert np.all(np.isfinite(current))


# ---------------------------------------------------------------------------
# 12. FELsim adapter
# ---------------------------------------------------------------------------

def test_felsim_adapter():
    """FELsimAdapter instantiates and runs a simulation."""
    if not EXCEL_PATH.exists():
        pytest.skip(f"Excel file not found: {EXCEL_PATH}")

    from felsimAdapter import FELsimAdapter

    adapter = FELsimAdapter(excel_path=str(EXCEL_PATH))
    particles = adapter.generate_particles(num_particles=50)
    result = adapter.simulate(particles)

    assert result.success
    assert result.final_particles.shape == particles.shape
    assert np.all(np.isfinite(result.final_particles))


# ---------------------------------------------------------------------------
# 13. Adapter multi-format loading
# ---------------------------------------------------------------------------

def test_felsim_adapter_json():
    """FELsimAdapter loads from JSON lattice path."""
    if not JSON_PATH.exists():
        pytest.skip(f"JSON file not found: {JSON_PATH}")

    from felsimAdapter import FELsimAdapter

    adapter = FELsimAdapter(lattice_path=str(JSON_PATH))
    assert len(adapter._native_beamline) > 0
    particles = adapter.generate_particles(num_particles=50)
    result = adapter.simulate(particles)
    assert result.success


def test_felsim_adapter_yaml():
    """FELsimAdapter loads from YAML lattice path."""
    if not YAML_PATH.exists():
        pytest.skip(f"YAML file not found: {YAML_PATH}")

    from felsimAdapter import FELsimAdapter

    adapter = FELsimAdapter(lattice_path=str(YAML_PATH))
    assert len(adapter._native_beamline) > 0
    particles = adapter.generate_particles(num_particles=50)
    result = adapter.simulate(particles)
    assert result.success


def test_lattice_loader():
    """latticeLoader.create_beamline() works for all supported formats."""
    import latticeLoader

    for path in [EXCEL_PATH, JSON_PATH, YAML_PATH]:
        if not path.exists():
            continue
        bl = latticeLoader.create_beamline(str(path))
        assert len(bl) > 0
        for elem in bl:
            assert elem.length >= 0


def test_lattice_loader_parse():
    """latticeLoader.parse_beamline() works for all supported formats."""
    import latticeLoader

    for path in [EXCEL_PATH, JSON_PATH, YAML_PATH]:
        if not path.exists():
            continue
        parsed = latticeLoader.parse_beamline(str(path))
        assert len(parsed) > 0
        for elem in parsed:
            assert "type" in elem
            assert "length" in elem


# ---------------------------------------------------------------------------
# 15. Simulator factory
# ---------------------------------------------------------------------------

def test_simulator_factory():
    """SimulatorFactory lists available simulators and creates FELsim."""
    from simulatorFactory import SimulatorFactory

    available = SimulatorFactory.get_available_simulators()
    assert "felsim" in available

    sim = SimulatorFactory.create("felsim")
    assert sim.name == "Python"
