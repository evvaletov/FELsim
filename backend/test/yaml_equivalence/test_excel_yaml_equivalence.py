"""
Equivalence tests: verify that YAML lattice loading produces the same
beamline as Excel loading.

Mirrors the JSON equivalence tests in test_excel_json_equivalence.py.

Run directly:   python test_excel_yaml_equivalence.py
Run via pytest:  pytest test_excel_yaml_equivalence.py -v

Author: Eremey Valetov
"""

import sys
import os
import tempfile
import numpy as np
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

EXCEL_PATH = Path(__file__).resolve().parent.parent.parent.parent / "beam_excel" / "Beamline_elements.xlsx"
YAML_PATH = Path(__file__).resolve().parent.parent.parent.parent / "var" / "UH_FEL_beamline.yaml"

from excelElements import ExcelElements
from beamlineBuilder import BeamlineBuilder
from yamlLatticeLoader import YamlLatticeLoader
from excelToYaml import convert
from beamline import driftLattice, qpfLattice, qpdLattice, dipole, dipole_wedge

TOL = 1e-10
MATRIX_TOL = 1e-12


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def merge_consecutive_drifts(beamline):
    """Merge adjacent driftLattice elements into single drifts."""
    merged = []
    for elem in beamline:
        if isinstance(elem, driftLattice) and merged and isinstance(merged[-1], driftLattice):
            merged[-1] = driftLattice(merged[-1].length + elem.length)
        else:
            merged.append(elem)
    return merged


def merge_consecutive_drift_dicts(beamline):
    """Merge adjacent DRIFT dicts into single drifts."""
    merged = []
    for elem in beamline:
        if elem["type"] == "DRIFT" and merged and merged[-1]["type"] == "DRIFT":
            merged[-1] = {"type": "DRIFT", "length": merged[-1]["length"] + elem["length"]}
        else:
            merged.append(elem)
    return merged


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_yaml_schema_validation():
    """Verify that the YAML beamline passes schema validation."""
    loader = YamlLatticeLoader(str(YAML_PATH), validate_schema=True)
    assert len(loader.create_beamline()) > 0


def test_yaml_create_beamline():
    """Compare beamline.py class instances from Excel and YAML loading."""
    ee = ExcelElements(str(EXCEL_PATH))
    excel_bl = merge_consecutive_drifts(ee.create_beamline())

    loader = YamlLatticeLoader(str(YAML_PATH), validate_schema=False)
    yaml_bl = merge_consecutive_drifts(loader.create_beamline())

    assert len(excel_bl) == len(yaml_bl), \
        f"Element count: Excel={len(excel_bl)} vs YAML={len(yaml_bl)}"

    for i, (ex, ym) in enumerate(zip(excel_bl, yaml_bl)):
        ex_type = type(ex).__name__
        ym_type = type(ym).__name__
        assert ex_type == ym_type, f"[{i}] type: {ex_type} vs {ym_type}"
        assert abs(ex.length - ym.length) < TOL, \
            f"[{i}] {ex_type} length: {ex.length} vs {ym.length}"

        if isinstance(ex, (qpfLattice, qpdLattice)):
            assert abs(ex.current - ym.current) < TOL, \
                f"[{i}] {ex_type} current: {ex.current} vs {ym.current}"
        elif isinstance(ex, dipole):
            assert abs(ex.angle - ym.angle) < TOL, \
                f"[{i}] {ex_type} angle: {ex.angle} vs {ym.angle}"
        elif isinstance(ex, dipole_wedge):
            for attr in ("angle", "dipole_length", "dipole_angle", "pole_gap"):
                assert abs(getattr(ex, attr) - getattr(ym, attr)) < TOL, \
                    f"[{i}] {ex_type} {attr}: {getattr(ex, attr)} vs {getattr(ym, attr)}"


def test_yaml_parse_beamline():
    """Compare BeamlineBuilder-style dicts from Excel and YAML loading."""
    bb = BeamlineBuilder(str(EXCEL_PATH))
    excel_bl = merge_consecutive_drift_dicts(bb.parse_beamline())

    loader = YamlLatticeLoader(str(YAML_PATH), validate_schema=False)
    yaml_bl_raw = loader.parse_beamline()
    yaml_bl = merge_consecutive_drift_dicts(
        [e for e in yaml_bl_raw if e["length"] > 0 or e["type"] == "DRIFT"]
    )

    assert len(excel_bl) == len(yaml_bl), \
        f"Element count: Excel={len(excel_bl)} vs YAML={len(yaml_bl)}"

    for i, (ex, ym) in enumerate(zip(excel_bl, yaml_bl)):
        assert ex["type"] == ym["type"], \
            f"[{i}] type: {ex['type']} vs {ym['type']}"
        assert abs(ex["length"] - ym["length"]) < TOL, \
            f"[{i}] {ex['type']} length: {ex['length']} vs {ym['length']}"

        if ex["type"] in ("QPF", "QPD"):
            assert abs(ex["current"] - ym["current"]) < TOL, \
                f"[{i}] {ex['type']} current"
        elif ex["type"] == "DPH":
            assert abs(ex["angle"] - ym["angle"]) < TOL, \
                f"[{i}] {ex['type']} angle"
        elif ex["type"] == "DPW":
            for key in ("angle", "wedge_angle", "gap_wedge", "pole_gap"):
                assert abs(ex[key] - ym[key]) < TOL, \
                    f"[{i}] {ex['type']} {key}: {ex[key]} vs {ym[key]}"


def test_yaml_transfer_matrices():
    """Compare numeric transfer matrices from Excel and YAML beamlines."""
    ee = ExcelElements(str(EXCEL_PATH))
    excel_bl = merge_consecutive_drifts(ee.create_beamline())

    loader = YamlLatticeLoader(str(YAML_PATH), validate_schema=False)
    yaml_bl = merge_consecutive_drifts(loader.create_beamline())

    assert len(excel_bl) == len(yaml_bl)

    for i, (ex, ym) in enumerate(zip(excel_bl, yaml_bl)):
        ex_mat = ex._compute_numeric_matrix()
        ym_mat = ym._compute_numeric_matrix()
        diff = np.max(np.abs(ex_mat - ym_mat))
        assert diff < MATRIX_TOL, \
            f"[{i}] {type(ex).__name__}: max matrix diff = {diff:.2e}"


def test_yaml_full_propagation():
    """Propagate a test particle through both beamlines and compare output."""
    ee = ExcelElements(str(EXCEL_PATH))
    excel_bl = merge_consecutive_drifts(ee.create_beamline())

    loader = YamlLatticeLoader(str(YAML_PATH), validate_schema=False)
    yaml_bl = merge_consecutive_drifts(loader.create_beamline())

    particle = [[1.0, 0.1, 0.5, 0.05, 0.2, 0.01]]

    excel_out = list(particle)
    yaml_out = list(particle)
    for elem in excel_bl:
        excel_out = elem.useMatrice(excel_out)
    for elem in yaml_bl:
        yaml_out = elem.useMatrice(yaml_out)

    diff = np.max(np.abs(np.array(excel_out[0]) - np.array(yaml_out[0])))
    assert diff < 1e-10, f"Propagation max diff = {diff:.2e}"


def test_yaml_round_trip():
    """Convert Excel→YAML→beamline and verify it matches Excel→beamline."""
    with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False, mode="w") as f:
        tmp_path = f.name

    try:
        convert(str(EXCEL_PATH), tmp_path, name="round_trip_test")
        loader = YamlLatticeLoader(tmp_path, validate_schema=False)
        yaml_bl = merge_consecutive_drifts(loader.create_beamline())

        ee = ExcelElements(str(EXCEL_PATH))
        excel_bl = merge_consecutive_drifts(ee.create_beamline())

        assert len(excel_bl) == len(yaml_bl), \
            f"Element count: {len(excel_bl)} vs {len(yaml_bl)}"

        for i, (ex, ym) in enumerate(zip(excel_bl, yaml_bl)):
            assert type(ex).__name__ == type(ym).__name__, \
                f"[{i}] type mismatch"
            assert abs(ex.length - ym.length) < TOL, \
                f"[{i}] length mismatch"
    finally:
        os.unlink(tmp_path)


def test_yaml_kind_only_inline():
    """Inline v2 YAML with kind-only elements loads correctly."""
    import yaml

    lattice = {
        "beamline": {
            "metadata": {
                "format_version": 2,
                "name": "test",
                "version": "1.0",
                "reference_energy_mev": 45.0,
                "particle_type": "electron",
            },
            "beam_parameters": {
                "particle": {
                    "type": "electron",
                    "kinetic_energy_mev": 45.0,
                    "mass_mev": 0.51099895,
                    "charge_e": -1,
                },
                "rf_frequency_hz": 2.856e9,
            },
            "elements": [
                {
                    "name": "D1",
                    "kind": "Drift",
                    "s_start_m": 0.0,
                    "s_end_m": 0.5,
                    "length_m": 0.5,
                    "parameters": {},
                },
                {
                    "name": "Q1",
                    "kind": "Quadrupole",
                    "polarity": "focusing",
                    "s_start_m": 0.5,
                    "s_end_m": 0.6,
                    "length_m": 0.1,
                    "parameters": {"current_a": 2.0},
                },
                {
                    "name": "B1",
                    "kind": "SBend",
                    "s_start_m": 0.6,
                    "s_end_m": 0.7,
                    "length_m": 0.1,
                    "parameters": {
                        "bending_angle_deg": 15.0,
                        "dipole_length_m": 0.1,
                    },
                },
                {
                    "name": "K1",
                    "kind": "Kicker",
                    "plane": "horizontal",
                    "s_start_m": 0.7,
                    "s_end_m": 0.7,
                    "length_m": 0.0,
                    "parameters": {},
                },
                {
                    "name": "BPM1",
                    "kind": "Instrument",
                    "instrument_type": "BPM",
                    "s_start_m": 0.7,
                    "s_end_m": 0.7,
                    "length_m": 0.0,
                    "parameters": {},
                },
                {
                    "name": "M1",
                    "kind": "Marker",
                    "s_start_m": 0.7,
                    "s_end_m": 0.7,
                    "length_m": 0.0,
                    "parameters": {},
                },
            ],
        }
    }

    with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False, mode="w") as f:
        yaml.dump(lattice, f, sort_keys=False, default_flow_style=False)
        tmp_path = f.name

    try:
        loader = YamlLatticeLoader(tmp_path, validate_schema=False)
        assert loader._format_version == 2

        bl = loader.create_beamline()
        assert len(bl) == 3
        assert type(bl[0]).__name__ == "driftLattice"
        assert type(bl[1]).__name__ == "qpfLattice"
        assert type(bl[2]).__name__ == "dipole"

        parsed = loader.parse_beamline()
        type_list = [d["type"] for d in parsed]
        assert "DRIFT" in type_list
        assert "QPF" in type_list
        assert "DPH" in type_list
        assert "STH" in type_list
        assert "BPM" in type_list
    finally:
        os.unlink(tmp_path)


# ---------------------------------------------------------------------------
# Standalone mode
# ---------------------------------------------------------------------------

def main():
    print(f"Excel: {EXCEL_PATH}")
    print(f"YAML:  {YAML_PATH}")

    if not EXCEL_PATH.exists():
        print(f"ERROR: Excel file not found: {EXCEL_PATH}")
        sys.exit(1)
    if not YAML_PATH.exists():
        print(f"ERROR: YAML file not found: {YAML_PATH}")
        sys.exit(1)

    tests = [
        ("yaml_schema_validation", test_yaml_schema_validation),
        ("yaml_create_beamline", test_yaml_create_beamline),
        ("yaml_parse_beamline", test_yaml_parse_beamline),
        ("yaml_transfer_matrices", test_yaml_transfer_matrices),
        ("yaml_full_propagation", test_yaml_full_propagation),
        ("yaml_round_trip", test_yaml_round_trip),
        ("yaml_kind_only_inline", test_yaml_kind_only_inline),
    ]

    results = {}
    for name, func in tests:
        print(f"\n--- {name} ---")
        try:
            func()
            results[name] = True
            print(f"  PASSED")
        except Exception as e:
            results[name] = False
            print(f"  FAILED: {e}")

    print("\n" + "=" * 50)
    all_passed = all(results.values())
    for name, passed in results.items():
        print(f"  {name:25s} {'PASSED' if passed else 'FAILED'}")
    print()
    print("All tests passed." if all_passed else "Some tests FAILED.")
    if not all_passed:
        sys.exit(1)


if __name__ == "__main__":
    main()
