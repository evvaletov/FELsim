"""
Equivalence tests: verify that JSON lattice loading produces the same
beamline as Excel loading.

Tests compare:
  1. create_beamline() — element types, lengths, and constructor parameters
  2. parse_beamline()  — BeamlineBuilder-compatible dict fields
  3. Transfer matrices — numeric matrices for all active elements
  4. Full propagation  — end-to-end particle tracking
  5. Schema validation
  6. Round-trip conversion (Excel → JSON → beamline)

Run directly:   python test_excel_json_equivalence.py
Run via pytest:  pytest test_excel_json_equivalence.py -v

Author: Eremey Valetov
"""

import sys
import os
import tempfile
import numpy as np
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

EXCEL_PATH = Path(__file__).resolve().parent.parent.parent.parent / "beam_excel" / "Beamline_elements.xlsx"
JSON_PATH = Path(__file__).resolve().parent.parent.parent.parent / "var" / "UH_FEL_beamline.json"

from excelElements import ExcelElements
from beamlineBuilder import BeamlineBuilder
from jsonLatticeLoader import JsonLatticeLoader
from excelToJson import convert
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

def test_schema_validation():
    """Verify that the converted JSON passes schema validation."""
    loader = JsonLatticeLoader(str(JSON_PATH), validate_schema=True)
    assert len(loader.create_beamline()) > 0


def test_create_beamline():
    """Compare beamline.py class instances from Excel and JSON loading."""
    ee = ExcelElements(str(EXCEL_PATH))
    excel_bl = merge_consecutive_drifts(ee.create_beamline())

    loader = JsonLatticeLoader(str(JSON_PATH), validate_schema=False)
    json_bl = merge_consecutive_drifts(loader.create_beamline())

    assert len(excel_bl) == len(json_bl), \
        f"Element count: Excel={len(excel_bl)} vs JSON={len(json_bl)}"

    for i, (ex, js) in enumerate(zip(excel_bl, json_bl)):
        ex_type = type(ex).__name__
        js_type = type(js).__name__
        assert ex_type == js_type, f"[{i}] type: {ex_type} vs {js_type}"
        assert abs(ex.length - js.length) < TOL, \
            f"[{i}] {ex_type} length: {ex.length} vs {js.length}"

        if isinstance(ex, (qpfLattice, qpdLattice)):
            assert abs(ex.current - js.current) < TOL, \
                f"[{i}] {ex_type} current: {ex.current} vs {js.current}"

        elif isinstance(ex, dipole):
            assert abs(ex.angle - js.angle) < TOL, \
                f"[{i}] {ex_type} angle: {ex.angle} vs {js.angle}"

        elif isinstance(ex, dipole_wedge):
            for attr in ("angle", "dipole_length", "dipole_angle", "pole_gap"):
                assert abs(getattr(ex, attr) - getattr(js, attr)) < TOL, \
                    f"[{i}] {ex_type} {attr}: {getattr(ex, attr)} vs {getattr(js, attr)}"


def test_parse_beamline():
    """Compare BeamlineBuilder-style dicts from Excel and JSON loading."""
    bb = BeamlineBuilder(str(EXCEL_PATH))
    excel_bl = merge_consecutive_drift_dicts(bb.parse_beamline())

    loader = JsonLatticeLoader(str(JSON_PATH), validate_schema=False)
    json_bl_raw = loader.parse_beamline()
    json_bl = merge_consecutive_drift_dicts(
        [e for e in json_bl_raw if e["length"] > 0 or e["type"] == "DRIFT"]
    )

    assert len(excel_bl) == len(json_bl), \
        f"Element count: Excel={len(excel_bl)} vs JSON={len(json_bl)}"

    for i, (ex, js) in enumerate(zip(excel_bl, json_bl)):
        assert ex["type"] == js["type"], \
            f"[{i}] type: {ex['type']} vs {js['type']}"
        assert abs(ex["length"] - js["length"]) < TOL, \
            f"[{i}] {ex['type']} length: {ex['length']} vs {js['length']}"

        if ex["type"] in ("QPF", "QPD"):
            assert abs(ex["current"] - js["current"]) < TOL, \
                f"[{i}] {ex['type']} current"

        elif ex["type"] == "DPH":
            assert abs(ex["angle"] - js["angle"]) < TOL, \
                f"[{i}] {ex['type']} angle"

        elif ex["type"] == "DPW":
            for key in ("angle", "wedge_angle", "gap_wedge", "pole_gap"):
                assert abs(ex[key] - js[key]) < TOL, \
                    f"[{i}] {ex['type']} {key}: {ex[key]} vs {js[key]}"


def test_transfer_matrices():
    """Compare numeric transfer matrices from Excel and JSON beamlines."""
    ee = ExcelElements(str(EXCEL_PATH))
    excel_bl = merge_consecutive_drifts(ee.create_beamline())

    loader = JsonLatticeLoader(str(JSON_PATH), validate_schema=False)
    json_bl = merge_consecutive_drifts(loader.create_beamline())

    assert len(excel_bl) == len(json_bl)

    for i, (ex, js) in enumerate(zip(excel_bl, json_bl)):
        ex_mat = ex._compute_numeric_matrix()
        js_mat = js._compute_numeric_matrix()
        diff = np.max(np.abs(ex_mat - js_mat))
        assert diff < MATRIX_TOL, \
            f"[{i}] {type(ex).__name__}: max matrix diff = {diff:.2e}"


def test_full_propagation():
    """Propagate a test particle through both beamlines and compare output."""
    ee = ExcelElements(str(EXCEL_PATH))
    excel_bl = merge_consecutive_drifts(ee.create_beamline())

    loader = JsonLatticeLoader(str(JSON_PATH), validate_schema=False)
    json_bl = merge_consecutive_drifts(loader.create_beamline())

    particle = [[1.0, 0.1, 0.5, 0.05, 0.2, 0.01]]

    excel_out = list(particle)
    json_out = list(particle)
    for elem in excel_bl:
        excel_out = elem.useMatrice(excel_out)
    for elem in json_bl:
        json_out = elem.useMatrice(json_out)

    diff = np.max(np.abs(np.array(excel_out[0]) - np.array(json_out[0])))
    assert diff < 1e-10, f"Propagation max diff = {diff:.2e}"


def test_round_trip_conversion():
    """Convert Excel→JSON→beamline and verify it matches Excel→beamline."""
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False, mode="w") as f:
        tmp_path = f.name

    try:
        convert(str(EXCEL_PATH), tmp_path, name="round_trip_test")
        loader = JsonLatticeLoader(tmp_path, validate_schema=False)
        json_bl = merge_consecutive_drifts(loader.create_beamline())

        ee = ExcelElements(str(EXCEL_PATH))
        excel_bl = merge_consecutive_drifts(ee.create_beamline())

        assert len(excel_bl) == len(json_bl), \
            f"Element count: {len(excel_bl)} vs {len(json_bl)}"

        for i, (ex, js) in enumerate(zip(excel_bl, json_bl)):
            assert type(ex).__name__ == type(js).__name__, \
                f"[{i}] type mismatch"
            assert abs(ex.length - js.length) < TOL, \
                f"[{i}] length mismatch"
    finally:
        os.unlink(tmp_path)


# ---------------------------------------------------------------------------
# Standalone mode
# ---------------------------------------------------------------------------

def main():
    print(f"Excel: {EXCEL_PATH}")
    print(f"JSON:  {JSON_PATH}")

    if not EXCEL_PATH.exists():
        print(f"ERROR: Excel file not found: {EXCEL_PATH}")
        sys.exit(1)
    if not JSON_PATH.exists():
        print(f"ERROR: JSON file not found: {JSON_PATH}")
        sys.exit(1)

    tests = [
        ("schema_validation", test_schema_validation),
        ("create_beamline", test_create_beamline),
        ("parse_beamline", test_parse_beamline),
        ("transfer_matrices", test_transfer_matrices),
        ("full_propagation", test_full_propagation),
        ("round_trip", test_round_trip_conversion),
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
