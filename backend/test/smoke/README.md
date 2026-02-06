# FELsim Smoke Tests

Fast build-verification tests that exercise core FELsim code paths. Run them
after installation, dependency updates, or major refactors to catch breakage
early.

## Running

```bash
cd backend

# Dynamic CLI runner (spinner + coloured output)
python test/smoke/run_smoke.py

# Or via pytest
python -m pytest test/smoke/ -v
```

## What each test verifies

| Test | What it checks |
|------|----------------|
| `test_imports` | All core backend modules import without errors |
| `test_physical_constants` | CODATA 2018 values and relativistic parameter formulas |
| `test_beam_generation` | `ebeam.gen_6d_gaussian` produces correct shape and finite values |
| `test_drift_matrix` | Drift transfer matrix structure and unit determinant |
| `test_quadrupole_matrix` | QPF/QPD symplecticity; zero-current reduces to drift |
| `test_dipole_matrix` | Dipole and dipole_wedge matrices are finite with det ≈ 1 |
| `test_excel_loading` | `ExcelElements` loads the standard Excel beamline |
| `test_json_loading` | `JsonLatticeLoader` loads the JSON beamline |
| `test_yaml_loading` | `YamlLatticeLoader` loads the YAML beamline |
| `test_excel_to_json` | Excel → JSON round-trip preserves total beamline length |
| `test_beamline_propagation` | Particles through full beamline stay finite |
| `test_felsim_adapter` | `FELsimAdapter` instantiates and simulates |
| `test_felsim_adapter_json` | `FELsimAdapter` loads from JSON lattice path |
| `test_felsim_adapter_yaml` | `FELsimAdapter` loads from YAML lattice path |
| `test_lattice_loader` | `latticeLoader.create_beamline()` works for all formats |
| `test_lattice_loader_parse` | `latticeLoader.parse_beamline()` works for all formats |
| `test_simulator_factory` | `SimulatorFactory` lists and creates simulators |

Tests that require data files (`Beamline_elements.xlsx`, `UH_FEL_beamline.json`,
`UH_FEL_beamline.yaml`) or optional dependencies (RF-Track, COSY) are skipped
gracefully when those resources are unavailable.

## Adding a new test

1. Add a `def test_<name>():` function to `test_smoke.py`.
2. Use `pytest.skip(reason)` for optional-dependency guards.
3. Keep tests fast (< 2 s each) — these are smoke tests, not benchmarks.
