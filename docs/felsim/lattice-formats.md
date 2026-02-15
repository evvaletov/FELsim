# Lattice Formats

The UH FEL beamline can be described in three data formats. The Excel
spreadsheet is the primary source of truth; JSON and YAML provide
machine-readable alternatives with richer metadata.

## Format Summary

| Format | File | Version | Notes |
|--------|------|---------|-------|
| Excel | `beam_excel/Beamline_elements.xlsx` | — | Source of truth |
| JSON | `var/UH_FEL_beamline.json` | v1 | FELsim-native type names |
| JSON | `var/lattice_specification.json` | v2 | PALS-aligned, specification + examples |
| YAML | `var/UH_FEL_beamline.yaml` | v2 | PALS-aligned |

## Format Versions

### Version 1

Uses FELsim-native uppercase type names:

```json
{
  "type": "QUADRUPOLE",
  "polarity": "focusing",
  "length_m": 0.076,
  "parameters": {"current_a": 3.5}
}
```

### Version 2 (PALS-Aligned)

Uses [PALS](https://pals-project.readthedocs.io/) CamelCase type names
with `kind` as an alias for `type`:

```yaml
- kind: Quadrupole
  polarity: focusing
  length_m: 0.076
  parameters:
    current_a: 3.5
```

Key v2 conventions:

| PALS Name | Resolves To | Rule |
|-----------|-------------|------|
| `Quadrupole` | QPF / QPD | Via `polarity` field |
| `SBend` / `RBend` | DPH | Sector / rectangular bend |
| `Wiggler` | UND | Undulator / wiggler |
| `Kicker` | STV / STH | Via `plane` field |
| `Instrument` | BPM / OTR / SPC | Via `instrument_type` field |
| `Marker` | DRIFT | Zero-length drift |

## Loader Architecture

```
latticeLoader.py          ← unified entry point
├── jsonLatticeLoader.py  ← JSON file I/O
├── yamlLatticeLoader.py  ← YAML file I/O
└── latticeLoaderBase.py  ← format-independent parsing & type resolution
```

`latticeLoader.create_beamline(path)` returns a list of `beamline.py` class
instances. `latticeLoader.parse_beamline(path)` returns a list of dicts
compatible with `BeamlineBuilder` (used by the COSY adapter).

The loaders wrap data in `TrackedDict` objects to detect unused fields,
ensuring that the lattice file doesn't contain data that silently goes
unread.

## Excel Format

The Excel spreadsheet (`Beamline_elements.xlsx`) contains one row per
beamline element with columns for type code, position, length, current,
dipole parameters, and Enge coefficients. It is parsed by
`excelElements.py` → `ExcelElements` class, which produces `beamline.py`
objects directly.

## Converting Formats

`excelToYaml.py` converts an Excel beamline to YAML v2 format:

```bash
cd backend
python excelToYaml.py ../beam_excel/Beamline_elements.xlsx -o ../var/UH_FEL_beamline.yaml
```

## Specification

Full format specification with field descriptions and examples:
`manuals/lattice_specification.md`
