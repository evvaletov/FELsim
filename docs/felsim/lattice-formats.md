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
| JSON Schema | `var/lattice_schema_v3.json` | v3 | Extends v2 with Bn1/BendP |

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

### Version 3

Extends v2 with optional fields from the
[official PALS specification](https://pals-project.readthedocs.io/):

- `MagneticMultipoleP.Bn1` — quadrupole pole-tip field (Tesla), alternative to `current_a`
- `BendP.g_ref` — bend strength (1/m), alternative to `bending_angle_deg`
- `BendP.e1`/`e2` — edge angles (radians), alternative to `entrance_edge_angle_deg`/`exit_edge_angle_deg`

If both old and new fields are present, the new PALS-aligned fields take
precedence. Full specification: `manuals/lattice_specification_v3.md`.

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

## COSY INFINITY Conversion

[pals2cosy](https://github.com/evvaletov/pals2cosy) is a standalone tool that
converts both official PALS and FELsim v2/v3 lattice files to COSY INFINITY
FOX code. It auto-detects the format and supports `MagneticMultipoleP.Bn1`,
`BendP`, and `DIPOLE_WEDGE` triplet consolidation.

## PALS Use Case

FELsim is registered as a use case for the
[Particle Accelerator Lattice Standard](https://pals-project.readthedocs.io/)
(PALS).  The v2/v3 lattice format independently adopted PALS-compatible type
names, and the v3 extension adds optional PALS parameter groups
(`MagneticMultipoleP.Bn1`, `BendP`).  Submission materials are in
`reports/pals_submission/`.

## Specification

Full format specification with field descriptions and examples:

- v2: `manuals/lattice_specification.md`
- v3: `manuals/lattice_specification_v3.md`
