# FELsim Lattice Specification v3

## Overview

Version 3 extends the v2 lattice format with optional fields aligned to the
[official PALS specification](https://pals-project.readthedocs.io/) while
maintaining full backward compatibility with v1 and v2 files.

A v3 loader reads all v1/v2 files unchanged. New optional fields provide
alternative parameter representations that coexist with existing ones.

## Changes from v2

| Area | v2 | v3 | Migration |
|------|----|----|-----------|
| `format_version` | `2` | `3` | Bump; v1/v2 files still load |
| Quad strength | `current_a` only | `current_a` or `MagneticMultipoleP.Bn1` | If Bn1 present, takes precedence |
| Dipole edges | `entrance_edge_angle_deg`/`exit_edge_angle_deg` in parameters | Also `BendP.e1`/`BendP.e2` (radians) | If BendP present, converted to degrees |
| Dipole bend | `bending_angle_deg` | Also `BendP.g_ref` (1/m) | If g_ref present, angle = g_ref × L × 180/π |
| Angle unit | Degrees (implicit) | `global_settings.angle_unit`: `deg` (default) or `rad` | Old files unchanged |
| Positioning | Requires `s_start_m`/`s_end_m` | Optional; inferred from lengths if absent | Flat arrays still need positions |

## Quadrupole: `MagneticMultipoleP.Bn1`

PALS defines quadrupole strength via the normal magnetic multipole component
Bn1 — the pole-tip field in Tesla.

```yaml
- name: Q1
  kind: Quadrupole
  polarity: focusing
  s_start_m: 0.358775
  s_end_m: 0.447675
  length_m: 0.0889
  parameters:
    current_a: 0.885719           # FELsim v2: current in Amperes
  MagneticMultipoleP:
    Bn1: -0.03221                 # v3: pole-tip field in Tesla
```

**Precedence:** If `MagneticMultipoleP.Bn1` is present, it takes precedence
over `current_a`. The loader computes `current_a = Bn1 / (G × r)` for backward
compatibility with adapters that expect current.

**Sign convention (PALS):** Positive Bn1 = horizontally focusing for positive
charge. For electrons (negative charge), positive Bn1 = horizontally defocusing.

**Relation to current:** Bn1 = sign × G × I × r, where G is the gradient
coefficient (T/A/m) and r is the bore radius (m). Note that `aperture_m` is
the bore **diameter** (not radius); the loader divides by 2 internally.

## Dipole: `BendP`

PALS uses the `BendP` parameter group for bend geometry:

```yaml
- name: B1
  kind: SBend
  s_start_m: 1.0
  s_end_m: 1.089
  length_m: 0.0889
  parameters:
    bending_angle_deg: 1.5        # FELsim v2: angle in degrees
    dipole_length_m: 0.0889
    pole_gap_m: 0.014478
    entrance_edge_angle_deg: 0.0
    exit_edge_angle_deg: 1.5
  BendP:
    g_ref: 0.29476                # v3: 1/m → angle = g_ref × L
    e1: 0.0                       # v3: entrance edge (radians)
    e2: 0.02618                   # v3: exit edge (radians)
```

**Precedence:** If `BendP` is present, its fields take precedence over the
corresponding v2 `parameters` fields. The loader converts:

- `angle_deg = BendP.g_ref × length_m × 180/π`
- `entrance_edge_angle_deg = BendP.e1 × 180/π`
- `exit_edge_angle_deg = BendP.e2 × 180/π`

## Angle Unit (reserved for future use)

A new optional `angle_unit` in `global_settings` is reserved for future use:

```yaml
global_settings:
  angle_unit: rad  # or deg (default)
```

**Status:** Not yet implemented. If `angle_unit: rad` is present, the loader
emits a warning and continues to interpret all `*_deg` parameter values as
degrees. The `BendP` fields (`g_ref`, `e1`, `e2`) are always in radians
regardless of this setting.

## JSON Schema

A v3 JSON schema is provided at `var/lattice_schema_v3.json`. It extends the v2
schema with:

- `MagneticMultipoleP` object on Quadrupole elements (optional)
- `BendP` object on SBend/RBend elements (optional)
- `global_settings.angle_unit` enum (optional)
- `format_version: 3`

## Backward Compatibility

| File format_version | v3 loader behavior |
|---------------------|-------------------|
| 1 | Full support (unchanged) |
| 2 | Full support (unchanged) |
| 3 | New fields recognized |

No existing field is removed or reinterpreted. Files using only v2 fields work
identically in the v3 loader.

## Conversion Tools

[pals2cosy](https://github.com/evvaletov/pals2cosy) converts both official PALS
(`PALS:` root) and FELsim v2/v3 (`beamline:` root) lattice files to COSY
INFINITY FOX code. It supports `MagneticMultipoleP.Bn1` and `BendP` natively.
