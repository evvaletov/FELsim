# FELsim JSON Lattice Specification

Version 1.0

## Overview

The JSON lattice format provides a human-readable, version-controllable alternative to the Excel-based beamline specification (`Beamline_elements.xlsx`). It contains all information needed to define a beamline: element sequence, physical parameters, beam properties, and simulator-specific settings.

The reference example is located at `var/lattice_specification.json`.


## Top-Level Structure

```json
{
  "beamline": {
    "metadata":           {},
    "beam_parameters":    {},
    "elements":           [],
    "lattice_structure":  {},
    "global_settings":    {},
    "simulator_specific": {}
  }
}
```

All data lives under the top-level `beamline` key.


## Versioning

The specification separates two concerns: the **format** evolving (structural changes to the JSON schema) and the **lattice data** changing (different magnet currents, added elements, etc.).

| Field            | Type    | Meaning                                                   |
|------------------|---------|-----------------------------------------------------------|
| `format_version` | integer | Schema version. Incremented when the JSON structure changes. The parser checks this value and rejects files it cannot handle. |
| `version`        | string  | User-managed label for the lattice data (e.g. `"2.3-chicane-opt"`). Purely informational; the parser does not interpret it. |

### Compatibility rules

- A parser that supports format_version *N* must accept any file with `format_version <= N`.
- Additions of new **optional** fields do not increment `format_version` (backward-compatible).
- Additions of new **required** fields, renaming of existing fields, or structural changes increment `format_version`.
- The JSON Schema file `var/lattice_schema_v1.json` defines the validation rules for format_version 1. Future versions get their own schema file (`lattice_schema_v2.json`, etc.).

### Validation

Structural validation uses JSON Schema (`var/lattice_schema_v1.json`). This catches missing required fields, wrong types, and invalid enum values at load time before any simulation code runs. To validate:

```python
import json, jsonschema

with open('var/lattice_specification.json') as f:
    data = json.load(f)
with open('var/lattice_schema_v1.json') as f:
    schema = json.load(f)

jsonschema.validate(data, schema)
```

The schema uses `additionalProperties: true` everywhere — extra fields are allowed, never rejected. This means the schema validates that required structure is present without being brittle to user extensions.

### Tracking unconsumed data

Structural validation catches malformed files; it does not catch data the parser silently ignores (misspelled parameter names, fields added in a newer lattice file than the parser expects, etc.). For this, the parser wraps the loaded JSON in a `TrackedDict` (`backend/trackedMapping.py`) — a dict-like object that records which keys are actually read. After parsing completes, calling `unaccessed()` returns every dotted path that was never touched:

```python
from trackedMapping import TrackedDict

data = TrackedDict(json.load(f))
# ... parser reads what it needs via normal dict access ...
unused = data.unaccessed()
# → ['beamline.elements[3].parameters.some_typo',
#    'beamline.simulator_specific.xsuite']
```

This is zero-maintenance: there is no list of "expected keys" to keep in sync with the schema. It catches both parser omissions (forgot to read a field) and user extras (added a field the parser doesn't handle). Unaccessed paths are logged as informational messages, not errors — extra data is tolerated, but the user is told about it.


## Sections

### `metadata`

| Field                  | Type    | Required | Description                               |
|------------------------|---------|----------|-------------------------------------------|
| `format_version`       | integer | yes      | Schema version (currently `1`)            |
| `name`                 | string  | yes      | Beamline identifier                       |
| `version`              | string  | yes      | Lattice data version (user-managed)       |
| `description`          | string  | no       | Human-readable description                |
| `author`               | string  | no       | Author or group                           |
| `date`                 | string  | no       | Date of last modification (ISO 8601)      |
| `reference_energy_mev` | float   | yes      | Reference kinetic energy in MeV           |
| `particle_type`        | string  | yes      | Particle species (`"electron"`, `"proton"`, or isotope string) |


### `beam_parameters`

Defines the reference particle and beam RF structure.

```json
"beam_parameters": {
  "particle": {
    "type": "electron",
    "kinetic_energy_mev": 45.0,
    "mass_mev": 0.51099895,
    "charge_e": -1
  },
  "rf_frequency_hz": 2.856e9
}
```

| Field              | Type   | Required | Description                                    |
|--------------------|--------|----------|------------------------------------------------|
| `particle.type`    | string | yes      | Particle species                               |
| `particle.kinetic_energy_mev` | float | yes | Kinetic energy in MeV                   |
| `particle.mass_mev`| float  | yes      | Rest mass in MeV/c²                           |
| `particle.charge_e`| int    | yes      | Charge in units of elementary charge           |
| `rf_frequency_hz`  | float  | yes      | RF frequency in Hz (used in M56 matrix terms)  |

The RF frequency is critical: it enters every transfer matrix calculation for the longitudinal plane (the M56 element).


### `elements`

Ordered array of beamline elements. Each element has:

| Field           | Type   | Required | Description                                |
|-----------------|--------|----------|--------------------------------------------|
| `name`          | string | yes      | Unique element identifier                  |
| `type`          | string | yes      | Element type (see table below)             |
| `s_start_m`     | float  | yes*     | Longitudinal start position in metres      |
| `s_end_m`       | float  | yes*     | Longitudinal end position in metres        |
| `length_m`      | float  | yes      | Physical length in metres (= s_end − s_start) |
| `parameters`    | object | yes      | Type-specific parameters (see below)       |
| `polarity`      | string | cond.    | `"focusing"` or `"defocusing"` (required for type `QUADRUPOLE`) |
| `aperture_m`    | float  | no       | Element aperture in metres                 |
| `fringe_fields` | object | no       | Fringe field specification                 |
| `optimization`  | object | no       | Optimization variable definition           |
| `metadata`      | object | no       | Descriptive metadata                       |

\* `s_start_m` and `s_end_m` may be `null` for placeholder/example elements not placed in the lattice.

**Drift insertion:** The parser automatically inserts drift spaces between elements when `s_start_m` of element N+1 exceeds `s_end_m` of element N. Explicit DRIFT elements are optional but recommended for clarity.


## Element Types

### Type names and aliases

Each element type has a canonical name and optional short aliases matching the internal code representation.

| Canonical Name  | Alias(es)      | Description                              |
|-----------------|----------------|------------------------------------------|
| `DRIFT`         | —              | Empty drift space                        |
| `QUADRUPOLE`    | `QPF`, `QPD`   | Quadrupole magnet (requires `polarity` unless alias used) |
| `DIPOLE`        | `DPH`          | Horizontal bending dipole                |
| `DIPOLE_WEDGE`  | `DPW`          | Dipole wedge (pole face rotation)        |
| `SOLENOID`      | `SOL`          | Solenoid magnet                          |
| `RF_CAVITY`     | `RFC`          | RF accelerating cavity                   |
| `SEXTUPOLE`     | `SXT`          | Sextupole magnet                         |
| `UNDULATOR`     | `UND`          | Undulator / wiggler                      |
| `BPM`           | —              | Beam position monitor                    |
| `OTR`           | —              | OTR screen                               |
| `CORRECTOR_V`   | `STV`          | Vertical corrector                       |
| `CORRECTOR_H`   | `STH`          | Horizontal corrector                     |
| `SPECTROMETER`  | `SPC`          | Spectrometer                             |

When using the short aliases `QPF` or `QPD`, the `polarity` field is not required (the polarity is implied by the type name). When using `QUADRUPOLE`, the `polarity` field is required.

**Zero-length elements** (`BPM`, `OTR`, `CORRECTOR_V`, `CORRECTOR_H`, `SPECTROMETER`) have `length_m = 0` and `s_start_m == s_end_m`. They are treated as pass-through markers: they do not affect beam transport but carry metadata (channel numbers, labels) used by the control system and diagnostics.


### Type-specific parameters

#### DRIFT

No parameters required. The `parameters` object should be empty (`{}`).

#### QUADRUPOLE

| Parameter    | Type  | Required | Description                                     |
|--------------|-------|----------|-------------------------------------------------|
| `current_a`  | float | yes      | Excitation current in Amperes                   |

The field gradient is computed as: `G × current`, where G is the gradient coefficient from `global_settings.quadrupole_gradient_coefficient_t_per_a_per_m` (default 2.694 T/A/m for UH FEL quadrupoles).

To override the gradient coefficient for a specific element, include `gradient_coefficient_t_per_a_per_m` in `parameters`.

#### DIPOLE

| Parameter          | Type  | Required | Description                                     |
|--------------------|-------|----------|-------------------------------------------------|
| `bending_angle_deg`| float | yes      | Bending angle in degrees (signed: positive = bend right) |
| `dipole_length_m`  | float | yes      | Effective magnetic length in metres              |
| `pole_gap_m`       | float | no       | Pole gap in metres (used for fringe field calculations) |

`dipole_length_m` is the effective magnetic length used in transfer matrix calculations. It typically equals `length_m` but is stated explicitly because the two can differ in principle.

#### DIPOLE_WEDGE

| Parameter          | Type  | Required | Description                                     |
|--------------------|-------|----------|-------------------------------------------------|
| `wedge_angle_deg`  | float | yes      | Wedge (pole face rotation) angle in degrees      |
| `dipole_angle_deg` | float | yes      | Bending angle of the associated dipole in degrees |
| `dipole_length_m`  | float | yes      | Magnetic length of the associated dipole in metres |
| `pole_gap_m`       | float | yes      | Pole gap in metres                               |

`length_m` for a DIPOLE_WEDGE corresponds to the wedge gap length (the `Gap wedge (m)` column in the Excel format). This is the thin-lens physical extent of the wedge element.

A dipole with wedge pole faces is represented as three consecutive elements: DIPOLE_WEDGE (entrance) → DIPOLE → DIPOLE_WEDGE (exit). The wedge elements carry the associated dipole's angle and length because they need these values for their transfer matrix calculation.

#### SOLENOID

| Parameter  | Type  | Required | Description                  |
|------------|-------|----------|------------------------------|
| `field_t`  | float | yes      | Axial magnetic field Bz in Tesla |

#### RF_CAVITY

| Parameter      | Type  | Required | Description                    |
|----------------|-------|----------|--------------------------------|
| `voltage_mv`   | float | yes      | Peak voltage in MV             |
| `frequency_hz` | float | yes      | RF frequency in Hz             |
| `phase_deg`    | float | yes      | RF phase in degrees            |

#### SEXTUPOLE

| Parameter  | Type  | Required | Description                       |
|------------|-------|----------|-----------------------------------|
| `strength` | float | yes      | Integrated sextupole strength     |

#### UNDULATOR

| Parameter      | Type  | Required | Description                       |
|----------------|-------|----------|-----------------------------------|
| `period_m`     | float | no       | Undulator period in metres        |
| `num_periods`  | int   | no       | Number of periods                 |
| `K_parameter`  | float | no       | Undulator K parameter             |
| `peak_field_t` | float | no       | Peak magnetic field in Tesla      |

The undulator is currently treated as a drift space for beam transport. The undulator-specific parameters are metadata for FEL calculations.


### Fringe fields

```json
"fringe_fields": {
  "enge_coefficients": [56.49, -50.79, 19.32, -3.621, 0.3315, -0.01193]
}
```

| Field               | Type       | Description                                        |
|---------------------|------------|----------------------------------------------------|
| `enge_coefficients` | float[] or null | Enge function coefficients for fringe field model |

Set to `null` or omit the `fringe_fields` object entirely when no Enge coefficients are available. The Enge function models the fringe field falloff at dipole edges:

$$B(s) = \frac{B_0}{1 + \exp\left(\sum_{i} c_i \left(\frac{s}{D}\right)^i\right)}$$

where $c_i$ are the Enge coefficients and $D$ is the pole gap.

In the UH FEL beamline, only the MkIII chicane (FC1/FC2) dipoles carry measured Enge coefficients. Other dipoles use the analytic triangle fringe field model built into the FELsim transfer matrix code.


### Optimization

```json
"optimization": {
  "variable": "I_QF1",
  "bounds": [0.0, 10.0]
}
```

| Field      | Type      | Description                                      |
|------------|-----------|--------------------------------------------------|
| `variable` | string    | Symbolic variable name for the optimizer          |
| `bounds`   | [min,max] | Allowed range for the variable                    |

This defines a named optimization variable associated with the element's primary parameter (current for quadrupoles, angle for dipoles). The optimizer can reference this variable by name.


### Element metadata

```json
"metadata": {
  "nomenclature": "DC1.QPF.021",
  "element_name": "Chromacity quad",
  "channel": 28,
  "label": "DPHQ",
  "sector": "DC1"
}
```

| Field          | Type   | Description                                           |
|----------------|--------|-------------------------------------------------------|
| `nomenclature` | string | Formal element designation (e.g. `"DC1.QPF.021"`)    |
| `element_name` | string | Descriptive name (e.g. `"Quad"`, `"OTR Screen"`)     |
| `channel`      | int    | Control system channel number                         |
| `label`        | string | Short label (e.g. `"DPHQ"`, `"OS1"`)                 |
| `sector`       | string | Beamline sector (e.g. `"LIN"`, `"DC1"`, `"FC1"`, `"FEL"`) |

All metadata fields are optional.


### `lattice_structure`

Groups elements into named sectors for organisational purposes.

```json
"lattice_structure": {
  "sectors": [
    {
      "name": "LIN",
      "description": "Linac section",
      "element_names": ["D_LIN_01", "LQ1", "BPM_LIN", "LQ2", "VC1"]
    }
  ]
}
```

Sectors reference elements by their `name` field. This grouping is purely organisational and does not affect simulation.


### `global_settings`

```json
"global_settings": {
  "coordinate_system": "felsim",
  "length_unit": "m",
  "angle_unit": "deg",
  "field_unit": "T",
  "current_unit": "A",
  "energy_unit": "MeV",
  "quadrupole_gradient_coefficient_t_per_a_per_m": 2.694
}
```

| Field                | Type   | Description                                           |
|----------------------|--------|-------------------------------------------------------|
| `coordinate_system`  | string | Default coordinate system: `felsim`, `cosy`, `rftrack`, `elegant` |
| `length_unit`        | string | Length unit (always `"m"`)                            |
| `angle_unit`         | string | Angle unit (always `"deg"`)                           |
| `field_unit`         | string | Magnetic field unit (always `"T"`)                    |
| `current_unit`       | string | Current unit (always `"A"`)                           |
| `energy_unit`        | string | Energy unit (always `"MeV"`)                          |
| `quadrupole_gradient_coefficient_t_per_a_per_m` | float | Default gradient calibration for quadrupoles (T/A/m) |

The gradient coefficient converts quadrupole excitation current to field gradient: `gradient [T/m] = G × I [A]`. The default value 2.694 T/A/m is the measured calibration for UH FEL quadrupoles.


### `simulator_specific`

Optional per-simulator configuration that does not belong in the generic lattice description.

```json
"simulator_specific": {
  "cosy": {
    "order": 3,
    "dimensions": 2,
    "use_enge_coefficients": true,
    "dipole_aperture_m": 0.0127,
    "quad_aperture_m": 0.027
  },
  "rftrack": {
    "default_aperture_m": 0.05,
    "space_charge": false
  },
  "felsim": {
    "fringe_delta_z_m": 0.01,
    "fringe_origin_factor": 0.99
  }
}
```


## Mapping to Internal Representations

When the JSON is parsed, element types and parameters are mapped to the internal dictionary format used by `BeamlineBuilder` and the simulator adapters.

### Type mapping

| JSON type       | JSON polarity    | Internal type | beamline.py class |
|-----------------|------------------|---------------|-------------------|
| `DRIFT`         | —                | `DRIFT`       | `driftLattice`    |
| `QUADRUPOLE`    | `"focusing"`     | `QPF`         | `qpfLattice`      |
| `QUADRUPOLE`    | `"defocusing"`   | `QPD`         | `qpdLattice`      |
| `QPF`           | —                | `QPF`         | `qpfLattice`      |
| `QPD`           | —                | `QPD`         | `qpdLattice`      |
| `DIPOLE`        | —                | `DPH`         | `dipole`          |
| `DPH`           | —                | `DPH`         | `dipole`          |
| `DIPOLE_WEDGE`  | —                | `DPW`         | `dipole_wedge`    |
| `DPW`           | —                | `DPW`         | `dipole_wedge`    |
| `SOLENOID`      | —                | `SOLENOID`    | —                 |
| `RF_CAVITY`     | —                | `RF_CAVITY`   | —                 |
| `SEXTUPOLE`     | —                | `SEXTUPOLE`   | —                 |
| `UNDULATOR`     | —                | `UND`         | `driftLattice`*   |

\* Undulators are treated as drift spaces for beam transport purposes.

### Parameter mapping to `BeamlineBuilder` dict

| JSON field                    | Internal dict key    | Applies to        |
|-------------------------------|----------------------|--------------------|
| `type` (mapped)               | `type`               | all                |
| `length_m`                    | `length`             | all                |
| `s_start_m`                   | `z_start`            | all                |
| `s_end_m`                     | `z_end`              | all                |
| `parameters.current_a`        | `current`            | QUADRUPOLE         |
| `parameters.bending_angle_deg`| `angle`              | DIPOLE             |
| `parameters.dipole_length_m`  | (used as `length` for dipole class) | DIPOLE |
| `parameters.wedge_angle_deg`  | `wedge_angle`        | DIPOLE_WEDGE       |
| `parameters.dipole_angle_deg` | `angle`              | DIPOLE_WEDGE       |
| `parameters.dipole_length_m`  | (dipole_length)      | DIPOLE_WEDGE       |
| `parameters.pole_gap_m`       | `pole_gap`           | DIPOLE, DIPOLE_WEDGE |
| `fringe_fields.enge_coefficients` | `enge_fct`       | DIPOLE, DIPOLE_WEDGE |

### Parameter mapping to Excel columns

For reference, here is how JSON fields correspond to the Excel columns in `Beamline_elements.xlsx`:

| JSON field                     | Excel column                          |
|--------------------------------|---------------------------------------|
| `s_start_m`                    | `z_start`                             |
| `s_end_m`                      | `z_end`                               |
| `parameters.current_a`         | `Current (A)`                         |
| `parameters.bending_angle_deg` | `Dipole Angle (deg)`                  |
| `parameters.dipole_length_m`   | `Dipole length (m)`                   |
| `parameters.wedge_angle_deg`   | `Dipole wedge (deg)`                  |
| `length_m` (for DIPOLE_WEDGE)  | `Gap wedge (m)`                       |
| `parameters.pole_gap_m`        | `Pole gap (m)`                        |
| `fringe_fields.enge_coefficients` | `Fringe Field Enge coefficients`   |
| `metadata.nomenclature`        | `Nomenclature`                        |
| `metadata.element_name`        | `Element name`                        |
| `metadata.channel`             | `Channel`                             |
| `metadata.label`               | `Label`                               |
| `metadata.sector`              | `Sector`                              |
| `type` (mapped)                | `Element`                             |


## Coordinate Systems

The `global_settings.coordinate_system` field specifies the default coordinate system. FELsim supports four systems, each with different conventions:

| System    | Coordinates                                                            |
|-----------|------------------------------------------------------------------------|
| `felsim`  | [x(mm), x'(mrad), y(mm), y'(mrad), ΔToF/T×10³, δW/W×10³]           |
| `cosy`    | [x(m), a, y(m), b, l(m), δK]                                         |
| `rftrack` | [x(mm), x'(mrad), y(mm), y'(mrad), t(mm/c), P(MeV/c)]              |
| `elegant` | [x(m), x'(rad), y(m), y'(rad), t(s), δ]                              |

The lattice file itself always uses SI units (metres, degrees, Tesla, Amperes). Coordinate system conversions are handled by the simulator adapters at runtime.


## Design Notes

**Why separate `length_m` and `dipole_length_m`?**
For a DIPOLE, `length_m` is the physical extent along the beamline (s_end − s_start). `dipole_length_m` is the effective magnetic length used in transfer matrix calculations. They are typically equal but can differ (e.g. if the magnetic field extends beyond the physical yoke). Both are stated explicitly to avoid ambiguity.

**Why does DIPOLE_WEDGE carry dipole parameters?**
A wedge element represents the pole-face rotation at the entrance or exit of a dipole. Its transfer matrix depends on the bending radius of the associated dipole, which is computed from `dipole_angle_deg` and `dipole_length_m`. Carrying these parameters directly on the wedge element avoids fragile cross-references between elements.

**Drift insertion:**
Drifts between elements are inferred from position gaps: if element N ends at s=1.0 and element N+1 starts at s=1.5, a 0.5 m drift is automatically inserted. Explicit DRIFT elements may be included for clarity or to attach metadata to specific drift sections. Including them is recommended for a complete, self-documenting lattice file.

**Zero-length elements:**
Diagnostic and corrector elements (BPM, OTR, STV, STH, SPC) have zero physical length. They appear in the element list with `length_m = 0` and `s_start_m == s_end_m`. They do not contribute to beam transport but are preserved for control system integration and lattice documentation.
