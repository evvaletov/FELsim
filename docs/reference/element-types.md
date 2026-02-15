# Beamline Element Types

Catalog of beamline element types supported by the lattice loaders. The
internal short names are used throughout FELsim; the v1 and v2 columns
show the names accepted in JSON/YAML lattice files.

## Active Elements

| Internal | v1 Name | v2 (PALS) Name | Key Parameters | Transfer Matrix |
|----------|---------|----------------|----------------|-----------------|
| QPF | `QUADRUPOLE` (polarity: focusing) | `Quadrupole` | `current_a` | Thick quad (focusing) |
| QPD | `QUADRUPOLE` (polarity: defocusing) | `Quadrupole` | `current_a` | Thick quad (defocusing) |
| DPH | `DIPOLE` / `DPH` | `SBend` / `RBend` | `bending_angle_deg`, `dipole_length_m`, `pole_gap_m` | Sector dipole |
| DPW | `DIPOLE_WEDGE` / `DPW` | — | `wedge_angle_deg`, `dipole_angle_deg`, `dipole_length_m`, `pole_gap_m`, Enge coefficients | Wedge dipole with fringe fields |
| SOL | `SOLENOID` / `SOL` | `Solenoid` | `current_a` | — |
| RFC | `RF_CAVITY` / `RFC` | `RFCavity` | — | — |
| SXT | `SEXTUPOLE` / `SXT` | `Sextupole` | — | — |
| UND | `UNDULATOR` / `UND` | `Wiggler` | — | Drift (in transport) |

## Correctors

| Internal | v1 Name | v2 (PALS) Name | Resolution |
|----------|---------|----------------|------------|
| STV | `CORRECTOR_V` / `STV` | `Kicker` (plane: vertical) | Via `plane` field |
| STH | `CORRECTOR_H` / `STH` | `Kicker` (plane: horizontal) | Via `plane` field |

## Diagnostics

| Internal | v1 Name | v2 (PALS) Name | Resolution |
|----------|---------|----------------|------------|
| BPM | `BPM` | `Instrument` (instrument_type: BPM) | Via `instrument_type` |
| OTR | `OTR` | `Instrument` (instrument_type: OTR) | Via `instrument_type` |
| SPC | `SPECTROMETER` / `SPC` | `Instrument` (instrument_type: SPECTROMETER) | Via `instrument_type` |

## Passive / Structural

| Internal | v1 Name | v2 (PALS) Name | Notes |
|----------|---------|----------------|-------|
| DRIFT | `DRIFT` | `Drift` | Free space |
| — | — | `Marker` | Zero-length, treated as drift |
| XRS | `XRS` | — | |
| BSW | `BSW` | — | |

## Type Resolution

The `_TYPE_ALIASES` mapping in `latticeLoaderBase.py` resolves all accepted
names to internal short names. Some v2 types require additional fields:

- **Quadrupole** → QPF or QPD based on `polarity` (focusing/defocusing)
- **Kicker** → STV or STH based on `plane` (vertical/horizontal)
- **Instrument** → BPM, OTR, or SPC based on `instrument_type`
- **Marker** → DRIFT (zero-length)

## Dipole Parameters

### DPH (Sector Dipole)

`dipole_length_m` is the **effective magnetic length**, which may differ
from the physical element length. The bending angle is specified in degrees.

### DPW (Wedge Dipole)

The wedge dipole represents the entrance/exit face of a dipole with a
non-normal beam angle. `gap_wedge` is the physical wedge length (= element
length). `dipole_length_m` is the associated sector dipole's magnetic
length. Enge coefficients specify the fringe field falloff; only the FC1
and FC2 chicane dipoles have measured values.

## Gradient Convention

The quadrupole gradient $G = 2.694$ T/(A·m) is hardcoded in `beamline.py`.
The RF-Track adapter allows configuring this value independently.
