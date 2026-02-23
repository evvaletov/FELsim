# FELsim Transport Line

FELsim models beam transport from the linac exit to the undulator entrance.
The primary goal is to optimize quadrupole currents so that the beam Twiss
parameters at the undulator match the targets required for FEL lasing.

## Architecture

FELsim integrates three simulation codes through a unified adapter interface.
All adapters accept a `lattice_path=` parameter that can point to an Excel,
JSON, or YAML beamline description.

| Code | Adapter | Method |
|------|---------|--------|
| FELsim (1st order) | `felsimAdapter.py` | 6×6 transfer matrices |
| COSY INFINITY | `cosyAdapter.py` | Differential algebra maps |
| RF-Track | `rftrackAdapter.py` | Full particle tracking |

The first-order code is used for rapid optimization; COSY and RF-Track
provide cross-validation with higher-fidelity physics.

## Beamline

The transport line consists of 118 elements over approximately 12.8 m:

- **23 quadrupoles** (QPF/QPD) — main matching elements
- **4 dipoles** — MkIII diagnostic chicane
- **Correctors** — vertical and horizontal steering
- **Diagnostics** — BPMs, OTR screens, spectrometer

Beamline data is maintained in three formats (see
[Lattice Formats](lattice-formats.md)).

## Key Files

| File | Purpose |
|------|---------|
| `backend/felsimAdapter.py` | First-order transfer matrix adapter |
| `backend/cosyAdapter.py` | COSY INFINITY adapter |
| `backend/rftrackAdapter.py` | RF-Track particle tracking adapter |
| `backend/beamline.py` | Transfer matrix element classes |
| `backend/beamOptimizer.py` | Nelder-Mead optimization wrapper |
| `backend/latticeLoader.py` | Unified lattice loading (Excel/JSON/YAML) |
| `backend/test/UHM_beamline_opt_*.py` | Optimization study scripts |

```{toctree}
:maxdepth: 1

simulation-codes
lattice-formats
optimization
beam-generation
fieldmap
```
