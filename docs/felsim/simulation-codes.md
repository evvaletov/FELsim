# Simulation Codes

FELsim uses three independent simulation codes for beam transport. All share
a unified lattice loading interface via `latticeLoader.py`.

## First-Order (FELsim Native)

**Adapter:** `felsimAdapter.py`

The built-in transport model uses 6×6 transfer matrices (`beamline.py`).
Each element has an analytic matrix: drifts, thin/thick quadrupoles, sector
dipoles, and wedge dipoles with Enge fringe fields.

Advantages:
- Fast (~milliseconds per beamline evaluation)
- Suitable for optimization (thousands of evaluations per study)
- Direct access to Twiss parameters at every element

Limitations:
- Linear optics only — no higher-order aberrations
- No space charge
- Fringe fields modelled analytically (Enge function), not from measured data

**Usage:**

```python
from felsimAdapter import FELsimAdapter

adapter = FELsimAdapter(lattice_path="var/UH_FEL_beamline.yaml")
adapter.create_beamline()
twiss = adapter.calculate_twiss(quad_currents)
```

## COSY INFINITY

**Adapter:** `cosyAdapter.py`

COSY INFINITY uses differential algebra (DA) to compute transfer maps to
arbitrary order. The Python adapter generates FOX code, executes COSY, and
parses the output.

Advantages:
- Higher-order maps capture nonlinear effects (chromaticity, geometric aberrations)
- Supports measured fringe fields via MGE (Mid-plane Generalized Enge) elements
- DA framework provides exact derivatives

Limitations:
- Requires COSY INFINITY binary
- Slower than first-order (~seconds per evaluation)
- FOX code generation adds complexity

The adapter writes a complete FOX program, runs the COSY binary, and reads
back transfer map coefficients and beam envelopes from the output files.

## RF-Track

**Adapter:** `rftrackAdapter.py`

RF-Track (A. Latina, CERN) performs full 6D particle tracking with optional
space charge. The adapter uses COMSOL field maps for dipoles and measured
gradient data for quadrupoles.

Advantages:
- Full particle tracking — captures all nonlinear effects
- Space charge modelling
- COMSOL field maps for realistic field distributions

Limitations:
- Slowest of the three codes (~minutes per evaluation)
- Requires RF-Track installation and field map files
- Not practical for optimization loops

**Usage:**

```python
from rftrackAdapter import RFTrackAdapter

adapter = RFTrackAdapter(lattice_path="var/UH_FEL_beamline.yaml")
adapter.create_beamline()
result = adapter.track(beam)
```

## When to Use Each Code

| Scenario | Recommended Code |
|----------|-----------------|
| Quadrupole optimization | First-order |
| Parameter sensitivity scans | First-order |
| Chromaticity / aberration checks | COSY INFINITY |
| Final validation with space charge | RF-Track |
| Cross-validation of linear results | COSY + RF-Track |

## Unified Interface

All three adapters accept the same `lattice_path=` parameter:

```python
adapter = FELsimAdapter(lattice_path="beam_excel/Beamline_elements.xlsx")
adapter = FELsimAdapter(lattice_path="var/UH_FEL_beamline.json")
adapter = FELsimAdapter(lattice_path="var/UH_FEL_beamline.yaml")
```

The `latticeLoader` module auto-detects the format and returns beamline
objects or dicts as appropriate for each adapter.
