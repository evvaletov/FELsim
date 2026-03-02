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
- Internal FIT optimizer for multi-stage Twiss matching

Limitations:
- Requires COSY INFINITY binary
- Slower than first-order (~seconds per evaluation)
- FOX code generation adds complexity
The adapter writes a complete FOX program, runs the COSY binary, and reads
back transfer map coefficients and beam envelopes from the output files.

**FIT optimization objectives:**

The COSY adapter supports both transverse and longitudinal FIT objectives:

| Axis | Parameter | COSY Expression | Notes |
|------|-----------|----------------|-------|
| x/y | alpha, beta, gamma | A0, B0, G0 | Twiss from transfer matrix + initial conditions |
| x/y | dispersion | ME(1,6), ME(3,6) | Chromatic dispersion |
| x/y | envelope | sqrt(epsilon * beta) | Requires `set_geometric_emittance()` |
| l | r56 | ME(5,6) | Path length vs energy; requires `dimensions=3` |
| l | r51, r52 | ME(5,1), ME(5,2) | Longitudinal coupling; requires `dimensions=3` |

Longitudinal objectives (R56 etc.) are needed for non-zero chirp or bunch
compression studies. For zero-chirp transverse matching, the longitudinal
phase space decouples and transverse-only objectives suffice.

**Aperture cuts:**

Element apertures are set via constructor parameters `quad_aperture` (default
0.027 m) and `dipole_aperture` (default 0.0127 m). These affect the
quadrupole `B_pole` computation for the MQ command. In particle tracking
mode, `enable_aperture_cuts()` generates COSY AP commands after each element
to kill rays exceeding the physical aperture:

- Quadrupoles: elliptic `AP r r 1` with `r = quad_aperture / 2`
- Dipoles: rectangular `AP w h 2` with `h = pole_gap / 2`,
  `w = dipole_half_width` (configurable, default 0.050 m)

The pipeline handles particle loss gracefully — zero survivors, partial loss,
and NaN coordinates all produce well-defined results with transmission logging.

## RF-Track

**Adapter:** `rftrackAdapter.py`

RF-Track (A. Latina, CERN) performs full 6D particle tracking with optional
space charge.

Advantages:
- Full particle tracking — captures all nonlinear effects
- Space charge modelling
- COMSOL field maps for realistic field distributions

Limitations:
- Slowest of the three codes (~minutes per evaluation)
- Requires RF-Track installation and field map files
- Not practical for optimization loops (but see prefix caching below)

**Usage:**

```python
from rftrackAdapter import RFTrackAdapter

adapter = RFTrackAdapter(lattice_path="var/UH_FEL_beamline.yaml")
adapter.create_beamline()
result = adapter.track(beam)
```

### SBend P/δ bug and analytical correction

RF-Track v2.5.5 has an upstream bug where `SBend` interprets Bunch6d's
absolute momentum (MeV/c) as momentum deviation δ = ΔP/P₀, producing
catastrophically wrong tracking output (~900 mm displacement for a 1 mm
off-axis particle).  See [SBend P/δ bug report](bugs/rftrack-sbend-p-delta.md).

**Workaround:** The adapter tracks dipole bodies as `Drift(L)` and applies an
analytical correction $M_\text{sector} \times M_\text{drift}^{-1}$ to the
$(x, x')$ coordinates, plus dispersion ($R_{16}$, $R_{26}$) and path-length
($R_{56}$) terms.  This is exact for a sector-bend transfer matrix.

### Segmented tracking architecture

`track_elements(beam, start, end)` tracks a sub-lattice with analytical dipole
corrections.  Internally, `_track_segmented()` groups consecutive non-dipole
elements into RF-Track `Lattice` objects for native tracking, and applies the
analytical correction at each DPH (dipole body) element.  DPW (wedge) elements
are treated as thin-lens edge kicks via `Quadrupole(L=1e-10, K1L)`.

This enables prefix caching: the static upstream portion of the beamline
(elements 0–86) is tracked once, and only the suffix (elements 87–117) is
re-tracked per optimizer evaluation.

### DPW edge kick workaround

FELsim models dipole wedges as thin-lens edge kicks with no physical drift
propagation.  The RF-Track adapter maps each DPW to a thin `Quadrupole`
with $K_{1L} = -|K_0| \tan\eta$, where $K_0 = |\theta|/L$ uses the unsigned
curvature.  See [edge kick sign bug report](bugs/rftrack-edge-kick-sign.md)
for why `abs()` is required.

### Known limitation

The thin-quad DPW model does not include the triangle-model fringe correction
$\phi$ in the y-plane edge kick formula $-\tan(\eta - \phi)/\rho$.  This
produces a small residual in $\beta_y$ compared to FELsim (which includes
$\phi$ analytically).

## When to Use Each Code

| Scenario | Recommended Code |
|----------|-----------------|
| Quadrupole optimization | First-order |
| Parameter sensitivity scans | First-order |
| Chromaticity / aberration checks | COSY INFINITY |
| R56 / bunch compression optimization | COSY INFINITY (dimensions=3) |
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
