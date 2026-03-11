# BUG: DPW Edge Kick Sign Inversion for Negative-Angle Dipoles

| Field | Value |
|-------|-------|
| **Date found** | 2026-03-01 |
| **Date fixed** | 2026-03-01 |
| **Severity** | Critical — wrong physics for all negative-angle dipoles |
| **Affected file** | `backend/beamline.py`, class `dipole_wedge` |
| **Affected methods** | `_compute_numeric_matrix()`, `_compute_symbolic_matrix()` |
| **Affects** | All FELsim simulations involving chicane (C1–C4) or negative-angle spectrometer dipoles |
| **Does NOT affect** | COSY adapter (correct by construction: CB + positive angle), RF-Track adapter |

## Symptom

FELsim transfer matrices disagreed with COSY INFINITY by orders of magnitude
for beamlines containing chicane dipoles. COSY particle tracking gave
$\sigma_z \approx 31$ ps for a 2 ps input beam (beam blowup due to
transverse instability), while FELsim's first-order optimizer reported a
stable solution.

Specifically:

| Code | $|\operatorname{Tr}_x/2|$ | $|\operatorname{Tr}_y/2|$ |
|------|---------------------------|---------------------------|
| **COSY INFINITY** | 4.17 (unstable) | 21.0 (unstable) |
| **FELsim (buggy)** | 0.14 (stable) | 2.69 (unstable) |
| **FELsim (fixed)** | 4.17 (unstable) | 23.2 (unstable) |

After the fix, FELsim's X-plane trace matches COSY exactly. The residual
Y-plane difference (21.0 vs 23.2) is due to the fringe-field model
(FELsim uses a triangle model; COSY uses hard-edge at `FR 0`).

## Root Cause

The `dipole_wedge` class computes the bending radius $R = \rho$ from the
associated dipole's signed angle:

```python
# BEFORE (buggy)
By = (M*c*beta*gamma / Q) * (dipole_angle * pi/180 / dipole_length)
R  = M*c*beta*gamma / (Q * By)
```

This simplifies to $R = L_\text{dipole} / (\theta_\text{dipole})$, where
$\theta$ carries the sign of the bending angle. For chicane dipoles with
$\theta = -4°$, $R < 0$.

The edge-kick matrix elements are:

$$
M_{21} = \frac{\tan\eta}{R}, \qquad M_{43} = -\frac{\tan(\eta - \phi)}{R}
$$

With $R < 0$, both kicks flip sign. This is physically wrong: the edge-kick
direction depends on the pole-face geometry (wedge angle $\eta$), not on the
bending direction. The standard beam optics convention (TRANSPORT, MAD-X,
COSY) uses $|\rho|$ in the edge-kick formula.

### Affected elements

| Dipole | Bending angle | Wedge angle | $R$ (buggy) | $T_x/R$ (buggy) | $T_x/R$ (correct) |
|--------|---------------|-------------|-------------|------------------|--------------------|
| C1–C4 (chicane) | $-4.0°$ | $2.018°$ | $-0.582$ m | $-0.0605$ | $+0.0605$ |
| SP3, SP4 (spectrometer) | $-11.25°$ | $11.25°$ | $-0.190$ m | $-1.045$ | $+1.045$ |
| FC1, FC2 (LINAC) | $+1.5°$ | various | $+3.396$ m | correct | correct |
| BC1, BC2 (bunch compressor) | $+5.0°$ | $2.536°$ | $+0.466$ m | correct | correct |

Positive-angle dipoles were unaffected because $R > 0$ already.

### Consequence for optimization

FELsim's Nelder-Mead optimizer found quad currents that stabilized the
**wrong** transfer matrix. When these currents were applied to COSY
(which correctly handles edge kicks via `CB` + positive angle), the
resulting beamline was unstable in both planes, causing massive beam blowup.

## Fix

Replace the signed-angle bending-radius computation with the absolute value
(lines 686–687 in `_compute_numeric_matrix`, lines 741–742 in
`_compute_symbolic_matrix`):

```python
# AFTER (fixed)
R = dipole_length / (abs(dipole_angle) * pi / 180)
```

This ensures $R = |\rho| > 0$ regardless of bending direction.

## Verification

### Element-by-element comparison (diagnostic script `diag_element_matrices.py`)

After the fix, the COSY analytical and FELsim element matrices agree:

- **Quadrupoles**: max $|\Delta M| \sim 10^{-4}$ (M56 drift term only)
- **Dipoles (transverse)**: max $|\Delta M| \sim 3 \times 10^{-4}$ (fringe model difference)
- **Dipoles (M50/M51)**: $\sim 0.22$ (FELsim DPW has M56 from physical wedge
  length; COSY DIL edge kick does not — a model difference, not a bug)

### Transfer matrix trace

With re-optimized quad currents (corrected DPW):
- $|\operatorname{Tr}_x/2| = 0.39$ (stable)
- RMS $= 8.2 \times 10^{-2}$

## Diagnostic scripts

All in `backend/test/results/W10/`:

| Script | Purpose |
|--------|---------|
| `diag_element_matrices.py` | Element-by-element COSY vs FELsim matrix comparison |
| `diag_transfer_map.py` | Full transfer matrix comparison (runs actual COSY) |
| `diag_beamline_compare.py` | Element types, lengths, and quad parameters |
| `diag_cumulative_map.py` | Cumulative FELsim trace evolution |
| `diag_edge_angles.py` | COSY tracking with/without edge angles |

## Timeline

1. Bug present since the original `beamline.py` commit (`98b1639`, initial import)
2. Never caught because:
   - FC dipoles ($+1.5°$) were unaffected
   - Previous W4 cross-validation used COSY optimization (self-consistent)
   - Y-plane was already marginally unstable, masking the additional error
3. Exposed during W10 study when FELsim-optimized currents were used in COSY
   particle tracking, revealing $\sigma_z$ blowup (31 ps vs expected 2 ps)
