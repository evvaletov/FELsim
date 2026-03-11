# BUG: RF-Track Adapter Edge Kick Sign Inversion

| Field | Value |
|-------|-------|
| **Date found** | 2026-03-02 |
| **Date fixed** | 2026-03-02 |
| **Severity** | Critical — chicane and spectrometer edge kicks inverted |
| **Affected file** | `backend/rftrackAdapter.py`, method `_annotate_dipole_edges()` |
| **Affects** | All RF-Track simulations with negative-angle dipoles (chicane C1–C4, spectrometer SP3–SP4) |
| **Does NOT affect** | FELsim (separate bug fixed in `beamline.py`, see [DPW edge kick sign](dpw-edge-kick-sign.md)), COSY adapter |

## Symptom

After implementing the analytical sector-bend correction (see
[SBend P/δ bug](rftrack-sbend-p-delta.md)), the RF-Track adapter produced
severely incorrect y-plane optics through the chicane:

| Parameter | Expected | With signed K0 | With abs(K0) |
|-----------|----------|----------------|--------------|
| $\beta_y$ at undulator | 0.242 m | ~87 m | ~0.055 m |
| RMS | ~1.4 | 39.8 | 1.41 |

The x-plane was relatively unaffected because the chicane DPW wedge angles
are small ($\eta \approx 2°$), making the x-kick $\tan\eta / R$ small
in absolute terms.  The y-plane kick $-\tan(\eta - \phi) / R$ is more
sensitive because $\phi$ (the fringe-field integral) is comparable to $\eta$.

## Root Cause

The `_annotate_dipole_edges()` method computed the curvature parameter for
each DPW element as:

```python
K0 = np.radians(dph_angle) / dph_length  # signed!
```

For chicane dipoles with $\theta = -4°$, this gives $K0 < 0$.  The thin-lens
edge kick uses:

$$
K_{1L} = -|K_0| \tan(\eta)
$$

but the code was using $-K_0 \cdot \tan(\eta)$.  With $K_0 < 0$, the kick
flips sign: defocusing edges become focusing and vice versa.

This is the RF-Track adapter's analogue of the
[DPW edge kick sign](dpw-edge-kick-sign.md) bug in `beamline.py`.  The root
cause is the same: using signed curvature where the standard beam optics
convention requires $|\rho|$.

### Affected elements

| Dipole | Bending angle | K0 (signed) | K0 (unsigned) |
|--------|---------------|-------------|---------------|
| C1–C4 (chicane) | $-4.0°$ | $-0.349$ m$^{-1}$ | $+0.349$ m$^{-1}$ |
| SP3, SP4 (spectrometer) | $-11.25°$ | $-0.982$ m$^{-1}$ | $+0.982$ m$^{-1}$ |
| FC1, FC2 (LINAC) | $+1.5°$ | $+0.131$ m$^{-1}$ | $+0.131$ m$^{-1}$ |
| BC1, BC2 (bunch compressor) | $+5.0°$ | $+0.436$ m$^{-1}$ | $+0.436$ m$^{-1}$ |

Positive-angle dipoles were unaffected (signed = unsigned for $\theta > 0$).

## Fix

```python
# BEFORE (buggy)
K0 = np.radians(dph_angle) / dph_length

# AFTER (fixed)
K0 = abs(np.radians(dph_angle) / dph_length)
```

This matches FELsim's `beamline.py` convention: $R = L / |\theta|$.

## Verification

After the fix, the full-beamline RFT-val (FELsim currents evaluated in
RF-Track) gives:

| Metric | Before fix | After fix |
|--------|-----------|-----------|
| RMS | 39.8 | 1.41 |
| $\beta_y$ | ~87 m | ~0.29 m |
| $\beta_x$ | ~0.5 m | ~5.1 m |

The remaining RMS = 1.41 is due to genuine model differences between FELsim
(transfer matrices with triangle-rule fringe) and RF-Track (analytical
sector-bend correction with thin-quad edge kicks), not a sign error.

## Relationship to DPW edge kick sign bug

This bug is conceptually identical to the [DPW edge kick sign](dpw-edge-kick-sign.md)
bug in `beamline.py`, but in a different file (`rftrackAdapter.py`).  Both
use signed curvature where unsigned is required.  The FELsim bug was found
on 2026-03-01; this adapter bug was found the following day when implementing
the RF-Track dipole workaround.
