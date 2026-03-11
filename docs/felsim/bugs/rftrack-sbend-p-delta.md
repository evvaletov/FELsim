# BUG: RF-Track SBend P/δ Confusion (Upstream)

| Field | Value |
|-------|-------|
| **Date found** | 2026-03-02 |
| **Status** | Open (upstream); workaround implemented |
| **Severity** | Critical — all RF-Track dipole body tracking broken |
| **Affected code** | RF-Track v2.5.5 (proprietary, A. Latina, CERN) |
| **Workaround file** | `backend/rftrackAdapter.py` |
| **Affects** | All RF-Track simulations involving `SBend` or `RBend` elements |
| **Does NOT affect** | FELsim (first-order), COSY INFINITY |

## Symptom

RF-Track's `SBend` element produces catastrophically wrong tracking output:

- **Without `P_Q` argument:** `SBend(L, angle)` acts as identity — no bending
  at all.  The `set_K0()` method accepts a value silently but has no effect.
- **With `P_Q` argument:** `SBend(L, angle, P_Q)` produces ~910 mm horizontal
  displacement for a 1 mm off-axis particle at 40 MeV (expected: sub-mm).

A single-particle test tracking through a 22.5° sector bend:

| Configuration | $x_\text{out}$ (mm) | Expected $x_\text{out}$ (mm) |
|---------------|---------------------|------------------------------|
| `SBend(L, θ)` | 1.000 (identity) | 0.924 |
| `SBend(L, θ, P_Q)` | 910.6 | 0.924 |
| Analytical $\cos\theta\cdot x$ | 0.924 | — |

## Root Cause

RF-Track's `SBend` transfer matrix reads Bunch6d's 6th coordinate — which
contains the particle's **absolute momentum** $P$ in MeV/c — as if it were
the **momentum deviation** $\delta = \Delta P / P_0$.  For a 40 MeV electron,
$P_0 \approx 40.5$ MeV/c, so the code interprets $\delta \approx 40.5$
(a 4050% deviation), producing enormous dispersive kicks.

Evidence: the displacement scales linearly with $P_0$:
- 40 MeV electron ($P_0 = 40.5$): $\Delta x \approx 910$ mm
- If $P_0$ were 4 MeV/c: $\Delta x$ drops proportionally

The `set_K0()` setter (without `P_Q` in the constructor) does not activate
the bending field at all — the element reduces to a drift.  All other API
methods (`set_h`, `set_Bfield`, `set_E1`, `set_E2`) accept values silently
but have no effect on tracking output.

## Workaround

The RF-Track adapter implements an **analytical sector-bend correction** that
recovers correct dipole physics without using `SBend`:

1. **Track dipole body as `Drift(L)`:** Preserves the y-plane and path length
   correctly (dipole body is straight-through in thin-element limit).

2. **Apply analytical correction:** After drift tracking, apply
   $M_\text{correction} = M_\text{sector} \times M_\text{drift}^{-1}$
   to the $(x, x')$ coordinates of each particle.

3. **Add dispersion and path-length terms:**
   - $R_{16} = \rho(1 - \cos\theta)$ and $R_{26} = \sin\theta / \rho$ for
     energy-dependent orbit shift
   - $R_{56} = -L(1 - \sin\theta / \theta)$ for path-length vs energy

The correction is exact to all significant figures for the sector-bend
transfer matrix.  It is applied in `_apply_sector_bend_correction()`, called
from `track_elements()` and `_track_segmented()`.

### Implementation details

```python
# In rftrackAdapter.py
def _apply_sector_bend_correction(self, particles, angle_rad, length):
    """Apply M_sector × M_drift⁻¹ correction to (x, x') after drift tracking."""
    theta = angle_rad
    rho = length / theta
    ct, st = np.cos(theta), np.sin(theta)

    # M_sector (x, x' block)
    M_sec = np.array([[ct, rho * st], [-st / rho, ct]])
    # M_drift (x, x' block)
    M_drift = np.array([[1, length], [0, 1]])
    # Correction
    M_corr = M_sec @ np.linalg.inv(M_drift)

    x = particles[:, 0]
    xp = particles[:, 1]
    particles[:, 0] = M_corr[0, 0] * x + M_corr[0, 1] * xp
    particles[:, 1] = M_corr[1, 0] * x + M_corr[1, 1] * xp

    # Dispersion (δ = ΔP/P₀ from RF-Track 6th coordinate)
    delta = particles[:, 5] / self._P0 - 1  # convert absolute P to δ
    particles[:, 0] += rho * (1 - ct) * delta
    particles[:, 1] += st * delta

    # R56 path-length correction
    particles[:, 4] += -length * (1 - st / theta) * delta
```

### Segmented tracking

For beamlines with interleaved dipoles and other elements, `_track_segmented()`
groups consecutive non-dipole elements into sub-lattices, tracks them with
RF-Track's native tracking, and applies the analytical correction at each
dipole.  DPW (wedge) elements are handled separately as thin-lens edge kicks.

## Verification

### Single-particle test

A 1 mm off-axis particle tracked through a 22.5° sector bend:

| Coordinate | Analytical | RF-Track + correction | Difference |
|------------|-----------|----------------------|------------|
| $x$ (mm) | 0.92388 | 0.92388 | $< 10^{-6}$ |
| $x'$ (mrad) | −3.82683 | −3.82683 | $< 10^{-6}$ |

### Full beamline validation

With the analytical correction (ε_n = 8, 3 restarts):

| Metric | Before correction | After correction |
|--------|-------------------|------------------|
| RFT-val RMS | 2.00 | 1.41 |
| RFT-opt RMS | 8.4e-2 | 8.4e-2 |

The RFT-val improvement confirms that the correction removes
the dominant source of model disagreement.  The residual RMS reflects genuine
model differences (FELsim uses triangle-rule fringe; RF-Track + correction
uses pure sector bend).

## MWE

`backend/test/rftrack_sbend_bug_mwe.py` — standalone script demonstrating
the bug with RF-Track v2.5.5.  Tracks on-axis and off-axis particles through
`SBend` with and without `P_Q`, compares with analytical sector-bend matrix.

## Timeline

1. **2026-02-24 (W8):** RF-Track adapter first used for Stage 11 optimization.
   SBend element mapped to `Drift(0)` as a temporary workaround —
   RFT-val RMS = 50.5 (no body focusing, no dispersion).
2. **2026-03-02 (C5):** P/δ confusion identified. Analytical sector-bend
   correction implemented.  RFT-val RMS drops to 1.41.
3. **TODO (C5-BUG):** File bug report with RF-Track maintainer (CERN GitLab
   or email to A. Latina).
