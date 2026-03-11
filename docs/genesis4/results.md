# Oscillator Results

Summary of Genesis4 oscillator simulation results as of February 2026.

## Superradiance Breakthrough

The key finding is that Genesis4's `SAMPLE` parameter controls whether the
simulation resolves superradiant pulse dynamics. With `SAMPLE=1`
(single-wavelength time slices), the oscillator reaches the superradiant
regime; with the default `SAMPLE=3`, it reaches only conventional
saturation.

## Parameter Sweep

All runs use shot-noise self-start, $d = 0$ (zero desynchronization),
`NGRID=129`, `DGRID=8mm`.

| Configuration | SAMPLE | SLEN | Peak Power | Pulse FWHM | Status |
|--------------|--------|------|------------|------------|--------|
| Baseline | 3 | 2.5 mm | 14 kW | ~1.2 ps | Saturated (pass ~150) |
| SLEN only | 3 | 10 mm | 10.8 kW | ~1.2 ps | Saturated |
| **SAMPLE only** | **1** | **2.5 mm** | **597 MW** | **0.15 ps** | **Near saturation (pass 400)** |
| Both | 1 | 10 mm | 184 MW+ | TBD | Running (pass 210+) |
| Paper (GINGER-3D) | — | — | 80 MW | 0.31 ps | Published |

**Conclusion**: `SAMPLE=1` is necessary and sufficient for superradiance.
Increasing `SLEN` alone has no effect.

## Physics Summary

### Normal Saturation (SAMPLE=3)

With 3λ time slices (~33 fs resolution), Genesis4 correctly computes
conventional FEL saturation. The optical pulse is comparable in duration
to the electron bunch (~2 ps), and peak power is limited by electron
trapping in the ponderomotive potential. This produces ~14 kW — a correct
result for the resolved physics, but not superradiance.

### Superradiance (SAMPLE=1)

With 1λ time slices (~11 fs resolution), the simulation resolves
sub-wavelength dynamics that drive superradiant pulse narrowing:

1. **Slippage-driven amplification**: The leading edge of the optical
   pulse overlaps with fresh, unmodulated electrons each pass.
2. **Pulse narrowing**: Only the leading edge is efficiently amplified;
   the trailing portion overlaps with depleted electrons.
3. **Power scaling**: Peak power grows much faster than pulse energy
   because the pulse is compressing.

The result: 597 MW peak power in a 0.15 ps pulse (vs. 14 kW in ~1.2 ps).

### Discrepancy with Paper

Our Genesis4 result overshoots the published GINGER-3D value (80 MW) by
~7.5×. Possible causes:

- **Not fully saturated**: Gain is still ~0.1%/pass at pass 400.
- **Circular-shift wrapping**: Our FFT circular shift wraps the temporal
  field periodically; GINGER-3D absorbs field at the boundaries.
- **SVEA vs. ADI**: Different physics approximations in the deep
  nonlinear regime.

Under investigation via `SAMPLE=1` + `SLEN=10mm` run and planned absorbing
boundary tests.

## Resolved Bugs

| # | Issue | Impact | Fix |
|---|-------|--------|-----|
| 1 | Missing slippage compensation | Field exits window in ~16 passes | Shift field by $S$ per pass |
| 2 | `interp1d` zero-fill boundary loss | 6.2% energy loss/pass, shot noise below threshold | FFT circular shift |
| 3 | Carrier at resonant λ (3.229 μm) | Gain peak outside SVEA bandwidth | Use gain peak λ (3.29 μm) |
| 4 | Mirror phase aliasing | Excess cavity loss | Exact Helmholtz kernel (no aperture needed) |

## Koa Run Log

| Job | Passes | Config | Time | Result |
|-----|--------|--------|------|--------|
| 10753099 | 500 | SAMPLE=3, baseline | 2.9h | 14 kW saturation |
| 10753100 | 300 | SAMPLE=3, 10W seed | 1.6h | 50 kW limit cycle |
| 10756522 | 400 | SAMPLE=1, SLEN=2.5mm | 11h | **597 MW** |
| 10756523 | 308 | SAMPLE=3, SLEN=10mm | 12h (timeout) | 10.8 kW |
| 10756524 | 400 | SAMPLE=1, SLEN=10mm | ~30h | Running |
