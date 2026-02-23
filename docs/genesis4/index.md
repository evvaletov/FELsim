# Genesis4 FEL Oscillator

Multi-pass FEL oscillator simulation for the UH MkV FEL, reproducing and
extending results from arXiv:2510.14061v1 (Weinberg, Fisher, Li).

## Simulation Approach

Each oscillator pass consists of two steps:

1. **Genesis4** simulates a single undulator pass: the electron beam
   interacts with the radiation field over 47 undulator periods.
2. **Fourier-optics cavity propagation** (Python) transports the radiation
   field back through the optical cavity — two confocal mirrors with
   free-space propagation between them.

These steps repeat for $N$ passes (typically 300–400) until the radiation
field reaches saturation. The cavity model uses thin-lens Fourier-optics
propagation; Genesis4's built-in cavity model is not used.

## Key Files

| File | Purpose |
|------|---------|
| `genesis4/oscillator.py` | Multi-pass oscillator driver |
| `genesis4/uhfel.lat` | Genesis4 lattice file (undulator) |
| `genesis4/plot_diagnostics.py` | Diagnostic plotting (power, profiles, etc.) |
| `genesis4/uhfel_ss.in` | Steady-state single-pass input |
| `genesis4/uhfel_td.in` | Time-dependent single-pass input |
| `src/Main.py` | Original oscillator code (Levi Fisher, GINGER-3D) |

## Status (February 2026)

**Superradiance achieved** with `SAMPLE=1` (single-wavelength time slices):

| Configuration | Peak Power | Pulse FWHM | Regime |
|--------------|------------|------------|--------|
| `SAMPLE=3` (default) | 14 kW | ~1.2 ps | Normal saturation |
| `SAMPLE=1` | 597 MW | 0.15 ps | Superradiance |
| Paper (GINGER-3D) | 80 MW | 0.31 ps | Superradiance |

The `SAMPLE` parameter is critical: it controls whether Genesis4 can
resolve the sub-wavelength pulse dynamics that drive superradiant
narrowing. With 3λ slices (`SAMPLE=3`), the FEL reaches only conventional
saturation. With 1λ slices (`SAMPLE=1`), superradiance is unlocked.

The 7.5× power overshoot vs. the paper is under investigation (likely
related to periodic boundary conditions in the slippage compensation).

### Resolved Issues

1. **Slippage walk-off**: Field exits simulation window without compensation.
   Fixed by shifting the field backward by $S = N_u \lambda$ per pass.
2. **Shift boundary loss**: `interp1d` with zero-fill lost 6.2% of energy per
   pass, keeping shot noise below threshold. Fixed with FFT circular shift.
3. **Carrier wavelength**: Genesis4 SVEA requires $\lambda_0$ at the gain peak
   (3.29 μm), not the resonant wavelength (3.229 μm).
4. **Cavity kernel**: Exact Helmholtz kernel matches GINGER-3D (vs. Fresnel).

### Active Work

- `SAMPLE=1` + `SLEN=10mm` run on Koa (testing temporal window effects)
- Absorbing boundary conditions (to match GINGER-3D behavior)
- Desynchronization scans ($d$ parameter, paper Figs. 4–6)

```{toctree}
:maxdepth: 1

parameters
running
results
```
