# UH Beam Dynamics — Initial Parameters & Code References

Source: Email thread "Re: UH beam dynamics", Niels Bidault ↔ Eremey Valetov,
Aug 7 – Sep 14, 2025.

## Initial Beam Parameters (from Niels, Aug 7 2025)

From `beamline_optimization.ipynb`:

```python
Energy = 40             # MeV
f = 2856e6              # Hz
bunch_spread = 2        # ps
energy_std_percent = 0.5  # 0.3% from Hadmack RSI 84, 063302 (2013); conservative 0.5%
h = 5e9                 # 1/s energy vs ToF correlations, to be assessed
epsilon_n = 8           # pi.mm.mrad transverse normalized emittance
x_std = 0.8             # mm
y_std = 0.8             # mm
nb_particles = 10000
```

**Key note:** h = 5e9 was explicitly marked "to be assessed" from the very beginning.

## Code References

- Original FELsim: https://github.com/k0m0code/FELsim/tree/main
- `excelElements.py` — recreates beamline from Excel spreadsheet
- `beamline_optimization.ipynb` — usage examples and beam input properties

## Enge Function Discussion (Aug 26 – Sep 14)

- MkIII chicane dipole Enge coefficients: `56.49, -50.79, 19.32, -3.621, 0.3315, -0.01193`
- These were **test values**, not fully adjusted to fringe field
- Mike Hadmack measured fringe field: B_max ~ 5330 gauss, dipole width L = 1.466" = 3.724 cm
- Chicane calibration factor = 1.510
- ∫B dx = 29973 gauss·cm; (B_max)·L = 19849 gauss·cm
- For θ = 11.25°: B_max(gauss) = 116.472 × [(K.E./MeV) + 0.511]
- Fringe field falloff too long for simple Enge function fitting
- COSY INFINITY MGE element suggested for direct field data use
- Displacement "Displ (cm)" — Niels to double-check if z-axis or arc length s
- Reference: UH-044 by Mike Hadmack, pp. 48–81 (MkIII dipole modelling)

## Framework Architecture (Sep 14)

Five simulation blocks:
1. **Gun cavity** — Analytical, GPT, other codes
2. **Injection line** — quads, alpha-magnet with knife edge (GPT, COSY, 1st order, Xsuite)
3. **Linac** — analytical, GPT, other codes
4. **Beam diagnostic chicane** — GPT, COSY, 1st order, Xsuite
5. **FEL** — Ginger3D (Siqi's group)

Framework capabilities (existing and planned):
- 6D beam distribution input (Gaussian, KV, user-defined)
- Beamline from Excel/CSV/JSON, web-app constructor
- 6D phase space visualization, figures of merit at every z
- Optimization of magnet values (Nelder-Mead, algebraic; consider Xopt)
- Future: virtual accelerator for ML/AI training
