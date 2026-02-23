# UH FEL Simulation

Simulation tools for the University of Hawai'i at Mānoa Mark V Free Electron
Laser. The project covers the full electron beam path from injector to
undulator, with three main simulation domains:

- **Transport line** (FELsim) — beam dynamics from linac exit to undulator
  entrance, including diagnostic chicane, quadrupole matching, and cross-validation
  across multiple simulation codes
- **Undulator / FEL** (Genesis4) — multi-pass oscillator simulation using
  Genesis4 for undulator physics and Fourier-optics for optical cavity propagation
- **Injector** — thermionic RF gun, alpha-magnet, and linac modelling (early stage)

## Status

| Domain | Codes | Maturity |
|--------|-------|----------|
| Transport line | FELsim (1st order), COSY INFINITY, RF-Track | Production — optimization studies ongoing |
| FEL oscillator | Genesis4 + Fourier-optics | Active — superradiance achieved (SAMPLE=1), quantitative comparison ongoing |
| Injector | RF-Track (gun only) | Early — gun cavity simulation by Niels Bidault |

## Repository

- **FELsim**: [evvaletov/FELsim](https://github.com/evvaletov/FELsim)
  (upstream: [komochristian/FELsim](https://github.com/komochristian/FELsim))
- **Genesis4 oscillator**: `/home/evaletov/UH/GENESIS/UHFEL_undulator` (local)
- **RF gun**: [Ioniels/RF-TRACK_UH_ThermionicRFgun](https://github.com/Ioniels/RF-TRACK_UH_ThermionicRFgun)

```{toctree}
:maxdepth: 2
:caption: Contents

accelerator
felsim/index
genesis4/index
injector/index
reference/element-types
reference/people
reference/bibliography
```
