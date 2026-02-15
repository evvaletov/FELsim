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
| `genesis4/uhfel_ss.in` | Steady-state single-pass input |
| `genesis4/uhfel_td.in` | Time-dependent single-pass input |
| `src/Main.py` | Original oscillator code (Levi Fisher, GINGER-3D) |

## Current Work

- Reproducing Fig. 2 of arXiv:2510.14061v1 (peak power vs pass number)
- Single-pass Genesis4 vs GINGER-3D comparison
- Desynchronization scans ($d$ parameter)
- Runtime optimization (~32 s/pass with 4 MPI ranks on Koa)

```{toctree}
:maxdepth: 1

parameters
running
```
