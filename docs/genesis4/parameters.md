# Machine and Simulation Parameters

Parameters for the Genesis4 oscillator simulation, from arXiv:2510.14061v1
and the `oscillator.py` configuration.

## Electron Beam

| Parameter | Value | Notes |
|-----------|-------|-------|
| Beam energy | 40 MeV | $\gamma \approx 78.3$ |
| Energy spread | 0.5% | Fractional |
| Normalized emittance | 8 μm | = 8 π·mm·mrad |
| Peak current | 30 A | |
| Bunch FWHM | 2 ps | |

## Twiss at Undulator Entrance

| Plane | $\beta$ (m) | $\alpha$ |
|-------|-------------|---------|
| $x$ | 1.4 | 0.4714 |
| $y$ | 0.24 | 0.0 |

These are the same targets used in the FELsim transport line optimization.

## Undulator

| Parameter | Value |
|-----------|-------|
| Period $\lambda_u$ | 2.3 cm |
| $K$ | 1.2 |
| $a_w = K/\sqrt{2}$ | 0.849 |
| Number of periods | 47 |
| Undulator length | 1.081 m |
| Radiation wavelength | 3.29 μm (gain peak) |

The resonant wavelength for 40 MeV / $K = 1.2$ is ~3.229 μm. The FEL gain
peak is red-shifted to ~3.29 μm due to energy spread effects.

## Optical Cavity

| Parameter | Value |
|-----------|-------|
| Cavity length | 2.0469 m |
| Mirror radius of curvature | 1.3 m |
| Configuration | Confocal pair |
| Focal length | 0.65 m |
| Power reflectivity | 93% |

## Genesis4 Simulation Grid

| Parameter | Value |
|-----------|-------|
| `ngrid` | 129 |
| `dgrid` | 8 mm (half-grid-size) |
| `sample` | 3 |
| `slen` | 2.5 mm (time window) |
| Number of passes | 400 |
| Desynchronization $d$ | 0.0 (default) |

The desynchronization parameter $d = 2 \Delta L / S$ where $S$ is the
slippage length and $\Delta L$ is the cavity shortening.
