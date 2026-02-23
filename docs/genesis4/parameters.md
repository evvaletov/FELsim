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

| Parameter | Value | Notes |
|-----------|-------|-------|
| `ngrid` | 129 | Converged (tested 65–257) |
| `dgrid` | 8 mm (half-grid-size) | Converged (tested 3–8 mm) |
| `sample` | **1** | **Critical for superradiance** (default 3 insufficient) |
| `slen` | 2.5 mm (time window) | 10 mm under investigation |
| Number of passes | 400 | |
| Desynchronization $d$ | 0.0 (default) | |

The desynchronization parameter $d = 2 \Delta L / S$ where $S$ is the
slippage length and $\Delta L$ is the cavity shortening.

### The SAMPLE Parameter

`SAMPLE` controls the number of radiation wavelengths per time slice:
$\Delta s = \texttt{SAMPLE} \times \lambda_0$. This is the most important
parameter for resolving superradiant dynamics:

- **`SAMPLE=3`** (default): $\Delta s \approx 10\,\mu\text{m}$ (33 fs).
  Too coarse for sub-picosecond pulse narrowing → normal saturation at 14 kW.
- **`SAMPLE=1`**: $\Delta s \approx 3.3\,\mu\text{m}$ (11 fs).
  Resolves single-wavelength structures → superradiance at 597 MW.

GINGER-3D's ADI solver naturally resolves sub-wavelength dynamics (no
SAMPLE equivalent), which is why the published results show superradiance
without special parameter tuning.

### Carrier Wavelength

Genesis4 uses the slowly-varying envelope approximation (SVEA). The carrier
wavelength $\lambda_0$ must be set near the **gain peak** (3.29 μm), not the
resonant wavelength (3.229 μm). With the resonant wavelength as carrier, the
gain peak falls outside the SVEA bandwidth and the FEL does not lase.

GINGER-3D (non-SVEA) has no such constraint and uses the resonant wavelength.
