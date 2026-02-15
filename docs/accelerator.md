# UH MkV FEL Accelerator

Overview of the full accelerator layout from electron source to undulator.

## Beamline Layout

```
TM010 RF cavity (thermionic cathode, 2856 MHz)
  → 2 quads
  → Alpha-Magnet (electrostatic kickers, knife-edge slit)
  → 2 quads
  → LINAC (S-band, 2856 MHz, → ~40 MeV)
  → Transport line (diagnostic chicane, 23 quads, 4 dipoles)
  → FEL (undulator, 47 periods, λ_u = 2.3 cm)
```

## Machine Parameters

### Electron Source

| Parameter | Value | Notes |
|-----------|-------|-------|
| Cavity mode | TM010 | Quarter-wave (λ/4) geometry |
| RF frequency | 2856 MHz | S-band |
| Cathode | Thermionic | Richardson-Dushman emission |
| Cathode radius | 1.57 mm | |
| Exit energy | ~0.6–1.1 MeV | Distribution; nominal 1 MeV |
| Accelerating voltage | V_peak = 0.937 MV | From COMSOL field map |

### Alpha-Magnet

| Parameter | Value | Source |
|-----------|-------|--------|
| Gradient $G$ | 0.103754 T/(m·A) | Design |
| Typical current | 17.5 A | Logbook |
| Slit width | 16 mm | Logbook |

The alpha-magnet provides energy selection (knife-edge slit) and bunch
compression. The positive energy chirp from the gun is compensated by the
alpha-magnet's negative $R_{56}$.

### Linac

S-band (2856 MHz) accelerating structure. Accelerates from ~1 MeV to the
nominal 40 MeV operating energy.

### Transport Line

The diagnostic chicane and beam transport from linac exit to undulator
entrance:

- 118 elements over ~12.8 m
- 23 quadrupoles (QPF/QPD)
- 4 dipoles (MkIII chicane)
- Correctors (vertical/horizontal steering)
- Diagnostics (BPM, OTR screens, spectrometer)

### Undulator

| Parameter | Value |
|-----------|-------|
| Period $\lambda_u$ | 2.3 cm |
| $K$ parameter | 1.2 |
| $a_w = K/\sqrt{2}$ | 0.849 |
| Number of periods | 47 |
| Undulator length | 1.081 m |
| Radiation wavelength | 3.29 μm (gain peak) |

### Optical Cavity

| Parameter | Value |
|-----------|-------|
| Cavity length | 2.0469 m |
| Mirror radius of curvature | 1.3 m |
| Configuration | Confocal pair |
| Power reflectivity | 93% |

## Beam Parameters at Linac Exit

These are the nominal parameters used in transport line simulations:

| Parameter | Value | Notes |
|-----------|-------|-------|
| Kinetic energy | 40 MeV | $\gamma \approx 78.3$ |
| Normalized emittance $\varepsilon_n$ | 8 π·mm·mrad | Baseline |
| Energy spread $\sigma_E/E$ | 0.5% | Hadmack RSI 84, 063302 |
| Bunch length | 0.5–2 ps | 2 ps nominal, 0.5 ps for short-bunch studies |
| Peak current | 30 A | |
| Energy chirp $h$ | 0 | Default; $h = 5 \times 10^9$ /s exploratory |

### Twiss Targets at Undulator Entrance

From arXiv:2510.14061v1 Table I:

| Plane | $\beta$ (m) | $\alpha$ |
|-------|-------------|---------|
| $x$ | 1.4 | 0.4714 |
| $y$ | 0.24 | 0.0 |

## References

- Wiedemann, *Particle Accelerator Physics* (2015) — pp. 449 (alpha-magnet), 628 (linac)
- Hadmack, Rev. Sci. Instrum. **84**, 063302 (2013) — energy spread measurement
- arXiv:2510.14061v1 — UH FEL undulator Twiss targets (Table I)
- UH-044 by Mike Hadmack — MkIII dipole modelling (pp. 48–81)
