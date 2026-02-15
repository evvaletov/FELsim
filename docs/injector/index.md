# Injector Simulation

The injector chain — RF gun, alpha-magnet, and linac — sits upstream of the
FELsim transport line. Injector simulations are at an early stage, with
the RF gun being the most developed component.

## RF Gun (Niels Bidault)

Niels Bidault is simulating the thermionic RF gun cavity using RF-Track
with COMSOL field maps.

**Repository:** [Ioniels/RF-TRACK_UH_ThermionicRFgun](https://github.com/Ioniels/RF-TRACK_UH_ThermionicRFgun)

### Gun Parameters

| Parameter | Value |
|-----------|-------|
| Cavity mode | TM010 |
| RF frequency | 2856 MHz |
| Physical length | λ/4 = 26.24 mm |
| Peak $|E_z|$ | ~3.8 × 10⁷ V/m (COMSOL) |
| Accelerating voltage | $V_\text{peak}$ = 0.937 MV |
| Cathode type | Thermionic (Richardson-Dushman) |
| Cathode radius | 1.57 mm |
| Work function | 2.1 eV |
| Temperature | 1700 K |

### Gun Exit Beam

At the gun exit (no space charge, $\phi = 100°$):

| Parameter | Value |
|-----------|-------|
| Mean kinetic energy | ~356 keV |
| $\sigma_x$, $\sigma_y$ | ~2.85 mm |
| Energy spread | ~15% |
| Temporal spread | ~350 ps |
| Chirp | Strong positive |

### Gap to FELsim Parameters

The gun output is far from the linac-exit parameters used in FELsim
transport simulations. The alpha-magnet and linac perform major
transformations:

| Parameter | Gun Exit | FELsim (linac exit) | Transformation |
|-----------|----------|--------------------:|----------------|
| Energy | 0.356 MeV | 40 MeV | Linac accelerates (112×) |
| $\varepsilon_n$ | 0.08–0.84 π·mm·mrad | 8 | Emittance growth in chain |
| $\sigma_t$ | ~350 ps | 2 ps | Alpha-magnet compresses |
| $\sigma_E/E$ | ~15% | 0.5% | Alpha-magnet selects + linac reduces |
| Chirp | Strong positive | ~0 | Alpha-magnet compensates |

## Alpha-Magnet (Planned)

The alpha-magnet provides energy selection and bunch compression. Mason
McMahon is modelling it as a series of dipoles.

| Parameter | Value | Source |
|-----------|-------|--------|
| Gradient $G$ | 0.103754 T/(m·A) | Design |
| Typical current | 17.5 A | Logbook |
| Slit width | 16 mm | Logbook |

## Linac (Planned)

S-band (2856 MHz) accelerating structure, ~1 MeV → ~40 MeV. Niels Bidault
is working on analytical calculations. A simulation model does not yet exist.

## References

- Wiedemann, *Particle Accelerator Physics* (2015) — pp. 449 (alpha-magnet), 628 (linac)
