# Beam Parameters at Linac Exit — Email 2026-02-04

Source: Email from Niels Bidault, 4 Feb 2026, regarding initial Twiss parameters.

## Transverse Planes

| Parameter | Value | Notes |
|-----------|-------|-------|
| ε_n | 8 π·mm·mrad | Normalized emittance |
| σ_x | 0.8 mm | |
| σ_y | 0.8 mm | |
| α_x, α_y | 0 | No x-y correlations |

Assumption based on previous papers and the previous group's documentation.

## Longitudinal

| Parameter | Value | Notes |
|-----------|-------|-------|
| Energy | 40 MeV | Range: 20–45 MeV |
| f_RF | 2856 MHz | S-band |
| Bunch length | 2 ps | RMS |
| σ_E | 0.5% | Conservative; 0.3% from Hadmack RSI 84, 063302 (2013) without RF feedforward |

Reference: M. Hadmack, Rev. Sci. Instrum. **84**, 063302 (2013); doi:10.1063/1.4809938.

## Longitudinal Correlations (Chirp)

> "I ran some tests with different values of the correlation between energy and
> ToF, for example, using a chirp h = 5e9 (1/s), which remains to be determined
> by injector simulations. You can assume different correlations between dE and
> dT, or fix alpha_z = 0 for the moment. We will have to finalize the injector
> simulations to see how the cavity, alpha-magnet, and linac transform and rotate
> the longitudinal dynamics, in order to obtain a better estimate of this phase
> space."

**Conclusion:** h = 5e9 /s was exploratory; the actual chirp remains to be
determined by injector simulations. Default: h = 0 (α_z = 0).
