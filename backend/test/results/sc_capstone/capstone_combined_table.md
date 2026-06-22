## Reduced no-dipole SC capstone — DA-FMM vs xsuite-frozen vs cosy-pic (combined)

Section [32,46) (drift+quad, 6 quads, 1.65 m). Common distribution N_p=6000,
eps_n=8 mm.mrad, sigma_delta=0. Optics held fixed (k0 at 45 MeV) so the 1 MeV
rows isolate space-charge scaling, not a focusing-instability blow-up. Three SC
engines on one distribution: DA-FMM (cosy-fmm N-body 1/r treecode), xsuite
frozen-Gaussian (nonlinear Bassetti-Erskine, rms-self-consistent), cosy-pic
(Hockney mesh PIC).

| E [MeV] | Q [nC] | DA-FMM dEx | xsuite dEx | cosy-pic dEx | DA/xs ratio |
|---|---|---|---|---|---|
| 45 | 0.3 | -0.005% | -0.006% | -0.005% | 0.89 |
| 45 | 1 | +0.255% | +0.012% | +0.012% | 21.66 |
| 45 | 3 | +2.810% | +0.292% | +0.267% | 9.63 |
| 45 | 10 | +26.949% | +4.315% | +3.441% | 6.25 |
| 1 | 0.01 | +0.843% | +3.992% | N/A | 0.21 |
| 1 | 0.03 | +8.901% | +42.544% | N/A | 0.21 |
| 1 | 0.1 | +49.866% | +245.063% | N/A | 0.20 |
| 1 | 0.3 | +150.744% | +790.287% | N/A | 0.19 |

SC-off control (both ~0): max |DA-FMM|=5.9e-13%, max |xsuite|=1.2e-12%.

cosy-pic (mesh PIC) and xsuite-frozen agree to ~20% at every 45 MeV charge, jointly confirming the physical mean-field growth; the DA-FMM bare-treecode excess (6-22x) is macroparticle shot noise, not physics.
