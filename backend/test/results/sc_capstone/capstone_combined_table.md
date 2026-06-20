## Reduced no-dipole SC capstone — DA-FMM vs xsuite-frozen (combined)

Section [32,46) (drift+quad, 6 quads, 1.65 m). Common distribution N_p=6000,
eps_n=8 mm.mrad, sigma_delta=0. Optics held fixed (k0 at 45 MeV) so the 1 MeV
rows isolate space-charge scaling, not a focusing-instability blow-up.
SC: xsuite frozen-Gaussian (nonlinear Bassetti-Erskine, rms-self-consistent) vs DA-FMM N-body 1/r treecode (resolves the actual profile).

| E [MeV] | Q [nC] | DA-FMM dEx | xsuite dEx | DA/xs ratio |
|---|---|---|---|---|
| 45 | 0.1 | -0.017% | -0.003% | 6.19 |
| 45 | 0.3 | -0.005% | -0.006% | 0.89 |
| 45 | 1 | +0.255% | +0.012% | 21.66 |
| 45 | 3 | +2.810% | +0.292% | 9.63 |
| 45 | 10 | +26.949% | +4.315% | 6.25 |
| 1 | 0.01 | +0.843% | +3.992% | 0.21 |
| 1 | 0.03 | +8.901% | +42.544% | 0.21 |
| 1 | 0.1 | +49.866% | +245.063% | 0.20 |
| 1 | 0.3 | +150.744% | +790.287% | 0.19 |

SC-off control (both ~0): max |DA-FMM|=5.9e-13%, max |xsuite|=1.2e-12%.
