# Reduced no-dipole SC capstone — DA-FMM vs xsuite-frozen

Section (32, 46) (drift+quad, 6 quads, 1.65 m) of the FELsim line.
Common distribution: N_p=6000, eps_n=8.0 mm.mrad, sigma_delta=0, seed=20260619. Space charge: frozen-Gaussian (xsuite) vs full N-body 1/r treecode (DA-FMM).

| E [MeV] | Q [nC] | DA-FMM dEx | DA-FMM dEy | xsuite dEx | xsuite dEy | DA-vs-xs (x) |
|---|---|---|---|---|---|---|
| 45 | 0.1 | -0.017% | +0.006% | -0.003% | -0.003% | 0.014 pp |
| 45 | 0.3 | -0.005% | +0.061% | -0.006% | -0.007% | 0.001 pp |
| 45 | 1 | +0.255% | +0.310% | +0.012% | +0.001% | 0.243 pp |
| 45 | 3 | +2.810% | +3.215% | +0.292% | +0.208% | 2.518 pp |
| 45 | 10 | +26.949% | +26.391% | +4.315% | +3.988% | 22.633 pp |
| 1 | 0.1 | +8221.117% | -100.000% | +nan% | +nan% | nan pp |
| 1 | 0.3 | -100.000% | +24747185.628% | +nan% | +nan% | nan pp |
| 1 | 1 | -100.000% | -100.000% | +nan% | +nan% | nan pp |
| 1 | 3 | -100.000% | -100.000% | +nan% | +nan% | nan pp |
| 1 | 10 | -100.000% | +559967051.470% | +nan% | +nan% | nan pp |

SC-off control (both engines should be ~0%): max |DA-FMM off| = 1.00e+02%, max |xsuite off| = 9.79e-13%.
