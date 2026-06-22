# Reduced no-dipole SC capstone — DA-FMM vs xsuite-frozen

Section (32, 46) (drift+quad, 6 quads, 1.65 m) of the FELsim line.
Common distribution: N_p=6000, eps_n=8.0 mm.mrad, sigma_delta=0, seed=20260619. Space charge: frozen-Gaussian (xsuite) vs full N-body 1/r treecode (DA-FMM).

| E [MeV] | Q [nC] | DA-FMM dEx | DA-FMM dEy | xsuite dEx | xsuite dEy | cosy-pic dEx | cosy-pic dEy | DA-vs-xs (x) |
|---|---|---|---|---|---|---|---|---|
| 45 | 0.3 | -0.005% | +0.061% | -0.006% | -0.007% | -0.005% | -0.004% | 0.001 pp |
| 45 | 1 | +0.255% | +0.310% | +0.012% | +0.001% | +0.012% | +0.009% | 0.243 pp |
| 45 | 3 | +2.810% | +3.215% | +0.292% | +0.208% | +0.267% | +0.193% | 2.518 pp |
| 45 | 10 | +26.949% | +26.391% | +4.315% | +3.988% | +3.441% | +2.850% | 22.633 pp |

SC-off control (both engines should be ~0%): max |DA-FMM off| = 4.57e-13%, max |xsuite off| = 9.79e-13%.
