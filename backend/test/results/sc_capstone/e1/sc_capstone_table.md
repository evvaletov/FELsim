# Reduced no-dipole SC capstone — DA-FMM vs xsuite-frozen

Section (32, 46) (drift+quad, 6 quads, 1.65 m) of the FELsim line.
Common distribution: N_p=6000, eps_n=8.0 mm.mrad, sigma_delta=0, seed=20260619. Space charge: frozen-Gaussian (xsuite) vs full N-body 1/r treecode (DA-FMM).

| E [MeV] | Q [nC] | DA-FMM dEx | DA-FMM dEy | xsuite dEx | xsuite dEy | DA-vs-xs (x) |
|---|---|---|---|---|---|---|
| 1 | 0.01 | +0.843% | +1.191% | +3.992% | +3.656% | 3.149 pp |
| 1 | 0.03 | +8.901% | +8.590% | +42.544% | +40.223% | 33.643 pp |
| 1 | 0.1 | +49.866% | +40.382% | +245.063% | +237.569% | 195.197 pp |
| 1 | 0.3 | +150.744% | +120.519% | +790.287% | +785.048% | 639.543 pp |

SC-off control (both engines should be ~0%): max |DA-FMM off| = 5.87e-13%, max |xsuite off| = 1.15e-12%.
