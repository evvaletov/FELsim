# Beam-loading refinement - cell-resolved TW linac

Author: Eremey Valetov, 2026-06-20. Steady-state constant-gradient (CG)
fundamental-mode beam loading on the validated cell-resolved TW model.
Formulas panel-verified (Wangler / SLAC two-mile, S.Y. Lee).

| Quantity | Value |
|---|---:|
| CG beam-loading factor k(tau=0.57) | 0.2320 |
| Shunt impedance r | 53 Mohm/m |
| Energy loss slope (on crest) | 37.00 MeV/A |
| Average current for 3% droop | 32 mA |
| Loaded gradient sag E_b(L) at 100 mA | 3.021 MV/m |

Self-checks:
- numeric integral(E_b dz) = 3.6998 MeV vs closed-form I r L k = 3.6998 MeV (0.00%).
- group velocity vg(L)/vg(0): analytic CG e^-2tau = 0.320, aperture-geometry estimate 0.319 (consistent).

| I_avg (mA) | K_out model (MeV) | K_out analytic (MeV) |
|---:|---:|---:|
| 0 | 40.953 | 40.953 |
| 10 | 40.583 | 40.583 |
| 20 | 40.213 | 40.213 |
| 30 | 39.843 | 39.843 |
| 50 | 39.103 | 39.103 |
| 75 | 38.177 | 38.178 |
| 100 | 37.252 | 37.253 |

The gradient droops toward the output end (E_b grows with z); the energy
droop is linear in current and matches the closed-form CG result. For a
low-current FEL macropulse this is a few-percent correction; it becomes
the dominant gradient effect only at high average current. The per-cell
vg profile from the iris taper sets where the sag concentrates.