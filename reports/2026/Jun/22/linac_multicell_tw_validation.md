# Cell-resolved TW linac model - validation

Author: Eremey Valetov, 2026-06-19. Built from the production tau=0.57
cell geometry (SLAC-75 Table 6-6, 86 cells).

| Quantity | This model | elegant ref | agreement |
|---|---:|---:|---:|
| Peak K_out (MeV) | 40.953 | 41.442 | 1.18% |
| Optimal phase (deg) | 339 | 70 | conv. offset |
| Adiabatic damping det(Rx)=p_in/p_out | 0.0343 | 0.0359 | 4.5% |

L_cell = 34.99 mm (2pi/3 at 2.856 GHz), L_total = 3.0091 m, E0 = 13.3 MV/m (constant gradient).
Iris taper 2a 1.032 -> 0.757 in; estimated vg/c 0.0204 -> 0.0065, fill time 0.81 us.

The autophasing integration reproduces the reference peak energy gain,
confirming the beta<1 phase slippage is handled (a naive beta=1 cavity
chain collapses to ~7 MeV). The residual 1.2% is the structure-length effect: 86 cells x 34.99 mm = 3.009 m vs the 3.048 m nominal (the known fractional-cell tail),
not a model error. The optimal-phase offset vs elegant is a convention
difference. Next: per-cell gradient sag from beam loading (uses the vg
profile above), and wiring as an xsuite/xtrack element with an
energy-ramp reference.