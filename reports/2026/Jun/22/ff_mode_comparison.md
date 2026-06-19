# FF-mode comparison - UH MkV FEL transport line

Author: Eremey Valetov. FF = fringe field. Undulator target beta_y = 0.242 m, beta_x = 1.40 m.

## Table A - COSY fringe-order convergence (computed 2026-06-19)

Each fringe order is warm-started to the undulator match. All reach the target with a different current set; the last column is the RMS quad-current shift relative to FR0.

| COSY mode | beta_x (m) | beta_y (m) | match MSE | RMS |dI| vs FR0 (A) |
|---|---:|---:|---:|---:|
| FR0 (hard-edge) | 1.400 | 0.241 | 2.27e-07 | 0.000 |
| FR1 (Enge 1st, warm) | 1.400 | 0.242 | 3.65e-09 | 0.193 |
| FR2 (Enge 2nd, warm) | 1.400 | 0.242 | 7.49e-09 | 0.654 |
| FR3 (RK fringe, warm) | 1.400 | 0.242 | 4.39e-09 | 1.084 |

## Cross-code fringe status

A cross-CODE fringe comparison must re-optimise each code to the target (applying one code's currents to another does not match: R2 established different currents -> same Twiss). Established results:

- FELsim (triangle-rule DPW edge) and COSY (FR0/FR1) agree on matched Twiss (R2: MSE 1.2e-6 / 2.3e-7 / 3.7e-9 at eps_n=8).
- RF-Track: the R2-era vertical beta_y deficit (0.055 vs 0.242 m) came from a missing triangle-rule edge phi. The triangle-phi correction was since added (rftrackAdapter._annotate_dipole_edges); RF-Track now re-optimises to the target (RFT-opt MSE 7.0e-3, eps_n=8), i.e. the deficit is closed.
- xsuite: no dipole edge/fringe model (drift edges) -> the no-fringe baseline. Its bend body now builds after the xsuiteAdapter h->angle fix (2026-06-19); quantifying its matched beta_y needs a re-optimisation run (remaining quick task).

REMAINING (post-Monday): re-optimise RF-Track and xsuite to the undulator target and tabulate achieved beta_y alongside COSY FR0-FR3 and FELsim for a single 5-code FF figure.
