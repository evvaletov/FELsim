# Quadrupole Optimization

The transport line has 23 quadrupoles whose currents must be optimized so
that the Twiss parameters at the undulator entrance match the targets from
arXiv:2510.14061v1 Table I.

## Optimization Strategy

The optimizer uses an 11-stage sequential Nelder-Mead approach
(`beamOptimizer.py`). Each stage optimizes a small group of quadrupoles
while holding the others fixed, progressively refining the match from
upstream to downstream.

**Objective function:** mean squared error (MSE) between computed and target
Twiss parameters ($\beta_x$, $\alpha_x$, $\beta_y$, $\alpha_y$) at the
undulator entrance.

**Quality thresholds:**

| MSE | Rating |
|-----|--------|
| < $10^{-3}$ | Excellent |
| < $10^{-2}$ | Acceptable |
| < $10^{-1}$ | Marginal |
| ≥ $10^{-1}$ | Failed |

The final stage (Stage 11) jointly optimizes 4 variables: the chromaticity
section quadrupole (q5) plus the final triplet (q21, q22, q23). A
multi-start variant retries with 5 random initial conditions if the
single-start result exceeds the MSE threshold.

## Completed Studies

### S1 — 2 ps Baseline

Paper-aligned asymmetric Twiss targets ($\beta_x = 1.4$ m,
$\alpha_x = 0.47$, $\beta_y = 0.24$ m, $\alpha_y = 0$).
Joint 4-variable final stage.

- Script: `UHM_beamline_opt_v2.py`
- Report: `UHM_beamline_opt_v2_report.md`

### S2 — 0.5 ps Longitudinal Emittance Conservation

$\sigma_E = 2\%$, $h = 20 \times 10^9$ /s from longitudinal emittance
conservation. Symmetric Twiss targets ($\beta = 0.24$ m, $\alpha = 0$).

- Script: `UHM_beamline_opt_05ps.py`

### S3 — 0.5 ps Paper-Aligned

$\sigma_E = 0.5\%$, $h = 5 \times 10^9$ /s unchanged from baseline.
Asymmetric Twiss targets from Table I.

- Script: `UHM_beamline_opt_05ps_v2.py`

### S4 — 0.5 ps Parameter Sensitivity

Three 1D sweeps: energy spread (0.1–5%), chirp ($0$–$40 \times 10^9$ /s),
emittance (1–20 π·mm·mrad). 500 particles per point.

- Script: `UHM_beamline_opt_05ps_params.py`
- Results: `results/params_05ps/`

### W1 — Chirp Effect on Twiss Matching

Confirmed chirp has negligible effect on Twiss matching. Both $h = 0$ and
$h = 5 \times 10^9$ /s achieve MSE < $3 \times 10^{-5}$.

### W2 — Emittance Scan with Multi-Start

Re-ran emittance scan ($\varepsilon_n = 1$–$20$, 20 points) with multi-start
fallback. Resolved dips at $\varepsilon_n = 14$–$16$; all now excellent
(MSE ~ $10^{-6}$–$10^{-9}$).

## Planned Work

See `backend/test/PRIORITIES.md` for the full roadmap. Key upcoming items:

- **S5** — 2D coupled parameter scans ($\sigma_E$ vs $h$, $\sigma_E$ vs
  $\varepsilon_n$, $h$ vs $\varepsilon_n$)
- **W4** — COSY INFINITY cross-validation of the Python results
- **C1** — RF-Track cross-validation at key parameter points
- **O3** — Evolutionary optimization for high-dimensional problems
