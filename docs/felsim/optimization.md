# Quadrupole Optimization

The transport line has 23 quadrupoles whose currents must be optimized so
that the Twiss parameters at the undulator entrance match the targets from
arXiv:2510.14061v1 Table I.

## Optimization Strategy

The optimizer (`beamOptimizer.py`) supports two methods:

- **Nelder-Mead** (default): 11-stage sequential approach. Each stage
  optimizes a small group of quadrupoles while holding the others fixed,
  progressively refining the match from upstream to downstream.
- **Glyfada** (`method='glyfada'`): evolutionary optimizer
  (`glyfadaAdapter.py`) using the glyfada C++ MPI binary. Uses
  population-based search to explore multiple basins simultaneously,
  useful for high-dimensional problems or when Nelder-Mead gets trapped in
  local minima.

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

### W4 — COSY INFINITY Cross-Validation

COSY's gradient-based FIT reproduces the 11-stage optimization.
FR 0 (hard-edge): MSE = $4.5 \times 10^{-9}$.
FR 1 (1st-order fringe, warm-started from FR 0): MSE = $7.9 \times 10^{-8}$.
Cold-starting FR 1 fails (MSE ~ 0.2) due to changed edge kicks creating
incompatible local minima. Stage 5 consistently finds negative-polarity
currents (valid solution inaccessible to FELsim's bounded NM).

- Script: `UHM_beamline_opt_cosy.py`
- Report: `reports/2026/Mar/04/R2_unified_comparison_report.pdf` (§4–5)

### W6/W7 — Glyfada Optimizer Benchmarks

W6 benchmarked glyfada (ULS algorithm, 600 evals, uniform random init)
against NM for Stage 11. NM outperformed by 3–6 orders of magnitude at
$\varepsilon_n = 5, 8, 14$. W7 re-benchmarked with CMA-ES, warm-starting,
tight bounds (±3 A), and feasibility-rules constraint handling. CMA-ES
still failed at all points (MSE $7.7 \times 10^5$ at $\varepsilon_n = 5$,
MSE $17.8$ at $\varepsilon_n = 8$).

The FELsim MSE landscape has an extremely narrow feasibility basin;
evolutionary search cannot navigate it efficiently.

- Scripts: `UHM_beamline_opt_05ps_params.py --w6` / `--w7`
- Results: `results/params_05ps/W6/`, `results/params_05ps/W7/`
- Report: `reports/2026/Mar/04/R2_unified_comparison_report.pdf` (§7)

### W8 — RF-Track Stage 11 Optimization

Hybrid FELsim/RF-Track: stages 1–10 use FELsim, stage 11 uses RF-Track
particle tracking with prefix caching. At $\varepsilon_n = 8$:
RFT-opt MSE = $7.0 \times 10^{-3}$ (limited by $\beta_y$ residual from
missing triangle-rule fringe correction). At $\varepsilon_n = 5$:
RFT-opt MSE = $2.6 \times 10^{-7}$, $110\times$ better than FELsim.

- Script: `UHM_rftrack_opt.py`
- Results: `results/rftrack_opt/`
- Report: `reports/2026/Mar/04/R2_unified_comparison_report.pdf` (§4)

### W12 — Bunch Compression Feasibility

Can the transport line compress 2 ps bunches to 0.5 ps?  Chirp sweep
(analytical + COSY map propagation), RF-Track validation (post-C7 fix),
extended current bounds (15 A), $T_{566}$ assessment.

- $R_{56} = 27.09$ mm (geometry-locked), $T_{566} = 0$
- Compression floor $\approx 0.45$ ps ($R_{56} \times \sigma_\delta$)
- RF-Track Part B: at $C = 4$ chirp, RF-Track gives $\sigma_z = 1.94$ ps
  vs COSY map 0.67 ps (aperture losses + nonlinear effects)
- Transport line is not a compressor; compression should occur upstream
- Script: `W12_compression_feasibility.py`
- Results: `results/W12/`
- Report: `reports/2026/Mar/04/R3_longitudinal_report.pdf` (§5)

### R2 — Cross-Code Validation Report

Unified report documenting 3-code agreement (FELsim, COSY, RF-Track),
the 9 bug fixes enabling it, quad current comparison, parameter
sensitivity, and optimizer benchmarks.

- Report: `reports/2026/Mar/04/R2_unified_comparison_report.pdf`

### R3 — Longitudinal & Compression Report

Combined report merging W9 (0.5 ps / 2 ps transfer map, bunch
propagation) and W12 (compression feasibility with RF-Track validation).

- Report: `reports/2026/Mar/04/R3_longitudinal_report.pdf`

## Planned Work

See `backend/test/PRIORITIES.md` for the full roadmap. Key upcoming items:

- **S5** — 2D coupled parameter scans ($\sigma_E$ vs $h$, $\sigma_E$ vs
  $\varepsilon_n$, $h$ vs $\varepsilon_n$) — smoke test done, full scans pending
