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
undulator entrance. Results are reported as **RMS** $= \sqrt{\text{MSE}}$,
which has the same units as the Twiss residuals and is the standard
figure of merit in beam physics matching.

**Quality thresholds:**

| RMS | MSE | Rating |
|-----|-----|--------|
| < $3.2 \times 10^{-2}$ | < $10^{-3}$ | Excellent |
| < $10^{-1}$ | < $10^{-2}$ | Acceptable |
| < $3.2 \times 10^{-1}$ | < $10^{-1}$ | Marginal |
| ≥ $3.2 \times 10^{-1}$ | ≥ $10^{-1}$ | Failed |

The final stage (Stage 11) jointly optimizes 4 variables: the chromaticity
section quadrupole (q5) plus the final triplet (q21, q22, q23). A
multi-start variant retries with 5 random initial conditions if the
single-start result exceeds the threshold.

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
$h = 5 \times 10^9$ /s achieve RMS $< 5.5 \times 10^{-3}$.

### W2 — Emittance Scan with Multi-Start

Re-ran emittance scan ($\varepsilon_n = 1$–$20$, 20 points) with multi-start
fallback and CMA-ES polishing.  Current results: 14 Excellent, 3 Acceptable,
3 Marginal, 0 Failed.  Both original Failed points ($\varepsilon_n = 2, 5$)
eliminated.

### W4 — COSY INFINITY Cross-Validation

COSY's gradient-based FIT reproduces the 11-stage optimization.
FR 0 (hard-edge): RMS $= 6.7 \times 10^{-5}$.
FR 1 (1st-order fringe, warm-started from FR 0): RMS $= 2.8 \times 10^{-4}$.
Cold-starting FR 1 fails (RMS ~ 0.45) due to changed edge kicks creating
incompatible local minima. Stage 5 consistently finds negative-polarity
currents (valid solution inaccessible to FELsim's bounded NM).

- Script: `UHM_beamline_opt_cosy.py`
- Report: `reports/2026/Mar/04/R2_unified_comparison_report.pdf` (§4–5)

### W6/W7 — Glyfada Optimizer Benchmarks

W6 benchmarked glyfada (ULS algorithm, 600 evals, uniform random init)
against NM for Stage 11. NM outperformed by 3–6 orders of magnitude at
$\varepsilon_n = 5, 8, 14$. W7 re-benchmarked with CMA-ES, warm-starting,
tight bounds (±3 A), and feasibility-rules constraint handling. CMA-ES
still failed at all points (RMS $= 878$ at $\varepsilon_n = 5$,
RMS $= 4.2$ at $\varepsilon_n = 8$).

The objective landscape has an extremely narrow feasibility basin;
evolutionary search cannot navigate it efficiently.

- Scripts: `UHM_beamline_opt_05ps_params.py --w6` / `--w7`
- Results: `results/params_05ps/W6/`, `results/params_05ps/W7/`
- Report: `reports/2026/Mar/04/R2_unified_comparison_report.pdf` (§7)

### W8 — RF-Track Stage 11 Optimization

Hybrid FELsim/RF-Track: stages 1–10 use FELsim, stage 11 uses RF-Track
particle tracking with prefix caching. At $\varepsilon_n = 8$:
RFT-opt RMS $= 8.4 \times 10^{-2}$ (limited by $\beta_y$ residual from
missing triangle-rule fringe correction). At $\varepsilon_n = 5$:
RFT-opt RMS $= 5.1 \times 10^{-4}$, $110\times$ better than FELsim.

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

### S5 — 2D Coupled Parameter Scans

Three 10×10 grids exploring pairwise coupling between beam parameters,
500 particles per point:

| Scan | Parameters | Fixed | Excellent | Acceptable | Marginal | Failed |
|------|-----------|-------|-----------|------------|----------|--------|
| S5a | $\sigma_E \times h$ | $\varepsilon_n = 8$ | 92 | 11 | 2 | 0 |
| S5b | $\sigma_E \times \varepsilon_n$ | $h = 5 \times 10^9$ | 38 | 27 | 17 | 18 |
| S5c | $h \times \varepsilon_n$ | $\sigma_E = 0.5\%$ | 54 | 12 | 14 | 20 |

**Key finding:** The feasibility boundary is dominated by emittance —
$\sigma_E$ and $h$ have negligible impact on matchability at the baseline
$\varepsilon_n = 8$.  Failures cluster at $\varepsilon_n < 5$ and
$\varepsilon_n > 17$, consistent with S4 and W2.

- Script: `S5_2d_parameter_scans.py`, `S5_analysis.py`
- Results: `results/params_05ps_2d/`, `results/S5_analysis/`

### S6 — Bunch Length Sensitivity

Swept bunch length from 0.1 to 2.0 ps (15 points) at two parameter sets:

| Config | $\sigma_E$ | $h$ | Excellent | Acceptable | Marginal |
|--------|-----------|-----|-----------|------------|----------|
| Baseline | 0.5% | $5 \times 10^9$ | 15/15 | — | — |
| Emittance-conserved | 2% | $20 \times 10^9$ | 11/15 | 3 | 1 |

**Result:** Transverse Twiss matching is independent of bunch length.
FELsim's linear transfer matrices decouple transverse and longitudinal
planes — the 4×4 transverse block and dispersion column do not depend on
the column-5 (time-of-flight) distribution.  This confirms S9's analytical
prediction.  The small dips in the emittance-conserved config are NM noise
(not bunch-length-correlated).

- Script: `S6_bunch_length_sensitivity.py`
- Results: `results/S6/`

### S7 — Verification Runs at Key Transition Points

Re-ran S4 transition points with $N = 500, 1000, 2000$ particles to test
statistical robustness of the 500-particle results.

| Sweep | Points | Consistent? | Notes |
|-------|--------|-------------|-------|
| Energy spread | $\sigma_E = 0.4, 0.55, 0.7\%$ | All consistent | All Excellent at every $N$ |
| Emittance | $\varepsilon_n = 1, 3, 5, 8, 14, 16, 20$ | Only $\varepsilon_n = 8$ | Quality varies with $N$ at extremes |

**Key finding:** S4/S5 emittance results at extreme values are **not
statistically robust** — they depend on the specific random particle
realization, not solely on physics.  S8 later showed that Stage 11 is
unimodal, so the variability originates from stages 1–10 (upstream beam
sampling), not from multiple local minima in the final stage.

- Script: `S7_verification_runs.py`
- Results: `results/S7/`

### S8 — Multi-Start Robustness (Stage 11)

10 random Stage 11 starting points at 5 extreme emittance values.

| $\varepsilon_n$ | RMS | Quality | $\beta_y$ (m) | 10/10 identical? |
|-----------------|-----|---------|---------------|-----------------|
| 1 | $1.08 \times 10^{-1}$ | Marginal | $\approx 0$ | Yes |
| 3 | $138.8$ | Failed | 264 | Yes |
| 5 | $5.58 \times 10^{-3}$ | Excellent | 0.242 | Yes |
| 14 | $9.53 \times 10^{-2}$ | Acceptable | 0.030 | Yes |
| 16 | $1.06 \times 10^{-1}$ | Marginal | 0.005 | Yes |

**Key finding:** Stage 11's objective landscape is **unimodal** — all 10 random
starts converge to exactly the same solution (identical MSE, quad currents, and
Twiss to machine precision).  The S7 seed-dependence at extreme emittances
originates from stages 1–10 (random beam generation), not Stage 11's search.
Multi-start Stage 11 cannot improve results.

- Script: `S8_multistart_robustness.py`
- Results: `results/S8/`

### O1 — Warm-Starting from Neighbors (Negative Result)

Sequential warm-starting: pass previous scan point's optimized currents as
initial guess for all 11 stages.

**Result:** Warm start wins only 4/20 emittance points.  The $\varepsilon_n = 1$
basin (Marginal) traps the optimizer — once locked in, the Marginal solution
propagates through $\varepsilon_n = 2$–$13$.  Cold start independently finds
better basins at each point.

- Script: `O1_warm_start_validation.py`
- Results: `results/O1/`

### R1 — Interactive Parameter Explorer

Plotly/Dash dashboard for S4 1D scans and S5 2D heatmaps.  Displays RMS with
IQR-based robust y-axis limits.

- Script: `R1_parameter_explorer.py`
- Run: `python R1_parameter_explorer.py` → http://localhost:8050

### P8 — Order-by-Order DA Convergence

COSY INFINITY DA orders 1, 2, 3, 5 with fixed quad currents at FR=0 and FR=3.

**Hard-edge (FR=0):** Linear transfer map and Twiss parameters are rigorously
identical across all DA orders.  FELsim's first-order model is exact.
Key chromatic aberrations: $T_{116} = 50.6$, $T_{126} = 624$.
RMS $= 4.8 \times 10^{-4}$ at all orders.

**3rd-order fringe (FR=3):** Linear map elements change by up to 0.3\% with DA
order (fringe-field nonlinearities feed into linear map).  Twiss match best at
$O = 3$ (RMS $= 5.9 \times 10^{-5}$) where currents were optimized; $O = 1$
gives RMS $= 6.4 \times 10^{-4}$ (still Excellent).  Large geometric
aberration $U_{1111} = 21577$.

**Finding:** First-order optics is sufficient for Twiss matching.  Higher-order
DA terms do not affect the COSY adapter's particle tracking (first-order
element maps used regardless of DA computation order).

- Script: `P8_order_convergence.py`
- Results: `results/P8/`

### P9 — Chromaticity Analysis

Swept energy deviation $\delta$ from $-3\%$ to $+3\%$ using near-mono-energetic
beams ($\sigma_E = 0.05\%$) with chromatic transport (per-particle
momentum-dependent matrices).

At $\delta = 0$: chromatic transport matches achromatic (RMS $\approx 0.09$).
At $|\delta| = 0.5\%$: RMS jumps to $1.6$–$5.4$ (Failed).
Chromaticity: $\mathrm{d}\beta_x/\mathrm{d}\delta \approx -0.9$ m/\%,
$\mathrm{d}\beta_y/\mathrm{d}\delta \approx -1.0$ m/\%.

**Finding:** Acceptance bandwidth is $|\delta| < \sim 0.3\%$ for Acceptable
Twiss matching.  The transport line is highly chromatic; beam energy stability
must be controlled to $\sim 0.3\%$.

- Script: `P9_chromaticity_analysis.py`
- Results: `results/P9/`

### P10 — Emittance Preservation

Tracked $\varepsilon_n(s)$ element-by-element through the transport line.

**Achromatic transport:** Dispersion-corrected $\varepsilon_n$ conserved to
machine precision ($\Delta = 0.00\%$) in both planes.  Raw
$\varepsilon_{n,x}$ grows $16\%$ through the chicane (x-$\delta$ coupling),
raw $\varepsilon_{n,y}$ is perfectly conserved.

**Chromatic transport:** Dispersion-corrected $\varepsilon_{n,x}$ grows by
$591\%$ — emittance is **not conserved** due to chromatic filamentation
(energy-dependent optics creates irreversible phase space distortion in the
$\sigma$-matrix sense).

- Script: `P10_emittance_evolution.py`
- Results: `results/P10/`

### P11 — Fringe Field Impact

Quantified the effect of DPW triangle-model fringe correction ($\phi$) on
Twiss matching within FELsim's first-order model.

**FELsim fringe architecture:**
- DPW edge kick always includes $\phi = (l_e/6) \cdot h \cdot (1 + \sin^2\eta) / \cos\eta$
- `fringeType` parameter affects field profile only (drift-space matrix)
- Quadrupole fringe: not modeled

**Results with fixed (warm-start) currents:**

| Mode | RMS | $\beta_y$ (m) | $\alpha_y$ |
|------|-----|---------------|------------|
| With $\phi$ | $7.8 \times 10^{-2}$ | 0.232 | 0.013 |
| Without $\phi$ | $5.7 \times 10^{-1}$ | 0.241 | $-1.12$ |

The $\phi$ correction modifies $M_{43}$ by 2–8% (transport dipoles) and
5% (FC1 chicane).  Removing it degrades RMS by $7\times$, primarily
through $\alpha_y$.  Zero-angle DPW faces acquire "pure" fringe kicks
(0.046 m$^{-1}$) that have no hard-edge equivalent.

- Script: `P11_fringe_field_impact.py`
- Results: `results/P11/`

### C4 — CI Pipeline Expansion

Expanded the GitHub Actions test suite from 8 to 16 test files.  Added
`test_chromatic_dipole`, `test_chromatic_quad`, `test_optimizer`,
`test_edge_cases`, `test_rftrack`, `test_felsim_unified`,
`test_cosy_unified`, and `test_rftrack_unified`.  Tests requiring external
tools (COSY, RF-Track) are gracefully skipped via `pytest.skip()` at module
level.  Visual tests are excluded with `-m "not visual"`.

- Config: `.github/workflows/tests.yml`
- Markers: `conftest.py` (`visual`, `cosy`, `rftrack`)

### Bug Fixes — AlgebraicOptimization Variable Ordering

Fixed non-deterministic `set.pop()` variable extraction in
`AlgebraicOptimization.py` — replaced with `sorted()` by symbol name to
ensure reproducible x/y assignment across Python versions and runs.

### Visualization Decoupling

Added `--plots-only` flag to P9, P10, and P11 study scripts.  Regenerates
figures from cached `summary.json` without re-running computation.  P10
updated to save full element-by-element evolution data (was 5 sample
points).  S4/S5/R2 already had this capability.

## Planned Work

See `backend/test/PRIORITIES.md` for the full roadmap.
