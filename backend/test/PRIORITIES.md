# UH MkV FEL Beamline Optimization — Priorities & Roadmap

**Date:** 2026-02-11 (updated 2026-03-03)
**Scripts:** `backend/test/UHM_beamline_opt_*.py`

---

## Due by 2026-02-18 (Wednesday)

### W1. Table I Optimisations Without Chirp [DONE 2026-02-12]
- Ran 0.5 ps and 2 ps Table I studies with h=0 and h=5e9 (1000 particles each).
- Confirmed: chirp has negligible effect on Twiss matching.
  Both configurations achieve MSE < 3e-5.
  h=0 actually gives slightly lower MSE (1.27e-6 vs 2.6e-5).
- Results: `results/params_05ps/w1_chirp_comparison.csv`
- Implementation: `--w1` flag in `UHM_beamline_opt_05ps_params.py`

### W2. Emittance Scan — Multi-Start + Raised Bounds [DONE 2026-02-12]
- Re-ran emittance scan (ε_n = 1–20, 20 points) with fallback multi-start:
  - Default 10 A bounds first; if MSE > 1e-3, retry with 5 Stage-11 restarts
    and optionally raised chromaticity bounds (15 A)
- **Primary goal achieved:** ε_n=14–16 dips resolved — all now Excellent
  (MSE ≈ 10⁻⁶–10⁻⁹) via multi-start Stage 11
- Remaining issues (candidates for evolutionary optimizer):
  - ε_n=2: Failed (physics limit — beam too small to match)
  - ε_n=5: Failed (Nelder-Mead local minimum trap)
  - ε_n=17: Acceptable (q10 needed 15 A)
- Results: `results/params_05ps/scan_emittance_w2.csv`,
  `mse_vs_emittance_w2_comparison.eps`
- Implementation: `--w2` flag in `UHM_beamline_opt_05ps_params.py`

### W3. Confirm ε_n = 8 μm Baseline [DONE 2026-02-11]
- All scripts confirmed using `epsilon_n = 8` with comment `# pi.mm.mrad`.
  This is 8 π·mm·mrad normalised emittance — matches the expected photocathode
  emittance for the UH MkV FEL.

### W5. Chirp Value Assessment [DONE 2026-02-12]
- **Assessment:** h = 5×10⁹ /s is **low-to-moderate** for an S-band FEL at 40 MeV.
  Corresponds to ~16° off-crest. Adds 0.25%/σ_t to energy spread (total σ_E
  grows from 0.50% to 0.56%, +12%). Comparable to SDL/BNL residual chirp.
- **h = 20×10⁹ /s** (emittance-conserved) is moderate-to-high. Adds 1.0%/σ_t,
  more than doubling total σ_E to 1.12%. Comparable to SPARC high-chirp.
- **Provenance:** h = 5e9 originates from Niels Bidault's initial beam parameters
  (Aug 7 2025 email: "to be assessed"; Feb 4 2026 email: "remains to be determined
  by injector simulations"). Heritage from legacy FELsim `goldTwiss.py`.
  Not from arXiv:2510.14061v1. Not hardcoded in FELsim core (defaults to h = 0).
- **Plots:** `results/w5_chirp/chirp_phase_space_comparison.eps`,
  `chirp_energy_spread_growth.eps`, `chirp_context_comparison.eps`
- **Script:** `w5_chirp_assessment.py`

### W6. Glyfada vs Nelder-Mead Benchmark [DONE 2026-02-23]
- Benchmarked glyfada evolutionary optimizer against NM for Stage 11 at
  ε_n = 5, 8, 14.
- **Result:** Glyfada fails at all three points. NM outperforms by 3–6
  orders of magnitude in MSE. At ε_n=5, all 600 glyfada evaluations hit
  the penalty value (unstable optics).
- **Root cause:** Uniform random initialization over wide bounds wastes
  evaluations in infeasible regions; 600 evals insufficient for a 4D
  landscape with large singular regions.
- **Conclusion:** DH evaluator protocol not suitable for FELsim's
  fast-evaluation Stage 11. An in-process global optimizer
  (`scipy.optimize.differential_evolution`) with warm-starting is the
  recommended path forward.
- Report: `W6_glyfada_benchmark_report.tex`
- Results: `results/params_05ps/W6/benchmark_results.csv`,
  `mse_comparison.eps`, `time_comparison.eps`
- Implementation: `--w6` flag in `UHM_beamline_opt_05ps_params.py`

### W7. Glyfada Config Optimization & Re-Benchmark [IN PROGRESS]
- Re-benchmarks Glyfada against NM with optimized configurations leveraging
  CMA-ES, constraint handling (feasibility_rules), warm-starting from NM
  solution, and tighter bounds (±3A around NM result).
- **Config A (CMA-ES):** pop_size=20, max_gen=150, initial_sigma=0.3,
  feasibility_rules constraint handling. 3000 evals.
- **Config B (NSGA-II Phased):** pop_size=30, max_gen=100, two phases
  (σ=0.05 → σ=0.01 at gen 50). 3000 evals.
- **Key improvements over W6:** warm-starting, CMA-ES algorithm, constraint
  handling replaces blunt 1e6 penalty, informed bounds, 5× more evals.
- glyfadaAdapter updated with extra_config pass-through for new Glyfada
  features (CMA-ES, phased optimization, constraints).
- glyfada_eval updated with stability constraint output.
- Results: `results/params_05ps/W7/`
- Implementation: `--w7` flag in `UHM_beamline_opt_05ps_params.py`

### W9. COSY Longitudinal Study [DONE 2026-02-25]
- Full 3D (6D phase space) COSY simulation with longitudinal diagnostics
- Extracts R56, T566, coupling terms from optimised beamline
- Propagates 6D bunches for 0.5 ps and 2 ps modes
- Tests adding R56=0 as an optimisation objective
- Reviewer verdict on what changes when switching bunch length
- Script: `W9_cosy_longitudinal_study.py`
- Results: `results/W9/`

### I6. COSY σ_z Blowup Investigation [DONE 2026-03-02]
- **Root cause found:** The 60–100× blowup was NOT reproduced with W9-optimized
  currents (MSE ≈ 2.3e-7). With properly matched optics:
  - ORDER 1: σ_z = 2.264 ps (1.13× growth, consistent with R56 × σ_δ)
  - ORDER 2: σ_z = 3.740 ps (1.87× growth)
  - ORDER 3: σ_z = 4.158 ps (2.08× growth)
- **Conclusion:** The original blowup occurred with poorly matched currents
  (MSE ≈ 7.1e-3). Higher-order growth (ORDER 3 > ORDER 1) is correct physics
  from nonlinear path-length terms, not a bug.
- **Element-by-element evolution:** σ_l grows smoothly through the beamline with
  a sharp jump at the chicane (elements 60–75) due to R56 coupling.
- **Single particle test (Part C):** COSY RRAY with N=1 loses the particle at
  every checkpoint (known RRAY limitation with single data particles).
- Script: `I6_sigma_z_diagnostic.py`
- Results: `results/I6/`

### C5. RF-Track SBend Bug & Analytical Workaround [DONE]
- **Root cause:** RF-Track v2.5.5 SBend body tracking interprets Bunch6d's
  absolute momentum P [MeV/c] (6th column) as momentum deviation δ = ΔP/P₀.
  Constructor `SBend(L, angle, P_Q)` produces ~910 mm displacement. Setter-only
  `set_K0` gives identity (= Drift). `RBend` has the same issue.
- **Workaround implemented (2026-03-02):**
  1. **Analytical sector-bend correction (DPH):** Track dipole body as Drift
     (preserves y-plane and path length), then apply M_correction = M_sector × M_drift⁻¹
     to (x, x') plus dispersion (R₁₆, R₂₆) and R₅₆ corrections. Exact to 6 decimal
     places vs analytical sector-bend matrix. Implemented in `_apply_sector_bend_correction()`,
     `_track_segmented()`, `track_elements()`, and `collect_evolution()`.
  2. **Edge kicks (DPW):** Thin-lens `Quadrupole(L=1e-10)` with
     `K1L = -|K0| * tan(wedge_angle)`. Uses unsigned |K0| to match FELsim's
     `R = L/|θ|` convention — signed K0 inverts chicane edge kicks.
  3. **DPW-DPH-DPW triplet detection:** `_annotate_dipole_edges()` scans
     the beamline and writes `dipole_K0` into DPW parameters.
  4. **Corrector normalization:** `Kx = BdL/P₀` (same as Quadrupole set_strength).
- **Key implementation bugs found and fixed:**
  - R-matrix unit conversion: R12/R21 do NOT need ×1000/÷1000 scaling in (mm,mrad)
    coordinates — both x and x' scale by same factor, so R-matrix is invariant.
  - Edge kick sign: must use `abs(K0)` for DPW, matching FELsim's `R = L/|θ|`.
    Signed K0 inverted chicane edge kicks → β_y blew up (MSE=1583).
- **Validation results (ε_n=8, 3 restarts):**
  - FELsim MSE = 6.13e-6 (baseline)
  - RFT-val (FELsim currents → RF-Track): MSE = 1.985 (model differences shift optimum)
  - **RFT-opt (RF-Track optimized): MSE = 7.0e-3** (matches FELsim quality)
  - RF-Track can now independently optimize to good Twiss match
- **Remaining discrepancy:** DPW thin-quad doesn't include triangle-model fringe
  correction φ → β_y slightly off (0.055 vs target 0.2418)
- **Bug report:** File with RF-Track maintainer. See C5-BUG below.
- **Files:** `rftrackAdapter.py`, `test/rftrack_sbend_bug_mwe.py`

### C5-BUG. RF-Track SBend Bug Report [TODO]
- **Prepare a minimal working example demonstrating that `SBend` has no
  transverse tracking effect in RF-Track v2.5.5.**
- **MWE structure:**
  1. Create SBend via constructor `SBend(L=0.2, angle=0.393, P_Q=40.5)`
  2. Track on-axis and 1mm off-axis particles (Bunch6d)
  3. Show output = input (identity), compare with expected cos θ / sin θ
  4. Show constructor gives ~910 mm displacement (P/δ confusion)
  5. Show `set_E1`/`set_E2` has no effect
  6. Show RF-Track version, platform, compilation info
- **Context:** SBend is documented in the RF-Track reference manual §4.6
  with constructor `SBend(L, angle, P_Q, E1, E2, K1)`. The API methods
  (`set_K0`, `set_h`, `set_Bfield`, etc.) all accept values silently
  but none affect the tracking output.
- **File with:** GitLab issue on https://gitlab.cern.ch/rf-track or email

### O4. Glyfada 26D Feasibility [FINDING — 2026-03-02]
- Ran glyfada on Koa with 26 quad currents and FELsim transfer-matrix
  objective. Both wide bounds [0, 10 A] and tight bounds [NM ± 2 A] tried.
- **Result:** ~27,000 evaluations, ~99% returned penalty (unstable optics).
  Best solution = NM starting point itself. The FELsim MSE landscape has an
  extremely narrow feasibility basin in 26D — evolutionary search cannot
  navigate it.
- **Auto-algorithm selected SA** (Simulated Annealing) over CMA-ES, SHADE,
  NSGA-II. Landscape characterised as rugged=1.0 with 7 modes.
- **Conclusion:** Glyfada cannot improve on FELsim NM with the same objective.
  Value of glyfada requires either: (a) RF-Track particle-tracking objective
  (install RF-Track on Koa), or (b) softer penalty function providing gradient
  info for unstable solutions.
- **Config:** `parameters_tight.json`, `parameters_wide.json` on Koa scratch

### W11. Throughput Optimization — Maximize Peak Current [DONE 2026-03-02]
- Extends Twiss-only Stage 11 optimization to include transmission and bunch
  length objectives via weighted scalar cost function
- Two scenarios: 2 ps → 2 ps (transport) and 2 ps → 0.5 ps (compression)
- RF-Track particle tracking with physical apertures for Stage 11 NM optimization
- Objective: w_t × MSE_Twiss + w_T × (1-T)² + w_σ × (σ_t/σ_target - 1)²
- **Results (post-C5 re-run with analytical dipole correction):**
  - 2ps transport: MSE=0.175, T=42.8%, σ_t=1.83 ps, I_peak=5.6 A
  - 0.5ps compress: MSE=0.071, T=28.4%, σ_t=1.31 ps (target 0.5), I_peak=5.2 A
  - Compression scenario failed — σ_t target unreachable with Stage 11 quads alone
- Script: `W11_throughput_opt.py`
- Results: `results/W11/`

### W12. Bunch Compression Feasibility Study [DONE 2026-03-02]
- Can the transport line compress 2 ps → 0.5 ps?
- Part A: Analytical compression + chirp sweep via COSY map propagation (41 points)
  - Compression floor ≈ 0.45 ps (R56 × σ_δ), C=4 chirp gives ~0.67 ps
  - T566 = 0 confirmed — no second-order effect
- Part B: RF-Track validation with C7 coord5 fix (re-run of W10 Part B scenarios)
- Part C: Extended bounds (15 A, 5 restarts) — improves transverse MSE, σ_z unchanged
- Part D: Feasibility summary — transport line is not a compressor
- **Conclusion:** compression should occur upstream (velocity bunching / dedicated compressor)
- Script: `W12_compression_feasibility.py`
- Report: `W12_compression_feasibility_report.tex`
- Results: `results/W12/`

### W10. Beam Losses & Bunch Compression Study [DONE 2026-03-02]
- Quantifies particle losses through the full transport line with physical apertures
- Part A: Transmission baseline at 2 ps and 0.5 ps (COSY + RF-Track)
- Part B: Bunch compression via negative chirp — demonstrates chirp required,
  energy spread alone does NOT compress (R56 ≈ +27 mm elongates unchirped beams)
- **NOTE:** W10 Part B RF-Track results are invalidated by the C7 coord5 bug
  (pass-through inflated initial σ(ct) by 9.5×). See W12 Part B for corrected results.
- Part C: Charge scan (20–300 pC) at both operating modes, RF-Track SC ≥ 100 pC
- RF-Track adapter extended with per-element physical apertures
  (`enable_physical_apertures()`)
- Apertures: quad bore 27 mm, dipole gap 14.5 mm, dipole width 50 mm (placeholder)
- Script: `W10_beam_losses_compression.py`
- Results: `results/W10/`

### W4. COSY INFINITY Cross-Validation [DONE 2026-02-15]
- COSY's internal FIT reproduces the 11-stage optimisation in a single run.
- **FR 0 (no fringe):** MSE = 4.5e-9 — converges from default starting points.
- **FR 1 (1st-order fringe):** MSE = 7.9e-8 — requires warm-start from FR 0;
  cold-start fails (MSE ~0.2) due to local minima from changed edge kicks.
- **FR 3 (3rd-order fringe):** Too slow — each lattice evaluation takes >15 min
  due to Runge-Kutta fringe field integration. Single-dipole comparison available.
- **Fieldmap (MGE):** Blocked — DELTAS=0 and field values ~50× too small (see C3).
- **Key finding:** Warm-starting is essential for fringe field optimisation.
  The sequential 11-stage approach is sensitive to the dipole model; different
  fringe settings change edge kicks enough to create incompatible local minima.
- Stage 5 consistently converges to negative-polarity currents in COSY (valid
  solution, inaccessible to FELsim's bounded Nelder-Mead).
- S1 and S3 produce identical transverse optimisation results.
- Report: `W4_cosy_xval_report.tex` (7 pages)
- Results: `results/cosy_s1_fr0.json`, `results/cosy_s1_fr1_warm.json`
- Script: `UHM_beamline_opt_cosy.py` (supports `--fr`, `--order`, `--warm-start`)

---

## Completed Studies

### S1. 2 ps Baseline Optimization [DONE]
- Script: `UHM_beamline_opt_v2.py`, report: `UHM_beamline_opt_v2_report.md`
- Paper-aligned asymmetric Twiss targets (β_x=1.4 m, α_x=0.47)
- Joint 4-variable final stage (chromaticity 5 + triplet)

### S2. 0.5 ps Fixed Parameters — Longitudinal Emittance Conservation [DONE]
- Script: `UHM_beamline_opt_05ps.py`, report: `UHM_beamline_opt_05ps_report.md`
- σ_E=2%, h=20e9 from longitudinal emittance conservation
- Symmetric Twiss targets (β=0.24 m, α=0 both planes)

### S3. 0.5 ps Fixed Parameters — Paper-Aligned [DONE]
- Script: `UHM_beamline_opt_05ps_v2.py`, report: `UHM_beamline_opt_05ps_v2_report.md`
- σ_E=0.5%, h=5e9 unchanged from baseline per arXiv:2510.14061v1 §III
- Asymmetric Twiss targets from Table I
- Exploration: `UHM_beamline_opt_05ps_v2_explore.py`

### S4. 0.5 ps 1D Parameter Sensitivity [DONE 2026-02-11]
- Script: `UHM_beamline_opt_05ps_params.py`
- Report: `UHM_beamline_opt_05ps_params_report.md`
- Results: `results/params_05ps/`
- Three 1D sweeps at 0.5 ps: energy spread (0.1–5%), chirp (0–40e9), emittance (1–20)
- MSE quality thresholds: <1e-3 excellent, <0.01 acceptable, <0.1 marginal

---

## Category S: Parameter Studies

### S5. 0.5 ps 2D Coupled Scans [IN PROGRESS — script done, full scans pending]
- **Motivation:** 1D scans hold other parameters fixed; realistic operation involves
  correlated changes (e.g., shorter bunch → larger energy spread). 2D scans map
  the feasibility surface.
- **Scans:**
  - S5a: (σ_E, h) grid at ε_n=8 — energy spread vs chirp coupling
  - S5b: (σ_E, ε_n) grid at h=5e9 — degradation interaction
  - S5c: (h, ε_n) grid — chirp compensation vs emittance
- **Design:** 10×10 grids (100 points each), 500 particles, ~5 hours per scan.
  Checkpoint/resume via CSV. Contour/heatmap plots (MSE LogNorm, Twiss deviation,
  feasibility bands). CLI: `--s5a/--s5b/--s5c/--all/--plots-only/--grid N`.
- **Status:** 3×3 smoke test (S5a) completed, all 9 points converge.
  Full 10×10 scans not yet launched.
- Script: `S5_2d_parameter_scans.py`
- **Output:** `results/params_05ps_2d/`, contour plots, feasibility boundary curves
- **Prerequisite:** S4 results to identify interesting regions

### S6. Bunch Length Sensitivity (0.1–2 ps) [MEDIUM PRIORITY]
- **Motivation:** The FELsim request asks about 0.5 ps specifically, but understanding
  the full bunch length range is valuable context.
- **Design:** Sweep bunch_spread from 0.1 to 2.0 ps (15 points) at two parameter
  sets: (a) baseline (σ_E=0.5%, h=5e9) and (b) emittance-conservation scaled.
- **Key question:** At what bunch length does the optimizer start to degrade?

### S7. Verification Runs at Key Points [MEDIUM PRIORITY]
- **Motivation:** The 500-particle sweeps trade accuracy for speed. Key points
  (boundaries, transitions) need 1000–2000 particle confirmation.
- **Design:** From S4 results, identify 5–8 key points per scan where MSE crosses
  thresholds. Re-run with 1000 and 2000 particles, compare.
- **Output:** Verification table in report, error bars on sensitivity plots

### S8. Multi-Start Robustness Study [LOW PRIORITY]
- **Motivation:** Nelder-Mead is a local optimizer. Different starting points may
  find different local minima, especially at extreme parameter values.
- **Design:** At 5 extreme parameter points, run 10 random starts per point.
  Report: best/worst/median MSE, current spread.
- **Key question:** Is the optimization landscape convex enough for single-start?

### S9. Bunch Length Independence Study [DONE 2026-02-22]
- **Motivation:** S1 (2 ps) and S3 (0.5 ps) produce identical quad currents.
  A reviewer asks why, whether this is correct, and how to switch in practice.
- **Root cause:** Transverse Twiss depends on ε, σ_δ, not σ_z; bunch length
  enters only column 4 (time-of-flight), which `cal_twiss()` ignores.
- **Part A:** Analytic estimates — R56 (FELsim + COSY vs analytic), compression
  chirp table, CSR (Derbenev–Saldin), LSC, resistive wall wakefields.
  Summary table classifying each effect at 2 ps and 0.5 ps.
- **Part B1:** Pre-compressed beam: re-run optimization at 0.5 ps with varied σ_E
  (0.5%, 2%, 3%) and h=0 to confirm transverse decoupling.
- **Part B3:** σ_E sensitivity scan (0.1–3%) at 0.5 ps, h=0. Identify threshold
  where matching degrades.
- **Part B4:** Cross-validation with COSY INFINITY (transfer maps, R56, fringe
  fields) and RF-Track (particle tracking with space charge / CSR).
- **Part C:** Report section (LaTeX): linear decoupling explanation, compression
  options for UH MkV, collective effects at 0.5 ps.
- **Scripts:** `S9_bunch_length_study.py` (Parts A + B1–B3),
  `S9_rftrack_cosy_validation.py` (Part B4)
- **Results:** `results/S9/`

---

## Category O: Optimizer Improvements

### O1. Warm-Starting from Neighboring Points [MEDIUM PRIORITY]
- **Motivation:** Each scan point starts from the same fixed initial guess. Using
  the optimized currents from the previous (neighboring) scan point as the start
  could improve convergence speed and robustness at extreme parameter values.
- **Implementation:** Pass previous `quad_currents` dict as initial guess to
  `run_optimization()`. Add `warm_start=` parameter.
- **Expected benefit:** 2–5× faster scans, better convergence at boundaries.

### O2. Adaptive Scan Resolution [LOW PRIORITY]
- **Motivation:** Uniform spacing wastes points in flat regions and undersamples
  transition regions.
- **Design:** Bisection refinement: run coarse scan (5 points), identify where
  MSE changes rapidly, insert midpoints. Iterate to target resolution.

### O3. Evolutionary Optimisation [DONE 2026-02-22]
- Glyfada adapter implemented: `backend/glyfadaAdapter.py` (GlyfadaOptimizer class)
  and `backend/glyfada_eval.py` (DH evaluator script).
- Integrated into `beamOptimizer.py` via `method='glyfada'` in `calc()`.
- Uses the DH evaluator protocol: pickles objective function, spawns
  `mpirun -np N optimiser --config=parameters.json`.
- Defaults tuned for FELsim: pop_size=50, max_gen=100, sigma=0.05, multistart mode.
- Supports `n_processes`, `pop_size`, `max_gen`, `sigma`, `algorithm` kwargs.
- **Benchmark (W6):** Glyfada underperforms NM at all tested emittance points
  (see W6). Root causes: insufficient eval budget, no warm-starting, wide bounds.
- **Re-benchmark (W7):** Optimized configs with CMA-ES, warm-starting from NM,
  ±3A bounds, constraint handling. Tests whether Glyfada can match NM with
  proper configuration.

---

## Category R: Reporting & Visualization

### R1. Interactive Parameter Explorer [LOW PRIORITY]
- **Motivation:** Static EPS plots are good for papers; interactive plots help
  with exploration.
- **Design:** Plotly/Dash dashboard reading CSV data. Sliders for parameter
  selection, hover for quad currents.

### R2. Comparison Table Across All Studies [DONE 2026-03-02]
- Aggregates data from W4, S4, W8, W9, W10, W11, W12 into 5 cross-code tables
  and 3 summary plots.
- **Tables:** (1) Baseline cross-code optimization, (2) Parameter sensitivity summary,
  (3) Bunch length & transmission, (4) Compression feasibility, (5) Quad currents.
- **Plots:** 3-panel MSE vs parameter, cross-code Twiss bar chart, compression curve.
- Script: `R2_unified_comparison.py`
- Report: `R2_unified_comparison_report.tex`
- Results: `results/R2/`

---

## Category I: Infrastructure

### I1. Propagate Element Labels Through the Pipeline [DONE 2026-02-22]
- Added `name` attribute to `lattice` base class and all subclasses.
- Excel path: Label column (Nomenclature fallback) propagated through
  `excelElements.create_beamline()` and `beamlineBuilder.parse_beamline()`.
- JSON/YAML path: `name` field propagated through `latticeLoaderBase` to both
  beamline objects (`create_beamline()`) and dict output (`parse_beamline()`).
- Also fixed `tracked_mapping` → `tracked_dict` import (stale module name)
  in 4 files, resolving 6 smoke test failures.

### I3. PALS / v2 Lattice → COSY INFINITY Converter [v0.3.0 DONE 2026-03-03]
- **Goal:** Standalone tool that converts PALS-aligned lattice files
  (JSON or YAML) into COSY INFINITY FOX input decks.
- **v0.3.0 changes:**
  1. **Official PALS format support:** Parses `PALS:` root key, `facility:` element
     definitions, `BeamLine` composition with `line:` references, `inherit:` overrides,
     `repeat:` expansion. New `pals_parser.py` module (~230 lines).
  2. **MagneticMultipoleP.Bn1:** Quadrupole pole-tip field (Tesla) emitted directly
     as COSY MQ `b_pole` — no gradient coefficient needed.
  3. **BendP support:** `g_ref` (1/m) → angle (deg), `e1`/`e2` (rad) → edge angles.
  4. **CLI auto-detection:** `--mode auto` (new default) routes to the correct parser
     based on root key (`PALS:` vs `beamline:`). New `--beamline NAME` flag.
  5. **Example lattices:** `fodo.pals.yaml` (official PALS FODO cell),
     `uhfel_excerpt.pals.yaml` (UH FEL first section with Bn1 quads + BendP dipole).
  6. **Test suite:** 55 passing tests (was 34).
  7. **Sphinx documentation** updated with both formats, architecture, element mapping.
- **v0.2.0 changes:** Particle type, element comments, strict PALS mode, FC
  suppression, 34 passing tests.
- **Known issue:** YAML has 8 spectrometer dipoles vs 4 in Excel → COSY round-trip
  fails with "no fixed point". Requires reconciling YAML beamline endpoint.
- **Deferred:** Optimization (FIT blocks), MGE/fieldmap, particle tracking,
  solenoid/RF cavity/sextupole physics
- **Future work — 3-way MAD↔PALS↔COSY converter umbrella project:**
  1. Build a 3-way converter supporting MAD-X, PALS, and COSY INFINITY
     formats (any→any direction). MAD-X is the de facto standard for
     accelerator lattice interchange.
  2. Validate via round-trip QA loop: MAD → PALS → COSY → MAD. Use the
     MAD-X lattice library on this machine as a test suite of real-world
     lattices covering diverse element types and beamline topologies.
  3. The public `pals2cosy` repo remains a focused, standalone PALS→COSY
     converter — simple and easy to use. The 3-way umbrella project
     informs its design and provides QA, but `pals2cosy` does not depend
     on the umbrella. Shared internal representations and conversion
     logic can be extracted from the umbrella into pals2cosy as needed.
- **Repo:** `~/COSY/PALS2COSY/` (git: evvaletov/pals2cosy)

### I4. COSY Aperture Commands for Particle Tracking [DONE 2026-02-22]
- AP commands generated after each element when `enable_aperture_cuts()` is called.
  Quads: elliptic `AP r r 1` (r = quad_aperture/2).
  Dipoles: rectangular `AP w h 2` (h = pole_gap/2, w = dipole_half_width).
- Opt-in via `sim.enable_aperture_cuts(dipole_half_width=0.050)`.
- Forwarded through `COSYAdapter.enable_aperture_cuts()`.
- Particle loss robustness: 0-ray handling in `_read_rray_format()`,
  N<2 guard in `calculate_twiss_from_particles()`, transmission logging in
  `read_checkpoints()` and `collect_evolution()`, graceful `all_particles_lost`
  metadata in `simulate()`.
- **TODO:** Determine actual UH MkV dipole pole face width (currently 50 mm placeholder).

### I6. MCNP-Style Robustness & Foolproofness [IN PROGRESS]
- **Motivation:** MCNP is a gold standard for production code robustness: every
  input is validated, edge cases are caught with clear diagnostics, defaults are
  sensible, and the code never silently produces wrong results. FELsim should
  adopt this level of rigour.
- **Actions:**
  1. Input validation at system boundaries: lattice files, API payloads,
     CLI arguments, Excel data. Fail loudly with descriptive errors.
  2. Guard against silent numerical failures: NaN/Inf propagation,
     singular matrices, zero-length elements, particle loss without warning.
  3. Consistent error handling: no bare `except:`, no swallowed exceptions.
     Every failure path either recovers correctly or raises with context.
  4. Default values must be physically sensible (not 0 or 1 by convenience).
  5. Audit all `setattr`/`getattr` patterns for typo-resilience (consider
     `__slots__` or property validation on beamline element classes).
  6. Configuration validation: warn on unused/unknown keys, reject
     contradictory settings.
- **Progress (2026-03-10):**
  - `beamline.py`: `setE()` and `setMQE()` validate inputs (positive, finite).
    `dipole_wedge` guards zero-angle case in both `_compute_numeric_matrix`
    and chromatic `useMatrice`.
  - `ebeam.py`: `ellipse_sym()` sqrt guard, `cal_twiss()` epsilon threshold
    uses `np.finfo(float).tiny` instead of arbitrary 1e-30.
  - `beamOptimizer.py`: NaN/Inf guard after MSE computation (`_optiSpeed`).
  - `cosyParticleSimulator.py`: bare `except:` narrowed to
    `except (ValueError, IndexError):`.
  - `felAPI.py`: Silently caught exceptions now logged.
  - **Broad except audit (8 handlers narrowed):**
    `simulatorBase.py` (2), `cosySimulator.py` (1), `cosyAdapter.py` (3),
    `beamlineBuilder.py` (2). All changed from `except Exception` to specific
    exception tuples.
  - Configuration validation: `_report_unaccessed()` elevated from `info`
    to `warning` level for unused lattice keys.
  - Attribute typo guard: `test_attribute_guard.py` (13 tests) — AST-based
    scan of all backend setattr calls, validates targets against element
    class attribute whitelist, checks for suspicious near-duplicate names.
  - **Current test suite:** 176 pass, 7 skip, 0 fail across 7 test modules.

### I7. Multi-Code Simulation Framework [DONE 2026-03-10]
- **Motivation:** Different simulation codes have different strengths: RF-Track
  excels at 3D space charge, COSY INFINITY excels at high-order DA maps and
  fringe fields. A production beamline study should be able to use RF-Track
  for one section and COSY for another — seamlessly and configurably.
- **Implementation:**
  - `CoordinateTransformer` (simulatorFactory.py): all 6 pairwise transforms
    (FELSIM↔COSY, FELSIM↔RFTRACK, COSY↔RFTRACK) as static methods
  - `MultiCodeSimulator` (multiCodeSimulator.py): orchestrator that chains
    multiple SimulatorBase instances on contiguous beamline sections.
    All adapters use FELsim coordinates as I/O format (each handles its
    own internal transforms), so no inter-section coordinate conversion
    is needed.
  - `_felsim_to_generic()`: converts FELsim native elements to generic
    `BeamlineElement` for non-FELsim adapters, preserving all parameters
    including DPW pole_gap, dipole_angle, dipole_length
  - `SimSection` dataclass: (name, simulator_key, element_range, config)
  - `SimulatorFactory.create('multicode', ...)`: lazy-imported registration
  - `from_config()`: dict-based construction for YAML/JSON configuration
  - Test suite (test_multicode.py): 22 tests — SimSection, init, coord
    roundtrips (all 3 pairs), FELsim split equivalence (2/3-section),
    element conversion (drift, quad, DPW params), factory registration,
    FELsim→RF-Track hybrid (successful run, cross-validation, DPW params)
  - CI: test_multicode.py + test_attribute_guard.py added to pipeline
- **Production validated:** FELsim→RF-Track hybrid at Stage 11 boundary
  (element 87) runs successfully. Hybrid vs full RF-Track shows qualitatively
  similar results (transverse RMS within order of magnitude) with expected
  differences from dipole model (transfer matrix vs analytical sector-bend).
  - COSY→FELsim handoff testing with real DA map tracking
  - Per-section config passthrough (space_charge, fringe fields, etc.)

### I5. T566 Objective via 2nd-Order DA Map [LOW PRIORITY — NOT NEEDED FOR UH FEL]
- **Status:** `("l", "t566")` is in MEASURE_MAP but raises NotImplementedError.
- **Goal:** Extract T566 = (∂²l/∂δ²)/2 from the COSY DA polynomial map
  and use it as a FIT objective. Requires `transfer_matrix_order >= 2`.
- **Use case:** Bunch compression optimization where both R56 and T566 matter.
- **W12 finding:** T566 = 0 for the UH FEL transport line (W9 Part A).
  A T566 FIT objective is redundant for this beamline. Implementation
  may still be useful for other beamlines with non-zero T566.

### I8. FELsim v3 Lattice Format — PALS Alignment [DONE 2026-03-03]
- **Goal:** Extend the v2 lattice format with optional PALS-aligned fields while
  maintaining full backward compatibility with v1/v2 files.
- **Changes:**
  1. **Specification:** `manuals/lattice_specification_v3.md` — documents
     `MagneticMultipoleP.Bn1`, `BendP` (g_ref, e1, e2), `angle_unit`, precedence
     rules, backward compatibility matrix.
  2. **JSON Schema:** `var/lattice_schema_v3.json` — extends v2 schema with
     `MagneticMultipoleP`, `BendP` objects, `angle_unit` enum, `format_version` enum [1,2,3].
  3. **Loader:** `latticeLoaderBase.py` updated — `SUPPORTED_FORMAT_VERSIONS = [1, 2, 3]`,
     Bn1 → current conversion for quads (`current = Bn1 / (G × r)`), BendP override
     for dipole angle and edge angles. Both `parse_beamline()` and `create_beamline()`
     paths handle v3 fields.
  4. **Documentation:** `docs/felsim/lattice-formats.md` updated with v3 section.
- **Backward compatibility:** v1/v2 files load identically. Verified with existing
  v2 YAML (137 elements from `create_beamline()`).
- **Deferred:** Full PALS root key (`PALS:`) support in FELsim loaders, `line:`
  composition, implicit positioning from lengths.

### I2. Engage with PALS as a Real-World Use Case [HIGH PRIORITY]
- **Context:** The PALS (Particle Accelerator Lattice Standard) group is looking
  for ~10 real-world use cases of their evolving standard before proceeding to
  the next stage. FELsim's v2 lattice format already uses PALS-aligned type
  names (Quadrupole, SBend, Wiggler, Kicker, Instrument, etc.).
- **Goal:** Register UH MkV FEL / FELsim as a PALS use case. Contribute
  feedback on gaps (e.g., DIPOLE_WEDGE has no PALS equivalent, label/name
  conventions, polarity handling).
- **Actions:**
  1. Identify the PALS working group contact and submission process
  2. Prepare a short description of the UH FEL lattice and how FELsim uses
     the PALS-aligned format (JSON + YAML, 118 elements, 23 quads, dipole
     sandwiches, chromaticity sections)
  3. Document feedback from our implementation experience (e.g., label/name
     conventions, polarity handling). Note: DIPOLE_WEDGE is not a real gap —
     it is an internal FELsim modeling artifact (thin-lens edge kick); in the
     interchange format, edge angles and Enge coefficients are attributes of
     the dipole element (SBend/RBend), which is the standard approach.
     However, DIPOLE_WEDGE could be worth mentioning to PALS as a pattern
     for ingesting legacy or substandard lattice specifications where dipole
     edges are defined as separate elements. If the standard incorporates
     this, it should be a distinct category — not deprecated in the standard
     itself, but marked as a "legacy compatibility" or "suboptimal" element
     type, signaling that individual lattice files using it should migrate
     toward folding edge parameters into the parent dipole element
  4. Submit as a use case

---

## Category C: Cross-Validation

### C1. RF-Track Cross-Validation & Optimisation [DONE 2026-02-24]
- **Motivation:** The FELsim optimizer uses transfer matrices; RF-Track uses
  full particle tracking. Comparing at key parameter points validates the
  transfer matrix model.
- **Part A (S9):** RF-Track validation with FELsim-optimised currents. DONE.
- **Part B (W8):** Hybrid FELsim→RF-Track Stage 11 optimisation. DONE.
  - FELsim runs stages 1–10 (fast transfer matrices), then RF-Track
    particle tracking optimises Stage 11 (4 quads → undulator Twiss match).
  - Prefix caching: elements 0:87 pre-tracked once (~0.2 s), suffix trackings
    (87:118 + 87:93) per NM eval. Full optimisation ~80 s per emittance point.
  - **Key findings (ε_n=8 smoke test):**
    - FELsim currents in RF-Track give MSE=2552 (β_x=86 vs target 1.4).
      Root cause: accumulated dipole edge-kick model differences through 118
      elements shift the beam state enough that FELsim-optimal currents are
      far from RF-Track-optimal.
    - RF-Track optimiser finds its own solution with MSE=6.4e-6 (125× better
      than FELsim's 8.1e-4), but with very different Stage 11 currents.
    - This validates that RF-Track particle tracking can match the undulator
      Twiss targets more precisely than FELsim's transfer matrices.
  - **Adapter fixes applied:**
    - `DIPOLE_WEDGE → Drift(0)`: FELsim models wedges as thin-lens edge kicks
      (no drift propagation); RF-Track was adding spurious 10 mm drifts.
    - `Quadrupole.set_strength(k1*L)`: Verified correct — `set_strength(S)`
      is internally used as k1*L, not the manual's nominal P/q·k1·L formula.
  - Script: `UHM_rftrack_opt.py` (`--smoke`, `--emittance`, `--space-charge`)
  - Results: `results/rftrack_opt/`
  - TODO: Run full comparison at ε_n = 5, 8, 14 with 5 restarts.

### C3. FR3+MGE Optimization [IN PROGRESS — CMA-ES v2 on Koa]
- **Fieldmap fix (2026-02-22):** DELTAS 0→0.001, removed 0.835× scaling. Peak field
  now 0.5307 T matching source CSV.
- **COSY FIT attempts (2026-03-07/08):**
  1. Graduated chain FR1→FR2→FR3→FR3+MGE: Steps 1-2 succeed (MSE<1e-6), Step 3
     (FR3+MGE) fails — Stage 11 diverges (MSE=8.9e10).
  2. Direct FR1→FR3+MGE: same failure (MSE=6.4e12).
  3. Stage-11-only with FR3 upstream: 2D parameter scan shows NO stable solutions
     exist with fixed Stages 1-10 under MGE.
  4. Koa warm-start (MSE=1030→FIT): Stages 1-10 FIT re-optimized to different values,
     destroying Koa's global consistency → Stage 11 diverged (S11_Ic=-8.67 A,
     MSE=6.0e12).
- **Root cause:** Sequential FIT cannot solve FR3+MGE. The stability island for all
  23 variables is far from the per-stage local optima. DA gradients work within a
  stage but the sequential decomposition loses global coupling.
- **CMA-ES v1 (Koa job 11427816, 2026-03-07):** 10,020 evals, sigma=2→0.5,
  bounds 0-10 for some vars. MSE=1030 (cos_mu≈1.03, barely unstable). Warm-start
  refinement at any sigma (0.01-0.1) fails — stability boundary too steep.
  Result: `test/results/koa_cosy_mge_result.json`
- **CMA-ES v2 (Koa job 11451841, submitted 2026-03-08):** 50,000 evals, sigma=0.5,
  all bounds [-10,10], BIPOP restarts ×9, warm-start from v1 MSE=1030.
  Est. runtime: ~89 hours (~4 days). Script: `test/koa_cosy_mge_opt.py`
  - **Status (2026-03-10, 44.5 hrs):** 15,392/50,000 evals (31%), MSE=1000.005
    (best, barely unstable). Sigma collapsed to 4.7e-05 (from 0.5), axis ratio
    530:1. First BIPOP restart imminent — sigma exhausted around current best.
    No stable solution found in this basin; restart will explore larger volume.
- **Files:** `fields/chicane_dipole_fieldmap.dat`, `test/koa_cosy_mge_opt.py`,
  `test/koa_cosy_mge_opt.slurm`, `test/results/koa_cosy_mge_result.json`,
  `test/results/koa_cosy_mge_result_indexed.json`

### C4. Systematic Testing, Validation & Verification [IN PROGRESS]
- **Motivation:** FELsim currently relies on ad-hoc cross-validation studies.
  A systematic V&V programme is needed for production confidence.
- **Actions:**
  1. **Unit tests:** Core physics routines (transfer matrices, Twiss
     computation, dispersion, coordinate transforms) need pytest coverage
     with known analytic results (e.g., thin-lens quad, drift, FODO).
  2. **Regression tests:** Each optimization study should produce a frozen
     reference result. CI runs confirm that code changes don't alter results
     beyond numerical noise.
  3. **Cross-code benchmarks:** Extend S9/C1/C2 pattern — for each major
     beamline section, compare FELsim, RF-Track, and COSY Twiss functions
     element-by-element. Automate as a benchmark suite.
  4. **Edge case testing:** ε_n → 0, σ_E → 0, single particle, 10⁵ particles,
     zero-length elements, degenerate optics (β → ∞).
  5. **Adapter round-trip tests:** Load lattice in all three formats
     (Excel/JSON/YAML), verify identical beamline objects.
  6. **CI pipeline:** Automated test runs on commit (at minimum: unit tests
     + adapter round-trip + one optimization smoke test).
- **Output:** `backend/test/test_*.py` files, CI config, benchmark report.
- **Progress (2026-03-10):**
  - Created `test_chromatic_physics.py`: 30 tests covering chromatic quads,
    dipoles, wedges, apertures, and FODO integration — all pass.
  - Fixed `test_transfer_matrices.py`: quad known-element test now uses
    element's own constants (was using truncated test-file constants).
  - Fixed `test_fieldmap_validation.py`: `np.trapezoid` → `np.trapz` for
    numpy 1.x compatibility.
  - Fixed `test_adapter_roundtrip.py`: YAML tests skip gracefully when
    `tracked_dict` is unavailable (requires Python ≥3.10).
  - Fixed `latticeLoaderBase.py`, `jsonLatticeLoader.py`: `TrackedDict`
    import made optional with recursive fallback class.
  - Added frozen regression tests to `test_adapter_roundtrip.py`: cumulative
    matrix entries, R56, symplecticity, and particle output at 40 MeV — all
    frozen from validated runs to catch silent physics changes.
  - **Current test suite:** 156 pass, 7 skip, 0 fail across 6 test modules
    (`test_chromatic_physics`, `test_transfer_matrices`, `test_twiss`,
    `test_adapter_roundtrip`, `test_optimizer`, `test_fieldmap_validation`).
  - **Remaining:** CI pipeline, cross-code benchmark automation.

### C2. COSY INFINITY Cross-Validation [DONE — see W4]
- **Motivation:** Independent DA-based simulation. Particularly valuable for
  verifying dispersion and chromatic behavior at high σ_E.
- **Design:** Run COSY optimisation directly (W4), then compare element-by-element
  Twiss functions with the Python results.
- **Due:** 2026-02-18 (W4)

### C8. Chromatic Quadrupole Matrices [DONE 2026-03-10]
- **Problem:** RF-Track applies momentum-dependent focusing (k_eff = k₀ × P₀/P),
  while FELsim uses fixed k₀ for all particles. At σ_p = 0.5%, this causes
  β_y to diverge by 65–1100% between codes.
- **Solution:** Added `chromatic` flag to `lattice` base class (default False).
  When enabled, `qpfLattice.useMatrice()` and `qpdLattice.useMatrice()` compute
  per-particle k from coord6 (ΔK/K₀ × 10³) using vectorized numpy operations.
  Formula: k = k₀ × βγ₀/βγ(δ) ≡ k₀ × P₀/P.
- **Validation (D15):** Three-way comparison (linear FELsim, chromatic FELsim,
  RF-Track) at σ_p = 0%, 0.1%, 0.25%, 0.5%, 1.0%:
  - Chromatic FELsim ≡ RF-Track to 0.00% across all Twiss parameters (β, α, ε, σ)
  - Linear FELsim errors: β up to 1466%, ε up to 98%, σ up to 90%
- **Properties:** Opt-in, backward-compatible, symplectic (det M = 1.0 to
  machine precision), no Python loops.
- **Files:** `beamline.py` (qpfLattice, qpdLattice), `test/test_chromatic_quad.py`,
  `test/D15_chromatic_comparison.py`

---

## Category P: Physics Model Upgrades

Items are ordered by estimated impact on the UH MkV FEL beamline
(40 MeV, σ_p = 0.5%, 2 ps bunch, chicane + 26 quads, 118 elements).

### P1. Aperture Loss Tracking in FELsim [DONE 2026-03-10]
- **Motivation:** RF-Track shows 40–50% particle loss at the 0.5 ps compression
  operating point (W10, W11). FELsim never drops particles, so its optimizer
  can find "optimal" solutions that are physically unrealizable (beam clips
  apertures, especially in chicane dipoles with 14.5 mm pole gap).
- **Impact:** >5%. Largest unmodeled effect. Directly changes optimization outcomes.
- **Implementation:**
  1. Added `aperture_x`, `aperture_y` attributes to `lattice` base class (default None).
  2. Added `apply_aperture(particles)` method: masks particles exceeding the aperture,
     returns surviving subset. `useMatrice` unchanged (no API break).
  3. Default apertures set automatically:
     - QPF/QPD: ±13.5 mm (27 mm bore, class attribute `BORE_RADIUS_MM`)
     - DPH: ±pole_gap/2 in y (from Excel data, passed via new `pole_gap` parameter)
     - DPW: ±pole_gap/2 in y (from existing `self.pole_gap` attribute)
     - Drift: no aperture (None)
  4. Tracking loop calls `apply_aperture` after each `useMatrice`.
- **Validation (FELsim-optimized currents, 2000 particles, σ_p=0.5%):**
  - 64.5% transmission (35.5% loss = 710 particles)
  - Dominant loss: FC1.DPW.111 (element 100, ±6.3 mm) — 395 particles (56% of losses)
  - Chicane exit dipoles (72–73, 89–90): 281 particles
  - Quad losses minor: 34 particles total
- **Optimizer integration:** `calc(..., use_apertures=True)` enables aperture cuts
  in `_optiSpeed`. Smooth penalty: returns `1e6 × (1-T)` if <2 particles survive
  (can't compute Twiss). Opt-in, backward-compatible (default: `use_apertures=False`).
- **Transmission objective [DONE 2026-03-10]:** `calc(..., transmission_weight=W)` adds
  `W × (1-T)²` to MSE, where T = n_surviving/n_initial. Tracks `transmission` in
  `trackGoals` for plotting. Backward-compatible (default weight = 0).
- **Chromatic re-optimization [DONE 2026-03-10]:** Multi-pass coordinate-descent
  optimization with `chromatic=True` + `use_apertures=True` + `transmission_weight=5`.
  3 passes over 11-stage NM sequence. Best result (Pass 2):
  - MSE = 4.0e-4 (β_x=1.40, α_x=0.44, β_y=0.23, α_y=0.03)
  - Transmission = 98.1% (up from 64.5% with linear-optimized currents)
  - Significant current changes in chicane region (elements 33-43) and
    downstream quads (70, 87, 93-97)
  - Script: `test/UHM_beamline_opt_chromatic.py`, results: `results/felsim_chromatic_warm.json`
- **Files:** `beamline.py` (lattice base, qpfLattice, qpdLattice, dipole, dipole_wedge),
  `excelElements.py` (pole_gap → dipole constructor)

### P2. Chromatic DPW Edge Kicks [DONE 2026-03-10]
- **Motivation:** The DPW edge kick is tan(η)/R where R = ρ = L/|θ|. For
  off-momentum particles, ρ ∝ P (magnetic rigidity), so the kick strength
  scales as P₀/P — same mechanism as quad chromaticity.
- **Impact:** ~2% on β_y per 1‰ energy offset in chicane regions.
- **Implementation:**
  1. Override `useMatrice` in `dipole_wedge` (same pattern as chromatic quads).
  2. Per-particle R: R_p = R₀ × βγ_p/βγ₀.
  3. Both h = 1/R and the fringe correction φ = K·g·h·(1+sin²η)/cos(η) are
     momentum-dependent (h_p = h₀ × P₀/P).
  4. Full chromatic DPW kick: M21_p = Tx/R_p, M43_p = -Ty_p/R_p where
     Ty_p = tan(η - φ_p). Vectorized over all particles.
- **Validation:**
  - On-momentum: max|Δ| < 10⁻¹² vs linear (exact consistency)
  - Off-momentum: Δyp ≈ 1.4e-6/mrad per ‰ of δ (correctly momentum-scaled)
  - Symplecticity: det(M_x) = det(M_y) = 1.0 to machine precision
- **Note:** FELsim chromatic DPW is more accurate than RF-Track's DPW model
  (RF-Track uses K1L = -|K0|·tan(η) without fringe correction φ). COSY's
  DIL includes all edge effects via DA — inherently chromatic.
- **Files:** `beamline.py` (dipole_wedge.useMatrice)

### P3. Chromatic Sector-Bend Body [DONE 2026-03-10]
- **Motivation:** Dipole body focusing (cos θ, sin θ) uses reference ρ.
  Off-momentum particles see different bending radii (ρ ∝ P = magnetic rigidity).
  Per-particle θ = L/ρ(δ) produces momentum-dependent bending, dispersion, and R56.
- **Impact:** <1% for σ_p = 0.5% (dominated by existing M16/M26 terms), but
  needed for model consistency with chromatic quads and edges.
- **Implementation:**
  1. Override `useMatrice` in `dipole` class (same pattern as chromatic quads).
  2. Per-particle ρ_p = ρ₀ × βγ_p/βγ₀, θ_p = L/ρ_p.
  3. Full sector-bend body: M11 = cos(θ_p), M12 = ρ_p·sin(θ_p), M16 = ρ_p·(1-cos θ_p)·gfac·δ,
     M21 = -sin(θ_p)/ρ_p, M22 = cos(θ_p), M26 = sin(θ_p)·gfac·δ. Y-plane: drift.
  4. Longitudinal: M51, M52, R56 all per-particle (momentum-dependent dispersion).
- **Validation:**
  - On-momentum: max|Δ| < 10⁻¹⁰ vs linear matrix (exact consistency)
  - Off-momentum: chromaticity effect = 9.8e-2 (correct direction)
  - Symplecticity: det(Mx) = det(My) = 1.0 at δ = 0, ±5, +10
  - ρ scaling: verified ρ_p/ρ₀ = βγ_p/βγ₀ at δ = 0.1–2.0%
- **Files:** `beamline.py` (dipole.useMatrice), `test/test_chromatic_dipole.py`

### P4. Fringe Field Sextupole Terms [LOW PRIORITY]
- **Motivation:** COSY with FR ≥ 1 includes Enge sextupole contributions from
  fringe fields. Neither FELsim nor RF-Track model these.
- **Impact:** ~2–5% on β_y for high-order dipole fringe fields. Negligible for
  linear Enge model (FR 0).
- **Design:** Would require extracting sextupole kick from Enge function gradient.
  More naturally handled by COSY (already does this via DA).
- **Effort:** ~1 week. Medium difficulty (Enge coefficient extraction and validation).

### P5. Chromatic M56 [LOW PRIORITY]
- **Motivation:** M56 = -L·f/(c·β·γ·(γ+1)) uses reference β, γ. Per-particle
  M56 would capture velocity-dependent time-of-flight.
- **Impact:** <0.1% for σ_p = 0.5% at 40 MeV. Negligible.
- **Effort:** Trivial (same per-particle γ calculation already in chromatic quads).

### P6. CSR Module [LOW PRIORITY — NOT NEEDED FOR UH FEL]
- **Motivation:** Coherent Synchrotron Radiation causes energy loss and emittance
  growth in short bunches traversing bends. Derbenev-Saldin formula gives
  Δσ_E ≈ 0.5% for 0.5 ps bunch at UH FEL chicane parameters.
- **Impact:** Negligible at 2 ps operating point. Potentially significant at <1 ps,
  but the beamline cannot compress to that regime anyway (W12 conclusion).
- **Design:** 1D CSR solver (steady-state or transient). External module.
- **Effort:** ~2 weeks.

### P7. COSY Model Assessment [DONE 2026-03-10]
- **Finding:** The COSY model uses NO kicks. Both element types are proper
  thick-element representations:
  - Quadrupoles: `MQ L B_pole R` — thick magnetic quadrupole (full DA map)
  - Dipoles: `DIL L θ g/2 e₁ 0 e₂ 0` — thick dipole with integrated edge angles
  - DPW-DPH-DPW triplets are consolidated into single DIPOLE_CONSOLIDATED
    with entrance/exit edge angles passed to DIL
  - Fringe fields via FC/FD (Enge) or MGE (fieldmap)
- **Conclusion:** No element replacement needed in COSY. The DIL + FC/FD approach
  is the most accurate available representation.

---

## Workflow

- **Deadlines:** Complete full beamline (to-undulator) studies relevant to FEL
  operations/physics by **Monday** each week. The real presentation target is
  Wednesday, but Monday gives buffer time for analysis and review. This deadline
  does not apply to simulation R&D (optimizer improvements, code work, etc.).
- **Reports:** LaTeX format with figures and tables. The user reads the LaTeX
  report and makes their own PowerPoint from it. Do not generate PPTX for weekly
  reports — only when specifically requested.
- **Report revisions:** Keep notes about requested revisions to optimise reports
  over time (see `report_style.md` if it exists).
- **Units:** Always use correct, consistent units in reports and presentations.
  Key conventions:
  - Normalised emittance: π·mm·mrad (= μm = μm·rad)
  - Energy spread: % of total energy (σ_E / E₀ × 100)
  - Chirp: 1/s (energy-time correlation coefficient)
  - Beta function: m
  - Current: A
  - Bunch length: ps

## Known Issues (from QA Passes 1–4)

- **CRITICAL (pre-existing):** `/plot-parameters` endpoint (felAPI.py:580) references non-existent
  `Beamline` class and `findSegmentAtPos` method — endpoint always crashes with NameError.
  Part of the frontend-facing API (Christian's responsibility to fix or coordinate).
- **MODERATE:** Mutable default `shape={}` in `ebeam.plotXYZ` (ebeam.py:568) and
  `draw_beamline.plotBeamPositionTransform` (schematic.py:273). Shared across calls.
  Low risk since the dict is only read, not mutated, but violates Python best practice.
- **MODERATE (needs domain knowledge):** `dipole_wedge` fringe field integral
  (beamline.py:707) uses `le = self.length` which comes from `gap_wedge` in the Excel
  lattice. Needs clarification whether `gap_wedge` is the fringe field extent or the
  inter-dipole drift gap — if the latter, the K integral is wrong.
- **MODERATE:** `ebeam.ellipse_sym` (ebeam.py:79–82) divides by Twiss beta or gamma
  without zero guard. Crashes on degenerate beam distributions.
- **MODERATE:** `AlgebraicOptimization.py` (lines 276–277) uses `set.pop()` for variable
  extraction — non-deterministic ordering can silently swap x/y solutions between runs.

## Longer-Term Improvements (from multi-AI review 2026-03-10)

Source: 4-perspective expert review (FEL scientist, Berz-style computational physicist, SWE, UX/UI)

### Physics & Validation
- [ ] **Order-by-order convergence study**: Run COSY at orders 1, 2, 3, 5 — quantify how optimized currents and final phase space change with map order
- [ ] **Emittance preservation plot**: ε_n(s) through the chicane to quantify CSR-driven emittance growth
- [ ] **Chromaticity analysis**: Twiss parameters and beam size as a function of energy deviation δ
- [ ] **Fringe field treatment in FELsim**: Currently fringe_field_order=0; currents are optimized for the wrong model. Add fringe field support or at least quantify the impact
- [ ] **Sensitivity / error analysis**: Magnet errors, misalignments, power supply ripple — DA methods can compute high-order sensitivities directly
- [ ] **Multi-seed robustness study**: Run optimization with seeds 42, 137, 2023+ and compare results to confirm global minimum was found

### Code Quality
- [ ] **Decouple optimization from visualization**: Cache optimization results (HDF5/NPZ) so figure generation doesn't require re-running the 2-minute optimization
- [ ] **Add pytest test suite**: Unit tests for Twiss computation, integration test for simplified optimization, visual regression with pytest-mpl
- [ ] **Extract rcParams to .mplstyle file**: Reusable seminar style across projects
- [ ] **Make output formats configurable**: argparse for PDF/PNG selection

## Notes

- **Frontend ownership:** The frontend (`fel-app/`) is developed exclusively by Christian Komo. QA and code changes should focus on the backend. Minor frontend improvements are acceptable, but avoid making extensive changes that step on his work.
- All optimization scripts use `seed=42` for reproducibility
- NewFELsim conda environment required: `/home/evaletov/.conda/envs/NewFELsim/bin/python`
- Run commands from `backend/` with `MPLBACKEND=Agg PYTHONPATH=$(pwd)`
- **Review methodology**: Apply Berz-style computational physics perspective regularly when reviewing simulation results. Use multi-AI collab for second opinions.
