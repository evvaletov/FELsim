# UH MkV FEL Beamline Optimization — Priorities & Roadmap

**Date:** 2026-02-11 (updated 2026-02-24)
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

### W9. COSY Longitudinal Study [IN PROGRESS]
- Full 3D (6D phase space) COSY simulation with longitudinal diagnostics
- Extracts R56, T566, coupling terms from optimised beamline
- Propagates 6D bunches for 0.5 ps and 2 ps modes
- Tests adding R56=0 as an optimisation objective
- Reviewer verdict on what changes when switching bunch length
- Script: `W9_cosy_longitudinal_study.py`
- Results: `results/W9/`

### I6. COSY σ_z Blowup Investigation [HIGH PRIORITY — MANUAL]
- COSY particle tracking shows σ_z ≈ 92–233 ps (60–100× blowup) for 2 ps input.
  RF-Track gives correct σ_z ≈ 2 ps for the same currents.
- Both 2 ps and 0.5 ps inputs produce identical output σ_z, confirming the output
  is dominated by energy-spread × R56 coupling, not initial bunch length.
- **Not a coordinate conversion bug:** verified FELsim ↔ COSY round-trip is exact.
  The FELsim col 4 std genuinely grows 60× during COSY tracking.
- **Possible causes:** (a) COSY higher-order tracking amplifies dispersion errors
  from the poorly matched optics (MSE = 7.1e-3); (b) COSY's l₀ coordinate picks
  up path-length errors in the chicane that FELsim transfer matrices miss;
  (c) the COSY FOX RP procedure may define the 5th coordinate differently than
  assumed in `cosyParticleSimulator.transform_from_cosy_coordinates()`.
- **Action:** Manually inspect COSY RRAY checkpoint files, trace l₀ evolution
  element-by-element through the chicane, compare with FELsim prediction.
  Check RP settings for coordinate convention.
- **Expected resolution:** Either fix a unit convention mismatch or confirm that
  better optics (lower MSE from glyfada) eliminate the blowup.

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
- **Results (pre-C5 fix, with missing edge angles):**
  - 2ps transport: MSE=0.175, T=42.8%, σ_t=1.83 ps, I_peak=5.6 A
  - 0.5ps compress: MSE=0.071, T=28.4%, σ_t=1.31 ps (target 0.5), I_peak=5.2 A
  - Compression scenario failed — σ_t target unreachable with Stage 11 quads alone
- **Note:** C5 workaround now implemented (analytical dipole correction). Re-run
  to get updated results with correct dipole physics.
- Script: `W11_throughput_opt.py`
- Results: `results/W11/`

### W10. Beam Losses & Bunch Compression Study [IN PROGRESS]
- Quantifies particle losses through the full transport line with physical apertures
- Part A: Transmission baseline at 2 ps and 0.5 ps (COSY + RF-Track)
- Part B: Bunch compression via negative chirp — demonstrates chirp required,
  energy spread alone does NOT compress (R56 ≈ +27 mm elongates unchirped beams)
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

### S5. 0.5 ps 2D Coupled Scans [HIGH PRIORITY]
- **Motivation:** 1D scans hold other parameters fixed; realistic operation involves
  correlated changes (e.g., shorter bunch → larger energy spread). 2D scans map
  the feasibility surface.
- **Scans:**
  - S5a: (σ_E, h) grid at ε_n=8 — energy spread vs chirp coupling
  - S5b: (σ_E, ε_n) grid at h=5e9 — degradation interaction
  - S5c: (h, ε_n) grid — chirp compensation vs emittance
- **Design:** 10×10 grids (100 points each), 500 particles, ~5 hours per scan.
  Use contour/heatmap plots of MSE and Twiss deviation.
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

### R2. Comparison Table Across All Studies [MEDIUM PRIORITY]
- **Motivation:** S1–S4 each have separate reports. A unified comparison table
  showing how the same beamline responds to different parameter regimes.
- **Output:** Standalone comparison document or section in a master report.

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

### I3. PALS / v2 Lattice → COSY INFINITY Converter [HIGH PRIORITY — due ~2026-03-08]
- **Goal:** Standalone tool that converts a v2 (PALS-aligned) lattice file
  (JSON or YAML) into a COSY INFINITY FOX input deck. Intended for the
  official COSY INFINITY website as a community resource.
- **Design considerations:**
  - Mode switch: strict FELsim v2 (our superset with DIPOLE_WEDGE, etc.)
    vs generic PALS (broader compatibility for general PALS-to-COSY use)
  - Output: complete FOX procedure with element definitions, beamline
    sequence, and initial beam parameters (energy, particle type)
  - Element mapping: Quadrupole→MQ, SBend/RBend→DI/MC, Drift→DL,
    Solenoid→SOLND, Sextupole→MH, RFCavity→RFC, diagnostics→DL(0)
  - Fringe field options: FR 0/1/3 selectable; Enge coefficients carried
    through where available
  - Should work independently of FELsim runtime (no FastAPI dependency)
- **Scope:**
  1. Core converter module (reads JSON/YAML, emits FOX)
  2. CLI entry point (`pals2cosy` or `lattice2cosy`)
  3. Validation against our UH FEL beamline (round-trip: YAML → FOX → COSY
     run → compare Twiss with existing W4 results)
  4. Minimal documentation / README for the COSY website

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

### I6. MCNP-Style Robustness & Foolproofness [HIGH PRIORITY]
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

### I7. Multi-Code Simulation Framework [HIGH PRIORITY]
- **Motivation:** Different simulation codes have different strengths: RF-Track
  excels at 3D space charge, COSY INFINITY excels at high-order DA maps and
  fringe fields. A production beamline study should be able to use RF-Track
  for one section and COSY for another — seamlessly and configurably.
- **Design goals:**
  1. Per-section code assignment: lattice file or config specifies which
     simulator handles each beamline segment (e.g., elements 0–86 with COSY,
     87–117 with RF-Track).
  2. Beam state handoff: well-defined coordinate transforms between codes
     at junction points. Currently `transform_coordinates()` handles
     FELsim↔RF-Track; extend to include COSY particle format.
  3. Unified result format: `SimulationResult` already provides this;
     ensure all adapters populate it consistently.
  4. Configuration: YAML/JSON config with per-section `simulator` key
     (e.g., `{simulator: rftrack, elements: [87, 117], space_charge: true}`).
  5. Prototype: the hybrid FELsim→RF-Track Stage 11 optimisation (C1/W8) is
     the first instance of this pattern. Generalise from there.
- **Prerequisite:** C1 Part B (RF-Track optimisation) validates the handoff approach.

### I5. T566 Objective via 2nd-Order DA Map [LOW PRIORITY]
- **Status:** `("l", "t566")` is in MEASURE_MAP but raises NotImplementedError.
- **Goal:** Extract T566 = (∂²l/∂δ²)/2 from the COSY DA polynomial map
  and use it as a FIT objective. Requires `transfer_matrix_order >= 2`.
- **Use case:** Bunch compression optimization where both R56 and T566 matter.

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

### C3. Field Map Scaling for MGE Dipoles [FIX APPLIED 2026-02-22]
- **Root cause:** The Mathematica notebook (`fields/calculation/UH_chicane_fringe.nb`)
  applied an erroneous 0.835× momentum scaling (P/P_45) during fieldmap generation,
  and set DELTAS=0.0 instead of 0.001 m.
- **Fix:** DELTAS corrected to 0.001 m; all 201 field values multiplied by
  1/0.8351818473537908. Peak field now 0.5307 T, matching the source CSV
  (`fields/calculation/UH_chicane_fringe.csv`). Fieldmap length = 200 × 0.001 = 0.2 m.
- **TODO:** Re-run COSY optimization with MGE (FR 3) using the corrected fieldmap.
  This is the most physically accurate dipole model and was previously blocked.
- **Files:** `fields/chicane_dipole_fieldmap.dat`, `fields/calculation/chicane_dipole_fieldmap.dat`,
  `cosySimulator.py:113-117` (MGE parameters), `cosySimulator.py:1254-1279` (MGE FOX generation)

### C4. Systematic Testing, Validation & Verification [HIGH PRIORITY]
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

### C2. COSY INFINITY Cross-Validation [DONE — see W4]
- **Motivation:** Independent DA-based simulation. Particularly valuable for
  verifying dispersion and chromatic behavior at high σ_E.
- **Design:** Run COSY optimisation directly (W4), then compare element-by-element
  Twiss functions with the Python results.
- **Due:** 2026-02-18 (W4)

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

## Notes

- All optimization scripts use `seed=42` for reproducibility
- NewFELsim conda environment required: `/home/evaletov/.conda/envs/NewFELsim/bin/python`
- Run commands from `backend/` with `MPLBACKEND=Agg PYTHONPATH=$(pwd)`
