# UH MkV FEL Beamline Optimization — Priorities & Roadmap

**Date:** 2026-02-11 (updated 2026-02-12)
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

### O3. Evolutionary Optimisation [MEDIUM PRIORITY]
- **Motivation:** Nelder-Mead is a local optimizer; at extreme parameters, multiple
  solution families exist (e.g., the ε_n=14–16 vs ε_n=18–20 quad patterns).
  Evolutionary optimisation (EO) can explore multiple basins simultaneously.
- **Options:**
  - User has a highly effective custom EO optimizer to integrate later
  - scipy `differential_evolution` as interim (v2 exploration showed it matches NM
    for the 4-variable problem, but at 10× cost)
- **Additional consideration:** optimise more quad currents simultaneously (beyond
  the current 4 Stage-11 variables). EO handles high-dimensional problems better
  than NM.

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

## Category C: Cross-Validation

### C1. RF-Track Cross-Validation [HIGH PRIORITY]
- **Motivation:** The FELsim optimizer uses transfer matrices; RF-Track uses
  full particle tracking. Comparing at key parameter points validates the
  transfer matrix model.
- **Design:** Run RF-Track at baseline + 3–4 extreme S4 points. Compare
  final Twiss at undulator entrance.
- **Prerequisite:** RF-Track adapter (`rftrackAdapter.py`) with lattice_path support

### C3. Field Map Scaling for MGE Dipoles [HIGH PRIORITY — BLOCKED on user input]
- **Status:** The chicane dipole fieldmap (`fields/chicane_dipole_fieldmap.dat`) has
  DELTAS=0.0 and field values (~5 mT) that are ~50× too small for the expected
  dipole fields (~250 mT). The `mge_scaling = P/P_45` only applies momentum
  scaling (~0.89), not the missing absolute scaling factor.
- **Action needed:** User to provide the scripts/files used to generate the field map
  so the correct scaling coefficient can be recovered. Once the fieldmap is fixed,
  re-run COSY optimization with MGE + FR 3 for the most physically accurate model.
- **Files:** `fields/chicane_dipole_fieldmap.dat`, `cosySimulator.py:113-117`
  (MGE parameters), `cosySimulator.py:1254-1279` (MGE FOX generation)

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
