# UH MkV FEL Beamline Optimization — Priorities & Roadmap

**Date:** 2026-02-11 (updated 2026-06-19)
**Scripts:** `backend/test/UHM_beamline_opt_*.py`

---

## Action items from 2026-05-04 meeting with Niels [NEW — git-bug umbrella]

Verbatim items captured 2026-05-29; cross-referenced to existing roadmap/git-bug.

### Linac
- [ ] **Xsuite — more detail.** Extend the xsuite linac model. ↔ L3.3 / git-bug `343b42f` (RF-Track ↔ xsuite linac SC comparison).
- [ ] **TW (travelling-wave).** COSY TW extension (N28 Phase 2, designed not coded — see `meeting_slides.txt`). RF-Track TW_Structure already done in L1.

### COSY space charge (`~/COSY/cosy-fmm/demo/spch_demo/`)
- [x] **Analyse as a function of equivalent macroparticle charge** (Q/N_p). DONE — F1: excess growth ∝ q_mp^0.46 (≈√, = 1/√N_p shot noise); threshold q_mp\*≈0.037 pC.
- [x] **Check SC charge-density profile + potential** — DONE — F5: ρ(r), Φ(r), E_r(r); exposed φ via `phi.dat` in `spch_kick.f90`; E_r matches analytic 2D Gaussian.
- [x] **Run at 1 MeV** — DONE (F6, matched). Matched beam on a STABLE FODO (1/f=1, μ=60°, 4 cells); the demo 1/f=2 cell is at μ=180° (no matched optics). 45 MeV @1 nC = 0.31% (matched) vs 1 MeV @1 nC = 4.6e4% (SC-limit blowup); ~7-order gap is pure SC since both start matched. Risk R1 resolved.
- [ ] **Check longitudinal SC effects.** ↔ L3.5 (was deferred — now active per Niels). Phase 4 (Fortran treecode edit).

### COSY ↔ RF-Track ↔ Xsuite three-code SC comparison (1 MeV and 45 MeV)
- [x] **16 ps bunch and 1 ps bunch** — COSY FODO done (F7): 1 ps=3.99%, 16 ps=0.58%. Three-code/line version = Phase 3.
- [x] **Small spot size** — COSY FODO done (F7): 0.25 mm=26%, 4 mm=0.01%. Three-code/line version = Phase 3.
- [x] **Enable SC in RF-Track and in Xsuite** — Xsuite frozen-Gaussian SC exercised in the capstone (`sc_capstone_run.py`); XsuiteAdapter has frozen+PIC. RF-Track PIC works on FODO; full-line PIC core-dumps (`af9d56c`).
- [ ] **Add BEAMPATH ("Beampass") simulation — for the injector** (confirmed by Niels 2026-05-29; NOT the linac). Batygin code; Niels already uses it for the injector (gun + α-magnet + matching quads, tickler 2026-04-05). Scope = injector cross-check.
- [ ] **Later: full PIC in COSY.** ↔ git-bug `1cbd8ea` (SC Option 1).
- [ ] **DA-FMM + RF-Track + Xsuite on the focusing/transport line, at 45 MeV and 1 MeV:** — groundwork DONE 2026-06-19, see `backend/test/sc_capstone/` + `results/sc_capstone/CAPSTONE_REPORT.md`. Common-distribution generator + reproducibility manifest, no-SC cross-code handoff, and the DA-FMM-vs-xsuite-frozen SC run all built + verified (11/11 tests).
  - [ ] **with dipoles** — BLOCKED: xsuite has no dipole edge/fringe model (full-line no-SC handoff: xsuite fails, FELsim-vs-RF-Track 18× on σ_x) + SC-inside-magnets slicing needed. The binding prerequisite.
  - [x] **without dipoles** — DONE: DA-FMM vs xsuite-frozen on section [32,46) at 45 & 1 MeV. Cold-beam agrees <0.1%; 1 MeV (fair, q_mp≪q_mp*) shows the frozen Bassetti-Erskine model over-predicts ε_n growth ~5× vs the N-body treecode in strong SC; 45 MeV high-charge DA-FMM excess re-demonstrates the Phase-1 q_mp* shot-noise threshold on a real section.

---

## Action items from early-June 2026 meeting with Niels [NEW — git-bug umbrella `1aa5c44`]

Captured 2026-06-19 from handwritten notes; cross-referenced to repo/git-bug/memory via a per-item investigation pass. **Target: complete by Mon 2026-06-22. Delivery: results emailed to Niels, off the regular Monday-meeting cadence.** Exact meeting date TBC (Monday cadence → 2026-06-01 or 2026-06-08; recalled "~June 3").

1. **elegant scope / channel** *(largely settled — confirm scope with Niels in the email).* Note: *"Don't need elegant for the diag[nostic] chicane — maybe for RF gun."* (Channel read by Eremey 2026-06-19 as the **diagnostic chicane**.) Reading: elegant is **not** needed for the diagnostic chicane — its optics are already covered by COSY / Xsuite / FELsim-1st-order, and elegant's chicane/CSR strengths aren't required at 40 MeV / this charge; elegant stays scoped to **RF acceleration**. *"maybe for RF gun"* = the gun (RF-Track-only today) is the candidate next elegant cross-check. Fact check: elegant is currently wired ONLY as the linac reference (RFCA+TWLA, `elegant_linac/`); no `ElegantAdapter` for the chicane/gun/α-magnet/transport. Action: confirm in the email + offer an elegant gun cross-check. ↔ git-bug `35e31e6` (L1; Phase 5 report + Phase 6 polish still open).
2. **Linac: codes × physics-feature table** *(in-progress / mostly new).* Existing P1 capability matrix (`generate_ipac_figures.py`, `results/ipac/P1_capability_matrix.*`) is codes × beamline-BLOCKS (0/1/2), not codes × PHYSICS-FEATURES, not linac-scoped, no implemented-vs-potential split. Build a new linac table (rows = codes, cols = TW/SW, field maps, SC, beam loading, wakefields/CSR, longitudinal; cells = implemented vs potential).
3. **Linac: more complex / refined model** *(in-progress — continuation of 2026-05-04 (a)+(b)).* xsuite linac is a lumped single `Cavity` (`compare_xsuite.py`) and `XsuiteAdapter` treats `RF_CAVITY` as a drift; COSY TW designed-not-coded; RF-Track `TW_Structure` is the detailed one. **First cell-resolved TW model built + validated (2026-06-19):** `backend/test/rftrack_linac/linac_multicell_tw.py` reads the production τ=0.57 geometry (84 cells + 2 couplers, SLAC-75 Table 6-6), autophases the 1 MeV electron (handles the β<1 slippage), and reproduces the elegant peak energy gain to **1.18%** (40.95 vs 41.44 MeV; residual = 3.009-vs-3.048 m fractional-cell tail), det(Rx)=0.0343, fill time 0.81 µs. **Beam loading added (2026-06-20):** steady-state CG model (k(0.57)=0.232, 37 MeV/A, ~3% droop at 32 mA), loaded-gradient sag toward output, droop matches closed-form I·r·L·k to <0.01 MeV; analytic vg ratio e^{-2τ}=0.32 confirms aperture-derived 0.319 (`linac_beam_loading.md`). **xtrack-native element DONE (2026-06-20):** `xsuite_linac_tw.py` drives an xtrack Line via per-cell Cavity + ReferenceEnergyIncrease (energy-ramp reference); on-axis K_out matches the integrator to 0.000% (δ~1e-16), transverse det(R_x)=0.0343=p_in/p_out to 0.00% — a real xtrack Line ready for the xsuite SC engines. **XsuiteAdapter integration DONE (2026-06-20):** `RF_CAVITY` now builds the multi-cell TW (per-cell Cavity + ReferenceEnergyIncrease) instead of a drift, with energy-aware quad k1 (running energy threaded through `_build_line`); `slac_linac.json` accelerates 1→41.4 MeV (δ~1e-15, ~0.06% vs elegant/RF-Track), transport-only lines unchanged (`test_xsuite_linac_adapter.py`, `xsuite_adapter_linac.md`). **Per-cell SC inside the linac DONE (2026-06-20):** the adapter SC path interleaves one frozen/PIC SC kick per TW cell, each after the cell's reference ramp, so each acts at the local energy; validated (`validate_sc_linac.py`): 87 SC kicks one-per-cell, a 4000-particle bunch reaches 41.4 MeV with all surviving, and the per-cell SC defocus drops 24265× from 1→41.5 MeV (matches 1/(βγ)³, predicted 25782×) — `sc_linac_validation.md`. **σ(s) envelope INVESTIGATED (2026-06-20):** the SC σ(s) is already self-consistent — for a tracked bunch the frozen SC runs `update_on_track`, recomputing σ per cell from the actual beam (not the fixed injection σ). The beam stays ~1 mm here (focusing-free section), so a naive adiabatic 1/√(βγ) law would be WRONG. Added an OPTIONAL prescribed-envelope path (`sig_env` + `sc_envelope_prepass`) for matched focusing lattices / deterministic+PIC-mesh use, but it is NOT the default: a fixed σ over-defocuses an unmatched beam (ε_n 5.8 vs 2.8) because it does not self-correct (`sc_envelope_validation.md`). **simulate() exit-energy frame DONE (2026-06-20):** the output coordinate transform now uses the EXIT reference energy (read from `p.energy0`), not the injection energy — so an accelerating line returns correct FELsim δ/angles (linac on-axis output |δ|≈0.004, was ~40; transport lines unchanged). **Prescribed-envelope regime characterized (2026-06-20):** charge scan on the linac shows the prescribed path agrees with self-consistent at low charge (0% at 1 pC) and diverges as SC grows (16% at 0.1 nC, 551% at 1 nC) — valid for stable/matched beams only; multi-cavity `sig_env` generalized + verified for contiguous cavities (`validate_sc_envelope_regime.py`; interspersed non-cavity SC = future). The full matched-FODO focusing test is deferred (a quick synthetic FODO was not stable). Remaining: COSY TW DA map (the third linac code); full matched-focusing accelerating lattice on the real injector+linac. ↔ `343b42f`, `35e31e6`, `f59e144`.
4. **PALS beam-file standard?** *(not-started).* **Answer: NO** — PALS is a lattice standard only (whole submission package + I2/I3/I8 are lattice-only). Beam/particle-distribution interchange is undecided: git-bug `8d0d833` (P3) still open; FELsim writes only COSY rray/ascii (`cosyParticleSimulator.py`), no openPMD/SDDS/HDF5. Decide a format (recommend **openPMD-beamphysics**, ASCII fallback); ax3l/Axel Huebl (openPMD, a PALS participant) reachable via PALS Issue #176. Investigate Genesis4 native particle-input path (gates the choice).
5. **Fringe-field (FF) mode comparison across codes** *(largely-done).* FF = fringe field. Already covered by R2 (FELsim / COSY FR0 / COSY FR1 / RF-Track) + W4 (FR0/FR1/FR3 convergence) + P8 (FR0-vs-FR3 across DA orders) + P11 (FELsim triangle-φ impact); tooling `--fr {0,1,2,3} --warm-start`, 18 `cosy_s1_fr*` JSONs. Remaining: add xsuite (no-fringe baseline — adapter treats edges as drift), **refresh the RF-Track column** (triangle-φ correction has since been added → R2-era β_y deficit narrative is stale), one consolidated full-lattice FR0–FR3 table. No git-bug yet.
6. **Cross-calibration (between codes)** *(in-progress).* Backbone DONE: single anchor G=2.694 T/A/m, identical k1 formula in FELsim/RF-Track/xsuite, COSY same gradient; δ-convention conversions in `CoordinateTransformer` (COSY ΔK/K₀ vs xsuite Δp/p₀), verified vs MAD-X ~1e-10; xsuite tests pass. Remaining (= live ask): the three-code SC-on transport-line capstone at 45 & 1 MeV ±dipoles; xsuite still lacks dipole-edge/fringe + RF models. ↔ `f59e144` item l, `343b42f`.
7. **FEL objectives A/B — continue iterations** *(in-progress).* Ablation done twice: `results/ablation/` (A 10% / B 15% / C 75% NM-fail) and corrected `results/ablation_MOP6318/` (A 35% / B 35% / C 80%). "Aware of diff between objectives" = at MOP6318 targets A and B BOTH fail ~35% (not robust as the abstract framed) → narrative shift to confirm (git-bug `cf11eb7`). Remaining: **S6 BO baseline** (P0 paper deliverable; needs Niels's xopt hyperparameters — requested in `email_niels_ipac.txt`, confirm if replied), write the MOP6318 `summary.md` interpretation (placeholder), reconcile undulator-Twiss provenance (β_x=1.267/α_x=0.560 from a code comment), S4/S5/S7. ↔ `7f690aa`, `cf11eb7`.
8. **Try optimization with glyfada** *(largely-done).* Already tried + characterised: W6 (fails at ε_n=5/8/14, NM beats by 3–6 orders), W7 (conclusion = two-phase NM→pycma-CMA-ES; distributed glyfada not needed for this 4-var narrow-basin problem), O4 (26D Koa ~27k evals ~99% penalty), O5 (pycma adopted into production). Remaining (pick target with Niels): if the IPAC objective matrix, run glyfada on **Config C** (stress case) + report head-to-head with NM/BO (overlaps `7f690aa` S4/S6); likely a negative-result-report deliverable. Use the `NewFELsim` conda env (system python3 can't import glyfada).
9. **COSY DA-FMM SC + PIC — continue work** *(in-progress).* Phase 1 DONE 2026-05-29 (FODO physics F1–F7; XsuiteAdapter enabler). **Phase-3 groundwork DONE 2026-06-19** (`backend/test/sc_capstone/`, `results/sc_capstone/CAPSTONE_REPORT.md`): (i) common-distribution generator + reproducibility manifest; (ii) no-SC cross-code handoff — codes agree <0.1% on the no-dipole section [32,46), full line NOT comparable (xsuite fails on dipoles, FELsim-vs-RF-Track 18×) = sequencing trap quantified; (iii) DA-FMM-vs-xsuite-frozen SC on that section at 45 & 1 MeV (1 MeV fair: frozen over-predicts ~5×; 45 MeV high-charge = q_mp* shot noise). 11/11 tests. Remaining: **with-dipoles** needs xsuite dipole edge/fringe + SC-inside-magnets slicing (binding); RF-Track 2.5.5 full-line PIC core-dump (`af9d56c` / Risk R5); DA kick-composition (`fmm_eval_treecode_da`, Niels's route); Phase 2 linac TW; Phase 4 longitudinal SC, BEAMPATH injector (LANL access), full PIC (cosy-pic Phase 3 native SPCKICK). ↔ `f59e144`, `35e03f2`, `1cbd8ea`, `343b42f`.

---

## Due by 2026-04-13 (Monday — next meeting with Niels)

### L1. RF-Track Linac Model vs elegant Benchmark [IN PROGRESS — Phase 0-2 DONE 2026-04-05]
- **Git-bug:** 35e31e6 (P1-high) — Big item for Monday
- **Scope:** RF-Track TW_Structure model of SLAC 3-m S-band linac, cell-by-cell
  geometry from SLAC-two-mile-report10.pdf, benchmark vs elegant RFCA+TWLA
- **Standalone:** `backend/test/rftrack_linac/linac_standalone.py` — 73-point
  phase scan, peak K_out=41.468 MeV at phid=0 (autophased)
- **Benchmark headline:** **0.06% agreement at peak** (RFT 41.468 vs elegant
  RFCA-optimal 41.442 MeV at 1 MeV injection, 13.3 MV/m, 3.048 m)
- **Phase 0 DONE:** RF-Track API reconnaissance, directory scaffolding,
  phase-convention understanding (autophase → phid=0 ≡ optimal phase)
- **Phase 1 DONE:** Standalone script, 73-point scan, CSV + PDF output,
  loss handling (11/73 points in deep-decel zones phid ∈ [+85°, +155°])
- **Phase 2 DONE:** `rfCavityLattice` class in beamline.py; RFC branch in
  latticeLoaderBase._element_to_object; _build_rf_cavity() in rftrackAdapter
  dispatching to TW/SW/RFCA; linac-only JSON `var/slac_linac.json`;
  integration test passes (0.005% standalone↔adapter agreement); no
  regressions (32 smoke + 22 JSON/YAML equivalence tests pass)
- **Phase 3 DONE (2026-04-05):** `benchmark_vs_elegant.py` adapter-path
  benchmark producing 3 figures in `reports/2026/Apr/13/figures/`:
  phase_vs_Eout.pdf, detRx_vs_phase.pdf, twiss_evolution.pdf. Agreement:
  energy gain 0.063% at peak, det(R_x) 5.7% at peak, exit Twiss β_x 2.8%
  / α_x 1.4%. R-matrix extraction via 3-particle unit-perturbation;
  Twiss evolution via 30 progressively-longer sub-lattices with their
  own autophase (approximation, small systematic error)
- **Phase 4 DONE (2026-04-05):** Cell data extracted from
  SLAC-two-mile-report10.pdf Tables I+II (86 cells, all 4 columns,
  both disk thicknesses). **Critical finding:** these tables are for a
  tau=0.40 prototype that was REJECTED — production SLAC structure uses
  tau=0.57. Tables I+II are the tau=0.40 prototype; archival reference in
  `docs/sbend_linac/slac_cells_tau040_table{I,II}.csv`. Our analytical
  single-Fourier-coefficient model is the CORRECT approach for the
  constant-gradient TW benchmark.
  - **CORRECTION (2026-06-19):** the earlier claim that production tau=0.57
    cell data "is not in our documents" is STALE. The production tau=0.57
    cell-by-cell geometry WAS located the next day (2026-04-06) in a different
    source — SLAC-75 "The Stanford Two-Mile Accelerator," Table 6-6 (p. 145) —
    and extracted to `docs/sbend_linac/slac_cells_tau057_table66.csv` (100
    cells, 2a/2b, 2pi/3 mode, rho=0.1215 in, t=0.23 in), cross-checked vs
    SLAC-75 Fig 6-22 (`plot_cell_geometry.py`). So a cell-by-cell multi-cell TW
    model is a BUILD task, not blocked on external input. A measured/simulated
    on-axis field map would be a later refinement (derivable from the geometry
    via constant-gradient relations).
- **Remaining:** Phase 5 (LaTeX report at reports/2026/Apr/13/),
  Phase 6 (polish)
- **Open TODOs (from 2026-04-20 review):**
  - **n_cells fractional:** initial coupler hypothesis was wrong
    (tested empirically 2026-04-20). The 0.11-cell tail is only
    3.89 mm, too small to be 2 couplers (~35 mm each). Real cause:
    (i) 3.048 m is nominal/round (10 ft) vs. synchronous 87 × 34.99
    = 3.0441 m, and (ii) production tau=0.57 is constant-gradient
    with slightly varying cell lengths along the structure — our
    uniform-L_cell_sync model absorbs both into a fractional tail.
    Document in the Phase 5 report; no code change needed.
  - **β_x = β_y = 1 [RESOLVED 2026-04-20]:** not a flat line —
    both codes evolve from β=1 m → ~14 m. The "1" is just the
    placeholder initial Twiss (`linac_twiss.ele` line 17 and
    `rft_twiss_evolution` default). Substitute RF-gun-matched
    injector β once Niels's repo run finishes; benchmark script
    now accepts `--beta0 / --alpha0` overrides.
  - **Niels's repo reviewed [DONE 2026-04-20]:** identified
    `ode_epsabs` as the "to.l." tolerance knob; `dt_mm` equivalent
    on Structure is `set_nsteps`. Documented in `linac_model.md`.
  - **Δt/to.l. convergence study [DONE 2026-04-20]:**
    `backend/test/rftrack_linac/convergence_study.py`. K_out
    converges ~1/n² under RK2; default nsteps=872 is conservative
    (nsteps=200 gives ~8e-6 MeV error at 0.012 s vs 0.033 s).
    `epsabs` has negligible effect above 1e-5. All GSL integrators
    agree to 1e-5 MeV — no reason to switch from rk2.
  - **β_x = β_y = 1 in Twiss plot:** investigate why the adapter-path
    Twiss evolution flattens at β=1. Suspected uninitialised Twiss
    input (default identity) vs. a physical injector β; confirm
    whether `Bunch6d`/sub-lattice autophase wipes the initial Twiss,
    and pass a matched β from the injector design instead.
  - **Niels's GH repo:** review Ioniels/RF-TRACK_UH_ThermionicRFgun
    modelling patterns (field-map handling, autophase workflow,
    COMSOL integration) for ideas applicable to our linac model.
  - **RF-Track convergence study:** profile how quickly RF-Track
    converges as integration step Δt is refined. Use Δt as a
    metaparameter (sweep, compare K_out / Twiss / R-matrix). Also
    clarify `to.l.` option from the manual — likely the length-step
    / longitudinal tolerance control — and sweep it in the same way.

### L3. Linac + COSY space-charge follow-ups [DONE 2026-05-04 — see reports/2026/May/04/L3_sc_studies.pdf, 14 pages, includes L3.1-L3.7]
The "COSY SC demo" referenced here is **`~/COSY/cosy-fmm/demo/spch_demo/`**
(DA-FMM kicks driven from a COSY INFINITY FOX program via Fortran-bridge
file handshake; 45 MeV / 1 nC / 1k macroparticles / 2 m FODO; presented
2026-04-22). It is the shipped Demo-A-flavored deliverable, not Demo B
(0a70783, analytical linear Gaussian, still TODO).

- **L3.1 — N_p convergence (spch_demo).** Sweep
  N_p ∈ {500, 1k, 2k, 5k, 10k, 20k, 50k} at fixed
  N_slice=20, Q=1 nC, θ=0.3. Diagnostics: ε_n,x/y(2 m), σ_x/y(s),
  wall-time. Identify N_p* where moments stabilise <1% relative.
- **L3.2 — N_slice convergence (spch_demo).** Sweep
  N_slice ∈ {5, 10, 20, 40, 80, 160, 320} at N_p=N_p*, Q=1 nC, θ=0.3.
  Verify second-order convergence in Δs. Identify N_slice*. Bonus axis:
  θ ∈ {0.1, 0.3, 0.5} (DA-FMM MAC).
- **L3.3 — RF-Track ↔ xsuite linac comparison.** Git-bug 343b42f
  (P2-medium). Same input bunch + cavity-chain xsuite approximation of
  the SLAC TW structure; compare moments along s. Distinct from the
  COSY spch_demo work — this is the linac SC validation axis.
- **L3.4 — Bunch charge sweep (spch_demo).** Q logarithmically spaced:
  Q ∈ {0.1, 0.3, 1, 3, 10, 30} nC at converged (N_p*, N_slice*).
  Watch ε_n growth ∝ Q in the linear regime; identify the charge at
  which DA-FMM nonlinear corrections become non-negligible. Maps to
  README extension #1 in spch_demo.
- **L3.5 — Longitudinal SC.** Deferred. Re-evaluate after L3.1–L3.4
  if longitudinal coupling shows up in σ_z growth or energy chirp drift.

### L2. COSY SC interspersed transfer-map demo (Demo B) [DONE 2026-05-04]
- **Git-bug:** 0a70783 (P2-medium) — closed
- **Implementation:** `~/COSY/cosy-fmm/demo/demob_gaussian/demob_gaussian.fox`
  (pure FOX, linearised Bassetti-Erskine kick, no DA-FMM, no Fortran bridge)
- **Result:** RMS emittance preserved to 7 sig figs at all Q ∈ {0.1…30} nC
  (linear forces ⇒ symplectic ⇒ second moments invariant). Cross-validates
  with L3.6/L3.7: the ~0.5% smooth-field growth in xsuite/L3.7-extrap is
  the nonlinear-BE contribution that this demo deliberately drops.
- **Write-up:** L2 section in `reports/2026/May/04/L3_sc_studies.pdf` (15 pages now).

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
- **Re-run with CMA-ES polishing (2026-03-11):** 20 emittance points re-run with
  current code (post-QA fixes) + CMA-ES polishing. Key changes vs original (Feb 15):
  - **ε_n = 2:** Failed → Marginal (5.9× improvement)
  - **ε_n = 5:** Failed → Excellent (114× improvement — the NM trap is resolved)
  - **ε_n = 14, 16:** Excellent → Acceptable (regressions from QA code fixes
    changing the objective landscape, not CMA-ES)
  - Quality distribution: 14 Excellent, 3 Acceptable, 3 Marginal, **0 Failed**
    (was 15E, 1A, 2M, 2F). Both Failed points eliminated.
  - Backup: `results/params_05ps/pre_cmaes_backup/scan_emittance_w2.csv`

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

### W7. Glyfada Config Optimization & Re-Benchmark [DONE 2026-03-11]
- Re-benchmarks Glyfada against NM with corrected CMA-ES configurations.
- **Bugs found in original W7 (2026-03-11):**
  1. `cma_es.initial_mean` was never set — CMA-ES used population centroid
     instead of NM warm-start point
  2. `feasibility_rules` constraint handling silently ignored for CMA-ES
     (only works with NSGA-II). Changed to `death_penalty`.
  3. Evaluator constraint was binary (0/1) — changed to continuous for
     gradient information
- **Re-benchmark results (2026-03-11):** pycma CMA-ES warm-started from NM
  dramatically outperforms NM alone at all emittance points:
  - ε_n=5: pycma MSE=4.5e-16 vs NM 8.2e-06 (~10¹⁰× better, 11s vs 32s)
  - ε_n=8: pycma MSE=1.2e-15 vs NM 6.1e-06 (~10⁹× better, 11s vs 33s)
  - ε_n=14: pycma MSE=6.3e-03 vs NM 8.9e-03 (~1.4× better, 25s vs 30s)
  - Glyfada binary configs (G-A/B/C) not tested — binary not built locally.
    Key result: CMA-ES algorithm is highly effective; glyfada's distributed
    infrastructure is not needed for this 4-variable problem.
- **Conclusion:** Two-phase NM→CMA-ES is the recommended Stage 11 strategy.
  NM finds the right basin; CMA-ES polishes to near-machine-precision.
- glyfadaAdapter rewritten to use `glyfada.optimize()` Python API (was subprocess+CSV)
- glyfada_eval: continuous constraint values, multi-objective support
- Results: `results/params_05ps/W7/`
- Script: `W7_glyfada_rebenchmark.py`

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
- Quality thresholds (MSE): <1e-3 excellent, <0.01 acceptable, <0.1 marginal
  (RMS equivalents: <3.2e-2, <1e-1, <3.2e-1)

### O-Abl. Objective-design ablation for Niels IPAC paper [DONE 2026-05-10]
- Scripts: `ablation_run.py`, `ablation_analyze.py`, `diagnose_seed91.py`
- Results: `results/ablation/` (60 random + 3 Sobol JSONs, plots, summary.md)
- Configs: **A** verbatim, **B** mild rescale (dispersion ref=0.5 m), **C** B + Stage 1 weight=0→0.5 typo fix + Stage 7 envelope 0→1.5 mm
- 20 seeds + 1 Sobol per config; ~5 min on starfield01
- Headline: A fails 10% (2/20), B fails 15% (3/20), **C fails 75% (15/20)**
- Failure mechanism: Stage 11 NM gets stuck in local min after Stage 10 shifts the starting basin
- Sobol deterministic beam also lands in the bad basin for C (consistent with multi-seed result, not an independent confirmation)
- Implication for IPAC paper: Config C is a useful stress case for global-search optimisers (BO, CMA-ES). Whether BO actually recovers NM's 75% failure rate is the open empirical question; A is the correct NM baseline.
- Open follow-ups: S4 CMA-ES drop-in, S5 joint 26-var, S6 BO baseline, S7 warm-chain starting points

---

## Category S: Parameter Studies

### S5. 0.5 ps 2D Coupled Scans [DONE 2026-03-11]
- **Motivation:** 1D scans hold other parameters fixed; realistic operation involves
  correlated changes (e.g., shorter bunch → larger energy spread). 2D scans map
  the feasibility surface.
- **Scans:**
  - S5a: (σ_E, h) grid at ε_n=8 — energy spread vs chirp coupling
  - S5b: (σ_E, ε_n) grid at h=5e9 — degradation interaction
  - S5c: (h, ε_n) grid — chirp compensation vs emittance
- **Design:** 10×10 grids (100 points each), 500 particles. Checkpoint/resume via
  CSV. Contour/heatmap plots (MSE LogNorm, Twiss deviation, feasibility bands).
  CLI: `--s5a/--s5b/--s5c/--all/--plots-only/--grid N`.
- **Results (2026-03-11):**
  - S5a: 100/100, **0 failures**, 92 Excellent, 11 Acceptable, 2 Marginal.
    Energy spread and chirp do not compromise matching at ε_n=8.
  - S5b: 100/100, **18 failures** at extreme emittances (ε_n<5 or ε_n>17),
    38 Excellent, 27 Acceptable, 17 Marginal. Consistent with W2 findings.
  - S5c: 100/100, **20 failures** at extreme emittances, 54 Excellent,
    12 Acceptable, 14 Marginal. Chirp does not rescue high-emittance failure.
  - **Key finding:** The feasibility boundary is dominated by emittance —
    σ_E and h have negligible impact on matchability at baseline ε_n=8.
- Script: `S5_2d_parameter_scans.py`
- **Output:** `results/params_05ps_2d/`, contour plots, feasibility boundary curves
- **Prerequisite:** S4 results to identify interesting regions

### S6. Bunch Length Sensitivity (0.1–2 ps) [DONE 2026-03-11]
- **Motivation:** The FELsim request asks about 0.5 ps specifically, but understanding
  the full bunch length range is valuable context.
- **Design:** Sweep bunch_spread from 0.1 to 2.0 ps (15 points) at two parameter
  sets: (a) baseline (σ_E=0.5%, h=5e9) and (b) emittance-conservation scaled.
- **Results:**
  - Baseline: 15/15 Excellent. Bunch length has **no effect** on transverse matching.
    Confirms S9's linear decoupling prediction: transverse Twiss depends only on
    the 4×4 block + dispersion column, not the longitudinal (column 5) distribution.
  - Emittance-conserved (σ_E=2%, h=20e9): 11/15 Excellent, 3 Acceptable, 1 Marginal.
    Dips are NM noise (not bunch-length-correlated), consistent with S4's finding
    that σ_E=2% introduces some optimizer sensitivity.
  - **Conclusion:** No bunch length threshold exists — the optimizer does not degrade.
- Script: `S6_bunch_length_sensitivity.py`
- Results: `results/S6/`

### S7. Verification Runs at Key Points [DONE 2026-03-11]
- **Motivation:** The 500-particle sweeps trade accuracy for speed. Key points
  (boundaries, transitions) need 1000–2000 particle confirmation.
- **Design:** From S4 results, re-run 7 emittance points (ε_n = 1,3,5,8,14,16,20)
  and 3 energy spread points (σ_E = 0.4,0.55,0.7%) at N=500, 1000, 2000.
- **Results:**
  - **Energy spread: fully consistent.** All 9 runs (3 points × 3 particle counts)
    are Excellent. The σ_E=0.55% Acceptable dip in S4 was a statistical artifact
    (different seed or NM basin).
  - **Emittance: highly inconsistent at extremes.** Only ε_n=8 (baseline) is
    consistent across particle counts. At ε_n=1,3,5,14,16,20, quality classification
    varies wildly between N=500/1000/2000 — same emittance, different particle
    samples → different NM local minima.
  - **Key finding:** S4/S5 emittance results at extreme values are **not robust** —
    they depend on the specific random particle realization, not solely on physics.
    The optimizer landscape has multiple local minima that the beam sample selects.
    Publication-quality claims require multi-start + multi-seed statistics.
  - ε_n=14,16: N=2000 achieves Excellent (MSE~2e-5) while N=500,1000 give
    Acceptable/Marginal — suggests the problem is NM basin selection, not physics.
- Script: `S7_verification_runs.py`
- Results: `results/S7/`

### S8. Multi-Start Robustness Study [DONE 2026-03-11]
- **Motivation:** S7 showed that single-start emittance results at extremes are not
  statistically robust. This study directly characterizes the optimizer landscape.
- **Design:** 5 emittance points (ε_n = 1, 3, 5, 14, 16) × 10 random Stage 11 starts.
  Reports: best/worst/median RMS, quality classification histogram, box-and-whisker plot.
- **Implementation:** `S8_multistart_robustness.py`. Uses `stage11_startPoint` to inject
  random starting currents for the 4 Stage 11 quads (q87, q93, q95, q97). Checkpoint/resume
  via CSV. Runtime ~17 min.
- **Key finding: Stage 11 is unimodal.** All 10 random starts at each ε_n converge to
  exactly the same solution (identical MSE, quad currents, and Twiss to machine precision).
  The Stage 11 objective landscape has a single global basin — multi-start cannot improve it.
- **Consequence:** The S7 seed-dependence at extreme emittances originates entirely from
  stages 1–10 (random beam generation affecting upstream matching), not from Stage 11's
  Nelder-Mead search. Improving extreme-emittance robustness requires addressing the
  upstream stages (e.g., multi-seed averaging, deterministic beam generation).
- **Results:**

  | ε_n | RMS | Quality | β_y (m) |
  |-----|-----|---------|---------|
  | 1 | 1.08e-1 | Marginal | ~0 |
  | 3 | 138.8 | Failed | 264 |
  | 5 | 5.58e-3 | Excellent | 0.242 |
  | 14 | 9.53e-2 | Acceptable | 0.030 |
  | 16 | 1.06e-1 | Marginal | 0.005 |

- **Results dir:** `results/S8/`

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

### O1. Warm-Starting from Neighboring Points [DONE 2026-03-11 — negative result]
- **Motivation:** Each scan point starts from the same fixed initial guess. Using
  the optimized currents from the previous (neighboring) scan point as the start
  could improve convergence speed and robustness at extreme parameter values.
- **Implementation:** `warm_start_currents` parameter added to `run_optimization()`,
  `warm_start=True` flag added to `run_scan()`. When enabled, each point passes the
  previous point's `quad_currents` result as starting guess for all 11 stages
  (clamped to bounds). Failed points do not propagate warm-start.
- **Validation:** `O1_warm_start_validation.py` — re-runs S4 emittance sweep with
  warm_start=True and compares RMS, convergence speed, and current continuity.
- **Result: Naive sequential warm-starting is counterproductive.**
  Warm start wins only 4/20 emittance points. The ε_n=1 basin (Marginal, RMS≈0.108)
  traps the optimizer: once locked into those currents, stages 1–10 propagate the
  Marginal solution through ε_n=2–13. Cold start independently finds better basins
  at each point.
- **Possible improvements (not yet implemented):**
  - Run both cold and warm starts, keep the best
  - Warm-start only Stage 11 (not stages 1–10)
  - Bidirectional scan (low→high and high→low, merge best)
- **Results dir:** `results/O1/`

### O5. CMA-ES Polishing in Production Optimizer [DONE 2026-03-11]
- **Motivation:** W7 showed pycma CMA-ES warm-started from NM dramatically outperforms
  NM alone. Integrated as automatic post-NM polishing step in `run_optimization()`.
- **Implementation:** After Stage 11 NM multi-start finds the best result, pycma
  CMA-ES refines from that point (σ=0.1, popsize=20, maxfevals=3000). Falls back
  gracefully if pycma not installed (`except ImportError: pass`).
- **Validation (ε_n=5,8,14, 5 restarts):**
  - ε_n=5: MSE 8.18e-06 → 7.76e-06 (1.1× better)
  - ε_n=8: MSE 6.13e-06 → 6.13e-06 (dispersion-limited, NM already optimal)
  - ε_n=14: MSE 8.93e-03 → 5.05e-03 (1.8× better)
- **Note:** The composite 5-goal objective (4 Twiss + dispersion) has a dispersion
  floor ~2.5e-6 that stage 11 quads cannot reduce. Pure 4-Twiss MSE reaches 1e-15
  (see W7). CMA-ES provides most value at difficult emittances (ε_n=14).
- **Full re-scan validation (2026-03-11):** All parameter scans (S4, S5, S7) re-run
  with CMA-ES polishing active. Key results:
  - **S7 verification:** 8/21 points improved (up to 31.6×), 0 regressions. Notable:
    ε_n=1,3 at N=1000 improved 25–32×; ε_n=14 at N=500,1000 improved 1.4–1.9×.
  - **S5a (σ_E × h):** 18/100 improved (13 by >2×), 0 worse.
    Quality: 87→96 Excellent, 11→3 Acceptable, 2→1 Marginal.
  - **S5b (σ_E × ε_n):** 50/100 improved (34 by >2×), 0 worse.
    Quality: 38→60 Excellent, 27→19 Acceptable, 17→7 Marginal, 18→14 Failed.
  - **S5c (h × ε_n):** 27/100 improved (18 by >2×), 0 worse.
    Quality: 54→65 Excellent, 12→11 Acceptable, 14→8 Marginal, 20→16 Failed.
  - **S4 emittance:** All identical (confirms Stage 11 unimodality with fixed seed).
  - **S4 energy spread:** 7/15 better, 2 worse (NM basin variability, all Excellent).
  - **S4 chirp:** 5/12 better, 4 worse (NM basin variability, all Excellent).
  - **Conclusion:** CMA-ES polishing provides zero regressions and substantial
    improvements in the emittance-varying parameter space. The largest gains are
    in 2D scans involving emittance (S5b: 50% of points improved). The S4 1D
    emittance scan is fully deterministic (no change), confirming S8's unimodality
    finding. Pre-CMA-ES results backed up to `pre_cmaes_backup/` directories.

### O2. Adaptive Scan Resolution [DONE 2026-03-11]
- **Motivation:** Uniform spacing wastes points in flat regions and undersamples
  transition regions.
- **Implementation:** `run_adaptive_scan()` in `UHM_beamline_opt_05ps_params.py`.
  Starts with a coarse uniform grid (default 5 pts), computes RMS slope between
  adjacent points, inserts midpoints where slope exceeds `slope_multiplier × median`
  (default 1.5×). Iterates up to `max_iter` passes or `max_points` total.
  Supports warm-start (uses nearest evaluated neighbor's quad currents).
- **Parameters:** `param_range=(lo,hi)`, `n_coarse=5`, `max_points=30`,
  `max_iter=4`, `slope_multiplier=1.5`.

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

### R1. Interactive Parameter Explorer [IMPLEMENTED]
- **Motivation:** Static EPS plots are good for papers; interactive plots help
  with exploration.
- **Design:** Plotly/Dash dashboard reading CSV data. Two tabs: S4 1D scans
  (RMS + Twiss subplots with hover showing quad currents), S5 2D heatmaps
  (log₁₀ RMS with threshold contours). IQR-based robust y-axis limits for
  outlier handling.
- **Script:** `R1_parameter_explorer.py`. Run: `python R1_parameter_explorer.py`
  → opens at http://localhost:8050.
- **Dependencies:** `pip install dash plotly`

### R2. Comparison Table Across All Studies [DONE 2026-03-02]
- Aggregates data from W4, S4, W8, W9, W10, W11, W12 into 5 cross-code tables
  and 3 summary plots.
- **Tables:** (1) Baseline cross-code optimization, (2) Parameter sensitivity summary,
  (3) Bunch length & transmission, (4) Compression feasibility, (5) Quad currents.
- **Plots:** 3-panel RMS vs parameter, cross-code Twiss bar chart, compression curve.
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

### I6. MCNP-Style Robustness & Foolproofness [DONE 2026-03-11]
- **Motivation:** MCNP is a gold standard for production code robustness: every
  input is validated, edge cases are caught with clear diagnostics, defaults are
  sensible, and the code never silently produces wrong results. FELsim should
  adopt this level of rigour.
- **Actions:**
  1. Input validation at system boundaries: lattice files, API payloads,
     CLI arguments, Excel data. Fail loudly with descriptive errors. [DONE]
  2. Guard against silent numerical failures: NaN/Inf propagation,
     singular matrices, zero-length elements, particle loss without warning. [DONE]
  3. Consistent error handling: no bare `except:`, no swallowed exceptions.
     Every failure path either recovers correctly or raises with context. [DONE]
  4. Default values must be physically sensible (not 0 or 1 by convenience). [DONE]
  5. Audit all `setattr`/`getattr` patterns for typo-resilience (consider
     `__slots__` or property validation on beamline element classes). [DONE]
  6. Configuration validation: warn on unused/unknown keys, reject
     contradictory settings. [DONE]
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
- **Progress (2026-03-11):**
  - **`__slots__` enforcement (item 5):** Added `__slots__` to all beamline
    element classes (`lattice`, `driftLattice`, `qpfLattice`, `qpdLattice`,
    `dipole`, `dipole_wedge`, `fringeField`). Typo attributes now raise
    `AttributeError` at runtime. 6 new runtime enforcement tests in
    `test_attribute_guard.py`.
  - **Configuration validation (item 6):** `SimulatorFactory.create()` now
    validates kwargs against per-simulator known key sets. Warns on unknown
    keys and on features passed to unsupported simulators (e.g., `space_charge`
    on FELsim). `MultiCodeSimulator._init_simulators()` validates per-section
    config keys against known set.
  - **Edge case test suite:** 48 new tests in `test_edge_cases.py` covering
    degenerate distributions, aperture edge cases, matrix composition precision,
    extreme parameters, chromatic tracking, and input validation.
  - **Current test suite:** 278+ pass, 0 fail.

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
  - Test suite (test_multicode.py): 38 tests — SimSection, init, coord
    roundtrips (all 3 pairs), FELsim split equivalence (2/3-section),
    element conversion (drift, quad, DPW params), factory registration,
    FELsim→RF-Track hybrid (successful run, cross-validation, DPW params,
    physical apertures config), per-section config (runtime key filtering,
    separate vs shared instances), COSY adapter set_beamline (generic elements,
    dicts, DPW conversion, FELsim slice, transfer matrix), hybrid FELsim+COSY,
    COSY→RF-Track (chain, cross-validation), 3-way FELsim+COSY+RF-Track
    (chain, metadata, cross-validation)
  - CI: test_multicode.py + test_attribute_guard.py added to pipeline
- **Per-section config passthrough:** Runtime keys (space_charge, sc_mesh,
  physical_apertures, aperture) applied via setter methods before each
  section's tracking. Creation-time keys (G_quad, particle_mass) used to
  key simulator instances — different creation configs get separate instances,
  same creation config shares one instance.
- **Production validated:** FELsim→RF-Track hybrid at Stage 11 boundary
  (element 87) runs successfully. Hybrid vs full RF-Track shows qualitatively
  similar results (transverse RMS within order of magnitude) with expected
  differences from dipole model (transfer matrix vs analytical sector-bend).
- **COSY adapter integration (2026-03-10):**
  - `COSYAdapter.set_beamline()`: converts generic `BeamlineElement` objects
    to COSY beamline dict format, sets `_native_sim.beamline` directly.
    Handles DRIFT, QPF, QPD, DPH, DPW with full parameter passthrough.
  - `BeamlineBuilder.__init__` updated to accept `excel_path=None` (skips
    file validation when beamline will be set programmatically).
  - Auto-enables particle tracking at all elements when in particle_tracking
    mode, so `final_particles` is available for section handoff.
  - DPW triplet consolidation handled correctly: checkpoint at all elements
    avoids index mismatch from `_detect_dipole_triplets()`.
  - Verified: `MultiCode(felsim:0-10 + cosy:10-20)` runs to completion
    with particle handoff through COSY particle tracking.
  - COSY→RF-Track hybrid: `MultiCode(cosy:0-60 + rftrack:60-137)` runs to
    completion. Cross-validation vs full RF-Track: transverse RMS within
    order of magnitude — expected differences from COSY (DA maps + fringe
    fields) vs RF-Track (analytical sector-bend).
  - 3-way hybrid: `MultiCode(felsim:0-10 + cosy:10-60 + rftrack:60-137)`
    runs to completion. Metadata correctly records all 3 sections. Particle
    count non-increasing across handoffs. Cross-validation vs full RF-Track
    passes (transverse RMS ratios 0.1–10×).

### I5. T566 Objective via 2nd-Order DA Map [DONE — NOT NEEDED FOR UH FEL]
- **Status:** Fully implemented. `("l", "t566"): "ME(5,66)"` in MEASURE_MAP
  (`cosySimulator.py`). Validation enforces `transfer_matrix_order >= 2` and
  `dimensions >= 3`. FOX code generation maps to `ME(5,66)` in FIT blocks.
- **Use case:** Bunch compression optimization where both R56 and T566 matter.
- **W12 finding:** T566 = 0 for the UH FEL transport line (W9 Part A).
  A T566 FIT objective is redundant for this beamline but the implementation
  is ready for other beamlines with non-zero T566.

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

### I2. Engage with PALS as a Real-World Use Case [IN PROGRESS — brief drafted 2026-03-10]
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
- **Progress (2026-03-10):**
  - Use case brief drafted: `reports/pals_use_case_brief.md`
  - PALS contact: Jean-Luc Vay (@jlvay on GitHub)
  - Submission channel: GitHub discussion on pals-project/pals + email
  - Currently 3 participating codes (BMAD, BLAST ImpactX, BLAST WarpX);
    FELsim would be the first multi-code adapter-based participant
  - PALS confirms: no DIPOLE_WEDGE element, edge angles are BendP.e1/e2
    (matches our v3 implementation). Fringe via edge1_int/edge2_int.
  - **Remaining:** Open GitHub discussion, provide example v3 lattice file,
    join weekly meetings

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
- **Part C (C1B, 2026-03-10):** Four-way hybrid comparison via MultiCodeSimulator.
  Script: `C1B_hybrid_comparison.py` (`--smoke`, `--emittance`)
  - Methods: FELsim / MC-val / RFT-val / RFT-opt
  - MC-val: MultiCodeSimulator with FELsim(0:87) + RF-Track(87:137),
    FELsim-optimised currents — validates production hybrid architecture
  - **Full comparison (ε_n=5,8,14, 5 restarts):**
    | ε_n | FELsim | MC-val | RFT-val | RFT-opt |
    |-----|--------|--------|---------|---------|
    | 5 | 5.6e-3 | 1.16 | 1.84 | **6.2e-4** |
    | 8 | **2.5e-3** | 1.16 | 1.37 | 7.0e-2 |
    | 14 | 1.0e-1 | 12.7 | 2.84 | **9.9e-2** |
  - MC-val/RFT-val confirm large model differences (RMS ~1-13) when using
    FELsim-optimised currents in RF-Track.
  - RFT-opt recovers excellent matching at ε_n=5 (better than FELsim)
    but both methods struggle with β_y at ε_n=14.
- **Part D (C1C, 2026-03-10):** Five-way comparison adding MC-opt.
  Script: `C1C_multicode_optimization.py` (`--smoke`, `--emittance`)
  - Methods: FELsim / MC-val / MC-opt / RFT-val / RFT-opt
  - MC-opt: MultiCodeSimulator with FELsim(0:87) + RF-Track(87:137),
    NM optimization of Stage 11 quad currents [87, 93, 95, 97]
  - Architecture: mc_full and mc_disp share `_master_beamline` by reference;
    NM mutates `.current` on shared beamline, both simulators see updates
  - **Full comparison (ε_n=5,8,14, 5 restarts, SC=OFF):**
    | ε_n | FELsim | MC-val | MC-opt | RFT-val | RFT-opt | Winner |
    |-----|--------|--------|--------|---------|---------|--------|
    | 5 | 5.6e-3 | 1.16 | 6.2e-3 | 1.84 | **6.2e-4** | RFT-opt |
    | 8 | **2.5e-3** | 1.16 | 1.0e-1 | 1.37 | 7.0e-2 | FELsim |
    | 14 | 1.0e-1 | 12.7 | **9.0e-2** | 2.84 | 9.9e-2 | MC-opt |
  - No method dominates across all emittances. RFT-opt excels at ε_n=5
    (deep NM basin), FELsim wins at ε_n=8 (linear model adequate),
    MC-opt marginal winner at ε_n=14 (all struggle with β_y).
  - All optimizers struggle with β_y at ε_n=14 (0.01-0.08 vs target 0.24),
    indicating a physical limitation, not a code artifact.
  - **Space charge (ε_n=8 smoke test, SC=ON):** Identical MSE to SC=OFF.
    SC forces negligible at 40 MeV / 500 particles. Timing confirms SC solver
    active (MC-val: 1.1s SC=ON vs 0.3s SC=OFF).
  - SC wiring added: `--space-charge` now applies to MC-val/MC-opt via
    per-section config on RF-Track suffix (FELsim prefix unaffected).
  - Results: `results/C1C/`

### C3. FR3+MGE Optimization [IN PROGRESS — C3v4 on Koa]
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
  - **Cancelled (2026-03-10):** ~16k evals, best MSE=1000.005 (barely unstable).
    Sigma collapsed to 1.3e-05 — search exhausted, no stable solution found.
  - **C3v3 (Koa job 11471707):** maxfevals=8000/restart, BIPOP×9, tolx=1e-4.
    Killed after 247 evals. No improvement over warm start.
- **Root cause analysis (2026-03-11):** The objective function allows unstable
  solutions (penalty ~1000) to beat barely-stable solutions with bad Twiss
  (MSE >> 1000). CMA-ES *prefers* the unstable side of the boundary because
  the Twiss MSE for random stable solutions can be worse than the instability
  penalty. ~26k total evals across v1-v3, zero stable solutions.
- **C3v4 (Koa job 11473642, submitted 2026-03-11):** Three approaches:
  - **A (most likely):** Capped objective (stable MSE capped at 999, ensuring
    ANY stable solution beats ANY unstable one) + cold start [0]*23 + σ=2.0 +
    BIPOP×15 (up to 80k evals). Addresses the root cause directly.
  - **B:** Capped objective + warm start from v1 + σ=1.0 + BIPOP×15
  - **C:** Two-phase (stability-only first half → Twiss MSE second half)
  - Script: `C3v4_cosy_mge_opt.py`, `C3v4_cosy_mge_opt.slurm`
  - **Interim (2026-03-11, ~10h running):** Approach A completed 5006 evals —
    still unstable (cos_mu_x=1.2, cos_mu_y=-3.4, capped obj=2230). Now in
    BIPOP restart at ~1131 evals, best obj=6924 (unstable, improving).
    Zero stable solutions found in ~6000+ total evals. The FR3+MGE stability
    basin in 23D appears to be extremely narrow or nonexistent with the
    current fieldmap parameters.
- **C3v5: Fringe-shape homotopy + direct FR3 warm-start (2026-03-12):**
  After ~37k evals across v1-v4, zero stable solutions found. v5 introduces
  two fundamentally new ideas:
  1. **Fieldmap shape homotopy:** Continuously deform MGE field from uniform
     (alpha=0, hard-edge-like) to measured (alpha=1), optimizing at each step.
     `B_alpha(j) = B_mean + alpha*(B_measured(j) - B_mean)` — field integral
     preserved at all alpha. CMA-ES warm-starts from previous step. Bisection
     down to delta_alpha=0.025 if stability breaks.
  2. **Direct FR3 warm-start (Phase D):** CMA-ES from proven FR3 currents
     (MSE=2.5e-9) with full MGE — the obvious approach never tried in v1-v4.
  - Phases: 0 (verify FR0/FR3 at alpha=0,1), 1 (optimize alpha=0),
    2 (sweep alpha 0→1, 10 steps), 3 (polish), D (direct, parallel)
  - Est. ~27k evals (~14h) without bisection, up to ~100k with bisection
  - Scripts: `C3v5_homotopy_opt.py`, `C3v5_homotopy_opt.slurm`,
    `C3v5_direct_opt.slurm`
- **Files:** `fields/chicane_dipole_fieldmap.dat`, `test/koa_cosy_mge_opt.py`,
  `test/C3v4_cosy_mge_opt.py`, `test/C3v5_homotopy_opt.py`,
  `test/results/koa_cosy_mge_result.json`

### C4. Systematic Testing, Validation & Verification [IN PROGRESS]
- **Motivation:** FELsim currently relies on ad-hoc cross-validation studies.
  A systematic V&V programme is needed for production confidence.
- **Actions:**
  1. **Unit tests:** Core physics routines (transfer matrices, Twiss
     computation, dispersion, coordinate transforms) need pytest coverage
     with known analytic results (e.g., thin-lens quad, drift, FODO). [DONE]
  2. **Regression tests:** Each optimization study should produce a frozen
     reference result. CI runs confirm that code changes don't alter results
     beyond numerical noise. [DONE]
  3. **Cross-code benchmarks:** Extend S9/C1/C2 pattern — for each major
     beamline section, compare FELsim, RF-Track, and COSY Twiss functions
     element-by-element. Automate as a benchmark suite. [DONE — Tier 1+2]
  4. **Edge case testing:** ε_n → 0, σ_E → 0, single particle, 10⁵ particles,
     zero-length elements, degenerate optics (β → ∞).
  5. **Adapter round-trip tests:** Load lattice in all three formats
     (Excel/JSON/YAML), verify identical beamline objects. [DONE]
  6. **CI pipeline:** Automated test runs on commit (at minimum: unit tests
     + adapter round-trip + one optimization smoke test). [DONE]
- **test_cross_code_benchmark.py (2026-03-10):**
  - Tier 1 (CI, no external deps): frozen Twiss regression at 10 key
    checkpoints (500 particles, seed=42, 40 MeV, rtol=1e-10); y-emittance
    conservation; x-emittance pre-dipole; beamline geometry; determinism.
  - Tier 2 (RF-Track): drift/quad beta agreement within 10% (excluding
    dipole neighborhoods); full-beamline RMS envelope within 2 OOM;
    y-emittance conservation in both codes.
  - Physics note: x-emittance shows apparent growth through dispersive
    chicane regions (x-δ coupling). This is expected — the raw (x, x')
    emittance is not conserved in the presence of dispersion.
- **Current CI suite:** 210+ tests, 16 test files (expanded 2026-03-11 from 8).
  Added: `test_chromatic_dipole`, `test_chromatic_quad`, `test_optimizer`,
  `test_edge_cases`, `test_rftrack`, `test_felsim_unified`, `test_cosy_unified`,
  `test_rftrack_unified`. Tests requiring COSY/RF-Track gracefully skip.
  Visual tests excluded via `-m "not visual"`. Added `cloudpickle` to CI deps.
  Registered `cosy` and `rftrack` pytest markers in `conftest.py`.
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

### P8. Order-by-Order DA Convergence Study [DONE 2026-03-11]
- Runs COSY at DA orders 1, 2, 3, 5 with fixed quad currents.
  Two current sets: FR0 (hard-edge) and FR3 (3rd-order fringe).
- **FR=0 (hard-edge):**
  - Linear map and Twiss rigorously identical across all DA orders
  - FELsim's first-order model is exact for hard-edge elements
  - Aberration coefficients well-converged by O2 (no change at O3, O5)
  - Key chromatic terms: T_116 = 50.6, T_126 = 624, T_336 = -77.7
  - RMS = 4.8e-4 at all orders (Excellent)
- **FR=3 (3rd-order fringe):**
  - Linear map elements change by up to 0.3% with DA order
    (fringe-field nonlinearities feed into the linear map)
  - Twiss match best at O3 (RMS=5.9e-5) where currents were optimized;
    O1 gives RMS=6.4e-4 (still Excellent)
  - Key chromatic terms: T_116 = -205.5, T_126 = -1800,
    T_336 = 24.3, T_346 = -330.7
  - Large geometric aberration U_1111 = 21577
- **Cross-cutting findings:**
  - COSY particle tracker (via adapter) is order-independent — uses
    first-order element maps for tracking regardless of DA computation order
  - To assess nonlinear beam dynamics, must apply DA map to particles
    manually or use RF-Track
  - Map complexity: 15 coefficients (O1) → ~594 (O5)
- Script: `P8_order_convergence.py`
- Results: `results/P8/`

### P9. Chromaticity Analysis [DONE 2026-03-11]
- Swept energy deviation δ from -3% to +3% (13 points, 500 particles)
  with σ_E = 0.05% (near mono-energetic) to isolate chromaticity effects.
- **Achromatic transport:** Twiss constant across all δ (trivially, since
  matrices are momentum-independent).
- **Chromatic transport (achromatic currents):**
  - At δ=0: RMS = 0.087 (≈ achromatic, 0.075) — chromatic effects
    negligible for on-energy beam with narrow spread
  - At δ=±0.5%: RMS = 1.6–5.4 — dramatic degradation
  - At δ=±1%: RMS = 3.9–17 — catastrophic
  - Chromaticity: dβ_x/dδ ≈ -0.9 m/%, dβ_y/dδ ≈ -1.0 m/%
- **Acceptance bandwidth:** Only |δ| < ~0.3% achieves Acceptable matching.
  The transport line is highly chromatic; beam energy stability must be
  controlled to ~0.3% for adequate Twiss matching.
- Script: `P9_chromaticity_analysis.py`
- Results: `results/P9/`

### P12. Multi-Seed Robustness Study [DONE 2026-03-12]
- **Motivation:** S7 tested particle count sensitivity (fixed seed=42), S8 tested
  Stage 11 optimizer landscape (fixed seed, varied starting point). Neither varied
  the random beam realization itself. P12 fills this gap: does the 11-stage
  optimizer produce consistent results across different random seeds?
- **Design:** 20 seeds × 3 emittance points (5, 8, 14), 500 particles, cold-start.
  Seeds: [42 + 100·i for i in range(20)]. 60 total runs.
- **Results:**

  | ε_n | Exc | Acc | Mar | Fail | Best RMS | Median RMS | CV% |
  |-----|-----|-----|-----|------|----------|------------|-----|
  | 5 | 2 | 4 | 4 | 10 | 5.58e-3 | 3.16e-1 | 407% |
  | 8 | 13 | 3 | 1 | 3 | 4.57e-3 | 5.34e-3 | 203% |
  | 14 | 12 | 4 | 1 | 3 | 4.29e-3 | 4.81e-3 | 209% |

- **Key finding: Single-seed results are not representative.** At ε_n=5,
  50% of seeds fail entirely (RMS > 0.32). At ε_n=8 and 14, the majority
  (65–60%) achieve Excellent but ~15% still fail. The variability originates
  from upstream stages (beam generation), confirming S7/S8's findings that
  Stage 11 is unimodal but stages 1–10 are seed-sensitive.
- **Quad current variability:** Mean CV = 56% (ε_n=5), 25% (ε_n=8), 27% (ε_n=14).
  Worst: q16 (408% CV at ε_n=5), q80 (~130–193% CV at ε_n=8,14).
- **Implication:** Publication-quality results require multi-seed statistics
  (report median ± IQR, not single-seed values). Improving robustness
  requires upstream fixes (deterministic beam generation, multi-seed averaging,
  or ensemble optimization).
- **Plots:** RMS box plot, quad current CV% bars, quality histogram, Twiss scatter.
- **CLI:** `--seeds N`, `--particles N`, `--emittances ...`, `--plots-only`.
- Script: `P12_multi_seed_robustness.py`
- Results: `results/P12/`

### P10. Emittance Preservation Along Transport Line [DONE 2026-03-11]
- Tracked ε_n(s) element-by-element through 118-element beamline.
- **Achromatic transport:**
  - Dispersion-corrected ε_n conserved to -0.00% (both planes)
  - Raw ε_n,x grows 16% (x-δ dispersion coupling in chicane)
  - ε_n,y perfectly conserved (no vertical dispersion)
  - Peak η_x = 98 mm in chicane, peak raw ε_n,x = 55.5 π·mm·mrad
- **Chromatic transport:**
  - Dispersion-corrected ε_n,x grows **591%** — NOT conserved
  - Raw ε_n,x grows 659%, ε_n,y grows 387%
  - Cause: energy-dependent optics creates phase space filamentation
    (chromatic emittance dilution), appearing as irreversible growth
    in the sigma-matrix emittance
  - This is a known accelerator physics effect, not a code bug
- Script: `P10_emittance_evolution.py`
- Results: `results/P10/achromatic/`, `results/P10/chromatic/`

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
- **MODERATE (FIXED 2026-03-11):** Mutable default `shape={}` in `ebeam.plotXYZ` and
  `schematic.plotBeamPositionTransform` changed to `shape=None`. Callers already
  checked `if shape is not None` — no behavior change.
- **MODERATE (needs domain knowledge):** `dipole_wedge` fringe field integral
  (beamline.py:707) uses `le = self.length` which comes from `gap_wedge` in the Excel
  lattice. Needs clarification whether `gap_wedge` is the fringe field extent or the
  inter-dipole drift gap — if the latter, the K integral is wrong.
- **MODERATE (FIXED 2026-03-11):** `ebeam.ellipse_sym` (ebeam.py:79–82) — zero guards
  already present (`abs(beta) > 1e-30`, `denom > 0` checks). No change needed.
- **MODERATE (FIXED 2026-03-11):** `AlgebraicOptimization.py` (line 276) — `set.pop()` replaced
  with `sorted(sett, key=lambda s: s.name)` for deterministic variable ordering.

## Longer-Term Improvements (from multi-AI review 2026-03-10)

Source: 4-perspective expert review (FEL scientist, Berz-style computational physicist, SWE, UX/UI)

### Physics & Validation
- [x] **Order-by-order convergence study**: See P8 — first-order sufficient for hard-edge, ~0.3% linear map variation with fringe fields
- [x] **Emittance preservation plot**: See P10 — achromatic: conserved to 0.00%; chromatic: 591% growth (filamentation)
- [x] **Chromaticity analysis**: See P9 — acceptance bandwidth |δ| < ~0.3%, dβ/dδ ≈ 1 m/%
- [x] **Fringe field treatment in FELsim**: See P11 — DPW φ correction modifies M43 by 2–8%, removing it degrades RMS 7×. fringeType parameter is field-profile-only (drift matrix). Quad fringe not modeled.
- [ ] **Sensitivity / error analysis**: Magnet errors, misalignments, power supply ripple — DA methods can compute high-order sensitivities directly
- [x] **Multi-seed robustness study**: See P12 — 50% failure at ε_n=5, 65% Excellent at ε_n=8, 60% at ε_n=14. Single-seed not representative; upstream stages are seed-sensitive

### Code Quality
- [x] **Decouple optimization from visualization**: P9, P10, P11 now support `--plots-only` to regenerate plots from cached summary.json without re-running computation. S5, R2, S4 parameter scans already had this. P10 updated to save full evolution data.
- [ ] **Add pytest test suite**: Unit tests for Twiss computation, integration test for simplified optimization, visual regression with pytest-mpl
- [x] **Extract rcParams to .mplstyle file**: `felsim.mplstyle` created with shared settings (font, grid, DPI). P9-P11 and generate_seminar_figures.py updated to use it.
- [ ] **Make output formats configurable**: argparse for PDF/PNG selection

## Notes

- **Frontend ownership:** The frontend (`fel-app/`) is developed exclusively by Christian Komo. QA and code changes should focus on the backend. Minor frontend improvements are acceptable, but avoid making extensive changes that step on his work.
- All optimization scripts use `seed=42` for reproducibility
- NewFELsim conda environment required: `/home/evaletov/.conda/envs/NewFELsim/bin/python`
- Run commands from `backend/` with `MPLBACKEND=Agg PYTHONPATH=$(pwd)`
- **Review methodology**: Apply Berz-style computational physics perspective regularly when reviewing simulation results. Use multi-AI collab for second opinions.
