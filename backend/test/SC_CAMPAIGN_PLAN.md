# Space-Charge Study Campaign — Execution Plan

**Source:** 2026-05-04 meeting with Niels Bidault (git-bug `f59e144`).
**Authored:** 2026-05-29 (Eremey Valetov). Planning input: 3-way with codex (gpt-5.5) + agy (Gemini 3.1 Pro).
**Canonical checklist:** `PRIORITIES.md` → "Action items from 2026-05-04 meeting with Niels".

## Item index

| # | Item | Phase | Status |
|---|------|-------|--------|
| a | Linac: Xsuite — more detail (multi-cell TW, currently lumped Cavity) | 2 | Blocked on cavity geometry from Niels |
| b | Linac: TW (COSY TW extension N28; RF-Track TW done) | 2/4 | COSY TW = hard DA dev |
| c | COSY SC vs equivalent macroparticle charge (Q/N_p) | 1 | **In progress (plots)** |
| d | COSY SC charge-density profile ρ(r) + potential Φ(r) | 1 | Cheap (expose φ from treecode) |
| e | COSY SC at 1 MeV | 1 | Needs re-match first (see Risk R1) |
| f | Longitudinal SC effects | 4 | Hard — Fortran treecode core edit |
| g | 3-code: 16 ps + 1 ps bunch | 1→3 | FODO now; line later |
| h | 3-code: small spot size | 1→3 | FODO now; line later |
| i | 3-code: RF-Track / Xsuite with SC | 1→3 | Works on FODO; line = item l |
| j | BEAMPATH — injector cross-check | 4 | Blocked on LANL access + deck |
| k | Later: full PIC in COSY (git-bug `1cbd8ea`) | 4 | Long-term dev |
| l | DA-FMM+RFT+Xsuite on focusing/transport line, 45 & 1 MeV, ± dipoles | 3 | **Integrative capstone** |

## Phases

**Phase 1 — Close FODO space-charge physics (quick wins, unattended now).**
Isolate SC physics on the existing 2 m FODO before adding lattice complexity. Covers c, d, e, and the FODO-level slices of g/h/i. Most of the raw data already exists from the L3 studies; the work is re-analysis + a few cheap new runs.

**Phase 2 — Make the linac models comparable.**
SC comparison through accelerating structures is meaningless until RF optics align. Upgrade the xsuite linac from a lumped Cavity to a multi-cell TW model (a); implement the COSY TW extension (b); re-confirm energy gain vs elegant and transverse focusing vs RF-Track TW.

**Phase 3 — Three-code transport-line benchmark (the capstone, item l).**
DA-FMM (COSY) + RF-Track + Xsuite on the FELsim focusing/transport line at 45 MeV and 1 MeV, with and without dipoles. Build order: no-SC handoff verification → linear/frozen SC → PIC/FMM. Requires an **Xsuite adapter** in the multi-code framework (new) and resolving the RF-Track full-line PIC core-dump (`af9d56c`).

**Phase 4 — Injector + advanced SC.**
BEAMPATH injector cross-check (j, needs access), longitudinal SC in COSY (f, core treecode edit + RF-bucket normalization), full PIC in COSY (k).

## Quick wins (doable now, no Niels input) — Phase 1
- **(c)** Re-plot L3 sweeps as emittance growth vs equivalent macroparticle charge q_mp = Q/N_p; extract the threshold q_mp\* where numerical collisionality overtakes physical mean-field SC. *(Pure analysis — done this session.)*
- **(d)** Expose Φ at macroparticle positions (one line in `spch_kick.f90`; treecode already computes it), add a gridded Φ(x,y)/ρ(x,y) diagnostic.
- **(e)** 1 MeV FODO run — but re-match the lattice first (Risk R1).
- **(g/h)** 1 ps / 16 ps and small-spot FODO variants (config + 1-line FOX edits).
- Three-code FODO SC comparison plots from existing `sc_compare_output/` data.

## Blocked / needs input
- **(a) xsuite multi-cell TW** — needs cavity geometry / field-map + phase convention from Niels.
- **(b) COSY TW** — symplectic DA map for a time-dependent TW is non-trivial dev.
- **(f) longitudinal SC** — edit the Fortran treecode leapfrog (z-kick) + define RF-bucket normalization; not a FOX-only change.
- **(j) BEAMPATH** — LANL license/access + reference input deck; agree on authoritative output.
- **RF-Track full-line PIC** — prior core-dump; needs a bisected minimal-lattice reproducer.

## Risks & mitigations
- **R1 — 1 MeV ≠ same FODO.** Matched Twiss changes drastically at low β with SC; naively reusing the 45 MeV FODO gives meaningless emittance metrics. *Mitigate:* re-match (adjust 1/f) before every 1 MeV run; add an automated Twiss re-match step to the framework.
- **R2 — false code disagreement from coordinate/convention mismatch.** *Mitigate:* no-SC single-particle + beam-moment benchmark at every code handoff and section boundary, before turning SC on.
- **R3 — SC shot noise mistaken for physical growth.** *Mitigate:* always pair runs with N_p convergence + a fixed-q_mp family and the linear-Gaussian baseline (demob_gaussian).
- **R4 — dipoles dominate discrepancies.** *Mitigate:* run the line without dipoles first, then one dipole, then the full set.
- **R5 — RF-Track PIC stalls the capstone.** *Mitigate:* keep DA-FMM/COSY + xsuite-frozen as the primary comparison path; isolate the PIC crash separately.

## Cross-cutting infrastructure (raised by codex + agy — adopt before Phase 3)
- **Common initial-distribution generator** — identical macroparticles fed to COSY / RF-Track / Xsuite wherever possible.
- **Explicit agreement tolerances** — % envelope error, emittance-growth tolerance, centroid, energy-spread bands defined up front.
- **Reproducibility manifest** per run — code versions, seed, N_p, N_slice, mesh, Δs/Δt, SC solver settings.
- **Wall-clock tracking** — accuracy-per-cost is a first-class result (DA-FMM is ~100× slower than xsuite-frozen on the FODO).

## On item (c) specifically
The L3 finding (DA-FMM emittance growth ∝ macroparticle shot noise ~1/√N_p) is the *foundation*, not the full answer. Niels likely wants the **rule-of-thumb threshold**: the q_mp = Q/N_p at which numerical collisionality equals the physical mean-field growth. Deliverable: excess growth (DA-FMM − xsuite) vs q_mp, collapsing the charge-sweep and N_p-sweep onto one trend, with the crossover q_mp\* annotated.

---

## Execution log
- **2026-05-29** — Plan authored + 3-way reviewed. Phase 1 result plots generated from existing L3 data (`results/sc_campaign/`): see `sc_campaign_plots.py`. Items addressed: (c) headline + N_p convergence, three-code FODO comparison (i/l preview), N_slice robustness.
- **2026-05-29** — Phase 1 extended with new COSY DA-FMM runs (`cosy-fmm`: `spch_kick.f90` now dumps `phi.dat`; `spch_demo.py` parameterised by energy/spot/bunch via `sc_field_profile.py`, `sc_energy_scaling.py`):
  - **F5** (d): ρ(r)/Φ(r)/E_r(r) profiles, 1 nC/45 MeV; E_r matches the analytic 2D-Gaussian field (peak Φ≈609 V, E_r≈5.9e4 V/m).
  - **F6** (e): energy scaling — 1 MeV @0.001 nC = 123% vs 45 MeV @1 nC = 1.8%; ~6 decades over 1→45 MeV (1/β²γ³). *Fixed-beam; SC-matched comparison is Phase 3.*
  - **F7** (g/h): bunch length 1→16 ps (3.99→0.58%) and spot 0.25→4 mm (26→0.01%).
  - Full montage: `results/sc_campaign/overview.png`.
- **2026-05-29** — Phase 2/3 enabler + F6 matched fix:
  - **XsuiteAdapter** added (`backend/xsuiteAdapter.py`) and wired into the
    multi-code framework (CoordinateSystem.XSUITE + FELSIM↔XSUITE transforms +
    factory registration). Verified: roundtrip 5.5e-13, drift vs FELsim 7.8e-10,
    stable FODO vs FELsim 7.5e-5, SC frozen smoke test; test_multicode 38/38.
    Unblocks the three-code transport-line capstone (item l).
  - **F6 re-done as a matched comparison** (Risk R1 resolved): the demo FODO
    (1/f=2) is at μ=180° with no matched optics, so the study moved to a stable
    FODO (1/f=1) with a beam matched to its bare optics over 4 cells. 45 MeV @1 nC
    = 0.31% vs 1 MeV @1 nC = 4.6e4% (SC-limit blowup). Tracker gained
    `bare_matched_twiss` / `make_matched_bunch` and k_focus/n_cells params.
  - Results mirrored to `/mnt/hgfs/Documents/003_Niels/sc_campaign/` (PNG + ASCII).
