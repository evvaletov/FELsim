# cosy-fmm space charge in FELsim — capstone groundwork + first cross-code SC

Eremey Valetov, 2026-06-19. Tracker: git-bug `f59e144` item l (2026-05-04 Niels list);
PRIORITIES.md item 9 ("COSY DA-FMM SC + PIC"). Code: `backend/test/sc_capstone/`.

## Context

The 2026-05-04 ask was a three-code space-charge benchmark — cosy-fmm DA-FMM (the
COSY-side Coulomb treecode), RF-Track, and xsuite — on the FELsim transport line at
45 MeV and 1 MeV, with and without dipoles. A 2-model design panel (codex gpt-5.5 +
agy Gemini 3.1 Pro) flagged the sequencing trap up front: **you cannot benchmark
full-line space charge while the codes still disagree on linear transport**, because
xsuite has no dipole edge/fringe model and the DA-FMM kicks only happen in drifts. So
this work builds the prerequisite layer first, then runs the SC comparison where it is
currently meaningful.

Three components:

1. **Common distribution + reproducibility manifest** (`common_distribution.py`) — one
   macroparticle array in FELsim coordinates from a manifest that records every
   parameter, the code versions, the FELsim + cosy-fmm git hashes, and the array
   SHA-256. Every engine starts from the identical particles.
2. **No-SC cross-code handoff** (`nosc_handoff.py`) — the cold-beam agreement table.
3. **Reduced SC capstone** (`sc_section.py`, `sc_capstone_run.py`) — cosy-fmm DA-FMM vs
   xsuite frozen-Gaussian on the real no-dipole section, same beam, SC on.

All numbers below: N_p = 6000, ε_n = 8 mm·mrad, β_x = β_y = 3 m, seed 20260619.

## Result 1 — cold-beam agreement (no space charge)

**No-dipole section, elements [32, 46)** (6 quads, 1.65 m — the longest drift+quad run
in the line). FELsim / xsuite / RF-Track agree to **< 0.1 %**:

| quantity | FELsim | xsuite | RF-Track | worst vs FELsim |
|---|---|---|---|---|
| σ_x [mm] | 0.27581 | 0.27593 | 0.27593 | 0.041 % |
| σ_y [mm] | 0.24857 | 0.24853 | 0.24853 | 0.018 % |
| ε_{n,x} [mm·mrad] | 8.1645 | 8.1702 | 8.1702 | 0.070 % |
| ε_{n,y} [mm·mrad] | 7.9511 | 7.9505 | 7.9505 | 0.007 % |

→ on this section the codes are calibrated and space charge can be compared.

**Full line, elements [0, 137)** — not yet comparable, exactly as the panel predicted.
xsuite **fails outright** (no dipole edge/fringe model), and FELsim vs RF-Track disagree
by **~18×** on σ_x (1279 mm vs 27 mm; ε_{n,x} 321 vs 72 mm·mrad). The chicane/spectrometer
dipole models, the DPW edge-kick workarounds, and dispersion × σ_δ diverge between codes.
Turning SC on over the full line now would report dipole-model differences as
space-charge disagreement. **The +dipole headline needs the xsuite dipole edge/fringe
model (and SC-inside-magnets element slicing) first.**

## Result 2 — space charge on the no-dipole section (DA-FMM vs xsuite-frozen)

Same common distribution, σ_δ = 0 (isolates SC from the ~2 % chromatic emittance growth
both chromatic codes otherwise show). Split-operator slicing, ds ≈ 20 mm. The DA-FMM
tracker's SC-off transport reproduces FELsim **exactly** (0.000 %, emittance conserved to
1e-13), so the SC-on growth is trustworthy. For the 1 MeV rows the quad currents are
scaled to preserve the 45-MeV optics (analytic re-match), so the energy comparison
isolates space-charge scaling rather than a focusing instability (the 45-MeV-tuned quads
are over-focused ~32× at 1 MeV and would blow up the beam on their own — Risk R1).

ε_{n,x} growth (figure: `capstone_combined_growth.png`):

| E [MeV] | Q [nC] | DA-FMM | xsuite-frozen | DA/xs |
|---|---|---|---|---|
| 45 | 0.1 | −0.017 % | −0.003 % | (noise) |
| 45 | 0.3 | −0.005 % | −0.006 % | (noise) |
| 45 | 1 | +0.255 % | +0.012 % | 21.7 |
| 45 | 3 | +2.810 % | +0.292 % | 9.6 |
| 45 | 10 | +26.95 % | +4.315 % | 6.2 |
| 1 | 0.01 | +0.843 % | +3.992 % | 0.21 |
| 1 | 0.03 | +8.901 % | +42.54 % | 0.21 |
| 1 | 0.1 | +49.87 % | +245.1 % | 0.20 |
| 1 | 0.3 | +150.7 % | +790.3 % | 0.19 |

SC-off control (both engines, every point): |growth| < 1e-12 %.

Two regimes, two distinct physics findings:

- **45 MeV (weak SC).** The DA-FMM excess over xsuite (6–22×) is *not* a physical
  nonlinear effect — at N_p = 6000 the macroparticle charge q_mp = Q/N_p exceeds the
  Phase-1 threshold q_mp* ≈ 0.037 pC for Q ≥ 1 nC, so the DA-FMM excess is **numerical
  collisionality (shot noise ∝ q_mp^0.46)**, which the analytic frozen model does not
  have. An N_p sweep at 45 MeV/1 nC (`np_convergence.png`) confirms the physical growth
  here is only ≲ 0.05 % — below the DA-FMM single-realization shot-noise floor, so 45 MeV
  is too weak-SC to resolve a clean code-vs-code physical difference at single seed. This
  **re-demonstrates the Phase-1 q_mp* threshold on a real transport section** (previously
  only the toy FODO).

- **1 MeV (strong SC), fair comparison.** At 1 MeV/0.01 nC, q_mp = 0.0017 pC ≪ q_mp*, so
  there is no shot-noise inflation and the comparison is physical. Here the ordering
  **reverses**: the xsuite frozen-Gaussian over-predicts emittance growth by **~5×**
  relative to the DA-FMM N-body treecode (0.84 % vs 3.99 %; 151 % vs 790 % at 0.3 nC).
  Mechanism (note: `SpaceChargeBiGaussian` is *not* a linearized model — it applies the
  exact **nonlinear** Bassetti–Erskine field of a Gaussian; with σ_δ = 0 a purely linear
  kick would give zero rms-emittance growth by Liouville, so the growth in both codes is
  genuine nonlinear filamentation). The frozen model re-imposes a Gaussian
  (Bassetti–Erskine) field reconstructed from the *evolving rms* at every step, so it
  keeps injecting the Gaussian profile's nonlinear field — and spurious electrostatic free
  energy — even after the real beam core has flattened/hollowed, whose true field is more
  linear and whose free energy is depleted (the physical growth saturates). The N-body
  treecode tracks the actual non-Gaussian distribution and saturates; the frozen model
  over-drives filamentation. **This is the headline DA-FMM-vs-xsuite physics result: where
  space charge is strong enough to matter, the rms-self-consistent frozen model
  over-estimates growth, and the cosy-fmm treecode is the corrective.** (Direction +
  mechanism cross-checked against an independent 2-model review, 2026-06-19.)

## Method validation

- Common distribution: deterministic (manifest → identical array + SHA-256), emittance
  and Twiss recovered to < 2 %, FELsim↔SI round-trip exact. 7/7 unit tests.
- DA-FMM section tracker: SC-off reproduces FELsim transport to 0.000 % with emittance
  conserved to 1e-13; xsuite-frozen SC-off control 1e-12 % (k1 convention verified). 4/4
  smoke tests. Total suite 11/11.

## How to reproduce

```
cd backend/   # env: /home/evaletov/.conda/envs/NewFELsim/bin/python
export MPLBACKEND=Agg PYTHONPATH=$(pwd):$(pwd)/test/sc_capstone
python test/sc_capstone/common_distribution.py --out test/results/sc_capstone/dist_default
python test/sc_capstone/nosc_handoff.py --codes felsim xsuite rftrack --cases nodip full
python test/sc_capstone/sc_capstone_run.py --charges-nc 0.1 0.3 1 3 10 --energies-mev 45 --optics-energy 45 --outdir test/results/sc_capstone/e45
python test/sc_capstone/sc_capstone_run.py --charges-nc 0.01 0.03 0.1 0.3 --energies-mev 1 --optics-energy 45 --outdir test/results/sc_capstone/e1
python test/sc_capstone/combine_capstone.py
python test/sc_capstone/np_convergence.py
pytest test/sc_capstone/ -q
```

## What is next (panel-ranked)

1. **xsuite dipole edge/fringe model + SC-inside-magnets element slicing** → unlocks the
   full-line and +dipole comparison (currently the binding gap). Tracked: git-bug `660da00`
   (P1-high).
2. **RF-Track full-line PIC core-dump (v2.5.5)** isolated as a parallel debug task, not a
   capstone blocker.
3. **DA-level kick composition** (`fmm_eval_treecode_da`) — COSY's unique high-order-DA
   differentiator, the route Niels asked about; deferred until the baseline capstone
   defines the reference.
4. **Clean weak-SC number at 45 MeV** — now supported in `sc_capstone_run.py` via two
   routes (added 2026-06-20): `--seeds s1 s2 …` averages DA-FMM over realizations
   (reports mean ± std; at 45 MeV/1 nC, 3 seeds give 0.19 % ± 0.13 % — the physical value
   is buried in sampling variance), and `--softening <eps_m | auto>` applies Plummer
   softening to suppress the close-encounter collisional heating (validated: at 10 nC,
   ε = 3e-4 m pulls DA-FMM 26.9 % → 18.5 % toward xsuite 4.3 %; ε must be a meaningful
   fraction of σ_x — the `auto` = σ_x/√N_p ≈ 6.7 µm heuristic is too small to suppress
   much). Multi-seed handles low-charge sampling variance; softening handles high-charge
   collisional excess.
