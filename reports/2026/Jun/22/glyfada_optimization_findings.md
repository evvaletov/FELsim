# glyfada optimization on the UH MkV FEL transport line: findings

Eremey Valetov, 2026-06-19
Re: Niels's action item "Try optimization using glyfada"

## Bottom line

glyfada has already been tried thoroughly and characterised. On FELsim's
transfer-matrix Twiss-matching objective its distributed evolutionary
algorithms underperform the production optimizer by 3 to 10 orders of
magnitude, and they cannot navigate the narrow feasibility basin in the
high-dimensional (26-current) variant at all. The production approach,
adopted as a result of these studies, is two-phase: Nelder-Mead (NM) to find
the basin, then in-process pycma CMA-ES to polish to near machine precision.
The integration is mature and remains wired in (`method='glyfada'`), so
re-running it is a single call. The one open use worth the effort is to run
glyfada or another global optimiser on the objective-design Config C stress
case, where sequential NM fails 75-80% of the time, as a head-to-head for the
IPAC paper. That is the only setting where glyfada might earn its keep, and
even there it is a paper deliverable, not a production need.

## Evidence

Four completed studies, all reproducible from `backend/test/`:

| Study | Setup | Result |
| --- | --- | --- |
| W6 (2026-02-23) | Stage-11 4-var match, eps_n = 5, 8, 14; glyfada pop=30 x 20 gen = 600 evals, uniform-random init over wide bounds, vs NM | glyfada loses at all three points. eps_n=5: all 600 evals hit the 1e6 penalty (no feasible point sampled). eps_n=8: glyfada MSE=2.96 vs NM 1.2e-6. eps_n=14: glyfada 0.31 vs NM 1.1e-3. NM beats by 3-6 orders. |
| W7 (2026-03-11) | Re-benchmark after fixing 3 glyfada config bugs (unset `initial_mean`, ignored `feasibility_rules` for CMA-ES, binary 0/1 constraint); NM vs warm-started pycma CMA-ES, eps_n = 5, 8, 14 | Two-phase NM->pycma wins decisively. eps_n=5: pycma 4.5e-16 vs NM 8.2e-6 (~1e10x, 10.5s vs 31.7s). eps_n=8: pycma 1.2e-15 vs NM 6.1e-6 (~1e9x). eps_n=14: pycma 6.3e-3 vs NM 8.9e-3 (~1.4x). Conclusion: the CMA-ES algorithm is what helps; glyfada's distributed infrastructure is not needed for a 4-variable problem. |
| O4 (2026-03-02) | glyfada on Koa, full 26-current single-shot objective; wide bounds [0,10 A] and tight bounds [NM +/- 2 A]; ~27,000 evals | ~99% of evaluations returned the penalty (unstable optics). Best solution found = the NM starting point itself. glyfada's auto-selector chose Simulated Annealing over CMA-ES/SHADE/NSGA-II; it characterised the landscape as rugged=1.0 with 7 modes. glyfada cannot improve on NM here. |
| O5 (2026-03-11) | Adopt the W7 winner into production | pycma CMA-ES warm-started from NM added as an automatic post-NM polishing step in `run_optimization()` (sigma=0.1, popsize=20, maxfevals=3000), with graceful `ImportError` fallback. Re-validated across the S4/S5/S7 scans: zero regressions, gains up to ~32x at difficult emittances. |

W7 numbers are taken from `backend/test/results/params_05ps/W7/rebenchmark_results.json`;
W6 numbers from `backend/test/results/params_05ps/W6/benchmark_results.csv`;
W6/O4 wording from `backend/test/W6_glyfada_benchmark_report.pdf` and
`backend/test/PRIORITIES.md` sections W6/W7/O4/O5.

## Why glyfada underperforms on this objective

The objective is a smooth, low-dimensional, transfer-matrix Twiss MSE with a
narrow feasible basin surrounded by a hard penalty wall. Three properties make
distributed evolutionary search a poor fit:

1. Narrow feasible basin. Outside a small region of current space the linear
   optics go unstable, Twiss invariants become NaN/inf, and the evaluator
   returns a 1e6 (DH protocol) or `-1e6` (fitness) penalty. In 26D the basin
   is narrow enough that ~99% of evaluations land on the penalty wall (O4).

2. Uniform-random initialization. glyfada seeds its population uniformly over
   the wide current bounds, so most initial members and most subsequent
   samples are infeasible. The 600-evaluation budget in W6 is spent almost
   entirely outside the basin; at eps_n=5 not a single feasible point was
   sampled.

3. Penalty wall destroys gradient. Collapsing every infeasible point to a
   single constant penalty erases the slope that would guide search toward the
   basin. The fitness landscape becomes a flat plateau with a thin spike, the
   worst case for population methods that rely on selection pressure across a
   smooth surface.

By contrast NM, started from physically motivated default currents, begins
inside or adjacent to the basin and descends the smooth interior. pycma
CMA-ES, warm-started from the NM point, adapts a covariance to the local
geometry and reaches machine precision. The distributed island/MPI machinery
glyfada provides buys nothing for a 4-variable smooth problem and is
counterproductive when the global landscape is mostly penalty.

## The one genuinely open use

The objective-design ablation (memory `objective_design_ablation.md`,
git-bug `7f690aa`) produced Config C: a variant of the production objective
with three plausible "fixes" (Stage-1 x.beta weight 0->0.5, Stage-7 envelope
target 0->1.5 mm, dispersion rescaling). Config C makes sequential per-stage
NM fail 75% of the time (20 seeds) in the original run and 80% at the MOP6318
targets, because Stage 10 shifts the basin and Stage 11 NM gets stuck. This
is exactly the regime where a global optimiser could matter: the question of
whether BO/CMA-ES/glyfada recovers Config C's NM failures is the open one the
IPAC paper is meant to answer. It overlaps the planned S4 (CMA-ES drop-in)
and S6 (BO baseline) deliverables under `7f690aa`.

Running glyfada on Config C is therefore a legitimate, scoped experiment, but
it is a negative-or-positive-result paper deliverable, not a production need.
Expect it to either (a) confirm that even global search struggles on Config
C's pathological basin, strengthening the "objective design matters" message,
or (b) recover where NM fails, which is the interesting positive result. Pair
it head-to-head with NM and the BO baseline on the same A/B/C matrix and fixed
eval budget so the comparison is apples-to-apples.

## Recommendation

1. Treat the glyfada action item as largely done. The honest answer to "try
   glyfada" is "tried, characterised, and superseded by NM->pycma-CMA-ES,
   which is now the production optimizer."

2. Do not adopt distributed glyfada for the standard transport-line match.
   It is strictly dominated by the two-phase in-process optimizer on this
   objective, at lower wall-time and far higher precision.

3. If a glyfada run is wanted for the paper, target Config C only, as the
   global-search stress test, and report it alongside NM and the S6 BO
   baseline. Run it in the `NewFELsim` conda env; the `glyfada` module is not
   importable from system `python3`.

## Reproduction notes

- Integration: `backend/glyfadaAdapter.py` (delegates to `glyfada.optimize()`)
  and `backend/glyfada_eval.py`, dispatched from `backend/beamOptimizer.py`
  via `method == 'glyfada'`. Binary at `~/ML/paradiseo/glyfada/build/optimiser`.
- Environment: `glyfada` imports only in the `NewFELsim` conda env
  (`/home/evaletov/.conda/envs/NewFELsim/bin/python`), not in system python3.
- Re-run W6/W7: `backend/test/W7_glyfada_rebenchmark.py`; W6 via the `--w6`
  flag in `UHM_beamline_opt_05ps_params.py`.
- Config C ablation: `backend/test/ablation_run.py`, `ablation_analyze.py`,
  `diagnose_seed91.py`; results in `backend/test/results/ablation/` and
  `results/ablation_MOP6318/`.
