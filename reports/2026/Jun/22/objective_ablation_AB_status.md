# Objective-design ablation: A/B robustness status

Eremey Valetov, 2026-06-22
For: Niels Bidault

## Bottom line

At the MOP6318 undulator-match targets (beta_x = 1.267 m, alpha_x = 0.560),
the original objective (Config A) and the mildly-rescaled objective (Config B)
each fail 35% of the time (7/20 seeds) by the undulator-RMS > 1e-2 criterion.
This is not robust, contrary to the seminar/abstract framing. The defensible
"fixed" objective (Config C) fails 80% (16/20). The narrative shift: sequential
Nelder-Mead is not a reliable matcher for this lattice at these targets, which
strengthens the case for the BO baseline as the paper's deliverable rather than
weakening it.

## Evidence

Two ablations, 20 random-seed beams per config, 11-stage sequential NM, identical
setup except the Stage-11 targets.

| Targets | Source | A fail | B fail | C fail |
| --- | --- | --- | --- | --- |
| (1.4 m, 0.47) | arXiv:2510.14061 radiation-mode | 2/20 (10%) | 3/20 (15%) | 15/20 (75%) |
| (1.267 m, 0.560) | MOP6318 Table 1 | 7/20 (35%) | 7/20 (35%) | 16/20 (80%) |

Counts re-derived from the JSONs at threshold RMS > 1e-2 (the pass/fail rule in
`ablation_analyze.py`). MOP6318 results: `backend/test/results/ablation_MOP6318/`.

Mechanism (diagnose_seed91.py, confirmed on the MOP6318 data): failures are a
vertical-plane collapse at the final undulator match. Across the seven failing
A seeds, beta_x stays on target (1.24-1.27 m) while beta_y collapses to
0.01-0.12 m against the 0.242 m target. Stage 10 (Triplet C) reaches a lower
local MSE on the failing seeds but shifts the alpha/beta basin so Stage 11 then
lands ~300x worse (seed 91 Stage-11 MSE 7.1e-3 vs passing seed 3 at 2.4e-5).
The distribution is bimodal: seeds either converge to RMS ~ 1e-4 or land near
1e-1, nothing in between. Iteration counts are flat across configs (medians
1654-1808), so this is basin selection, not budget.

## Done

- Ablation study run twice (old radiation-mode targets and MOP6318 targets),
  20 seeds x 3 configs each. Logic: `backend/test/ablation_run.py`. Analysis +
  threshold: `backend/test/ablation_analyze.py`. Interpretation written into
  `backend/test/results/ablation_MOP6318/analysis/summary.md`.
- Paper figures F2 (A/B/C fail-rate, 35/35/80% at MOP6318 targets) and F6
  (GD-vs-BO template) built. F4 shows the per-stage Stage-11 divergence; F5
  shows the beta_y collapse. Figures in `backend/test/results/ipac/`
  (git-bug cf11eb7).

## Open

- S6 BO baseline (git-bug 7f690aa) is the actual paper deliverable and is
  blocked on your xopt hyperparameters. F6 renders the NM bars as real and the
  BO bars as a placeholder until S6 runs on the same A/B/C matrix.
- Twiss provenance to reconcile: the figures use beta_x = 1.267 m,
  alpha_x = 0.560, which I took from a code comment. These match the MOP6318
  Table 1 values, but I want to confirm them against the abstract before the
  proceedings write-up.
- Remaining ablation follow-ups (git-bug 7f690aa): S4 CMA-ES drop-in, S5 joint
  26-variable single-shot optimization (tests whether the basin escape is an
  artifact of the sequential chain), S7 warm-chained stage starting points.

## Ask

1. Your xopt hyperparameters (acquisition, kernel, eval budget per stage) so
   S6 runs on the same setup as your linac optimization.
2. Confirm the narrative shift for the proceedings: at the MOP6318 targets,
   A and B both fail ~35%, so the framing moves from "the original objective is
   robust" to "sequential NM is unreliable here, motivating the BO baseline."
   This changes the abstract's framing and needs your sign-off.
