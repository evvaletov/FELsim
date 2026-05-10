# 11-stage NM objective-design ablation
Three objective configurations were run on the UH FEL transport line with 20 random-seed beams each, plus a single deterministic Sobol beam per config. NM is sequential per-stage; starting points are reset per stage to the configured defaults.
## Configurations
- **A**: Verbatim original (`branch_niels/UHM_beamline_opt.py`).
- **B**: A + mild per-measure-type rescaling. Each squared residual is divided by `MEASURE_REF**2`; only dispersion is actually rescaled (`ref = 0.5 m`, weight x4). Other measures have `ref = 1` so their numerical contribution is unchanged.
- **C**: B + Stage 1 `x.beta` weight `0.0 -> 0.5` (typo fix) + Stage 7 `envelope` goal `0.0 mm -> 1.5 mm` (finite physical target).

## Final undulator RMS (20 seeds per config)
| Config | Median | Min | Max | Mean | Std | Fail rate (RMS>1e-2) | Sobol |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | 4.436e-05 | 3.502e-05 | 1.186e-01 | 1.106e-02 | 3.393e-02 | 2/20 (10%) | 4.131e-05 |
| B | 1.804e-04 | 1.492e-04 | 1.125e-01 | 1.581e-02 | 3.825e-02 | 3/20 (15%) | 1.560e-04 |
| C | 1.122e-01 | 1.134e-04 | 1.814e-01 | 8.672e-02 | 5.437e-02 | 15/20 (75%) | 8.652e-02 |

## Total NM iterations (20 seeds per config)
| Config | Median | Min | Max | Mean |
| --- | --- | --- | --- | --- |
| A | 1816 | 1285 | 2381 | 1790 |
| B | 1808 | 1000 | 2381 | 1791 |
| C | 1647 | 1206 | 2455 | 1648 |

## Per-seed detail
| Config | Seed | RMS | Total iters | Runtime (s) |
| --- | --- | --- | --- | --- |
| A | 1 | 1.186e-01 | 1285 | 2.9 |
| A | 2 | 5.673e-05 | 1901 | 4.4 |
| A | 3 | 4.076e-05 | 1608 | 3.7 |
| A | 5 | 4.132e-05 | 1886 | 4.3 |
| A | 7 | 5.603e-05 | 1783 | 4.1 |
| A | 11 | 3.502e-05 | 2097 | 4.8 |
| A | 13 | 6.982e-05 | 2381 | 5.6 |
| A | 17 | 4.072e-05 | 1485 | 3.4 |
| A | 19 | 3.967e-05 | 1613 | 3.7 |
| A | 23 | 4.288e-05 | 1917 | 4.3 |
| A | 29 | 5.028e-05 | 1621 | 3.7 |
| A | 31 | 3.828e-05 | 2096 | 4.8 |
| A | 37 | 3.919e-05 | 1849 | 4.2 |
| A | 41 | 4.139e-05 | 1512 | 3.5 |
| A | 42 | 5.642e-05 | 1897 | 4.3 |
| A | 43 | 3.513e-05 | 2096 | 4.8 |
| A | 47 | 4.584e-05 | 1498 | 3.4 |
| A | 53 | 5.161e-05 | 1709 | 3.9 |
| A | 91 | 1.013e-01 | 1577 | 3.7 |
| A | 137 | 6.468e-04 | 1986 | 4.6 |
| B | 1 | 1.125e-01 | 1781 | 4.1 |
| B | 2 | 2.268e-04 | 1901 | 4.4 |
| B | 3 | 1.783e-04 | 1602 | 3.7 |
| B | 5 | 1.521e-04 | 1885 | 4.3 |
| B | 7 | 2.094e-04 | 1728 | 4.0 |
| B | 11 | 1.492e-04 | 2102 | 4.7 |
| B | 13 | 2.565e-04 | 2381 | 5.4 |
| B | 17 | 9.353e-02 | 1000 | 2.2 |
| B | 19 | 1.756e-04 | 1730 | 4.0 |
| B | 23 | 1.750e-04 | 1935 | 4.4 |
| B | 29 | 1.742e-04 | 1667 | 3.8 |
| B | 31 | 1.519e-04 | 1997 | 4.6 |
| B | 37 | 1.667e-04 | 1897 | 4.4 |
| B | 41 | 1.746e-04 | 1557 | 3.6 |
| B | 42 | 2.147e-04 | 1834 | 4.2 |
| B | 43 | 1.523e-04 | 2060 | 4.6 |
| B | 47 | 1.825e-04 | 1455 | 3.3 |
| B | 53 | 4.760e-04 | 1742 | 4.0 |
| B | 91 | 1.068e-01 | 1577 | 3.6 |
| B | 137 | 2.281e-04 | 1986 | 4.6 |
| C | 1 | 1.168e-01 | 1648 | 3.7 |
| C | 2 | 1.101e-01 | 1206 | 2.7 |
| C | 3 | 1.166e-01 | 1341 | 3.0 |
| C | 5 | 1.190e-01 | 1661 | 3.9 |
| C | 7 | 2.239e-04 | 2093 | 4.8 |
| C | 11 | 8.548e-02 | 2034 | 4.8 |
| C | 13 | 1.144e-01 | 1582 | 3.6 |
| C | 17 | 1.814e-01 | 1646 | 3.9 |
| C | 19 | 1.261e-04 | 1639 | 3.8 |
| C | 23 | 1.251e-04 | 1659 | 3.8 |
| C | 29 | 2.290e-04 | 2079 | 4.7 |
| C | 31 | 1.177e-01 | 1335 | 3.1 |
| C | 37 | 1.192e-01 | 2455 | 5.7 |
| C | 41 | 1.126e-01 | 1223 | 2.8 |
| C | 42 | 1.134e-04 | 1892 | 4.3 |
| C | 43 | 1.182e-01 | 1707 | 4.0 |
| C | 47 | 1.091e-01 | 1457 | 3.3 |
| C | 53 | 1.145e-01 | 1317 | 3.0 |
| C | 91 | 1.118e-01 | 1261 | 2.8 |
| C | 137 | 8.674e-02 | 1724 | 4.0 |

<!-- INTERPRETATION (preserve when regenerating) -->

## Interpretation

**Headline result, 20 seeds + 1 Sobol per config:**

| Config | Fail rate | When successful, RMS | Sobol RMS |
| --- | --- | --- | --- |
| A | 2/20 (10%) | 3.5e-5 to 7e-5 (median 4.4e-5) | 4.1e-5 (success) |
| B | 3/20 (15%) | 1.5e-4 to 5e-4 (median 1.8e-4) | 1.6e-4 (success) |
| C | **15/20 (75%)** | 1.1e-4 to 2.2e-4 when good | 8.7e-2 (**fail**) |

Config C's "fixes" (typo + envelope=0 replaced, on top of mild dispersion
rescaling) introduce a **pathological basin structure** for sequential
per-stage NM. 75% of seeds get stuck in a bad basin. The deterministic
Sobol beam also lands there, which is consistent with the multi-seed
finding rather than confirming it independently (it is a single
data point).

Configs A and B remain robust (10-15% fail rate). Both Sobol runs for
A and B land within the typical success-case RMS range, consistent
with the multi-seed result.

**Where the failure happens (per-stage diagnostic; see
`diagnose_seed91.py --target 91 --compare-seed 42`):**

Per-stage MSE comparison between a working seed (42) and a failing
seed (91) shows that:
- Stages 1-9 converge to similar MSE for both seeds.
- Stage 10 (Triplet C) shows a 28x worse final MSE for the failing
  seed (3.8e-5 vs 1.3e-6), still small in absolute terms.
- Stage 11 (Undulator match) is where the divergence becomes severe:
  290x worse final MSE for the failing seed (8.3e-3 vs 2.9e-5).

The mechanism: Stage 10 leaves the beam in a slightly different
state (different alpha/beta at segs 85/86), which shifts Stage 11's
starting basin. NM's local-only search at Stage 11 then gets stuck
in a local minimum. The "fixes" in Config C amplify this effect
across the upstream stages — more seeds land in shifted basins.

**Implications for the IPAC paper:**

1. **Single-seed comparisons mislead.** With 5 seeds we estimated
   Config C at 60% fail; with 20 seeds the actual rate is 75%.
   Multi-seed (>=20) is essential for the GD-vs-BO benchmark.

2. **Config A is the right NM baseline.** Lowest median RMS, lowest
   fail rate, and it matches what was used to produce the seminar
   figures.

3. **The "fixes" we identified in the email (typo, envelope=0,
   rescaling) make NM dramatically less robust** on this lattice,
   even though they make the objective design more defensible. For
   NM, the original placeholder-laden objectives are well-tuned to
   the optimization trajectory; the fixes shift the basin structure
   into a regime NM cannot handle.

4. **Config C is a useful stress case for global-search optimisers.**
   Stage 11 NM gets stuck in a local minimum on 75% of seeds; whether
   BO recovers that depends on its setup and is not demonstrated by
   this ablation alone. The natural follow-up is to run BO on the
   same A/B/C matrix with a fixed evaluation budget (see open
   git-bug 7f690aa S6) and report the comparison.

5. **Sobol deterministic beam is a useful single-shot probe.** It
   gives ~ median-of-successful-cases for A and B. For C it falls
   in the bad basin, consistent with the multi-seed result. Could
   be a cheap pre-flight check for objective designs before
   committing to multi-seed runs.

**Caveats.**

- Rescaling design (Config B) is one specific choice (only dispersion
  rescaled, ref=0.5 m, weight*4). A different rescaling could change
  B's behaviour. The C result (the typo+envelope fix combination)
  is the more robust finding.
- All starting points are the per-stage defaults from the original
  `UHM_beamline_opt.py`. Warm-chaining stage outputs as next-stage
  starting points may shift the failure rates.
