# UH MkV FEL Beamline Optimization — 2 ps Bunch, Paper-Aligned Twiss Targets

Author: Eremey Valetov
Date: 2026-02-06
Script: `backend/test/UHM_beamline_opt_v2.py`

## Motivation

This study revises the baseline 2 ps beamline optimization to use the undulator
Twiss matching conditions from Weinberg, Fisher & Li (arXiv:2510.14061v1,
Table I).  All beam parameters are identical to the original `UHM_beamline_opt.py`;
only the undulator matching targets change.

The original script targeted symmetric β = 0.2418 m and α = 0 in both planes,
using the natural undulator focusing formula β = γλ_u/(2πK) for both x and y.
The paper specifies plane-dependent matching:

| Plane | α | β | Matching condition |
|-------|---|---|-------------------|
| Vertical | 0 | 0.24 m (= γλ_u/(2πK)) | Beam at waist, collimated through undulator |
| Horizontal | 0.47 rad | 1.4 m | Beam matched to radiation mode, waist at undulator center |

The horizontal matching accounts for the radiation fundamental mode, with the
electron beam size matched to the radiation spot size: σ_x = σ_r = √(λ_r Z_R / 4π).


## Beam Parameters (unchanged from original)

| Parameter | Value |
|-----------|-------|
| Energy | 40 MeV |
| RF frequency | 2.856 GHz |
| Bunch spread | 2 ps |
| Energy spread | 0.5% |
| Chirp h | 5e9 /s |
| ε_n | 8 π·mm·mrad |
| x_std, y_std | 0.8 mm |
| N particles | 1000 (seed=42) |


## Optimization Results

Beamline truncated to 118 elements.  All stages use Nelder-Mead.
The final stage uses joint optimization of chromaticity quad 5 (index 87) and
the final triplet (indices 93, 95, 97) — 4 variables for 4 asymmetric Twiss
objectives plus 1 secondary dispersion goal.

### Optimized Quadrupole Currents

| Section | Element Index | Description | Current (A) |
|---------|--------------|-------------|-------------|
| First doublet | 1 | LQ1 (QPF) | 0.8218 |
| First doublet | 3 | LQ2 (QPD) | 1.0430 |
| Chromaticity 1 | 10 | DC1 DPHQ | 3.8834 |
| Triplet 1 | 16 | QPF | 2.2396 |
| Triplet 1 | 18 | QPD | 4.9532 |
| Triplet 1 | 20 | QPF | 3.4258 |
| Chromaticity 2 | 27 | DC1 DPHQ | 4.6657 |
| Double triplet | 33 | QPF | 2.6942 |
| Double triplet | 35 | QPD | 2.6523 |
| Double triplet | 37 | QPF | 0.2768 |
| Mirror triplet | 39 | = elem 37 | 0.2768 |
| Mirror triplet | 41 | = elem 35 | 2.6523 |
| Mirror triplet | 43 | = elem 33 | 2.6942 |
| Chromaticity 3 | 50 | FC1 DPHQ | 4.6739 |
| IP doublet | 56 | QPF | 3.1219 |
| IP doublet | 58 | QPD | 3.3129 |
| Post-IP doublet | 61 | QPF | 5.1775 |
| Post-IP doublet | 63 | QPD | 4.0434 |
| Chromaticity 4 | 70 | DPHQ | 4.6818 |
| Triplet 2 | 76 | QPD | 3.9336 |
| Triplet 2 | 78 | QPF | 4.0787 |
| Triplet 2 | 80 | QPD | 0.0139 |
| Chromaticity 5 | 87 | DPHQ (joint) | 1.3624 |
| Final triplet | 93 | QPD | 0.9452 |
| Final triplet | 95 | QPF | 2.8851 |
| Final triplet | 97 | QPD | 2.1921 |

### Undulator Matching (joint, 4 variables)

| Twiss parameter | Target | Achieved |
|-----------------|--------|----------|
| α_x (elem 117) | 0.47 | 0.4700 |
| α_y (elem 117) | 0 | -3.22e-5 |
| β_x (elem 117) | 1.4 m | 1.3999 m |
| β_y (elem 117) | 0.2418 m | 0.2419 m |
| x dispersion (elem 92) | 0 | 0.0170 |

Convergence: 334 iterations, MSE = 2.89e-5.

**All four Twiss parameters match the paper targets to within 0.1%.**

### Chromaticity Quad Results

| Section | Element | Current (A) | Residual dispersion |
|---------|---------|-------------|---------------------|
| DC1 (1st) | 15 | 3.8834 | -3.88e-7 |
| DC1 (2nd) | 32 | 4.6657 | -1.94e-7 |
| FC1 (3rd) | 55 | 4.6739 | 7.46e-9 |
| FC2 (4th) | 75 | 4.6818 | 1.06e-6 |
| FC2 (5th, joint) | 92 | 1.3624 | 1.70e-2 |

### Interaction Point (z = 7.11 m)

| Twiss parameter | Value |
|-----------------|-------|
| x envelope (elem 59) | 0.0190 mm |
| y envelope (elem 59) | 0.0560 mm |


## Comparison with Original Baseline

| Quantity | Original | This study |
|----------|----------|------------|
| β_x target | 0.2418 m | 1.4 m |
| α_x target | 0 | 0.47 |
| β_y target | 0.2418 m | 0.2418 m |
| Final stage | chrom5 separate + 3-var triplet | Joint 4-var |

The original script used the natural undulator focusing formula for both planes,
giving symmetric targets.  With 3 variables and 4 objectives (3 for symmetric
case where α_x ≈ α_y), the original approach was marginally constrained.
The asymmetric targets from the paper require 4 independent objectives, making
the joint approach necessary.


## Reproducibility

```bash
cd backend
MPLBACKEND=Agg PYTHONPATH=$(pwd) python test/UHM_beamline_opt_v2.py
```

Results are deterministic (seed=42).


## References

- Weinberg, Fisher & Li, arXiv:2510.14061v1, Table I and §III.
