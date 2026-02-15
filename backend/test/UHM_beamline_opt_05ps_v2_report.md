# UH MkV FEL Beamline Optimization — 0.5 ps Bunch, Paper-Aligned Twiss Targets

Author: Eremey Valetov
Date: 2026-02-06
Script: `backend/test/UHM_beamline_opt_05ps_v2.py`
Exploration: `backend/test/UHM_beamline_opt_05ps_v2_explore.py`

## Motivation

This study revises the 0.5 ps bunch spread optimization to align with the
Twiss matching conditions from Weinberg, Fisher & Li (arXiv:2510.14061v1).
Two corrections are made relative to the previous 0.5 ps study
(`UHM_beamline_opt_05ps.py`):

1. **Energy spread and chirp are unchanged from the 2 ps baseline.**
   The paper (§III, p.5) explicitly states: *"we only change the bunch length
   and keep the other machine and beam parameters the same as the nominal case
   in our simulation to isolate the effect of bunch length shortening."*
   The previous study scaled energy spread from 0.5% to 2.0% by longitudinal
   emittance conservation; this revision keeps it at 0.5%.

2. **Undulator Twiss targets are asymmetric.**
   The paragraph below Table I specifies plane-dependent matching:
   - *Vertical*: beam at waist, collimated through undulator — α_y = 0,
     β_y = γλ_u/(2πK) = 0.24 m
   - *Horizontal*: waist at undulator center, beam size matched to radiation
     spot size — α_x = 0.47, β_x = 1.4 m
   The previous study used symmetric targets (β = 0.24 m, α = 0 in both planes).


## Parameter Changes from Previous 0.5 ps Study

| Parameter | Previous 0.5 ps | This study | Rationale |
|-----------|----------------|------------|-----------|
| `energy_std_percent` | 2.0% | 0.5% | Paper: unchanged from baseline |
| `h` (chirp) | 20e9 /s | 5e9 /s | Paper: unchanged from baseline |
| α_x target | 0 | 0.47 rad | Table I: waist at undulator center |
| β_x target | 0.2418 m | 1.4 m | Table I: radiation mode matching |
| α_y target | 0 | 0 | Unchanged |
| β_y target | 0.2418 m | 0.2418 m | Unchanged (same formula) |


## Beam Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| Energy | 40 MeV | Set by linac |
| RF frequency | 2.856 GHz | |
| Bunch spread | 0.5 ps | Study parameter |
| Energy spread | 0.5% | Unchanged from 2 ps baseline |
| Chirp h | 5e9 /s | Unchanged from 2 ps baseline |
| ε_n | 8 π·mm·mrad | |
| x_std, y_std | 0.8 mm | |
| N particles | 1000 | seed=42 |

### Undulator Matching Targets (Table I)

| Parameter | Vertical | Horizontal | Source |
|-----------|----------|------------|--------|
| β (m) | 0.2418 | 1.4 | γλ_u/(2πK) / radiation mode |
| α (rad) | 0 | 0.47 | waist at entrance / center |
| Matching | Natural undulator focusing | σ_x = σ_r, waist at L_u/2 | Ref. [36] / §III |

The horizontal β_x is determined by matching the electron beam size to the
radiation fundamental mode: σ_x ~ σ_r = √(λ_r Z_R / 4π), β_x = σ_r²/ε_x.
The waist at the undulator center (z = L_u/2 = 0.540 m from entrance) gives
α_x = 0.47 at the entrance.  Consistency check: z_waist = α_x / γ_Twiss =
0.47 / 0.872 = 0.539 m ≈ L_u/2. ✓

### Derived quantities

- γ = 79.279, β_rel = 0.99992
- Geometric emittance ε = 0.1009 π·mm·mrad
- tof_std = 0.5e-12 × 2.856e9 = 1.428e-3
- energy_std = 0.5% × 10 = 5.0


## Optimization Strategy

### Why joint optimization is needed

The undulator matching requires satisfying 4 Twiss objectives (α_x, α_y, β_x,
β_y) at element 117.  With only 3 quadrupole currents in the final triplet
(indices 93, 95, 97), the system is overconstrained.

This overconstraint is **worse for asymmetric targets** than for symmetric ones.
Systematic exploration (see `UHM_beamline_opt_05ps_v2_explore.py`) shows:

| # | Strategy | Variables | MSE | β_x (m) | β_y (m) | α_x | α_y |
|---|----------|-----------|-----|---------|---------|------|------|
| 1 | NM 3-var, symmetric targets | 3 | 5.80e-3 | 0.10 | 0.18 | +0.015 | -0.010 |
| 2 | NM 3-var, asymmetric targets | 3 | 8.59e-2 | 0.97 | 0.51 | +0.185 | -0.037 |
| **3** | **NM 4-var joint, asymmetric** | **4** | **4.26e-9** | **1.40** | **0.24** | **+0.470** | **-0.000** |
| 4 | DiffEvo 3-var, asymmetric | 3 | 8.59e-2 | 0.97 | 0.51 | +0.185 | -0.037 |
| 5 | DiffEvo 4-var joint, asymmetric | 4 | 1.21e-14 | 1.40 | 0.24 | +0.470 | -0.000 |

Key findings:

1. **3-variable approaches are fundamentally overconstrained** regardless of
   optimizer.  Differential evolution (global optimizer) finds the same MSE
   floor as Nelder-Mead (8.59e-2 for asymmetric targets).

2. **Asymmetric targets are harder** than symmetric ones — MSE increases 15×
   (5.8e-3 → 8.6e-2) going from symmetric to asymmetric with 3 variables.

3. **Adding chromaticity quad 5 as a 4th variable resolves the problem.**
   The Nelder-Mead joint approach achieves MSE = 4.3e-9 in 42 seconds.
   DiffEvo confirms this is at the global optimum (MSE = 1.2e-14).


### Adopted approach

Joint Nelder-Mead optimization of chromaticity quad 5 (index 87) and the final
triplet (indices 93, 95, 97), with a secondary dispersion objective at
element 92 (weight 0.5).  All upstream stages use the same settings as the
baseline script.


## Optimization Results

Beamline truncated to 118 elements (MkIII undulator entrance at z = 12.389 m).
All stages use Nelder-Mead.  1000 particles, seed=42.

### Optimized Quadrupole Currents

| Section | Element Index | Description | Current (A) |
|---------|--------------|-------------|-------------|
| First doublet | 1 | LQ1 (QPF) | 0.8226 |
| First doublet | 3 | LQ2 (QPD) | 1.0433 |
| Chromaticity 1 | 10 | DC1 DPHQ | 4.9243 |
| Triplet 1 | 16 | QPF | 2.7164 |
| Triplet 1 | 18 | QPD | 5.0253 |
| Triplet 1 | 20 | QPF | 3.0592 |
| Chromaticity 2 | 27 | DC1 DPHQ | 4.9311 |
| Double triplet | 33 | QPF | 2.4455 |
| Double triplet | 35 | QPD | 3.3353 |
| Double triplet | 37 | QPF | 1.1979 |
| Mirror triplet | 39 | = elem 37 | 1.1979 |
| Mirror triplet | 41 | = elem 35 | 3.3353 |
| Mirror triplet | 43 | = elem 33 | 2.4455 |
| Chromaticity 3 | 50 | FC1 DPHQ | 4.7764 |
| IP doublet | 56 | QPF | 3.1237 |
| IP doublet | 58 | QPD | 3.3253 |
| Post-IP doublet | 61 | QPF | 5.1800 |
| Post-IP doublet | 63 | QPD | 4.0408 |
| Chromaticity 4 | 70 | DPHQ | 4.5854 |
| Triplet 2 | 76 | QPD | 3.9602 |
| Triplet 2 | 78 | QPF | 4.1619 |
| Triplet 2 | 80 | QPD | 0.0000 |
| Chromaticity 5 | 87 | DPHQ (joint) | 1.4668 |
| Final triplet | 93 | QPD | 0.8235 |
| Final triplet | 95 | QPF | 2.8419 |
| Final triplet | 97 | QPD | 2.2091 |

### Twiss Functions at Key Points

#### Chromaticity Quad Results (target: dispersion = 0)

| Section | Element | Current (A) | Residual dispersion |
|---------|---------|-------------|---------------------|
| DC1 (1st) | 15 | 4.9243 | -8.97e-8 |
| DC1 (2nd) | 32 | 4.9311 | -7.58e-7 |
| FC1 (3rd) | 55 | 4.7764 | -1.01e-6 |
| FC2 (4th) | 75 | 4.5854 | -1.79e-8 |
| FC2 (5th, joint) | 92 | 1.4668 | 1.61e-2 |

Chromaticity quads 1–4 converge to effectively zero dispersion.
Chromaticity quad 5 is optimized jointly with the final triplet, trading a
small residual dispersion (0.016) for perfect undulator matching.

#### Interaction Point (z = 7.11 m)

| Twiss parameter | Value |
|-----------------|-------|
| x envelope (elem 59) | 0.0141 mm |
| y envelope (elem 59) | 0.0624 mm |

Both envelopes effectively focused to a waist.

#### Chromaticity 5 + Final Triplet — Undulator Matching (joint, 4 variables)

| Twiss parameter | Target | Achieved |
|-----------------|--------|----------|
| α_x (elem 117) | 0.47 | 0.4700 |
| α_y (elem 117) | 0 | -1.83e-5 |
| β_x (elem 117) | 1.4 m | 1.4000 m |
| β_y (elem 117) | 0.2418 m | 0.2420 m |
| x dispersion (elem 92) | 0 | 0.0161 |

Convergence: 351 iterations, MSE = 2.60e-5.

**All four Twiss parameters match the paper targets to within 0.1%.**


## Comparison with Previous 0.5 ps Study

| Quantity | Previous study | This study |
|----------|---------------|------------|
| Energy spread | 2.0% | 0.5% |
| Chirp | 20e9 /s | 5e9 /s |
| β_x target | 0.2418 m | 1.4 m |
| α_x target | 0 | 0.47 |
| β_x achieved | 0.2419 m | 1.4000 m |
| α_x achieved | -3.49e-5 | 0.4700 |
| β_y achieved | 0.2419 m | 0.2420 m |
| α_y achieved | +3.89e-6 | -1.83e-5 |
| Chrom5 current (A) | 2.1919 | 1.4668 |
| Final triplet (A) | 1.15 / 3.03 / 2.17 | 0.82 / 2.84 / 2.21 |

Both studies achieve excellent undulator matching using the joint 4-variable
approach.  The current study uses the physically correct asymmetric Twiss
targets from the paper, with unchanged energy spread and chirp.


## Comparison with 2 ps Baseline (v2)

The companion script `UHM_beamline_opt_v2.py` optimizes the 2 ps baseline
with the same asymmetric Twiss targets.  Results are very similar:

| Quantity | 2 ps v2 | 0.5 ps v2 |
|----------|---------|-----------|
| β_x achieved | 1.3999 m | 1.4000 m |
| α_x achieved | 0.4700 | 0.4700 |
| β_y achieved | 0.2419 m | 0.2420 m |
| α_y achieved | -3.22e-5 | -1.83e-5 |
| Chrom5 current (A) | 1.3624 | 1.4668 |

The upstream stages are essentially identical since the energy spread is the
same.  Only the final stage shows small current differences, reflecting the
slightly different longitudinal distribution (4× shorter bunch).


## Reproducibility

```bash
cd backend
MPLBACKEND=Agg PYTHONPATH=$(pwd) python test/UHM_beamline_opt_05ps_v2.py
```

Output figures saved to `backend/`:
- `dynamics_plot_z_7.11.eps` — beam envelope evolution
- `phase_space_z_7.11.eps` — phase space at the interaction point

Full strategy exploration:
```bash
MPLBACKEND=Agg PYTHONPATH=$(pwd) python test/UHM_beamline_opt_05ps_v2_explore.py
```

Results are deterministic (seed=42).


## References

- Weinberg, Fisher & Li, *"Three-Dimensional Simulation of the University of
  Hawai'i FEL Oscillator: Superradiant Emission and Cavity Desynchronization"*,
  arXiv:2510.14061v1, Table I and §III.
