# UH MkV FEL Beamline Optimization — 0.5 ps Bunch Spread Study

Author: Eremey Valetov
Date: 2026-02-06
Script: `backend/test/UHM_beamline_opt_05ps.py`
Exploration: `backend/test/UHM_beamline_opt_05ps_explore.py`

## Motivation

This study replicates the baseline UH MkV FEL beamline optimization
(`UHM_beamline_opt.py`, 2 ps bunch spread) for a 4x shorter bunch at 0.5 ps.
Shorter bunches are relevant for improved FEL gain and higher peak current
but come with larger energy spread, making chromatic control more demanding.

## Parameter Changes

| Parameter | Baseline (2 ps) | This study (0.5 ps) | Rationale |
|-----------|-----------------|---------------------|-----------|
| `bunch_spread` | 2 ps | 0.5 ps | Study parameter |
| `energy_std_percent` | 0.5% | 2.0% | Longitudinal emittance conservation: sigma_t * sigma_E ~ const; 4x compression -> 4x energy spread |
| `h` (chirp) | 5e9 /s | 20e9 /s | Energy-time correlation steepens proportionally with compression ratio |
| `Energy` | 40 MeV | 40 MeV | Unchanged — set by linac |
| `epsilon_n` | 8 pi.mm.mrad | 8 pi.mm.mrad | Unchanged — transverse emittance decouples from longitudinal |
| `x_std`, `y_std` | 0.8 mm | 0.8 mm | Unchanged |
| Undulator K | 1.2 | 1.2 | Unchanged — hardware parameter |
| Undulator period | 2.3 cm | 2.3 cm | Unchanged — hardware parameter |

### Why energy but not energy spread?

The beam energy (40 MeV) is set by the linac accelerating gradient and is
independent of the bunch length. The energy *spread*, however, is coupled to
the bunch length through longitudinal phase space conservation. When a bunch
is compressed (e.g. by velocity bunching or magnetic compression), the
longitudinal emittance sigma_t * sigma_E is approximately conserved. A 4x
shorter bunch therefore requires a 4x larger energy spread.

The chirp parameter `h` represents the linear energy-time correlation within
the bunch. For a more strongly compressed bunch, this correlation steepens
proportionally.

### Derived quantities (unchanged)

- gamma = 79.279
- beta = 0.99992
- Geometric emittance epsilon = 0.1009 pi.mm.mrad
- Matched beta_y at undulator = 0.2418 m
- Matched y_std at undulator = 0.156 um
- Matched y'_std at undulator = 0.646 mrad

### Derived quantities (changed by bunch parameters)

- tof_std = bunch_spread * 1e-9 * f = 0.5e-12 * 2.856e9 = 1.428e-3 (in RF period units)
- energy_std = energy_std_percent * 10 = 20.0 (in delta_W/W * 1e3 units)

For the baseline: tof_std = 5.712e-3, energy_std = 5.0.
The 4x larger energy spread is the dominant change affecting optics.


## Optimization Strategy

### Initial approach: stage-by-stage Nelder-Mead

The baseline script optimizes the beamline sequentially, stage by stage. The
final stage optimizes 3 quadrupole currents (indices 93, 95, 97) to match 4
Twiss parameters (x-alpha, y-alpha, x-beta, y-beta) at the undulator entrance.
This is an overconstrained problem (3 variables, 4 objectives), and with the
0.5 ps energy spread, the initial Nelder-Mead approach achieved poor x-beta
matching (0.083 m vs 0.242 m target).

### Strategy exploration

To resolve this, 13 alternative strategies were systematically evaluated
(see `UHM_beamline_opt_05ps_explore.py`), including:

| # | Strategy | Variables | MSE | x beta (m) | y beta (m) |
|---|----------|-----------|-----|------------|------------|
| 1 | Nelder-Mead baseline (3 var) | 3 | 6.45e-3 | 0.083 | 0.223 |
| 2 | Powell (3 var) | 3 | 2.58e-1 | 0.052 | 1.219 |
| 3 | L-BFGS-B (3 var) | 3 | 9.68e-3 | 0.228 | 0.046 |
| 4 | COBYLA (3 var) | 3 | 5.90e-2 | 0.614 | 0.146 |
| 5 | Multi-start NM, 20 trials (3 var) | 3 | 6.45e-3 | 0.083 | 0.223 |
| 6 | NM with x-beta weight=3 (3 var) | 3 | 9.78e-3 | 0.239 | 0.047 |
| **7** | **Joint chrom5 + final triplet NM (4 var)** | **4** | **1.56e-9** | **0.242** | **0.242** |
| 8 | NM wide bounds 0–15 A (3 var) | 3 | 6.45e-3 | 0.083 | 0.223 |
| 9 | Joint triplet2 + chrom5 + final NM (7 var) | 7 | 6.93e-10 | 0.242 | 0.242 |
| 10 | Differential evolution (3 var) | 3 | 8.75e-3 | 0.055 | 0.251 |
| 11 | DiffEvo joint chrom5 + triplet (4 var) | 4 | 6.44e-21 | 0.242 | 0.242 |

Key findings:

1. **The 3-variable problem is fundamentally overconstrained.** No 3-variable
   method — regardless of algorithm, starting point, or bounds — achieves
   simultaneous x-beta and y-beta matching. Strategies 1–6, 8, 10 all show
   one plane matched at the expense of the other.

2. **Adding chromaticity quad 5 as a 4th variable resolves the problem.**
   Strategy 7 (Nelder-Mead, 4 variables: chrom5 + final triplet) achieves
   essentially perfect matching (MSE ~ 1e-9) in only 41 seconds. The extra
   degree of freedom from the chromaticity quad provides the knob needed to
   balance both planes simultaneously.

3. **Global optimization confirms the result.** Strategy 11 (differential
   evolution, 4 variables) converges to MSE = 6.4e-21, confirming the global
   optimum. The NM solution (strategy 7) is already at this optimum.

### Adopted approach

The final script uses joint optimization of chromaticity quad 5 (index 87) and
the final triplet (indices 93, 95, 97) with Nelder-Mead. The 5th chromaticity
section's dispersion objective is included with reduced weight (0.5) to allow
the optimizer to trade off a small amount of residual dispersion for better
undulator matching.


## Optimization Results

The beamline was truncated to 118 elements (up to the MkIII undulator entrance
at z = 12.389 m). All stages used the Nelder-Mead algorithm.
1000 particles were tracked.

### Optimized Quadrupole Currents

| Section | Element Index | Description | Current (A) |
|---------|--------------|-------------|-------------|
| First doublet | 1 | LQ1 (QPF) | 0.8373 |
| First doublet | 3 | LQ2 (QPD) | 1.0513 |
| Chromaticity 1 | 10 | DC1 DPHQ | 3.6517 |
| Triplet 1 | 16 | QPF | 2.0777 |
| Triplet 1 | 18 | QPD | 4.9116 |
| Triplet 1 | 20 | QPF | 3.5105 |
| Chromaticity 2 | 27 | DC1 DPHQ | 4.5599 |
| Double triplet | 33 | QPF | 2.6824 |
| Double triplet | 35 | QPD | 2.3708 |
| Double triplet | 37 | QPF | 0.0000 |
| Mirror triplet | 39 | = elem 37 | 0.0000 |
| Mirror triplet | 41 | = elem 35 | 2.3708 |
| Mirror triplet | 43 | = elem 33 | 2.6824 |
| Chromaticity 3 | 50 | FC1 DPHQ | 4.6328 |
| IP doublet | 56 | QPF | 3.1212 |
| IP doublet | 58 | QPD | 3.3086 |
| Post-IP doublet | 61 | QPF | 5.1735 |
| Post-IP doublet | 63 | QPD | 4.0436 |
| Chromaticity 4 | 70 | DPHQ | 4.7364 |
| Triplet 2 | 76 | QPD | 3.9270 |
| Triplet 2 | 78 | QPF | 4.1579 |
| Triplet 2 | 80 | QPD | 0.1858 |
| Chromaticity 5 | 87 | DPHQ (joint) | 2.1919 |
| Final triplet | 93 | QPD | 1.1539 |
| Final triplet | 95 | QPF | 3.0329 |
| Final triplet | 97 | QPD | 2.1724 |

### Twiss Functions at Optimization Points

#### First Quadrupole Doublet (target: alpha = 0)

| Twiss parameter | Value |
|-----------------|-------|
| x alpha (elem 8) | -7.08e-4 |
| x beta (elem 8) | 1.004 m |
| y alpha (elem 9) | -5.00e-4 |
| y beta (elem 9) | 0.127 m |

Convergence: 61 iterations, residual 9.42e-5.

#### Chromaticity Quad Results (target: dispersion = 0)

| Chromaticity section | Element | Current (A) | Residual dispersion |
|---------------------|---------|-------------|---------------------|
| DC1 (1st) | 15 | 3.6517 | -6.95e-8 |
| DC1 (2nd) | 32 | 4.5599 | -6.25e-7 |
| FC1 (3rd) | 55 | 4.6328 | -1.36e-6 |
| FC2 (4th) | 75 | 4.7364 | 4.42e-7 |
| FC2 (5th, joint) | 92 | 2.1919 | 1.15e-2 |

Chromaticity quads 1–4 converge to effectively zero dispersion. Chromaticity
quad 5 (index 87) is optimized jointly with the final triplet, so it trades
a small residual dispersion (0.012) for the undulator matching improvement.

#### Quadrupole Triplet 1 (target: alpha = 0, beta ~ 0.1 m)

| Twiss parameter | Value |
|-----------------|-------|
| x alpha (elem 25) | -2.63e-4 |
| x beta (elem 25) | 0.0155 m |
| y alpha (elem 26) | -9.61e-5 |
| y beta (elem 26) | 0.0121 m |

Convergence: 96 iterations, residual 1.86e-3.

#### Double Triplet (target: alpha = 0, envelope ~ 2 mm)

| Twiss parameter | Value |
|-----------------|-------|
| x alpha (elem 37) | -0.00862 |
| y alpha (elem 37) | -0.00528 |
| x envelope (elem 37) | 1.290 mm |
| y envelope (elem 37) | 1.478 mm |

Convergence: 265 iterations, residual 1.94e-1. Envelopes below
the 2.0 mm target due to the larger energy spread.

#### Interaction Point Doublet (target: envelope -> 0)

| Twiss parameter | Value |
|-----------------|-------|
| x envelope (elem 59) | 0.0228 mm |
| y envelope (elem 59) | 0.0511 mm |

Convergence: 83 iterations, residual 1.56e-3.
Both envelopes effectively focused to a waist at z = 7.11 m.

#### Post-IP Doublet (target: alpha = 0, beta ~ 0.1 m)

| Twiss parameter | Value |
|-----------------|-------|
| x alpha (elem 68) | 0.00108 |
| x beta (elem 68) | 0.0628 m |
| y alpha (elem 69) | 2.13e-5 |
| y beta (elem 69) | 0.0150 m |

Convergence: 96 iterations, residual 1.08e-3.

#### Quadrupole Triplet 2 (target: alpha = 0, beta ~ 0.1 m)

| Twiss parameter | Value |
|-----------------|-------|
| x alpha (elem 85) | 2.12e-4 |
| x beta (elem 85) | 0.0774 m |
| y alpha (elem 86) | 2.45e-4 |
| y beta (elem 86) | 0.0839 m |

Convergence: 157 iterations, residual 9.62e-5.

#### Chromaticity 5 + Final Triplet — Undulator Matching (joint, 4 variables)

Target: alpha = 0 in both planes, beta = 0.2418 m in both planes,
with secondary dispersion goal at element 92.

| Twiss parameter | Target | Achieved |
|-----------------|--------|----------|
| x alpha (elem 117) | 0 | -3.49e-5 |
| y alpha (elem 117) | 0 | +3.89e-6 |
| x beta (elem 117) | 0.2418 m | 0.2419 m |
| y beta (elem 117) | 0.2418 m | 0.2419 m |
| x dispersion (elem 92) | 0 | 0.0115 |

Convergence: 232 iterations, residual 1.32e-5.

**Both x and y beta functions match the undulator target to within 0.01%.**
Alpha functions are effectively zero in both planes. This is a dramatic
improvement over the initial 3-variable approach, which achieved x-beta of
only 0.083 m (66% below target).


## Discussion

1. **Chromaticity correction is robust.** The first four chromaticity quads
   drive dispersion to effectively zero with currents in the 3.7–4.7 A range.
   The 5th chromaticity quad is partially repurposed for undulator matching
   (current reduced from ~3.9 A to 2.2 A), resulting in a small residual
   dispersion of 0.012 at element 92. This trade-off is acceptable because
   the undulator is downstream and the FEL interaction depends primarily on
   the beam Twiss parameters at entry, not dispersion at element 92.

2. **Interaction point focus is preserved.** The IP doublet achieves
   sub-0.06 mm envelopes in both planes — the waist at z = 7.11 m is
   essentially unaffected by the energy spread change.

3. **Undulator matching is achieved** by jointly optimizing chromaticity
   quad 5 with the final triplet, providing 4 variables for 4 objectives
   (+ 1 secondary objective). This resolves the overconstrained 3-variable
   limitation of the stage-by-stage approach. The key insight is that the
   chromaticity quad provides an additional knob that couples into the
   transverse optics downstream.

4. **Total optimization time:** ~136 s (1148 total iterations across 10
   optimization stages), faster than the original sequential approach.


## Reproducibility

To reproduce these results:

```bash
cd backend
MPLBACKEND=Agg PYTHONPATH=$(pwd) python test/UHM_beamline_opt_05ps.py
```

Output figures are saved to `backend/`:
- `dynamics_plot_z_7.11.eps` — beam envelope evolution along the beamline
- `phase_space_z_7.11.eps` — phase space at the interaction point

To reproduce the full strategy exploration:

```bash
cd backend
MPLBACKEND=Agg PYTHONPATH=$(pwd) python test/UHM_beamline_opt_05ps_explore.py
```

Note: results depend on the random seed for particle generation (1000 particles,
6D Gaussian). The exploration script uses `np.random.seed(42)` for deterministic
reproducibility.
