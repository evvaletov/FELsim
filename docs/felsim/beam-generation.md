# Beam Generation

FELsim generates 6D Gaussian beam distributions at the linac exit for
transport through the beamline. The beam represents a single electron bunch.

## Coordinate System

The 6D phase space vector is $(x, x', y, y', \Delta t, \Delta E / E)$:

| Coordinate | Description | Units |
|------------|-------------|-------|
| $x$ | Horizontal position | m |
| $x'$ | Horizontal divergence | rad |
| $y$ | Vertical position | m |
| $y'$ | Vertical divergence | rad |
| $\Delta t$ | Time offset from bunch center | s |
| $\Delta E / E$ | Fractional energy deviation | — |

## Configurable Parameters

| Parameter | Symbol | Baseline | Range |
|-----------|--------|----------|-------|
| Kinetic energy | $E_k$ | 40 MeV | 20–45 MeV |
| Normalized emittance | $\varepsilon_n$ | 8 π·mm·mrad | 1–20 |
| Bunch length | $\sigma_t$ | 2 ps | 0.5–2 ps |
| Energy spread | $\sigma_E / E$ | 0.5% | 0.1–5% |
| Energy chirp | $h$ | 0 | 0–$40 \times 10^9$ /s |
| Number of particles | $N$ | 1000 | 500–2000 |

## Energy Chirp

The chirp parameter $h$ (units: 1/s) introduces a linear correlation
between time offset and energy deviation:

$$
\frac{\Delta E}{E} = h \cdot \Delta t + \text{uncorrelated spread}
$$

At $h = 5 \times 10^9$ /s the chirp adds 0.25%/$\sigma_t$ to the energy
spread (total $\sigma_E$ grows from 0.50% to 0.56%, a +12% increase).
At $h = 20 \times 10^9$ /s (emittance-conservation scaled), the chirp
more than doubles $\sigma_E$ to 1.12%.

The chirp value $h = 5 \times 10^9$ /s originates from Niels Bidault's
initial beam parameters. It is **not** from the arXiv paper and remains
to be determined by injector simulations. FELsim defaults to $h = 0$.

## Twiss Initialization

The beam is generated with matched Twiss parameters at the linac exit.
These are typically obtained from the first stage of the optimizer or
set to design values. The Twiss parameters at the undulator entrance
(the optimization targets) are:

| Plane | $\beta$ (m) | $\alpha$ |
|-------|-------------|---------|
| $x$ | 1.4 | 0.4714 |
| $y$ | 0.24 | 0.0 |

## Bunch Length Independence

The bunch length parameter does not affect transverse Twiss matching (S6).
This is because FELsim's 6×6 transfer matrices decouple transverse and
longitudinal phase space: the 4×4 transverse block and the dispersion column
(column 6) are independent of the column-5 ($\Delta t$) distribution.
Changing $\sigma_t$ from 0.1 to 2.0 ps produces identical optimization
results.  See also S9 (analytical derivation) and the
[optimization studies](optimization.md#s6--bunch-length-sensitivity).

## Sampling Method

Two sampling methods are supported via `gen_6d_gaussian(..., method=)`:

- **`'random'`** (default): Pseudo-random via `np.random.normal`. Fast, but
  results depend on the random seed.  P12 showed that different seeds produce
  highly variable optimization results at extreme emittances (CV up to 407%).
- **`'sobol'`**: Sobol quasi-random sequence with inverse normal CDF
  (`scipy.stats.qmc.Sobol`).  Deterministic — produces identical output
  regardless of seed.  Requires `num_particles` to be a power of 2.
  Better space-filling properties eliminate the upstream seed variability
  identified in P12/S7/S8. See [P13](optimization.md#p13--deterministic-beam-generation-sobol).

## Random Seed

All optimization scripts use `seed=42` for reproducibility.  With
`method='sobol'`, the seed only affects the Nelder-Mead / CMA-ES optimizer
internals (Stage 11 multi-start), not the beam distribution itself.
