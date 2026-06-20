# XsuiteAdapter linac integration

Author: Eremey Valetov, 2026-06-20

The multi-code framework's `XsuiteAdapter` previously treated `RF_CAVITY` as a
drift, so xsuite could not run the linac. It now builds the accelerating
structure, so the framework (and the space-charge engines) can act through the
S-band linac.

## What changed (`backend/xsuiteAdapter.py`)

- `RF_CAVITY` builds a multi-cell travelling-wave chain (per-cell `xt.Cavity`
  on crest + `xt.ReferenceEnergyIncrease`) from the element's
  `gradient_mv_per_m` / `frequency_hz` / `phase_advance_deg`. The per-cell energy
  gains come from the autophasing model (`_tw_synchronous_profile`), so the
  reference momentum ramps with the synchronous particle.
- Quad `k1` is computed at the local reference energy: a magnet after the linac
  is weaker (`k1 ~ 1/p`), so the running energy is threaded through `_build_line`.
- RF-cavity parameters now survive `_convert_element_from_native` (they were
  dropped before).

## Validation (`backend/test/rftrack_linac/test_xsuite_linac_adapter.py`)

| Check | Result |
|---|---|
| `slac_linac.json` acceleration | 1 MeV -> 41.469 MeV (delta ~ 1e-15); matches RF-Track 41.468 to 0.001 MeV, elegant 41.442 to 0.06% |
| energy-aware k1 | quad after cavity uses k1=96 (post-linac), not the 1 MeV value 2840 |
| element length | TW chain spans 3.048 m exactly (L_cell = L/n_cells, no downstream displacement) |
| RF param guard | cavity missing frequency_hz falls back to a drift |
| transport-only line | unchanged vs the original fixed-energy build (backward compatible) |

## Remaining

- Space charge currently treats the cavity as a block; per-cell SC inside the
  linac (with local gamma) is the next step for the SC-on-linac capstone.
- `simulate()` uses one reference energy for the in/out coordinate transform;
  a linac that changes energy needs the exit energy for the output frame.
- COSY TW DA map (the third linac code).
