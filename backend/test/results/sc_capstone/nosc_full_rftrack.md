# No-SC cross-code handoff — cold-beam agreement

- Distribution: N_p=6000, E=45.0 MeV, eps_n=8.0 mm.mrad, seed=20260619, sha256=b02d509d4426c771
- Codes: felsim, rftrack (reference: felsim)
- Space charge: OFF

### Full transport line (0, 137) (cold beam; xsuite has no dipole edge/fringe -> divergence expected)

| quantity | initial | felsim | rftrack |
|---|---|---|---|
| sig_x [mm] | 0.52598 | 1278.6 | 26.969 |
| sig_y [mm] | 0.51812 | 2.5627 | 1.2081 |
| eps_n,x [mm.mrad] | 8.1645 | 321.32 | 72.195 |
| eps_n,y [mm.mrad] | 7.9511 | 7.9511 | 61.769 |
| centroid_x [mm] | 0.010504 | -6.1055 | 0.49997 |
| centroid_y [mm] | 0.0025607 | 0.0051428 | -0.11342 |

**Agreement vs felsim (relative):**

| quantity | rftrack |
|---|---|
| sig_x [mm] | 97.891% |
| sig_y [mm] | 52.858% |
| eps_n,x [mm.mrad] | 77.532% |
| eps_n,y [mm.mrad] | 676.860% |

Worst sigma/emittance disagreement vs felsim: **676.860%** (tolerance 2% -> EXCEEDS)

Wall clock: felsim 0.01s, rftrack 1.09s
