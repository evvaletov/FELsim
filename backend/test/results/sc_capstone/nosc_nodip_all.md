# No-SC cross-code handoff — cold-beam agreement

- Distribution: N_p=6000, E=45.0 MeV, eps_n=8.0 mm.mrad, seed=20260619, sha256=b02d509d4426c771
- Codes: felsim, xsuite, rftrack (reference: felsim)
- Space charge: OFF

### No-dipole subset (32, 46) (drift+quad, cold beam)

| quantity | initial | felsim | xsuite | rftrack |
|---|---|---|---|---|
| sig_x [mm] | 0.52598 | 0.27581 | 0.27593 | 0.27593 |
| sig_y [mm] | 0.51812 | 0.24857 | 0.24853 | 0.24853 |
| eps_n,x [mm.mrad] | 8.1645 | 8.1645 | 8.1702 | 8.1702 |
| eps_n,y [mm.mrad] | 7.9511 | 7.9511 | 7.9505 | 7.9505 |
| centroid_x [mm] | 0.010504 | -0.0062432 | -0.0064016 | -0.0064016 |
| centroid_y [mm] | 0.0025607 | -0.001627 | -0.0016813 | -0.0016813 |

**Agreement vs felsim (relative):**

| quantity | xsuite | rftrack |
|---|---|---|
| sig_x [mm] | 0.041% | 0.041% |
| sig_y [mm] | 0.018% | 0.018% |
| eps_n,x [mm.mrad] | 0.070% | 0.070% |
| eps_n,y [mm.mrad] | 0.007% | 0.007% |

Worst sigma/emittance disagreement vs felsim: **0.070%** (tolerance 2% -> PASS)

Wall clock: felsim 0.00s, xsuite 0.68s, rftrack 0.23s
