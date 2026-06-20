# No-SC cross-code handoff — cold-beam agreement

- Distribution: N_p=6000, E=45.0 MeV, eps_n=8.0 mm.mrad, seed=20260619, sha256=b02d509d4426c771
- Codes: felsim, xsuite (reference: felsim)
- Space charge: OFF

### Full transport line (0, 137) (cold beam; xsuite has no dipole edge/fringe -> divergence expected)

| quantity | initial | felsim |
|---|---|---|
| sig_x [mm] | 0.52598 | 1278.6 |
| sig_y [mm] | 0.51812 | 2.5627 |
| eps_n,x [mm.mrad] | 8.1645 | 321.32 |
| eps_n,y [mm.mrad] | 7.9511 | 7.9511 |
| centroid_x [mm] | 0.010504 | -6.1055 |
| centroid_y [mm] | 0.0025607 | 0.0051428 |


**Codes that did not complete:**
- `xsuite`: {'failed_section': 'xsuite', 'section_index': 0}

Wall clock: felsim 0.03s
