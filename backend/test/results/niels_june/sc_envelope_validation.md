# sigma(s) envelope for per-cell linac space charge -- finding

Author: Eremey Valetov, 2026-06-20.

Investigated the roadmap item 'SC uses a fixed input sigma; feed the matched envelope'. Two findings change the picture:

1. **The per-cell SC sigma(s) is already self-consistent.** For a tracked bunch the frozen SC runs with `update_on_track=True`, so each cell recomputes sigma from the actual beam -- it is not fixed at the injection value. The build-time sigma is only an initial/mesh value.
2. **The beam does not adiabatically shrink here.** Pre-pass envelope (87 cells): sigma_x cell1=1.000, mid=1.037, last=1.058 mm -- roughly constant. This SLAC section is focusing-free, so the naive adiabatic 1/sqrt(beta gamma) law (which would predict 0.19 mm at the end) is WRONG; the envelope is set by the optics + SC, which the self-consistent track captures.

A new OPTIONAL prescribed-envelope path was added (`sig_env` + `sc_envelope_prepass`) for matched lattices and deterministic/PIC-mesh use. It is verified to apply the per-cell sigma correctly. But it is NOT the default and is LESS accurate for an unmatched beam: prescribed (fixed sigma) gives eps_n 5.80 mm.mrad / sigma_x 1.82 mm vs self-consistent 2.78 / 1.06 mm. A fixed sigma does not self-correct, so once the unmatched beam grows it over-defocuses; `update_on_track` self-limits and is the correct default. For a matched focusing lattice the two converge.

**Conclusion:** keep `update_on_track` (self-consistent) as the default -- the SC sigma(s) is already correct. The prescribed envelope is a documented option for matched lattices / reproducible studies. Remaining: `simulate()` exit-energy frame; revisit the prescribed path on the real injector+linac (with matching quads), where a matched envelope exists.