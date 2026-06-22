# Per-cell space charge in the xsuite TW linac

Author: Eremey Valetov, 2026-06-20.

1. **Structure:** 87 SC kicks interleaved, one per TW cell, each right after that cell's ReferenceEnergyIncrease (so it acts at the local energy).
2. **Acceleration:** with SC on, a 4000-particle bunch still reaches 41.4 MeV and all 4000 particles survive.
3. **Local-gamma law:** the frozen SC defocus dpx/x falls from 1.07e+01 1/m at 1 MeV to 4.42e-04 1/m at 41.5 MeV (**x24265**), tracking the 1/(beta^2 gamma^3) fixed-charge relativistic shielding (predicted x24265). A single fixed-energy SC block would apply the 1 MeV strength throughout and over-estimate SC by this factor.

Frozen SC reads the ramped reference automatically; PIC gets the local gamma. Emittance is a poor probe here (the field is nearly linear across a round Gaussian -> envelope perturbation, not eps_n growth, in one short pass). Remaining: sigma(s) envelope evolution (fixed input sigma here); simulate() exit-energy coordinate frame.