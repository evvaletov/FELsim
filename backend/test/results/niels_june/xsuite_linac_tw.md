# xtrack-native multi-cell TW linac

Author: Eremey Valetov, 2026-06-20. The validated autophasing TW model
(linac_multicell_tw.py) drives an xtrack Line via per-cell Cavity +
ReferenceEnergyIncrease (energy-ramp reference), so xtrack can track the
1->40 MeV linac it otherwise cannot, and provides the 6D transport for the
xsuite space-charge path.

| Check | xtrack | reference | agreement |
|---|---:|---:|---:|
| on-axis K_out (MeV) | 40.953 | 40.953 (integrator) | 0.000% |
| transverse det(R_x) | 0.0343 | 0.0343 (p_in/p_out) | 0.00% |
| final delta (synchronous) | -7.2e-16 | 0 | - |

86 cells, optimal injection phase 339 deg. The reference
ramps with the synchronous energy (delta stays ~0), and the transverse map
shows the correct adiabatic damping det(R_x)=p_in/p_out. This Line is a
drop-in for the xsuite SC engines (frozen / PIC). Beam loading is modelled
in this standalone validation script via synchronous_profile(..., I_amp=...);
the production XsuiteAdapter cavity is currently unloaded (zero beam current).
Remaining: true-phase cavities for self-consistent longitudinal bunching
(this build is on-crest, energy-exact).