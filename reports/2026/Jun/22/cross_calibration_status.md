# Cross-calibration of the FELsim multi-code framework

Author: Eremey Valetov, 2026-06-19

## Bottom line

The calibration backbone is in place and verified across all four codes
(FELsim, COSY INFINITY, RF-Track, xsuite): they consume the same magnet
calibration constant, the same current->strength formula, and consistent
energy/dispersion conventions, and they agree on matched Twiss to the tolerances
below. The remaining cross-calibration work is the space-charge-on, full
transport-line three-code run (the 2026-05-04 capstone), not the linear optics.

## What is calibrated and verified

**Magnet strength.** One anchor, `G_quad_default = 2.694 T/A/m`
(`physicalConstants.py`). Every adapter maps current to normalised gradient with
the identical formula

    k1 = |Q . G . I| / (m . c . beta . gamma)

in FELsim (`beamline.py`), RF-Track (`rftrackAdapter.py`), and xsuite
(`xsuiteAdapter.py`); COSY consumes the same 2.694 T/A/m as a field gradient in
its `MQ` commands. G is configurable per adapter for sensitivity studies.

**Energy / dispersion convention.** Handled explicitly by the
`CoordinateTransformer` (`simulatorFactory.py`):
- FELsim <-> COSY uses delta = (K - K0)/K0 (kinetic) with the gamma/(1+gamma)
  longitudinal factor;
- FELsim <-> xsuite uses delta = (p - p0)/p0 (no gamma factor).
- Dispersion converts as D_standard = D_cosy . (p0c / K0); verified against MAD-X
  to ~1e-10 (2026-02-24).

**Agreement (verification numbers).**
- xsuite adapter calibration (`test_xsuite_adapter.py`): coordinate round-trip
  5.5e-13, drift vs FELsim 7.8e-10, stable FODO beta vs FELsim 7.5e-5;
  `test_multicode` 38/38.
- FELsim vs RF-Track (`test_cross_code_benchmark.py`): beta within 10% in
  drift/quad regions, y-emittance conserved both codes.
- Three-code matched Twiss at eps_n = 8 (R2, 2026-03-04): FELsim NM MSE 1.2e-6,
  COSY FR0 2.3e-7, COSY FR1 3.7e-9, RF-Track opt 7.0e-3 (after 9 cross-code bug
  fixes). Each code reaches the same Twiss with its own currents.
- Linac energy gain RF-Track vs elegant: 0.06% at the peak.

## Remaining

1. Space-charge-on three-code run (DA-FMM + RF-Track + xsuite) on the
   focusing/transport line at 45 MeV and 1 MeV, with and without dipoles
   (2026-05-04 item l). FODO-level previews exist (F1-F7); the full line run does
   not. This is where calibration is the prerequisite, not the deliverable.
2. xsuite has no dipole edge/fringe and no RF-cavity model, so any +dipole or
   linac comparison diverges there until those are added (the bend body now
   builds after the 2026-06-19 h->angle fix).

## See also

- `reports/2026/Jun/22/FF1_fringe_mode_comparison.md` (fringe-mode side).
- git-bug `f59e144` item l, `343b42f` (RF-Track <-> xsuite linac SC).
