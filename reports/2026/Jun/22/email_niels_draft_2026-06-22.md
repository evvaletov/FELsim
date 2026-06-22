DRAFT - review before sending. Outward-facing; not sent autonomously.

To: Niels Bidault <nbidault@hawaii.edu>
Subject: FEL status - early-June action items

Hi Niels,

Status on the items from our early-June meeting. Detail is in the repo under
reports/2026/Jun/22/, the key plots are attached; the short version:

Settled / answered:
- elegant stays scoped to the RF/linac section. It is the linac reference only
  (0.06% vs RF-Track) and is not needed for the diagnostic chicane section, which
  COSY/xsuite/FELsim already cover. For an independent gun cross-check a dedicated
  photoinjector code (ASTRA or GPT) is probably a better fit than elegant, since
  the low-energy, space-charge-dominated gun is outside the regime elegant is
  built for; the gun is RF-Track-only today.
- PALS does not define a beam-file standard; it is lattice-only. For the
  particle-distribution interchange I recommend openPMD-beamphysics (HDF5), which
  reads ASTRA, elegant/SDDS, and Genesis4, with plain ASCII as a fallback.
  Genesis4 is HDF5 on both ends anyway.
- Linac codes-vs-physics-features table: built (TW/SW, field maps, multi-cell,
  space charge, beam loading, wakefields, longitudinal, autophasing, det(R),
  implemented vs potential per code).
- Fringe-field modes: COSY matches the undulator optics in every mode FR0-FR3,
  each with its own current set. The RF-Track vertical-beta deficit we saw in
  the spring is closed by the edge-fringe correction added since.
- The cross-calibration backbone across FELsim/COSY/RF-Track/xsuite is verified
  for the linear optics and RF (shared gradient calibration, consistent
  energy/dispersion conventions, agreement to the test tolerances); the
  with-dipoles space-charge line is the one piece still open (see below).
- glyfada: following your suggestion I re-checked the benchmark: its evolutionary
  search underperforms the current Nelder-Mead + CMA-ES workflow by several orders
  of magnitude on this objective, so it is not the right tool for the production
  matching. I would still rerun the Config-C stress case for the optimisation paper.

Recent work, and where I need input from you:
- Objective design: in these seed tests both objective A and B fail ~35% of seeds
  at the MOP6318 targets (C fails ~80%), which points to a robustness limitation
  there. I think this makes the Bayesian-optimisation baseline the right
  comparison point for the paper. Two things would help: (i) the xopt
  hyperparameters you use (acquisition, kernel, per-stage eval budget) so I can
  run the BO baseline on an identical setup, and (ii) how you would prefer to
  frame the A/B robustness. Could you also confirm the MOP6318 Twiss targets (I
  have beta_x = 1.267 m, alpha_x = 0.560 from a code comment and want to verify
  them against the abstract)?
- Linac model refinement: done. On top of the RF-Track model I built a
  cell-by-cell multi-cell travelling-wave linac (tau = 0.57, SLAC-75 Table 6-6),
  validated to 1.18% against elegant for the standalone cell model (the xsuite
  production adapter, integrating the full structure, matches to 0.06%), and added
  steady-state beam loading. The main later refinement I see is a
  measured/simulated field map.
- COSY space charge: the comparison on the no-dipole section of the transport line
  is set up, with the first cross-code results there (earlier checks were on a toy
  FODO). The codes now share one macroparticle distribution (with a
  reproducibility manifest), and at zero current they agree to <0.1% on that
  focusing section, so they are calibrated. The full line is not yet comparable:
  xsuite has no dipole edge/fringe model and the dipole/dispersion handling
  diverges between codes, so the with-dipoles run needs that edge/fringe model
  plus SC-inside-magnet element slicing first (the DA-FMM kicks currently act only
  in drifts). With space charge on the no-dipole section I have three engines side
  by side: the cosy-fmm DA-FMM treecode, xsuite frozen-Gaussian, and the new
  cosy-pic mesh PIC solver (RF-Track's PIC would be the fourth, but it is blocked
  by a 2.5.5 core-dump, so it is out of this comparison for now). Two results
  stand out. At 45 MeV the cosy-pic mesh PIC and xsuite frozen-Gaussian agree
  within ~20% across the charge scan, while the bare treecode over-predicts the
  emittance growth by 6-22x over the 1-10 nC scan (largest at low charge; at the
  nominal 0.3 nC all three show sub-1% growth, in the noise floor). That excess is
  consistent with shot noise from too few macroparticles (the threshold we flagged
  in the spring) rather than physics, so cosy-pic and xsuite frozen-Gaussian give
  the best physics estimate here. At
  1 MeV the ordering reverses: there the frozen-Gaussian over-predicts emittance
  growth ~5x versus the DA-FMM treecode (cosy-pic was run only at 45 MeV),
  consistent with the frozen model re-imposing a Gaussian field instead of letting
  the distribution relax. Write-up and figures are in
  backend/test/results/sc_capstone/.

Needed from you:
- the gun cross-check preference (ASTRA/GPT, or run elegant anyway);
- the xopt hyperparameters for the BO baseline, plus how you would like to frame
  the A/B robustness and a quick confirm of the MOP6318 Twiss targets.

Best,
Eremey
