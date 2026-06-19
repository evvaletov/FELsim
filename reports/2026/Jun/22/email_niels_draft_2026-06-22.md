DRAFT - review before sending. Outward-facing; not sent autonomously.

To: Niels Bidault <nbidault@hawaii.edu>
Subject: FEL status - early-June action items

Hi Niels,

Status on the items from our early-June meeting. Detail and figures are in the
repo under reports/2026/Jun/22/; the short version:

Settled / answered:
- elegant stays scoped to the RF/linac section. It is our linac reference only
  (0.06% vs RF-Track) and is not needed for the diagnostic chicane, which
  COSY/xsuite/FELsim already cover. I'd suggest pointing elegant at the RF gun
  next as an independent cross-check (the gun is RF-Track-only today) -- let me
  know if you want that.
- PALS does not define a beam-file standard; it is lattice-only. For the
  particle-distribution interchange I recommend openPMD-beamphysics (HDF5), which
  reads ASTRA, elegant/SDDS, and Genesis4, with plain ASCII as a fallback.
  Genesis4 is HDF5 on both ends, so this fits cleanly.
- Linac codes-vs-physics-features table: built (TW/SW, field maps, multi-cell,
  space charge, beam loading, wakefields, longitudinal, autophasing, det(R),
  implemented vs potential per code).
- Fringe-field modes: COSY matches the undulator optics in every mode FR0-FR3,
  each with its own current set. The RF-Track vertical-beta deficit we saw in
  the spring is closed by the edge-fringe correction added since.
- Cross-calibration across FELsim/COSY/RF-Track/xsuite is verified (shared
  gradient calibration, consistent energy/dispersion conventions, agreement to
  the test tolerances).
- glyfada: I had already benchmarked it; its evolutionary search underperforms
  our Nelder-Mead + CMA-ES by several orders on this problem, so it is not the
  production path. The one place it is worth a fresh run is the Config-C stress
  case for the optimisation paper.

In progress, and where I need input from you:
- Objective design: at the MOP6318 targets both objective A and B fail ~35% of
  seeds (C ~80%), so they are not as robust as the abstract framing implied.
  That strengthens the case for the Bayesian-optimisation baseline as the paper
  deliverable. Could you send your xopt hyperparameters so I run the BO baseline
  on the same setup? And confirm you are comfortable with the A/B-not-robust
  framing.
- Linac model refinement: RF-Track is the detailed model (single-Fourier
  uniform-cell). Next step is a cell-by-cell multi-cell travelling-wave model;
  the production tau = 0.57 cell geometry (SLAC-75 Table 6-6) is already in hand,
  so this is a build task on my side. I'll fold a field map in later if needed.
- COSY space charge: Phase 1 is done; next is the three-code DA-FMM + RF-Track +
  xsuite run on the transport line at 45 and 1 MeV, with and without dipoles.

Best,
Eremey
