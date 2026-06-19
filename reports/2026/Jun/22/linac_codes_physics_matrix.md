# Linac simulation codes x physics features

Eremey Valetov, 2026-06-19

RF-Track is the only code in the FELsim stack with a detailed, cross-validated
linac model: an analytical single-Fourier `TW_Structure` benchmarked to 0.06%
against elegant at the energy-gain peak, with adiabatic-damping det(R) confirmed
against the elegant RFCA reference. elegant serves as the independent RFCA+TWLA
reference. COSY INFINITY has a Maxwell-consistent standing-wave cavity
(`RFCAVX_TM`, TM010 pillbox, 491-line passing test suite) but it is not yet wired
into FELsim's COSYAdapter, and its travelling-wave extension is designed, not
coded (N28). xsuite models the linac only as a lumped single `Cavity` in a
standalone comparison; its production adapter treats `RF_CAVITY` as a drift.
FELsim's 1st-order transfer-matrix core has no RF element at all.

This table is codes x linac PHYSICS FEATURES. It is distinct from the P1
capability matrix (`backend/test/generate_ipac_figures.py:p1_capability_matrix`),
which is codes x beamline BLOCKS. Each cell is tagged by adoption status in the
FELsim stack, not by what the upstream code could do in principle.

## Matrix

| Feature \ Code                     | RF-Track | elegant | COSY INFINITY | xsuite | FELsim 1st-order | analytical |
|------------------------------------|:--------:|:-------:|:-------------:|:------:|:----------------:|:----------:|
| TW structure                       |   [x]    |   [x]   |      [o]      |  [o]   |       [-]        |    [~]     |
| SW structure                       |   [~]    |   [~]   |      [~]      |  [~]   |       [-]        |    [o]     |
| RF field maps                      |   [~]    |   [o]   |      [o]      |  [o]   |       [-]        |    [-]     |
| Multi-cell / cell-by-cell geometry |   [~]    |   [~]   |      [o]      |  [o]   |       [-]        |    [o]     |
| Space charge                       |   [~]    |   [o]   |      [x]      |  [~]   |       [-]        |    [~]     |
| Beam loading                       |   [~]    |   [o]   |      [o]      |  [o]   |       [-]        |    [o]     |
| Wakefields / CSR                   |   [~]    |   [~]   |      [o]      |  [o]   |       [-]        |    [~]     |
| Longitudinal dynamics, energy spread |  [x]   |   [x]   |      [x]      |  [~]   |       [o]        |    [~]     |
| Autophasing                        |   [x]    |   [x]   |      [o]      |  [o]   |       [-]        |    [~]     |
| Adiabatic damping det(R)           |   [x]    |   [x]   |      [~]      |  [o]   |       [o]        |    [~]     |

Legend:
- `[x]` implemented and cross-validated
- `[~]` implemented (code-supported and wired, not independently cross-validated in this stack)
- `[o]` possible, not done
- `[-]` not applicable

## Evidence per cell

Verified by reading the cited files on 2026-06-19.

RF-Track (`backend/rftrackAdapter.py`, `backend/test/rftrack_linac/`):
- TW [x]: `_build_rf_cavity` (lines 405-491) builds `rft.TW_Structure` from a
  single Fourier coefficient (constant-gradient peak-field model). Benchmarked
  0.06% vs elegant at the energy-gain peak (41.468 vs 41.442 MeV;
  `benchmark_vs_elegant.py`, PRIORITIES L1, git-bug 35e31e6).
- SW [~]: `_build_rf_cavity` builds `rft.SW_Structure` for `structure_type='SW'`;
  no SW cross-validation run exists.
- RF field maps [~]: RF-Track supports COMSOL field maps (adapter docstring line
  6; Niels's `Ioniels/RF-TRACK_UH_ThermionicRFgun` gun workflow). Not used in the
  linac model, which is the analytical single-Fourier TW.
- Multi-cell [~]: `n_cells` is derived from length and phase advance or supplied
  explicitly (lines 446-463); the SLAC structure runs as ~87 synchronous cells.
  This is the analytical constant-gradient model, not extracted production
  cell-by-cell geometry (production tau=0.57 geometry is not in our documents;
  PRIORITIES L1 Phase 4).
- Space charge [~]: `set_space_charge` / `_setup_space_charge` (lines 1112-1207)
  wire `SpaceCharge_P2P` or `SpaceCharge_PIC_FreeSpace` with per-element
  `set_sc_nsteps`. Full-line PIC currently blocked by an RF-Track 2.5.5 core-dump
  (git-bug af9d56c / Risk R5).
- Beam loading [~]: advertised by RF-Track (manual), not exercised in any FELsim
  study.
- Wakefields / CSR [~]: advertised (adapter docstring lines 6, 40, 69), not
  exercised.
- Longitudinal / energy spread [x]: phase scan (`linac_standalone.py`,
  `rftrack_linac_phase_scan.csv`) and energy-gain-vs-phase agreement with elegant
  validate longitudinal energy dynamics. Off-crest chirp covered in W5.
- Autophasing [x]: RF-Track's internal `autophase()` sets phid=0 to on-crest
  (`linac_standalone.py` lines 150-161); validated against the elegant +70 deg
  phase-slippage offset at 1 MeV.
- Adiabatic damping det(R) [x]: `benchmark_vs_elegant.py` extracts det(R_x) via
  unit-perturbation tracking (lines 139-177); the elegant reference reports
  det(R_x)=0.034=p_in/p_out (`meeting_slides.txt` line 20).

elegant (`backend/test/elegant_linac/`):
- The independent reference. RFCA+TWLA model (`slac_linac.ele`, `linac_twiss.ele`,
  `phase_scan.py`); on-crest 39.6 MeV, optimal-phase 41.4 MeV at 1 MeV injection;
  adiabatic damping det(R_x)=0.034 verified; full 6x6 transfer matrix (1st and 2nd
  order) extracted (`meeting_slides.txt`). TW [x] and SW [~] (RFCA is the
  lumped-cavity / SW-like reference; TWLA is the TW reference). CSR/wakefields are
  elegant core strengths [~] (available, not run here for the linac). Space charge
  [o] (elegant SC not wired in this stack).

COSY INFINITY (`backend/cosySimulator.py`, `meeting_slides.txt`):
- SW [~]: `RFCAVX_TM` is a Maxwell-consistent TM010 pillbox cavity with DA Bessel
  functions, transit-time factor, and RF defocusing, implemented in COSY with a
  491-line passing test suite. It is not yet wired into FELsim's COSYAdapter, so
  it is implemented upstream but not in the FELsim linac pipeline.
- TW [o]: travelling-wave extension designed (N28) but not coded.
- Space charge [x]: COSY DA-FMM space charge is the most-validated SC engine in
  the stack (P3 emittance-growth sweep vs xsuite frozen and PIC3D; FODO physics
  F1-F7; PRIORITIES items L3, N9).
- Longitudinal / energy spread [x]: 6D DA maps, R56/T566 (W9, I5, I6).
- det(R) [~]: DA transfer maps make det(R) and adiabatic damping a direct map
  property; not separately cross-validated for an RF section here.
- RF field maps / multi-cell / beam loading / wakefields / autophasing [o]:
  reachable through COSY but not implemented for the linac.

xsuite (`backend/xsuiteAdapter.py`, `backend/test/rftrack_linac/compare_xsuite.py`):
- TW [o] / multi-cell [o]: `compare_xsuite.py` models the SLAC linac as a single
  lumped `Cavity` with V_total = E0 * L (RFCA-equivalent), explicitly noting that
  a per-cell `Cavity` chain is unfaithful at 1 MeV injection without per-cell
  autophasing. A genuine multi-cell TW model is not implemented.
- SW [~]: the lumped `Cavity` is the RFCA/SW-flavored model used in the standalone
  comparison.
- Space charge [~]: `xfields` SpaceCharge3D (PIC) and SpaceChargeBiGaussian
  (frozen) wired via split-operator (`xsuiteAdapter.py` lines 11-12, 171-178);
  used in the P3 SC comparison.
- `XsuiteAdapter` treats `RF_CAVITY`, dipole edges, and fringe as drift (lines
  14, 162-164): autophasing [o], field maps [o], det(R) [o], longitudinal [~]
  (tracking carries longitudinal coordinates but no RF model).

FELsim 1st-order (`backend/beamline.py`):
- Transfer-matrix framework with no RF element (`meeting_slides.txt` line 39: no
  `rfcavityLattice` class). All RF features [-]. Longitudinal [o] via R56/T566 map
  terms; det(R) [o] as an algebraic map property, not exercised for RF.

analytical:
- Closed-form energy-gain / chirp expressions (dE/dphi off-crest, W5). Useful for
  TW [~], longitudinal/energy-spread [~], CSR estimates [~] (Derbenev-Saldin, S9),
  and as a sanity check on det(R) [~]. Not a tracking code.

## Gaps and next steps

The two highest-value linac gaps are both partly waiting on input from Niels.

1. xsuite multi-cell TW. The lumped `Cavity` captures integrated energy gain but
   not the TW spatial extension or longitudinal phase slippage along the
   structure. A faithful per-cell or true TW model is partly blocked on the
   production cavity geometry (tau=0.57 cell-by-cell data is not in our documents).
   If Niels can supply the production geometry or a COMSOL/SuperFish field map,
   xsuite can move from lumped-RFCA to a cross-validatable multi-cell model.
   (git-bug 343b42f, 35e31e6.)

2. COSY travelling-wave DA map. COSY's standing-wave `RFCAVX_TM` is implemented
   and tested; the N28 travelling-wave superposition is designed, not coded.
   Coding it would give a high-order DA TW map cross-checkable against RF-Track's
   analytical TW and elegant's TWLA, and would let COSY cover the linac at the
   same fidelity it already covers the transport line and space charge.

Two smaller items: wire COSY `RFCAVX_TM` into the FELsim COSYAdapter so the
standing-wave cavity is usable from the multi-code framework (it is currently
implemented only inside COSY); and exercise RF-Track's beam loading and
wakefield/CSR paths against the elegant TWLA/CSR reference, which would upgrade
those three RF-Track cells from [~] to [x].
