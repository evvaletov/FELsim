# Decision Memo: Does PALS Provide a Beam-Files Standard?

**To:** Niels Bidault
**From:** Eremey Valetov
**Date:** 2026-06-22
**Re:** Particle-distribution interchange format for the UH MkV FEL toolchain

## Answer

No. PALS (Particle Accelerator Lattice Standard) standardizes the lattice and
optics description only. It does not define a beam-files standard for particle
distributions or phase-space dumps. The choice of beam-file format is ours to
make, independently of PALS.

Recommendation: adopt openPMD-beamphysics (HDF5) as the primary particle
interchange format across FELsim, the linac/transport codes, and Genesis4,
with plain ASCII as the human-readable fallback.

## Evidence

### PALS covers lattices, not beams

Every PALS artifact in this repo is lattice-scoped. The submission package
(`reports/pals_submission/`) and the use-case brief
(`reports/pals_use_case_brief.md`) describe element types (`Quadrupole`,
`SBend`, `Drift`, `Wiggler`, `Kicker`, `Instrument`), magnet strengths
(`MagneticMultipoleP.Bn1`), and bend geometry (`BendP` with `g_ref`, `e1`,
`e2`). Nothing in the package mentions particle distributions, phase-space
files, charge, or bunch profiles.

The three roadmap items that engage PALS are all lattice work:

- I2 (`PRIORITIES.md`): register UH MkV FEL as a PALS use case; feedback on
  element-type gaps (DIPOLE_WEDGE, edge angles).
- I3: the `pals2cosy` converter, lattice-to-COSY only.
- I8: FELsim v3 lattice format, PALS-aligned fields.

The PALS reference implementations listed in the submission README are BMAD,
BLAST ImpactX, and BLAST WarpX, joined by lattice exchange. The standard's
scope is the machine description, not the beam in it.

### FELsim has no portable beam-file I/O today

`backend/cosyParticleSimulator.py` writes particles only in COSY-native and
plain-text forms (`write_particle_file`, around lines 941-999):

- `rray` / `sray`: COSY's 8-coordinate native format (X, A, Y, B, T, D, G, Z),
  with G and Z zero-padded.
- `ascii` / `ascii_simple`: 6-column whitespace text, optionally with a count
  header.
- `binary`: `NotImplementedError`.

There is no openPMD, SDDS, ASTRA, or HDF5 writer. `load_particles_from_file`
(line 576) is also a stub (`NotImplementedError`). `backend/requirements.txt`
carries no `h5py`, `openpmd-api`, `pmd_beamphysics`, or SDDS dependency (it
does carry `xsuite`, which bundles HDF5 access transitively, but FELsim does
not use it for beam I/O). So the beam-file question is genuinely open: we are
not locked into anything yet.

### The downstream handoff is the real driver

The open tracker is git-bug `8d0d833` (P3): FELsim end-of-beamline 6D phase
space into Genesis4 for the FEL stage. It explicitly lists ASTRA, SDDS, and
openPMD-beamphysics as candidates and flags the Genesis4 particle-input path
as "TBD, needs investigation."

I investigated it. Current state of `~/UH/GENESIS/UHFEL_undulator/`:

- The beam is specified analytically in every active input deck via the
  `&beam` namelist (`current`, `delgam`, `ex`, `ey`, `betax/y`, `alphax/y`)
  plus `&profile_gauss` for the current profile. No external particle file is
  imported in any committed `.in` file.
- Genesis4 already writes 6D macroparticle dumps as HDF5 (`*_beam.par.h5`:
  gamma, theta, x, y, px, py), consumed by the project's post-processing
  (`oscillator.py`, `plot_phase_space_fig4cd.py`).
- Genesis4 v4 has a native HDF5 particle-import path (`importbeam` /
  `importdistribution`), not yet exercised here.

The implication: Genesis4 is HDF5-native on both ends. A handoff format that
maps cleanly onto Genesis4's HDF5 particle layout minimizes friction. That
points at openPMD-beamphysics, whose `pmd_beamphysics` reader/writer set
already includes Genesis4.

The other end of the toolchain, the elegant linac reference at
`backend/test/elegant_linac/`, is SDDS-native (output files carry the `SDDS1`
magic header). So SDDS is the incumbent on the linac side and cannot be
ignored; openPMD-beamphysics reads it.

## Format comparison

| Format | Producers in our stack | Consumers / readers | Metadata | Verdict |
|---|---|---|---|---|
| openPMD-beamphysics (HDF5, `pmd_beamphysics`) | proposed FELsim writer | Genesis4, elegant, ASTRA, Bmad, Impact-T, GPT (via `pmd_beamphysics`); self-describing for any HDF5 reader | rich, self-describing: species, charge, units, reference momentum, t/z slicing, metadata dict | Primary. Single hub format, native units and charge, Genesis4 path supported, reads the SDDS and ASTRA we already have. |
| elegant SDDS | elegant (incumbent at `backend/test/elegant_linac/`) | elegant, `pmd_beamphysics`, SDDS toolkit | strong, column-typed with units and parameters | Keep as elegant's native I/O; bridge via `pmd_beamphysics`. Not the cross-code hub (ecosystem narrower than openPMD). |
| Genesis4 native dist (HDF5) | Genesis4 (`*_beam.par.h5`) | Genesis4, `pmd_beamphysics` | Genesis4-specific layout (gamma, theta, x, y, px, py) | Endpoint format only; subsumed by reading/writing through openPMD-beamphysics. |
| ASTRA | Niels's injector (BEAMPATH/ASTRA chain) | ASTRA, `pmd_beamphysics`, elegant import | fixed-column text, charge per particle, status flags | Keep as the injector's native output; convert into the hub via `pmd_beamphysics`. Not the hub itself (text, no rich metadata). |
| plain ASCII | FELsim (`write_particle_file` ascii) | anything; trivial to parse | none (bare columns, conventions in code only) | Fallback / debug only. No units, no charge, no provenance. Retain for inspection and quick diffs. |
| COSY rray/sray | FELsim COSY backend (`cosyParticleSimulator.py`) | COSY INFINITY only | 8 coords, COSY units, no charge | Internal to the COSY backend. Not an interchange format. |

## Recommendation

1. openPMD-beamphysics (HDF5) as the primary interchange format for particle
   distributions across FELsim, the linac/injector codes, and Genesis4. It is
   self-describing (species, charge, units, reference momentum), and its
   `pmd_beamphysics` package already reads and writes the formats we touch:
   ASTRA (injector), SDDS/elegant (linac), and Genesis4 (FEL). One hub format,
   converters to the rest for free.

2. Plain ASCII as the fallback for human inspection and quick cross-checks.
   Keep FELsim's existing ASCII writer; do not invest further in it.

3. Keep each code's native format at its own boundary (COSY rray/sray inside
   the COSY backend, SDDS at elegant, ASTRA at the injector) and bridge through
   openPMD-beamphysics rather than forcing every code onto one on-disk format.

This decision is orthogonal to the PALS work. We continue to align the lattice
format with PALS (I2/I3/I8) and separately adopt openPMD-beamphysics for beams.

## Next steps

1. Confirm the Genesis4 particle-input path. Exercise the v4
   `importbeam`/`importdistribution` namelist on a `pmd_beamphysics`-written
   HDF5 file end to end in `~/UH/GENESIS/UHFEL_undulator/` (the input decks
   currently use analytic `&beam` only). This validates the handoff before we
   commit a FELsim writer.

2. Prototype a `pmd_beamphysics` writer in FELsim. Add an export path in
   `cosyParticleSimulator.py` (or a thin sibling) that dumps the
   end-of-beamline 6D distribution as an openPMD-beamphysics `ParticleGroup`,
   adding `pmd_beamphysics` (and its `h5py` dependency) to
   `backend/requirements.txt`. Round-trip test against the existing COSY rray
   output. Closes git-bug `8d0d833`.

3. One line to the openPMD-beamphysics maintainers via PALS. Christopher Mayes
   (ChristopherMayes) maintains `pmd_beamphysics`; Axel Huebl (ax3l) leads the
   openPMD standard and participates in PALS through BLAST. A short note on
   PALS Issue #176 confirms the Genesis4 reader status and surfaces our use
   case to the right people. Low cost, useful signal.

Reviewed evidence: `reports/pals_submission/*`, `reports/pals_use_case_brief.md`,
`backend/cosyParticleSimulator.py`, `backend/requirements.txt`,
`backend/test/PRIORITIES.md`, `backend/test/elegant_linac/`, git-bug `8d0d833`,
`~/UH/GENESIS/UHFEL_undulator/genesis4/*.in`.
