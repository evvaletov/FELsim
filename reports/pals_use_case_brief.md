# PALS Use Case: University of Hawaii Mark V Free Electron Laser

**Prepared by:** Eremey Valetov, University of Hawaii at Manoa
**Date:** 2026-03-10
**Contact:** evaletov@hawaii.edu

## Facility

The University of Hawaii Mark V Free Electron Laser (UH MkV FEL) is a
40 MeV electron-beam transport line and oscillator FEL under design at
UH Manoa.  The beamline comprises 137 elements including 23 quadrupoles
(in 11 matching stages), four chicane dipole assemblies, corrector magnets,
diagnostic stations, and a planar undulator.

## FELsim and the Lattice Standard

FELsim is a web-based simulation and optimisation platform for the UH FEL
beamline.  It supports three simulation backends through a unified adapter
layer:

- **FELsim** (in-house) -- linear transfer matrices, fast Stage-11
  Nelder-Mead optimisation
- **COSY INFINITY** -- high-order differential-algebraic transfer maps,
  fringe fields via Enge functions
- **RF-Track** (CERN) -- fully relativistic particle tracking with space
  charge, wakefields, and field maps

A multi-code orchestrator (`MultiCodeSimulator`) chains these backends on
contiguous beamline sections, converting between native coordinate systems
transparently.  The lattice description is the single point of truth for all
three codes.

### Current lattice format (v2/v3)

FELsim's lattice format is a flat YAML (or JSON) array of elements with
PALS-aligned type names adopted independently before we became aware of
the PALS effort:

| FELsim `kind` | PALS equivalent | Notes |
|---------------|----------------|-------|
| `Quadrupole`  | `Quadrupole`   | polarity: focusing / defocusing |
| `SBend`       | `SBend`        | sector dipoles in chicanes |
| `Drift`       | `Drift`        | inter-element spacing |
| `Wiggler`     | `Wiggler`      | undulator section |
| `Kicker`      | `Kicker`       | corrector magnets |
| `Instrument`  | `Instrument`   | BPMs, screens, current monitors |

Version 3 (2026-03) adds optional PALS parameter groups:

- **`MagneticMultipoleP.Bn1`** for quadrupole strength (pole-tip field in
  Tesla), with precedence over the legacy `current_a` parameter.
- **`BendP`** with `g_ref`, `e1`, `e2` for dipole geometry (SI/radians),
  with precedence over legacy degree-based parameters.

A JSON Schema (`lattice_schema_v3.json`) validates all three format versions.

## Value as a PALS Use Case

1. **Multi-code interoperability.** The same lattice file drives three
   independent simulation codes with different coordinate conventions,
   element representations, and physics models.  This is precisely the
   use case PALS targets.

2. **Independent convergence.** FELsim's type vocabulary (`Quadrupole`,
   `SBend`, `Wiggler`, etc.) was chosen before learning of PALS, confirming
   that the naming conventions are natural for the accelerator community.

3. **Real optimisation workload.** The lattice is actively used for
   23-quadrupole matching optimisation across 100+ parameter-space points,
   with results cross-validated between all three codes.

4. **Edge-element pattern (feedback).** FELsim's internal representation
   uses `DIPOLE_WEDGE` thin elements for dipole edge kicks -- a legacy
   artefact from the original beamline spreadsheet.  In the interchange
   format, edge angles are correctly folded into the parent `SBend` as
   `BendP.e1`/`BendP.e2`.  However, we note that some legacy lattice files
   in the community do represent edges as separate elements; PALS could
   consider documenting a recommended migration path for such files.

5. **Strength representation (feedback).** Our quadrupoles are driven by
   current (Amperes), not by normalised gradient K1.  The `Bn1` intermediate
   representation is a practical bridge, but the standard could benefit
   from documenting the `current -> Bn1 -> K1` conversion chain explicitly,
   including the sign conventions for electron beams (negative charge).

## Lattice Summary

| Property | Value |
|----------|-------|
| Beam energy | 40 MeV (electrons) |
| Total elements | 137 (incl. implicit drifts) |
| Quadrupoles | 23 (11 matching stages) |
| Dipole assemblies | 4 chicanes (DPW-DPH-DPW sandwiches) |
| Correctors | Horizontal + vertical kickers |
| Undulator | Planar wiggler |
| Diagnostics | BPMs, OTR screens, current transformers |
| Beamline length | ~15 m |
| Lattice format | YAML/JSON, v1-v3 compatible |

## Next Steps

- Share this brief with the PALS working group via GitHub discussion
- Provide a representative v3 lattice file as an example
- Participate in weekly PALS meetings to discuss feedback items
