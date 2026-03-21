# PALS Use Case: University of Hawaii Mark V Free Electron Laser

**Prepared by:** Eremey Valetov, University of Hawaii at Manoa
**Contact:** evaletov@hawaii.edu
**Date:** 2026-03-21

## Facility

The University of Hawaii Mark V FEL is a 40 MeV electron-beam transport line
and oscillator FEL under design at UH Manoa. The beamline has 137 elements:
23 quadrupoles in 11 matching stages, four chicane dipole assemblies, corrector
magnets, diagnostics (BPMs, OTR screens), and a planar undulator.

## Multi-code simulation

FELsim drives three simulation backends from a single lattice file:

| Backend | Method | Role |
|---------|--------|------|
| FELsim (in-house) | Linear transfer matrices | Fast matching optimisation |
| COSY INFINITY | Arbitrary-order DA maps | Fringe fields, aberrations |
| RF-Track (CERN) | Particle tracking | Space charge, wakefields |

A multi-code orchestrator chains backends on contiguous beamline sections,
converting between coordinate systems transparently. The lattice file is
the single source of truth for all three codes.

## PALS alignment

FELsim's lattice format (YAML/JSON) independently adopted PALS-compatible
type names: `Quadrupole`, `SBend`, `Drift`, `Wiggler`, `Kicker`, `Instrument`.
A v3 extension adds optional PALS parameter groups:

- **`MagneticMultipoleP.Bn1`** — quadrupole pole-tip field (Tesla)
- **`BendP`** with `g_ref` (1/m), `e1`, `e2` (radians)

PALS fields coexist with legacy parameters and take precedence when present.
A JSON Schema validates all format versions.

The [pals2cosy](https://github.com/evvaletov/pals2cosy) converter handles both
official PALS and FELsim lattice formats.

## Value as a use case

1. **Multi-code interoperability.** One lattice file drives three codes with
   different coordinate conventions, element models, and physics — the central
   PALS use case.

2. **Independent convergence.** FELsim's type names were chosen before we
   learned of PALS, confirming the conventions are natural.

3. **Production workload.** The lattice is used for 23-quadrupole matching
   optimisation across 100+ parameter points, cross-validated between all
   three codes.

## Feedback

1. **Edge elements.** Our internal model uses separate `DIPOLE_WEDGE` thin
   elements for dipole edge kicks (legacy artefact). In the interchange format
   we fold edge angles into `BendP.e1`/`e2`. Some community lattice files also
   represent edges as separate elements; a documented migration path would help.

2. **Strength representation.** Our quadrupoles are driven by current (Amperes),
   not K1. `Bn1` bridges this, but the standard could document the full
   `current → Bn1 → K1` conversion chain, including sign conventions for
   electron beams.

## Lattice summary

| Property | Value |
|----------|-------|
| Beam energy | 40 MeV electrons |
| Elements | 137 |
| Quadrupoles | 23 (11 matching stages) |
| Dipoles | 4 chicane assemblies |
| Length | ~15 m |
| Format | YAML/JSON, v1–v3 compatible |

## Attached

- `uhfel_excerpt_v3.yaml` — 15-element excerpt with PALS fields
- Full 137-element lattice available on request
