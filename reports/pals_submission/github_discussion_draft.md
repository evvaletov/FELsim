# GitHub Discussion Draft

**Repository:** https://github.com/pals-project/pals
**Category:** Use Cases

---

**Title:** Use Case: UH Mark V FEL — Multi-Code Transport Line Simulation

**Body:**

We would like to register the University of Hawaii Mark V Free Electron Laser
as a PALS use case.

**Facility.** The UH MkV FEL is a 40 MeV electron transport line (~15 m,
137 elements) under design at UH Manoa: 23 quadrupoles, four chicane dipole
assemblies, correctors, diagnostics, and a planar undulator.

**Multi-code architecture.** Our simulation platform (FELsim) uses three
backends from a single lattice description:

- **FELsim** — linear transfer matrices for fast matching optimisation
- **COSY INFINITY** — arbitrary-order DA transfer maps with fringe fields
- **RF-Track** (CERN) — fully relativistic particle tracking with space charge

A multi-code orchestrator chains these on contiguous beamline sections. The
lattice file is the single source of truth for all three codes — precisely
the use case PALS targets.

**PALS alignment.** Our YAML/JSON lattice format independently adopted
PALS-compatible type names (`Quadrupole`, `SBend`, `Drift`, `Wiggler`,
`Kicker`, `Instrument`). A v3 extension adds optional PALS parameter groups:

- `MagneticMultipoleP.Bn1` — quadrupole pole-tip field (Tesla)
- `BendP` with `g_ref` (1/m), `e1`, `e2` (radians)

Both coexist with legacy parameters; PALS fields take precedence. The
[pals2cosy](https://github.com/evvaletov/pals2cosy) converter handles both
official PALS and FELsim formats.

**Example.** An excerpt of the UH FEL beamline in v3 format is attached below.
It shows quadrupoles with `Bn1`, sector bends with `BendP`, and the edge-kick
pattern used for chicane dipoles. The full 137-element lattice is available on
request.

**Feedback.** Two items from our implementation experience:

1. **Edge elements.** We internally use separate `DIPOLE_WEDGE` elements for
   dipole edge kicks (legacy artefact). In the interchange format, edge angles
   fold into `BendP.e1`/`e2`. Some community lattice files also represent edges
   separately — a documented migration path would help.

2. **Strength chain.** Our quadrupoles are driven by current (Amperes), not K1.
   `Bn1` bridges this, but the standard could document the full
   `current → Bn1 → K1` conversion, including sign conventions for electrons.

---

**Attachment:** `uhfel_excerpt_v3.yaml`
