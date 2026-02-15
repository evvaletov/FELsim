# Injector and Linac Parameters — Presentation Summary

Source: Email from Niels Bidault, Sep 26, 2025, cc Mason McMahon.

## Beamline Layout (Gun to FEL)

```
TM010 RF cavity (thermionic cathode)
  → 2 quads
  → Alpha-Magnet (with electrostatic kickers, knife-edge slit)
  → 2 quads
  → LINAC (S-band, 2856 MHz)
  → Transport line (diagnostic chicane)
  → FEL (undulator)
```

## RF Gun

- TM010 mode, 2856 MHz
- Thermionic cathode
- CAD model and field map being prepared from recovered drawings
- Energy distribution at exit: ~0.6–1.1 MeV
- Nominal energy for linac injection: 1 MeV

## Alpha-Magnet

| Parameter | Value | Source |
|-----------|-------|--------|
| Gradient G | 0.103754 T/(m·A) | Design |
| Typical current | 17.5 A | Logbook |
| Slit width | 16 mm | Logbook |

- Provides energy selection (knife-edge slit) and bunch compression
- Mason McMahon modelling it as a series of dipoles for beam dynamics
- Reference: Wiedemann, Particle Accelerator Physics (2015), p. 449

## Linac

- S-band (2856 MHz), accelerates to ~40 MeV
- Niels working on analytical calculations
- Reference: Wiedemann (2015), p. 628

## Action Items (from Sep 26)

- Check COSY Infinity models for linac / alpha-magnet / gun
- Look into General Particle Tracer (GPT) for cavity with space charge
- Mason to provide dipole parameters for alpha-magnet emulation
- Niels: alpha-magnet Python calculations in `branch_niels`

## References

- Wiedemann, *Particle Accelerator Physics* (2015) — pp. 449 (alpha-magnet), 628 (linac)
- Hadmack, Rev. Sci. Instrum. **84**, 063302 (2013) — energy spread measurement
- UH-044 by Mike Hadmack — MkIII dipole modelling (pp. 48–81)
