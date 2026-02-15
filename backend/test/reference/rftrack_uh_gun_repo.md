# RF-TRACK UH Thermionic RF Gun — Repository Analysis

**Repo:** https://github.com/Ioniels/RF-TRACK_UH_ThermionicRFgun.git
**Author:** Niels Bidault (GitHub: Ioniels)
**Analyzed:** 2026-02-12
**Status:** Active development (5 commits, Feb 9–11 2026)

## What It Simulates

**Only the thermionic RF gun cavity** — the first element in the UH FEL injector.
Does NOT include alpha magnet, linac, or downstream beamline.

- TM010 mode, quarter-wave (λ/4) geometry at f = 2.856 GHz
- Physical length: λ/4 = 26.24 mm
- Field from COMSOL simulation: on-axis peak |Ez| ~ 3.8×10⁷ V/m
- Accelerating voltage: V_peak = 0.937 MV
- Richardson-Dushman thermionic emission with Schottky barrier lowering

## Key Files

| File | Purpose |
|------|---------|
| `UH_gun_tracking_demo.ipynb` | Main simulation notebook |
| `utils.py` | All physics: emission model, I/Q phasor, tracking wrappers |
| `load_fieldmap_mat.py` | COMSOL .mat field map loading |
| `config.py` | RF-Track import wrapper |
| `field_maps/*.mat` | COMSOL XY/YZ plane field data (~65 MB) |

## Beam Parameters

### Cathode / Emission

| Parameter | Initial Commit | Latest Commit |
|-----------|---------------|--------------|
| Cathode radius | 1.57 mm | 1.57 mm |
| Temperature | 1450 K | 1700 K |
| Work function | 2.1 eV | 2.1 eV |
| Initial pz | 4 keV/c | 4 keV/c |
| Target charge | 50 pC | 200 pC |
| Space charge | off | **enabled** |
| N particles | 10,000 | 100,000 |

### Thermal Emittance at Cathode

| Temperature | Effective R (0.157 mm) | Full R (1.57 mm) |
|-------------|----------------------|-------------------|
| 1450 K | 0.078 π·mm·mrad | 0.78 π·mm·mrad |
| 1700 K | 0.084 π·mm·mrad | 0.84 π·mm·mrad |

### Gun Exit (φ = 100°, no space charge, initial commit)

| Parameter | Value |
|-----------|-------|
| Mean kinetic energy | ~356 keV |
| Mean pz | ~0.70 MeV/c |
| σ_x, σ_y | ~2.85 mm (expanding) |
| Energy spread | ~15% (huge) |
| Temporal spread | ~350 ps (full RF period) |
| Chirp | Strong positive (higher E at tail) |

## Relevance to FELsim

| Parameter | FELsim (linac exit) | Gun Exit | Gap |
|-----------|--------------------:|:---------|-----|
| Energy | 40 MeV | 0.356 MeV | 112× — linac accelerates |
| ε_n | 8 π·mm·mrad | 0.08–0.84 | Emittance growth in gun, α-magnet, linac |
| σ_x, σ_y | 0.8 mm | 2.85 mm | Linac + matching optics refocus |
| σ_t | 2 ps | ~350 ps | α-magnet compresses |
| σ_E/E | 0.5% | ~15% | α-magnet selects + linac reduces |
| Chirp h | 0 | Strong + | α-magnet compensates; linac phase sets final |

**Key takeaways:**
1. Gun does not directly provide linac-exit parameters — major transformation
   chain (α-magnet → linac) intervenes
2. Thermal emittance (~0.84 π·mm·mrad) is 10× smaller than FELsim's 8;
   the difference is emittance growth through the injector chain
3. Positive chirp from gun is compressed by α-magnet's negative R56
4. 200 pC at 2 ps → ~100 A peak current, consistent with FEL requirements
