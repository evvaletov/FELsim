# Chicane Dipole Fieldmap

## Overview

The MkIII diagnostic chicane contains four dipoles that bend the beam through a
dogleg path for energy spread measurement. Each dipole's fringe field profile is
represented as a 1D fieldmap used by COSY INFINITY's midplane-based geometric
element (MGE) via the `FR 3` fringe field model.

The fieldmap was derived from an OPERA-3D simulation of the dipole:

1. **OPERA-3D** — 3D magnet simulation produced $B_y(z)$ along the midplane
2. **CSV export** — 132 non-uniformly spaced points over $\pm 10$ mm (`UH_chicane_fringe.csv`)
3. **Mathematica Enge fit** — analytical fit to the CSV data (`UH_chicane_fringe.nb`)
4. **Uniform resampling** — 201 points at 1 mm spacing over $\pm 100$ mm (`chicane_dipole_fieldmap.dat`)

## Fieldmap Parameters

| Parameter | Value |
|-----------|-------|
| Number of points ($N$) | 201 |
| Step size (DELTAS) | 0.001 m (1 mm) |
| Range | $\pm 100$ mm |
| Peak $B_y$ | 0.5307 T |
| FWHM | 47.2 mm |
| Format | COSY 1D MGE (header + values) |

## Corrected Fieldmap Profile

The corrected fieldmap profile with the OPERA-3D source data overlaid in the
$\pm 10$ mm region where measurements exist:

```{figure} images/fieldmap_profile.png
:width: 100%
:alt: Corrected chicane dipole fieldmap profile with OPERA-3D source data overlay

Corrected $B_y(z)$ profile (blue) with OPERA-3D source data (red scatter).
Peak field 0.5307 T matches the source data.
```

## Before/After Correction

The fieldmap was corrected on 2026-02-22. The old version had an erroneous
$\times 0.835$ momentum scaling factor applied to all field values, reducing the
peak from 0.5307 T to 0.4433 T — a 19.7% underestimate:

```{figure} images/fieldmap_correction.png
:width: 100%
:alt: Before and after fieldmap correction comparison

Old (dashed red) vs corrected (solid blue) fieldmap. The shaded region shows
the field that was missing due to the erroneous scaling.
```

## Bug History

Two bugs were present in the original fieldmap file:

1. **DELTAS = 0** — The step size was set to 0.0 instead of 0.001 m. COSY
   interpreted this as zero spacing between field samples, collapsing the
   entire fringe field to a point.

2. **Erroneous $\times 0.835$ scaling** — All field values were multiplied by
   $P/P_{45°} = 0.8352$, a momentum ratio that had no physical basis in this
   context. This reduced the peak field from 0.5307 T to 0.4433 T.

Both bugs were discovered during the W4 cross-validation study
(`backend/test/W4_cosy_xval_report.tex`), which compared COSY chicane tracking
against the first-order FELsim model. The COSY results showed anomalously small
chicane dispersion, prompting inspection of the fieldmap.

These bugs blocked the COSY MGE (`FR 3`) fringe field optimization, which
requires a correct fieldmap as input. With the corrected fieldmap, the MGE
optimization can proceed.

## Source Files

| File | Description |
|------|-------------|
| `fields/chicane_dipole_fieldmap.dat` | Corrected fieldmap (201 points, COSY format) |
| `fields/calculation/UH_chicane_fringe.csv` | OPERA-3D source data (132 points, mm/Gauss) |
| `fields/calculation/UH_chicane_fringe.nb` | Mathematica notebook: Enge fit + resampling |
| `fields/calculation/UH_chicane_dipole.txt` | Dipole specification notes |
