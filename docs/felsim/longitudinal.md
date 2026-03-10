# Longitudinal Dynamics

This page summarizes the longitudinal beam dynamics in the UH MkV FEL
transport line: coordinate conventions, momentum compaction, bunch
compression, and key results from the W9 and W12 studies.

## Coordinate Conventions

| Code | coord 5 | coord 6 | Notes |
|------|---------|---------|-------|
| FELsim | $\Delta\text{ToF}/T \times 10^3$ | $\Delta W/W \times 10^3$ | ToF = time of flight, $T = 1/f_\text{RF}$ |
| COSY | $l$ (m) | $\delta = \Delta K/K_0$ | $l$ = path length deviation |
| RF-Track | $ct$ (mm) | $P$ (MeV/c) | Absolute momentum; C7 fix corrects FELsim↔RFT conversion |

**COSY energy convention:** $\delta = \Delta K / K_0$ (kinetic energy
deviation), not $\Delta p / p_0$.  Conversion:
$(l|\delta_p) = (l|\delta_K) \times p_0 c / K_0$.  For 40 MeV electron:
$p_0 c / K_0 = 40.50/40 = 1.012$.

## R56 and T566

**$R_{56}$ (momentum compaction):** the path-length dependence on energy
deviation.  In COSY notation, $R_{56} = (l|\delta) = \text{ME}(5,6)$.

$$
\Delta l = R_{56} \cdot \delta + T_{566} \cdot \delta^2 + \cdots
$$

For the UH FEL transport line (W9 Part A):
- $R_{56} = 27.09$ mm (COSY, $\delta = \Delta K/K_0$)
- $T_{566} = 0$ (no second-order effect)
- $R_{56}$ is **geometry-locked**: determined by dipole angles and spacings,
  not by quadrupole currents (confirmed in W9 Part C)

$R_{56}$ is given by the dipole geometry:

$$
R_{56} = \sum_i \rho_i (\sin\theta_i - \theta_i \cos\theta_i)
$$

where $\rho_i = L_i / \theta_i$ is the bending radius.

## Bunch Compression

The compression factor for a chirped beam ($\delta = h \cdot \Delta t$)
traversing a beamline with momentum compaction $R_{56}$:

$$
C(h) = \frac{1}{|1 + h \cdot R_{56} / (\beta c)|}
$$

Full compression at $h_\text{opt} = -\beta c / R_{56}$.

**Compression floor:** even at optimal chirp, the output bunch length is
bounded by the uncorrelated energy spread:

$$
\sigma_{z,\text{min}} \approx R_{56} \cdot \sigma_{\delta,0}
$$

For $\sigma_\delta = 0.5\%$: $\sigma_{z,\text{min}} = 135\ \mu\text{m}
\approx 0.45$ ps.

### Chirp ↔ Off-crest Angle

For S-band ($f_\text{RF} = 2856$ MHz), assuming single linac section
with $eV_0 \approx K_0$:

$$
\varphi = \arcsin\left(\frac{-h}{2\pi f_\text{RF}}\right)
$$

| Chirp $h$ (s$^{-1}$) | Off-crest $\varphi$ | Application |
|-----------------------|--------------------|-------------|
| $-8.3 \times 10^9$ | 27.5° | $C = 4$ target |
| $-1.1 \times 10^{10}$ | 38.1° | Full compression |
| $\pm 1.8 \times 10^{10}$ | $\pm 90°$ | Single-pass limit |

### Energy Spread Effect

Larger $\sigma_\delta$ always elongates the bunch through
$R_{56} \times \sigma_\delta$, regardless of chirp.  This is why
scenarios with increased energy spread (no chirp) show bunch
*elongation*, not compression.

## Summary of Studies

### W9: COSY Longitudinal Study

- Extracted full 6D transfer map from COSY INFINITY
- Confirmed transverse-longitudinal decoupling: $(\delta|x_j) = 0$
- $R_{56} = 27.09$ mm, $T_{566} = 0$
- Adding $R_{56} = 0$ objective has zero effect (geometry-locked)
- Bunch lengthening: 35% for 0.5 ps, 2.5% for 2 ps (at $\sigma_\delta = 0.5\%$)

### W12: Bunch Compression Feasibility

- Chirp sweep: analytical + COSY map propagation (41 chirp values)
- Compression floor $\approx 0.45$ ps prevents reaching 0.5 ps target
- At $C = 4$ chirp: COSY map $\sigma_{z,\text{out}} \approx 0.67$ ps
- Extended quad bounds (15 A) do not affect $\sigma_z$
- **Conclusion:** transport line is not a bunch compressor; compression
  should occur upstream

#### RF-Track Validation (Part B, post-C7)

Six compression scenarios tracked with the C7-fixed RF-Track adapter.
RF-Track and the COSY map **disagree significantly** for chirped beams:

| Scenario | $\sigma_z$ RF-Track (ps) | $\sigma_z$ COSY map (ps) | Transmission |
|----------|--------------------------|--------------------------|--------------|
| B1: baseline ($h = 0$) | 2.36 | 2.04 | 83.5% |
| B3: $C = 4$ chirp | 1.94 | 0.67 | 58.8% |
| B4: $C = 6$ chirp | 1.81 | 0.52 | 52.5% |
| B5: $\sigma_E = 2\%$ | 2.32 | 2.68 | 56.2% |

The discrepancy grows with chirp ($2.9\times$ at $C = 4$) due to:
aperture losses removing extreme-$\delta$ particles that drive
compression, higher-order chromatic aberrations, and model differences
(analytical sector-bend correction vs triangle-rule fringe).

### R3: Combined Report

The R3 report consolidates W9 and W12 into a single document covering
transfer map analysis, 6D bunch propagation, compression feasibility,
and the RF-Track validation findings.

- Report: `reports/2026/Mar/04/R3_longitudinal_report.pdf`
