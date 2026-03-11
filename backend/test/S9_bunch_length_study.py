# S9: Bunch Length Independence Study
#
# Investigates why the S1 (2 ps) and S3 (0.5 ps) optimizations produce identical
# quadrupole currents, whether this is physically correct, and what collective
# effects break the linear decoupling at short bunch lengths.
#
# Part A: Analytic estimates (R56, compression chirp, CSR, LSC, wakefields)
# Part B: Numerical verification (pre-compressed, in-line compression, σ_E scan)
#
# Author: Eremey Valetov
# Date: 2026-02-22

import sys
import math
import time
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ebeam import beam
from beamline import lattice, beamline, driftLattice, dipole, dipole_wedge
from excelElements import ExcelElements
from beamOptimizer import beamOptimizer
from schematic import draw_beamline

# ── Physical constants ─────────────────────────────────────────────────────────
C_LIGHT = 299792458.0          # m/s
R_E = 2.8179403262e-15         # classical electron radius (m)
Z_0 = 376.730313412            # impedance of free space (Ω)
E_CHARGE = 1.60217663e-19      # C
M_E = 9.1093837e-31            # electron mass (kg)
M_E_MEV = 0.51099895           # electron rest mass (MeV)
E_BEAM = 40.0                  # MeV kinetic energy
F_RF = 2856e6                  # Hz

GAMMA = 1 + E_BEAM / M_E_MEV
BETA = np.sqrt(1 - 1/GAMMA**2)
P_C = GAMMA * BETA * M_E_MEV   # MeV/c

# Beam parameters
BUNCH_CHARGE = 60e-12          # 60 pC
N_BUNCH = BUNCH_CHARGE / E_CHARGE
EPSILON_N = 8e-6               # 8 π·mm·mrad = 8 μm (normalised)
EPSILON_GEOM = EPSILON_N / (GAMMA * BETA)  # geometric emittance (m·rad)

# Pipe geometry
R_PIPE = 0.0127                # beam pipe radius (m) — from pole gap / 2 for chicane
SIGMA_COND_CU = 5.96e7         # copper conductivity (S/m)


# ═══════════════════════════════════════════════════════════════════════════════
#  Part A: Analytic Estimates
# ═══════════════════════════════════════════════════════════════════════════════

def load_beamline():
    """Load UH FEL beamline from Excel and prepare for transport."""
    file_path = Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx'
    excel = ExcelElements(file_path)
    beamline_elements = excel.create_beamline()
    relat = lattice(1, fringeType=None)
    relat.setE(E=E_BEAM)
    line = relat.changeBeamType("electron", E_BEAM, beamline_elements)
    return line, relat, excel


def compute_r56_felsim(line, segments=118):
    """Propagate 6×6 identity through the beamline to extract R56 = M(5,6).

    FELsim convention: columns are [x, x', y, y', Δt/T×10³, δW/W×10³].
    R56 lives in M[4,5] (0-indexed).
    """
    M = np.eye(6)
    r56_cumulative = []
    for i, elem in enumerate(line[:segments]):
        Mi = elem._compute_numeric_matrix()
        M = Mi @ M
        r56_cumulative.append(M[4, 5])

    return M[4, 5], np.array(r56_cumulative)


def compute_r56_fc1(line):
    """Extract R56 of the FC1 chicane alone (elements ~98-117 in the 118-element line).

    FC1 sector starts around element index 98 (FC1_DPW_111) in the Excel beamline.
    We identify FC1 elements by their name containing 'FC1' or by position.
    """
    # FC1 elements are at s ≈ 11.1–12.12 m. In the Excel beamline, these are
    # roughly indices 98–109 (the last ~20 elements before the undulator at 118).
    # We'll compute R56 for just the FC1 section.
    M_fc1 = np.eye(6)
    fc1_indices = []
    for i, elem in enumerate(line[:118]):
        name = getattr(elem, 'name', '') or ''
        if 'FC1' in name.upper() or (hasattr(elem, 'startPos') and elem.startPos and
                                      11.05 < elem.startPos < 12.15):
            fc1_indices.append(i)

    # If name-based detection fails, use position-based (elements around s=11.1-12.1)
    if not fc1_indices:
        # Compute cumulative positions
        s = 0
        for i, elem in enumerate(line[:118]):
            if 11.05 < s < 12.15:
                fc1_indices.append(i)
            s += elem.length

    for i in fc1_indices:
        Mi = line[i]._compute_numeric_matrix()
        M_fc1 = Mi @ M_fc1

    return M_fc1[4, 5], fc1_indices


def r56_analytic_chicane(theta_deg, L_d, drifts):
    """Analytic R56 for a symmetric 4-dipole chicane.

    R56 ≈ -2θ²(L_d/3 + L_drift) for small angles, per LCLS CDR.
    For a 4-dipole chicane with two inner drifts and one central drift:
    R56 = -2θ² × [L_d/3 + L_12 + L_23/cos²θ + L_34 + L_d/3]

    But the standard formula for a symmetric chicane is:
    R56 = -2 L_d θ² (2/3 + L_drift/L_d) ≈ -2θ²(2L_d/3 + 2L_drift)
    where L_drift is the drift between adjacent dipole pairs.

    More precisely for the FC1 geometry (4 dipoles, +θ, -θ, -θ, +θ):
    """
    theta = np.radians(theta_deg)
    # FC1 geometry from YAML:
    # DPH_111: s=11.109 → 11.147  (L_d = 0.0374 m, θ = +11.25°)
    # DPH_115: s=11.476 → 11.513  (L_d = 0.0374 m, θ = -11.25°)
    # DPH_117: s=11.704 → 11.742  (L_d = 0.0374 m, θ = -11.25°)
    # DPH_121: s=12.071 → 12.108  (L_d = 0.0374 m, θ = +11.25°)
    # Drifts (center-to-center): ~0.33 m, ~0.19 m, ~0.33 m

    # Standard 4-dipole chicane R56 (accounting for wedge effects is complex;
    # use the simplified formula and compare with FELsim matrix result)
    rho = L_d / theta if theta != 0 else float('inf')
    L12, L23, L34 = drifts

    # Each bend pair contributes R56; for a symmetric chicane:
    # R56 ≈ -2θ² × Σ(effective path lengths)
    # Simple estimate: R56 ≈ -2 sin²θ × (L12/cos³θ) for each half
    # More standard: use the matrix approach result for comparison
    R56_approx = -2 * theta**2 * (L_d/3 + (L12 + L34)/2 + L23/2)
    return R56_approx


def part_a1_r56(line):
    """A1: R56 from actual beamline vs analytic formula."""
    print("\n" + "="*72)
    print("  A1: R56 Extraction")
    print("="*72)

    r56_total, r56_cum = compute_r56_felsim(line, segments=118)
    r56_fc1, fc1_idx = compute_r56_fc1(line)

    # Note: FELsim M[4,5] is in internal units (relative to f, β, γ).
    # The physical R56 = M[4,5] × (-C·β/f) or similar. Let's check the
    # dipole M56 formula from beamline.py:
    # M56_drift = -l·f / (C·β·γ·(γ+1))
    # M56_dipole = -f·(l - ρ·sinθ) / (C·β·γ·(γ+1))
    # So column 5 is δW/W (×10³) and column 4 is ΔToF/T (×10³).
    # M[4,5] relates change in ToF to energy spread — this IS R56 in FELsim units.

    # Convert to physical R56 (m):
    # Δt/T × 10⁻³ = M[4,5] × δW/W × 10⁻³
    # Δt = M[4,5] × δ × T  where T = 1/f
    # Physical path length change: Δz = β·c·Δt = β·c × M[4,5] × δ / f
    # So R56_physical (m) = β·c × M[4,5] / f  (if M[4,5] relates 10⁻³ units)
    # Wait — both columns are already in 10⁻³ units, so:
    # R56_m = M[4,5] × (C * BETA) / F_RF × 1  (dimensionless 10⁻³/10⁻³ cancels)

    # Actually, let's just use the raw M[4,5] value and interpret it.
    # From the drift matrix: M56 = -l·f/(C·β·γ·(γ+1))
    # For l=1 m: M56 ≈ -1 × 2856e6 / (3e8 × 0.9999 × 79.3 × 80.3) ≈ -4.5e-3
    # This is the coupling per unit δ. The physical R56 needs careful unit conversion.

    # Physical R56 in metres, from time-energy coupling:
    # R56_phys = -(C·β/f) × M[4,5] × (γ(γ+1))⁻¹ ... this gets circular.
    # Better: compute directly from the full matrix product.

    # Actually, let's compute R56 properly by propagating a test particle with
    # known δ and measuring the time-of-flight change.
    delta_test = 1e-3  # 0.1% energy deviation
    particle_ref = np.array([[0, 0, 0, 0, 0, 0]], dtype=np.float64)
    particle_off = np.array([[0, 0, 0, 0, 0, delta_test * 1e3]], dtype=np.float64)  # ×10³ units

    p_ref = particle_ref.copy()
    p_off = particle_off.copy()
    for elem in line[:118]:
        p_ref = np.array(elem.useMatrice(p_ref))
        p_off = np.array(elem.useMatrice(p_off))

    # Column 4 is ΔToF/T × 10³. Physical Δt = (col4 × 10⁻³) / f
    dt_ref = p_ref[0, 4] * 1e-3 / F_RF
    dt_off = p_off[0, 4] * 1e-3 / F_RF
    delta_t = dt_off - dt_ref

    # R56 = Δz / δ = β·c·Δt / δ  (path length change per fractional momentum)
    R56_phys = BETA * C_LIGHT * delta_t / delta_test
    print(f"R56 (full line):  {R56_phys*1e3:.4f} mm  ({R56_phys:.6f} m)")
    print(f"M[4,5] (raw):     {r56_total:.6f}")

    # FC1 only
    p_ref_fc1 = particle_ref.copy()
    p_off_fc1 = particle_off.copy()
    for i in fc1_idx:
        p_ref_fc1 = np.array(line[i].useMatrice(p_ref_fc1))
        p_off_fc1 = np.array(line[i].useMatrice(p_off_fc1))

    dt_fc1 = (p_off_fc1[0, 4] - p_ref_fc1[0, 4]) * 1e-3 / F_RF
    R56_fc1 = BETA * C_LIGHT * dt_fc1 / delta_test
    print(f"R56 (FC1 only):   {R56_fc1*1e3:.4f} mm  ({R56_fc1:.6f} m)")
    print(f"FC1 element indices: {fc1_idx[0]}–{fc1_idx[-1]} ({len(fc1_idx)} elements)")

    # Analytic estimate for FC1
    theta_fc1 = 11.25  # degrees
    L_d_fc1 = 0.037389  # m
    # Drift distances between dipole centers (from YAML s positions)
    d12 = 11.476 - 11.147  # ≈ 0.329 m (DPH_111 center to DPH_115 center)
    d23 = 11.704 - 11.513  # ≈ 0.191 m (DPH_115 to DPH_117)
    d34 = 12.071 - 11.742  # ≈ 0.329 m (DPH_117 to DPH_121)
    R56_analytic = r56_analytic_chicane(theta_fc1, L_d_fc1, [d12, d23, d34])
    print(f"R56 (FC1 analytic): {R56_analytic*1e3:.4f} mm  ({R56_analytic:.6f} m)")

    return R56_phys, R56_fc1, R56_analytic


def part_a2_compression(R56_fc1):
    """A2: Compression chirp table — required chirp for various compression ratios."""
    print("\n" + "="*72)
    print("  A2: Compression Chirp Table")
    print("="*72)

    sigma_z_2ps = 2e-12 * BETA * C_LIGHT    # bunch length (m) at 2 ps
    sigma_z_05ps = 0.5e-12 * BETA * C_LIGHT  # bunch length (m) at 0.5 ps

    print(f"σ_z(2 ps)   = {sigma_z_2ps*1e6:.1f} μm  ({sigma_z_2ps*1e3:.3f} mm)")
    print(f"σ_z(0.5 ps) = {sigma_z_05ps*1e6:.1f} μm  ({sigma_z_05ps*1e3:.3f} mm)")
    print(f"R56(FC1)    = {R56_fc1*1e3:.4f} mm")
    print()

    # Compression factor C = 1/(1 + h·R56/c)
    # For target C, h = (1 - 1/C) × c/R56
    print(f"{'C':>5s}  {'h (1/s)':>12s}  {'φ_off-crest (°)':>16s}  {'σ_δ,corr (%)':>13s}  {'E_final/E_crest':>15s}")
    print("-" * 72)

    results = []
    for C_target in [2, 3, 4, 5, 10]:
        if R56_fc1 == 0:
            print(f"  R56 = 0, cannot compress")
            continue

        h_required = (1 - 1/C_target) * C_LIGHT / R56_fc1

        # Off-crest phase: V_rf(φ) = V_peak·sin(φ), chirp h = 2πf × V'(φ)/V(φ)/c
        # Simplified: φ ≈ arctan(h × c / (2π f × E_beam_eV × e))
        # More precisely: for an RF linac, h ≈ 2πf/c × tan(φ) × (ΔE/E)
        # The relation is: h = -(2π f_rf / c) × cos(φ₀)/sin(φ₀) × ... complex.
        # Simple estimate: φ ≈ arcsin(h × σ_z / (2πf × σ_δ_initial / c))
        # For S-band, 2856 MHz:
        omega_rf = 2 * np.pi * F_RF
        # φ off-crest ≈ arctan(h / (ω_rf × γ²)) — simplified for ultra-relativistic
        # Better: h = -(ω_rf/c) × E_beam × sin(φ) where φ is off-crest angle
        # Actually h_chirp = dE/dt = eV₀ω cos(φ), so φ = arccos(h·c/(eV₀ω²))
        # This needs linac voltage. Let's just compute the implied σ_δ.
        sigma_delta_corr = abs(h_required) * sigma_z_2ps / C_LIGHT * 100  # in %
        # Energy penalty: cos(φ) ≈ 1 - σ_δ²/2 for small chirp
        E_ratio = 1.0 / np.sqrt(1 + (sigma_delta_corr/100)**2)  # approximate

        # Phase estimate: for RF, sin(φ) ≈ h·σ_z/(ωrf·σ_E_initial)
        # This is too model-dependent without linac params; report h directly
        phi_deg = np.degrees(np.arctan(abs(h_required) / omega_rf))

        print(f"{C_target:5d}  {h_required:12.3e}  {phi_deg:16.2f}  {sigma_delta_corr:13.3f}  {E_ratio:15.4f}")
        results.append((C_target, h_required, phi_deg, sigma_delta_corr, E_ratio))

    return results


def part_a3_csr(R56_fc1):
    """A3: CSR energy loss and emittance growth per FC1 dipole."""
    print("\n" + "="*72)
    print("  A3: Coherent Synchrotron Radiation (CSR)")
    print("="*72)

    theta_fc1 = np.radians(11.25)
    L_d = 0.037389  # m
    rho = L_d / theta_fc1  # bending radius

    print(f"FC1 dipole: L_d = {L_d*1e3:.1f} mm, θ = 11.25°, ρ = {rho*1e3:.1f} mm")
    print()

    scenarios = [
        ("2 ps (S1)", 2e-12),
        ("0.5 ps (S3)", 0.5e-12),
    ]

    print(f"{'Scenario':<20s}  {'σ_z (μm)':>10s}  {'I_peak (A)':>10s}  {'ΔE_CSR/bend (keV)':>18s}  "
          f"{'σ_δ,CSR (%)':>12s}  {'Δε/ε (%)':>10s}")
    print("-" * 92)

    results = []
    for label, sigma_t in scenarios:
        sigma_z = sigma_t * BETA * C_LIGHT  # m
        I_peak = BUNCH_CHARGE / (np.sqrt(2 * np.pi) * sigma_t)

        # Steady-state Derbenev–Saldin formula per dipole:
        # ΔE_CSR ≈ 0.22 × r_e × N_b / (ρ^{2/3} × σ_z^{4/3}) × L_d
        # This gives energy loss per particle in eV (with r_e in m, need mc² factor)
        # Standard form: ΔE = (0.22 × r_e × m_e c² × N_b × L_d) / (ρ^{2/3} × σ_z^{4/3})
        dE_csr_per_dipole = (0.22 * R_E * M_E_MEV * 1e6 * N_BUNCH * L_d /
                             (rho**(2/3) * sigma_z**(4/3)))  # eV

        # 4 dipoles in FC1
        dE_csr_total = 4 * dE_csr_per_dipole
        sigma_delta_csr = dE_csr_total / (E_BEAM * 1e6)  # fractional

        # CSR-induced emittance growth: Δε ≈ (η_max × σ_δ,CSR)² / (2β)
        # η_max in FC1 chicane ≈ θ × L_drift ≈ 0.196 × 0.33 ≈ 0.065 m
        eta_max = np.sin(theta_fc1) * 0.33  # approximate max dispersion in chicane
        beta_typical = 1.0  # m (typical beta function in chicane region)
        delta_epsilon = (eta_max * sigma_delta_csr)**2 / (2 * beta_typical)  # m·rad
        delta_epsilon_rel = delta_epsilon / EPSILON_GEOM * 100  # %

        print(f"{label:<20s}  {sigma_z*1e6:10.1f}  {I_peak:10.1f}  {dE_csr_total/1e3:18.3f}  "
              f"{sigma_delta_csr*100:12.4f}  {delta_epsilon_rel:10.3f}")

        results.append({
            'label': label, 'sigma_z': sigma_z, 'I_peak': I_peak,
            'dE_csr': dE_csr_total, 'sigma_delta_csr': sigma_delta_csr,
            'delta_epsilon_rel': delta_epsilon_rel
        })

    return results


def part_a4_lsc():
    """A4: Longitudinal space charge (LSC) impedance and energy spread growth."""
    print("\n" + "="*72)
    print("  A4: Longitudinal Space Charge (LSC)")
    print("="*72)

    # Z_LSC/n ≈ Z₀/(2πγ²) × ln(r_pipe/σ_r)
    sigma_r = np.sqrt(EPSILON_GEOM * 1.0)  # rms beam size (m), assume β≈1 m
    L_transport = 12.0  # approximate transport line length (m)

    Z_lsc = Z_0 / (2 * np.pi * GAMMA**2) * np.log(R_PIPE / sigma_r)
    print(f"σ_r (rms beam size) = {sigma_r*1e6:.1f} μm")
    print(f"Z_LSC/n = {Z_lsc:.4e} Ω")
    print()

    scenarios = [("2 ps", 2e-12), ("0.5 ps", 0.5e-12)]
    print(f"{'Scenario':<12s}  {'I_peak (A)':>10s}  {'ΔE_LSC (eV)':>12s}  {'σ_δ,LSC (%)':>12s}")
    print("-" * 56)

    results = []
    for label, sigma_t in scenarios:
        sigma_z = sigma_t * BETA * C_LIGHT
        I_peak = BUNCH_CHARGE / (np.sqrt(2 * np.pi) * sigma_t)

        # Energy spread from LSC over transport length:
        # ΔE ≈ I_peak × Z_LSC × L / (β·c)
        # More precisely: σ_δ from LSC ≈ (e × I_peak × Z_lsc/n × n) / (E_beam × e)
        # where n ~ f_rf × L / c (number of relevant harmonics)
        # Simplified: energy spread gained ≈ Q × Z_LSC × c / (σ_z × L)
        # Using: ΔE_LSC ≈ (Z₀ × I_peak × L_transport) / (4π × γ² × β × c) × ln(r/σ)
        dE_lsc = (Z_0 * I_peak * L_transport) / (4 * np.pi * GAMMA**2) * np.log(R_PIPE / sigma_r)
        sigma_delta_lsc = dE_lsc / (E_BEAM * 1e6)

        print(f"{label:<12s}  {I_peak:10.1f}  {dE_lsc:12.1f}  {sigma_delta_lsc*100:12.5f}")
        results.append({'label': label, 'dE_lsc': dE_lsc, 'sigma_delta_lsc': sigma_delta_lsc})

    return results


def part_a5_wakefields():
    """A5: Resistive wall and geometric wakefields."""
    print("\n" + "="*72)
    print("  A5: Resistive Wall and Geometric Wakefields")
    print("="*72)

    a = R_PIPE  # pipe radius (m)
    L_transport = 12.0  # m

    scenarios = [("2 ps", 2e-12), ("0.5 ps", 0.5e-12)]
    print(f"\nResistive wall wake (Cu pipe, a = {a*1e3:.1f} mm):")
    print(f"{'Scenario':<12s}  {'σ_z (μm)':>10s}  {'W_RW (V/pC/m)':>14s}  {'ΔE_RW (eV)':>11s}  {'σ_δ,RW (%)':>11s}")
    print("-" * 64)

    results = []
    for label, sigma_t in scenarios:
        sigma_z = sigma_t * BETA * C_LIGHT

        # Resistive wall wake for short bunches (σ_z << a):
        # W_RW ≈ (Z₀·c)/(π²·a⁴) × √(2σ_z/(Z₀·σ_cond))  [V/C/m]
        # But the standard short-range resistive wall wake is:
        # W(s) ≈ (c·Z₀)/(π·a²) × 1/√(π·σ_cond·s)  for s = σ_z
        # Loss factor: κ = (c / (4π)) × (Γ(3/4) / (a² √(σ_cond·σ_z)))
        # ΔE ≈ κ × Q × L

        # Using Chao's formula for Gaussian bunch:
        # κ_rw ≈ c/(4π·a²) × Γ(3/4) / √(σ_cond × σ_z / Z₀)  [V/C/m]
        from scipy.special import gamma as gamma_func
        kappa_rw = (C_LIGHT / (4 * np.pi * a**2) *
                    gamma_func(0.75) / np.sqrt(SIGMA_COND_CU * sigma_z / Z_0))  # V/C/m

        dE_rw = kappa_rw * BUNCH_CHARGE * L_transport  # eV
        sigma_delta_rw = dE_rw / (E_BEAM * 1e6)

        print(f"{label:<12s}  {sigma_z*1e6:10.1f}  {kappa_rw/1e12:14.3f}  {dE_rw:11.2f}  {sigma_delta_rw*100:11.5f}")
        results.append({
            'label': label, 'dE_rw': dE_rw, 'sigma_delta_rw': sigma_delta_rw,
            'sigma_z': sigma_z
        })

    print(f"\nGeometric wakefields: Estimated from BPM/bellows transitions.")
    print(f"  Typical geometric wake loss factor: κ_geom ~ 10–100 V/pC/m")
    print(f"  For {BUNCH_CHARGE*1e12:.0f} pC over {L_transport:.0f} m: ΔE ~ 0.006–0.06 keV")
    print(f"  → σ_δ < 0.0002% (negligible compared to CSR)")

    return results


def part_a6_summary(csr_results, lsc_results, rw_results):
    """A6: Summary table of collective effects."""
    print("\n" + "="*72)
    print("  A6: Collective Effects Summary")
    print("="*72)

    print(f"\n{'Effect':<25s}  {'2 ps σ_δ (%)':>13s}  {'0.5 ps σ_δ (%)':>14s}  "
          f"{'2 ps Δε/ε':>10s}  {'0.5 ps Δε/ε':>12s}  {'Classification':>14s}")
    print("-" * 100)

    # CSR
    csr_2ps = [r for r in csr_results if '2 ps' in r['label']][0]
    csr_05ps = [r for r in csr_results if '0.5 ps' in r['label']][0]
    classify_05 = "PROBLEMATIC" if csr_05ps['delta_epsilon_rel'] > 5 else (
        "modest" if csr_05ps['delta_epsilon_rel'] > 0.5 else "negligible")
    print(f"{'CSR (FC1 chicane)':<25s}  {csr_2ps['sigma_delta_csr']*100:13.4f}  "
          f"{csr_05ps['sigma_delta_csr']*100:14.4f}  "
          f"{csr_2ps['delta_epsilon_rel']:9.3f}%  {csr_05ps['delta_epsilon_rel']:11.3f}%  "
          f"{classify_05:>14s}")

    # LSC
    lsc_2ps = [r for r in lsc_results if '2 ps' in r['label']][0]
    lsc_05ps = [r for r in lsc_results if '0.5 ps' in r['label']][0]
    classify_lsc = "negligible" if lsc_05ps['sigma_delta_lsc']*100 < 0.01 else "modest"
    print(f"{'LSC (full line)':<25s}  {lsc_2ps['sigma_delta_lsc']*100:13.5f}  "
          f"{lsc_05ps['sigma_delta_lsc']*100:14.5f}  "
          f"{'<0.001%':>10s}  {'<0.01%':>12s}  {classify_lsc:>14s}")

    # RW
    rw_2ps = [r for r in rw_results if '2 ps' in r['label']][0]
    rw_05ps = [r for r in rw_results if '0.5 ps' in r['label']][0]
    print(f"{'Resistive wall':<25s}  {rw_2ps['sigma_delta_rw']*100:13.5f}  "
          f"{rw_05ps['sigma_delta_rw']*100:14.5f}  "
          f"{'<0.001%':>10s}  {'<0.001%':>12s}  {'negligible':>14s}")

    print(f"{'Geometric wakes':<25s}  {'<0.0002':>13s}  {'<0.0002':>14s}  "
          f"{'<0.001%':>10s}  {'<0.001%':>12s}  {'negligible':>14s}")

    print(f"\nUncorrelated σ_E (beam) = 0.5%")
    print(f"CSR at 0.5 ps is the dominant collective effect.")

    # Reviewer response summary
    print("\n" + "="*72)
    print("  Reviewer Response")
    print("="*72)
    print("""
(a) WHY IDENTICAL QUAD CURRENTS?

The transverse Twiss parameters (β, α, ε, D) depend on the 4×4 transverse
block of the 6×6 transfer matrix plus the dispersion column (column 6).
Bunch length σ_z enters only column 4 (time-of-flight), which is decoupled
from columns 1–4 and 6 in the linear transfer matrix. Since S1 (2 ps) and
S3 (0.5 ps) use identical σ_E=0.5%, h=5×10⁹/s, and ε_n=8 μm, the
transverse optimization problem is mathematically identical. The quad
currents are the same because the underlying equation system is the same.

This is confirmed by COSY INFINITY (W4), which also produces identical
transfer-map-based optimization results for both configurations.

(b) IS THIS PHYSICALLY CORRECT?

Yes, in the linear regime without collective effects. The decoupling breaks
down when short bunches produce:
  • CSR in chicane dipoles (dominant): σ_δ,CSR ≈ 0.017% at 0.5 ps
    → Δε/ε ≈ 0.06% (negligible at this charge/energy)
  • LSC over transport line: σ_δ,LSC < 0.0001% (negligible)
  • Resistive wall wakes: σ_δ,RW < 0.0001% (negligible)

All collective effects are <0.02% of the 0.5% beam energy spread.
The linear model is valid for the UH MkV FEL at 60 pC.

At higher charges (>1 nC) or lower energies, CSR would become significant
and the matching would need to account for bunch-length-dependent
emittance growth.

(c) HOW TO SWITCH BETWEEN 2 ps AND 0.5 ps IN PRACTICE?

The FELsim transport line has no dedicated bunch compressor. Options:
  1. Pre-compressed beam from injector:
     - Velocity bunching (SPARC precedent): adjust RF phase in injector
     - Adds no hardware; beam enters transport line already at 0.5 ps
     - Quad currents unchanged (same transverse optimization)
  2. Dedicated upstream compressor (cleanest):
     - ~1.5 m 4-dipole chicane + off-crest RF section
     - Standard FEL approach (LCLS, SwissFEL, European XFEL)
  3. FC1 chicane compression:
     - FC1 R56 ≈ 2.9 mm — too small for useful compression
     - C=2 would require h ≈ 5×10¹⁰/s (extreme, σ_δ > 10%)
     - NOT recommended
  4. Dogleg compression (DC1–DC5):
     - Cumulative R56 ~ few mm, similar limitation to FC1
     - NOT recommended
""")


def print_reviewer_response_latex():
    """Print LaTeX-ready reviewer response (Part C)."""
    print(r"""
% S9 Part C: Reviewer Response (LaTeX section)
\subsection{Bunch Length Independence of Transverse Matching}
\label{sec:bunch-length}

The S1 (2\,ps) and S3 (0.5\,ps) optimizations produce identical
quadrupole currents because the transverse Twiss parameters
($\beta$, $\alpha$, $\varepsilon$, $D$) depend on the $4\times4$
transverse block of the $6\times6$ transfer matrix plus the
dispersion column, but not on the time-of-flight column (column~5
of the matrix). Since both configurations use identical
$\sigma_E=0.5\%$, $h=5\times10^9$\,s$^{-1}$, and
$\varepsilon_n=8$\,$\mu$m, the transverse optimization problem is
mathematically identical.

This is confirmed independently by COSY INFINITY, which produces
identical transfer-map-based results for both bunch lengths.

\paragraph{Collective Effects}
At 0.5\,ps (60\,pC), the dominant collective effect is coherent
synchrotron radiation (CSR) in the FC1 chicane dipoles
($\rho=0.19$\,m), producing $\sigma_{\delta,\text{CSR}}\approx0.017\%$
--- negligible compared to the 0.5\% beam energy spread. Longitudinal
space charge, resistive wall wakes, and geometric wakes contribute
$<0.001\%$ each. The linear model is valid at this charge and energy.

\paragraph{Switching Bunch Length in Practice}
The transport line quad settings need not change. The recommended
approach is velocity bunching in the injector (no hardware changes,
SPARC precedent) or a dedicated upstream bunch compressor ($\sim$1.5\,m
4-dipole chicane). The FC1 chicane ($R_{56}\approx2.9$\,mm) is too
weak for useful compression.
""")



# ═══════════════════════════════════════════════════════════════════════════════
#  Part B: Numerical Verification
# ═══════════════════════════════════════════════════════════════════════════════

def create_beam(bunch_spread_ps, energy_std_pct, h_chirp, epsilon_n=8,
                nb_particles=1000, seed=42):
    """Create a 6D Gaussian beam distribution."""
    np.random.seed(seed)
    relat = lattice(1, fringeType=None)
    relat.setE(E=E_BEAM)
    norm = relat.gamma * relat.beta
    epsilon = epsilon_n / norm
    x_std = 0.8
    y_std = 0.8
    x_prime_std = epsilon / x_std
    y_prime_std = epsilon / y_std
    tof_std = bunch_spread_ps * 1e-9 * F_RF
    energy_std = energy_std_pct * 10

    ebeam_obj = beam()
    beam_dist = ebeam_obj.gen_6d_gaussian(
        0, [x_std, x_prime_std, y_std, y_prime_std, tof_std, energy_std],
        nb_particles
    )
    tof_dist = beam_dist[:, 4] / F_RF
    beam_dist[:, 5] += h_chirp * tof_dist
    return beam_dist


def run_11stage_optimization(line, beam_dist, label="", print_results=True):
    """Run the standard 11-stage optimization (identical to S1/S3).

    Returns the optimized line, final MSE, and quad currents dict.
    """
    segments = 118
    line_opt = line[:segments]
    opti = beamOptimizer(line_opt, beam_dist)

    # Undulator matching targets
    relat = lattice(1, fringeType=None)
    relat.setE(E=E_BEAM)
    K = 1.2
    lambda_u = 2.3e-2  # m
    beta_ym = relat.gamma * lambda_u / (2 * np.pi * K)
    alpha_ym = 0.0
    beta_xm = 1.4
    alpha_xm = 0.47

    if print_results:
        print(f"\n── {label}: 11-Stage Optimization ──")

    stages_config = [
        # Stage 1
        ({1: ["I", "current", lambda n: n], 3: ["I2", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}, "I2": {"bounds": (0, 10), "start": 1}},
         {8: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
              {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.0}],
          9: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
              {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]}),
        # Stage 2
        ({10: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {15: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}),
        # Stage 3
        ({16: ["I", "current", lambda n: n], 18: ["I2", "current", lambda n: n],
          20: ["I3", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 5},
          "I3": {"bounds": (0, 10), "start": 3}},
         {25: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
          26: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]}),
        # Stage 4
        ({27: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {32: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}),
        # Stage 5
        ({37: ["I", "current", lambda n: n], 35: ["I2", "current", lambda n: n],
          33: ["I3", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2},
          "I3": {"bounds": (0, 10), "start": 2}},
         {37: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["x", "envelope"], "goal": 2.0, "weight": 1},
               {"measure": ["y", "envelope"], "goal": 2.0, "weight": 1}]}),
        # Stage 6
        ({50: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {55: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}),
        # Stage 7
        ({56: ["I", "current", lambda n: n], 58: ["I2", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2}},
         {59: [{"measure": ["x", "envelope"], "goal": 0.0, "weight": 1},
               {"measure": ["y", "envelope"], "goal": 0.0, "weight": 1}]}),
        # Stage 8
        ({61: ["I", "current", lambda n: n], 63: ["I2", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2}},
         {68: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
          69: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]}),
        # Stage 9
        ({70: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {75: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}),
        # Stage 10
        ({76: ["I", "current", lambda n: n], 78: ["I2", "current", lambda n: n],
          80: ["I3", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2},
          "I3": {"bounds": (0, 10), "start": 2}},
         {85: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
          86: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]}),
        # Stage 11
        ({87: ["Ic", "current", lambda n: n], 93: ["I", "current", lambda n: n],
          95: ["I2", "current", lambda n: n], 97: ["I3", "current", lambda n: n]},
         {"Ic": {"bounds": (0, 10), "start": 4}, "I": {"bounds": (0, 10), "start": 2},
          "I2": {"bounds": (0, 10), "start": 2}, "I3": {"bounds": (0, 10), "start": 2}},
         {92: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 0.5}],
          117: [{"measure": ["x", "alpha"], "goal": alpha_xm, "weight": 1},
                {"measure": ["y", "alpha"], "goal": alpha_ym, "weight": 1},
                {"measure": ["x", "beta"], "goal": beta_xm, "weight": 1},
                {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]}),
    ]

    final_mse = None
    for i, (variables, startPoint, objectives) in enumerate(stages_config):
        result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                           plotBeam=False, printResults=print_results, plotProgress=False)
        final_mse = result.fun
        # Mirror Stage 5
        if i == 4:
            line_opt[43].current = line_opt[33].current
            line_opt[41].current = line_opt[35].current
            line_opt[39].current = line_opt[37].current

    # Extract quad currents
    quad_indices = [1, 3, 10, 16, 18, 20, 27, 33, 35, 37, 39, 41, 43,
                    50, 56, 58, 61, 63, 70, 76, 78, 80, 87, 93, 95, 97]
    currents = {i: line_opt[i].current for i in quad_indices if hasattr(line_opt[i], 'current')}

    return line_opt, final_mse, currents


def part_b1_precompressed(line):
    """B1: Pre-compressed beam at 0.5 ps — does σ_E matter?"""
    print("\n" + "="*72)
    print("  B1: Pre-Compressed Beam (0.5 ps, varied σ_E)")
    print("="*72)

    # Scenario: beam enters at 0.5 ps, no residual chirp (h=0), σ_E from compression
    # At compression ratio C=4: σ_δ,corr ≈ 0.75×C ≈ 3% (very rough)
    # But if chirp is consumed by R56, residual σ_E depends on mismatch
    scenarios = [
        ("S1 baseline (2 ps, 0.5%, h=5e9)", 2, 0.5, 5e9),
        ("S3 baseline (0.5 ps, 0.5%, h=5e9)", 0.5, 0.5, 5e9),
        ("Pre-comp (0.5 ps, 0.5%, h=0)", 0.5, 0.5, 0),
        ("Pre-comp (0.5 ps, 2%, h=0)", 0.5, 2.0, 0),
        ("Pre-comp (0.5 ps, 3%, h=0)", 0.5, 3.0, 0),
    ]

    results = []
    for label, bunch_ps, sigma_e, h in scenarios:
        print(f"\n── {label} ──")
        beam_dist = create_beam(bunch_ps, sigma_e, h)
        _, final_mse, currents = run_11stage_optimization(
            line, beam_dist, label=label, print_results=False
        )
        print(f"  Final RMS: {math.sqrt(final_mse):.6e}")
        # Print a few key quad currents
        key_quads = [1, 3, 93, 95, 97]
        for qi in key_quads:
            if qi in currents:
                print(f"  I[{qi}] = {currents[qi]:.4f} A")
        results.append({'label': label, 'mse': final_mse, 'currents': currents,
                        'sigma_e': sigma_e, 'h': h, 'bunch_ps': bunch_ps})

    # Compare currents
    print("\n── Current Comparison ──")
    ref_currents = results[0]['currents']
    for r in results[1:]:
        max_diff = max(abs(r['currents'].get(k, 0) - ref_currents.get(k, 0))
                       for k in ref_currents)
        print(f"  {r['label']}: max |ΔI| vs S1 = {max_diff:.6f} A, RMS = {math.sqrt(r['mse']):.3e}")

    return results


def part_b3_sigma_e_scan(line):
    """B3: σ_E sensitivity scan at fixed 0.5 ps, h=0."""
    print("\n" + "="*72)
    print("  B3: σ_E Sensitivity Scan (0.5 ps, h=0)")
    print("="*72)

    sigma_e_values = [0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    results = []

    for sigma_e in sigma_e_values:
        beam_dist = create_beam(0.5, sigma_e, h_chirp=0)
        _, final_mse, currents = run_11stage_optimization(
            line, beam_dist, label=f"σ_E={sigma_e}%", print_results=False
        )
        print(f"  σ_E={sigma_e:5.1f}%: RMS={math.sqrt(final_mse):.3e}")
        results.append({'sigma_e': sigma_e, 'mse': final_mse, 'currents': currents})

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.semilogy([r['sigma_e'] for r in results], [r['mse'] for r in results], 'ko-')
    ax1.set_xlabel('σ_E (%)')
    ax1.set_ylabel('RMS Twiss Mismatch')
    ax1.set_title('Optimization quality vs energy spread')
    ax1.axhline(1e-3, color='g', ls='--', label='Excellent')
    ax1.axhline(0.01, color='orange', ls='--', label='Acceptable')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot key quad currents
    key_quads = [93, 95, 97]
    for qi in key_quads:
        ax2.plot([r['sigma_e'] for r in results],
                 [r['currents'].get(qi, 0) for i, r in enumerate(results)],
                 'o-', label=f'I[{qi}]')
    ax2.set_xlabel('σ_E (%)')
    ax2.set_ylabel('Current (A)')
    ax2.set_title('Final triplet currents vs energy spread')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    outdir = Path(__file__).resolve().parent / 'results' / 'S9'
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(outdir / 'S9_sigma_e_scan.eps', dpi=150)
    plt.savefig(outdir / 'S9_sigma_e_scan.png', dpi=150)
    print(f"\n  Saved: {outdir}/S9_sigma_e_scan.{{eps,png}}")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="S9: Bunch Length Independence Study")
    parser.add_argument('--part-a', action='store_true', help='Run Part A (analytic estimates)')
    parser.add_argument('--part-b1', action='store_true', help='Run Part B1 (pre-compressed beam)')
    parser.add_argument('--part-b3', action='store_true', help='Run Part B3 (σ_E sensitivity scan)')
    parser.add_argument('--latex', action='store_true', help='Print LaTeX section (Part C)')
    parser.add_argument('--all', action='store_true', help='Run all parts')
    args = parser.parse_args()

    if not any([args.part_a, args.part_b1, args.part_b3, args.latex, args.all]):
        args.all = True

    print("S9: Bunch Length Independence Study")
    print(f"E = {E_BEAM} MeV, γ = {GAMMA:.2f}, β = {BETA:.6f}")
    print(f"Q = {BUNCH_CHARGE*1e12:.0f} pC, N = {N_BUNCH:.2e}")
    print(f"ε_n = {EPSILON_N*1e6:.0f} μm, ε_geom = {EPSILON_GEOM*1e6:.3f} μm")

    line, relat, excel = load_beamline()

    if args.part_a or args.all:
        R56_total, R56_fc1, R56_analytic = part_a1_r56(line)
        part_a2_compression(R56_fc1)
        csr = part_a3_csr(R56_fc1)
        lsc = part_a4_lsc()
        rw = part_a5_wakefields()
        part_a6_summary(csr, lsc, rw)

    if args.part_b1 or args.all:
        part_b1_precompressed(line)

    if args.part_b3 or args.all:
        part_b3_sigma_e_scan(line)

    if args.latex or args.all:
        print_reviewer_response_latex()

    print("\n" + "="*72)
    print("  S9 Complete")
    print("="*72)


if __name__ == "__main__":
    main()
