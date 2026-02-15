# W5: Chirp Value Assessment for UH MkV FEL
#
# Assesses whether h = 5e9 /s is high, low, or moderate for this beamline.
# Produces phase-space plots and comparison with h = 20e9 /s.
#
# Author: Eremey Valetov
# Date: 2026-02-12

import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ebeam import beam
from beamline import lattice

# ── Constants ─────────────────────────────────────────────────────────────────
ENERGY = 40             # MeV
RF_FREQ = 2856e6        # Hz
REST_MASS = 0.511       # MeV/c²

relat = lattice(1, fringeType=None)
relat.setE(E=ENERGY)
GAMMA = relat.gamma
BETA = relat.beta

# Baseline beam parameters (0.5 ps bunch)
BUNCH_SPREAD = 0.5      # ps
ENERGY_STD_PCT = 0.5    # %
EPSILON_N = 8           # π·mm·mrad
X_STD = 0.8             # mm
NB_PARTICLES = 5000


def generate_beam(h=0, seed=42):
    """Generate 6D Gaussian beam with optional chirp, return beam array."""
    norm = GAMMA * BETA
    epsilon = EPSILON_N / norm
    x_prime_std = epsilon / X_STD
    tof_std = BUNCH_SPREAD * 1e-9 * RF_FREQ
    energy_std = ENERGY_STD_PCT * 10  # δW/W × 10³

    np.random.seed(seed)
    ebeam_gen = beam()
    dist = ebeam_gen.gen_6d_gaussian(
        0, [X_STD, x_prime_std, X_STD, x_prime_std, tof_std, energy_std],
        NB_PARTICLES)

    # Apply chirp
    if h != 0:
        tof_seconds = dist[:, 4] / RF_FREQ  # column 4 → 10³ × Δt/T, ÷ f → 10³ × Δt(s)
        dist[:, 5] += h * tof_seconds  # adds to δW/W × 10³

    return dist


def beam_to_physical(dist):
    """Convert FELsim coordinates to physical units for plotting.

    Returns (Δt in ps, ΔE/E in %).
    """
    # Column 4: Δt/T_RF × 10³ → Δt(s) = col4 × 10⁻³ / f → ps
    dt_ps = dist[:, 4] * 1e-3 / RF_FREQ * 1e12
    # Column 5: δW/W × 10³ → %
    dE_pct = dist[:, 5] * 1e-3 * 100
    return dt_ps, dE_pct


def plot_phase_space_comparison(outdir):
    """Plot energy-time phase space at h=0, 5e9, 20e9 side by side."""
    chirp_values = [0, 5e9, 20e9]
    labels = [r'$h = 0$', r'$h = 5\times10^9$ /s', r'$h = 20\times10^9$ /s']

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)

    for ax, h, label in zip(axes, chirp_values, labels):
        dist = generate_beam(h=h)
        dt_ps, dE_pct = beam_to_physical(dist)

        ax.scatter(dt_ps, dE_pct, s=0.5, alpha=0.3, c='steelblue', rasterized=True)
        ax.set_xlabel(r'$\Delta t$ (ps)')
        ax.set_title(label, fontsize=12)
        ax.set_xlim(-2.0, 2.0)
        ax.set_ylim(-4.0, 4.0)
        ax.axhline(0, color='grey', lw=0.5, ls='--')
        ax.axvline(0, color='grey', lw=0.5, ls='--')

        # Show RMS energy spread
        sigma_E = np.std(dE_pct)
        ax.text(0.05, 0.95, rf'$\sigma_{{E}}$ = {sigma_E:.2f}%',
                transform=ax.transAxes, va='top', fontsize=9,
                bbox=dict(boxstyle='round', fc='white', alpha=0.8))

        # Show correlation slope
        if h != 0:
            dt_range = np.array([-1.5, 1.5])
            # slope in %/ps: h × 1e-12 × 100 = h × 1e-10
            slope_pct_per_ps = h * 1e-12 * 100
            ax.plot(dt_range, slope_pct_per_ps * dt_range,
                    'r-', lw=1.5, label=f'slope = {slope_pct_per_ps:.2f} %/ps')
            ax.legend(loc='lower right', fontsize=8)

    axes[0].set_ylabel(r'$\Delta E / E$ (%)')

    fig.suptitle(r'Energy–Time Phase Space at 0.5 ps Bunch Length ($\varepsilon_n$ = 8 $\pi\cdot$mm$\cdot$mrad)',
                 fontsize=12, y=1.02)
    fig.tight_layout()

    path = outdir / 'chirp_phase_space_comparison.eps'
    fig.savefig(path, bbox_inches='tight', dpi=150)
    print(f'  Saved {path}')
    plt.close(fig)


def plot_energy_spread_vs_chirp(outdir):
    """Plot total RMS energy spread as a function of chirp."""
    h_values = np.linspace(0, 40e9, 50)
    sigma_E_total = []

    for h in h_values:
        dist = generate_beam(h=h)
        _, dE_pct = beam_to_physical(dist)
        sigma_E_total.append(np.std(dE_pct))

    # Analytic prediction: σ_E_total = sqrt(σ_E0² + (h × σ_t)²)
    sigma_t_s = BUNCH_SPREAD * 1e-12
    sigma_E0_pct = ENERGY_STD_PCT
    h_analytic = np.linspace(0, 40e9, 200)
    sigma_analytic = np.sqrt(sigma_E0_pct**2 + (h_analytic * sigma_t_s * 100)**2)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(np.array(h_values) / 1e9, sigma_E_total, 'o', ms=3, label='Simulation')
    ax.plot(h_analytic / 1e9, sigma_analytic, 'r-', lw=1.5,
            label=r'$\sqrt{\sigma_{E,0}^2 + (h\,\sigma_t)^2}$')

    # Mark key values
    for h_mark, color, ls in [(5, 'green', '--'), (20, 'orange', '--')]:
        sigma_mark = np.sqrt(sigma_E0_pct**2 + (h_mark * 1e9 * sigma_t_s * 100)**2)
        ax.axvline(h_mark, color=color, ls=ls, lw=1, alpha=0.7)
        ax.annotate(f'h = {h_mark}e9\n({sigma_mark:.2f}%)',
                    xy=(h_mark, sigma_mark), fontsize=8,
                    xytext=(h_mark + 2, sigma_mark + 0.2),
                    arrowprops=dict(arrowstyle='->', color=color),
                    color=color)

    ax.set_xlabel(r'Chirp $h$ ($10^9$ /s)')
    ax.set_ylabel(r'Total RMS $\sigma_E$ (%)')
    ax.set_title('Energy Spread Growth with Chirp (0.5 ps bunch)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    path = outdir / 'chirp_energy_spread_growth.eps'
    fig.savefig(path, bbox_inches='tight', dpi=150)
    print(f'  Saved {path}')
    plt.close(fig)


def plot_chirp_in_context(outdir):
    """Plot energy-time tilt per σ_t for context."""
    fig, ax = plt.subplots(figsize=(7, 4.5))

    # ΔE per σ_t at different chirp values
    sigma_t_s = BUNCH_SPREAD * 1e-12  # 0.5 ps
    h_values = np.linspace(0, 40e9, 200)
    dE_per_sigma = h_values * sigma_t_s * ENERGY  # MeV per σ_t

    ax.plot(h_values / 1e9, dE_per_sigma, 'b-', lw=2)

    # Mark key values
    for h_mark, label, color in [
        (5, 'UH baseline (5e9)', 'green'),
        (20, 'Emittance-conserved (20e9)', 'orange'),
        (1.36, 'SDL/BNL residual (1.36e9)', 'purple'),
    ]:
        dE = h_mark * 1e9 * sigma_t_s * ENERGY
        ax.plot(h_mark, dE, 'o', color=color, ms=8, zorder=5)
        ax.annotate(f'{label}\n{dE:.3f} MeV/σ_t',
                    xy=(h_mark, dE), fontsize=8,
                    xytext=(h_mark + 3, dE + 0.02),
                    arrowprops=dict(arrowstyle='->', color=color),
                    color=color)

    # Off-crest angle axis (h ≈ ω_RF × sin(φ) for single-section approximation)
    omega_RF = 2 * np.pi * RF_FREQ
    valid = h_values <= omega_RF
    ax2 = ax.twiny()
    if np.any(valid):
        phi_max = np.degrees(np.arcsin(min(h_values[valid].max() / omega_RF, 1.0)))
        ax2.set_xlim(0, phi_max)
    ax2.set_xlabel(r'Off-crest angle $\phi$ (degrees, single-section approx.)')

    ax.set_xlabel(r'Chirp $h$ ($10^9$ /s)')
    ax.set_ylabel(r'$\Delta E$ per $\sigma_t$ (MeV)')
    ax.set_title(r'Energy Tilt at 0.5 ps, $E_0$ = 40 MeV')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    path = outdir / 'chirp_context_comparison.eps'
    fig.savefig(path, bbox_inches='tight', dpi=150)
    print(f'  Saved {path}')
    plt.close(fig)


def off_crest_str(h, omega_RF):
    """Format off-crest angle, or 'N/A' if h > ω_RF (multi-section regime)."""
    ratio = h / omega_RF
    if abs(ratio) <= 1:
        return f"{np.degrees(np.arcsin(ratio)):.1f}°"
    return "N/A (multi-section)"


def print_assessment():
    """Print the chirp assessment summary."""
    sigma_t_s = BUNCH_SPREAD * 1e-12
    omega_RF = 2 * np.pi * RF_FREQ

    print("\n" + "=" * 70)
    print("W5: Chirp Value Assessment for UH MkV FEL")
    print("=" * 70)

    print(f"\nBeamline: UH MkV FEL, E = {ENERGY} MeV, f_RF = {RF_FREQ/1e6:.0f} MHz")
    print(f"Bunch: σ_t = {BUNCH_SPREAD} ps, σ_E = {ENERGY_STD_PCT}%")
    print(f"ω_RF = 2π × {RF_FREQ/1e6:.0f} MHz = {omega_RF:.3e} rad/s")

    print("\n--- Chirp values under assessment ---")
    for h, label in [(5e9, "UH baseline"), (20e9, "Emittance-conserved (S2)")]:
        dE_pct = h * sigma_t_s * 100
        dE_MeV = h * sigma_t_s * ENERGY
        sigma_total = np.sqrt(ENERGY_STD_PCT**2 + dE_pct**2)
        print(f"\n  h = {h:.0e} /s ({label})")
        print(f"    ΔE/E per σ_t  = {dE_pct:.3f}%  ({dE_MeV:.3f} MeV)")
        print(f"    Off-crest      ≈ {off_crest_str(h, omega_RF)}")
        print(f"    Total σ_E      = {sigma_total:.3f}%  (vs {ENERGY_STD_PCT}% uncorrelated)")
        print(f"    Energy spread increase = {(sigma_total/ENERGY_STD_PCT - 1)*100:.1f}%")

    print("\n--- Comparison with other S-band FELs ---")
    comparisons = [
        ("SDL/BNL (70 MeV, residual)", 1.36e9),
        ("SPARC (150 MeV, low)", 6.6e9),
        ("SPARC (150 MeV, high)", 22.7e9),
        ("LCLS (typical post-BC1)", 30e9),
    ]
    for name, h in comparisons:
        print(f"  {name}: h = {h:.2e} /s  ({off_crest_str(h, omega_RF)} at 2856 MHz)")

    print("\n--- Assessment ---")
    print("  h = 5e9 /s is LOW-TO-MODERATE for an S-band FEL at 40 MeV.")
    print("  It corresponds to ~16° off-crest, adding only 0.25%/σ_t to the")
    print("  energy spread (total σ_E grows from 0.50% to 0.56%, a 12% increase).")
    print("  This is comparable to the SDL/BNL residual chirp level.")
    print()
    print("  h = 20e9 /s is MODERATE-TO-HIGH. It adds 1.0%/σ_t, doubling")
    print("  the total energy spread to 1.12%. This is comparable to SPARC")
    print("  high-chirp operation.")
    print()
    print("  PROVENANCE: h = 5e9 /s is NOT from arXiv:2510.14061v1 (paper")
    print("  does not mention chirp). Earliest codebase appearance: goldTwiss.py.")
    print("  Not hardcoded in FELsim core (which defaults to h = 0).")


if __name__ == "__main__":
    outdir = Path(__file__).resolve().parent / 'results' / 'w5_chirp'
    outdir.mkdir(parents=True, exist_ok=True)

    print("W5: Chirp Value Assessment — Generating plots...")

    print("\n[1/3] Phase-space comparison (h = 0, 5e9, 20e9)")
    plot_phase_space_comparison(outdir)

    print("\n[2/3] Energy spread growth vs chirp")
    plot_energy_spread_vs_chirp(outdir)

    print("\n[3/3] Chirp in context (ΔE per σ_t)")
    plot_chirp_in_context(outdir)

    print_assessment()
    print(f"\nPlots saved to: {outdir}/")
