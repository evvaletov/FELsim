"""
RF-Track vs xsuite linac comparison (L3.3).

Phase 1 (this script): deterministic transport -- validate that the xsuite
cavity-chain approximation of the SLAC 3-m S-band TW structure agrees with
the elegant RFCA reference and quantifies the TW-vs-RFCA gap to RF-Track.

  - Single-particle phase scan: K_out vs phase, three codes overlaid
  - Gaussian bunch tracking, no SC: σ_x(s), σ_y(s), σ_z(s) at on-crest

Phase 2 (TODO; out of scope for first pass): enable SC in both codes
(RF-Track SpaceCharge_PIC_FreeSpace, xfields SpaceChargeBiGaussian)
and compare moment evolution. Requires Phase 1 agreement first.

Eremey Valetov, 2026-05-04
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import RF_Track as rft
import xtrack as xt
import xpart as xp

# Physical constants (SI; energies in MeV elsewhere).
MC2_MEV = 0.510998950
CLIGHT = 299792458.0
QE = 1.602176634e-19

# SLAC 3-m S-band TW structure (UH FEL operating point).
FREQ_HZ = 2.856e9
PHASE_ADV = 2.0 * np.pi / 3.0
PEAK_GRADIENT = 13.3e6  # V/m
L_TARGET = 3.048
K_INJECT = 1.0  # MeV

# Cell length from synchronous β=1 phase advance.
L_CELL = CLIGHT * PHASE_ADV / (2.0 * np.pi * FREQ_HZ)  # ≈ 34.99 mm
N_CELLS_NOMINAL = int(round(L_TARGET / L_CELL))  # 87

# Phase shift mapping RF-Track phid (autophased to crest) to elegant
# convention (peak ~70° at 1 MeV injection).
PHASE_SHIFT_DEG = 70.0

WORK_DIR = Path(__file__).parent
ELEG_CSV = WORK_DIR.parent / "elegant_linac" / "phase_scan_results.csv"


# ----------------------------------------------------------------------
# Lattice builders
# ----------------------------------------------------------------------


def rft_structure(phi_deg: float = 0.0) -> rft.TW_Structure:
    """RF-Track TW_Structure, single Fourier coefficient, peak-field model."""
    n_cells = L_TARGET / L_CELL
    tw = rft.TW_Structure(PEAK_GRADIENT, 0, FREQ_HZ, PHASE_ADV, n_cells)
    tw.set_phid(phi_deg)
    return tw


def xsuite_cavity_chain(phi_deg: float, n_cells: int = N_CELLS_NOMINAL) -> xt.Line:
    """Build an xsuite Line modelling the SLAC linac as a single lumped Cavity.

    A multi-cavity chain (one xtrack.Cavity per cell) is unfaithful at 1 MeV
    injection: each Cavity is static-phase, so low-β phase slippage between
    cells (β_inj=0.94, ΔΦ_slip ≈ 8°/cell vs the 120° synchronous advance)
    sends most cells off-crest and the integrated K_out collapses to ~7 MeV.
    Fixing this would require per-cell autophasing (what RF-Track's
    TW_Structure does internally), which defeats the purpose of using
    xsuite as an *independent* check.

    Lumped-cavity model: one Cavity with V_total = E_0 · L_target
    sandwiched by L_target/2 drifts. This is equivalent to elegant's RFCA
    element. It captures the integrated energy gain accurately; what it
    *does not* capture is the TW spatial extension or longitudinal phase
    slippage along the structure. Disagreement with RF-Track is the TW
    correction to RFCA.

    Args:
        phi_deg: lag angle for the lumped cavity, in degrees, elegant
            convention (peak ≈ 70° at 1 MeV injection in this model).
        n_cells: kept for API; unused (lumped model).

    Returns:
        xtrack Line: drift_half, cavity, drift_half.
    """
    v_total = PEAK_GRADIENT * L_TARGET  # V

    elements = [
        xt.Drift(length=L_TARGET / 2),
        xt.Cavity(voltage=v_total, frequency=FREQ_HZ, lag=phi_deg),
        xt.Drift(length=L_TARGET / 2),
    ]
    line = xt.Line(elements=elements)
    line.particle_ref = xp.Particles(
        mass0=MC2_MEV * 1e6,  # eV
        q0=-1.0,
        kinetic_energy0=K_INJECT * 1e6,  # eV
    )
    return line


# ----------------------------------------------------------------------
# Single-particle phase scan
# ----------------------------------------------------------------------


def rft_track_one(phi_deg: float) -> float:
    tw = rft_structure(phi_deg=phi_deg)
    lat = rft.Lattice()
    lat.append(tw)
    P_in = np.sqrt((K_INJECT + MC2_MEV) ** 2 - MC2_MEV**2)
    bunch = rft.Bunch6d(
        MC2_MEV, 1.0, -1.0, np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P_in]])
    )
    bunch_out = lat.track(bunch)
    M = bunch_out.get_phase_space("%x %xp %y %yp %t %Pc")
    if M.shape[0] == 0:
        return float("nan")
    P_out = M[0, 5]
    return float(np.sqrt(P_out**2 + MC2_MEV**2) - MC2_MEV)


def xsuite_track_one(phi_deg_elegant: float) -> float:
    """Single-particle xsuite scan; phi_deg_elegant is in elegant convention."""
    line = xsuite_cavity_chain(phi_deg=phi_deg_elegant)
    p = line.particle_ref.copy()
    line.track(p)
    P_out_eV = float(p.energy[0]) - line.particle_ref.energy0[0]
    # `p.energy` is total energy (eV); kinetic = total - rest mass × c²
    K_out_eV = float(p.energy[0]) - MC2_MEV * 1e6
    return K_out_eV / 1e6


def phase_scan_compare(
    n_phases: int = 73, out_pdf: Path | None = None
) -> dict[str, np.ndarray]:
    """Three-code phase scan: RF-Track, xsuite, elegant (overlay)."""
    phi_rft = np.linspace(-180, 180, n_phases)
    phi_eleg = (phi_rft + PHASE_SHIFT_DEG) % 360.0  # elegant-convention phase

    print(f"RF-Track phase scan, {n_phases} points...", flush=True)
    K_rft = np.array([rft_track_one(p) for p in phi_rft])

    print(f"xsuite cavity-chain phase scan, {n_phases} points...", flush=True)
    K_xs = np.array([xsuite_track_one(p) for p in phi_eleg])

    K_eleg = None
    if ELEG_CSV.exists():
        eleg = np.loadtxt(ELEG_CSV, delimiter=",", skiprows=1)
        K_eleg = (eleg[:, 0], eleg[:, 1])  # (phase, K_out)

    if out_pdf is not None:
        fig, ax = plt.subplots(figsize=(10, 6))
        ordr = np.argsort(phi_eleg)
        ax.plot(
            phi_eleg[ordr],
            K_rft[ordr],
            "b.-",
            ms=5,
            lw=1.2,
            label="RF-Track TW_Structure (TW model)",
        )
        ax.plot(
            phi_eleg[ordr],
            K_xs[ordr],
            "g.--",
            ms=5,
            lw=1.2,
            label="xsuite lumped Cavity (RFCA-equivalent)",
        )
        if K_eleg is not None:
            ax.plot(
                K_eleg[0],
                K_eleg[1],
                "r-",
                lw=1.5,
                alpha=0.7,
                label="elegant RFCA reference",
            )
        ax.set_xlabel("Phase [deg], elegant convention")
        ax.set_ylabel("K_out [MeV]")
        ax.set_xlim(0, 360)
        ax.set_title(
            "SLAC 3-m TW Linac single-particle phase scan\n"
            f"f={FREQ_HZ/1e9:.3f} GHz, E_0={PEAK_GRADIENT/1e6:.1f} MV/m, "
            f"K_inj={K_INJECT} MeV"
        )
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf, dpi=150, bbox_inches="tight")
        fig.savefig(out_pdf.with_suffix(".png"), dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"  -> {out_pdf}")

    valid_rft = ~np.isnan(K_rft)
    valid_xs = ~np.isnan(K_xs)
    print(f"\nPeaks:")
    if valid_rft.any():
        i = int(np.nanargmax(K_rft))
        print(f"  RF-Track: K_max={K_rft[i]:.4f} MeV at phi_eleg={phi_eleg[i]:.2f} deg")
    if valid_xs.any():
        i = int(np.nanargmax(K_xs))
        print(f"  xsuite  : K_max={K_xs[i]:.4f} MeV at phi_eleg={phi_eleg[i]:.2f} deg")
    if K_eleg is not None:
        i = int(np.argmax(K_eleg[1]))
        print(f"  elegant : K_max={K_eleg[1][i]:.4f} MeV at phi={K_eleg[0][i]:.2f} deg")

    return {
        "phase_eleg": phi_eleg,
        "K_rft": K_rft,
        "K_xs": K_xs,
    }


# ----------------------------------------------------------------------
# Multi-particle Gaussian bunch -- no SC, sigma evolution
# ----------------------------------------------------------------------


def make_bunch_arrays(
    n_part: int,
    sig_x: float,
    sig_xp: float,
    sig_y: float,
    sig_yp: float,
    sig_z: float,
    sig_delta: float,
    seed: int,
) -> dict[str, np.ndarray]:
    """Generate a fixed-seed Gaussian bunch in canonical (x, xp, y, yp, z, delta)
    where delta = (P - P_ref)/P_ref (xsuite convention; not dE/E_total --
    relativistic conversion factor 1/beta^2 differs by ~13% for 1 MeV electrons)."""
    rng = np.random.default_rng(seed)
    return {
        "x": sig_x * rng.standard_normal(n_part),
        "xp": sig_xp * rng.standard_normal(n_part),
        "y": sig_y * rng.standard_normal(n_part),
        "yp": sig_yp * rng.standard_normal(n_part),
        "z": sig_z * rng.standard_normal(n_part),
        "delta": sig_delta * rng.standard_normal(n_part),
    }


def rft_track_bunch_traces(
    bunch_arrays: dict[str, np.ndarray], n_samples: int, phi_deg: float = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    """Track a Gaussian bunch through RF-Track TW_Structure, sampled at n_samples
    sub-lattice exit points along s. Returns (s_array, sigmas_array of shape
    [n_samples+1, 3] for σ_x, σ_y, σ_z in mm)."""
    K_total = K_INJECT + MC2_MEV
    P_ref = np.sqrt(K_total**2 - MC2_MEV**2)

    # RF-Track Bunch6d: [X XP Y YP T P]  with units [mm mrad mm mrad mm/c MeV/c]
    # Convert: x→mm, xp→mrad, y→mm, yp→mrad, z→T as t = -z/c (with t in mm/c
    # = mm of length / c → set as -z (mm), since RFT's t-coordinate units are mm/c
    # and "mm/c" multiplied by c gives mm, so numerically t = -z when both in mm).
    # δE/E maps to (P-P_ref)/P_ref to first order; convert to P:
    n_part = len(bunch_arrays["x"])
    P_per = P_ref * (1.0 + bunch_arrays["delta"])
    phase_space = np.column_stack(
        [
            bunch_arrays["x"] * 1e3,  # mm
            bunch_arrays["xp"] * 1e3,  # mrad
            bunch_arrays["y"] * 1e3,  # mm
            bunch_arrays["yp"] * 1e3,  # mrad
            # RF-Track t [mm/c] = arrival time at S; for a leading particle
            # (z > 0) t < 0. Strict relation t = -z/(beta*c), so the [mm/c]
            # numerical value is -z[mm]/beta. We drop the 1/beta factor here
            # (~6% systematic at 1 MeV electrons, beta=0.94); the no-SC
            # comparison is qualitative and the bunch is small around the
            # reference, so this is acceptable for a sanity check.
            -bunch_arrays["z"] * 1e3,  # mm  (approximate; see comment)
            P_per,  # MeV/c
        ]
    )
    bunch_in = rft.Bunch6d(MC2_MEV, 1.0, -1.0, phase_space)

    # Sample at n_samples sub-lattice exit points by tracking through
    # progressively-longer fractions of one TW structure (matches the
    # benchmark_vs_elegant.py technique).
    s_arr = np.linspace(L_TARGET / n_samples, L_TARGET, n_samples)
    sigmas = np.zeros((n_samples + 1, 3))
    sigmas[0] = (
        np.std(phase_space[:, 0]),
        np.std(phase_space[:, 2]),
        np.std(phase_space[:, 4]),
    )

    for i, s_end in enumerate(s_arr, start=1):
        n_cells_partial = s_end / L_CELL
        tw = rft.TW_Structure(
            PEAK_GRADIENT, 0, FREQ_HZ, PHASE_ADV, n_cells_partial
        )
        tw.set_phid(phi_deg)
        lat = rft.Lattice()
        lat.append(tw)
        bunch_partial = rft.Bunch6d(MC2_MEV, 1.0, -1.0, phase_space)
        bunch_out = lat.track(bunch_partial)
        M = bunch_out.get_phase_space("%x %xp %y %yp %t %Pc")
        if M.shape[0] == 0:
            sigmas[i] = (np.nan, np.nan, np.nan)
        else:
            sigmas[i] = (np.std(M[:, 0]), np.std(M[:, 2]), np.std(M[:, 4]))

    s_full = np.concatenate(([0.0], s_arr))
    return s_full, sigmas


def xsuite_track_bunch_traces(
    bunch_arrays: dict[str, np.ndarray], phi_deg_elegant: float
) -> tuple[np.ndarray, np.ndarray]:
    """Track a Gaussian bunch through the xsuite cavity-chain, capturing σ
    after each element. Returns (s_array, sigmas in mm)."""
    line = xsuite_cavity_chain(phi_deg=phi_deg_elegant)
    n_part = len(bunch_arrays["x"])

    # Build particles from the canonical Gaussian arrays.
    # xsuite δ = (P - P_ref)/P_ref; to first order δ ≈ δE/(β² E_total) -- use
    # `delta=` to set directly. For 1 MeV electrons, β=0.94, β²·E_total ≈ 1.34 MeV.
    # Better: convert δE/E_total → δ via xpart factor; for now use delta = δE/E.
    p = xp.Particles(
        _context=line._context,
        mass0=MC2_MEV * 1e6,
        q0=-1.0,
        kinetic_energy0=K_INJECT * 1e6,
        x=bunch_arrays["x"],
        px=bunch_arrays["xp"],  # px ≈ xp at this energy
        y=bunch_arrays["y"],
        py=bunch_arrays["yp"],
        zeta=bunch_arrays["z"],
        delta=bunch_arrays["delta"],
    )

    # Insert monitors after each Drift+Cavity pair to sample σ.
    # Simpler: track step-by-step using `ele_start`/`ele_stop` slicing.
    s_positions = []
    sigmas = []
    s_cum = 0.0
    s_positions.append(s_cum)
    sigmas.append((np.std(p.x), np.std(p.y), np.std(p.zeta)))

    for i_el, el in enumerate(line.elements):
        line.track(p, ele_start=i_el, num_elements=1)
        s_cum += el.length
        s_positions.append(s_cum)
        sigmas.append((np.std(p.x), np.std(p.y), np.std(p.zeta)))

    return np.asarray(s_positions), np.asarray(sigmas)


def bunch_compare(out_pdf: Path) -> None:
    """Compare RF-Track vs xsuite σ evolution along the linac."""
    n_part = 1000
    sig_x = sig_y = 0.5e-3  # m  (nominal injector σ; placeholder)
    sig_xp = sig_yp = 0.1e-3  # rad
    sig_z = 0.6e-3  # m
    sig_delta = 5e-3  # σ_p/p (xsuite convention); spch_demo SD0 nominally 5e-3
    seed = 20260504

    bunch = make_bunch_arrays(
        n_part, sig_x, sig_xp, sig_y, sig_yp, sig_z, sig_delta, seed
    )

    print("\nBunch tracking -- no SC")
    print("  RF-Track (sub-lattice exit points), phid=0 (autophased)...", flush=True)
    s_rft, sig_rft_mm = rft_track_bunch_traces(bunch, n_samples=15, phi_deg=0.0)

    # xsuite lumped Cavity: peak at lag=90° in xsuite convention.
    print("  xsuite lumped Cavity, lag=90 (xsuite peak)...", flush=True)
    s_xs, sig_xs_m = xsuite_track_bunch_traces(bunch, phi_deg_elegant=90.0)

    # RF-Track sigmas come back in mm; xsuite returns m. Normalize to mm.
    sig_xs_mm = sig_xs_m * 1e3

    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    labels = [r"$\sigma_x$ [mm]", r"$\sigma_y$ [mm]", r"$\sigma_z$ [mm]"]
    for k in range(3):
        ax = axes[k]
        ax.plot(s_rft, sig_rft_mm[:, k], "b.-", lw=1.2, label="RF-Track TW_Structure")
        ax.plot(
            s_xs, sig_xs_mm[:, k], "g.--", lw=1.2, label="xsuite cavity-chain"
        )
        ax.set_ylabel(labels[k])
        ax.grid(alpha=0.3)
        if k == 0:
            ax.legend(fontsize=9, loc="best")
    axes[-1].set_xlabel("s [m]")
    fig.suptitle(
        "Gaussian bunch through SLAC linac -- RF-Track vs xsuite (no SC)\n"
        f"N={n_part}, σ_x=σ_y=0.5 mm, σ_z=0.6 mm, σ_p/p=5e-3, on-crest"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, dpi=150, bbox_inches="tight")
    fig.savefig(out_pdf.with_suffix(".png"), dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {out_pdf}")


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--what",
        choices=["phase", "bunch", "all"],
        default="all",
        help="Which comparison to run",
    )
    ap.add_argument(
        "--nphase",
        type=int,
        default=73,
        help="Number of single-particle phase-scan points",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=WORK_DIR / "xsuite_compare_output",
        help="Directory for figures + CSVs",
    )
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.what in ("phase", "all"):
        scan = phase_scan_compare(
            n_phases=args.nphase, out_pdf=args.out_dir / "phase_scan.pdf"
        )
        np.savetxt(
            args.out_dir / "phase_scan.csv",
            np.column_stack([scan["phase_eleg"], scan["K_rft"], scan["K_xs"]]),
            header="phase_eleg_deg,K_rft_MeV,K_xs_MeV",
            delimiter=",",
            comments="",
        )

    if args.what in ("bunch", "all"):
        bunch_compare(args.out_dir / "bunch_no_sc.pdf")

    print(f"\nAll outputs in: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
