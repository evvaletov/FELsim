#!/usr/bin/env python3
"""xtrack-native multi-cell TW linac, driven by the validated TW model.

xtrack/xsuite is fixed-reference (storage-ring) and cannot, on its own, track a
1->40 MeV linac: a cavity chain with a fixed reference mistimes the RF once the
particle accelerates, and the canonical (x', delta) coordinates blow up. The fix
is to ramp the reference momentum cell-by-cell with xt.ReferenceEnergyIncrease,
paired with an xt.Cavity that delivers the per-cell energy gain. Here the per-cell
gains come from the validated autophasing model (linac_multicell_tw.py), so the
xtrack line reproduces the energy profile and provides the 6D transport (transverse
adiabatic damping) that the standalone integrator does not -- and is ready for the
xsuite space-charge path.

Validation: (1) on-axis K_out matches the integrator; (2) the transverse map
det(R_x) equals p_in/p_out (adiabatic damping).

Author: Eremey Valetov
"""

import sys
from pathlib import Path
import numpy as np
import xtrack as xt
import xpart as xp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from linac_multicell_tw import (
    CLIGHT, MC2, OMEGA, FREQ, E0, K_INJECT, L_CELL, TAU, R_SHUNT,
    load_cells, Eb_profile,
)

OUT = HERE.parent / "results" / "niels_june"
EV = 1e6  # MeV -> eV


def synchronous_profile(phi_inj_deg, n_cells, I_amp=0.0, nsub=20):
    """Integrate the synchronous electron cell-by-cell (autophasing, optional
    beam loading); return K [MeV] at each cell boundary, length n_cells+1."""
    L_total = n_cells * L_CELL
    K = K_INJECT
    psi = np.deg2rad(phi_inj_deg)
    Ks = [K]

    def deriv(K, psi, z):
        g = 1 + K / MC2
        if g <= 1:
            return None
        b = np.sqrt(1 - 1 / g ** 2)
        Eb = Eb_profile(z / L_total, I_amp)
        return (E0 * np.cos(psi) - Eb) / 1e6, (OMEGA / CLIGHT) * (1 / b - 1)

    def pad(lst):
        return np.array(lst + [np.nan] * (n_cells + 1 - len(lst)))

    for n in range(n_cells):
        zs = np.linspace(n * L_CELL, (n + 1) * L_CELL, nsub + 1)
        for i in range(nsub):
            z, h = zs[i], zs[i + 1] - zs[i]
            k1 = deriv(K, psi, z)
            if k1 is None:
                return pad(Ks)
            k2 = deriv(K + .5 * h * k1[0], psi + .5 * h * k1[1], z + .5 * h)
            if k2 is None:
                return pad(Ks)
            k3 = deriv(K + .5 * h * k2[0], psi + .5 * h * k2[1], z + .5 * h)
            if k3 is None:
                return pad(Ks)
            k4 = deriv(K + h * k3[0], psi + h * k3[1], z + h)
            if k4 is None:
                return pad(Ks)
            K += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            psi += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        Ks.append(K)
    return np.array(Ks)


def momentum_ev(K_mev):
    return np.sqrt((K_mev + MC2) ** 2 - MC2 ** 2) * EV


def build_linac(Ks):
    """xtrack Line from per-cell synchronous energies Ks [MeV] (len n_cells+1).
    Each cell: half-drift, on-crest Cavity (V = dK), ReferenceEnergyIncrease
    (matched dp so the synchronous particle stays at delta~0), half-drift."""
    elems = []
    for n in range(len(Ks) - 1):
        dK = (Ks[n + 1] - Ks[n]) * EV                       # eV
        dp = momentum_ev(Ks[n + 1]) - momentum_ev(Ks[n])    # eV/c
        elems += [
            xt.Drift(length=L_CELL / 2),
            xt.Cavity(frequency=FREQ, voltage=dK, phase=np.pi / 2),
            xt.ReferenceEnergyIncrease(Delta_p0c=dp),
            xt.Drift(length=L_CELL / 2),
        ]
    line = xt.Line(elements=elems)
    line.particle_ref = xp.Particles(mass0=MC2 * EV, q0=-1.0,
                                     kinetic_energy0=K_INJECT * EV)
    line.build_tracker()
    return line


def transverse_R(line):
    """Extract the 2x2 horizontal map R_x from unit-perturbation tracking."""
    x0, a0 = 1e-6, 1e-6
    pa = line.build_particles(x=[x0, 0.0], px=[0.0, a0],
                              y=[0, 0], py=[0, 0], zeta=[0, 0], delta=[0, 0])
    line.track(pa)
    R = np.array([[pa.x[0] / x0, pa.x[1] / a0],
                  [pa.px[0] / x0, pa.px[1] / a0]])
    return R


def main():
    labels, _, _, _ = load_cells()
    n_cells = len(labels)

    # optimal injection phase from the validated model
    phis = np.arange(0, 360, 1.0)
    Kend = np.array([synchronous_profile(p, n_cells)[-1] for p in phis])
    phipk = phis[int(np.nanargmax(np.where(np.isfinite(Kend), Kend, -np.inf)))]
    Ks = synchronous_profile(phipk, n_cells)
    K_model = Ks[-1]

    line = build_linac(Ks)
    p = line.build_particles(x=[0.0], px=[0.0], y=[0.0], py=[0.0],
                             zeta=[0.0], delta=[0.0])
    line.track(p)
    K_xtrack = (p.energy[0] - MC2 * EV) / EV          # MeV
    delta_final = p.delta[0]

    R = transverse_R(line)
    detR = R[0, 0] * R[1, 1] - R[0, 1] * R[1, 0]
    detR_expected = momentum_ev(K_INJECT) / momentum_ev(K_model)

    print(f"optimal phase {phipk:.0f} deg, {n_cells} cells")
    print(f"on-axis K_out: integrator {K_model:.3f} MeV, "
          f"xtrack {K_xtrack:.3f} MeV "
          f"({abs(K_xtrack-K_model)/K_model*100:.3f}%); final delta {delta_final:.2e}")
    print(f"transverse det(R_x): xtrack {detR:.4f}, "
          f"expected p_in/p_out {detR_expected:.4f} "
          f"({abs(detR-detR_expected)/detR_expected*100:.2f}%)")

    md = [
        "# xtrack-native multi-cell TW linac",
        "",
        "Author: Eremey Valetov, 2026-06-20. The validated autophasing TW model",
        "(linac_multicell_tw.py) drives an xtrack Line via per-cell Cavity +",
        "ReferenceEnergyIncrease (energy-ramp reference), so xtrack can track the",
        "1->40 MeV linac it otherwise cannot, and provides the 6D transport for the",
        "xsuite space-charge path.",
        "",
        "| Check | xtrack | reference | agreement |",
        "|---|---:|---:|---:|",
        f"| on-axis K_out (MeV) | {K_xtrack:.3f} | {K_model:.3f} (integrator) | "
        f"{abs(K_xtrack-K_model)/K_model*100:.3f}% |",
        f"| transverse det(R_x) | {detR:.4f} | {detR_expected:.4f} (p_in/p_out) | "
        f"{abs(detR-detR_expected)/detR_expected*100:.2f}% |",
        f"| final delta (synchronous) | {delta_final:.1e} | 0 | - |",
        "",
        f"{n_cells} cells, optimal injection phase {phipk:.0f} deg. The reference",
        "ramps with the synchronous energy (delta stays ~0), and the transverse map",
        "shows the correct adiabatic damping det(R_x)=p_in/p_out. This Line is a",
        "drop-in for the xsuite SC engines (frozen / PIC). Beam loading is modelled",
        "in this standalone validation script via synchronous_profile(..., I_amp=...);",
        "the production XsuiteAdapter cavity is currently unloaded (zero beam current).",
        "Remaining: true-phase cavities for self-consistent longitudinal bunching",
        "(this build is on-crest, energy-exact).",
    ]
    (OUT / "xsuite_linac_tw.md").write_text("\n".join(md))

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    z = np.arange(len(Ks)) * L_CELL
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(z, Ks, 'C0-', lw=2)
    ax.set_xlabel('z (m)')
    ax.set_ylabel('synchronous kinetic energy (MeV)')
    ax.set_title(f'xtrack-native TW linac: 1 -> {K_xtrack:.1f} MeV '
                 f'(det Rx={detR:.4f}, delta~{delta_final:.0e})')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f"xsuite_linac_tw.{ext}", dpi=150)
    print(f"\nWrote {OUT}/xsuite_linac_tw.md and xsuite_linac_tw.png")


if __name__ == "__main__":
    main()
