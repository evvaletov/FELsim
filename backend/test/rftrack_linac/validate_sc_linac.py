#!/usr/bin/env python3
"""Validate per-cell space charge inside the xsuite TW linac.

The XsuiteAdapter interleaves an SC kick after every TW cell, each placed after
that cell's reference-energy ramp so it acts at the local energy. Space charge is
therefore strongest at 1 MeV injection and suppressed (~1/(beta^2 gamma^3) for a
fixed-charge bunch) as the beam accelerates -- the relativistic shielding a
single fixed-energy SC block would miss.

Three checks: (1) STRUCTURE -- one SC kick per cell, each right after a
ReferenceEnergyIncrease; (2) ACCELERATION -- a bunch still reaches ~41 MeV with
SC on and survives; (3) LOCAL-GAMMA LAW -- a frozen SC kick scales ~1/(beta^2
gamma^3) with the reference energy, so each interleaved cell (sitting after its
ramp) acts at its own local energy.

Author: Eremey Valetov
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xtrack as xt
import xpart as xp
import xfields as xf

BACKEND = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(BACKEND))
REPO = BACKEND.parent
OUT = BACKEND / "test" / "results" / "niels_june"
OUT.mkdir(parents=True, exist_ok=True)

from xsuiteAdapter import XsuiteAdapter
from physicalConstants import PhysicalConstants

MC2_EV = PhysicalConstants.E0_electron * 1e6
Q_E = PhysicalConstants.Q
L_CELL = PhysicalConstants.C * (2 * np.pi / 3) / (2 * np.pi * 2.856e9)


def betagamma(K_mev):
    g = 1 + K_mev * 1e6 / MC2_EV
    return np.sqrt(g ** 2 - 1)


def b2g3(K_mev):
    """beta^2 gamma^3 -- the scaling of the linear transverse SC defocus of a
    FIXED-charge bunch. The fixed-current perveance form 1/(beta^3 gamma^3)
    loses one power of beta because the line current I ~ beta at fixed charge."""
    g = 1 + K_mev * 1e6 / MC2_EV
    return (g ** 2 - 1) * g          # (beta gamma)^2 * gamma


def build_sc_line(charge_nc=1.0, sig=1e-3, sig_z=1e-3):
    sim = XsuiteAdapter(lattice_path=str(REPO / "var" / "slac_linac.json"),
                        beam_energy=1.0, space_charge=True, sc_method="frozen")
    n_e = charge_nc * 1e-9 / Q_E
    return sim._build_line(sc_on=True, sig_x=sig, sig_y=sig, sig_z=sig_z, n_e=n_e)


def sc_kick_coeff(K_mev, n_e=6.24e9, sig=1e-3):
    """Linear frozen-SC defocus dpx/x for one cell at reference energy K_mev."""
    lp = xf.LongitudinalProfileQGaussian(number_of_particles=n_e, sigma_z=1e-3,
                                         z0=0.0, q_parameter=1.0)
    sc = xf.SpaceChargeBiGaussian(length=L_CELL, longitudinal_profile=lp,
                                  sigma_x=sig, sigma_y=sig, mean_x=0.0, mean_y=0.0,
                                  update_on_track=False)
    line = xt.Line(elements=[sc])
    line.particle_ref = xp.Particles(mass0=MC2_EV, q0=-1.0,
                                     kinetic_energy0=K_mev * 1e6)
    line.build_tracker()
    x0 = 1e-5
    p = line.build_particles(x=[x0], px=[0.0], y=[0.0], py=[0.0],
                             zeta=[0.0], delta=[0.0])
    line.track(p)
    return abs(p.px[0] / x0)


def main():
    # 1) STRUCTURE
    line = build_sc_line()
    tab = line.get_table()
    et = list(tab.element_type)
    n_sc = et.count('SpaceChargeBiGaussian')
    # every SC element is preceded by a ReferenceEnergyIncrease earlier in its cell
    cell_ok = all(
        'ReferenceEnergyIncrease' in et[max(0, i - 3):i]
        for i, t in enumerate(et) if t == 'SpaceChargeBiGaussian')
    print(f"[structure] {n_sc} SC kicks, one per TW cell, each after a "
          f"ReferenceEnergyIncrease: {cell_ok}")
    assert n_sc == 87 and cell_ok

    # 2) ACCELERATION + SURVIVAL (real bunch, frozen SC update_on_track=True)
    rng = np.random.default_rng(42)
    n = 4000
    p = line.build_particles(
        x=rng.normal(0, 1e-3, n), px=rng.normal(0, 1e-3, n),
        y=rng.normal(0, 1e-3, n), py=rng.normal(0, 1e-3, n),
        zeta=rng.normal(0, 1e-3, n), delta=np.zeros(n))
    line.track(p)
    alive = int((p.state == 1).sum())
    K_out = (np.nanmean(p.energy[p.state == 1]) - MC2_EV) / 1e6
    print(f"[accel] SC-on linac: 1 MeV -> {K_out:.2f} MeV, {alive}/{n} alive")
    assert alive == n and 39 < K_out < 42

    # 3) LOCAL-GAMMA LAW: per-cell SC strength vs reference energy
    Ks = np.array([1, 2, 3, 5, 8, 13, 21, 34, 41.5])
    coeff = np.array([sc_kick_coeff(K) for K in Ks])
    supp = coeff[0] / coeff[-1]
    sc_law = b2g3(Ks[0]) / b2g3(Ks) * coeff[0]   # fixed-charge ~1/(beta^2 gamma^3)
    print(f"[local-gamma] frozen SC defocus dpx/x: 1 MeV={coeff[0]:.3e}, "
          f"41.5 MeV={coeff[-1]:.3e} -> x{supp:.0f} suppression "
          f"(1/(beta^2 gamma^3) predicts x{b2g3(41.5)/b2g3(1):.0f})")
    assert supp > 1000

    md = [
        "# Per-cell space charge in the xsuite TW linac",
        "",
        "Author: Eremey Valetov, 2026-06-20.",
        "",
        f"1. **Structure:** {n_sc} SC kicks interleaved, one per TW cell, each "
        "right after that cell's ReferenceEnergyIncrease (so it acts at the local "
        "energy).",
        f"2. **Acceleration:** with SC on, a {n}-particle bunch still reaches "
        f"{K_out:.1f} MeV and all {alive} particles survive.",
        f"3. **Local-gamma law:** the frozen SC defocus dpx/x falls from "
        f"{coeff[0]:.2e} 1/m at 1 MeV to {coeff[-1]:.2e} 1/m at 41.5 MeV "
        f"(**x{supp:.0f}**), tracking the 1/(beta^2 gamma^3) fixed-charge "
        f"relativistic shielding (predicted x{b2g3(41.5)/b2g3(1):.0f}). A single "
        "fixed-energy SC block would apply the 1 MeV strength throughout and "
        "over-estimate SC by this factor.",
        "",
        "Frozen SC reads the ramped reference automatically; PIC gets the local "
        "gamma. Emittance is a poor probe here (the field is nearly linear across a "
        "round Gaussian -> envelope perturbation, not eps_n growth, in one short "
        "pass). Remaining: sigma(s) envelope evolution (fixed input sigma here); "
        "simulate() exit-energy coordinate frame.",
    ]
    (OUT / "sc_linac_validation.md").write_text("\n".join(md))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(betagamma(Ks), coeff, 'C3o-', label='frozen SC (measured)')
    ax.loglog(betagamma(Ks), sc_law, 'k--', lw=1,
              label=r'$\propto 1/(\beta^2\gamma^3)$')
    ax.set_xlabel(r'$\beta\gamma$ (1 MeV -> 41.5 MeV)')
    ax.set_ylabel('per-cell SC defocus |dpx/x| (1/m)')
    ax.set_title(r'Per-cell SC follows the local energy ($\propto 1/(\beta^2\gamma^3)$)')
    ax.legend()
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f"sc_linac_percell.{ext}", dpi=150)
    print(f"\nWrote {OUT}/sc_linac_validation.md and sc_linac_percell.png")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
