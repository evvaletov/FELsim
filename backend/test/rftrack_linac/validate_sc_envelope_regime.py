#!/usr/bin/env python3
"""Regime of validity of the prescribed matched-envelope SC, + multi-cavity sig_env.

The optional prescribed-envelope SC path uses a fixed per-cell sigma from a
pre-pass, so it cannot self-correct. It is therefore accurate only when the beam
stays close to that envelope -- a stable / weakly-perturbed / matched beam -- and
diverges when strong unmatched space charge moves the beam away from it.

Two checks:
  1. CHARGE SCAN on the SLAC linac: prescribed vs self-consistent agreement as a
     function of bunch charge. At low charge (weak SC, stable beam) they agree; at
     high charge (strong unmatched SC) they diverge. This is the regime of
     validity. Transverse focusing on the real injector+linac would keep the beam
     near the envelope and extend the agreement -- that lattice-design validation
     is deferred (a quick synthetic accelerating FODO was not matched/stable).
  2. MULTI-CAVITY sig_env: a 2-cavity line builds with a flat envelope indexed by
     a running offset, each cavity's SC using the correct slice.

Author: Eremey Valetov
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BACKEND = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(BACKEND))
REPO = BACKEND.parent
OUT = BACKEND / "test" / "results" / "niels_june"
OUT.mkdir(parents=True, exist_ok=True)

from xsuiteAdapter import XsuiteAdapter
from simulatorBase import BeamlineElement
from physicalConstants import PhysicalConstants

MC2_EV = PhysicalConstants.E0_electron * 1e6
Q_E = PhysicalConstants.Q
SIG, SIGP, SIGZ = 1e-3, 1e-3, 1e-3


def adapter():
    return XsuiteAdapter(lattice_path=str(REPO / "var" / "slac_linac.json"),
                         beam_energy=1.0, space_charge=True, sc_method="frozen")


def eps_n(line, seed=1, n=3000):
    rng = np.random.default_rng(seed)
    p = line.build_particles(
        x=rng.normal(0, SIG, n), px=rng.normal(0, SIGP, n),
        y=rng.normal(0, SIG, n), py=rng.normal(0, SIGP, n),
        zeta=rng.normal(0, SIGZ, n), delta=np.zeros(n))
    line.track(p)
    m = p.state == 1
    bg = p.p0c[0] / MC2_EV
    eg = np.sqrt(max(np.var(p.x[m]) * np.var(p.px[m]) - np.cov(p.x[m], p.px[m])[0, 1] ** 2, 0))
    return bg * eg * 1e6


def main():
    sim = adapter()

    # 1) charge scan: prescribed vs self-consistent agreement
    charges = [0.001, 0.01, 0.03, 0.1, 0.3, 1.0]   # nC
    rel = []
    for q in charges:
        n_e = q * 1e-9 / Q_E
        line_sc = sim._build_line(sc_on=True, sig_x=SIG, sig_y=SIG, sig_z=SIGZ, n_e=n_e)
        en_sc = eps_n(line_sc)
        env = sim.sc_envelope_prepass(SIG, SIG, SIGP, SIGP, SIGZ, n_e)
        line_env = sim._build_line(sc_on=True, sig_x=SIG, sig_y=SIG, sig_z=SIGZ,
                                   n_e=n_e, sig_env=env)
        en_env = eps_n(line_env)
        r = abs(en_env - en_sc) / en_sc
        rel.append(r)
        print(f"[scan] q={q:>5} nC: prescribed eps_n={en_env:.3f}, "
              f"self-consistent={en_sc:.3f}, rel diff {r*100:5.1f}%")
    rel = np.array(rel)
    assert rel[0] < 0.05, "prescribed should agree at low charge (stable beam)"
    assert rel[-1] > rel[0] * 3, "prescribed should diverge as SC grows"
    print(f"[regime] agreement {rel[0]*100:.1f}% at {charges[0]} nC -> "
          f"{rel[-1]*100:.0f}% at {charges[-1]} nC: valid for stable beams, "
          "diverges for strong unmatched SC")

    # 2) multi-cavity sig_env: two CONTIGUOUS cavities, env spans both with the
    # running offset. (Contiguous so every SC cell is a cavity cell; a non-cavity
    # SC element between cavities would not be counted by the cavity-only offset --
    # that interspersed case is a documented limitation.)
    sim2 = XsuiteAdapter(beam_energy=1.0, space_charge=True, sc_method="frozen")
    sim2.beamline = [
        BeamlineElement('RF_CAVITY', 0.5, frequency_hz=2.856e9,
                        gradient_mv_per_m=13.3, phase_advance_deg=120.0),
        BeamlineElement('RF_CAVITY', 0.5, frequency_hz=2.856e9,
                        gradient_mv_per_m=13.3, phase_advance_deg=120.0)]
    n_e = 0.05 * 1e-9 / Q_E
    env2 = sim2.sc_envelope_prepass(SIG, SIG, SIGP, SIGP, SIGZ, n_e)
    line2 = sim2._build_line(sc_on=True, sig_x=SIG, sig_y=SIG, sig_z=SIGZ,
                             n_e=n_e, sig_env=env2)
    tab = line2.get_table()
    sc_names = [tab.name[i] for i, t in enumerate(tab.element_type)
                if t == 'SpaceChargeBiGaussian']
    used = [line2[nm].sigma_x for nm in sc_names]
    match = np.allclose(used, [e[0] for e in env2], rtol=1e-6)
    print(f"[multi-cavity] 2 RF_CAVITY -> {len(env2)} SC cells; per-cell sigma "
          f"matches the flat envelope with running offset: {match}")
    assert match and len(env2) == len(sc_names)

    md = [
        "# Prescribed-envelope SC: regime of validity + multi-cavity sig_env",
        "",
        "Author: Eremey Valetov, 2026-06-20.",
        "",
        "The optional prescribed-envelope SC uses a fixed per-cell sigma, so it is "
        "accurate only while the beam stays near that envelope.",
        "",
        "## Charge scan (SLAC linac)",
        "| charge (nC) | prescribed/self-consistent eps_n rel diff |",
        "|---:|---:|",
    ]
    for q, r in zip(charges, rel):
        md.append(f"| {q} | {r*100:.1f}% |")
    md += [
        "",
        f"At low charge ({charges[0]} nC, weak SC, stable beam) the prescribed "
        f"envelope agrees with self-consistent to {rel[0]*100:.1f}%; at {charges[-1]} "
        f"nC (strong unmatched SC) it diverges to {rel[-1]*100:.0f}% because a fixed "
        "sigma over-defocuses once the beam moves. Transverse focusing keeps the "
        "beam near the envelope and extends the agreement -- that focusing-lattice "
        "validation is deferred to the real injector+linac (a designed matched "
        "FODO; a quick synthetic one was not stable).",
        "",
        f"## Multi-cavity sig_env",
        f"Two contiguous cavities build with a flat {len(env2)}-cell envelope "
        "indexed by a running offset; each cavity's SC uses the correct slice "
        "(verified). Limitation: a non-cavity SC element (a quad/drift slice) "
        "between cavities is not counted by the cavity-only offset, so a fully "
        "interspersed accelerating-FODO would need the envelope to span all SC "
        "cells -- future work. update_on_track stays the default.",
    ]
    (OUT / "sc_envelope_regime.md").write_text("\n".join(md))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(charges, rel * 100, 'C3o-')
    ax.axhline(5, ls='--', color='gray', lw=1, label='5% agreement')
    ax.set_xlabel('bunch charge (nC)')
    ax.set_ylabel('prescribed vs self-consistent eps_n diff (%)')
    ax.set_title('Prescribed envelope: valid for stable beams, diverges with SC')
    ax.legend()
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f"sc_envelope_regime.{ext}", dpi=150)
    print(f"\nWrote {OUT}/sc_envelope_regime.md and sc_envelope_regime.png")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
