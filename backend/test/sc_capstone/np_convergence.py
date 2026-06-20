"""
N_p convergence of the DA-FMM space-charge growth on the no-dipole section, at the
headline point (45 MeV, 1 nC). The DA-FMM treecode carries a numerical-collisionality
(shot-noise) excess that scales as q_mp^0.46 with q_mp = Q/N_p; below the Phase-1
threshold q_mp* ~ 0.037 pC (N_p ~ 27000 at 1 nC) the shot noise drops below the physical
mean-field growth. This script sweeps N_p and shows the bare DA-FMM growth converging
toward the xsuite frozen-Gaussian (shot-noise-free) value as q_mp falls below q_mp*,
which is the fair physical comparison the per-charge sweep cannot make at fixed N_p.

Author: Eremey Valetov
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common_distribution as cd
import sc_section as scs
from nosc_handoff import LATTICE, NODIP_RANGE
from sc_capstone_run import make_beam, track_xsuite_frozen, growth_pct

Q_MP_STAR_PC = 0.037   # Phase-1 threshold (memory: sc_campaign)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--energy-mev", type=float, default=45.0)
    ap.add_argument("--q-nc", type=float, default=1.0)
    ap.add_argument("--np-list", nargs="+", type=int, default=[6000, 15000, 30000, 60000])
    ap.add_argument("--ds-target", type=float, default=0.02)
    ap.add_argument("--outdir", type=Path, default=Path("test/results/sc_capstone"))
    args = ap.parse_args()

    segs = scs.extract_segments(LATTICE, *NODIP_RANGE, args.energy_mev,
                                optics_energy_mev=args.energy_mev)
    rows = []
    for n_p in args.np_list:
        arr, _ = make_beam(args.energy_mev, n_p, 20260619, 8.0, args.q_nc)
        da = scs.track_dafmm(arr, segs, args.energy_mev, args.q_nc, args.ds_target, spch_on=True)
        xs = track_xsuite_frozen(arr, segs, args.energy_mev, args.q_nc, args.ds_target, spch_on=True)
        q_mp_pc = args.q_nc * 1e3 / n_p
        r = {"n_p": n_p, "q_mp_pc": q_mp_pc,
             "dafmm_gx": da["epsnx_growth"] * 100,
             "xsuite_gx": growth_pct(xs["epsnx_exit"], da["epsnx_init"])}
        rows.append(r)
        print(f"[npconv] N_p={n_p:6d}  q_mp={q_mp_pc:.4f} pC "
              f"({'>' if q_mp_pc > Q_MP_STAR_PC else '<'} q_mp*={Q_MP_STAR_PC}) "
              f"DA-FMM={r['dafmm_gx']:+.4f}%  xsuite={r['xsuite_gx']:+.4f}%", flush=True)

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "np_convergence.json").write_text(json.dumps(rows, indent=2))

    fig, ax = plt.subplots(figsize=(7, 5))
    nps = [r["n_p"] for r in rows]
    ax.plot(nps, [r["dafmm_gx"] for r in rows], "o-", label="DA-FMM (bare 1/r, with shot noise)")
    xs_mean = sum(r["xsuite_gx"] for r in rows) / len(rows)
    ax.axhline(xs_mean, ls="--", color="C1", label="xsuite frozen-Gaussian (shot-noise-free)")
    npstar = args.q_nc * 1e3 / Q_MP_STAR_PC
    ax.axvline(npstar, ls=":", color="gray", label=f"q_mp* (N_p~{npstar:.0f})")
    ax.set_xlabel("N_p (macroparticles)"); ax.set_ylabel(r"$\epsilon_{n,x}$ growth [%]")
    ax.set_xscale("log")
    ax.set_title(f"DA-FMM N_p convergence at {args.energy_mev:g} MeV, {args.q_nc:g} nC")
    ax.legend(); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.outdir / "np_convergence.png", dpi=140)
    print(f"[npconv] wrote {args.outdir}/np_convergence.json + .png")


if __name__ == "__main__":
    main()
