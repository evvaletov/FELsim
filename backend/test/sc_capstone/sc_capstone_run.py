"""
Reduced no-dipole space-charge capstone: cosy-fmm DA-FMM vs xsuite frozen-Gaussian
on the real FELsim transport section [32, 46), driven from ONE common distribution.

This is the SC-on counterpart of nosc_handoff.py. The no-SC handoff proved the codes
agree to <0.1% on this section at zero current (and that the full line is NOT yet
comparable because of dipole-model gaps), so the section is the right place to turn
space charge on. Both engines start from the identical common macroparticle array
(sigma_delta = 0 to isolate SC from chromatic growth), use the same split-operator
slicing, and report normalized-emittance growth vs charge and vs energy.

  * DA-FMM        — cosy-fmm Coulomb treecode (libspch_kick.so), full N-body 1/r.
  * xsuite frozen — xfields.SpaceChargeBiGaussian: the nonlinear Bassetti-Erskine field
    of an rms-self-consistent Gaussian (updated each step), NOT a linearized model.

Author: Eremey Valetov
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

# Point the cosy-pic adapter at the actual binary (its built-in default path is
# stale: build/pic-core-cli vs build/pic_core/pic-core-cli). Must be set BEFORE the
# simulatorFactory import chain so the 'pic' engine registers.
os.environ.setdefault("COSY_PIC_DIR", os.path.expanduser("~/COSY/cosy-pic"))
os.environ.setdefault(
    "PIC_CORE_CLI",
    os.path.expanduser("~/COSY/cosy-pic/build/pic_core/pic-core-cli"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import common_distribution as cd
import sc_section as scs
from nosc_handoff import LATTICE, NODIP_RANGE, run_code

QE = scs.QE


def track_pic(felsim_arr, energy_mev, q_nc, mesh=(64, 64, 64), optics_energy=None):
    """cosy-pic mesh-PIC space charge on the no-dipole section, same common
    distribution. Returns None if the 'pic' engine is unavailable (cosy-pic not
    built / env unset) or the run fails, so the column degrades gracefully.

    The PicAdapter tracks FELsim's native element maps at `energy_mev`, so it
    cannot apply the analytic optics re-match the DA-FMM/xsuite trackers use via
    `optics_energy` (it would track at the wrong, over-focused optics). Skip with
    a note when an optics re-match is requested, rather than report a mismatched
    column."""
    from simulatorFactory import SimulatorFactory
    import latticeLoader
    if "pic" not in getattr(SimulatorFactory, "_registry", {}):
        return None
    if optics_energy is not None and abs(optics_energy - energy_mev) > 1e-9:
        return {"error": f"optics re-match to {optics_energy} MeV not supported for pic"}
    try:
        sim = SimulatorFactory.create("pic", lattice_path=LATTICE,
                                      beam_energy=energy_mev,
                                      bunch_charge_nc=q_nc, sc_mesh=mesh)
        bl = latticeLoader.create_beamline(LATTICE)
        for e in bl:
            e.setE(energy_mev)
        sim.set_beamline(bl[NODIP_RANGE[0]:NODIP_RANGE[1]])
        sim.set_space_charge(True, mesh=mesh, bunch_charge_nc=q_nc)
        t0 = time.perf_counter()
        res = sim.simulate(particles=felsim_arr.copy())
        wall = time.perf_counter() - t0
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}
    if not res.success or res.final_particles is None:
        return {"error": "pic simulate failed"}
    p = res.final_particles
    if p.shape[0] == 0:
        return {"error": "no particles survived"}
    g, b = cd.gamma_beta(energy_mev)
    bg = b * g
    x, xp, y, yp = p[:, 0], p[:, 1], p[:, 2], p[:, 3]  # mm, mrad

    def epsn(u, up):
        du, dup = u - u.mean(), up - up.mean()
        return bg * np.sqrt(max(np.mean(du**2) * np.mean(dup**2) - np.mean(du * dup)**2, 0.0))
    return {"epsnx_exit": float(epsn(x, xp)), "epsny_exit": float(epsn(y, yp)),
            "wall": wall, "n_out": int(p.shape[0])}


def make_beam(energy_mev, n_p, seed, eps_n, q_nc):
    man = cd.BeamManifest(energy_mev=energy_mev, n_p=n_p, seed=seed,
                          sig_delta=0.0, eps_n_mm_mrad=eps_n, q_nc=q_nc,
                          betx_m=3.0, bety_m=3.0, label="sc_capstone")
    arr, man = cd.build(man)
    return arr, man


def track_xsuite_frozen(felsim_arr, segs, energy_mev, q_nc, ds_target=0.02,
                        spch_on=True):
    """Track the common distribution through `segs` in xsuite with frozen-Gaussian
    SC (or no SC). Mirrors the DA-FMM split-operator (half-map / SC / half-map)."""
    import xpart as xp
    import xtrack as xt
    import xfields as xf

    si = cd.to_physical_si(felsim_arr, energy_mev)
    g, b = cd.gamma_beta(energy_mev)
    bg = b * g
    n_p = len(si["x"])
    q_C = q_nc * 1e-9
    N_e = q_C / QE
    sig_x0 = float(si["x"].std()); sig_y0 = float(si["y"].std())
    sig_z0 = float(si["z"].std())

    p = xp.Particles(
        mass0=cd.E0_E * 1e6, q0=-1.0, kinetic_energy0=energy_mev * 1e6,
        x=si["x"], px=si["xp"], y=si["y"], py=si["yp"],
        zeta=si["z"], delta=si["delta"],
        weight=N_e / n_p,
    )

    long_profile = xf.LongitudinalProfileQGaussian(
        number_of_particles=N_e, sigma_z=sig_z0, z0=0.0, q_parameter=1.0)

    steps = scs._slice_schedule(segs, ds_target)
    elements = []
    for h, kx, ky in steps:
        # half-map
        if kx == 0.0:
            elements.append(xt.Drift(length=h / 2))
        else:
            elements.append(xt.Quadrupole(length=h / 2, k1=kx))
        if spch_on:
            elements.append(xf.SpaceChargeBiGaussian(
                length=h, longitudinal_profile=long_profile,
                sigma_x=sig_x0, sigma_y=sig_y0, mean_x=0.0, mean_y=0.0,
                update_on_track=True))
        if kx == 0.0:
            elements.append(xt.Drift(length=h / 2))
        else:
            elements.append(xt.Quadrupole(length=h / 2, k1=kx))

    line = xt.Line(elements=elements)
    line.particle_ref = xp.Particles(
        mass0=cd.E0_E * 1e6, q0=-1.0, kinetic_energy0=energy_mev * 1e6)
    line.build_tracker()
    t0 = time.perf_counter()
    line.track(p)
    wall = time.perf_counter() - t0

    alive = p.state > 0
    x, px = np.asarray(p.x[alive]), np.asarray(p.px[alive])
    y, py = np.asarray(p.y[alive]), np.asarray(p.py[alive])

    def epsn(u, up):
        du, dup = u - u.mean(), up - up.mean()
        e2 = np.mean(du**2) * np.mean(dup**2) - np.mean(du * dup)**2
        return bg * np.sqrt(max(e2, 0.0)) * 1e6
    return {
        "epsnx_exit": epsn(x, px), "epsny_exit": epsn(y, py),
        "sigx_exit_mm": float(x.std() * 1e3), "sigy_exit_mm": float(y.std() * 1e3),
        "n_alive": int(alive.sum()), "n_steps": len(steps), "wall": wall,
    }


def growth_pct(exit_val, init_val):
    return (exit_val - init_val) / init_val * 100.0


def run_point(energy_mev, q_nc, n_p, seed, eps_n, ds_target, optics_energy=None,
              softening="0", with_pic=False, pic_mesh=(64, 64, 64)):
    """softening: Plummer softening length [m] for the DA-FMM kernel, or 'auto' for
    the eps = sigma_x_init/sqrt(N_p) heuristic, or '0'/0 for the bare 1/r kernel.
    Softening suppresses the macroparticle-collision shot noise so the treecode
    approximates the smooth mean field (the alternative to brute-force high N_p)."""
    arr, man = make_beam(energy_mev, n_p, seed, eps_n, q_nc)
    segs = scs.extract_segments(LATTICE, *NODIP_RANGE, energy_mev,
                                optics_energy_mev=optics_energy)
    si = cd.to_physical_si(arr, energy_mev)
    if str(softening) == "auto":
        eps = float(si["x"].std()) / np.sqrt(n_p)
    else:
        eps = float(softening)
    # DA-FMM
    da_off = scs.track_dafmm(arr, segs, energy_mev, q_nc, ds_target, spch_on=False,
                             softening_eps=eps)
    da_on = scs.track_dafmm(arr, segs, energy_mev, q_nc, ds_target, spch_on=True,
                            softening_eps=eps)
    # xsuite frozen
    xs_off = track_xsuite_frozen(arr, segs, energy_mev, q_nc, ds_target, spch_on=False)
    xs_on = track_xsuite_frozen(arr, segs, energy_mev, q_nc, ds_target, spch_on=True)
    init_enx, init_eny = da_off["epsnx_init"], da_off["epsny_init"]
    out = {
        "energy_mev": energy_mev, "q_nc": q_nc, "n_p": n_p,
        "init_epsnx": init_enx, "init_epsny": init_eny,
        "dafmm_gx": growth_pct(da_on["epsnx_exit"], init_enx),
        "dafmm_gy": growth_pct(da_on["epsny_exit"], init_eny),
        "xsuite_gx": growth_pct(xs_on["epsnx_exit"], init_enx),
        "xsuite_gy": growth_pct(xs_on["epsny_exit"], init_eny),
        # SC-off sanity: both should be ~0
        "dafmm_off_gx": growth_pct(da_off["epsnx_exit"], init_enx),
        "xsuite_off_gx": growth_pct(xs_off["epsnx_exit"], init_enx),
        "xs_sigx_off": xs_off["sigx_exit_mm"], "da_sigx_off": da_off["sigx_exit_mm"],
        "wall_dafmm": da_on["wall_sc"], "wall_xsuite": xs_on["wall"],
        "seed": seed, "softening_eps": eps,
    }
    if with_pic:
        pic = track_pic(arr, energy_mev, q_nc, mesh=pic_mesh, optics_energy=optics_energy)
        if pic and "epsnx_exit" in pic:
            out["pic_gx"] = growth_pct(pic["epsnx_exit"], init_enx)
            out["pic_gy"] = growth_pct(pic["epsny_exit"], init_eny)
            out["wall_pic"] = pic["wall"]
        else:
            out["pic_gx"] = None
            out["pic_err"] = (pic or {}).get("error", "pic engine unavailable")
    return out


def _aggregate_seeds(per_seed):
    """Mean +/- std over seed realizations; used to average down DA-FMM shot noise."""
    out = dict(per_seed[0])
    keys = ["dafmm_gx", "dafmm_gy", "xsuite_gx", "xsuite_gy"]
    if all(r.get("pic_gx") is not None for r in per_seed):
        keys += ["pic_gx", "pic_gy"]
    elif "pic_gx" in out:
        out["pic_gx"] = None  # any seed missing pic -> mark the aggregate N/A
    for key in keys:
        vals = np.array([r[key] for r in per_seed])
        out[key] = float(vals.mean())
        out[key + "_std"] = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
    out["n_seeds"] = len(per_seed)
    out["dafmm_off_gx"] = max(abs(r["dafmm_off_gx"]) for r in per_seed)
    out["xsuite_off_gx"] = max(abs(r["xsuite_off_gx"]) for r in per_seed)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="No-dipole SC capstone: DA-FMM vs xsuite-frozen")
    ap.add_argument("--charges-nc", nargs="+", type=float, default=[0.1, 0.3, 1.0, 3.0])
    ap.add_argument("--energies-mev", nargs="+", type=float, default=[45.0])
    ap.add_argument("--n-p", type=int, default=6000)
    ap.add_argument("--seed", type=int, default=20260619)
    ap.add_argument("--eps-n", type=float, default=8.0)
    ap.add_argument("--ds-target", type=float, default=0.02)
    ap.add_argument("--optics-energy", type=float, default=None,
                    help="preserve the section optics (k0) at this energy across all "
                         "tracking energies; use 45 for the energy-scaling comparison "
                         "so a 1 MeV run isn't an over-focused blow-up (Risk R1)")
    ap.add_argument("--softening", default="0",
                    help="DA-FMM Plummer softening length [m], or 'auto' (sigma_x/sqrt(N_p)), "
                         "or '0' (bare 1/r). Softening suppresses the macroparticle-collision "
                         "shot noise -> the smooth mean-field comparison without huge N_p.")
    ap.add_argument("--seeds", nargs="+", type=int, default=None,
                    help="multiple seeds -> average DA-FMM growth over realizations "
                         "(mean +/- std), to drive down the single-seed shot-noise scatter. "
                         "Defaults to the single --seed.")
    ap.add_argument("--with-pic", action="store_true",
                    help="add the cosy-pic mesh-PIC engine as a 3rd column (needs "
                         "~/COSY/cosy-pic built; slower, spawns pic-core-cli per SC station)")
    ap.add_argument("--pic-mesh", nargs=3, type=int, default=[64, 64, 64],
                    help="cosy-pic Hockney mesh nx ny nz (default 64 64 64)")
    ap.add_argument("--outdir", type=Path, default=Path("test/results/sc_capstone"))
    args = ap.parse_args()

    seeds = args.seeds if args.seeds else [args.seed]
    pic_mesh = tuple(args.pic_mesh)
    rows = []
    for E in args.energies_mev:
        for q in args.charges_nc:
            print(f"[capstone] E={E} MeV  Q={q} nC  "
                  f"(soft={args.softening}, {len(seeds)} seed(s)"
                  f"{', +pic' if args.with_pic else ''}) ...", flush=True)
            per_seed = [run_point(E, q, args.n_p, s, args.eps_n, args.ds_target,
                                  optics_energy=args.optics_energy,
                                  softening=args.softening,
                                  with_pic=args.with_pic, pic_mesh=pic_mesh)
                        for s in seeds]
            r = per_seed[0] if len(seeds) == 1 else _aggregate_seeds(per_seed)
            rows.append(r)
            sd = (f" +/-{r['dafmm_gx_std']:.3f}" if "dafmm_gx_std" in r else "")
            picstr = ""
            if args.with_pic:
                picstr = (f"   cosy-pic x = {r['pic_gx']:+.3f}%"
                          if r.get("pic_gx") is not None
                          else f"   cosy-pic = N/A ({r.get('pic_err', '?')})")
            print(f"   DA-FMM  growth x = {r['dafmm_gx']:+.3f}%{sd}   "
                  f"xsuite-frozen x = {r['xsuite_gx']:+.3f}%{picstr}   "
                  f"(SC-off: DA {r['dafmm_off_gx']:+.2e}%, xs {r['xsuite_off_gx']:+.2e}%)")

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "sc_capstone_results.json").write_text(json.dumps(rows, indent=2))

    # Markdown table
    lines = ["# Reduced no-dipole SC capstone — DA-FMM vs xsuite-frozen",
             "",
             f"Section {NODIP_RANGE} (drift+quad, 6 quads, 1.65 m) of the FELsim line.",
             f"Common distribution: N_p={args.n_p}, eps_n={args.eps_n} mm.mrad, "
             f"sigma_delta=0, seed={args.seed}. Space charge: frozen-Gaussian (xsuite) "
             f"vs full N-body 1/r treecode (DA-FMM)"
             + (" vs Hockney mesh PIC (cosy-pic)." if args.with_pic else "."),
             "",
             ("| E [MeV] | Q [nC] | DA-FMM dEx | DA-FMM dEy | xsuite dEx | xsuite dEy | "
              + ("cosy-pic dEx | cosy-pic dEy | " if args.with_pic else "")
              + "DA-vs-xs (x) |"),
             "|---|---|---|---|---|---|" + ("---|---|" if args.with_pic else "") + "---|"]
    for r in rows:
        dvx = abs(r["dafmm_gx"] - r["xsuite_gx"])
        piccols = ""
        if args.with_pic:
            if r.get("pic_gx") is not None:
                piccols = f"{r['pic_gx']:+.3f}% | {r['pic_gy']:+.3f}% | "
            else:
                piccols = "N/A | N/A | "
        lines.append(f"| {r['energy_mev']:g} | {r['q_nc']:g} | "
                     f"{r['dafmm_gx']:+.3f}% | {r['dafmm_gy']:+.3f}% | "
                     f"{r['xsuite_gx']:+.3f}% | {r['xsuite_gy']:+.3f}% | "
                     f"{piccols}{dvx:.3f} pp |")
    lines += ["",
              "SC-off control (both engines should be ~0%): "
              f"max |DA-FMM off| = {max(abs(r['dafmm_off_gx']) for r in rows):.2e}%, "
              f"max |xsuite off| = {max(abs(r['xsuite_off_gx']) for r in rows):.2e}%.",
              ""]
    (args.outdir / "sc_capstone_table.md").write_text("\n".join(lines))

    # Plot: growth vs charge per energy (x-plane)
    fig, ax = plt.subplots(figsize=(7, 5))
    for E in args.energies_mev:
        sub = [r for r in rows if r["energy_mev"] == E]
        q = [r["q_nc"] for r in sub]
        ax.plot(q, [r["dafmm_gx"] for r in sub], "o-", label=f"DA-FMM {E:g} MeV")
        ax.plot(q, [r["xsuite_gx"] for r in sub], "s--", label=f"xsuite-frozen {E:g} MeV")
        if args.with_pic and all(r.get("pic_gx") is not None for r in sub):
            ax.plot(q, [r["pic_gx"] for r in sub], "^:", label=f"cosy-pic {E:g} MeV")
    ax.set_xlabel("bunch charge [nC]")
    ax.set_ylabel(r"$\epsilon_{n,x}$ growth [%]")
    ax.set_title(f"No-dipole section {NODIP_RANGE}: DA-FMM vs xsuite-frozen SC")
    ax.legend(); ax.grid(True, alpha=0.3)
    if max(r["q_nc"] for r in rows) / min(r["q_nc"] for r in rows) > 10:
        ax.set_xscale("log")
    fig.tight_layout()
    fig.savefig(args.outdir / "sc_capstone_growth.png", dpi=140)
    print(f"[capstone] wrote {args.outdir}/sc_capstone_table.md + .json + growth.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
