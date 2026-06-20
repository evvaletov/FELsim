"""
No-SC cross-code handoff harness (capstone prerequisite + cold-beam deliverable).

Tracks ONE common macroparticle distribution (from common_distribution) through a
chosen element range of the FELsim transport line with space charge OFF, in each
available code (FELsim / xsuite / RF-Track), and reports how the second moments
(sigma_x, sigma_y, eps_n_x, eps_n_y) and centroid agree at the exit. This isolates
code-physics differences (e.g. xsuite's missing dipole edge/fringe) from any later
space-charge disagreement: at zero current the codes MUST agree on linear transport,
so this table sets the "cold-beam" tolerance floor for the SC capstone.

Two standard cases:
  * no-dipole subset (elements [32, 46) -- the longest drift+quad run, 6 quads,
    1.65 m): all codes should agree to tight tolerance.
  * full line ([0, 137)): xsuite has no dipole edge/fringe and RF-Track has its
    own bend handling, so divergence here is expected and is the point.

FELsim (analytic first-order maps) is the reference column.

Author: Eremey Valetov
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from multiCodeSimulator import MultiCodeSimulator, SimSection
import common_distribution as cd

LATTICE = str(cd.FELSIM_REPO / "var" / "UH_FEL_beamline.json")
NODIP_RANGE = (32, 46)     # longest no-dipole drift+quad section
FULL_RANGE = (0, 137)

E0_E = cd.E0_E


def moments(particles: np.ndarray, energy_mev: float) -> dict:
    """Second-moment beam parameters from a FELsim-coordinate array
    [x mm, x' mrad, y mm, y' mrad, ...]."""
    g, b = cd.gamma_beta(energy_mev)
    bg = b * g
    x, xp, y, yp = (particles[:, 0], particles[:, 1],
                    particles[:, 2], particles[:, 3])
    dx, dxp = x - x.mean(), xp - xp.mean()
    dy, dyp = y - y.mean(), yp - yp.mean()
    eps_x = np.sqrt(max(np.mean(dx**2) * np.mean(dxp**2) - np.mean(dx * dxp)**2, 0.0))
    eps_y = np.sqrt(max(np.mean(dy**2) * np.mean(dyp**2) - np.mean(dy * dyp)**2, 0.0))
    return {
        "n": int(particles.shape[0]),
        "cx_mm": float(x.mean()), "cy_mm": float(y.mean()),
        "sigx_mm": float(x.std()), "sigy_mm": float(y.std()),
        "epsnx": float(bg * eps_x), "epsny": float(bg * eps_y),
    }


def run_code(code: str, particles: np.ndarray, elem_range: tuple,
             energy_mev: float, space_charge: bool = False) -> dict:
    """Track `particles` through `elem_range` in one code; return exit moments."""
    cfg = {"space_charge": space_charge} if space_charge else {}
    mc = MultiCodeSimulator(
        sections=[SimSection(name=code, simulator_key=code,
                             element_range=elem_range, config=cfg)],
        lattice_path=LATTICE, beam_energy=energy_mev,
    )
    t0 = time.perf_counter()
    res = mc.simulate(particles=particles.copy())
    wall = time.perf_counter() - t0
    if not res.success or res.final_particles is None:
        return {"ok": False, "wall_s": wall,
                "err": res.metadata.get("reason", res.metadata)}
    m = moments(res.final_particles, energy_mev)
    m["ok"] = True
    m["wall_s"] = wall
    return m


def _rel(a, b):
    """Relative difference |a-b|/max(|b|, floor)."""
    denom = max(abs(b), 1e-12)
    return abs(a - b) / denom


def compare(particles, elem_range, energy_mev, codes, space_charge=False):
    init = moments(particles, energy_mev)
    out = {"init": init, "codes": {}}
    for c in codes:
        try:
            out["codes"][c] = run_code(c, particles, elem_range, energy_mev, space_charge)
        except Exception as e:
            out["codes"][c] = {"ok": False, "err": f"{type(e).__name__}: {e}"}
    return out


def render_table(result: dict, ref: str, label: str, tol: dict) -> str:
    init = result["init"]
    codes = result["codes"]
    keys = ["sigx_mm", "sigy_mm", "epsnx", "epsny", "cx_mm", "cy_mm"]
    names = {"sigx_mm": "sig_x [mm]", "sigy_mm": "sig_y [mm]",
             "epsnx": "eps_n,x [mm.mrad]", "epsny": "eps_n,y [mm.mrad]",
             "cx_mm": "centroid_x [mm]", "cy_mm": "centroid_y [mm]"}
    lines = [f"### {label}", ""]
    ok_codes = [c for c in codes if codes[c].get("ok")]
    hdr = "| quantity | initial | " + " | ".join(ok_codes) + " |"
    sep = "|" + "---|" * (len(ok_codes) + 2)
    lines += [hdr, sep]
    for k in keys:
        row = [names[k], f"{init[k]:.5g}"]
        for c in ok_codes:
            row.append(f"{codes[c][k]:.5g}")
        lines.append("| " + " | ".join(row) + " |")
    # agreement vs reference
    lines.append("")
    if ref in codes and codes[ref].get("ok"):
        others = [c for c in ok_codes if c != ref]
        if others:
            lines.append(f"**Agreement vs {ref} (relative):**")
            lines.append("")
            ahdr = "| quantity | " + " | ".join(others) + " |"
            asep = "|" + "---|" * (len(others) + 1)
            lines += [ahdr, asep]
            worst = 0.0
            for k in keys[:4]:  # sigmas + emittances (centroids ~0, skip rel)
                row = [names[k]]
                for c in others:
                    r = _rel(codes[c][k], codes[ref][k])
                    worst = max(worst, r)
                    row.append(f"{r*100:.3f}%")
                lines.append("| " + " | ".join(row) + " |")
            band = tol.get("sigma_emit", 0.02)
            verdict = "PASS" if worst <= band else "EXCEEDS"
            lines.append("")
            lines.append(f"Worst sigma/emittance disagreement vs {ref}: "
                         f"**{worst*100:.3f}%** (tolerance {band*100:.0f}% -> {verdict})")
    # failures
    failed = [(c, codes[c].get("err")) for c in codes if not codes[c].get("ok")]
    if failed:
        lines.append("")
        lines.append("**Codes that did not complete:**")
        for c, err in failed:
            lines.append(f"- `{c}`: {err}")
    # wall clock
    walls = [(c, codes[c].get("wall_s")) for c in ok_codes if codes[c].get("wall_s")]
    if walls:
        lines.append("")
        lines.append("Wall clock: " + ", ".join(f"{c} {w:.2f}s" for c, w in walls))
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description="No-SC cross-code handoff / cold-beam table")
    ap.add_argument("--manifest", type=Path,
                    default=Path("test/results/sc_capstone/dist_default.json"))
    ap.add_argument("--codes", nargs="+", default=["felsim", "xsuite", "rftrack"])
    ap.add_argument("--ref", default="felsim")
    ap.add_argument("--energy-mev", type=float, default=45.0)
    ap.add_argument("--cases", nargs="+", default=["nodip", "full"],
                    choices=["nodip", "full"])
    ap.add_argument("--out", type=Path,
                    default=Path("test/results/sc_capstone/nosc_agreement.md"))
    args = ap.parse_args()

    if args.manifest.exists():
        man = cd.load_manifest(args.manifest)
        man.energy_mev = args.energy_mev
        particles = cd.make_felsim_distribution(man)
        sha = cd.array_sha256(particles)
    else:
        man = cd.BeamManifest(energy_mev=args.energy_mev)
        particles, man = cd.build(man)
        sha = man.array_sha256

    tol = {"sigma_emit": 0.02}
    ranges = {"nodip": NODIP_RANGE, "full": FULL_RANGE}
    labels = {"nodip": f"No-dipole subset {NODIP_RANGE} (drift+quad, cold beam)",
              "full": f"Full transport line {FULL_RANGE} (cold beam; xsuite has no "
                      "dipole edge/fringe -> divergence expected)"}

    doc = [
        f"# No-SC cross-code handoff — cold-beam agreement",
        "",
        f"- Distribution: N_p={man.n_p}, E={args.energy_mev} MeV, eps_n="
        f"{man.eps_n_mm_mrad} mm.mrad, seed={man.seed}, sha256={sha[:16]}",
        f"- Codes: {', '.join(args.codes)} (reference: {args.ref})",
        f"- Space charge: OFF",
        "",
    ]
    for case in args.cases:
        print(f"[nosc] case={case} range={ranges[case]} codes={args.codes}")
        res = compare(particles, ranges[case], args.energy_mev, args.codes, space_charge=False)
        doc.append(render_table(res, args.ref, labels[case], tol))
        # console summary
        for c, m in res["codes"].items():
            if m.get("ok"):
                print(f"   {c:8s} sigx={m['sigx_mm']:.4f} sigy={m['sigy_mm']:.4f} "
                      f"epsnx={m['epsnx']:.4f} epsny={m['epsny']:.4f} ({m['wall_s']:.2f}s)")
            else:
                print(f"   {c:8s} FAILED: {m.get('err')}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(doc))
    print(f"[nosc] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
