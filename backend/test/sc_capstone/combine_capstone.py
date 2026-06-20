"""Combine the 45 MeV and 1 MeV capstone sweeps into one report table + figure."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("test/results/sc_capstone")


def load(p):
    f = OUT / p / "sc_capstone_results.json"
    return json.loads(f.read_text()) if f.exists() else []


def main():
    rows = load("e45") + load("e1")
    rows.sort(key=lambda r: (-r["energy_mev"], r["q_nc"]))

    lines = [
        "## Reduced no-dipole SC capstone — DA-FMM vs xsuite-frozen (combined)",
        "",
        "Section [32,46) (drift+quad, 6 quads, 1.65 m). Common distribution N_p=6000,",
        "eps_n=8 mm.mrad, sigma_delta=0. Optics held fixed (k0 at 45 MeV) so the 1 MeV",
        "rows isolate space-charge scaling, not a focusing-instability blow-up.",
        "SC: xsuite frozen-Gaussian (nonlinear Bassetti-Erskine, rms-self-consistent) "
        "vs DA-FMM N-body 1/r treecode (resolves the actual profile).",
        "",
        "| E [MeV] | Q [nC] | DA-FMM dEx | xsuite dEx | DA/xs ratio |",
        "|---|---|---|---|---|",
    ]
    for r in rows:
        ratio = (r["dafmm_gx"] / r["xsuite_gx"]) if abs(r["xsuite_gx"]) > 1e-6 else float("nan")
        lines.append(f"| {r['energy_mev']:g} | {r['q_nc']:g} | {r['dafmm_gx']:+.3f}% | "
                     f"{r['xsuite_gx']:+.3f}% | {ratio:.2f} |")
    offmax = max(abs(r["dafmm_off_gx"]) for r in rows)
    xoffmax = max(abs(r["xsuite_off_gx"]) for r in rows)
    lines += ["",
              f"SC-off control (both ~0): max |DA-FMM|={offmax:.1e}%, max |xsuite|={xoffmax:.1e}%.",
              ""]
    (OUT / "capstone_combined_table.md").write_text("\n".join(lines))

    # 2-panel figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, E in zip(axes, (45.0, 1.0)):
        sub = sorted([r for r in rows if r["energy_mev"] == E], key=lambda r: r["q_nc"])
        if not sub:
            continue
        q = [r["q_nc"] for r in sub]
        ax.plot(q, [r["dafmm_gx"] for r in sub], "o-", label="DA-FMM (N-body)")
        ax.plot(q, [r["xsuite_gx"] for r in sub], "s--", label="xsuite frozen-Gaussian")
        ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=0.01)
        ax.set_xlabel("bunch charge [nC]")
        ax.set_ylabel(r"$\epsilon_{n,x}$ growth [%]")
        ax.set_title(f"{E:g} MeV")
        ax.legend(); ax.grid(True, which="both", alpha=0.3)
    fig.suptitle("No-dipole section [32,46): cosy-fmm DA-FMM vs xsuite-frozen space charge")
    fig.tight_layout()
    fig.savefig(OUT / "capstone_combined_growth.png", dpi=140)
    print(f"wrote {OUT}/capstone_combined_table.md + capstone_combined_growth.png")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
