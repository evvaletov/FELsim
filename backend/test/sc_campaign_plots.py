#!/usr/bin/env python3
"""Space-charge campaign — Phase 1 result plots (2026-05-04 Niels items).

Reads existing L3 sweep data (COSY DA-FMM demo + FELsim three-code SC compare)
and produces a small, focused figure set:

  F1  Emittance growth vs equivalent macroparticle charge q_mp = Q/N_p   [item (c)]
  F2  Three-code SC charge sweep at 45 MeV (DA-FMM / xsuite-frozen / PIC)  [items (i),(l) preview]
  F3  N_p convergence of DA-FMM shot noise, with 1/sqrt(N_p) extrapolation [supports (c)]
  F4  N_slice integrator robustness (the N_slice=10 blind spot)            [numerical method]
  overview  2x2 montage of the above

No new tracking runs — pure re-analysis of committed sweep CSVs.

Author: Eremey Valetov
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COSY = Path("/home/evaletov/COSY/cosy-fmm/demo/spch_demo")
FELSIM_SC = Path(__file__).resolve().parent / "rftrack_linac" / "sc_compare_output"
OUT = Path(__file__).resolve().parent / "results" / "sc_campaign"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 150, "font.size": 11,
    "axes.grid": True, "grid.alpha": 0.3, "axes.axisbelow": True,
    "legend.framealpha": 0.9, "lines.markersize": 7,
})
C = {"dafmm": "#1f77b4", "frozen": "#2ca02c", "pic": "#d62728",
     "np": "#9467bd", "highnp": "#ff7f0e", "phys": "#555555"}


def macroparticle_charge_pC(Q_nC, N_p):
    """Equivalent charge per macroparticle, in pC."""
    return np.asarray(Q_nC) * 1e3 / np.asarray(N_p)


# ---------- load ----------
three = {k: pd.read_csv(FELSIM_SC / f"sweep_{n}.csv")
         for k, n in [("dafmm", "dafmm"), ("frozen", "xsuite-frozen"), ("pic", "xsuite-pic3d")]}
NP_FIXED = int(three["dafmm"]["alive"].iloc[0])           # 6000
# physical baseline = xsuite frozen (no macroparticle shot noise), in %
phys = three["frozen"].set_index("Q_nC")["epsny_growth"] * 100.0
phys_at1 = float(phys.loc[1.0])

np_extrap = pd.read_csv(COSY / "np_extrap_output" / "np_extrap_Q1.0nC.csv")   # Q=1nC, %
highnp = pd.read_csv(COSY / "highnp_charge_output" / "highnp_sweep_N20000.csv")  # N_p=2e4, %
N_HIGH = 20000
nslice = pd.read_csv(COSY / "sweeps" / "nslice" / "results.csv")


def agg(df, key, val):
    g = df.groupby(key)[val]
    return g.mean(), g.std().fillna(0.0), np.array(sorted(df[key].unique()))


# =====================================================================
# F1 — excess growth vs equivalent macroparticle charge  [item (c)]
# =====================================================================
def fig1():
    fig, ax = plt.subplots(figsize=(7.2, 5.2))

    # series A: charge sweep at N_p=6000 (vary Q) -> excess over xsuite at same Q
    QA = three["dafmm"]["Q_nC"].values
    exA = (three["dafmm"]["epsny_growth"].values - three["frozen"]["epsny_growth"].values) * 100
    qmpA = macroparticle_charge_pC(QA, NP_FIXED)

    # series B: N_p sweep at Q=1nC (vary N_p) -> excess over xsuite at Q=1
    mB, sB, npsB = agg(np_extrap, "n_p", "epsny_growth_pct")
    exB = mB.values - phys_at1
    qmpB = macroparticle_charge_pC(1.0, npsB)

    # series C: charge sweep at N_p=20000 (vary Q) -> excess over xsuite at same Q
    mC, sC, qsC = agg(highnp, "q_nc", "epsny_growth_pct")
    physC = np.interp(qsC, phys.index.values, phys.values)
    exC = mC.values - physC
    qmpC = macroparticle_charge_pC(qsC, N_HIGH)

    for q, e, c, m, lab in [
        (qmpA, exA, C["dafmm"], "o", f"charge sweep, $N_p$={NP_FIXED}"),
        (qmpB, exB, C["np"], "^", "$N_p$ sweep, $Q$=1 nC"),
        (qmpC, exC, C["highnp"], "s", f"charge sweep, $N_p$={N_HIGH}")]:
        pos = e > 0
        ax.plot(np.asarray(q)[pos], np.asarray(e)[pos], m, color=c, label=lab,
                mec="k", mew=0.4, ls="none")

    # power-law fit on the clean fixed-Q (series B) excess
    pos = exB > 0
    slope, intc = np.polyfit(np.log10(qmpB[pos]), np.log10(exB[pos]), 1)
    xx = np.logspace(np.log10(qmpB.min()) - 0.2, np.log10(qmpA.max()) + 0.2, 50)
    yy = 10 ** (intc + slope * np.log10(xx))
    ax.plot(xx, yy, "--", color=C["np"], lw=1.4, alpha=0.8,
            label=fr"fit $\propto q_{{mp}}^{{{slope:.2f}}}$")

    # threshold: numerical excess == physical mean-field (at Q=1 nC)
    q_star = 10 ** ((np.log10(phys_at1) - intc) / slope)
    ax.axhline(phys_at1, color=C["phys"], ls=":", lw=1.2,
               label=f"physical SC ($Q$=1 nC): {phys_at1:.2f}%")
    ax.axvline(q_star, color="k", ls="-.", lw=1.0, alpha=0.6)
    ax.annotate(f"$q_{{mp}}^*$ ≈ {q_star:.3f} pC\n($N_p$≈{1e3/q_star:.0f} at 1 nC)",
                xy=(q_star, phys_at1), xytext=(q_star * 1.2, phys_at1 * 4),
                fontsize=9, arrowprops=dict(arrowstyle="->", color="k", alpha=0.6))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("equivalent macroparticle charge  $q_{mp}=Q/N_p$  [pC]")
    ax.set_ylabel(r"excess $\Delta\varepsilon_{n,y}$ growth  (DA-FMM $-$ xsuite)  [%]")
    ax.set_title("Numerical space-charge noise scales with macroparticle charge")
    ax.legend(fontsize=8.5, loc="lower right")
    fig.tight_layout()
    save(fig, "F1_macroparticle_charge")
    return dict(slope=slope, q_star_pC=q_star, phys_at1_pct=phys_at1)


# =====================================================================
# F2 — three-code charge sweep at 45 MeV  [items (i),(l) preview]
# =====================================================================
def fig2():
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    style = [("frozen", "xsuite frozen-Gaussian", "o-", C["frozen"]),
             ("pic", "xsuite 3D PIC", "s--", C["pic"]),
             ("dafmm", f"COSY DA-FMM ($N_p$={NP_FIXED})", "^-", C["dafmm"])]
    for k, lab, ls, c in style:
        d = three[k]
        ax.plot(d["Q_nC"], d["epsny_growth"] * 100, ls, color=c, label=lab,
                mec="k", mew=0.4)
    ax.set_xscale("log")
    ax.set_xlabel("bunch charge  $Q$  [nC]")
    ax.set_ylabel(r"$\Delta\varepsilon_{n,y}$ growth  [%]")
    ax.set_title("Three-code space-charge agreement — 2 m FODO, 45 MeV")
    ax.legend(fontsize=9, loc="upper left")
    ax.annotate("DA-FMM offset = macroparticle\nshot noise (see F1/F3)",
                xy=(1.0, float(three['dafmm'].set_index('Q_nC')['epsny_growth'].loc[1.0] * 100)),
                xytext=(0.15, 2.0), fontsize=8.5,
                arrowprops=dict(arrowstyle="->", color=C["dafmm"], alpha=0.7))
    fig.tight_layout()
    save(fig, "F2_three_code_charge")


# =====================================================================
# F3 — N_p convergence / shot-noise extrapolation  [supports (c)]
# =====================================================================
def fig3():
    m, s, nps = agg(np_extrap, "n_p", "epsny_growth_pct")
    x = nps.astype(float)
    y = m.reindex(nps).values
    e = s.reindex(nps).values
    # fit growth = a + b/sqrt(N_p)
    A = np.vstack([np.ones_like(x), 1 / np.sqrt(x)]).T
    (a, b), *_ = np.linalg.lstsq(A, y, rcond=None)

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.errorbar(x, y, yerr=e, fmt="o", color=C["dafmm"], mec="k", mew=0.4,
                capsize=3, label="DA-FMM (3 seeds, $Q$=1 nC)")
    xx = np.linspace(x.min(), x.max(), 200)
    ax.plot(xx, a + b / np.sqrt(xx), "--", color=C["dafmm"],
            label=fr"fit $\Delta\varepsilon = {a:.2f} + {b:.0f}/\sqrt{{N_p}}$  [%]")
    ax.axhline(a, color=C["np"], ls="-.", lw=1.2,
               label=fr"$N_p\!\to\!\infty$ limit: {a:.2f}%")
    ax.axhline(phys_at1, color=C["frozen"], ls=":", lw=1.4,
               label=f"xsuite physical: {phys_at1:.2f}%")
    ax.set_xlabel("macroparticles  $N_p$")
    ax.set_ylabel(r"$\Delta\varepsilon_{n,y}$ growth  [%]")
    ax.set_title(r"DA-FMM shot noise $\sim 1/\sqrt{N_p}$ converges to physical SC")
    ax.legend(fontsize=9)
    fig.tight_layout()
    save(fig, "F3_np_convergence")
    return dict(np_inf_pct=float(a), shot_coeff=float(b))


# =====================================================================
# F4 — N_slice integrator robustness
# =====================================================================
def fig4():
    d = nslice.sort_values("N_slice")
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(d["N_slice"], d["epsny_growth"] * 100, "o-", color=C["dafmm"],
            mec="k", mew=0.4, label=r"$\varepsilon_{n,y}$")
    ax.plot(d["N_slice"], d["epsnx_growth"] * 100, "s--", color=C["highnp"],
            mec="k", mew=0.4, label=r"$\varepsilon_{n,x}$")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("longitudinal slices  $N_{slice}$  (split-operator steps)")
    ax.set_ylabel(r"$\Delta\varepsilon_n$ growth  [%]")
    ax.set_title("Integrator convergence — plateau for $N_{slice}\\gtrsim 80$")
    # flag the N_slice=10 anomaly if present
    if (d["N_slice"] == 10).any():
        row = d[d["N_slice"] == 10].iloc[0]
        ax.annotate("$N_{slice}$=10 blind spot",
                    xy=(10, row["epsny_growth"] * 100), xytext=(20, row["epsny_growth"] * 100 + 0.4),
                    fontsize=9, arrowprops=dict(arrowstyle="->", color="k", alpha=0.6))
    ax.legend(fontsize=9)
    fig.tight_layout()
    save(fig, "F4_nslice_convergence")


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}.png / .pdf")


def overview():
    """2x2 montage by re-importing the saved PNGs."""
    import matplotlib.image as mpimg
    names = ["F1_macroparticle_charge", "F2_three_code_charge",
             "F3_np_convergence", "F4_nslice_convergence"]
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.5))
    for ax, n in zip(axs.ravel(), names):
        ax.imshow(mpimg.imread(OUT / f"{n}.png")); ax.axis("off")
    fig.suptitle("Space-charge campaign — Phase 1 (2 m FODO, 45 MeV; L3 data re-analysis)",
                 fontsize=13, y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "overview.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  wrote overview.png")


if __name__ == "__main__":
    print(f"output -> {OUT}")
    r1 = fig1()
    fig2()
    r3 = fig3()
    fig4()
    overview()
    print("\nKey numbers:")
    print(f"  physical SC (xsuite, Q=1nC):  {r1['phys_at1_pct']:.3f} %")
    print(f"  excess-vs-q_mp power-law slope: {r1['slope']:.3f}")
    print(f"  threshold q_mp* (excess=physical): {r1['q_star_pC']:.4f} pC  "
          f"(N_p≈{1e3/r1['q_star_pC']:.0f} at 1 nC)")
    print(f"  DA-FMM N_p->inf limit:        {r3['np_inf_pct']:.3f} %  "
          f"(shot coeff {r3['shot_coeff']:.1f})")
