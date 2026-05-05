"""
L3.6 -- four-engine space-charge comparison on the spch_demo lattice.

Same 2 m FODO + same Gaussian bunch as spch_demo, swept across charge.
DA-FMM results pulled from the existing L3.4 CSV (no need to re-run).

Engines:
  * COSY DA-FMM           (existing CSV; spch_demo)
  * xsuite frozen-Gaussian (xfields.SpaceChargeBiGaussian, linearised analytical)
  * xsuite 3D PIC         (xfields.SpaceCharge3D, full nonlinear)
  * RF-Track PIC          (rft.SpaceCharge_PIC_FreeSpace, full nonlinear)

Lattice: 0.5 m drift, thin QF (1/f=2), 1.0 m drift, thin QD (1/f=2), 0.5 m drift.
Bunch:   45 MeV electrons, sigma_xy=1 mm, sigma_xpyp=0.1 mrad, sigma_z=0.6 mm,
         sigma_delta=5e-3, fixed seed, N_p=6000.
"""
from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import RF_Track as rft
import xfields as xf
import xobjects as xo
import xpart as xp
import xtrack as xt

# Constants
MC2_MEV = 0.510998950
QE = 1.602176634e-19

# Lattice
L_DRIFT_END = 0.5
L_DRIFT_MID = 1.0
L_TOT = 2.0
S_QF = 0.5
S_QD = 1.5
K_FOCUS = 2.0  # 1/f in 1/m

# Bunch
KE_MEV = 45.0
SIG_X = SIG_Y = 1.0e-3        # m
SIG_XP = SIG_YP = 1.0e-4      # rad
SIG_Z = 6.0e-4                # m
# delta = 0 deliberately: spch_demo's HALFDRIFT/THINQUAD are non-chromatic
# (no (1+delta) factors), so its transport ignores momentum spread. xsuite
# and RF-Track are rigorously chromatic; running them with sigma_delta>0
# adds ~17% spurious emittance growth that has no counterpart in spch_demo.
# Setting delta=0 isolates the SC contribution across all four codes.
SIG_DELTA = 0.0               # dp/p (see note above)
N_P = 6000
SEED = 20260420

# Sweep
CHARGE_GRID_NC = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0]

# Slicing
N_SLICE_SC = 80               # matches L3.2 N_slice*

WORK_DIR = Path(__file__).parent
OUT_DIR = WORK_DIR / "sc_compare_output"
DAFMM_CSV = Path(
    "/home/evaletov/COSY/cosy-fmm/demo/spch_demo/sweeps/charge/results.csv"
)


# ----------------------------------------------------------------------
# Bunch generation (shared)
# ----------------------------------------------------------------------


def make_bunch(seed: int = SEED) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    return {
        "x":  SIG_X  * rng.standard_normal(N_P),
        "xp": SIG_XP * rng.standard_normal(N_P),
        "y":  SIG_Y  * rng.standard_normal(N_P),
        "yp": SIG_YP * rng.standard_normal(N_P),
        "z":  SIG_Z  * rng.standard_normal(N_P),
        "delta": SIG_DELTA * rng.standard_normal(N_P),
    }


def gamma_beta() -> tuple[float, float]:
    gamma = 1.0 + KE_MEV / MC2_MEV
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    return gamma, beta


def moments_from_arrays(x, y, z) -> tuple[float, float, float]:
    return float(np.std(x)), float(np.std(y)), float(np.std(z))


def emittance_n(x, xp, betagamma) -> float:
    sx2 = float(np.mean(x * x))
    sxp2 = float(np.mean(xp * xp))
    sxxp = float(np.mean(x * xp))
    eps2 = sx2 * sxp2 - sxxp * sxxp
    if eps2 < 0:
        eps2 = 0.0
    return betagamma * np.sqrt(eps2)


# ----------------------------------------------------------------------
# Engine: xsuite (frozen-Gaussian or 3D PIC)
# ----------------------------------------------------------------------


def _xsuite_line_with_sc(
    bunch: dict[str, np.ndarray],
    Q_C: float,
    sc_kind: str,
    n_slice: int = N_SLICE_SC,
) -> tuple[xt.Line, xp.Particles]:
    """Build an xsuite line: drift sliced into n_slice pieces with SC inserted
    between halves (split-operator), thin quads at s=L/4 and s=3L/4 of total.
    Both `frozen` and `pic3d` use the same lattice, just different SC element.
    """
    gamma, beta = gamma_beta()
    bg = beta * gamma

    # Number of real electrons in the bunch
    N_e = Q_C / QE

    # Build particles
    p = xp.Particles(
        mass0=MC2_MEV * 1e6,
        q0=-1.0,
        kinetic_energy0=KE_MEV * 1e6,
        x=bunch["x"], px=bunch["xp"],
        y=bunch["y"], py=bunch["yp"],
        zeta=bunch["z"],
        delta=bunch["delta"],
        weight=N_e / N_P,  # real electrons per macroparticle
    )

    # SC element factory — fresh element per slice with the correct length
    ds = L_TOT / n_slice
    long_profile = xf.LongitudinalProfileQGaussian(
        number_of_particles=N_e,
        sigma_z=SIG_Z,
        z0=0.0,
        q_parameter=1.0,  # Gaussian
    )

    def make_sc():
        if sc_kind == "frozen":
            return xf.SpaceChargeBiGaussian(
                length=ds,
                longitudinal_profile=long_profile,
                sigma_x=SIG_X,
                sigma_y=SIG_Y,
                mean_x=0.0, mean_y=0.0,
                update_on_track=True,
            )
        if sc_kind == "pic3d":
            return xf.SpaceCharge3D(
                length=ds,
                update_on_track=True,
                x_range=(-8 * SIG_X, 8 * SIG_X),
                y_range=(-8 * SIG_Y, 8 * SIG_Y),
                z_range=(-8 * SIG_Z, 8 * SIG_Z),
                nx=64, ny=64, nz=64,
                solver="FFTSolver2p5D",
                gamma0=gamma,
            )
        raise ValueError(sc_kind)

    # Build the slice-by-slice element list, dropping the QF/QD at the right s
    elements = []
    s = 0.0
    for i in range(n_slice):
        elements.append(xt.Drift(length=ds / 2))
        elements.append(make_sc())
        elements.append(xt.Drift(length=ds / 2))
        s_new = (i + 1) * ds
        # Insert thin quads at s_new closest to the design positions
        # spch_demo uses the (s>L_q-HDS)*(s<L_q+HDS) window test; here we apply
        # at the slice that strictly bounds the quad position.
        if s < S_QF <= s_new:
            elements.append(xt.Multipole(knl=[0.0, +K_FOCUS]))  # QF
        if s < S_QD <= s_new:
            elements.append(xt.Multipole(knl=[0.0, -K_FOCUS]))  # QD
        s = s_new

    line = xt.Line(elements=elements)
    line.particle_ref = xp.Particles(
        mass0=MC2_MEV * 1e6, q0=-1.0, kinetic_energy0=KE_MEV * 1e6
    )
    line.build_tracker()
    return line, p


def run_xsuite(bunch, Q_C, sc_kind):
    line, p = _xsuite_line_with_sc(bunch, Q_C, sc_kind)
    t0 = time.perf_counter()
    line.track(p)
    wall = time.perf_counter() - t0
    gamma, beta = gamma_beta()
    bg = beta * gamma

    # Get arrays back. xsuite particles may have lost some — check state.
    alive = p.state > 0
    x = np.asarray(p.x[alive])
    px = np.asarray(p.px[alive])
    y = np.asarray(p.y[alive])
    py = np.asarray(p.py[alive])
    zeta = np.asarray(p.zeta[alive])

    sigx, sigy, sigz = moments_from_arrays(x, y, zeta)
    epsnx = emittance_n(x, px, bg)
    epsny = emittance_n(y, py, bg)
    return dict(sigx=sigx, sigy=sigy, sigz=sigz,
                epsnx=epsnx, epsny=epsny, wall=wall, alive=int(alive.sum()))


# ----------------------------------------------------------------------
# Engine: RF-Track PIC
# ----------------------------------------------------------------------


def run_rftrack(bunch, Q_C):
    """RF-Track PIC SC adapter (L3.6 Phase 2).

    Correct Volume pattern: explicit per-element add(elem, x, y, z) with
    set_s0/set_s1 boundaries; SC_engine set globally via rft.cvar.SC_engine.
    The earlier vol.add(lattice, ...) attempt core-dumped on RF-Track 2.5.5.
    """
    gamma, beta = gamma_beta()
    bg = beta * gamma
    P_ref_MeVc = beta * gamma * MC2_MEV  # MeV/c

    n_part = N_P
    P_per = P_ref_MeVc * (1.0 + bunch["delta"])

    # RF-Track Bunch6d state: [x[mm] xp[mrad] y[mm] yp[mrad] t[mm/c] P[MeV/c]]
    phase_space = np.column_stack([
        bunch["x"]  * 1e3,
        bunch["xp"] * 1e3,
        bunch["y"]  * 1e3,
        bunch["yp"] * 1e3,
        -bunch["z"] * 1e3,
        P_per,
    ])
    n_e_per_macro = (Q_C / QE) / n_part
    bunch_in = rft.Bunch6d(MC2_MEV, n_e_per_macro, -1.0, phase_space)

    # Set the PIC SC engine globally before building the volume
    sc = rft.SpaceCharge_PIC_FreeSpace(64, 64, 64)
    rft.cvar.SC_engine = sc

    # Build a Volume bounded by [0, L_TOT] with two thin Multipoles at
    # the QF / QD positions. No drifts: Volume time-integrates between
    # elements, drifts are implicit.
    vol = rft.Volume()
    vol.set_s0(0.0)
    vol.set_s1(L_TOT)
    vol.add(rft.Multipole(0.0, [0.0, +K_FOCUS]), 0.0, 0.0, S_QF,
            reference="center")
    vol.add(rft.Multipole(0.0, [0.0, -K_FOCUS]), 0.0, 0.0, S_QD,
            reference="center")
    # Integration time step in mm/c. For beta=0.94 over 2 m, total time
    # is ~2128 mm/c. dt_mm=25 gives ~85 steps, comparable to N_slice=80.
    vol.dt_mm = 25.0
    vol.odeint_algorithm = "rk2"

    line = rft.Lattice()
    line.append(vol)

    t0 = time.perf_counter()
    bunch_out = line.track(bunch_in)
    wall = time.perf_counter() - t0

    M = bunch_out.get_phase_space("%x %xp %y %yp %t %Pc")
    if M.shape[0] == 0:
        return dict(sigx=np.nan, sigy=np.nan, sigz=np.nan,
                    epsnx=np.nan, epsny=np.nan, wall=wall, alive=0)

    x = M[:, 0] * 1e-3
    xp = M[:, 1] * 1e-3
    y = M[:, 2] * 1e-3
    yp = M[:, 3] * 1e-3
    z = -M[:, 4] * 1e-3

    sigx, sigy, sigz = moments_from_arrays(x, y, z)
    epsnx = emittance_n(x, xp, bg)
    epsny = emittance_n(y, yp, bg)
    return dict(sigx=sigx, sigy=sigy, sigz=sigz,
                epsnx=epsnx, epsny=epsny, wall=wall, alive=int(M.shape[0]))


# ----------------------------------------------------------------------
# Sweep + report
# ----------------------------------------------------------------------


def load_dafmm() -> dict[float, dict[str, float]]:
    """Pull the L3.4 charge-sweep CSV: Q[nC] -> dict of exit moments + the
    CSV-recorded growth (which is computed against spch_demo's *own* initial
    emittance, NOT the numpy bunch's). Mixing the two would inflate dafmm
    growth by ~4pp because the FOX Box-Muller bunch and the numpy bunch
    differ at the few-percent level."""
    out: dict[float, dict[str, float]] = {}
    with DAFMM_CSV.open() as f:
        r = csv.DictReader(f)
        for row in r:
            Q = float(row["Q_bunch [nC]"])
            out[round(Q, 3)] = dict(
                sigx=float(row["sigx_exit_m"]),
                sigy=float(row["sigy_exit_m"]),
                sigz=float(row["sigz_exit_m"]),
                epsnx=float(row["epsnx_exit"]),
                epsny=float(row["epsny_exit"]),
                # Growth columns come straight from the L3.4 CSV (computed
                # in spch_demo against its own initial bunch).
                epsnx_growth=float(row["epsnx_growth"]),
                epsny_growth=float(row["epsny_growth"]),
                wall=float(row["wall_on_s"]),
                alive=N_P,
            )
    return out


def initial_emittance(bunch) -> tuple[float, float]:
    gamma, beta = gamma_beta()
    bg = beta * gamma
    return emittance_n(bunch["x"], bunch["xp"], bg), emittance_n(bunch["y"], bunch["yp"], bg)


def sweep(engines: list[str]) -> dict[str, list[dict]]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    bunch = make_bunch()
    eps_init = initial_emittance(bunch)
    print(f"Initial bunch eps_nx={eps_init[0]:.4e} eps_ny={eps_init[1]:.4e}")

    results: dict[str, list[dict]] = {e: [] for e in engines}

    for Q_nC in CHARGE_GRID_NC:
        Q_C = Q_nC * 1e-9
        print(f"\n=== Q = {Q_nC} nC ===")
        if "dafmm" in engines:
            d = load_dafmm()[round(Q_nC, 3)]
            d["Q_nC"] = Q_nC
            # Keep the CSV's growth values (computed against spch_demo's own
            # initial bunch, not the numpy/xsuite bunch).
            results["dafmm"].append(d)
            print(f"  dafmm    : epsny_growth={d['epsny_growth']*100:7.3g}% (cached)")

        for kind, label in [("frozen", "xsuite-frozen"), ("pic3d", "xsuite-pic3d")]:
            if kind not in engines and label not in engines:
                continue
            try:
                d = run_xsuite(bunch, Q_C, kind)
                d["Q_nC"] = Q_nC
                d["epsnx_growth"] = (d["epsnx"] - eps_init[0]) / eps_init[0]
                d["epsny_growth"] = (d["epsny"] - eps_init[1]) / eps_init[1]
                results[label].append(d)
                print(f"  {label:14s}: epsny_growth={d['epsny_growth']*100:7.3g}% "
                      f"alive={d['alive']:4d}/{N_P} wall={d['wall']:.1f}s")
            except Exception as e:
                print(f"  {label:14s}: FAILED {type(e).__name__}: {e}")
                results[label].append(dict(Q_nC=Q_nC, error=str(e)))

        if "rftrack" in engines:
            try:
                d = run_rftrack(bunch, Q_C)
                d["Q_nC"] = Q_nC
                d["epsnx_growth"] = (d["epsnx"] - eps_init[0]) / eps_init[0]
                d["epsny_growth"] = (d["epsny"] - eps_init[1]) / eps_init[1]
                results["rftrack"].append(d)
                print(f"  {'rftrack':14s}: epsny_growth={d['epsny_growth']*100:7.3g}% "
                      f"alive={d['alive']:4d}/{N_P} wall={d['wall']:.1f}s")
            except Exception as e:
                print(f"  {'rftrack':14s}: FAILED {type(e).__name__}: {e}")
                results["rftrack"].append(dict(Q_nC=Q_nC, error=str(e)))

    # Write per-engine CSVs and a unified comparison CSV
    for engine, rows in results.items():
        if not rows:
            continue
        path = OUT_DIR / f"sweep_{engine}.csv"
        keys = sorted({k for r in rows for k in r.keys()})
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"  -> {path}")

    return results


def plot_compare(results: dict[str, list[dict]]) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    color_map = {
        "dafmm": ("k", "o", "DA-FMM (COSY)"),
        "xsuite-frozen": ("C0", "s", "xsuite frozen-Gaussian"),
        "xsuite-pic3d": ("C2", "^", "xsuite 3D PIC"),
        "rftrack": ("C3", "D", "RF-Track PIC"),
    }

    for engine, rows in results.items():
        if engine not in color_map:
            continue
        col, mk, lab = color_map[engine]
        valid = [r for r in rows if "error" not in r]
        if not valid:
            continue
        Q = np.array([r["Q_nC"] for r in valid])
        gx = np.array([r["epsnx_growth"] for r in valid]) * 100
        gy = np.array([r["epsny_growth"] for r in valid]) * 100
        axes[0].plot(Q, gx, marker=mk, color=col, label=lab, lw=1.4)
        axes[1].plot(Q, gy, marker=mk, color=col, label=lab, lw=1.4)

    for ax, title in zip(axes, [r"$\Delta\varepsilon_{n,x}$", r"$\Delta\varepsilon_{n,y}$"]):
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=0.1)
        ax.set_xlabel(r"$Q_\mathrm{bunch}$  [nC]")
        ax.set_ylabel(f"{title} growth at exit  [%]")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=9)
        ax.set_title(title)

    fig.suptitle(
        "L3.6 four-engine SC comparison on the spch_demo lattice "
        "(45 MeV, 2 m FODO, N_p=6000)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = OUT_DIR / "compare_sc_engines.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=140)
    plt.close(fig)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--engines",
        nargs="+",
        default=["dafmm", "xsuite-frozen", "xsuite-pic3d"],
        help=(
            "Which engines to run/plot. RF-Track Volume+SC pattern dumped "
            "core under the lat-in-Volume add() approach used here; "
            "deferred to L3.6-Phase-2."
        ),
    )
    args = ap.parse_args()

    results = sweep(args.engines)
    out = plot_compare(results)
    print(f"\nFigure: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
