#!/usr/bin/env python3
"""11-stage NM objective-design ablation for IPAC paper benchmarking.

Runs the seminar Twiss-matching with three objective configurations:
  A: verbatim original (Niels's branch_niels/UHM_beamline_opt.py)
  B: A + per-measure-type rescaling so squared residuals are
     dimensionally comparable across alpha / beta / envelope / dispersion
  C: B + typo fix (Stage 1 x.beta weight 0.0 -> 0.5)
       + envelope=0.0 -> 1.5 mm (Stage 7 finite physical target)

Outputs JSON with config, seed, currents, final Twiss, per-stage MSE
traces, total iterations, runtime.

Usage:
    python ablation_run.py --config A --seed 42 --out results/ablation/A_seed42.json

Author: Eremey Valetov
"""

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")

script_dir = Path(__file__).resolve().parent

if (script_dir / "ebeam.py").exists():
    backend_dir = script_dir
    EXCEL_PATH = script_dir / "Beamline_elements.xlsx"
else:
    backend_dir = script_dir.parent
    EXCEL_PATH = backend_dir.parent / "beam_excel" / "Beamline_elements.xlsx"

sys.path.insert(0, str(backend_dir))

from ebeam import beam
from beamline import lattice
from excelElements import ExcelElements
from beamOptimizer import beamOptimizer

# Beam parameters (match seminar reproducer)
ENERGY = 40
EPSILON_N = 8
X_STD = 0.8
Y_STD = 0.8
FREQ = 2856e6
BUNCH_SPREAD = 2
ENERGY_STD_PCT = 0.5
H_CHIRP = 5e9
NB_PARTICLES = 5000

K_UND = 1.2
LAMBDA_U = 2.3e-2
QUAD_INDICES = [
    1, 3, 10, 16, 18, 20, 27, 33, 35, 37, 39, 41, 43,
    50, 56, 58, 61, 63, 70, 76, 78, 80, 87, 93, 95, 97,
]

STAGE_LABELS = [
    "Doublet", "Chicane 1 corr.", "Triplet A",
    "Chicane 2 corr.", "Symm. triplet", "Chicane 3 corr.",
    "Doublet B", "Triplet B", "Chicane 4 corr.",
    "Triplet C", "Undulator match",
]

# Per-measure reference scales for Config B/C rescaling. Squared residuals
# are normalised by ref**2. Mild design: only dispersion is rescaled, with
# a 0.5 m reference (weight x 4) to give it a modest sensitivity boost
# without overwhelming the Twiss-matching objectives. Other measures keep
# the natural unit scale.
MEASURE_REF = {
    "alpha": 1.0,        # dimensionless, no rescale
    "beta": 1.0,         # m, no rescale (typical beta values O(1))
    "envelope": 1.0,     # mm, no rescale (typical 0-3 mm)
    "dispersion": 0.5,   # m, weight x 4 -- mild sensitivity boost
}


def compute_params(beta_xm=1.4, alpha_xm=0.47):
    """Stage 11 horizontal Twiss targets are CLI-overridable.

    Defaults (1.4 m, 0.47) come from Weinberg, Fisher & Li
    (arXiv:2510.14061v1) Table I -- radiation-mode matching.
    Bidault, Weinberg, Purwar & Li (MOP6318, IPAC 2026, Monday poster) Table 1
    uses (1.267 m, 0.560), derived from the Rayleigh-length
    formula beta_x,o = Z_R * eps_r / eps_x.
    Vertical-plane target beta_y = gamma*lambda_u/(2*pi*K) is
    fixed by FEL natural focusing.
    """
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    norm = relat.gamma * relat.beta
    epsilon = EPSILON_N / norm
    beta_ym = relat.gamma * LAMBDA_U / (2 * np.pi * K_UND)
    return {
        "beta_xm": beta_xm, "alpha_xm": alpha_xm,
        "beta_ym": beta_ym, "alpha_ym": 0.0,
        "epsilon": epsilon, "gamma": relat.gamma,
        "beta_rel": relat.beta, "norm": norm,
    }


def generate_beam(params, seed, method="random", n=NB_PARTICLES):
    np.random.seed(seed)
    eps = params["epsilon"]
    ebeam_gen = beam()
    dist = ebeam_gen.gen_6d_gaussian(
        0,
        [X_STD, eps / X_STD, Y_STD, eps / Y_STD,
         BUNCH_SPREAD * 1e-9 * FREQ, ENERGY_STD_PCT * 10],
        n,
        method=method,
    )
    dist[:, 5] += H_CHIRP * dist[:, 4] / FREQ
    return dist


def build_line():
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    excel = ExcelElements(str(EXCEL_PATH))
    line = relat.changeBeamType(
        "electron", ENERGY, excel.create_beamline()
    )[:118]
    return line


def make_obj(plane, measure_type, goal, weight, config):
    """Build an objective dict, applying config-specific rescaling.

    Configs B and C divide weight by MEASURE_REF[measure_type]**2 so the
    squared residual term ((stat - goal)/ref)**2 contributes comparable
    magnitude across measure types.
    """
    if config == "A":
        w = weight
    else:
        ref = MEASURE_REF.get(measure_type, 1.0)
        w = weight / (ref ** 2)
    return {"measure": [plane, measure_type], "goal": goal, "weight": w}


def build_stages(params, config):
    """Construct (stages_pre_mirror, stages_post_mirror) for the given config."""
    axm, aym = params["alpha_xm"], params["alpha_ym"]
    bxm, bym = params["beta_xm"], params["beta_ym"]

    # Config-C-specific overrides
    s1_xbeta_weight = 0.5 if config == "C" else 0.0
    s7_env_goal = 1.5 if config == "C" else 0.0

    o = lambda *a: make_obj(*a, config=config)

    stages_pre = [
        # Stage 1: Doublet (segs 8, 9)
        ({1: ["I", "current", lambda n: n], 3: ["I2", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1},
          "I2": {"bounds": (0, 10), "start": 1}},
         {8: [o("x", "alpha", 0, 1),
              o("x", "beta", 0.1, s1_xbeta_weight)],
          9: [o("y", "alpha", 0, 1),
              o("y", "beta", 0.1, 0.5)]}),
        # Stage 2: Chicane 1 dispersion (seg 15)
        ({10: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {15: [o("x", "dispersion", 0, 1)]}),
        # Stage 3: Triplet A (segs 25, 26)
        ({16: ["I", "current", lambda n: n], 18: ["I2", "current", lambda n: n],
          20: ["I3", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2},
          "I2": {"bounds": (0, 10), "start": 5},
          "I3": {"bounds": (0, 10), "start": 3}},
         {25: [o("x", "alpha", 0, 1),
               o("x", "beta", 0.1, 0.5)],
          26: [o("y", "alpha", 0, 1),
               o("y", "beta", 0.1, 0.5)]}),
        # Stage 4: Chicane 2 dispersion (seg 32)
        ({27: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {32: [o("x", "dispersion", 0, 1)]}),
        # Stage 5: Symm. triplet (seg 37)
        ({37: ["I", "current", lambda n: n], 35: ["I2", "current", lambda n: n],
          33: ["I3", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2},
          "I2": {"bounds": (0, 10), "start": 2},
          "I3": {"bounds": (0, 10), "start": 2}},
         {37: [o("x", "alpha", 0, 1),
               o("y", "alpha", 0, 1),
               o("x", "envelope", 2.0, 1),
               o("y", "envelope", 2.0, 1)]}),
    ]

    stages_post = [
        # Stage 6: Chicane 3 dispersion (seg 55)
        ({50: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {55: [o("x", "dispersion", 0, 1)]}),
        # Stage 7: Doublet B (seg 59) -- envelope target placeholder in A/B
        ({56: ["I", "current", lambda n: n], 58: ["I2", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2},
          "I2": {"bounds": (0, 10), "start": 2}},
         {59: [o("x", "envelope", s7_env_goal, 1),
               o("y", "envelope", s7_env_goal, 1)]}),
        # Stage 8: Triplet B (segs 68, 69)
        ({61: ["I", "current", lambda n: n], 63: ["I2", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2},
          "I2": {"bounds": (0, 10), "start": 2}},
         {68: [o("x", "alpha", 0, 1),
               o("x", "beta", 0.1, 0.5)],
          69: [o("y", "alpha", 0, 1),
               o("y", "beta", 0.1, 0.5)]}),
        # Stage 9: Chicane 4 dispersion (seg 75)
        ({70: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {75: [o("x", "dispersion", 0, 1)]}),
        # Stage 10: Triplet C (segs 85, 86)
        ({76: ["I", "current", lambda n: n], 78: ["I2", "current", lambda n: n],
          80: ["I3", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2},
          "I2": {"bounds": (0, 10), "start": 2},
          "I3": {"bounds": (0, 10), "start": 2}},
         {85: [o("x", "alpha", 0, 1),
               o("x", "beta", 0.1, 0.5)],
          86: [o("y", "alpha", 0, 1),
               o("y", "beta", 0.1, 0.5)]}),
        # Stage 11: Undulator match (segs 92, 117)
        ({87: ["Ic", "current", lambda n: n], 93: ["I", "current", lambda n: n],
          95: ["I2", "current", lambda n: n], 97: ["I3", "current", lambda n: n]},
         {"Ic": {"bounds": (0, 10), "start": 4},
          "I": {"bounds": (0, 10), "start": 2},
          "I2": {"bounds": (0, 10), "start": 2},
          "I3": {"bounds": (0, 10), "start": 2}},
         {92: [o("x", "dispersion", 0, 0.5)],
          117: [o("x", "alpha", axm, 1),
                o("y", "alpha", aym, 1),
                o("x", "beta", bxm, 1),
                o("y", "beta", bym, 1)]}),
    ]

    return stages_pre, stages_post


def run_ablation(line, beam_dist, params, config):
    opti = beamOptimizer(line, beam_dist)
    stages_pre, stages_post = build_stages(params, config)

    stage_traces = []

    def record_stage(stage_idx):
        stage_traces.append({
            "stage": stage_idx,
            "label": STAGE_LABELS[stage_idx - 1] if stage_idx <= len(STAGE_LABELS) else "?",
            "iter_count": len(opti.plotMSE),
            "mse_trace": [float(m) for m in opti.plotMSE],
            "final_mse": float(opti.plotMSE[-1]) if opti.plotMSE else None,
        })

    for i, (seg, sp, obj) in enumerate(stages_pre, 1):
        opti.calc("Nelder-Mead", seg, sp, obj)
        record_stage(i)

    line[43].current = line[33].current
    line[41].current = line[35].current
    line[39].current = line[37].current

    for i, (seg, sp, obj) in enumerate(stages_post, len(stages_pre) + 1):
        opti.calc("Nelder-Mead", seg, sp, obj)
        record_stage(i)

    currents = {idx: float(line[idx].current) for idx in QUAD_INDICES}
    return currents, stage_traces


def compute_final_twiss(line, beam_dist):
    ebeam_calc = beam()
    particles = beam_dist.copy()
    for elem in line:
        particles = np.array(elem.useMatrice(particles))
    _, _, tw = ebeam_calc.cal_twiss(particles, ddof=1)
    return {
        "beta_x": float(tw.loc["x", r"$\beta$ (m)"]),
        "beta_y": float(tw.loc["y", r"$\beta$ (m)"]),
        "alpha_x": float(tw.loc["x", r"$\alpha$"]),
        "alpha_y": float(tw.loc["y", r"$\alpha$"]),
    }


def compute_undulator_rms(final_twiss, params):
    targets = {
        "beta_x": params["beta_xm"], "alpha_x": params["alpha_xm"],
        "beta_y": params["beta_ym"], "alpha_y": params["alpha_ym"],
    }
    sse = sum((final_twiss[k] - targets[k]) ** 2 for k in targets)
    return math.sqrt(sse / 4)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", choices=["A", "B", "C"], required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--beam-method", choices=["random", "sobol"], default="random")
    p.add_argument("--beta-xm", type=float, default=1.4,
                   help="Stage 11 target beta_x in m (default 1.4 from "
                        "arXiv:2510.14061v1; use 1.267 for MOP6318)")
    p.add_argument("--alpha-xm", type=float, default=0.47,
                   help="Stage 11 target alpha_x (default 0.47 from "
                        "arXiv:2510.14061v1; use 0.560 for MOP6318)")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    print(f"=== ablation_run config={args.config} seed={args.seed} ===")
    t0 = time.perf_counter()

    params = compute_params(beta_xm=args.beta_xm, alpha_xm=args.alpha_xm)
    line = build_line()
    beam_dist = generate_beam(params, args.seed, method=args.beam_method)
    currents, stage_traces = run_ablation(line, beam_dist, params, args.config)
    final_twiss = compute_final_twiss(line, beam_dist)
    undulator_rms = compute_undulator_rms(final_twiss, params)
    runtime = time.perf_counter() - t0
    total_iters = sum(s["iter_count"] for s in stage_traces)

    output = {
        "config": args.config,
        "seed": args.seed,
        "beam_method": args.beam_method,
        "beta_xm": args.beta_xm,
        "alpha_xm": args.alpha_xm,
        "nb_particles": NB_PARTICLES,
        "currents": currents,
        "final_twiss": final_twiss,
        "targets": {
            "beta_x": params["beta_xm"], "alpha_x": params["alpha_xm"],
            "beta_y": float(params["beta_ym"]), "alpha_y": params["alpha_ym"],
        },
        "undulator_rms": undulator_rms,
        "total_iters": total_iters,
        "runtime_sec": runtime,
        "stage_traces": stage_traces,
    }

    with open(args.out, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\n  Final Twiss: bx={final_twiss['beta_x']:.4f}, "
          f"by={final_twiss['beta_y']:.4f}, "
          f"ax={final_twiss['alpha_x']:.4f}, "
          f"ay={final_twiss['alpha_y']:.4f}")
    print(f"  Undulator RMS: {undulator_rms:.3e}")
    print(f"  Total iterations: {total_iters}")
    print(f"  Runtime: {runtime:.1f} s")
    print(f"  Wrote: {args.out}")


if __name__ == "__main__":
    main()
