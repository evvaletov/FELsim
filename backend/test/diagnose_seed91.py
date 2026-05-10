#!/usr/bin/env python3
"""Diagnose why a particular seed (default 91) fails the 11-stage match.

Two diagnostics:
  1. Input-beam moments comparison: are the failing seed's input
     statistics anomalous (skew, kurtosis, max excursion) vs working
     seeds? If yes, the failure is beam-driven.
  2. Per-stage MSE comparison: load the ablation result JSONs and
     print final-MSE-per-stage for a working seed vs the failing
     seed, per config. Identifies which stage diverges.

Usage:
    python diagnose_seed91.py [--seeds ...] [--target 91]
        [--results-dir DIR] [--compare-seed 42]

Author: Eremey Valetov
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import skew, kurtosis

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent))

from ebeam import beam
from beamline import lattice

ENERGY = 40
EPSILON_N = 8
X_STD = 0.8
Y_STD = 0.8
FREQ = 2856e6
BUNCH_SPREAD = 2
ENERGY_STD_PCT = 0.5
H_CHIRP = 5e9
NB_PARTICLES = 5000

COORD_LABELS = ["x", "x'", "y", "y'", "t", "dK/K"]


def compute_params():
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    norm = relat.gamma * relat.beta
    return EPSILON_N / norm


def generate_beam(seed, eps):
    np.random.seed(seed)
    ebeam_gen = beam()
    dist = ebeam_gen.gen_6d_gaussian(
        0,
        [X_STD, eps / X_STD, Y_STD, eps / Y_STD,
         BUNCH_SPREAD * 1e-9 * FREQ, ENERGY_STD_PCT * 10],
        NB_PARTICLES,
    )
    dist[:, 5] += H_CHIRP * dist[:, 4] / FREQ
    return dist


def beam_stats(dist):
    """Return dict of per-coordinate moments + max excursions."""
    stats = {}
    for i, lbl in enumerate(COORD_LABELS):
        col = dist[:, i]
        stats[lbl] = {
            "mean": float(np.mean(col)),
            "std": float(np.std(col, ddof=1)),
            "skew": float(skew(col)),
            "kurt": float(kurtosis(col)),
            "max_abs": float(np.max(np.abs(col))),
            "max_abs_sigma": float(np.max(np.abs(col)) / np.std(col, ddof=1)),
        }
    return stats


def per_stage_mse_comparison(results_dir: Path, target_seed: int,
                              compare_seed: int):
    """Load ablation JSONs and print final MSE per stage for the
    comparison seed (working) vs target seed (failing), per config."""
    print(f"\n=== Per-stage final MSE: seed {compare_seed} vs {target_seed} ===\n")
    runs = {}
    for cfg in "ABC":
        for seed in (compare_seed, target_seed):
            fp = results_dir / f"{cfg}_seed{seed}.json"
            if not fp.exists():
                print(f"  missing: {fp}")
                return
            with open(fp) as f:
                runs[(cfg, seed)] = json.load(f)

    header = (f"{'Stage':>5} {'label':<18}"
              f" {'A_s'+str(compare_seed):>10} {'A_s'+str(target_seed):>10}"
              f" {'B_s'+str(compare_seed):>10} {'B_s'+str(target_seed):>10}"
              f" {'C_s'+str(compare_seed):>10} {'C_s'+str(target_seed):>10}")
    print(header)
    print("-" * len(header))
    for i in range(11):
        label = runs[("A", compare_seed)]["stage_traces"][i]["label"]
        cells = [runs[(cfg, s)]["stage_traces"][i]["final_mse"]
                 for cfg in "ABC" for s in (compare_seed, target_seed)]
        print(f"{i+1:>5} {label:<18}", " ".join(f"{c:>10.2e}" for c in cells))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seeds", type=int, nargs="+",
                   default=[1, 17, 23, 42, 91, 137])
    p.add_argument("--target", type=int, default=91)
    p.add_argument("--compare-seed", type=int, default=42,
                   help="working seed to compare against (for per-stage MSE)")
    p.add_argument("--results-dir", type=Path,
                   default=Path(__file__).parent / "results" / "ablation")
    args = p.parse_args()

    eps = compute_params()
    print(f"Beam stats for {NB_PARTICLES} particles, "
          f"epsilon = {eps:.6f}\n")

    all_stats = {}
    for seed in args.seeds:
        dist = generate_beam(seed, eps)
        all_stats[seed] = beam_stats(dist)

    print(f"{'Seed':>5}", end="")
    for lbl in COORD_LABELS:
        print(f"  {lbl+' max/sig':>11}", end="")
    print()
    for seed in args.seeds:
        marker = " <-- TARGET" if seed == args.target else ""
        print(f"{seed:>5}", end="")
        for lbl in COORD_LABELS:
            print(f"  {all_stats[seed][lbl]['max_abs_sigma']:>11.3f}", end="")
        print(marker)

    print()
    print("Per-coordinate std, mean, skew, kurt (each row a seed):")
    print()
    for lbl in COORD_LABELS:
        print(f"--- {lbl} ---")
        print(f"{'Seed':>5} {'mean':>11} {'std':>11} "
              f"{'skew':>9} {'kurt':>9} {'max_abs':>11} {'max/sig':>9}")
        for seed in args.seeds:
            s = all_stats[seed][lbl]
            marker = "<<" if seed == args.target else "  "
            print(f"{seed:>5} {s['mean']:>11.4e} {s['std']:>11.4e} "
                  f"{s['skew']:>9.4f} {s['kurt']:>9.4f} "
                  f"{s['max_abs']:>11.4e} {s['max_abs_sigma']:>9.3f} {marker}")
        print()

    if args.results_dir.exists():
        per_stage_mse_comparison(args.results_dir, args.target,
                                  args.compare_seed)


if __name__ == "__main__":
    main()
