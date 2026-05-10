#!/usr/bin/env python3
"""Analyze 11-stage NM ablation results: load JSONs, generate plots + report.

Inputs: results/ablation/{A,B,C}_seed{42,17,23,91,137}.json (15 files).
Outputs: results/ablation/analysis/
  - convergence_per_stage.{pdf,png}: 11 panels, MSE-vs-iteration overlay
  - final_rms_boxplot.{pdf,png}: distribution of final RMS per config
  - iters_boxplot.{pdf,png}: total iteration count per config
  - summary.md: markdown table + interpretation

Usage:
    python ablation_analyze.py [--results-dir DIR]

Author: Eremey Valetov
"""

import argparse
import json
import math
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CONFIGS = ["A", "B", "C"]
CONFIG_LABELS = {
    "A": "A: verbatim original",
    "B": "B: + mild rescaling",
    "C": "C: + typo + envelope fix",
}
CONFIG_COLORS = {"A": "#0077BB", "B": "#EE7733", "C": "#009988"}
SEEDS = [1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 42, 43, 47, 53, 91, 137]
STAGE_LABELS = [
    "Doublet", "Chicane 1", "Triplet A",
    "Chicane 2", "Symm. triplet", "Chicane 3",
    "Doublet B", "Triplet B", "Chicane 4",
    "Triplet C", "Undulator match",
]


def load_runs(results_dir: Path):
    """Load all random-seed runs + (optionally) a Sobol run per config."""
    runs = {}
    sobol = {}
    for cfg in CONFIGS:
        runs[cfg] = []
        for seed in SEEDS:
            fp = results_dir / f"{cfg}_seed{seed}.json"
            if not fp.exists():
                continue
            with open(fp) as f:
                runs[cfg].append(json.load(f))
        # Sobol variant (optional)
        sobol_fp = results_dir / f"{cfg}_sobol.json"
        if sobol_fp.exists():
            with open(sobol_fp) as f:
                sobol[cfg] = json.load(f)
    return runs, sobol


def plot_convergence(runs, sobol, out_path: Path):
    fig, axes = plt.subplots(3, 4, figsize=(16, 9), sharey=True)
    axes = axes.flatten()

    for stage_idx in range(11):
        ax = axes[stage_idx]
        for cfg in CONFIGS:
            for run in runs[cfg]:
                if stage_idx >= len(run["stage_traces"]):
                    continue
                trace = run["stage_traces"][stage_idx]
                rms = [math.sqrt(max(m, 0)) for m in trace["mse_trace"]]
                ax.semilogy(rms, color=CONFIG_COLORS[cfg], alpha=0.3, lw=0.6)
            # Sobol overlay (dashed, full opacity)
            if cfg in sobol and stage_idx < len(sobol[cfg]["stage_traces"]):
                trace = sobol[cfg]["stage_traces"][stage_idx]
                rms = [math.sqrt(max(m, 0)) for m in trace["mse_trace"]]
                ax.semilogy(rms, color=CONFIG_COLORS[cfg], alpha=1.0, lw=1.5,
                            linestyle="--")
        ax.set_title(f"Stage {stage_idx+1}: {STAGE_LABELS[stage_idx]}", fontsize=9)
        ax.tick_params(labelsize=7)
        if stage_idx % 4 == 0:
            ax.set_ylabel(r"$\sqrt{\mathrm{MSE}}$", fontsize=9)
        if stage_idx >= 7:
            ax.set_xlabel("iteration", fontsize=9)

    axes[11].axis("off")
    handles = [plt.Line2D([0], [0], color=CONFIG_COLORS[c], lw=2,
                          label=CONFIG_LABELS[c]) for c in CONFIGS]
    handles.append(plt.Line2D([0], [0], color="black", lw=1.5,
                              linestyle="--", label="Sobol (deterministic)"))
    axes[11].legend(handles=handles, loc="center", fontsize=10,
                    title=f"{len(SEEDS)} seeds + 1 Sobol per config")

    fig.suptitle("Per-stage NM convergence — A vs B vs C, multi-seed overlay",
                 fontsize=13, y=0.995)
    plt.tight_layout()
    for fmt in ("pdf", "png"):
        fig.savefig(out_path.with_suffix(f".{fmt}"))
    plt.close(fig)


def plot_final_rms_boxplot(runs, sobol, out_path: Path):
    fig, ax = plt.subplots(figsize=(9, 5))
    data = [[r["undulator_rms"] for r in runs[cfg]] for cfg in CONFIGS]
    positions = range(1, len(CONFIGS) + 1)
    bp = ax.boxplot(data, positions=positions, widths=0.5, patch_artist=True,
                    showmeans=True, meanline=True)
    for patch, cfg in zip(bp["boxes"], CONFIGS):
        patch.set_facecolor(CONFIG_COLORS[cfg])
        patch.set_alpha(0.5)
    ax.set_xticks(positions)
    ax.set_xticklabels([CONFIG_LABELS[c] for c in CONFIGS], fontsize=10)
    ax.set_ylabel("Undulator RMS at convergence", fontsize=11)
    ax.set_yscale("log")
    ax.grid(True, axis="y", alpha=0.3, which="both")
    ax.set_title(f"Final undulator RMS — {len(SEEDS)} seeds per config",
                 fontsize=12)
    for i, cfg in enumerate(CONFIGS, 1):
        for r in runs[cfg]:
            ax.scatter(i, r["undulator_rms"], color="black",
                       alpha=0.5, s=15, zorder=3)
        if cfg in sobol:
            ax.scatter(i, sobol[cfg]["undulator_rms"], color="red",
                       marker="*", s=150, zorder=4,
                       label="Sobol" if i == 1 else None)
    if any(c in sobol for c in CONFIGS):
        ax.legend(loc="upper left", fontsize=9)
    plt.tight_layout()
    for fmt in ("pdf", "png"):
        fig.savefig(out_path.with_suffix(f".{fmt}"))
    plt.close(fig)


def plot_iters_boxplot(runs, out_path: Path):
    fig, ax = plt.subplots(figsize=(9, 5))
    data = [[r["total_iters"] for r in runs[cfg]] for cfg in CONFIGS]
    positions = range(1, len(CONFIGS) + 1)
    bp = ax.boxplot(data, positions=positions, widths=0.5, patch_artist=True,
                    showmeans=True, meanline=True)
    for patch, cfg in zip(bp["boxes"], CONFIGS):
        patch.set_facecolor(CONFIG_COLORS[cfg])
        patch.set_alpha(0.5)
    ax.set_xticks(positions)
    ax.set_xticklabels([CONFIG_LABELS[c] for c in CONFIGS], fontsize=10)
    ax.set_ylabel("Total NM iterations (sum across 11 stages)", fontsize=11)
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_title(f"Total iterations — {len(SEEDS)} seeds per config",
                 fontsize=12)
    for i, cfg in enumerate(CONFIGS, 1):
        for r in runs[cfg]:
            ax.scatter(i, r["total_iters"], color="black",
                       alpha=0.5, s=15, zorder=3)
    plt.tight_layout()
    for fmt in ("pdf", "png"):
        fig.savefig(out_path.with_suffix(f".{fmt}"))
    plt.close(fig)


def write_summary(runs, sobol, out_path: Path):
    lines = []
    n_seeds_per_cfg = len(runs[CONFIGS[0]]) if runs[CONFIGS[0]] else 0
    lines.append("# 11-stage NM objective-design ablation\n")
    lines.append(
        f"Three objective configurations were run on the UH FEL transport "
        f"line with {n_seeds_per_cfg} random-seed beams each, plus a "
        f"single deterministic Sobol beam per config. NM is sequential "
        f"per-stage; starting points are reset per stage to the "
        f"configured defaults.\n"
    )
    lines.append("## Configurations\n")
    lines.append("- **A**: Verbatim original (`branch_niels/UHM_beamline_opt.py`).\n")
    lines.append(
        "- **B**: A + mild per-measure-type rescaling. Each squared "
        "residual is divided by `MEASURE_REF**2`; only dispersion is "
        "actually rescaled (`ref = 0.5 m`, weight x4). Other measures "
        "have `ref = 1` so their numerical contribution is unchanged.\n"
    )
    lines.append(
        "- **C**: B + Stage 1 `x.beta` weight `0.0 -> 0.5` (typo fix) + "
        "Stage 7 `envelope` goal `0.0 mm -> 1.5 mm` (finite physical "
        "target).\n"
    )

    n_seeds = len(runs[CONFIGS[0]]) if runs[CONFIGS[0]] else 0
    lines.append(f"\n## Final undulator RMS ({n_seeds} seeds per config)\n")
    lines.append("| Config | Median | Min | Max | Mean | Std | Fail rate (RMS>1e-2) | Sobol |\n")
    lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |\n")
    for cfg in CONFIGS:
        rmss = [r["undulator_rms"] for r in runs[cfg]]
        if not rmss:
            lines.append(f"| {cfg} | (no runs) |\n")
            continue
        n_fail = sum(1 for r in rmss if r > 1e-2)
        sobol_rms = sobol.get(cfg, {}).get("undulator_rms")
        sobol_str = f"{sobol_rms:.3e}" if sobol_rms is not None else "(none)"
        lines.append(
            f"| {cfg} | {statistics.median(rmss):.3e} | "
            f"{min(rmss):.3e} | {max(rmss):.3e} | "
            f"{statistics.mean(rmss):.3e} | "
            f"{statistics.stdev(rmss) if len(rmss) > 1 else 0:.3e} | "
            f"{n_fail}/{len(rmss)} ({100*n_fail/len(rmss):.0f}%) | "
            f"{sobol_str} |\n"
        )

    lines.append(f"\n## Total NM iterations ({n_seeds_per_cfg} seeds per config)\n")
    lines.append("| Config | Median | Min | Max | Mean |\n")
    lines.append("| --- | --- | --- | --- | --- |\n")
    for cfg in CONFIGS:
        iters = [r["total_iters"] for r in runs[cfg]]
        if not iters:
            lines.append(f"| {cfg} | (no runs) |\n")
            continue
        lines.append(
            f"| {cfg} | {statistics.median(iters):.0f} | "
            f"{min(iters)} | {max(iters)} | "
            f"{statistics.mean(iters):.0f} |\n"
        )

    lines.append("\n## Per-seed detail\n")
    lines.append("| Config | Seed | RMS | Total iters | Runtime (s) |\n")
    lines.append("| --- | --- | --- | --- | --- |\n")
    for cfg in CONFIGS:
        for r in runs[cfg]:
            lines.append(
                f"| {cfg} | {r['seed']} | {r['undulator_rms']:.3e} | "
                f"{r['total_iters']} | {r['runtime_sec']:.1f} |\n"
            )

    # Preserve any manual interpretation already in summary.md.
    # Marker delimits auto-generated tables from the manual narrative.
    MARKER = "<!-- INTERPRETATION (preserve when regenerating) -->"
    lines.append(f"\n{MARKER}\n\n")

    existing = ""
    if out_path.exists():
        prev = out_path.read_text()
        if MARKER in prev:
            existing = prev.split(MARKER, 1)[1].lstrip("\n")

    if not existing:
        existing = (
            "## Interpretation\n\n"
            "*To be written by hand after reviewing the numbers.*\n"
        )

    out_path.write_text("".join(lines) + existing)


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--results-dir",
        type=Path,
        default=Path(__file__).parent / "results" / "ablation",
    )
    args = p.parse_args()

    out_dir = args.results_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    runs, sobol = load_runs(args.results_dir)
    n_total = sum(len(v) for v in runs.values())
    print(f"Loaded {n_total} random-seed runs + {len(sobol)} Sobol runs "
          f"across {len(runs)} configs.")

    plot_convergence(runs, sobol, out_dir / "convergence_per_stage")
    plot_final_rms_boxplot(runs, sobol, out_dir / "final_rms_boxplot")
    plot_iters_boxplot(runs, out_dir / "iters_boxplot")
    write_summary(runs, sobol, out_dir / "summary.md")

    print(f"Wrote analysis to {out_dir}")


if __name__ == "__main__":
    main()
