#!/usr/bin/env python3
"""FF1: fringe-field (FF) mode comparison for the UH MkV FEL transport line.

Refreshes the R2 cross-code fringe study with the two updates Niels asked for.

Two valid comparisons are produced:

  Table A (computed here): COSY INFINITY at fringe orders FR0/FR1/FR2/FR3, each
  warm-started to the undulator match. All four reach the target beta_y with a
  DIFFERENT current set, quantified by the RMS current shift vs FR0. This is the
  clean intra-code fringe-order result.

  Cross-code status (cited): a cross-CODE fringe comparison cannot be done by
  applying one code's currents to another -- each code needs its own currents to
  match (R2: different currents -> same Twiss). The established results are
  summarised, with the one new enabler (xsuite dipole bend now builds after the
  xsuiteAdapter h->angle fix) and the remaining step (re-optimise RF-Track and
  xsuite to the target to quantify the now-closed beta_y deficit).

Author: Eremey Valetov
"""

import sys
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

ENERGY = 40.0
K = 1.2
LAMBDA_U = 2.3e-2
RESULTS = Path(__file__).resolve().parent / "results"
OUT = RESULTS / "niels_june"
OUT.mkdir(parents=True, exist_ok=True)

MODES = [
    ("cosy_s1_fr0.json", "FR0 (hard-edge)"),
    ("cosy_s1_fr1_warm.json", "FR1 (Enge 1st, warm)"),
    ("cosy_s1_fr2_warm.json", "FR2 (Enge 2nd, warm)"),
    ("cosy_s1_fr3_warm.json", "FR3 (RK fringe, warm)"),
]


def beta_ym():
    from beamline import lattice
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    return relat.gamma * LAMBDA_U / (2 * np.pi * K)


def load_mode(json_path):
    d = json.load(open(json_path))
    tw = d['twiss_undulator']
    currents = {int(k): float(v) for k, v in d['currents'].items()}
    return {'beta_x': tw['beta_x'], 'beta_y': tw['beta_y'],
            'mse': d.get('mse'), 'fr': d.get('fringe_field_order'),
            'currents': currents}


def rms_current_shift(ref, other):
    keys = sorted(set(ref) & set(other))
    if not keys:
        return float('nan')
    d = np.array([abs(other[k]) - abs(ref[k]) for k in keys])
    return float(np.sqrt(np.mean(d ** 2)))


def main():
    bym = beta_ym()
    print(f"Undulator target beta_y = {bym:.4f} m (gamma*lambda_u/2piK)")

    modes = []
    for jn, label in MODES:
        p = RESULTS / jn
        if p.exists():
            modes.append((label, load_mode(p)))
        else:
            print(f"  missing {jn}")
    ref_currents = modes[0][1]['currents']  # FR0

    md = []
    md.append("# FF-mode comparison - UH MkV FEL transport line\n")
    md.append(f"Author: Eremey Valetov. FF = fringe field. Undulator target "
              f"beta_y = {bym:.3f} m, beta_x = 1.40 m.\n")

    md.append("## Table A - COSY fringe-order convergence (computed 2026-06-19)\n")
    md.append("Each fringe order is warm-started to the undulator match. All reach "
              "the target with a different current set; the last column is the RMS "
              "quad-current shift relative to FR0.\n")
    md.append("| COSY mode | beta_x (m) | beta_y (m) | match MSE | RMS |dI| vs FR0 (A) |")
    md.append("|---|---:|---:|---:|---:|")
    for label, m in modes:
        shift = rms_current_shift(ref_currents, m['currents'])
        md.append(f"| {label} | {m['beta_x']:.3f} | {m['beta_y']:.3f} | "
                  f"{m['mse']:.3g} | {shift:.3f} |")
    md.append("")

    md.append("## Cross-code fringe status\n")
    md.append(
        "A cross-CODE fringe comparison must re-optimise each code to the target "
        "(applying one code's currents to another does not match: R2 established "
        "different currents -> same Twiss). Established results:\n")
    md.append(
        "- FELsim (triangle-rule DPW edge) and COSY (FR0/FR1) agree on matched "
        "Twiss (R2: MSE 1.2e-6 / 2.3e-7 / 3.7e-9 at eps_n=8).")
    md.append(
        "- RF-Track: the R2-era vertical beta_y deficit (0.055 vs 0.242 m) came "
        "from a missing triangle-rule edge phi. The triangle-phi correction was "
        "since added (rftrackAdapter._annotate_dipole_edges); RF-Track now "
        "re-optimises to the target (RFT-opt MSE 7.0e-3, eps_n=8), i.e. the "
        "deficit is closed.")
    md.append(
        "- xsuite: no dipole edge/fringe model (drift edges) -> the no-fringe "
        "baseline. Its bend body now builds after the xsuiteAdapter h->angle fix "
        "(2026-06-19); quantifying its matched beta_y needs a re-optimisation run "
        "(remaining quick task).")
    md.append("")
    md.append("REMAINING (post-Monday): re-optimise RF-Track and xsuite to the "
              "undulator target and tabulate achieved beta_y alongside COSY "
              "FR0-FR3 and FELsim for a single 5-code FF figure.\n")

    md_text = "\n".join(md)
    (OUT / "FF1_fringe_mode_comparison.md").write_text(md_text)
    print("\n" + md_text)

    # Plot: beta_y per COSY fringe mode vs target (all converge).
    labels = [m[0].split(" ")[0] for m in modes]
    byvals = [m[1]['beta_y'] for m in modes]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(labels, byvals, color='#1f77b4')
    ax.axhline(bym, ls='--', color='k', lw=1, label=f'target {bym:.3f} m')
    ax.set_ylabel(r'$\beta_y$ at undulator (m)')
    ax.set_ylim(0, max(byvals + [bym]) * 1.3)
    ax.set_title('COSY fringe-order convergence to the undulator match')
    ax.legend()
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f"FF1_cosy_fringe_convergence.{ext}", dpi=150)
    print(f"\nWrote {OUT}/FF1_fringe_mode_comparison.md "
          f"and FF1_cosy_fringe_convergence.png")


if __name__ == "__main__":
    main()
