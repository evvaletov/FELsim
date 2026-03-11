"""S6: Bunch Length Sensitivity Sweep (0.1–2 ps)

Tests whether transverse Twiss matching depends on bunch length.
S9 predicts no dependence (linear transfer matrices decouple transverse
and longitudinal planes). This script confirms numerically.

Two parameter sets:
  (a) baseline: σ_E=0.5%, h=5e9
  (b) emittance-conservation scaled: σ_E=2%, h=20e9

Author: Eremey Valetov
"""

import sys
import time
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from UHM_beamline_opt_05ps_params import (
    run_optimization, run_scan, BASELINE, QUAD_INDICES, MSE_THRESHOLDS,
)

OUTDIR = Path(__file__).resolve().parent / 'results' / 'S6'

BUNCH_LENGTHS = np.linspace(0.1, 2.0, 15)  # ps

CONFIGS = {
    'baseline': {
        'energy_std_percent': 0.5,
        'h': 5e9,
        'label': r'Baseline ($\sigma_E$=0.5%, $h$=5e9)',
    },
    'emittance_conserved': {
        'energy_std_percent': 2.0,
        'h': 20e9,
        'label': r'Emittance-conserved ($\sigma_E$=2%, $h$=20e9)',
    },
}


def run_bunch_sweep(config_name, nb_particles=500, seed=42):
    cfg = CONFIGS[config_name]
    csv_path = OUTDIR / f'scan_{config_name}.csv'

    overrides = {k: v for k, v in cfg.items() if k != 'label'}

    rows = run_scan(
        scan_name=f's6_{config_name}',
        param_name='bunch_spread',
        param_values=BUNCH_LENGTHS,
        outdir=OUTDIR,
        nb_particles=nb_particles,
        seed=seed,
        **overrides,
    )
    return rows


def plot_results():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for config_name, cfg in CONFIGS.items():
        csv_path = OUTDIR / f'scan_s6_{config_name}.csv'
        if not csv_path.exists():
            continue

        with open(csv_path) as f:
            rows = list(csv.DictReader(f))

        bl = [float(r['param_value']) for r in rows]
        mse = [float(r['mse']) for r in rows]
        beta_x = [float(r['beta_x']) for r in rows]
        beta_y = [float(r['beta_y']) for r in rows]

        ax1.semilogy(bl, mse, 'o-', label=cfg['label'], markersize=5)
        ax2.plot(bl, beta_x, 'o-', label=cfg['label'] + r' $\beta_x$', markersize=4)
        ax2.plot(bl, beta_y, 's--', label=cfg['label'] + r' $\beta_y$', markersize=4)

    # Threshold lines
    for name, thresh in MSE_THRESHOLDS.items():
        ax1.axhline(thresh, ls=':', color='gray', alpha=0.5)
        ax1.text(0.15, thresh * 1.3, name, fontsize=8, color='gray')

    ax1.set_xlabel('Bunch length (ps)')
    ax1.set_ylabel('MSE')
    ax1.set_title('S6: MSE vs Bunch Length')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    ax2.axhline(1.4, ls=':', color='blue', alpha=0.4, label=r'$\beta_x$ target')
    ax2.axhline(0.2418, ls=':', color='red', alpha=0.4, label=r'$\beta_y$ target')
    ax2.set_xlabel('Bunch length (ps)')
    ax2.set_ylabel(r'$\beta$ (m)')
    ax2.set_title('S6: Twiss Parameters vs Bunch Length')
    ax2.legend(fontsize=7, ncol=2)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'S6_bunch_length_sensitivity.{ext}', dpi=150)
    print(f"Saved: S6_bunch_length_sensitivity.{{eps,png}}")
    plt.close(fig)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    print("S6: Bunch Length Sensitivity (0.1–2 ps)")
    print(f"  15 points × 2 configs = 30 optimizations")

    for config_name in CONFIGS:
        print(f"\n── {config_name} ──")
        run_bunch_sweep(config_name)

    plot_results()

    print("\n" + "=" * 60)
    print("  S6 Complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
