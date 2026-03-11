"""P9: Chromaticity analysis — Twiss sensitivity to energy deviation.

Sweeps energy deviation δ from -3% to +3% and computes Twiss parameters at
the undulator entrance for each δ using FELsim's chromatic transport matrices.
Compares chromatic vs achromatic transport to quantify momentum-dependent
matching degradation.

Author: Eremey Valetov
"""

import sys
import json
import time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

from ebeam import beam
from beamline import lattice
from excelElements import ExcelElements

ENERGY = 40          # MeV
RF_FREQ = 2856e6     # Hz
SEGMENTS = 118
K_UND = 1.2
LAMBDA_U = 2.3e-2   # m

XLSX = Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx'
RESULTS_DIR = Path(__file__).resolve().parent / 'results' / 'P9'

BASELINE = {
    'bunch_spread': 0.5, 'energy_std_percent': 0.5,
    'h': 5e9, 'epsilon_n': 8, 'x_std': 0.8,
}

QUAD_INDICES = [1, 3, 10, 16, 18, 20, 27, 33, 35, 37, 39, 41, 43,
                50, 56, 58, 61, 63, 70, 76, 78, 80, 87, 93, 95, 97]


def compute_targets():
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    norm = relat.gamma * relat.beta
    beta_ym = relat.gamma * LAMBDA_U / (2 * np.pi * K_UND)
    epsilon = BASELINE['epsilon_n'] / norm
    return {
        'beta_xm': 1.4, 'alpha_xm': 0.47,
        'beta_ym': beta_ym, 'alpha_ym': 0.0,
        'epsilon': epsilon, 'norm': norm,
    }


def load_currents(json_path):
    with open(json_path) as f:
        data = json.load(f)
    return {int(k): v for k, v in data['currents'].items()}


def setup_beamline(currents, chromatic=False):
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    excel = ExcelElements(XLSX)
    beamlineUH = excel.create_beamline()
    line = relat.changeBeamType("electron", ENERGY, beamlineUH)[:SEGMENTS]

    for idx, current in currents.items():
        if 0 <= idx < len(line) and hasattr(line[idx], 'current'):
            line[idx].current = current

    if chromatic:
        for elem in line:
            if hasattr(elem, 'chromatic'):
                elem.chromatic = True

    return line, relat


def generate_beam(n=500, seed=42, delta_offset=0.0, mono_energetic=False):
    """Generate beam with optional energy offset δ (in %).

    delta_offset shifts all particles by that amount in ΔK/K₀.
    mono_energetic: if True, set σ_E ≈ 0 for clean chromaticity measurement.
    """
    targets = compute_targets()
    epsilon = targets['epsilon']
    x_std = BASELINE['x_std']
    x_prime_std = epsilon / x_std
    tof_std = BASELINE['bunch_spread'] * 1e-9 * RF_FREQ

    if mono_energetic:
        energy_std = 0.5  # 0.05% — small enough for chromaticity, large enough for stable Twiss
    else:
        energy_std = BASELINE['energy_std_percent'] * 10  # % → ‰

    np.random.seed(seed)
    eb = beam()
    dist = eb.gen_6d_gaussian(
        0, [x_std, x_prime_std, x_std, x_prime_std, tof_std, energy_std], n)

    if not mono_energetic:
        tof_s = dist[:, 4] / RF_FREQ
        dist[:, 5] += BASELINE['h'] * tof_s

    # Energy offset (‰ units)
    dist[:, 5] += delta_offset * 10

    return dist


def transport_and_twiss(line, particles):
    eb = beam()
    p = particles.copy()
    for elem in line:
        p = np.array(elem.useMatrice(p))
    _, _, twiss_df = eb.cal_twiss(p)
    return twiss_df, p


def run_delta_sweep(currents, deltas, chromatic=True, n_particles=500,
                    mono_energetic=True):
    """Sweep δ values, return list of result dicts."""
    targets = compute_targets()
    results = []

    for delta in deltas:
        line, _ = setup_beamline(currents, chromatic=chromatic)
        particles = generate_beam(n=n_particles, delta_offset=delta,
                                  mono_energetic=mono_energetic)
        twiss_df, final_p = transport_and_twiss(line, particles)

        bx = twiss_df.loc['x', twiss_df.columns[2]]  # β (m)
        ax = twiss_df.loc['x', twiss_df.columns[1]]   # α
        by = twiss_df.loc['y', twiss_df.columns[2]]
        ay = twiss_df.loc['y', twiss_df.columns[1]]
        ex = twiss_df.loc['x', twiss_df.columns[0]]   # ε
        ey = twiss_df.loc['y', twiss_df.columns[0]]
        dx = twiss_df.loc['x', twiss_df.columns[4]]   # D (m)

        mse = ((bx - targets['beta_xm'])**2 + (ax - targets['alpha_xm'])**2 +
               (by - targets['beta_ym'])**2 + (ay - targets['alpha_ym'])**2) / 4

        results.append({
            'delta_pct': delta,
            'beta_x': float(bx), 'alpha_x': float(ax),
            'beta_y': float(by), 'alpha_y': float(ay),
            'epsilon_x': float(ex), 'epsilon_y': float(ey),
            'eta_x': float(dx),
            'sigma_x': float(np.std(final_p[:, 0])),
            'sigma_y': float(np.std(final_p[:, 2])),
            'mse': float(mse), 'rms': float(np.sqrt(mse)),
            'n_survived': len(final_p),
        })

    return results


def run_aperture_sweep(currents, deltas, n_particles=500):
    """Sweep δ with apertures enabled, tracking transmission."""
    targets = compute_targets()
    results = []

    for delta in deltas:
        line, _ = setup_beamline(currents, chromatic=True)
        particles = generate_beam(n=n_particles, delta_offset=delta)
        n_initial = len(particles)
        p = particles.copy()
        for elem in line:
            p = np.array(elem.useMatrice(p))
            if hasattr(elem, 'apply_aperture'):
                p = elem.apply_aperture(p)

        eb = beam()
        if len(p) >= 2:
            _, _, twiss_df = eb.cal_twiss(p)
            bx = float(twiss_df.loc['x', twiss_df.columns[2]])
            by = float(twiss_df.loc['y', twiss_df.columns[2]])
        else:
            bx = by = float('nan')

        results.append({
            'delta_pct': delta,
            'n_survived': len(p),
            'transmission': len(p) / n_initial,
            'beta_x': bx, 'beta_y': by,
        })

    return results


def print_table(results, targets, label):
    print(f"\n{'=' * 90}")
    print(f"  {label}")
    print(f"{'=' * 90}")
    hdr = (f"{'δ (%)':>7}  {'β_x (m)':>9}  {'α_x':>9}  {'β_y (m)':>9}  "
           f"{'α_y':>9}  {'η_x (m)':>9}  {'σ_x (mm)':>9}  {'RMS':>10}")
    print(hdr)
    print("-" * len(hdr))

    for r in results:
        print(f"{r['delta_pct']:>7.1f}  {r['beta_x']:>9.4f}  {r['alpha_x']:>9.4f}  "
              f"{r['beta_y']:>9.4f}  {r['alpha_y']:>9.4f}  {r['eta_x']:>9.5f}  "
              f"{r['sigma_x']:>9.4f}  {r['rms']:>10.4e}")

    print(f"{'target':>7}  {targets['beta_xm']:>9.4f}  {targets['alpha_xm']:>9.4f}  "
          f"{targets['beta_ym']:>9.4f}  {targets['alpha_ym']:>9.4f}")


def plot_chromaticity(results_A, results_B, targets, outdir, results_C=None):
    deltas = [r['delta_pct'] for r in results_A]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    params = [('beta_x', r'$\beta_x$ (m)', targets['beta_xm']),
              ('alpha_x', r'$\alpha_x$', targets['alpha_xm']),
              ('beta_y', r'$\beta_y$ (m)', targets['beta_ym']),
              ('alpha_y', r'$\alpha_y$', targets['alpha_ym'])]

    for ax, (key, label, target) in zip(axes.flat, params):
        ax.plot(deltas, [r[key] for r in results_A], 'o-',
                label='A: chrom transport, achrom currents', markersize=4)
        ax.plot(deltas, [r[key] for r in results_B], 's--',
                label='B: achrom transport (ref)', markersize=4)
        if results_C:
            ax.plot(deltas, [r[key] for r in results_C], '^-',
                    label='C: chrom transport, chrom currents', markersize=4)
        ax.axhline(target, color='k', ls=':', alpha=0.5, label='Target')
        ax.axvline(0, color='gray', ls=':', alpha=0.3)
        ax.set_xlabel(r'$\delta$ (%)')
        ax.set_ylabel(label)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    fig.suptitle('P9: Twiss Sensitivity to Energy Deviation', fontsize=12)
    fig.tight_layout()

    for ext in ['pdf', 'png']:
        fig.savefig(outdir / f'twiss_vs_delta.{ext}', dpi=150)
    plt.close(fig)

    # RMS plot
    fig2, ax2 = plt.subplots(figsize=(7, 4))
    ax2.semilogy(deltas, [r['rms'] for r in results_A], 'o-',
                 label='A: chrom, achrom currents', markersize=5)
    ax2.semilogy(deltas, [r['rms'] for r in results_B], 's--',
                 label='B: achrom (ref)', markersize=5)
    if results_C:
        ax2.semilogy(deltas, [r['rms'] for r in results_C], '^-',
                     label='C: chrom, chrom currents', markersize=5)
    ax2.axhline(3.2e-2, color='green', ls=':', alpha=0.5, label='Excellent')
    ax2.axhline(1e-1, color='orange', ls=':', alpha=0.5, label='Acceptable')
    ax2.set_xlabel(r'Energy deviation $\delta$ (%)')
    ax2.set_ylabel('RMS Twiss residual')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_title('P9: Matching Quality vs Energy Deviation')
    fig2.tight_layout()

    for ext in ['pdf', 'png']:
        fig2.savefig(outdir / f'rms_vs_delta.{ext}', dpi=150)
    plt.close(fig2)


def main():
    import argparse
    parser = argparse.ArgumentParser(description='P9: Chromaticity analysis')
    parser.add_argument('--particles', type=int, default=500)
    parser.add_argument('--points', type=int, default=13,
                        help='Number of δ points')
    parser.add_argument('--delta-max', type=float, default=3.0,
                        help='Max δ in %%')
    parser.add_argument('--currents', type=str,
                        default='results/felsim_nm_warm.json',
                        help='Path to optimized currents JSON')
    parser.add_argument('--chromatic-currents', type=str,
                        default='results/felsim_chromatic_warm.json',
                        help='Path to chromatic-optimized currents JSON')
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    targets = compute_targets()
    base = Path(__file__).resolve().parent

    currents_path = base / args.currents
    if not currents_path.exists():
        print(f"Currents file not found: {currents_path}")
        sys.exit(1)

    print(f"Achromatic currents: {currents_path.name}")
    currents = load_currents(currents_path)

    chrom_currents_path = base / args.chromatic_currents
    has_chrom_currents = chrom_currents_path.exists()
    if has_chrom_currents:
        print(f"Chromatic  currents: {chrom_currents_path.name}")
        chrom_currents = load_currents(chrom_currents_path)

    deltas = np.linspace(-args.delta_max, args.delta_max, args.points)

    # All sweeps use mono-energetic beams for clean chromaticity measurement
    mono = True

    # ── A: Achromatic currents with chromatic transport ──
    print(f"\nA. Chromatic transport, achromatic currents "
          f"({len(deltas)} pts, {args.particles} particles, mono-energetic)...")
    t0 = time.time()
    results_A = run_delta_sweep(currents, deltas, chromatic=True,
                                n_particles=args.particles, mono_energetic=mono)
    t_A = time.time() - t0
    print_table(results_A, targets, "A: CHROMATIC TRANSPORT, ACHROMATIC CURRENTS")

    # ── B: Achromatic currents with achromatic transport (reference) ──
    print(f"\nB. Achromatic transport, achromatic currents...")
    t0 = time.time()
    results_B = run_delta_sweep(currents, deltas, chromatic=False,
                                n_particles=args.particles, mono_energetic=mono)
    t_B = time.time() - t0
    print_table(results_B, targets, "B: ACHROMATIC TRANSPORT (REFERENCE)")

    # ── C: Chromatic currents with chromatic transport ──
    results_C = None
    if has_chrom_currents:
        print(f"\nC. Chromatic transport, chromatic-optimized currents...")
        t0 = time.time()
        results_C = run_delta_sweep(chrom_currents, deltas, chromatic=True,
                                    n_particles=args.particles, mono_energetic=mono)
        t_C = time.time() - t0
        print_table(results_C, targets, "C: CHROMATIC TRANSPORT, CHROMATIC CURRENTS")

    # Summary
    print(f"\n{'=' * 90}")
    print("  CHROMATICITY SUMMARY")
    print(f"{'=' * 90}")

    idx_0 = len(deltas) // 2
    if idx_0 > 0 and idx_0 < len(deltas) - 1:
        dd = deltas[idx_0 + 1] - deltas[idx_0 - 1]
        for key, label in [('beta_x', 'dβ_x/dδ'), ('beta_y', 'dβ_y/dδ'),
                           ('alpha_x', 'dα_x/dδ'), ('alpha_y', 'dα_y/dδ')]:
            deriv_A = (results_A[idx_0+1][key] - results_A[idx_0-1][key]) / dd
            print(f"  {label:>10} (A, achrom currents) = {deriv_A:>10.4f} /%", end='')
            if results_C:
                deriv_C = (results_C[idx_0+1][key] - results_C[idx_0-1][key]) / dd
                print(f",  (C, chrom currents) = {deriv_C:>10.4f} /%")
            else:
                print()

    # Acceptance bandwidth
    excellent = 3.2e-2
    for label, results in [('A (achrom currents)', results_A),
                           ('B (achrom transport)', results_B),
                           ('C (chrom currents)', results_C)]:
        if results is None:
            continue
        bw = [r['delta_pct'] for r in results if r['rms'] < excellent]
        if bw:
            print(f"  {label} acceptance: δ ∈ [{min(bw):+.1f}%, {max(bw):+.1f}%]")
        else:
            print(f"  {label} acceptance: NONE (no Excellent points)")

    print(f"\n  δ=0 RMS: A={results_A[idx_0]['rms']:.4e}, "
          f"B={results_B[idx_0]['rms']:.4e}", end='')
    if results_C:
        print(f", C={results_C[idx_0]['rms']:.4e}")
    else:
        print()

    # Plots
    plot_chromaticity(results_A, results_B, targets, RESULTS_DIR,
                      results_C=results_C)
    print(f"  Plots saved to {RESULTS_DIR}/")

    # JSON summary
    summary = {
        'targets': targets,
        'currents_file': str(currents_path.name),
        'chromatic_currents_file': str(chrom_currents_path.name) if has_chrom_currents else None,
        'n_particles': args.particles,
        'delta_range_pct': [-args.delta_max, args.delta_max],
        'n_points': args.points,
        'A_chrom_achrom_currents': results_A,
        'B_achrom_reference': results_B,
        'C_chrom_chrom_currents': results_C,
    }
    with open(RESULTS_DIR / 'summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  Summary saved to {RESULTS_DIR / 'summary.json'}")


if __name__ == '__main__':
    main()
