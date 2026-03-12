"""P11: FELsim fringe field impact on Twiss matching.

Quantifies the effect of dipole wedge (DPW) fringe correction on beam optics
by comparing tracking with and without the triangle-model phi correction
in the y-plane edge kick.

FELsim fringe architecture:
  - DPW edge kick always includes phi = (le/6)·h·(1+sin²η)/cosη (triangle model)
  - fringeType parameter ('decay', [[x],[y]]) only affects field profile
    visualization — the fringeField class uses drift-space transfer matrices
  - Quadrupole fringe: not modeled

This study computes beam transport with the phi correction removed from DPW
elements, quantifying the fringe effect within FELsim's first-order model.

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
from beamline import lattice, dipole_wedge
from excelElements import ExcelElements

ENERGY = 40
RF_FREQ = 2856e6
SEGMENTS = 118

XLSX = Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx'
RESULTS_DIR = Path(__file__).resolve().parent / 'results' / 'P11'

BASELINE = {
    'bunch_spread': 0.5, 'energy_std_percent': 0.5,
    'h': 5e9, 'epsilon_n': 8, 'x_std': 0.8,
}

TWISS_TARGET = {
    'beta_x': 1.4, 'alpha_x': 0.47,
    'beta_y': 0.24, 'alpha_y': 0.0,
}


def load_currents(json_path):
    with open(json_path) as f:
        data = json.load(f)
    return {int(k): v for k, v in data['currents'].items()}


def setup_beamline(currents):
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    norm = relat.gamma * relat.beta
    excel = ExcelElements(XLSX)
    beamlineUH = excel.create_beamline()
    line = relat.changeBeamType("electron", ENERGY, beamlineUH)[:SEGMENTS]

    for idx, current in currents.items():
        if 0 <= idx < len(line) and hasattr(line[idx], 'current'):
            line[idx].current = current

    return line, relat, norm


def dpw_matrix_nophi(elem):
    """Compute DPW transfer matrix with phi=0."""
    l = elem.length
    a = elem.angle
    if abs(elem.dipole_angle) < 1e-14:
        R = np.inf
    else:
        R = elem.dipole_length / (abs(elem.dipole_angle) * np.pi / 180)
    eta = a * np.pi / 180
    Tx = np.tan(eta)
    Ty = np.tan(eta)  # phi=0: Ty = tan(eta) instead of tan(eta - phi)
    M56 = -elem.f * (l / (elem.C * elem.beta * elem.gamma * (elem.gamma + 1)))
    return np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [Tx / R, 1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -Ty / R, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 1.0, M56],
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    ], dtype=np.float64)


def track_nophi(line, particles):
    """Track particles, substituting DPW matrices with phi=0 versions."""
    p = particles.copy()
    for elem in line:
        if isinstance(elem, dipole_wedge):
            mat = dpw_matrix_nophi(elem)
            p = (mat @ p.T).T
        else:
            p = np.array(elem.useMatrice(p))
    return p


def generate_beam(n=500, seed=42):
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    norm = relat.gamma * relat.beta
    epsilon = BASELINE['epsilon_n'] / norm
    x_std = BASELINE['x_std']
    x_prime_std = epsilon / x_std
    tof_std = BASELINE['bunch_spread'] * 1e-9 * RF_FREQ
    energy_std = BASELINE['energy_std_percent'] * 10

    np.random.seed(seed)
    eb = beam()
    dist = eb.gen_6d_gaussian(
        0, [x_std, x_prime_std, x_std, x_prime_std, tof_std, energy_std], n)
    tof_s = dist[:, 4] / RF_FREQ
    dist[:, 5] += BASELINE['h'] * tof_s
    return dist


def extract_twiss(particles):
    eb = beam()
    _, _, twiss = eb.cal_twiss(particles)
    return {
        'beta_x': float(twiss.loc['x', twiss.columns[2]]),
        'alpha_x': float(twiss.loc['x', twiss.columns[1]]),
        'beta_y': float(twiss.loc['y', twiss.columns[2]]),
        'alpha_y': float(twiss.loc['y', twiss.columns[1]]),
    }


def compute_mse(twiss_computed, twiss_target):
    residuals = [twiss_computed[k] - twiss_target[k]
                 for k in ['beta_x', 'alpha_x', 'beta_y', 'alpha_y']]
    return np.mean(np.array(residuals) ** 2)


def analyze_dpw_elements(line):
    """Extract phi correction details for each DPW element."""
    dpw_info = []
    s = 0.0
    for i, elem in enumerate(line):
        if isinstance(elem, dipole_wedge):
            eta_deg = elem.angle
            eta = eta_deg * np.pi / 180
            if abs(elem.dipole_angle) < 1e-14:
                s += getattr(elem, 'length', 0)
                continue
            R = elem.dipole_length / (abs(elem.dipole_angle) * np.pi / 180)
            h = 1.0 / R
            le = elem.length
            phi = (le / 6.0) * h * (1 + np.sin(eta) ** 2) / np.cos(eta)
            M43_with = -np.tan(eta - phi) / R
            M43_without = -np.tan(eta) / R
            # For zero-angle entries, delta% is meaningless (0/0)
            if abs(M43_without) > 1e-12:
                delta_pct = (M43_with - M43_without) / abs(M43_without) * 100
            else:
                delta_pct = None  # pure fringe kick (no hard-edge equivalent)
            dpw_info.append({
                'index': i,
                'name': getattr(elem, 'name', f'DPW_{i}'),
                's': s,
                'eta_deg': eta_deg,
                'phi_rad': phi,
                'phi_deg': phi * 180 / np.pi,
                'M43_with_phi': M43_with,
                'M43_without_phi': M43_without,
                'M43_delta_pct': delta_pct,
                'R_m': R,
            })
        s += getattr(elem, 'length', 0)
    return dpw_info


def plot_comparison(dpw_info, twiss_phi, twiss_nophi, outdir):
    fig, axes = plt.subplots(2, 1, figsize=(10, 7))

    # Panel 1: M43 with vs without phi
    ax = axes[0]
    names = [d['name'] for d in dpw_info]
    m43_with = [d['M43_with_phi'] for d in dpw_info]
    m43_without = [d['M43_without_phi'] for d in dpw_info]
    x = np.arange(len(names))
    w = 0.35
    ax.bar(x - w/2, m43_with, w, label=r'With $\phi$ (standard)', color='steelblue')
    ax.bar(x + w/2, m43_without, w, label=r'Without $\phi$', color='coral')
    ax.set_ylabel(r'$M_{43}$ (m$^{-1}$)')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=7)
    ax.legend(fontsize=9)
    ax.set_title('DPW y-plane edge kick: effect of triangle-model fringe correction')
    ax.grid(True, alpha=0.3)

    # Panel 2: Twiss comparison
    ax = axes[1]
    params = ['beta_x', 'alpha_x', 'beta_y', 'alpha_y']
    labels = [r'$\beta_x$ (m)', r'$\alpha_x$', r'$\beta_y$ (m)', r'$\alpha_y$']
    targets = [TWISS_TARGET[p] for p in params]
    vals_phi = [twiss_phi[p] for p in params]
    vals_nophi = [twiss_nophi[p] for p in params]

    x = np.arange(len(params))
    w = 0.25
    ax.bar(x - w, targets, w, label='Target', color='gray', alpha=0.5)
    ax.bar(x, vals_phi, w, label=r'With $\phi$', color='steelblue', alpha=0.8)
    ax.bar(x + w, vals_nophi, w, label=r'Without $\phi$', color='coral', alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.legend(fontsize=9)
    ax.set_title('Final Twiss parameters (fixed currents)')
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(outdir / f'fringe_comparison.{ext}', dpi=150)
    plt.close(fig)


def main():
    import argparse
    parser = argparse.ArgumentParser(description='P11: Fringe field impact')
    parser.add_argument('--particles', type=int, default=500)
    parser.add_argument('--currents', type=str,
                        default='results/felsim_nm_warm.json')
    parser.add_argument('--plots-only', action='store_true',
                        help='Regenerate plots from cached summary.json')
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if args.plots_only:
        summary_path = RESULTS_DIR / 'summary.json'
        if not summary_path.exists():
            print(f"No cached results at {summary_path}")
            sys.exit(1)
        with open(summary_path) as f:
            summary = json.load(f)
        plot_comparison(
            summary['dpw_elements'],
            summary['with_phi']['twiss'],
            summary['without_phi']['twiss'],
            RESULTS_DIR,
        )
        print(f"Plots regenerated from {summary_path}")
        return

    currents_path = Path(__file__).resolve().parent / args.currents
    if not currents_path.exists():
        print(f"Currents file not found: {currents_path}")
        sys.exit(1)

    currents = load_currents(currents_path)

    # Part A: DPW element analysis
    print("=" * 80)
    print("  P11: FRINGE FIELD IMPACT ANALYSIS")
    print("=" * 80)

    print("\n  Part A: DPW phi correction analysis")
    print("  " + "-" * 60)

    line, _, _ = setup_beamline(currents)
    dpw_info = analyze_dpw_elements(line)

    print(f"\n  Found {len(dpw_info)} DPW elements (excluding zero-angle):\n")
    print(f"  {'Name':>16} {'η (°)':>8} {'φ (°)':>10} "
          f"{'M43 w/ φ':>12} {'M43 w/o φ':>12} {'Δ (%)':>8}")
    print(f"  {'-' * 70}")

    for d in dpw_info:
        delta_str = f"{d['M43_delta_pct']:>8.2f}" if d['M43_delta_pct'] is not None else "   pure"
        print(f"  {d['name']:>16} {d['eta_deg']:>8.2f} {d['phi_deg']:>10.4f} "
              f"{d['M43_with_phi']:>12.6f} {d['M43_without_phi']:>12.6f} "
              f"{delta_str}")

    # Part B: Twiss comparison with fixed currents
    print("\n\n  Part B: Twiss match with fixed currents")
    print("  " + "-" * 60)

    particles = generate_beam(n=args.particles)

    # Standard tracking (with phi)
    line_std, _, _ = setup_beamline(currents)
    p_std = particles.copy()
    for elem in line_std:
        p_std = np.array(elem.useMatrice(p_std))
    twiss_phi = extract_twiss(p_std)
    mse_phi = compute_mse(twiss_phi, TWISS_TARGET)

    # No-phi tracking
    line_nophi, _, _ = setup_beamline(currents)
    p_nophi = track_nophi(line_nophi, particles)
    twiss_nophi = extract_twiss(p_nophi)
    mse_nophi = compute_mse(twiss_nophi, TWISS_TARGET)

    print(f"\n  {'':>20} {'β_x (m)':>10} {'α_x':>10} "
          f"{'β_y (m)':>10} {'α_y':>10} {'RMS':>12}")
    print(f"  {'-' * 75}")
    print(f"  {'Target':>20} {TWISS_TARGET['beta_x']:>10.4f} "
          f"{TWISS_TARGET['alpha_x']:>10.4f} {TWISS_TARGET['beta_y']:>10.4f} "
          f"{TWISS_TARGET['alpha_y']:>10.4f} {'—':>12}")
    print(f"  {'With phi':>20} {twiss_phi['beta_x']:>10.4f} "
          f"{twiss_phi['alpha_x']:>10.4f} {twiss_phi['beta_y']:>10.4f} "
          f"{twiss_phi['alpha_y']:>10.4f} {np.sqrt(mse_phi):>12.4e}")
    print(f"  {'Without phi':>20} {twiss_nophi['beta_x']:>10.4f} "
          f"{twiss_nophi['alpha_x']:>10.4f} {twiss_nophi['beta_y']:>10.4f} "
          f"{twiss_nophi['alpha_y']:>10.4f} {np.sqrt(mse_nophi):>12.4e}")

    # Difference
    d_beta_y = twiss_nophi['beta_y'] - twiss_phi['beta_y']
    print(f"\n  Δβ_y (no-phi − phi) = {d_beta_y:+.6f} m")
    print(f"  The phi correction shifts β_y by {abs(d_beta_y)*1e3:.2f} mm")

    # Architecture summary
    print("\n\n  FELsim fringe field architecture:")
    print("  " + "-" * 60)
    print("  1. DPW edge kick: always includes triangle-model phi correction")
    print("     phi = (l_e/6) * h * (1 + sin^2(eta)) / cos(eta)")
    print("  2. fringeType parameter: field profile only (drift-space matrix)")
    print("  3. Quadrupole fringe: not modeled in FELsim")
    print("  4. Compare to COSY (P8): FR=3 captures full fringe including")
    print("     quad fringe and higher-order fringe contributions")

    # Plots
    plot_comparison(dpw_info, twiss_phi, twiss_nophi, RESULTS_DIR)
    print(f"\n  Plots saved to {RESULTS_DIR}/")

    # JSON summary
    summary = {
        'dpw_elements': dpw_info,
        'with_phi': {'twiss': twiss_phi, 'mse': mse_phi,
                     'rms': float(np.sqrt(mse_phi))},
        'without_phi': {'twiss': twiss_nophi, 'mse': mse_nophi,
                        'rms': float(np.sqrt(mse_nophi))},
        'delta_beta_y_m': d_beta_y,
        'n_particles': args.particles,
        'currents_file': str(currents_path.name),
        'architecture': {
            'dpw_phi': 'always active, triangle model',
            'fringeType_param': 'field profile only, drift matrix',
            'quad_fringe': 'not modeled',
        },
    }
    with open(RESULTS_DIR / 'summary.json', 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"  Summary saved to {RESULTS_DIR / 'summary.json'}")


if __name__ == '__main__':
    main()
