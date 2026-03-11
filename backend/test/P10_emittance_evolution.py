"""P10: Emittance preservation along the transport line.

Tracks ε_n(s) element-by-element, computing both raw and dispersion-corrected
emittance in x and y planes. Shows that geometric emittance is conserved in
drift/quad sections and that apparent growth in the chicane is due to x-δ
dispersion coupling, not real emittance growth.

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

ENERGY = 40
RF_FREQ = 2856e6
SEGMENTS = 118

XLSX = Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx'
RESULTS_DIR = Path(__file__).resolve().parent / 'results' / 'P10'

BASELINE = {
    'bunch_spread': 0.5, 'energy_std_percent': 0.5,
    'h': 5e9, 'epsilon_n': 8, 'x_std': 0.8,
}


def load_currents(json_path):
    with open(json_path) as f:
        data = json.load(f)
    return {int(k): v for k, v in data['currents'].items()}


def setup_beamline(currents, chromatic=False):
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    norm = relat.gamma * relat.beta
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

    return line, relat, norm


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


def compute_emittance(particles):
    """Compute raw and dispersion-corrected emittance for x and y planes."""
    eb = beam()
    _, _, twiss = eb.cal_twiss(particles)

    # Raw emittance (without dispersion correction)
    cov = np.cov(particles, rowvar=False, ddof=1)
    raw_eps = {}
    for i, plane in enumerate(['x', 'y']):
        idx, idx_p = 2*i, 2*i + 1
        var = cov[idx, idx]
        var_p = cov[idx_p, idx_p]
        covar = cov[idx, idx_p]
        eps_sq = var * var_p - covar**2
        raw_eps[plane] = np.sqrt(max(eps_sq, 0))

    return {
        'eps_x_raw': raw_eps['x'],
        'eps_y_raw': raw_eps['y'],
        'eps_x_corr': float(twiss.loc['x', twiss.columns[0]]),
        'eps_y_corr': float(twiss.loc['y', twiss.columns[0]]),
        'beta_x': float(twiss.loc['x', twiss.columns[2]]),
        'alpha_x': float(twiss.loc['x', twiss.columns[1]]),
        'beta_y': float(twiss.loc['y', twiss.columns[2]]),
        'alpha_y': float(twiss.loc['y', twiss.columns[1]]),
        'eta_x': float(twiss.loc['x', twiss.columns[4]]),
        'sigma_x': np.std(particles[:, 0]),
        'sigma_y': np.std(particles[:, 2]),
        'n_particles': len(particles),
    }


def track_evolution(line, particles, norm, use_apertures=False):
    """Track beam element-by-element, recording emittance at each boundary."""
    evolution = []
    s = 0.0
    p = particles.copy()

    # Initial
    em = compute_emittance(p)
    em['s'] = 0.0
    em['element'] = 'START'
    em['element_type'] = ''
    em['eps_nx_raw'] = em['eps_x_raw'] * norm
    em['eps_ny_raw'] = em['eps_y_raw'] * norm
    em['eps_nx_corr'] = em['eps_x_corr'] * norm
    em['eps_ny_corr'] = em['eps_y_corr'] * norm
    evolution.append(em)

    for i, elem in enumerate(line):
        p = np.array(elem.useMatrice(p))
        if use_apertures and hasattr(elem, 'apply_aperture'):
            p = elem.apply_aperture(p)
            if len(p) < 2:
                break

        s += getattr(elem, 'length', 0)

        elem_name = getattr(elem, 'name', f'E{i}')
        elem_type = type(elem).__name__

        em = compute_emittance(p)
        em['s'] = s
        em['element'] = elem_name
        em['element_type'] = elem_type
        em['eps_nx_raw'] = em['eps_x_raw'] * norm
        em['eps_ny_raw'] = em['eps_y_raw'] * norm
        em['eps_nx_corr'] = em['eps_x_corr'] * norm
        em['eps_ny_corr'] = em['eps_y_corr'] * norm
        evolution.append(em)

    return evolution


def identify_chicane_region(line):
    """Find s-range of chicane dipoles."""
    s = 0.0
    dipole_s = []
    for elem in line:
        length_m = getattr(elem, 'length', 0)
        etype = type(elem).__name__
        if 'dipole' in etype.lower():
            dipole_s.append((s, s + length_m, etype))
        s += length_m
    if dipole_s:
        return dipole_s[0][0], dipole_s[-1][1]
    return None, None


def plot_emittance(evolution, norm, outdir, chicane_s=None):
    s_arr = [e['s'] for e in evolution]

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

    # Panel 1: Normalized emittance (x)
    ax = axes[0]
    ax.plot(s_arr, [e['eps_nx_raw'] for e in evolution],
            'b-', label=r'$\varepsilon_{n,x}$ raw', alpha=0.8)
    ax.plot(s_arr, [e['eps_nx_corr'] for e in evolution],
            'r-', label=r'$\varepsilon_{n,x}$ dispersion-corrected', alpha=0.8)
    ax.set_ylabel(r'$\varepsilon_{n,x}$ ($\pi\cdot$mm$\cdot$mrad)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel 2: Normalized emittance (y)
    ax = axes[1]
    ax.plot(s_arr, [e['eps_ny_raw'] for e in evolution],
            'b-', label=r'$\varepsilon_{n,y}$ raw', alpha=0.8)
    ax.plot(s_arr, [e['eps_ny_corr'] for e in evolution],
            'r-', label=r'$\varepsilon_{n,y}$ dispersion-corrected', alpha=0.8)
    ax.set_ylabel(r'$\varepsilon_{n,y}$ ($\pi\cdot$mm$\cdot$mrad)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel 3: Dispersion and beam sizes
    ax = axes[2]
    ax.plot(s_arr, [e['eta_x'] * 1e3 for e in evolution],
            'g-', label=r'$\eta_x$ (mm)')
    ax2 = ax.twinx()
    ax2.plot(s_arr, [e['sigma_x'] for e in evolution],
             'b--', alpha=0.6, label=r'$\sigma_x$ (mm)')
    ax2.plot(s_arr, [e['sigma_y'] for e in evolution],
             'r--', alpha=0.6, label=r'$\sigma_y$ (mm)')
    ax.set_ylabel(r'$\eta_x$ (mm)')
    ax2.set_ylabel(r'$\sigma$ (mm)')
    ax.set_xlabel('s (m)')
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc='upper left')
    ax.grid(True, alpha=0.3)

    if chicane_s and chicane_s[0] is not None:
        for a in axes:
            a.axvspan(chicane_s[0], chicane_s[1], alpha=0.08, color='orange',
                      label='Chicane')

    fig.suptitle('P10: Emittance Evolution Along Transport Line', fontsize=12)
    fig.tight_layout()

    for ext in ['pdf', 'png']:
        fig.savefig(outdir / f'emittance_evolution.{ext}', dpi=150)
    plt.close(fig)


def print_summary(evolution, norm):
    print(f"\n{'=' * 95}")
    print("  EMITTANCE EVOLUTION SUMMARY")
    print(f"{'=' * 95}")

    e0 = evolution[0]
    ef = evolution[-1]

    print(f"\n  {'':>10} {'ε_n,x raw':>12} {'ε_n,x corr':>12} "
          f"{'ε_n,y raw':>12} {'ε_n,y corr':>12} {'η_x (mm)':>10} {'N':>6}")
    print(f"  {'-' * 75}")

    for label, e in [('Initial', e0), ('Final', ef)]:
        print(f"  {label:>10} {e['eps_nx_raw']:>12.4f} {e['eps_nx_corr']:>12.4f} "
              f"{e['eps_ny_raw']:>12.4f} {e['eps_ny_corr']:>12.4f} "
              f"{e['eta_x']*1e3:>10.3f} {e['n_particles']:>6}")

    # Conservation check
    dx_raw = (ef['eps_nx_raw'] - e0['eps_nx_raw']) / e0['eps_nx_raw'] * 100
    dx_corr = (ef['eps_nx_corr'] - e0['eps_nx_corr']) / e0['eps_nx_corr'] * 100
    dy_raw = (ef['eps_ny_raw'] - e0['eps_ny_raw']) / e0['eps_ny_raw'] * 100
    dy_corr = (ef['eps_ny_corr'] - e0['eps_ny_corr']) / e0['eps_ny_corr'] * 100

    print(f"\n  Relative change (initial → final):")
    print(f"    ε_n,x raw:  {dx_raw:+.2f}%")
    print(f"    ε_n,x corr: {dx_corr:+.2f}%")
    print(f"    ε_n,y raw:  {dy_raw:+.2f}%")
    print(f"    ε_n,y corr: {dy_corr:+.2f}%")

    # Peak dispersion
    max_eta = max(abs(e['eta_x']) for e in evolution)
    max_eta_s = [e['s'] for e in evolution if abs(e['eta_x']) == max_eta][0]
    print(f"\n  Peak |η_x| = {max_eta*1e3:.2f} mm at s = {max_eta_s:.3f} m")

    # Peak raw emittance (x)
    max_enx = max(e['eps_nx_raw'] for e in evolution)
    max_enx_s = [e['s'] for e in evolution if e['eps_nx_raw'] == max_enx][0]
    print(f"  Peak ε_n,x(raw) = {max_enx:.4f} at s = {max_enx_s:.3f} m")


def main():
    import argparse
    parser = argparse.ArgumentParser(description='P10: Emittance evolution')
    parser.add_argument('--particles', type=int, default=500)
    parser.add_argument('--currents', type=str,
                        default='results/felsim_nm_warm.json')
    parser.add_argument('--chromatic', action='store_true',
                        help='Use chromatic transport')
    parser.add_argument('--apertures', action='store_true',
                        help='Enable aperture tracking')
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    currents_path = Path(__file__).resolve().parent / args.currents
    if not currents_path.exists():
        alt = Path(__file__).resolve().parent / 'results' / 'seminar' / 'optimized_currents.json'
        if alt.exists():
            currents_path = alt
        else:
            print(f"Currents file not found: {currents_path}")
            sys.exit(1)

    print(f"Loading currents from: {currents_path}")
    currents = load_currents(currents_path)

    line, relat, norm = setup_beamline(currents, chromatic=args.chromatic)
    particles = generate_beam(n=args.particles)
    chicane_s = identify_chicane_region(line)

    mode = 'chromatic' if args.chromatic else 'achromatic'
    print(f"\nTracking {args.particles} particles through {len(line)} elements "
          f"({mode}, apertures={'ON' if args.apertures else 'OFF'})...")

    t0 = time.time()
    evolution = track_evolution(line, particles, norm,
                                use_apertures=args.apertures)
    wall_s = time.time() - t0
    print(f"  Done in {wall_s:.1f}s ({len(evolution)} checkpoints)")

    print_summary(evolution, norm)

    # Chicane info
    if chicane_s[0] is not None:
        print(f"\n  Chicane region: s ∈ [{chicane_s[0]:.3f}, {chicane_s[1]:.3f}] m")
        chicane_evo = [e for e in evolution
                       if chicane_s[0] <= e['s'] <= chicane_s[1]]
        if chicane_evo:
            max_raw = max(e['eps_nx_raw'] for e in chicane_evo)
            print(f"  Peak ε_n,x(raw) in chicane: {max_raw * norm:.4f} π·mm·mrad")

    # Output directory with mode suffix
    outdir = RESULTS_DIR / mode
    outdir.mkdir(parents=True, exist_ok=True)

    # Plots
    plot_emittance(evolution, norm, outdir, chicane_s)
    print(f"\n  Plots saved to {outdir}/")

    # JSON summary (sample key points, not all 118)
    sample_indices = [0, len(evolution)//4, len(evolution)//2,
                      3*len(evolution)//4, len(evolution)-1]
    summary = {
        'mode': mode,
        'apertures': args.apertures,
        'n_particles': args.particles,
        'currents_file': str(currents_path.name),
        'norm_factor': norm,
        'chicane_s': list(chicane_s) if chicane_s[0] else None,
        'n_checkpoints': len(evolution),
        'wall_s': wall_s,
        'samples': [evolution[i] for i in sample_indices
                     if i < len(evolution)],
        'initial': evolution[0],
        'final': evolution[-1],
    }
    with open(outdir / 'summary.json', 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"  Summary saved to {outdir / 'summary.json'}")


if __name__ == '__main__':
    main()
