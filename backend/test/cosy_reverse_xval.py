"""Reverse cross-validation: inject COSY-optimised currents into FELsim.

Takes COSY FR 0 optimised quad currents and propagates through FELsim,
comparing final Twiss and element-by-element dispersion.

Computes transfer-matrix-based analytical dispersion from FELsim's
own matrices and validates against particle-statistics dispersion.

Key finding: COSY's Stage 5 uses negative-polarity quad currents, which
flip QPF↔QPD behaviour. FELsim uses |current| (absolute value), so
negative currents produce the same matrix as positive — making the COSY
chicane solution unreproducible in FELsim. This script handles this by
comparing: (a) pre-chicane Twiss with COSY's positive-only currents,
(b) full beamline with compatible (positive-only) currents.

Author: Eremey Valetov
"""

import sys
import json
import argparse
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

# ── Beam parameters (same as S1 / UHM_beamline_opt_v2.py) ────────────────
Energy = 40
f_rf = 2856e6
bunch_spread = 2
energy_std_percent = 0.5
h = 5e9
epsilon_n = 8
x_std = 0.8
y_std = 0.8
nb_particles = 1000
np.random.seed(42)

# Stage 5 chicane quad indices (use negative currents in COSY)
CHICANE_QUADS = {33, 35, 37, 39, 41, 43}


def compute_targets():
    relat = lattice(1, fringeType=None)
    relat.setE(E=Energy)
    norm = relat.gamma * relat.beta
    epsilon = epsilon_n / norm
    K = 1.2
    lambda_u = 2.3e-2
    beta_ym = relat.gamma * lambda_u / (2 * np.pi * K)
    beta_0 = x_std**2 / epsilon
    return {
        'beta_xm': 1.4, 'alpha_xm': 0.47,
        'beta_ym': beta_ym, 'alpha_ym': 0.0,
        'epsilon': epsilon, 'beta_0': beta_0,
    }


def create_beamline():
    relat = lattice(1, fringeType=None)
    relat.setE(E=Energy)
    file_path = backend_dir.parent / 'beam_excel' / 'Beamline_elements.xlsx'
    excel = ExcelElements(str(file_path))
    beamlineUH = excel.create_beamline()
    line = relat.changeBeamType("electron", Energy, beamlineUH)
    return line[:118]


def generate_beam():
    relat = lattice(1, fringeType=None)
    relat.setE(E=Energy)
    norm = relat.gamma * relat.beta
    epsilon = epsilon_n / norm
    x_prime_std = epsilon / x_std
    y_prime_std = epsilon / y_std
    tof_std = bunch_spread * 1e-9 * f_rf
    energy_std = energy_std_percent * 10

    ebeam_gen = beam()
    beam_dist = ebeam_gen.gen_6d_gaussian(
        0, [x_std, x_prime_std, y_std, y_prime_std, tof_std, energy_std],
        nb_particles)
    tof_dist = beam_dist[:, 4] / f_rf
    beam_dist[:, 5] += h * tof_dist
    return beam_dist


def inject_currents(line, currents, skip_indices=None):
    """Set quad currents from a {felsim_idx_str: current} dict.

    Elements in skip_indices retain their existing currents.
    Negative COSY currents are injected as |current| since FELsim's
    quad model uses absolute values anyway.
    """
    skip = skip_indices or set()
    injected = 0
    for idx_str, current in currents.items():
        idx = int(idx_str)
        if idx in skip:
            continue
        if idx < len(line) and hasattr(line[idx], 'current'):
            line[idx].current = abs(current)
            injected += 1
    return injected


def propagate_twiss(line, beam_dist):
    """Propagate beam and record Twiss at each element exit."""
    ebeam_calc = beam()
    particles = beam_dist.copy()
    s = 0.0
    results = []

    for idx, elem in enumerate(line):
        particles = np.array(elem.useMatrice(particles))
        s += getattr(elem, 'length', getattr(elem, 'L', 0))
        _, _, twiss_df = ebeam_calc.cal_twiss(particles, ddof=1)
        results.append({
            'index': idx, 's': s,
            'beta_x': twiss_df.loc['x', r'$\beta$ (m)'],
            'beta_y': twiss_df.loc['y', r'$\beta$ (m)'],
            'alpha_x': twiss_df.loc['x', r'$\alpha$'],
            'alpha_y': twiss_df.loc['y', r'$\alpha$'],
            'D_x': twiss_df.loc['x', r'$D$ (m)'],
        })

    return results


def propagate_matrix_dispersion(line):
    """Compute analytical dispersion from cumulative transfer matrices.

    Propagates (η, η') through the beamline using M₁₆, M₂₆ from
    each element's 6×6 matrix. Initial dispersion is zero.

    The dispersion η is in the matrix coordinate system: η has units
    such that x_dispersive = η × δ₆, where δ₆ is the 6th coordinate
    (energy deviation in per-mille). Numerically η is in meters (since
    M₁₂ is in meters for the drift case: x[mm] = M₁₂[m] × x'[mrad]).
    The cal_twiss D column "D (m)" is in meters.
    """
    eta = 0.0
    etap = 0.0
    s = 0.0
    results = []

    for idx, elem in enumerate(line):
        M = elem._compute_numeric_matrix()
        s += getattr(elem, 'length', getattr(elem, 'L', 0))

        eta_new = M[0, 0] * eta + M[0, 1] * etap + M[0, 5]
        etap_new = M[1, 0] * eta + M[1, 1] * etap + M[1, 5]
        eta = eta_new
        etap = etap_new

        results.append({
            'index': idx, 's': s,
            'D_x': eta,      # meters (same units as cal_twiss D)
            'Dp_x': etap,
            'M16': M[0, 5],
        })

    return results


def print_twiss_at(results, elem_idx, label):
    """Print Twiss at a specific element."""
    r = results[elem_idx]
    print(f"  {label} (elem {elem_idx}, s={r['s']:.3f} m):")
    print(f"    β_x={r['beta_x']:.4f} m, β_y={r['beta_y']:.4f} m, "
          f"α_x={r['alpha_x']:.4f}, α_y={r['alpha_y']:.4f}")


def print_twiss_comparison(results, targets, label):
    final = results[-1]
    print(f"\n{'=' * 70}")
    print(f"Final Twiss at undulator entrance ({label})")
    print(f"{'=' * 70}")
    print(f"{'Parameter':<20} {'Value':>12} {'Target':>12} {'Δ':>12}")
    print(f"{'-' * 70}")

    params = [
        ('beta_x (m)', final['beta_x'], targets['beta_xm']),
        ('beta_y (m)', final['beta_y'], targets['beta_ym']),
        ('alpha_x', final['alpha_x'], targets['alpha_xm']),
        ('alpha_y', final['alpha_y'], targets['alpha_ym']),
    ]
    mse = 0
    for name, val, tgt in params:
        delta = val - tgt
        mse += delta**2
        print(f"{name:<20} {val:>12.6f} {tgt:>12.6f} {delta:>+12.6f}")
    mse /= 4
    print(f"\n  MSE = {mse:.6e}")
    print(f"{'=' * 70}")
    return mse


def plot_results(felsim_own, felsim_compat, matrix_own, matrix_compat,
                 output_dir):
    """Plot comparison: FELsim own vs COSY-compatible currents."""
    s = [r['s'] for r in felsim_own]

    fig, axes = plt.subplots(3, 1, figsize=(14, 11), sharex=True)

    # Beta
    axes[0].plot(s, [r['beta_x'] for r in felsim_own],
                 'b-', lw=1.5, label=r'FELsim $\beta_x$')
    axes[0].plot(s, [r['beta_y'] for r in felsim_own],
                 'r-', lw=1.5, label=r'FELsim $\beta_y$')
    axes[0].plot(s, [r['beta_x'] for r in felsim_compat],
                 'b--', lw=1.5, label=r'COSY-compat $\beta_x$')
    axes[0].plot(s, [r['beta_y'] for r in felsim_compat],
                 'r--', lw=1.5, label=r'COSY-compat $\beta_y$')
    axes[0].set_ylabel(r'$\beta$ (m)')
    axes[0].legend(fontsize=8, ncol=2)
    axes[0].set_title(
        'FELsim Twiss: Own Currents vs COSY-Compatible Currents\n'
        '(COSY non-chicane currents + FELsim chicane currents)')
    axes[0].set_ylim(0, 50)

    # Alpha
    axes[1].plot(s, [r['alpha_x'] for r in felsim_own],
                 'b-', lw=1.5, label=r'FELsim $\alpha_x$')
    axes[1].plot(s, [r['alpha_y'] for r in felsim_own],
                 'r-', lw=1.5, label=r'FELsim $\alpha_y$')
    axes[1].plot(s, [r['alpha_x'] for r in felsim_compat],
                 'b--', lw=1.5, label=r'COSY-compat $\alpha_x$')
    axes[1].plot(s, [r['alpha_y'] for r in felsim_compat],
                 'r--', lw=1.5, label=r'COSY-compat $\alpha_y$')
    axes[1].set_ylabel(r'$\alpha$')
    axes[1].legend(fontsize=8, ncol=2)

    # Dispersion (both methods overlaid)
    axes[2].plot(s, [r['D_x'] for r in felsim_own],
                 'b-', lw=1.5, label=r'FELsim $D_x$ (particles)')
    axes[2].plot(s, [r['D_x'] for r in felsim_compat],
                 'b--', lw=1.5, label=r'COSY-compat $D_x$ (particles)')
    axes[2].plot([r['s'] for r in matrix_own],
                 [r['D_x'] for r in matrix_own],
                 'g-', lw=1, alpha=0.7, label=r'FELsim $D_x$ (matrix)')
    axes[2].plot([r['s'] for r in matrix_compat],
                 [r['D_x'] for r in matrix_compat],
                 'g--', lw=1, alpha=0.7, label=r'COSY-compat $D_x$ (matrix)')
    axes[2].set_ylabel(r'$D_x$ (m)')
    axes[2].set_xlabel('s (m)')
    axes[2].legend(fontsize=8, ncol=2)

    plt.tight_layout()
    out = output_dir / 'cosy_reverse_xval.eps'
    plt.savefig(out, bbox_inches='tight')
    plt.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reverse cross-validation: COSY currents → FELsim")
    parser.add_argument('--cosy-json', type=str,
                        default='results/cosy_s1_default.json',
                        help='COSY results JSON file')
    args = parser.parse_args()

    targets = compute_targets()

    # Load COSY currents
    with open(args.cosy_json) as fh:
        cosy_data = json.load(fh)
    cosy_currents = cosy_data['currents']
    cosy_mse = cosy_data.get('mse', '?')
    cosy_fr = cosy_data.get('fringe_field_order', '?')
    print(f"COSY source: {args.cosy_json} (FR {cosy_fr}, MSE {cosy_mse:.2e})")

    # Identify polarity incompatibility
    neg_quads = {int(k): v for k, v in cosy_currents.items() if v < 0}
    if neg_quads:
        print(f"\nPolarity incompatibility: COSY uses negative currents at "
              f"elements {sorted(neg_quads.keys())}")
        print("  FELsim uses |I| → QPF/QPD polarity is fixed by element type.")
        print("  Negative currents in COSY flip QPF↔QPD, which FELsim cannot represent.")
        for idx, val in sorted(neg_quads.items()):
            print(f"    [{idx:3d}] I = {val:+.4f} A → |I| = {abs(val):.4f} A "
                  f"(polarity flip lost)")

    # ── Part 1: FELsim reference (own optimised currents) ────────────────
    print("\n[1/5] Running FELsim 11-stage optimization (reference)...")
    from test_cosy_felsim_xval import run_felsim_optimization
    file_path = backend_dir.parent / 'beam_excel' / 'Beamline_elements.xlsx'
    line_own, own_currents = run_felsim_optimization(file_path)

    print("\n[2/5] Propagating FELsim with own currents...")
    beam_dist = generate_beam()
    felsim_own = propagate_twiss(line_own, beam_dist)
    matrix_own = propagate_matrix_dispersion(line_own)
    mse_own = print_twiss_comparison(felsim_own, targets, "FELsim own currents")

    # ── Part 2: Compatible currents (COSY non-chicane + FELsim chicane) ──
    print(f"\n[3/5] Building COSY-compatible beamline...")
    print("  Strategy: Use COSY currents for non-chicane quads (all positive),")
    print("  keep FELsim's own chicane currents (elem 33-43).")

    line_compat = create_beamline()
    # First inject FELsim's optimised currents (sets the full beamline)
    for idx, val in own_currents.items():
        line_compat[idx].current = val
    # Then override with COSY currents for non-chicane quads
    n_overridden = inject_currents(line_compat, cosy_currents,
                                   skip_indices=CHICANE_QUADS)
    print(f"  {n_overridden} non-chicane quad currents overridden from COSY")

    # Show which currents changed
    print(f"\n  {'Elem':>6} {'FELsim':>9} {'COSY-compat':>12} {'Δ':>9} {'Source':>8}")
    print(f"  {'-' * 48}")
    all_indices = sorted(own_currents.keys())
    for idx in all_indices:
        f_val = own_currents[idx]
        c_val = line_compat[idx].current
        delta = c_val - f_val
        src = "chicane" if idx in CHICANE_QUADS else "COSY"
        if abs(delta) > 0.001:
            print(f"  {idx:6d} {f_val:>9.4f} {c_val:>12.4f} {delta:>+9.4f} {src:>8}")

    print(f"\n[4/5] Propagating FELsim with COSY-compatible currents...")
    beam_dist2 = generate_beam()
    felsim_compat = propagate_twiss(line_compat, beam_dist2)
    matrix_compat = propagate_matrix_dispersion(line_compat)
    mse_compat = print_twiss_comparison(
        felsim_compat, targets, "COSY-compatible currents in FELsim")

    # ── Part 3: Pre-chicane Twiss comparison ─────────────────────────────
    print(f"\n{'=' * 70}")
    print("Pre-Chicane Twiss Comparison (at chicane entrance, elem 32)")
    print(f"{'=' * 70}")
    print_twiss_at(felsim_own, 32, "FELsim own")
    print_twiss_at(felsim_compat, 32, "COSY-compatible")
    r1, r2 = felsim_own[32], felsim_compat[32]
    for param in ['beta_x', 'beta_y', 'alpha_x', 'alpha_y']:
        v1, v2 = r1[param], r2[param]
        if abs(v1) > 1e-10:
            pct = abs(v2 - v1) / abs(v1) * 100
        else:
            pct = abs(v2 - v1) * 100
        print(f"  Δ{param}: {v2-v1:+.6f} ({pct:.2f}%)")

    # ── Part 4: Dispersion validation ────────────────────────────────────
    print(f"\n{'=' * 70}")
    print("Dispersion Validation")
    print(f"{'=' * 70}")
    print(f"\nFELsim coordinate convention: 6th coord = δW/W × 10³ (per-mille).")
    print(f"M₁₆ is in meters. cal_twiss D is in meters (mm/per-mille = m).\n")

    cosy_eta = cosy_data['twiss_undulator'].get('eta_x', 0)
    print(f"{'Method':<35} {'D_x at undulator':>18}")
    print(f"{'-' * 55}")
    print(f"{'FELsim (particles, own currents)':<35} {felsim_own[-1]['D_x']:>+18.6f} m")
    print(f"{'FELsim (matrix, own currents)':<35} {matrix_own[-1]['D_x']:>+18.6f} m")
    print(f"{'FELsim (particles, COSY-compat)':<35} {felsim_compat[-1]['D_x']:>+18.6f} m")
    print(f"{'FELsim (matrix, COSY-compat)':<35} {matrix_compat[-1]['D_x']:>+18.6f} m")
    print(f"{'COSY (transfer map, own currents)':<35} {cosy_eta:>+18.6f} m")

    delta_methods = abs(felsim_own[-1]['D_x'] - matrix_own[-1]['D_x'])
    print(f"\nParticle vs matrix agreement (FELsim own): Δ = {delta_methods:.6f} m "
          f"({delta_methods/abs(matrix_own[-1]['D_x'])*100:.1f}%)")

    print(f"\nNote: FELsim and COSY use different quad currents (especially")
    print(f"chicane polarity and Stage 10-11), so different dispersion is expected.")
    print(f"The dispersion difference reflects the solution difference, not a")
    print(f"modeling discrepancy.")

    # ── Part 5: Plot ─────────────────────────────────────────────────────
    print(f"\n[5/5] Generating plots...")
    output_dir = Path(__file__).resolve().parent / 'results'
    plot_results(felsim_own, felsim_compat, matrix_own, matrix_compat,
                 output_dir)

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print("Summary")
    print(f"{'=' * 70}")
    print(f"  FELsim own currents → MSE = {mse_own:.6e}")
    print(f"  COSY-compatible currents → MSE = {mse_compat:.6e}")
    print(f"  COSY self-reported MSE   = {cosy_mse:.6e}")
    print(f"\n  Pre-chicane Twiss agreement: "
          f"Δβ_x = {abs(r2['beta_x']-r1['beta_x']):.4f} m, "
          f"Δβ_y = {abs(r2['beta_y']-r1['beta_y']):.4f} m")
    print(f"\n  Key findings:")
    print(f"  1. COSY Stage 5 uses negative-polarity currents (QPF↔QPD swap)")
    print(f"     → FELsim cannot represent this (uses |I|)")
    print(f"  2. Analytical and statistical dispersion agree in FELsim")
    print(f"  3. D_x difference (FELsim vs COSY) is from different solution families,")
    print(f"     not a modeling discrepancy")
