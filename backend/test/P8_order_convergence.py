"""P8: Order-by-order DA convergence study via COSY INFINITY.

Evaluates the UH FEL transport line transfer map at DA orders 1, 2, 3, 5
using fixed quad currents (FELsim-optimized and COSY FR3-optimized).
Quantifies how Twiss parameters, aberration coefficients, and particle
distributions change with map order.

Author: Eremey Valetov
"""

import sys
import json
import math
import time
from pathlib import Path
import numpy as np

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

from cosyAdapter import COSYAdapter
from cosyResultsReader import COSYResultsReader
from cosyOptHelper import parse_beamline_felsim_indexed

XLSX = Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx'
RESULTS_DIR = Path(__file__).resolve().parent / 'results' / 'P8'

Energy = 40     # MeV
epsilon_n = 8   # pi.mm.mrad
x_std = 0.8     # mm
N_PARTICLES = 500
SEED = 42

# Physics constants
MASS_E = 0.51099895  # MeV/c²
GAMMA = (Energy + MASS_E) / MASS_E
BETA_REL = np.sqrt(1 - 1/GAMMA**2)
C_LIGHT = 299792458.0
BETA_C = BETA_REL * C_LIGHT


def compute_targets():
    from beamline import lattice
    relat = lattice(1, fringeType=None)
    relat.setE(E=Energy)
    norm = relat.gamma * relat.beta
    epsilon = epsilon_n / norm

    K = 1.2
    lambda_u = 2.3e-2
    beta_ym = relat.gamma * lambda_u / (2 * np.pi * K)

    return {
        'beta_xm': 1.4, 'alpha_xm': 0.47,
        'beta_ym': beta_ym, 'alpha_ym': 0.0,
        'epsilon': epsilon, 'beta_0': x_std**2 / epsilon,
    }


def load_currents(json_path):
    with open(json_path) as f:
        data = json.load(f)
    return {int(k): v for k, v in data['currents'].items()}


def inject_currents(sim, currents):
    """Set quad currents on the COSY beamline.

    Truncate to [:118] to match UHM_beamline_opt_cosy.py — elements beyond
    118 are diagnostics/UND that break the lattice optics in COSY.
    """
    bl = parse_beamline_felsim_indexed(str(XLSX))[:118]
    for idx, current in currents.items():
        if idx < len(bl):
            bl[idx]['current'] = current
    sim._native_sim.beamline = bl


def generate_beam_felsim(n=N_PARTICLES, seed=SEED):
    """Generate 6D beam in FELsim coordinates."""
    from ebeam import beam
    from beamline import lattice
    relat = lattice(1, fringeType=None)
    relat.setE(E=Energy)
    norm = relat.gamma * relat.beta
    epsilon = epsilon_n / norm

    RF_FREQ = 2856e6
    x_prime_std = epsilon / x_std
    tof_std = 0.5e-12 * 1e-9 * RF_FREQ  # 0.5 ps bunch
    energy_std = 0.5 * 10  # 0.5% → FELsim units

    np.random.seed(seed)
    eb = beam()
    dist = eb.gen_6d_gaussian(
        0, [x_std, x_prime_std, x_std, x_prime_std, tof_std, energy_std], n)
    # Add chirp
    h = 5e9
    tof_s = dist[:, 4] / RF_FREQ  # relative ToF in seconds
    dist[:, 5] += h * tof_s

    return dist


def run_at_order(order, fringe_order, currents, targets, output_subdir,
                 n_particles=N_PARTICLES, seed=SEED):
    """Run COSY at a specific DA order. Returns dict of results."""
    out_dir = str(RESULTS_DIR / output_subdir / f'order_{order}')

    sim = COSYAdapter(
        excel_path=str(XLSX),
        mode='particle_tracking',
        config={'simulation': {
            'order': order, 'dimensions': 3, 'KE': Energy,
            'transfer_matrix_order': order,
        }},
        fringe_field_order=fringe_order,
        debug=False,
    )
    inject_currents(sim, currents)

    # Generate beam
    particles = generate_beam_felsim(n_particles, seed)

    # Track through beamline
    t0 = time.time()
    evolution = sim.collect_evolution(particles, checkpoint_elements='all')
    wall_s = time.time() - t0

    # Read transfer map from fort.99
    reader = COSYResultsReader('results')
    linear_map = reader.read_linear_transfer_map()
    all_orders = reader.read_transfer_map_all_orders(max_order=order)

    # Twiss from transfer map propagation
    beta_0 = targets['beta_0']
    twiss = reader.get_twiss_from_transfer_map(
        initial_twiss_x={'beta': beta_0, 'alpha': 0.0, 'eta': 0.0, 'etap': 0.0},
        initial_twiss_y={'beta': beta_0, 'alpha': 0.0, 'eta': 0.0, 'etap': 0.0},
        include_dispersion=True,
    )

    # Also get COSY's own Twiss (GT MAP output)
    result_file = Path('results') / 'result.txt'
    cosy_twiss = {}
    if result_file.exists():
        try:
            with open(result_file) as f:
                data = json.loads(f.read().replace("'", '"'))
            cosy_twiss = data.get('twiss', {})
        except Exception:
            pass

    # Final beam statistics
    if evolution.s_positions:
        s_final = evolution.s_positions[-1]
        p = evolution.get_particles_at(s_final)
    else:
        p = None
    if p is not None and len(p) >= 2:
        n_survived = p.shape[0]
        sigma_x = np.std(p[:, 0])
        sigma_xp = np.std(p[:, 1])
        sigma_y = np.std(p[:, 2])
        sigma_yp = np.std(p[:, 3])
        sigma_t = np.std(p[:, 4]) if p.shape[1] > 4 else None
    else:
        n_survived = 0 if p is None else len(p)
        sigma_x = sigma_xp = sigma_y = sigma_yp = sigma_t = None

    # MSE vs targets
    bx = twiss.get('beta_x', 0)
    ax = twiss.get('alpha_x', 0)
    by = twiss.get('beta_y', 0)
    ay = twiss.get('alpha_y', 0)
    mse = ((bx - targets['beta_xm'])**2 + (ax - targets['alpha_xm'])**2 +
           (by - targets['beta_ym'])**2 + (ay - targets['alpha_ym'])**2) / 4

    # Count coefficients per order
    n_coeffs = {}
    for o in sorted(all_orders.keys()):
        if o == 0:
            n_coeffs[o] = 1
        elif o == 1:
            n_coeffs[o] = int(np.count_nonzero(all_orders[1]))
        else:
            n_coeffs[o] = len(all_orders[o])

    # Key aberrations
    aberrations = extract_key_aberrations(all_orders)

    # Save fort.99 and result.txt to output directory
    import shutil
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    for fn in ['fort.99', 'result.txt', 'input.fox']:
        src = Path('results') / fn
        if src.exists():
            shutil.copy(src, Path(out_dir) / fn)

    return {
        'order': order,
        'fringe_order': fringe_order,
        'wall_s': wall_s,
        'linear_map': linear_map,
        'all_orders': all_orders,
        'twiss': twiss,
        'cosy_twiss': cosy_twiss,
        'mse': mse,
        'rms': math.sqrt(mse),
        'n_coeffs': n_coeffs,
        'n_particles_in': n_particles,
        'n_survived': n_survived,
        'sigma_x': sigma_x,
        'sigma_xp': sigma_xp,
        'sigma_y': sigma_y,
        'sigma_yp': sigma_yp,
        'sigma_t': sigma_t,
        'aberrations': aberrations,
    }


def extract_key_aberrations(all_orders):
    """Extract physically significant aberration coefficients."""
    ab = {}

    if 2 in all_orders:
        s = all_orders[2]
        for name, key in {
            'T_166': (0, 5, 5),   # x from δ²
            'T_266': (1, 5, 5),   # x' from δ²
            'T_366': (2, 5, 5),   # y from δ²
            'T_466': (3, 5, 5),   # y' from δ²
            'T_566': (4, 5, 5),   # l from δ² (= T566)
            'T_111': (0, 0, 0),   # geometric x³→x
            'T_116': (0, 0, 5),   # chromatic x·δ→x
            'T_126': (0, 1, 5),   # chromatic x'·δ→x
            'T_336': (2, 2, 5),   # chromatic y·δ→y
            'T_346': (2, 3, 5),   # chromatic y'·δ→y
        }.items():
            ab[name] = s.get(key, 0.0)

    if 3 in all_orders:
        s = all_orders[3]
        for name, key in {
            'U_1666': (0, 5, 5, 5),   # x from δ³
            'U_1111': (0, 0, 0, 0),   # geometric x⁴→x
            'U_1116': (0, 0, 0, 5),   # x²·δ→x
            'U_5666': (4, 5, 5, 5),   # l from δ³
        }.items():
            ab[name] = s.get(key, 0.0)

    if 5 in all_orders:
        s = all_orders[5]
        for name, key in {
            'V_166666': (0, 5, 5, 5, 5, 5),
        }.items():
            ab[name] = s.get(key, 0.0)

    return ab


def print_twiss_table(results, targets):
    print("\n" + "=" * 80)
    print("TWISS PARAMETERS AT UNDULATOR ENTRANCE (from linear map propagation)")
    print("=" * 80)

    hdr = f"{'Order':>5}  {'β_x (m)':>10}  {'α_x':>10}  {'β_y (m)':>10}  {'α_y':>10}  {'η_x (m)':>10}  {'RMS':>10}  {'Time':>6}"
    print(hdr)
    print("-" * len(hdr))

    for r in results:
        t = r['twiss']
        print(f"{r['order']:>5}  {t.get('beta_x',0):>10.6f}  {t.get('alpha_x',0):>10.6f}  "
              f"{t.get('beta_y',0):>10.6f}  {t.get('alpha_y',0):>10.6f}  "
              f"{t.get('eta_x',0):>10.6f}  {r['rms']:>10.4e}  {r['wall_s']:>5.1f}s")

    print(f"{'Target':>5}  {targets['beta_xm']:>10.6f}  {targets['alpha_xm']:>10.6f}  "
          f"{targets['beta_ym']:>10.6f}  {targets['alpha_ym']:>10.6f}")

    # Convergence vs order 1
    if len(results) > 1:
        print("\n  Convergence (vs Order 1):")
        r1 = results[0]
        for r in results[1:]:
            d = {k: abs(r['twiss'].get(k, 0) - r1['twiss'].get(k, 0))
                 for k in ['beta_x', 'alpha_x', 'beta_y', 'alpha_y']}
            print(f"    Order {r['order']}: Δβ_x={d['beta_x']:.2e}, Δα_x={d['alpha_x']:.2e}, "
                  f"Δβ_y={d['beta_y']:.2e}, Δα_y={d['alpha_y']:.2e}")


def print_linear_map_comparison(results):
    print("\n" + "=" * 80)
    print("LINEAR TRANSFER MAP COMPARISON")
    print("=" * 80)

    M_ref = results[0]['linear_map']
    for r in results[1:]:
        M = r['linear_map']
        diff = np.abs(M - M_ref)
        max_diff = np.max(diff)
        rel = np.where(np.abs(M_ref) > 1e-10, diff / np.abs(M_ref), 0)

        print(f"\n  Order {r['order']} vs Order {results[0]['order']}:")
        print(f"    Max |ΔM|    = {max_diff:.6e}")
        print(f"    Max rel |Δ| = {np.max(rel):.6e}")

        if max_diff > 1e-10:
            flat = diff.flatten()
            top3 = np.argsort(flat)[-3:][::-1]
            for idx in top3:
                i, j = divmod(idx, 6)
                if diff[i, j] > 1e-10:
                    print(f"    M[{i+1},{j+1}] = {M_ref[i,j]:.6e} → {M[i,j]:.6e}  (Δ={diff[i,j]:.2e})")


def print_aberration_table(results):
    print("\n" + "=" * 80)
    print("KEY ABERRATION COEFFICIENTS")
    print("=" * 80)

    all_names = sorted(set().union(*(r['aberrations'].keys() for r in results)))

    hdr = f"{'Coeff':>10}" + "".join(f"{'O'+str(r['order']):>14}" for r in results)
    print(hdr)
    print("-" * len(hdr))

    for name in all_names:
        row = f"{name:>10}"
        for r in results:
            v = r['aberrations'].get(name)
            if v is None:
                row += f"{'—':>14}"
            elif abs(v) < 1e-15:
                row += f"{'0':>14}"
            else:
                row += f"{v:>14.4e}"
        print(row)


def print_particle_stats(results):
    print("\n" + "=" * 80)
    print("FINAL BEAM STATISTICS")
    print("=" * 80)

    hdr = f"{'Order':>5}  {'N_surv':>6}  {'σ_x (mm)':>10}  {'σ_xp (mr)':>10}  {'σ_y (mm)':>10}  {'σ_yp (mr)':>10}"
    print(hdr)
    print("-" * len(hdr))

    for r in results:
        def fmt(v): return f"{v:.4f}" if v is not None else "—"
        print(f"{r['order']:>5}  {r['n_survived']:>6}  "
              f"{fmt(r['sigma_x']):>10}  {fmt(r['sigma_xp']):>10}  "
              f"{fmt(r['sigma_y']):>10}  {fmt(r['sigma_yp']):>10}")


def print_map_complexity(results):
    print("\n" + "=" * 80)
    print("TRANSFER MAP COMPLEXITY (nonzero coefficients)")
    print("=" * 80)
    for r in results:
        parts = [f"O{o}: {n}" for o, n in sorted(r['n_coeffs'].items())]
        print(f"  Order {r['order']:>1}: {', '.join(parts)}  (total: {sum(r['n_coeffs'].values())})")


def run_study(current_file, fringe_order, label, n_particles=N_PARTICLES):
    orders = [1, 2, 3, 5]
    targets = compute_targets()

    print(f"\n{'#' * 80}")
    print(f"# P8: Order-by-Order Convergence — {label}")
    print(f"# FR = {fringe_order}, N = {n_particles}")
    print(f"{'#' * 80}")

    currents = load_currents(
        Path(__file__).resolve().parent / 'results' / current_file)

    results = []
    for order in orders:
        print(f"\n  Order {order} (FR {fringe_order})...", end=" ", flush=True)
        try:
            r = run_at_order(order, fringe_order, currents, targets,
                             f'{label}/fr{fringe_order}', n_particles)
            results.append(r)
            print(f"done ({r['wall_s']:.1f}s, {r['n_survived']}/{n_particles} surv, "
                  f"RMS={r['rms']:.4e})")
        except Exception as e:
            print(f"FAILED: {e}")
            import traceback
            traceback.print_exc()

    if len(results) < 2:
        print("  Insufficient successful runs.")
        return results

    print_twiss_table(results, targets)
    print_linear_map_comparison(results)
    print_aberration_table(results)
    print_particle_stats(results)
    print_map_complexity(results)

    # Save JSON summary
    summary = {
        'label': label, 'fringe_order': fringe_order,
        'current_file': current_file, 'n_particles': n_particles,
        'targets': targets,
        'orders': [{
            'order': r['order'], 'wall_s': r['wall_s'],
            'twiss': r['twiss'], 'cosy_twiss': r['cosy_twiss'],
            'mse': r['mse'], 'rms': r['rms'],
            'n_survived': r['n_survived'],
            'sigma_x': r['sigma_x'], 'sigma_y': r['sigma_y'],
            'n_coeffs': {str(k): v for k, v in r['n_coeffs'].items()},
            'aberrations': r['aberrations'],
        } for r in results],
    }

    out_json = RESULTS_DIR / label / f'fr{fringe_order}' / 'summary.json'
    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Saved: {out_json}")

    return results


def main():
    import argparse
    parser = argparse.ArgumentParser(description='P8: DA order convergence')
    parser.add_argument('--felsim-only', '--fr0-only', action='store_true',
                        help='Only run FR0 currents')
    parser.add_argument('--cosy-only', '--fr3-only', action='store_true',
                        help='Only run FR3 currents')
    parser.add_argument('--particles', type=int, default=500)
    parser.add_argument('--fr', type=int, default=0)
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Use COSY-optimized currents (FELsim currents are unstable in COSY
    # due to different dipole edge models — see W4/R2)
    if not args.cosy_only:
        run_study('cosy_s1_fr0.json', args.fr,
                  'cosy_fr0_currents', args.particles)

    if not args.felsim_only:
        run_study('cosy_s1_fr3_warm.json', 3,
                  'cosy_fr3_currents', args.particles)


if __name__ == '__main__':
    main()
