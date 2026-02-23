"""S9 Part B4: RF-Track and COSY INFINITY validation with collective effects.

Compares FELsim linear results with:
  - RF-Track particle tracking (with optional space charge / CSR)
  - COSY INFINITY transfer maps (with fringe fields and longitudinal tracking)

Author: Eremey Valetov
"""

import sys
import json
import copy
import argparse
from pathlib import Path
import numpy as np

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ebeam import beam
from beamline import lattice

# ── Beam parameters ──────────────────────────────────────────────────────────
Energy = 40       # MeV (design energy)
# FELsim lattice defaults to E=45 in lattice.__init__, but the study scripts
# call changeBeamType("electron", 40, elements) to set all elements to E=40.
# RF-Track must use the same energy.
F_RF = 2856e6     # Hz
epsilon_n = 8     # pi.mm.mrad (normalized)
x_std = 0.8       # mm
C_LIGHT = 299792458.0

OUTDIR = Path(__file__).resolve().parent / 'results' / 'S9'
EXCEL_PATH = Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx'


def compute_targets():
    relat = lattice(1, fringeType=None)
    relat.setE(E=Energy)
    norm = relat.gamma * relat.beta
    epsilon = epsilon_n / norm

    K = 1.2
    lambda_u = 2.3e-2
    beta_ym = relat.gamma * lambda_u / (2 * np.pi * K)
    beta_xm = 1.4
    alpha_xm = 0.47
    alpha_ym = 0.0
    beta_0 = x_std**2 / epsilon

    return {
        'beta_xm': beta_xm, 'alpha_xm': alpha_xm,
        'beta_ym': beta_ym, 'alpha_ym': alpha_ym,
        'epsilon': epsilon, 'beta_0': beta_0,
        'gamma': relat.gamma, 'beta_rel': relat.beta,
    }


def create_beam(bunch_spread_ps, energy_std_pct, h_chirp, nb_particles=1000, seed=42):
    np.random.seed(seed)
    relat = lattice(1, fringeType=None)
    relat.setE(E=Energy)
    norm = relat.gamma * relat.beta
    epsilon = epsilon_n / norm
    x_prime_std = epsilon / x_std
    y_prime_std = epsilon / x_std
    tof_std = bunch_spread_ps * 1e-9 * F_RF
    energy_std = energy_std_pct * 10

    ebeam_obj = beam()
    beam_dist = ebeam_obj.gen_6d_gaussian(
        0, [x_std, x_prime_std, x_std, y_prime_std, tof_std, energy_std],
        nb_particles
    )
    tof_dist = beam_dist[:, 4] / F_RF
    beam_dist[:, 5] += h_chirp * tof_dist
    return beam_dist


# ═══════════════════════════════════════════════════════════════════════════════
#  COSY INFINITY Validation
# ═══════════════════════════════════════════════════════════════════════════════

def run_cosy_r56(file_path, targets, fringe_order=0, order=3):
    """Run COSY to extract full 6×6 transfer map and R56.

    Uses COSY transfer-matrix mode WITHOUT optimization (no FIT blocks)
    to get the linear transport matrix and compare R56 with FELsim.
    """
    from cosyAdapter import COSYAdapter
    from cosyOptHelper import parse_beamline_felsim_indexed

    adapter = COSYAdapter(
        lattice_path=str(file_path), mode='transfer_matrix',
        config={'simulation': {'KE': Energy, 'order': order, 'dimensions': 3}},
        fringe_field_order=fringe_order, debug=False
    )
    sim = adapter.get_native_simulator()
    sim.beamline = parse_beamline_felsim_indexed(str(file_path))[:118]

    sim.set_geometric_emittance(targets['epsilon'])
    beta_0 = targets['beta_0']
    sim.set_initial_twiss(beta_x=beta_0, alpha_x=0.0, beta_y=beta_0, alpha_y=0.0)

    # No optimization — just propagate
    sim.set_optimization_enabled(False)

    print(f"Running COSY (FR {fringe_order}, order {order}) for transfer map...")
    sim_result = sim.run_simulation()

    if sim_result.get('status') != 'success':
        print(f"  COSY FAILED: {sim_result}")
        return None

    reader = sim.analyze_results()
    M = reader.read_linear_transfer_map()

    twiss = reader.get_twiss_from_transfer_map(
        initial_twiss_x={'beta': beta_0, 'alpha': 0.0},
        initial_twiss_y={'beta': beta_0, 'alpha': 0.0},
    )

    R56_cosy = M[4, 5] if M is not None else None

    print(f"  COSY R56 (M[4,5]) = {R56_cosy:.6e}" if R56_cosy is not None else "  M unavailable")
    print(f"  COSY β_x = {twiss.get('beta_x', '?'):.4f} m" if twiss else "  Twiss unavailable")

    return {'M': M, 'twiss': twiss, 'R56': R56_cosy}


def run_cosy_optimization_s9(file_path, targets, bunch_spread_ps, sigma_e_pct, h_chirp,
                              fringe_order=0, order=3, nmax=1000, nalg=1):
    """Run COSY 11-stage optimization for a given beam configuration.

    Reuses the build_stages() from UHM_beamline_opt_cosy.py.
    """
    # Import the COSY optimization machinery
    from cosyAdapter import COSYAdapter
    from cosyOptHelper import add_stages, get_optimized_currents, parse_beamline_felsim_indexed

    # Build stages (same as W4)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from UHM_beamline_opt_cosy import build_stages, compute_mse

    stages = build_stages(targets)

    adapter = COSYAdapter(
        lattice_path=str(file_path), mode='transfer_matrix',
        config={'simulation': {'KE': Energy, 'order': order, 'dimensions': 3}},
        fringe_field_order=fringe_order, debug=False
    )
    sim = adapter.get_native_simulator()
    sim.beamline = parse_beamline_felsim_indexed(str(file_path))[:118]

    sim.set_geometric_emittance(targets['epsilon'])
    beta_0 = targets['beta_0']
    sim.set_initial_twiss(beta_x=beta_0, alpha_x=0.0, beta_y=beta_0, alpha_y=0.0)

    sim.fit_nmax = nmax
    sim.fit_eps = 1e-8
    sim.fit_nalgorithm = nalg
    sim.fit_combined_mse = True

    result = add_stages(sim, stages)

    print(f"\nRunning COSY FIT: σ_t={bunch_spread_ps} ps, σ_E={sigma_e_pct}%, "
          f"h={h_chirp:.0e}, FR {fringe_order}")
    sim_result = sim.run_simulation()

    if sim_result.get('status') != 'success':
        print(f"  COSY FAILED")
        return {'success': False}

    reader = sim.analyze_results()
    twiss = reader.get_twiss_from_transfer_map(
        initial_twiss_x={'beta': beta_0, 'alpha': 0.0},
        initial_twiss_y={'beta': beta_0, 'alpha': 0.0},
    )
    currents = get_optimized_currents(reader, stages)
    mse = compute_mse(twiss, targets)

    print(f"  MSE = {mse:.6e}")
    print(f"  β_x = {twiss.get('beta_x', '?'):.4f}, α_x = {twiss.get('alpha_x', '?'):.4f}")
    print(f"  β_y = {twiss.get('beta_y', '?'):.4f}, α_y = {twiss.get('alpha_y', '?'):.4f}")

    return {
        'success': True, 'twiss': twiss, 'currents': currents, 'mse': mse,
        'bunch_ps': bunch_spread_ps, 'sigma_e': sigma_e_pct, 'h': h_chirp,
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  RF-Track Validation
# ═══════════════════════════════════════════════════════════════════════════════

def check_rftrack():
    """Check if RF-Track is available."""
    try:
        import RF_Track as rft
        return True
    except ImportError:
        return False


def get_optimized_felsim_currents(beam_dist):
    """Run FELsim 11-stage optimization and return {element_index: current} dict."""
    from excelElements import ExcelElements
    from beamOptimizer import beamOptimizer
    from beamline import lattice as lattice_cls

    excel = ExcelElements(str(EXCEL_PATH))
    beamline_elements = excel.create_beamline()

    # Apply E=40 MeV to all elements (matching load_beamline() in the study script)
    relat = lattice_cls(1, fringeType=None)
    relat.setE(E=Energy)
    line = relat.changeBeamType("electron", Energy, beamline_elements)

    segments = 118
    line_opt = line[:segments]
    opti = beamOptimizer(line_opt, beam_dist)

    norm = relat.gamma * relat.beta
    epsilon = epsilon_n / norm
    K = 1.2
    lambda_u = 2.3e-2
    beta_ym = relat.gamma * lambda_u / (2 * np.pi * K)
    alpha_ym = 0.0
    beta_xm = 1.4
    alpha_xm = 0.47

    stages_config = [
        ({1: ["I", "current", lambda n: n], 3: ["I2", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}, "I2": {"bounds": (0, 10), "start": 1}},
         {8: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
              {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.0}],
          9: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
              {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]}),
        ({10: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {15: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}),
        ({16: ["I", "current", lambda n: n], 18: ["I2", "current", lambda n: n],
          20: ["I3", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 5},
          "I3": {"bounds": (0, 10), "start": 3}},
         {25: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
          26: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]}),
        ({27: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {32: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}),
        ({37: ["I", "current", lambda n: n], 35: ["I2", "current", lambda n: n],
          33: ["I3", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2},
          "I3": {"bounds": (0, 10), "start": 2}},
         {37: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["x", "envelope"], "goal": 2.0, "weight": 1},
               {"measure": ["y", "envelope"], "goal": 2.0, "weight": 1}]}),
        ({50: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {55: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}),
        ({56: ["I", "current", lambda n: n], 58: ["I2", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2}},
         {59: [{"measure": ["x", "envelope"], "goal": 0.0, "weight": 1},
               {"measure": ["y", "envelope"], "goal": 0.0, "weight": 1}]}),
        ({61: ["I", "current", lambda n: n], 63: ["I2", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2}},
         {68: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
          69: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]}),
        ({70: ["I", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 1}},
         {75: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}),
        ({76: ["I", "current", lambda n: n], 78: ["I2", "current", lambda n: n],
          80: ["I3", "current", lambda n: n]},
         {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2},
          "I3": {"bounds": (0, 10), "start": 2}},
         {85: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
          86: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
               {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]}),
        ({87: ["Ic", "current", lambda n: n], 93: ["I", "current", lambda n: n],
          95: ["I2", "current", lambda n: n], 97: ["I3", "current", lambda n: n]},
         {"Ic": {"bounds": (0, 10), "start": 4}, "I": {"bounds": (0, 10), "start": 2},
          "I2": {"bounds": (0, 10), "start": 2}, "I3": {"bounds": (0, 10), "start": 2}},
         {92: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 0.5}],
          117: [{"measure": ["x", "alpha"], "goal": alpha_xm, "weight": 1},
                {"measure": ["y", "alpha"], "goal": alpha_ym, "weight": 1},
                {"measure": ["x", "beta"], "goal": beta_xm, "weight": 1},
                {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]}),
    ]

    for i, (variables, startPoint, objectives) in enumerate(stages_config):
        opti.calc("Nelder-Mead", variables, startPoint, objectives,
                  plotBeam=False, printResults=False, plotProgress=False)
        if i == 4:
            line_opt[43].current = line_opt[33].current
            line_opt[41].current = line_opt[35].current
            line_opt[39].current = line_opt[37].current

    quad_indices = [1, 3, 10, 16, 18, 20, 27, 33, 35, 37, 39, 41, 43,
                    50, 56, 58, 61, 63, 70, 76, 78, 80, 87, 93, 95, 97]
    currents = {i: line_opt[i].current for i in quad_indices if hasattr(line_opt[i], 'current')}
    return currents


def run_rftrack_validation(beam_dist, label, optimized_currents=None,
                           space_charge=False, sc_mesh=(32, 32, 64),
                           aperture=0.5):
    """Run RF-Track particle tracking with optimized currents.

    If optimized_currents is provided, applies them to the RF-Track
    adapter before tracking. Otherwise uses default (unoptimized) lattice.
    Aperture defaults to 0.5 m (large) to avoid artificial particle loss.
    """
    try:
        from rftrackAdapter import RFTrackAdapter
    except ImportError:
        print("  RF-Track adapter not available")
        return None

    try:
        sim = RFTrackAdapter(
            lattice_path=str(EXCEL_PATH),
            beam_energy=Energy,
            space_charge=space_charge,
            sc_mesh=sc_mesh,
            aperture=aperture,
        )
    except ImportError:
        print("  RF-Track package not installed (pip install RF-Track)")
        return None

    # Truncate to 118 elements (matching FELsim optimizer)
    sim.beamline = sim.beamline[:118]
    sim._build_lattice()

    # Apply optimized quad currents
    if optimized_currents:
        for idx, current in optimized_currents.items():
            if idx < len(sim.beamline):
                sim._modify_element(idx, current=current)
        sim._build_lattice()

    print(f"\n── RF-Track: {label} (SC={'ON' if space_charge else 'OFF'}) ──")
    result = sim.simulate(beam_dist)

    if not result.success:
        print("  FAILED")
        return None

    twiss = result.twiss_parameters_statistical.get('final', {})
    n_good = result.metadata.get('num_good', 0)
    n_lost = result.metadata.get('num_lost', 0)

    print(f"  Particles: {n_good} good, {n_lost} lost")
    for plane in ['x', 'y']:
        if plane in twiss:
            t = twiss[plane]
            print(f"  {plane}: β={t.get('beta',0):.4f} m, α={t.get('alpha',0):.4f}, "
                  f"ε={t.get('emittance',0):.4f} π·mm·mrad")

    return {'twiss': twiss, 'metadata': result.metadata, 'particles': result.final_particles}


# ═══════════════════════════════════════════════════════════════════════════════
#  Comparison
# ═══════════════════════════════════════════════════════════════════════════════

def compare_results(felsim_twiss, cosy_results, rftrack_results, targets, label):
    """Print comparison table of Twiss at undulator entrance."""
    print(f"\n{'='*72}")
    print(f"  Comparison: {label}")
    print(f"{'='*72}")

    print(f"\n{'Parameter':<20s}  {'Target':>8s}  {'FELsim':>8s}  {'COSY':>8s}  {'RF-Track':>10s}")
    print("-" * 60)

    params = [
        ('β_x (m)', 'beta_x', targets['beta_xm']),
        ('α_x', 'alpha_x', targets['alpha_xm']),
        ('β_y (m)', 'beta_y', targets['beta_ym']),
        ('α_y', 'alpha_y', targets['alpha_ym']),
    ]

    for name, key, target in params:
        fs_val = felsim_twiss.get(key, '—') if felsim_twiss else '—'
        cosy_val = cosy_results.get('twiss', {}).get(key, '—') if cosy_results and cosy_results.get('success') else '—'
        rft_plane = key.split('_')[1]
        rft_param = key.split('_')[0]
        rft_val = '—'
        if rftrack_results and rftrack_results.get('twiss'):
            rft_twiss = rftrack_results['twiss'].get(rft_plane, {})
            rft_val = rft_twiss.get(rft_param, '—')

        fs_str = f"{fs_val:8.4f}" if isinstance(fs_val, (int, float)) else f"{fs_val:>8s}"
        cosy_str = f"{cosy_val:8.4f}" if isinstance(cosy_val, (int, float)) else f"{cosy_val:>8s}"
        rft_str = f"{rft_val:10.4f}" if isinstance(rft_val, (int, float)) else f"{rft_val:>10s}"

        print(f"{name:<20s}  {target:8.4f}  {fs_str}  {cosy_str}  {rft_str}")


# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="S9 B4: RF-Track + COSY validation")
    parser.add_argument('--cosy-r56', action='store_true',
                        help='Extract R56 from COSY transfer map')
    parser.add_argument('--cosy-opt', action='store_true',
                        help='Run COSY optimization for S1/S3 and compare')
    parser.add_argument('--rftrack', action='store_true',
                        help='Run RF-Track validation (requires RF-Track)')
    parser.add_argument('--rftrack-sc', action='store_true',
                        help='Run RF-Track with space charge')
    parser.add_argument('--fr', type=int, default=0,
                        help='COSY fringe field order (0, 1, 3)')
    parser.add_argument('--order', type=int, default=3,
                        help='COSY DA order')
    parser.add_argument('--all', action='store_true', help='Run all available')
    args = parser.parse_args()

    if not any([args.cosy_r56, args.cosy_opt, args.rftrack, args.rftrack_sc, args.all]):
        args.all = True

    OUTDIR.mkdir(parents=True, exist_ok=True)
    targets = compute_targets()

    print("S9 Part B4: Cross-Validation")
    print(f"E = {Energy} MeV, ε_n = {epsilon_n} μm")
    print(f"Targets: β_x = {targets['beta_xm']:.2f} m, α_x = {targets['alpha_xm']:.2f}, "
          f"β_y = {targets['beta_ym']:.4f} m, α_y = {targets['alpha_ym']:.2f}")

    # ── COSY R56 ──
    if args.cosy_r56 or args.all:
        print("\n" + "="*72)
        print("  COSY: R56 Extraction")
        print("="*72)
        for fr in [0, 1]:
            r56_result = run_cosy_r56(EXCEL_PATH, targets, fringe_order=fr, order=args.order)
            if r56_result and r56_result['M'] is not None:
                M = r56_result['M']
                print(f"\n  FR {fr} — Full 6×6 transfer map (selected elements):")
                print(f"    M[0,0]={M[0,0]:.6f}  M[0,1]={M[0,1]:.6f}  M[0,5]={M[0,5]:.6f}")
                print(f"    M[4,0]={M[4,0]:.6f}  M[4,1]={M[4,1]:.6f}  M[4,5]={M[4,5]:.6f}")

    # ── COSY Optimization ──
    cosy_results = {}
    if args.cosy_opt or args.all:
        print("\n" + "="*72)
        print("  COSY: 11-Stage Optimization Comparison")
        print("="*72)

        scenarios = [
            ("S1 (2 ps, 0.5%, h=5e9)", 2, 0.5, 5e9),
            ("S3 (0.5 ps, 0.5%, h=5e9)", 0.5, 0.5, 5e9),
            ("Pre-comp (0.5 ps, 2%, h=0)", 0.5, 2.0, 0),
        ]

        for label, bunch_ps, sigma_e, h in scenarios:
            result = run_cosy_optimization_s9(
                EXCEL_PATH, targets, bunch_ps, sigma_e, h,
                fringe_order=args.fr, order=args.order
            )
            cosy_results[label] = result

        # Check if S1 and S3 produce identical currents
        s1 = cosy_results.get("S1 (2 ps, 0.5%, h=5e9)")
        s3 = cosy_results.get("S3 (0.5 ps, 0.5%, h=5e9)")
        if s1 and s1.get('success') and s3 and s3.get('success'):
            print("\n── S1 vs S3 Current Comparison (COSY) ──")
            for k in sorted(s1['currents'].keys()):
                d = abs(s1['currents'][k] - s3['currents'].get(k, 0))
                print(f"  I[{k}]: S1={s1['currents'][k]:.4f}  S3={s3['currents'][k]:.4f}  "
                      f"Δ={d:.6f}")

    # ── RF-Track ──
    rftrack_results = {}
    if args.rftrack or args.rftrack_sc or args.all:
        if not check_rftrack():
            print("\nRF-Track not available. Install with: pip install RF-Track")
        else:
            scenarios = [
                ("2 ps baseline", 2, 0.5, 5e9),
                ("0.5 ps baseline", 0.5, 0.5, 5e9),
            ]

            # Run FELsim optimization once to get optimized currents
            print("\n" + "="*72)
            print("  RF-Track: Running FELsim optimization for quad currents...")
            print("="*72)
            ref_beam = create_beam(2, 0.5, 5e9)
            opt_currents = get_optimized_felsim_currents(ref_beam)
            print(f"  Got {len(opt_currents)} optimized quad currents")
            for idx in sorted(opt_currents)[:5]:
                print(f"    I[{idx}] = {opt_currents[idx]:.4f} A")
            print(f"    ...")

            for label, bunch_ps, sigma_e, h in scenarios:
                bd = create_beam(bunch_ps, sigma_e, h)

                # Without space charge
                if args.rftrack or args.all:
                    rft_result = run_rftrack_validation(
                        bd, f"{label} (no SC)", optimized_currents=opt_currents)
                    rftrack_results[f"{label}_noSC"] = rft_result

                # With space charge
                if args.rftrack_sc or args.all:
                    rft_result_sc = run_rftrack_validation(
                        bd, f"{label} (SC on)", optimized_currents=opt_currents,
                        space_charge=True)
                    rftrack_results[f"{label}_SC"] = rft_result_sc

    # ── Summary comparison ──
    print("\n" + "="*72)
    print("  Summary")
    print("="*72)

    # Save results
    summary = {
        'targets': {k: float(v) for k, v in targets.items()},
        'cosy': {},
        'rftrack': {},
    }
    for k, v in cosy_results.items():
        if v and v.get('success'):
            summary['cosy'][k] = {
                'mse': float(v['mse']),
                'twiss': {kk: float(vv) for kk, vv in v['twiss'].items()
                          if isinstance(vv, (int, float))},
            }
    for k, v in rftrack_results.items():
        if v and v.get('twiss'):
            summary['rftrack'][k] = {
                'twiss': {plane: {pp: float(pv) for pp, pv in pd.items()}
                          for plane, pd in v['twiss'].items()},
            }

    json_path = OUTDIR / 'S9_validation_results.json'
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nResults saved to {json_path}")


if __name__ == "__main__":
    main()
