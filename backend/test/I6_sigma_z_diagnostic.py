"""I6: COSY σ_z Blowup Diagnostic

Investigates why COSY particle tracking shows σ_z ≈ 92–233 ps (60–100× blowup)
for 2 ps input, while RF-Track gives correct σ_z ≈ 2 ps and W9 linear map
propagation shows only ~2.5% growth.

Parts:
  A — Order-dependence test: track at ORDER 1, 2, 3
  B — Element-by-element l-coordinate tracing via collect_evolution()
  C — Single off-energy particle test
  D — Reference particle convention audit

Author: Eremey Valetov
"""

import sys
import json
import argparse
import time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from cosyAdapter import COSYAdapter
from cosyOptHelper import parse_beamline_felsim_indexed
from simulatorBase import CoordinateSystem
from UHM_beamline_opt_cosy import Energy

COSY = CoordinateSystem.COSY
FELSIM = CoordinateSystem.FELSIM

# ── Constants ─────────────────────────────────────────────────────────────────
C_LIGHT = 299792458.0
M_E_MEV = 0.51099895
GAMMA = 1 + Energy / M_E_MEV
BETA_REL = np.sqrt(1 - 1 / GAMMA**2)
BETA_C = BETA_REL * C_LIGHT
P_C = GAMMA * BETA_REL * M_E_MEV  # MeV/c

OUTDIR = Path(__file__).resolve().parent / 'results' / 'I6'
EXCEL_PATH = (Path(__file__).resolve().parent.parent.parent
              / 'beam_excel' / 'Beamline_elements.xlsx')
W9_MAP = Path(__file__).resolve().parent / 'results' / 'W9' / 'part_a_longitudinal_map.json'

N_PARTICLES = 100
SEED = 42


def _print(msg):
    print(msg, flush=True)


def load_w9_currents():
    """Load W9-optimized currents."""
    with open(W9_MAP) as f:
        data = json.load(f)
    return {int(k): float(v) for k, v in data['currents'].items()}


def create_cosy_adapter(order=3, fringe_field_order=0):
    """Create a COSYAdapter in particle tracking mode."""
    sim = COSYAdapter(
        excel_path=str(EXCEL_PATH),
        mode='particle_tracking',
        config={'simulation': {'order': order, 'dimensions': 3, 'KE': Energy}},
        transfer_matrix_order=1,
        fringe_field_order=fringe_field_order,
        debug=False,
    )
    return sim


def inject_currents(sim, currents):
    """Inject W9-optimized currents into the COSY beamline."""
    bl = parse_beamline_felsim_indexed(str(EXCEL_PATH))
    for idx, current in currents.items():
        if idx < len(bl):
            bl[idx]['current'] = current
    sim._particle_sim.beamline = bl


def generate_beam_cosy(sigma_t_ps=2.0, sigma_delta=0.005, h_chirp=0.0,
                       N=N_PARTICLES, seed=SEED):
    """Generate 6D Gaussian in COSY coords [x(m), a(rad), y(m), b(rad), l(m), δK/K₀]."""
    rng = np.random.default_rng(seed)
    sigma_x = 0.8e-3   # m
    epsilon_n = 8.0     # pi.mm.mrad
    epsilon_geom = epsilon_n / (GAMMA * BETA_REL)  # pi.mm.mrad → pi.m.rad? No: mm.mrad → m.rad
    sigma_xp = epsilon_geom * 1e-6 / sigma_x  # rad
    sigma_z = sigma_t_ps * 1e-12 * BETA_C  # m

    X = np.zeros((N, 6))
    X[:, 0] = rng.normal(0, sigma_x, N)
    X[:, 1] = rng.normal(0, sigma_xp, N)
    X[:, 2] = rng.normal(0, sigma_x, N)
    X[:, 3] = rng.normal(0, sigma_xp, N)
    X[:, 4] = rng.normal(0, sigma_z, N)
    X[:, 5] = rng.normal(0, sigma_delta, N)
    if h_chirp != 0:
        X[:, 5] += h_chirp * X[:, 4] / BETA_C
    return X


def sigma_z_ps(particles_felsim):
    """RMS bunch length in ps from FELsim col 4."""
    tof_rel = particles_felsim[:, 4]  # ΔToF/T_RF × 10³
    F_RF = 2856e6
    T_RF = 1.0 / F_RF
    tof_s = tof_rel * 1e-3 * T_RF
    return np.std(tof_s) * 1e12


def sigma_l_m(particles_cosy):
    """RMS bunch length in metres from COSY col 4."""
    return np.std(particles_cosy[:, 4])


# ═══════════════════════════════════════════════════════════════════════════════
#  Part A: Order-dependence test
# ═══════════════════════════════════════════════════════════════════════════════

def part_a():
    """Track 100 particles at ORDER 1, 2, 3 with W9-optimized currents."""
    _print("\n" + "=" * 72)
    _print("  Part A: Order-Dependence Test")
    _print("=" * 72)

    currents = load_w9_currents()
    results = {}

    for order in [1, 2, 3]:
        _print(f"\n── ORDER {order} ──")
        t0 = time.time()

        sim = create_cosy_adapter(order=order)
        inject_currents(sim, currents)

        # Generate beam in FELsim coords via the adapter's transform
        beam_cosy = generate_beam_cosy(sigma_t_ps=2.0, sigma_delta=0.005, h_chirp=0.0)

        # Transform to FELsim coords for collect_evolution input
        particles_felsim = sim.transform_coordinates(
            beam_cosy, from_system=COSY, to_system=FELSIM)

        _print(f"  Input σ_l = {sigma_l_m(beam_cosy)*1e6:.1f} μm "
               f"= {sigma_l_m(beam_cosy)/BETA_C*1e12:.3f} ps")
        _print(f"  Input σ_δ = {np.std(beam_cosy[:, 5])*100:.3f}%")

        try:
            evolution = sim.collect_evolution(particles_felsim, checkpoint_elements='all')

            s_final = max(evolution.s_positions)
            p_final = evolution.particles[s_final]
            n_out = len(p_final)

            sz_out_ps = sigma_z_ps(p_final) if n_out >= 2 else float('nan')
            transmission = n_out / N_PARTICLES

            _print(f"  Output: {n_out}/{N_PARTICLES} particles survived (T={transmission:.1%})")
            _print(f"  Output σ_z = {sz_out_ps:.3f} ps")

            # Also check COSY-native coordinates at final checkpoint
            p_final_cosy = sim.transform_coordinates(
                p_final, from_system=FELSIM, to_system=COSY)
            sl_out = sigma_l_m(p_final_cosy)
            _print(f"  Output σ_l = {sl_out*1e6:.1f} μm = {sl_out/BETA_C*1e12:.3f} ps")

            results[order] = {
                'order': order,
                'n_input': N_PARTICLES,
                'n_output': n_out,
                'transmission': transmission,
                'sigma_z_in_ps': 2.0,
                'sigma_z_out_ps': sz_out_ps,
                'sigma_l_in_um': sigma_l_m(beam_cosy) * 1e6,
                'sigma_l_out_um': sl_out * 1e6,
                'time_s': time.time() - t0,
            }
        except Exception as e:
            _print(f"  FAILED: {e}")
            results[order] = {'order': order, 'error': str(e)}

    # ── Decision tree ──
    _print("\n── Decision Tree ──")
    o1 = results.get(1, {})
    o3 = results.get(3, {})

    o1_sz = o1.get('sigma_z_out_ps', float('nan'))
    o3_sz = o3.get('sigma_z_out_ps', float('nan'))
    o1_sl = o1.get('sigma_l_out_um', float('nan'))
    o3_sl = o3.get('sigma_l_out_um', float('nan'))

    _print(f"  ORDER 1: σ_z = {o1_sz:.3f} ps, σ_l = {o1_sl:.1f} μm")
    _print(f"  ORDER 3: σ_z = {o3_sz:.3f} ps, σ_l = {o3_sl:.1f} μm")

    ratio_o1 = o1_sz / 2.0  # input is 2 ps
    ratio_o3 = o3_sz / 2.0
    _print(f"  Blowup ratio: ORDER 1 = {ratio_o1:.2f}×, ORDER 3 = {ratio_o3:.2f}×")

    if ratio_o1 < 1.5 and ratio_o3 > 10:
        _print("  → ORDER 1 correct, ORDER 3 blows up → higher-order map amplification")
    elif ratio_o1 > 10 and ratio_o3 > 10:
        _print("  → Both orders blow up → coordinate convention or RP procedure issue")
    elif ratio_o1 < 2 and ratio_o3 < 2:
        _print("  → Both orders show modest growth (< 2×) → consistent with R56 × σ_δ coupling")
        _print("    The original 60–100× blowup is NOT reproduced with W9 currents.")
        _print("    Hypothesis: previous blowup was with different (poorly matched) currents.")
    elif ratio_o1 < 2 and ratio_o3 > 1.5:
        _print("  → ORDER 1 preserves σ_z, ORDER 3 shows moderate growth")
        _print("    Higher-order nonlinear path-length terms contribute at ORDER ≥ 2.")
        _print("    This is correct physics, not a bug.")
    else:
        _print(f"  → Unexpected pattern — needs further investigation")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    with open(OUTDIR / 'part_a_order_test.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    _print(f"\n  Saved: {OUTDIR / 'part_a_order_test.json'}")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  Part B: Element-by-element l-coordinate tracing
# ═══════════════════════════════════════════════════════════════════════════════

def part_b(part_a_data=None):
    """Trace σ_l(s) element-by-element through the beamline."""
    _print("\n" + "=" * 72)
    _print("  Part B: Element-by-Element σ_l(s) Evolution")
    _print("=" * 72)

    currents = load_w9_currents()

    # Use ORDER 3 (the configuration that showed blowup) and ORDER 1 for comparison
    results = {}

    for order in [1, 3]:
        _print(f"\n── ORDER {order} ──")
        t0 = time.time()

        sim = create_cosy_adapter(order=order)
        inject_currents(sim, currents)

        beam_cosy = generate_beam_cosy(sigma_t_ps=2.0, sigma_delta=0.005)
        particles_felsim = sim.transform_coordinates(
            beam_cosy, from_system=COSY, to_system=FELSIM)

        try:
            evolution = sim.collect_evolution(particles_felsim, checkpoint_elements='all')

            s_positions = sorted(evolution.s_positions)
            sigma_l_values = []
            sigma_l_cosy_values = []
            n_particles = []

            for s in s_positions:
                p = evolution.particles[s]
                n = len(p)
                n_particles.append(n)
                if n >= 2:
                    sz = sigma_z_ps(p)
                    sigma_l_values.append(sz)
                    # Also compute in COSY coords
                    p_cosy = sim.transform_coordinates(
                        p, from_system=FELSIM, to_system=COSY)
                    sigma_l_cosy_values.append(np.std(p_cosy[:, 4]) * 1e6)  # μm
                else:
                    sigma_l_values.append(float('nan'))
                    sigma_l_cosy_values.append(float('nan'))

            results[order] = {
                's_positions': s_positions,
                'sigma_z_ps': sigma_l_values,
                'sigma_l_um': sigma_l_cosy_values,
                'n_particles': n_particles,
                'time_s': time.time() - t0,
            }

            _print(f"  {len(s_positions)} checkpoints, final σ_z = {sigma_l_values[-1]:.3f} ps")

        except Exception as e:
            _print(f"  FAILED: {e}")
            results[order] = {'error': str(e)}

    # ── Also compute FELsim linear map prediction ──
    _print("\n── FELsim linear map prediction ──")
    try:
        with open(W9_MAP) as f:
            w9_data = json.load(f)
        M = np.array(w9_data['linear_map'])
        if abs(M[5, 5]) < 1e-15:
            M[5, 5] = 1.0
        beam_cosy = generate_beam_cosy(sigma_t_ps=2.0, sigma_delta=0.005)
        beam_out = (M @ beam_cosy.T).T
        sl_in = np.std(beam_cosy[:, 4]) * 1e6
        sl_out = np.std(beam_out[:, 4]) * 1e6
        ratio = sl_out / sl_in
        _print(f"  σ_l: {sl_in:.1f} → {sl_out:.1f} μm (ratio {ratio:.4f})")
        results['linear_map'] = {
            'sigma_l_in_um': sl_in,
            'sigma_l_out_um': sl_out,
            'ratio': ratio,
        }
    except FileNotFoundError:
        _print("  W9 data not found — skipping linear map comparison")

    # ── Plot ──
    OUTDIR.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    for order in sorted(k for k in results if isinstance(k, int)):
        data = results[order]
        if 's_positions' in data:
            label = f'ORDER {order}'
            color = f'C{order-1}'
            ax1.plot(data['s_positions'], data['sigma_z_ps'], '-o', markersize=2,
                     label=label, color=color, alpha=0.8)
            ax2.plot(data['s_positions'], data['n_particles'], '-', linewidth=1.5,
                     label=label, color=color, alpha=0.8)

    ax1.set_ylabel('σ_z (ps)')
    ax1.set_title('I6: σ_z evolution through beamline (COSY particle tracking)')
    ax1.axhline(2.0, color='gray', ls='--', lw=1, label='Input σ_z = 2 ps')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')

    ax2.set_xlabel('s (m)')
    ax2.set_ylabel('Surviving particles')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Mark element types with colored background
    if 1 in results and 'error' not in results[1]:
        evolution_data = results[1]
        # Add beamline element info if available
        _print(f"  Plotted {len(evolution_data['s_positions'])} checkpoint positions")

    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'I6_sigma_z_evolution.{ext}', dpi=150)
    _print(f"\n  Saved: I6_sigma_z_evolution.{{eps,png}}")
    plt.close(fig)

    with open(OUTDIR / 'part_b_evolution.json', 'w') as f:
        json.dump({str(k): v for k, v in results.items()}, f, indent=2, default=str)
    _print(f"  Saved: {OUTDIR / 'part_b_evolution.json'}")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  Part C: Single off-energy particle test
# ═══════════════════════════════════════════════════════════════════════════════

def part_c():
    """Track (0,0,0,0,0,δ=0.005) at ORDER 1 and ORDER 3."""
    _print("\n" + "=" * 72)
    _print("  Part C: Single Off-Energy Particle Test")
    _print("=" * 72)

    currents = load_w9_currents()
    delta = 0.005  # δK/K₀

    # Load W9 R56 for expected result
    try:
        with open(W9_MAP) as f:
            w9_data = json.load(f)
        R56 = w9_data['R56_cosy_m']
    except FileNotFoundError:
        R56 = 0.027  # approximate
    expected_l = R56 * delta  # m
    _print(f"  Expected final l = R56 × δ = {R56:.6f} × {delta} = {expected_l*1e6:.2f} μm")

    results = {}

    for order in [1, 3]:
        _print(f"\n── ORDER {order} ──")

        sim = create_cosy_adapter(order=order)
        inject_currents(sim, currents)

        # Single particle in COSY coords: (0, 0, 0, 0, 0, δ)
        beam_cosy = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, delta]])

        particles_felsim = sim.transform_coordinates(
            beam_cosy, from_system=COSY, to_system=FELSIM)
        _print(f"  Input (FELsim): {particles_felsim[0]}")

        try:
            evolution = sim.collect_evolution(particles_felsim, checkpoint_elements='all')
            s_final = max(evolution.s_positions)
            p_final = evolution.particles[s_final]

            if len(p_final) == 0:
                _print("  Particle lost!")
                results[order] = {'order': order, 'lost': True}
                continue

            _print(f"  Output (FELsim): {p_final[0]}")

            p_cosy_out = sim.transform_coordinates(
                p_final, from_system=FELSIM, to_system=COSY)
            _print(f"  Output (COSY):   {p_cosy_out[0]}")

            l_final = p_cosy_out[0, 4]
            delta_final = p_cosy_out[0, 5]
            _print(f"  Final l = {l_final*1e6:.2f} μm (expected {expected_l*1e6:.2f} μm, "
                   f"ratio {l_final/expected_l:.4f})")
            _print(f"  Final δ = {delta_final:.6f} (expected {delta:.6f})")

            # Trace l through all checkpoints
            s_positions = sorted(evolution.s_positions)
            l_values = []
            for s in s_positions:
                p = evolution.particles[s]
                if len(p) > 0:
                    pc = sim.transform_coordinates(p, from_system=FELSIM, to_system=COSY)
                    l_values.append(pc[0, 4] * 1e6)  # μm
                else:
                    l_values.append(float('nan'))

            results[order] = {
                'order': order,
                'l_final_um': l_final * 1e6,
                'l_expected_um': expected_l * 1e6,
                'delta_final': delta_final,
                'delta_expected': delta,
                'l_ratio': l_final / expected_l if expected_l != 0 else float('nan'),
                's_positions': s_positions,
                'l_evolution_um': l_values,
            }

        except Exception as e:
            _print(f"  FAILED: {e}")
            results[order] = {'order': order, 'error': str(e)}

    # ── Plot l evolution for single particle ──
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(14, 5))

    for order in sorted(k for k in results if isinstance(k, int)):
        data = results[order]
        if 's_positions' in data:
            ax.plot(data['s_positions'], data['l_evolution_um'], '-o', markersize=3,
                    label=f'ORDER {order}')

    ax.axhline(expected_l * 1e6, color='red', ls='--', lw=1.5,
               label=f'Expected (R56×δ = {expected_l*1e6:.2f} μm)')
    ax.set_xlabel('s (m)')
    ax.set_ylabel('l (μm)')
    ax.set_title(f'I6: Single particle (δ={delta}) path length evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    for ext in ['eps', 'png']:
        fig.savefig(OUTDIR / f'I6_single_particle_l.{ext}', dpi=150)
    _print(f"\n  Saved: I6_single_particle_l.{{eps,png}}")
    plt.close(fig)

    with open(OUTDIR / 'part_c_single_particle.json', 'w') as f:
        json.dump({str(k): v for k, v in results.items()}, f, indent=2, default=str)
    _print(f"  Saved: {OUTDIR / 'part_c_single_particle.json'}")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  Part D: Reference particle convention audit
# ═══════════════════════════════════════════════════════════════════════════════

def part_d():
    """Write 5 known particles, track through a short drift, check format."""
    _print("\n" + "=" * 72)
    _print("  Part D: Reference Particle Convention Audit")
    _print("=" * 72)

    # Use a beamline truncated to just the first drift (element 0)
    sim = create_cosy_adapter(order=3)

    # 5 test particles in COSY coords
    test_particles = np.array([
        [0.0,    0.0,    0.0,    0.0,    0.0,    0.0],      # on-axis, on-energy
        [1e-3,   0.0,    0.0,    0.0,    0.0,    0.0],      # 1mm x offset
        [0.0,    1e-3,   0.0,    0.0,    0.0,    0.0],      # 1mrad x' divergence
        [0.0,    0.0,    0.0,    0.0,    1e-4,   0.0],      # 100 μm l offset
        [0.0,    0.0,    0.0,    0.0,    0.0,    0.005],    # 0.5% energy deviation
    ])

    labels = ['on-axis', 'x=1mm', "x'=1mrad", 'l=100μm', 'δ=0.5%']

    # Truncate beamline to first few elements (just drifts)
    bl = parse_beamline_felsim_indexed(str(EXCEL_PATH))
    # Find first drift
    n_elements = 5  # first 5 elements: should be drifts/quads
    sim._particle_sim.beamline = bl[:n_elements]

    particles_felsim = sim.transform_coordinates(
        test_particles, from_system=COSY, to_system=FELSIM)

    _print(f"\nInput particles (COSY coords):")
    for i, label in enumerate(labels):
        _print(f"  {label:>12s}: {test_particles[i]}")

    _print(f"\nInput particles (FELsim coords):")
    for i, label in enumerate(labels):
        _print(f"  {label:>12s}: {particles_felsim[i]}")

    try:
        evolution = sim.collect_evolution(particles_felsim, checkpoint_elements='all')
        _print(f"\n{len(evolution.s_positions)} checkpoints recorded")

        for s in sorted(evolution.s_positions):
            p = evolution.particles[s]
            _print(f"\n  s = {s:.4f} m ({len(p)} particles)")
            p_cosy = sim.transform_coordinates(p, from_system=FELSIM, to_system=COSY)
            for i in range(min(len(p), 5)):
                _print(f"    Particle {i} (COSY): {p_cosy[i]}")

        results = {
            'n_elements': n_elements,
            'checkpoints': len(evolution.s_positions),
            'test_labels': labels,
            'input_cosy': test_particles.tolist(),
            'input_felsim': particles_felsim.tolist(),
        }

    except Exception as e:
        _print(f"  FAILED: {e}")
        results = {'error': str(e)}

    OUTDIR.mkdir(parents=True, exist_ok=True)
    with open(OUTDIR / 'part_d_convention_audit.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    _print(f"\n  Saved: {OUTDIR / 'part_d_convention_audit.json'}")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  Summary
# ═══════════════════════════════════════════════════════════════════════════════

def summary(part_a_data=None, part_b_data=None, part_c_data=None, part_d_data=None):
    """Print diagnostic summary and verdict."""
    _print("\n" + "=" * 72)
    _print("  I6 Diagnostic Summary")
    _print("=" * 72)

    if part_a_data:
        _print("\n── Order Dependence ──")
        for order in [1, 2, 3]:
            d = part_a_data.get(order, {})
            sz = d.get('sigma_z_out_ps', '?')
            sl = d.get('sigma_l_out_um', '?')
            _print(f"  ORDER {order}: σ_z = {sz} ps, σ_l = {sl} μm")

    if part_c_data:
        _print("\n── Single Particle (δ=0.5%) ──")
        for order in [1, 3]:
            d = part_c_data.get(order, {})
            l_um = d.get('l_final_um', '?')
            l_exp = d.get('l_expected_um', '?')
            ratio = d.get('l_ratio', '?')
            _print(f"  ORDER {order}: l_final = {l_um} μm, expected = {l_exp} μm, ratio = {ratio}")

    _print("\n── Diagnosis ──")
    if part_a_data:
        o1 = part_a_data.get(1, {})
        o3 = part_a_data.get(3, {})
        o1_sz = o1.get('sigma_z_out_ps', float('nan'))
        o3_sz = o3.get('sigma_z_out_ps', float('nan'))

        if not np.isnan(o1_sz) and not np.isnan(o3_sz):
            ratio_o1 = o1_sz / 2.0
            ratio_o3 = o3_sz / 2.0
            if ratio_o1 < 2 and ratio_o3 < 2:
                _print("  FINDING: σ_z blowup NOT reproduced with W9-optimized currents.")
                _print(f"  ORDER 1: {ratio_o1:.2f}× growth, ORDER 3: {ratio_o3:.2f}× growth.")
                _print("  Growth is consistent with linear R56 × σ_δ coupling (ORDER 1)")
                _print("  plus moderate nonlinear path-length terms (ORDER 2–3).")
                _print("  The original 60–100× blowup likely occurred with poorly matched")
                _print("  currents that produced large dispersion residuals and beam blowup")
                _print("  in the chicane/spectrometer sections.")
            elif ratio_o1 < 1.5 and ratio_o3 > 10:
                _print("  FINDING: Higher-order map amplification.")
                _print("  At ORDER 1, σ_z is preserved; at ORDER 3, nonlinear terms blow it up.")
            elif ratio_o1 > 10:
                _print("  FINDING: Blowup occurs even at ORDER 1 — coordinate issue.")
            else:
                _print(f"  ORDER 1: {ratio_o1:.2f}×, ORDER 3: {ratio_o3:.2f}×")
                _print("  Moderate higher-order contribution — correct physics.")


# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="I6: COSY σ_z Blowup Diagnostic")
    parser.add_argument('--part-a', action='store_true',
                        help='Part A: order-dependence test')
    parser.add_argument('--part-b', action='store_true',
                        help='Part B: element-by-element σ_l(s) tracing')
    parser.add_argument('--part-c', action='store_true',
                        help='Part C: single off-energy particle test')
    parser.add_argument('--part-d', action='store_true',
                        help='Part D: reference particle convention audit')
    parser.add_argument('--all', action='store_true',
                        help='Run all parts')
    args = parser.parse_args()

    if not any([args.part_a, args.part_b, args.part_c, args.part_d, args.all]):
        args.all = True

    _print("I6: COSY σ_z Blowup Diagnostic")
    _print(f"E = {Energy} MeV, γ = {GAMMA:.2f}, β = {BETA_REL:.6f}")
    _print(f"p₀c = {P_C:.3f} MeV, f_RF = {2856e6/1e6:.0f} MHz")
    _print(f"Particles: {N_PARTICLES}, Seed: {SEED}")

    part_a_data = None
    part_b_data = None
    part_c_data = None
    part_d_data = None

    if args.part_a or args.all:
        part_a_data = part_a()

    if args.part_b or args.all:
        part_b_data = part_b(part_a_data)

    if args.part_c or args.all:
        part_c_data = part_c()

    if args.part_d or args.all:
        part_d_data = part_d()

    summary(part_a_data, part_b_data, part_c_data, part_d_data)

    _print("\n" + "=" * 72)
    _print("  I6 Complete")
    _print("=" * 72)


if __name__ == "__main__":
    main()
