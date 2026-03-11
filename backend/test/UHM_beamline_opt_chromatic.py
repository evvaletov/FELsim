"""UH MkV FEL Beamline Optimization — Chromatic + Apertures.

Re-optimization of UHM_beamline_opt_v2.py with all chromatic physics enabled:
- Chromatic quadrupole matrices (k_eff = k₀ × P₀/P)
- Chromatic dipole sector-bend body (ρ ∝ P)
- Chromatic dipole wedge edge kicks (R ∝ P)
- Aperture loss tracking (quad bore ±13.5 mm, dipole pole gap)
- Transmission as optimizer objective (weighted MSE term)

Warm-starts from felsim_nm_warm.json (linear NM solution) for comparison.

Author: Eremey Valetov
"""

import sys
import json
import time
from pathlib import Path
import numpy as np

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

from ebeam import beam
from beamline import lattice
from excelElements import ExcelElements
from beamOptimizer import beamOptimizer

# ── Beam parameters ──────────────────────────────────────────────────────────
Energy = 40
epsilon_n = 8           # pi.mm.mrad
x_std = 0.8             # mm
y_std = 0.8             # mm
f = 2856e6              # Hz
bunch_spread = 2        # ps
energy_std_percent = 0.5
h = 5e9                 # 1/s
nb_particles = 1000
SEGMENTS = 118

WARM_START_PATH = backend_dir / 'test' / 'results' / 'felsim_nm_warm.json'
OUTPUT_PATH = backend_dir / 'test' / 'results' / 'felsim_chromatic_warm.json'

# Transmission weight: w_T × (1-T)² added to MSE.
# Keep moderate — too high dominates Twiss objectives.
TRANSMISSION_WEIGHT = 5.0
N_PASSES = 3  # Coordinate-descent passes over 11-stage sequence

np.random.seed(42)


def load_warm_start():
    if WARM_START_PATH.exists():
        with open(WARM_START_PATH) as fh:
            data = json.load(fh)
        return {int(k): v for k, v in data['currents'].items()}
    return {}


def build_beam():
    relat = lattice(1, fringeType=None)
    relat.setE(E=Energy)
    norm = relat.gamma * relat.beta
    epsilon = epsilon_n / norm

    xp_std = epsilon / x_std
    yp_std = epsilon / y_std
    tof_std = bunch_spread * 1e-9 * f
    energy_std = energy_std_percent * 10  # ΔK/K₀ × 10³

    eb = beam()
    dist = eb.gen_6d_gaussian(0, [x_std, xp_std, y_std, yp_std, tof_std, energy_std],
                              nb_particles)
    tof_dist = dist[:, 4] / f
    dist[:, 5] += h * tof_dist

    # Undulator targets
    K = 1.2
    lambda_u = 2.3e-2  # m
    beta_ym = relat.gamma * lambda_u / (2 * np.pi * K)
    targets = {
        'beta_xm': 1.4, 'alpha_xm': 0.47,
        'beta_ym': beta_ym, 'alpha_ym': 0.0,
        'epsilon': epsilon,
    }
    return dist, relat, targets


def build_beamline(relat, chromatic=True):
    file_path = backend_dir.parent / 'beam_excel' / 'Beamline_elements.xlsx'
    excel = ExcelElements(str(file_path))
    bl = excel.create_beamline()
    line = relat.changeBeamType("electron", Energy, bl)[:SEGMENTS]
    if chromatic:
        for elem in line:
            elem.chromatic = True
    return line


def apply_currents(line, currents):
    for idx, current in currents.items():
        if idx < len(line):
            line[idx].current = abs(current)


def run_optimization(line, beam_dist, targets, warm_currents=None,
                     use_apertures=True, transmission_weight=0.0):
    """Run the full 11-stage sequential optimization."""
    axm = targets['alpha_xm']
    aym = targets['alpha_ym']
    bxm = targets['beta_xm']
    bym = targets['beta_ym']

    opti = beamOptimizer(line, beam_dist)

    def warm(idx, default=2.0):
        return abs(warm_currents.get(idx, default)) if warm_currents else default

    t0 = time.perf_counter()

    # ── Stage 1 ──────────────────────────────────────────────────────────────
    print("\n── Stage 1: First Quadrupole Doublet ──")
    variables = {
        1: ["I", "current", lambda num: num],
        3: ["I2", "current", lambda num: num],
    }
    startPoint = {
        "I": {"bounds": (0, 10), "start": warm(1, 1)},
        "I2": {"bounds": (0, 10), "start": warm(3, 1)},
    }
    objectives = {
        8: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
            {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.0}],
        9: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
            {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}],
    }
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              printResults=True, use_apertures=use_apertures)

    # ── Stage 2 ──────────────────────────────────────────────────────────────
    print("\n── Stage 2: First Chromaticity Quad ──")
    variables = {10: ["I", "current", lambda num: num]}
    startPoint = {"I": {"bounds": (0, 10), "start": warm(10, 1)}}
    objectives = {15: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              printResults=True, use_apertures=use_apertures)

    # ── Stage 3 ──────────────────────────────────────────────────────────────
    print("\n── Stage 3: Quadrupole Triplet 1 ──")
    variables = {
        16: ["I", "current", lambda num: num],
        18: ["I2", "current", lambda num: num],
        20: ["I3", "current", lambda num: num],
    }
    startPoint = {
        "I": {"bounds": (0, 10), "start": warm(16, 2)},
        "I2": {"bounds": (0, 10), "start": warm(18, 5)},
        "I3": {"bounds": (0, 10), "start": warm(20, 3)},
    }
    objectives = {
        25: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
             {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
        26: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
             {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}],
    }
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              printResults=True, use_apertures=use_apertures)

    # ── Stage 4 ──────────────────────────────────────────────────────────────
    print("\n── Stage 4: Second Chromaticity Quad ──")
    variables = {27: ["I", "current", lambda num: num]}
    startPoint = {"I": {"bounds": (0, 10), "start": warm(27, 1)}}
    objectives = {32: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              printResults=True, use_apertures=use_apertures)

    # ── Stage 5 ──────────────────────────────────────────────────────────────
    print("\n── Stage 5: Double Quadrupole Triplet ──")
    variables = {
        37: ["I", "current", lambda num: num],
        35: ["I2", "current", lambda num: num],
        33: ["I3", "current", lambda num: num],
    }
    startPoint = {
        "I": {"bounds": (0, 10), "start": warm(37, 2)},
        "I2": {"bounds": (0, 10), "start": warm(35, 2)},
        "I3": {"bounds": (0, 10), "start": warm(33, 2)},
    }
    objectives = {
        37: [
            {"measure": ["x", "alpha"], "goal": 0, "weight": 1},
            {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
            {"measure": ["x", "envelope"], "goal": 2.0, "weight": 1},
            {"measure": ["y", "envelope"], "goal": 2.0, "weight": 1},
        ]
    }
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              printResults=True, use_apertures=use_apertures,
              transmission_weight=transmission_weight)

    # Mirror symmetry
    line[43].current = line[33].current
    line[41].current = line[35].current
    line[39].current = line[37].current

    # ── Stage 6 ──────────────────────────────────────────────────────────────
    print("\n── Stage 6: Third Chromaticity Quad ──")
    variables = {50: ["I", "current", lambda num: num]}
    startPoint = {"I": {"bounds": (0, 10), "start": warm(50, 1)}}
    objectives = {55: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              printResults=True, use_apertures=use_apertures,
              transmission_weight=transmission_weight)

    # ── Stage 7 ──────────────────────────────────────────────────────────────
    print("\n── Stage 7: IP Doublet ──")
    variables = {
        56: ["I", "current", lambda num: num],
        58: ["I2", "current", lambda num: num],
    }
    startPoint = {
        "I": {"bounds": (0, 10), "start": warm(56, 2)},
        "I2": {"bounds": (0, 10), "start": warm(58, 2)},
    }
    objectives = {
        59: [
            {"measure": ["x", "envelope"], "goal": 0.0, "weight": 1},
            {"measure": ["y", "envelope"], "goal": 0.0, "weight": 1},
        ]
    }
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              printResults=True, use_apertures=use_apertures,
              transmission_weight=transmission_weight)

    # ── Stage 8 ──────────────────────────────────────────────────────────────
    print("\n── Stage 8: Post-IP Doublet ──")
    variables = {
        61: ["I", "current", lambda num: num],
        63: ["I2", "current", lambda num: num],
    }
    startPoint = {
        "I": {"bounds": (0, 10), "start": warm(61, 2)},
        "I2": {"bounds": (0, 10), "start": warm(63, 2)},
    }
    objectives = {
        68: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
             {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
        69: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
             {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}],
    }
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              printResults=True, use_apertures=use_apertures,
              transmission_weight=transmission_weight)

    # ── Stage 9 ──────────────────────────────────────────────────────────────
    print("\n── Stage 9: Fourth Chromaticity Quad ──")
    variables = {70: ["I", "current", lambda num: num]}
    startPoint = {"I": {"bounds": (0, 10), "start": warm(70, 1)}}
    objectives = {75: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              printResults=True, use_apertures=use_apertures,
              transmission_weight=transmission_weight)

    # ── Stage 10 ─────────────────────────────────────────────────────────────
    print("\n── Stage 10: Quadrupole Triplet 2 ──")
    variables = {
        76: ["I", "current", lambda num: num],
        78: ["I2", "current", lambda num: num],
        80: ["I3", "current", lambda num: num],
    }
    startPoint = {
        "I": {"bounds": (0, 10), "start": warm(76, 2)},
        "I2": {"bounds": (0, 10), "start": warm(78, 2)},
        "I3": {"bounds": (0, 10), "start": warm(80, 2)},
    }
    objectives = {
        85: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
             {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
        86: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
             {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}],
    }
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              printResults=True, use_apertures=use_apertures,
              transmission_weight=transmission_weight)

    # ── Stage 11 ─────────────────────────────────────────────────────────────
    print("\n── Stage 11: Chromaticity 5 + Final Triplet ──")
    variables = {
        87: ["Ic", "current", lambda num: num],
        93: ["I", "current", lambda num: num],
        95: ["I2", "current", lambda num: num],
        97: ["I3", "current", lambda num: num],
    }
    startPoint = {
        "Ic": {"bounds": (0, 10), "start": warm(87, 4)},
        "I": {"bounds": (0, 10), "start": warm(93, 2)},
        "I2": {"bounds": (0, 10), "start": warm(95, 2)},
        "I3": {"bounds": (0, 10), "start": warm(97, 2)},
    }
    objectives = {
        92: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 0.5}],
        117: [
            {"measure": ["x", "alpha"], "goal": axm, "weight": 1},
            {"measure": ["y", "alpha"], "goal": aym, "weight": 1},
            {"measure": ["x", "beta"], "goal": bxm, "weight": 1},
            {"measure": ["y", "beta"], "goal": bym, "weight": 1},
        ]
    }
    result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                        printResults=True, use_apertures=use_apertures,
                        transmission_weight=transmission_weight)

    dt = time.perf_counter() - t0
    print(f"\nTotal optimization time: {dt:.1f} s")

    return opti, result


def extract_currents(line):
    quad_indices = [1, 3, 10, 16, 18, 20, 27, 33, 35, 37, 39, 41, 43,
                    50, 56, 58, 61, 63, 70, 76, 78, 80, 87, 93, 95, 97]
    return {idx: line[idx].current for idx in quad_indices if idx < len(line)}


def compute_final_twiss(line, beam_dist, targets):
    """Track beam and compute final Twiss at undulator entrance."""
    ps = beam_dist.copy()
    for elem in line:
        ps = np.array(elem.useMatrice(ps))
    eb = beam()
    _, _, tw = eb.cal_twiss(ps, ddof=1)
    bx = tw.loc['x'][r"$\beta$ (m)"]
    ax = tw.loc['x'][r"$\alpha$"]
    by = tw.loc['y'][r"$\beta$ (m)"]
    ay = tw.loc['y'][r"$\alpha$"]

    mse = ((bx - targets['beta_xm'])**2 + (by - targets['beta_ym'])**2 +
           (ax - targets['alpha_xm'])**2 + (ay - targets['alpha_ym'])**2) / 4
    return {'beta_x': bx, 'alpha_x': ax, 'beta_y': by, 'alpha_y': ay, 'mse': mse}


def compute_transmission(line, beam_dist):
    """Track with aperture cuts, return survival fraction."""
    ps = beam_dist.copy()
    n0 = len(ps)
    for elem in line:
        ps = np.array(elem.useMatrice(ps))
        ps = elem.apply_aperture(ps)
    return len(ps) / n0


def save_results(currents, twiss, transmission, output_path):
    data = {
        'config': 'chromatic_aperture_NM',
        'energy_MeV': Energy,
        'epsilon_n': epsilon_n,
        'energy_spread_pct': energy_std_percent,
        'bunch_length_ps': bunch_spread,
        'nb_particles': nb_particles,
        'chromatic': True,
        'use_apertures': True,
        'transmission_weight': TRANSMISSION_WEIGHT,
        'transmission': transmission,
        'twiss_undulator': twiss,
        'currents': {str(k): float(v) for k, v in sorted(currents.items())},
    }
    with open(output_path, 'w') as fh:
        json.dump(data, fh, indent=2)
    print(f"Results saved to {output_path}")


if __name__ == "__main__":
    beam_dist, relat, targets = build_beam()

    warm_currents = load_warm_start()
    if warm_currents:
        print(f"Warm start: {len(warm_currents)} currents from {WARM_START_PATH.name}")
    else:
        print("Cold start (no warm start file found)")

    print(f"\nPhysics: chromatic=True, apertures=True, "
          f"transmission_weight={TRANSMISSION_WEIGHT}, passes={N_PASSES}")
    print(f"Beam: {Energy} MeV, ε_n={epsilon_n} π.mm.mrad, "
          f"σ_p={energy_std_percent}%, {bunch_spread} ps, "
          f"N={nb_particles}")

    current_warm = warm_currents
    best_mse = float('inf')
    best_currents = None

    for p in range(N_PASSES):
        print(f"\n{'#' * 70}")
        print(f"# PASS {p + 1}/{N_PASSES}")
        print(f"{'#' * 70}")

        line = build_beamline(relat, chromatic=True)
        if current_warm:
            apply_currents(line, current_warm)

        opti, result = run_optimization(
            line, beam_dist, targets,
            warm_currents=current_warm,
            use_apertures=True,
            transmission_weight=TRANSMISSION_WEIGHT,
        )

        currents = extract_currents(line)
        twiss = compute_final_twiss(line, beam_dist, targets)
        transmission = compute_transmission(line, beam_dist)

        print(f"\n  Pass {p + 1} summary: MSE={twiss['mse']:.6e}, T={transmission*100:.1f}%")

        if twiss['mse'] < best_mse:
            best_mse = twiss['mse']
            best_currents = currents.copy()
            best_twiss = twiss.copy()
            best_transmission = transmission

        # Feed back for next pass
        current_warm = currents

    # Apply best result
    print("\n" + "=" * 70)
    print(f"BEST RESULT (from {N_PASSES} passes)")
    print("=" * 70)
    print(f"\nOptimized currents:")
    ref = load_warm_start()
    for idx in sorted(best_currents):
        w = f" (linear: {ref.get(idx, 0):.4f})" if ref else ""
        print(f"  [{idx:3d}] {best_currents[idx]:8.4f} A{w}")

    print(f"\nFinal Twiss at undulator entrance (element 117):")
    print(f"  beta_x  = {best_twiss['beta_x']:8.4f} m   (target: {targets['beta_xm']:.4f})")
    print(f"  beta_y  = {best_twiss['beta_y']:8.4f} m   (target: {targets['beta_ym']:.4f})")
    print(f"  alpha_x = {best_twiss['alpha_x']:8.4f}     (target: {targets['alpha_xm']:.4f})")
    print(f"  alpha_y = {best_twiss['alpha_y']:8.4f}     (target: {targets['alpha_ym']:.4f})")
    print(f"  MSE     = {best_mse:.6e}")
    print(f"\n  Transmission: {best_transmission*100:.1f}%")
    print("=" * 70)

    save_results(best_currents, best_twiss, best_transmission, OUTPUT_PATH)
