"""
Explore alternative optimization strategies for the 0.5 ps bunch spread study.

Focus: improve undulator matching (final triplet, element indices 93/95/97)
where the baseline Nelder-Mead achieved x-beta = 0.085 m vs target 0.242 m.

Strategies tried:
  1. Alternative scipy methods (Powell, L-BFGS-B, COBYLA)
  2. Multiple random starting points with Nelder-Mead
  3. Joint optimization of chromaticity quad 5 + final triplet (4 variables)
  4. Different objective weighting (emphasize x-beta)
  5. Wider current bounds (0-15 A)
  6. Global optimization via scipy.optimize.differential_evolution

Author: Eremey Valetov
"""

import sys
import copy
import math
import time
from pathlib import Path
import numpy as np
import scipy.optimize as spo

from ebeam import beam
from beamline import lattice
from excelElements import ExcelElements
from beamOptimizer import beamOptimizer

# ── Beam Parameters (0.5 ps study) ──────────────────────────────────────────

Energy = 40  # MeV
f = 2856e6  # Hz
bunch_spread = 0.5  # ps
energy_std_percent = 2.0
h = 20e9  # 1/s

epsilon_n = 8  # pi.mm.mrad
x_std = 0.8  # mm
y_std = 0.8  # mm
nb_particles = 1000

np.random.seed(42)  # Reproducible particle distribution

relat = lattice(1, fringeType=None)
relat.setE(E=Energy)
norm = relat.gamma * relat.beta
epsilon = epsilon_n / norm
x_prime_std = epsilon / x_std
y_prime_std = epsilon / y_std

tof_std = bunch_spread * 1e-9 * f
energy_std = energy_std_percent * 10

ebeam_gen = beam()
beam_dist = ebeam_gen.gen_6d_gaussian(
    0, [x_std, x_prime_std, y_std, y_prime_std, tof_std, energy_std], nb_particles
)
tof_dist = beam_dist[:, 4] / f
beam_dist[:, 5] += h * tof_dist

# MkV undulator matching target
K = 1.2
lambda_u = 2.3e-2  # m
beta_ym = relat.gamma / (K * (2 * np.pi / lambda_u))

print(f"Target beta_ym = {beta_ym:.4f} m")
print(f"gamma = {relat.gamma:.4f}, beta = {relat.beta:.6f}")

# ── Load Beamline ───────────────────────────────────────────────────────────

file_path = Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx'
excel = ExcelElements(file_path)
beamlineUH = excel.create_beamline()
line_template = relat.changeBeamType("electron", Energy, beamlineUH)

segments = 118
line_template = line_template[:segments]

# ── Run stages 1–11 (same as baseline) ──────────────────────────────────────

def run_upstream_stages(line, beam_dist):
    """Run optimization stages 1–11, returning the optimizer object."""
    opti = beamOptimizer(line, beam_dist)

    # Stage 1: First Quadrupole Doublet
    opti.calc("Nelder-Mead",
              {1: ["I", "current", lambda n: n], 3: ["I2", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}, "I2": {"bounds": (0, 10), "start": 1}},
              {8: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                   {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.0}],
               9: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                   {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]},
              printResults=False)

    # Stage 2: First Chromaticity Quad
    opti.calc("Nelder-Mead",
              {10: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {15: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]},
              printResults=False)

    # Stage 3: Quadrupole Triplet
    opti.calc("Nelder-Mead",
              {16: ["I", "current", lambda n: n], 18: ["I2", "current", lambda n: n],
               20: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 5},
               "I3": {"bounds": (0, 10), "start": 3}},
              {25: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
               26: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]},
              printResults=False)

    # Stage 4: Second Chromaticity Quad
    opti.calc("Nelder-Mead",
              {27: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {32: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]},
              printResults=False)

    # Stage 5: Double Triplet
    opti.calc("Nelder-Mead",
              {37: ["I", "current", lambda n: n], 35: ["I2", "current", lambda n: n],
               33: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2},
               "I3": {"bounds": (0, 10), "start": 2}},
              {37: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["x", "envelope"], "goal": 2.0, "weight": 1},
                    {"measure": ["y", "envelope"], "goal": 2.0, "weight": 1}]},
              printResults=False)

    line[43].current = line[33].current
    line[41].current = line[35].current
    line[39].current = line[37].current

    # Stage 6: Third Chromaticity Quad
    opti.calc("Nelder-Mead",
              {50: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {55: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]},
              printResults=False)

    # Stage 7: IP Doublet
    opti.calc("Nelder-Mead",
              {56: ["I", "current", lambda n: n], 58: ["I2", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2}},
              {59: [{"measure": ["x", "envelope"], "goal": 0.0, "weight": 1},
                    {"measure": ["y", "envelope"], "goal": 0.0, "weight": 1}]},
              printResults=False)

    # Stage 8: Post-IP Doublet
    opti.calc("Nelder-Mead",
              {61: ["I", "current", lambda n: n], 63: ["I2", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2}},
              {68: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
               69: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]},
              printResults=False)

    # Stage 9: Fourth Chromaticity Quad
    opti.calc("Nelder-Mead",
              {70: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {75: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]},
              printResults=False)

    # Stage 10: Quadrupole Triplet 2
    opti.calc("Nelder-Mead",
              {76: ["I", "current", lambda n: n], 78: ["I2", "current", lambda n: n],
               80: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2},
               "I3": {"bounds": (0, 10), "start": 2}},
              {85: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
               86: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]},
              printResults=False)

    # Stage 11: Fifth Chromaticity Quad
    opti.calc("Nelder-Mead",
              {87: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {92: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]},
              printResults=False)

    return opti


def propagate_to(line, beam_dist, end_idx):
    """Propagate particles through line[0:end_idx] and return them."""
    particles = beam_dist.copy()
    for seg in line[:end_idx]:
        particles = np.array(seg.useMatrice(particles))
    return particles


def eval_undulator_twiss(line, beam_dist, beta_target):
    """Evaluate Twiss parameters at element 117 and return a summary dict."""
    eb = beam()
    particles = beam_dist.copy()
    for seg in line[:118]:
        particles = np.array(seg.useMatrice(particles))

    x_alpha = eb.alpha(particles, "x")
    y_alpha = eb.alpha(particles, "y")
    x_beta = eb.beta(particles, "x")
    y_beta = eb.beta(particles, "y")

    mse = (x_alpha**2 + y_alpha**2 +
           (x_beta - beta_target)**2 +
           (y_beta - beta_target)**2) / 4

    return {
        "x_alpha": x_alpha, "y_alpha": y_alpha,
        "x_beta": x_beta, "y_beta": y_beta,
        "mse": mse,
        "currents": (line[93].current, line[95].current, line[97].current),
    }


def print_result(label, r):
    print(f"\n{'─'*60}")
    print(f"  {label}")
    print(f"{'─'*60}")
    print(f"  x alpha = {r['x_alpha']:+.6f}")
    print(f"  y alpha = {r['y_alpha']:+.6f}")
    print(f"  x beta  = {r['x_beta']:.4f} m  (target {beta_ym:.4f})")
    print(f"  y beta  = {r['y_beta']:.4f} m  (target {beta_ym:.4f})")
    print(f"  RMS     = {math.sqrt(r['mse']):.6e}")
    print(f"  Currents: I93={r['currents'][0]:.4f}, I95={r['currents'][1]:.4f}, I97={r['currents'][2]:.4f}")


# ═════════════════════════════════════════════════════════════════════════════
# Run upstream stages once, then snapshot the line state
# ═════════════════════════════════════════════════════════════════════════════

print("Running upstream stages 1-11...")
line_base = list(line_template)  # shallow copy of list
opti_base = run_upstream_stages(line_base, beam_dist)
print("Upstream stages complete.\n")

# Save the upstream-optimized currents so we can restore them for each trial
upstream_currents = {}
for i in range(len(line_base)):
    if hasattr(line_base[i], 'current'):
        upstream_currents[i] = line_base[i].current

best_result = None
best_label = None


def try_strategy(label, run_fn):
    """Run a strategy, evaluate, print, and track the best."""
    global best_result, best_label
    t0 = time.perf_counter()
    run_fn()
    dt = time.perf_counter() - t0
    r = eval_undulator_twiss(line_base, beam_dist, beta_ym)
    r["time"] = dt
    print_result(label, r)
    print(f"  Time: {dt:.1f} s")
    if best_result is None or r["mse"] < best_result["mse"]:
        best_result = r.copy()
        best_label = label
    return r


def restore_upstream():
    """Restore all quad currents to post-stage-11 values."""
    for i, cur in upstream_currents.items():
        line_base[i].current = cur


# ═════════════════════════════════════════════════════════════════════════════
# Strategy 1: Nelder-Mead baseline (reproduce for reference with seed=42)
# ═════════════════════════════════════════════════════════════════════════════

def strat_nelder_mead_baseline():
    restore_upstream()
    opti = beamOptimizer(line_base, beam_dist)
    opti.calc("Nelder-Mead",
              {93: ["I", "current", lambda n: n], 95: ["I2", "current", lambda n: n],
               97: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2},
               "I2": {"bounds": (0, 10), "start": 2},
               "I3": {"bounds": (0, 10), "start": 2}},
              {117: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["x", "beta"], "goal": beta_ym, "weight": 1},
                     {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]},
              printResults=False)

results = {}
results["1. Nelder-Mead baseline"] = try_strategy("1. Nelder-Mead baseline", strat_nelder_mead_baseline)

# ═════════════════════════════════════════════════════════════════════════════
# Strategy 2: Powell method
# ═════════════════════════════════════════════════════════════════════════════

def strat_powell():
    restore_upstream()
    opti = beamOptimizer(line_base, beam_dist)
    opti.calc("Powell",
              {93: ["I", "current", lambda n: n], 95: ["I2", "current", lambda n: n],
               97: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2},
               "I2": {"bounds": (0, 10), "start": 2},
               "I3": {"bounds": (0, 10), "start": 2}},
              {117: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["x", "beta"], "goal": beta_ym, "weight": 1},
                     {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]},
              printResults=False)

results["2. Powell"] = try_strategy("2. Powell", strat_powell)

# ═════════════════════════════════════════════════════════════════════════════
# Strategy 3: L-BFGS-B method
# ═════════════════════════════════════════════════════════════════════════════

def strat_lbfgsb():
    restore_upstream()
    opti = beamOptimizer(line_base, beam_dist)
    opti.calc("L-BFGS-B",
              {93: ["I", "current", lambda n: n], 95: ["I2", "current", lambda n: n],
               97: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2},
               "I2": {"bounds": (0, 10), "start": 2},
               "I3": {"bounds": (0, 10), "start": 2}},
              {117: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["x", "beta"], "goal": beta_ym, "weight": 1},
                     {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]},
              printResults=False)

results["3. L-BFGS-B"] = try_strategy("3. L-BFGS-B", strat_lbfgsb)

# ═════════════════════════════════════════════════════════════════════════════
# Strategy 4: COBYLA method
# ═════════════════════════════════════════════════════════════════════════════

def strat_cobyla():
    restore_upstream()
    opti = beamOptimizer(line_base, beam_dist)
    opti.calc("COBYLA",
              {93: ["I", "current", lambda n: n], 95: ["I2", "current", lambda n: n],
               97: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0.01, 10), "start": 2},
               "I2": {"bounds": (0.01, 10), "start": 2},
               "I3": {"bounds": (0.01, 10), "start": 2}},
              {117: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["x", "beta"], "goal": beta_ym, "weight": 1},
                     {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]},
              printResults=False)

results["4. COBYLA"] = try_strategy("4. COBYLA", strat_cobyla)

# ═════════════════════════════════════════════════════════════════════════════
# Strategy 5: Nelder-Mead with many random starting points
# ═════════════════════════════════════════════════════════════════════════════

def strat_multistart():
    restore_upstream()
    best_mse = float("inf")
    best_currents = None
    rng = np.random.RandomState(123)

    for trial in range(20):
        restore_upstream()
        s1, s2, s3 = rng.uniform(0.1, 9.9, 3)
        opti = beamOptimizer(line_base, beam_dist)
        opti.calc("Nelder-Mead",
                  {93: ["I", "current", lambda n: n], 95: ["I2", "current", lambda n: n],
                   97: ["I3", "current", lambda n: n]},
                  {"I": {"bounds": (0, 10), "start": float(s1)},
                   "I2": {"bounds": (0, 10), "start": float(s2)},
                   "I3": {"bounds": (0, 10), "start": float(s3)}},
                  {117: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                         {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                         {"measure": ["x", "beta"], "goal": beta_ym, "weight": 1},
                         {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]},
                  printResults=False)
        r = eval_undulator_twiss(line_base, beam_dist, beta_ym)
        if r["mse"] < best_mse:
            best_mse = r["mse"]
            best_currents = r["currents"]

    # Apply best
    restore_upstream()
    line_base[93].current, line_base[95].current, line_base[97].current = best_currents

results["5. Multi-start NM (20 trials)"] = try_strategy("5. Multi-start NM (20 trials)", strat_multistart)

# ═════════════════════════════════════════════════════════════════════════════
# Strategy 6: Nelder-Mead with higher x-beta weight
# ═════════════════════════════════════════════════════════════════════════════

def strat_xbeta_weight():
    restore_upstream()
    opti = beamOptimizer(line_base, beam_dist)
    opti.calc("Nelder-Mead",
              {93: ["I", "current", lambda n: n], 95: ["I2", "current", lambda n: n],
               97: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2},
               "I2": {"bounds": (0, 10), "start": 2},
               "I3": {"bounds": (0, 10), "start": 2}},
              {117: [{"measure": ["x", "alpha"], "goal": 0, "weight": 0.5},
                     {"measure": ["y", "alpha"], "goal": 0, "weight": 0.5},
                     {"measure": ["x", "beta"], "goal": beta_ym, "weight": 3},
                     {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]},
              printResults=False)

results["6. NM x-beta weight=3"] = try_strategy("6. NM x-beta weight=3", strat_xbeta_weight)

# ═════════════════════════════════════════════════════════════════════════════
# Strategy 7: Joint optimization of chromaticity quad 5 + final triplet
# ═════════════════════════════════════════════════════════════════════════════

def strat_joint_chrom5():
    restore_upstream()
    opti = beamOptimizer(line_base, beam_dist)
    opti.calc("Nelder-Mead",
              {87: ["Ic", "current", lambda n: n],
               93: ["I", "current", lambda n: n],
               95: ["I2", "current", lambda n: n],
               97: ["I3", "current", lambda n: n]},
              {"Ic": {"bounds": (0, 10), "start": 4},
               "I": {"bounds": (0, 10), "start": 2},
               "I2": {"bounds": (0, 10), "start": 2},
               "I3": {"bounds": (0, 10), "start": 2}},
              {92: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 0.5}],
               117: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["x", "beta"], "goal": beta_ym, "weight": 1},
                     {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]},
              printResults=False)

results["7. Joint chrom5+triplet (4 var)"] = try_strategy("7. Joint chrom5+triplet (4 var)", strat_joint_chrom5)

# ═════════════════════════════════════════════════════════════════════════════
# Strategy 8: Wider bounds (0-15 A) with Nelder-Mead
# ═════════════════════════════════════════════════════════════════════════════

def strat_wide_bounds():
    restore_upstream()
    opti = beamOptimizer(line_base, beam_dist)
    opti.calc("Nelder-Mead",
              {93: ["I", "current", lambda n: n], 95: ["I2", "current", lambda n: n],
               97: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 15), "start": 2},
               "I2": {"bounds": (0, 15), "start": 2},
               "I3": {"bounds": (0, 15), "start": 2}},
              {117: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["x", "beta"], "goal": beta_ym, "weight": 1},
                     {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]},
              printResults=False)

results["8. NM wide bounds (0-15 A)"] = try_strategy("8. NM wide bounds (0-15 A)", strat_wide_bounds)

# ═════════════════════════════════════════════════════════════════════════════
# Strategy 9: Joint triplet2 + chrom5 + final triplet (6 variables)
# ═════════════════════════════════════════════════════════════════════════════

def strat_joint_wide():
    restore_upstream()
    # Re-optimize triplet2 + chrom5 + final triplet together
    opti = beamOptimizer(line_base, beam_dist)
    opti.calc("Nelder-Mead",
              {76: ["Ia", "current", lambda n: n],
               78: ["Ib", "current", lambda n: n],
               80: ["Ic", "current", lambda n: n],
               87: ["Id", "current", lambda n: n],
               93: ["I", "current", lambda n: n],
               95: ["I2", "current", lambda n: n],
               97: ["I3", "current", lambda n: n]},
              {"Ia": {"bounds": (0, 10), "start": line_base[76].current},
               "Ib": {"bounds": (0, 10), "start": line_base[78].current},
               "Ic": {"bounds": (0, 10), "start": line_base[80].current},
               "Id": {"bounds": (0, 10), "start": line_base[87].current},
               "I": {"bounds": (0, 10), "start": 2},
               "I2": {"bounds": (0, 10), "start": 2},
               "I3": {"bounds": (0, 10), "start": 2}},
              {85: [{"measure": ["x", "alpha"], "goal": 0, "weight": 0.3},
                    {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.1}],
               86: [{"measure": ["y", "alpha"], "goal": 0, "weight": 0.3},
                    {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.1}],
               92: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 0.3}],
               117: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["x", "beta"], "goal": beta_ym, "weight": 1},
                     {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]},
              printResults=False)

results["9. Joint triplet2+chrom5+final (7 var)"] = try_strategy(
    "9. Joint triplet2+chrom5+final (7 var)", strat_joint_wide)

# ═════════════════════════════════════════════════════════════════════════════
# Strategy 10: Global optimization (differential evolution) on final triplet
# ═════════════════════════════════════════════════════════════════════════════

def strat_diffevo():
    restore_upstream()
    eb = beam()

    def objective(x):
        line_base[93].current = x[0]
        line_base[95].current = x[1]
        line_base[97].current = x[2]
        particles = beam_dist.copy()
        for seg in line_base[:118]:
            particles = np.array(seg.useMatrice(particles))
        xa = eb.alpha(particles, "x")
        ya = eb.alpha(particles, "y")
        xb = eb.beta(particles, "x")
        yb = eb.beta(particles, "y")
        return (xa**2 + ya**2 + (xb - beta_ym)**2 + (yb - beta_ym)**2) / 4

    result = spo.differential_evolution(
        objective, bounds=[(0.01, 10), (0.01, 10), (0.01, 10)],
        seed=42, maxiter=200, tol=1e-8, polish=True,
    )
    line_base[93].current = result.x[0]
    line_base[95].current = result.x[1]
    line_base[97].current = result.x[2]

results["10. Differential evolution"] = try_strategy("10. Differential evolution", strat_diffevo)

# ═════════════════════════════════════════════════════════════════════════════
# Strategy 11: Differential evolution on joint chrom5 + final triplet (4 var)
# ═════════════════════════════════════════════════════════════════════════════

def strat_diffevo_joint():
    restore_upstream()
    eb = beam()

    def objective(x):
        line_base[87].current = x[0]
        line_base[93].current = x[1]
        line_base[95].current = x[2]
        line_base[97].current = x[3]
        particles = beam_dist.copy()
        for seg in line_base[:118]:
            particles = np.array(seg.useMatrice(particles))
        xa = eb.alpha(particles, "x")
        ya = eb.alpha(particles, "y")
        xb = eb.beta(particles, "x")
        yb = eb.beta(particles, "y")
        return (xa**2 + ya**2 + (xb - beta_ym)**2 + (yb - beta_ym)**2) / 4

    result = spo.differential_evolution(
        objective, bounds=[(0.01, 10), (0.01, 10), (0.01, 10), (0.01, 10)],
        seed=42, maxiter=200, tol=1e-8, polish=True,
    )
    line_base[87].current = result.x[0]
    line_base[93].current = result.x[1]
    line_base[95].current = result.x[2]
    line_base[97].current = result.x[3]

results["11. DiffEvo joint chrom5+triplet"] = try_strategy(
    "11. DiffEvo joint chrom5+triplet", strat_diffevo_joint)

# ═════════════════════════════════════════════════════════════════════════════
# Strategy 12: Differential evolution on wide joint (triplet2+chrom5+final, 7 var)
# ═════════════════════════════════════════════════════════════════════════════

def strat_diffevo_wide():
    restore_upstream()
    eb = beam()

    def objective(x):
        line_base[76].current = x[0]
        line_base[78].current = x[1]
        line_base[80].current = x[2]
        line_base[87].current = x[3]
        line_base[93].current = x[4]
        line_base[95].current = x[5]
        line_base[97].current = x[6]
        particles = beam_dist.copy()
        for seg in line_base[:118]:
            particles = np.array(seg.useMatrice(particles))
        xa = eb.alpha(particles, "x")
        ya = eb.alpha(particles, "y")
        xb = eb.beta(particles, "x")
        yb = eb.beta(particles, "y")
        return (xa**2 + ya**2 + (xb - beta_ym)**2 + (yb - beta_ym)**2) / 4

    result = spo.differential_evolution(
        objective,
        bounds=[(0.01, 10), (0.01, 10), (0.01, 10), (0.01, 10),
                (0.01, 10), (0.01, 10), (0.01, 10)],
        seed=42, maxiter=300, tol=1e-8, polish=True,
    )
    line_base[76].current = result.x[0]
    line_base[78].current = result.x[1]
    line_base[80].current = result.x[2]
    line_base[87].current = result.x[3]
    line_base[93].current = result.x[4]
    line_base[95].current = result.x[5]
    line_base[97].current = result.x[6]

results["12. DiffEvo wide (7 var)"] = try_strategy("12. DiffEvo wide (7 var)", strat_diffevo_wide)


# ═════════════════════════════════════════════════════════════════════════════
# Strategy 13: DiffEvo on all downstream quads (post-IP through final, 10 var)
# Indices: 61,63 (post-IP) + 70 (chrom4) + 76,78,80 (triplet2) +
#          87 (chrom5) + 93,95,97 (final triplet)
# ═════════════════════════════════════════════════════════════════════════════

def strat_diffevo_all_downstream():
    restore_upstream()
    eb = beam()
    quad_idx = [61, 63, 70, 76, 78, 80, 87, 93, 95, 97]

    def objective(x):
        for j, idx in enumerate(quad_idx):
            line_base[idx].current = x[j]
        particles = beam_dist.copy()
        for seg in line_base[:118]:
            particles = np.array(seg.useMatrice(particles))
        xa = eb.alpha(particles, "x")
        ya = eb.alpha(particles, "y")
        xb = eb.beta(particles, "x")
        yb = eb.beta(particles, "y")
        return (xa**2 + ya**2 + (xb - beta_ym)**2 + (yb - beta_ym)**2) / 4

    starts = [upstream_currents[i] for i in quad_idx]
    result = spo.differential_evolution(
        objective, bounds=[(0.01, 10)] * len(quad_idx),
        seed=42, maxiter=500, tol=1e-8, polish=True,
    )
    for j, idx in enumerate(quad_idx):
        line_base[idx].current = result.x[j]

results["13. DiffEvo all downstream (10 var)"] = try_strategy(
    "13. DiffEvo all downstream (10 var)", strat_diffevo_all_downstream)


# ═════════════════════════════════════════════════════════════════════════════
# Summary
# ═════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  SUMMARY: All strategies ranked by RMS")
print("=" * 80)
print(f"  {'Strategy':<42s} {'RMS':>12s} {'x_beta':>8s} {'y_beta':>8s} {'x_alpha':>9s} {'y_alpha':>9s}")
print("-" * 80)

ranked = sorted(results.items(), key=lambda kv: kv[1]["mse"])
for label, r in ranked:
    print(f"  {label:<42s} {math.sqrt(r['mse']):12.6e} {r['x_beta']:8.4f} {r['y_beta']:8.4f} "
          f"{r['x_alpha']:+9.5f} {r['y_alpha']:+9.5f}")

print(f"\n  Target: x_beta = y_beta = {beta_ym:.4f} m, alpha = 0")
print(f"\n  BEST: {best_label}")
print(f"  Currents: I93={best_result['currents'][0]:.6f}, "
      f"I95={best_result['currents'][1]:.6f}, I97={best_result['currents'][2]:.6f}")
