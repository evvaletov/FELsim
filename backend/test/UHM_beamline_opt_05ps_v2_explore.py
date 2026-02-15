"""
Explore optimization strategies for 0.5 ps bunch spread with paper-aligned Twiss targets.

Compares strategies for the final undulator matching stage with the asymmetric
Twiss targets from Weinberg, Fisher & Li (arXiv:2510.14061v1, Table I):
    Horizontal:  beta_x = 1.4 m,   alpha_x = 0.47
    Vertical:    beta_y = 0.24 m,   alpha_y = 0

Strategies:
  1. NM 3-var (symmetric targets) — baseline reference
  2. NM 3-var (asymmetric targets) — demonstrates overconstraint
  3. NM 4-var joint chrom5+triplet (asymmetric) — expected solution
  4. DiffEvo 3-var (asymmetric) — global search, 3-var limit
  5. DiffEvo 4-var joint (asymmetric) — global verification

Author: Eremey Valetov
Date: 2026-02-06
"""

import sys
import time
from pathlib import Path
import numpy as np
import scipy.optimize as spo

from ebeam import beam
from beamline import lattice
from excelElements import ExcelElements
from beamOptimizer import beamOptimizer

# ── Beam Parameters ──────────────────────────────────────────────────────────

Energy = 40
f = 2856e6
bunch_spread = 0.5      # ps
energy_std_percent = 0.5  # % — unchanged from baseline per arXiv:2510.14061v1
h = 5e9                 # 1/s — unchanged from baseline

epsilon_n = 8
x_std = 0.8
y_std = 0.8
nb_particles = 1000

np.random.seed(42)

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

# ── Undulator matching targets ───────────────────────────────────────────────

K = 1.2
lambda_u = 2.3e-2  # m
N_u = 47

# Vertical (natural undulator focusing)
beta_ym = relat.gamma * lambda_u / (2 * np.pi * K)
alpha_ym = 0.0

# Horizontal (radiation mode matching, Table I)
beta_xm = 1.4
alpha_xm = 0.47

# Symmetric target (used in original script, for comparison)
beta_sym = beta_ym  # 0.2418 m
alpha_sym = 0.0

print(f"Symmetric target: beta = {beta_sym:.4f} m, alpha = 0")
print(f"Asymmetric targets: beta_x = {beta_xm}, alpha_x = {alpha_xm}, beta_y = {beta_ym:.4f}, alpha_y = 0")
print(f"gamma = {relat.gamma:.4f}, beta = {relat.beta:.6f}")

# ── Load beamline ────────────────────────────────────────────────────────────

file_path = Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx'
excel = ExcelElements(file_path)
beamlineUH = excel.create_beamline()
line_template = relat.changeBeamType("electron", Energy, beamlineUH)

segments = 118
line_template = line_template[:segments]


# ── Upstream optimization (stages 1–10) ─────────────────────────────────────

def run_upstream_stages(line, beam_dist):
    opti = beamOptimizer(line, beam_dist)

    opti.calc("Nelder-Mead",
              {1: ["I", "current", lambda n: n], 3: ["I2", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}, "I2": {"bounds": (0, 10), "start": 1}},
              {8: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                   {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.0}],
               9: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                   {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]},
              printResults=False)

    opti.calc("Nelder-Mead",
              {10: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {15: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]},
              printResults=False)

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

    opti.calc("Nelder-Mead",
              {27: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {32: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]},
              printResults=False)

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

    opti.calc("Nelder-Mead",
              {50: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {55: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]},
              printResults=False)

    opti.calc("Nelder-Mead",
              {56: ["I", "current", lambda n: n], 58: ["I2", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2}},
              {59: [{"measure": ["x", "envelope"], "goal": 0.0, "weight": 1},
                    {"measure": ["y", "envelope"], "goal": 0.0, "weight": 1}]},
              printResults=False)

    opti.calc("Nelder-Mead",
              {61: ["I", "current", lambda n: n], 63: ["I2", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2}, "I2": {"bounds": (0, 10), "start": 2}},
              {68: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
               69: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                    {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]},
              printResults=False)

    opti.calc("Nelder-Mead",
              {70: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {75: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]},
              printResults=False)

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

    return opti


def eval_twiss(line, beam_dist, bx_target, ax_target, by_target, ay_target):
    """Evaluate Twiss at element 117 and return summary dict."""
    eb = beam()
    particles = beam_dist.copy()
    for seg in line[:118]:
        particles = np.array(seg.useMatrice(particles))

    xa = eb.alpha(particles, "x")
    ya = eb.alpha(particles, "y")
    xb = eb.beta(particles, "x")
    yb = eb.beta(particles, "y")

    mse = ((xa - ax_target)**2 + (ya - ay_target)**2 +
           (xb - bx_target)**2 + (yb - by_target)**2) / 4

    return {
        "x_alpha": xa, "y_alpha": ya,
        "x_beta": xb, "y_beta": yb,
        "mse": mse,
        "currents": (line[93].current, line[95].current, line[97].current),
        "I87": line[87].current,
    }


def print_result(label, r, bx_t, by_t, ax_t, ay_t):
    print(f"\n{'─'*70}")
    print(f"  {label}")
    print(f"{'─'*70}")
    print(f"  x alpha = {r['x_alpha']:+.6f}  (target {ax_t:+.4f})")
    print(f"  y alpha = {r['y_alpha']:+.6f}  (target {ay_t:+.4f})")
    print(f"  x beta  = {r['x_beta']:.4f} m  (target {bx_t:.4f})")
    print(f"  y beta  = {r['y_beta']:.4f} m  (target {by_t:.4f})")
    print(f"  MSE     = {r['mse']:.6e}")
    print(f"  Currents: I87={r['I87']:.4f}, I93={r['currents'][0]:.4f}, "
          f"I95={r['currents'][1]:.4f}, I97={r['currents'][2]:.4f}")


# ── Run upstream stages ──────────────────────────────────────────────────────

print("Running upstream stages 1-10...")
line_base = list(line_template)
opti_base = run_upstream_stages(line_base, beam_dist)
print("Upstream stages complete.\n")

upstream_currents = {}
for i in range(len(line_base)):
    if hasattr(line_base[i], 'current'):
        upstream_currents[i] = line_base[i].current

results = {}


def restore_upstream():
    for i, cur in upstream_currents.items():
        line_base[i].current = cur


# ══════════════════════════════════════════════════════════════════════════════
# Strategy 1: NM 3-var, symmetric targets (baseline reference)
# ══════════════════════════════════════════════════════════════════════════════

def strat_1():
    restore_upstream()
    # Run chrom5 independently first
    opti = beamOptimizer(line_base, beam_dist)
    opti.calc("Nelder-Mead",
              {87: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {92: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]},
              printResults=False)
    # Then 3-var triplet with symmetric targets
    opti.calc("Nelder-Mead",
              {93: ["I", "current", lambda n: n], 95: ["I2", "current", lambda n: n],
               97: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2},
               "I2": {"bounds": (0, 10), "start": 2},
               "I3": {"bounds": (0, 10), "start": 2}},
              {117: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
                     {"measure": ["x", "beta"], "goal": beta_sym, "weight": 1},
                     {"measure": ["y", "beta"], "goal": beta_sym, "weight": 1}]},
              printResults=False)

t0 = time.perf_counter()
strat_1()
dt = time.perf_counter() - t0
r1 = eval_twiss(line_base, beam_dist, beta_sym, alpha_sym, beta_sym, alpha_sym)
r1["time"] = dt
print_result("1. NM 3-var, symmetric targets (reference)", r1,
             beta_sym, beta_sym, alpha_sym, alpha_sym)
print(f"  Time: {dt:.1f} s")
results["1"] = r1


# ══════════════════════════════════════════════════════════════════════════════
# Strategy 2: NM 3-var, asymmetric paper targets
# ══════════════════════════════════════════════════════════════════════════════

def strat_2():
    restore_upstream()
    opti = beamOptimizer(line_base, beam_dist)
    opti.calc("Nelder-Mead",
              {87: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {92: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]},
              printResults=False)
    opti.calc("Nelder-Mead",
              {93: ["I", "current", lambda n: n], 95: ["I2", "current", lambda n: n],
               97: ["I3", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 2},
               "I2": {"bounds": (0, 10), "start": 2},
               "I3": {"bounds": (0, 10), "start": 2}},
              {117: [{"measure": ["x", "alpha"], "goal": alpha_xm, "weight": 1},
                     {"measure": ["y", "alpha"], "goal": alpha_ym, "weight": 1},
                     {"measure": ["x", "beta"], "goal": beta_xm, "weight": 1},
                     {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]},
              printResults=False)

t0 = time.perf_counter()
strat_2()
dt = time.perf_counter() - t0
r2 = eval_twiss(line_base, beam_dist, beta_xm, alpha_xm, beta_ym, alpha_ym)
r2["time"] = dt
print_result("2. NM 3-var, asymmetric targets", r2,
             beta_xm, beta_ym, alpha_xm, alpha_ym)
print(f"  Time: {dt:.1f} s")
results["2"] = r2


# ══════════════════════════════════════════════════════════════════════════════
# Strategy 3: NM 4-var joint chrom5+triplet, asymmetric targets
# ══════════════════════════════════════════════════════════════════════════════

def strat_3():
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
               117: [{"measure": ["x", "alpha"], "goal": alpha_xm, "weight": 1},
                     {"measure": ["y", "alpha"], "goal": alpha_ym, "weight": 1},
                     {"measure": ["x", "beta"], "goal": beta_xm, "weight": 1},
                     {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}]},
              printResults=False)

t0 = time.perf_counter()
strat_3()
dt = time.perf_counter() - t0
r3 = eval_twiss(line_base, beam_dist, beta_xm, alpha_xm, beta_ym, alpha_ym)
r3["time"] = dt
print_result("3. NM 4-var joint, asymmetric targets", r3,
             beta_xm, beta_ym, alpha_xm, alpha_ym)
print(f"  Time: {dt:.1f} s")
results["3"] = r3


# ══════════════════════════════════════════════════════════════════════════════
# Strategy 4: DiffEvo 3-var, asymmetric targets
# ══════════════════════════════════════════════════════════════════════════════

def strat_4():
    restore_upstream()
    # Run chrom5 independently first
    opti = beamOptimizer(line_base, beam_dist)
    opti.calc("Nelder-Mead",
              {87: ["I", "current", lambda n: n]},
              {"I": {"bounds": (0, 10), "start": 1}},
              {92: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]},
              printResults=False)

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
        return ((xa - alpha_xm)**2 + (ya - alpha_ym)**2 +
                (xb - beta_xm)**2 + (yb - beta_ym)**2) / 4

    spo.differential_evolution(
        objective, bounds=[(0.01, 10), (0.01, 10), (0.01, 10)],
        seed=42, maxiter=200, tol=1e-10, polish=True)

t0 = time.perf_counter()
strat_4()
dt = time.perf_counter() - t0
r4 = eval_twiss(line_base, beam_dist, beta_xm, alpha_xm, beta_ym, alpha_ym)
r4["time"] = dt
print_result("4. DiffEvo 3-var, asymmetric targets", r4,
             beta_xm, beta_ym, alpha_xm, alpha_ym)
print(f"  Time: {dt:.1f} s")
results["4"] = r4


# ══════════════════════════════════════════════════════════════════════════════
# Strategy 5: DiffEvo 4-var joint chrom5+triplet, asymmetric targets
# ══════════════════════════════════════════════════════════════════════════════

def strat_5():
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
        return ((xa - alpha_xm)**2 + (ya - alpha_ym)**2 +
                (xb - beta_xm)**2 + (yb - beta_ym)**2) / 4

    spo.differential_evolution(
        objective, bounds=[(0.01, 10), (0.01, 10), (0.01, 10), (0.01, 10)],
        seed=42, maxiter=200, tol=1e-10, polish=True)

t0 = time.perf_counter()
strat_5()
dt = time.perf_counter() - t0
r5 = eval_twiss(line_base, beam_dist, beta_xm, alpha_xm, beta_ym, alpha_ym)
r5["time"] = dt
print_result("5. DiffEvo 4-var joint, asymmetric targets", r5,
             beta_xm, beta_ym, alpha_xm, alpha_ym)
print(f"  Time: {dt:.1f} s")
results["5"] = r5


# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 90)
print("  SUMMARY")
print("=" * 90)
print(f"  {'#':<3s} {'Strategy':<45s} {'MSE':>12s} {'x_beta':>8s} {'y_beta':>8s} "
      f"{'x_alpha':>9s} {'y_alpha':>9s} {'Time':>7s}")
print("-" * 90)

for key in sorted(results.keys()):
    r = results[key]
    lbl = {
        "1": "NM 3-var, symmetric (reference)",
        "2": "NM 3-var, asymmetric targets",
        "3": "NM 4-var joint, asymmetric targets",
        "4": "DiffEvo 3-var, asymmetric targets",
        "5": "DiffEvo 4-var joint, asymmetric targets",
    }[key]
    print(f"  {key:<3s} {lbl:<45s} {r['mse']:12.3e} {r['x_beta']:8.4f} {r['y_beta']:8.4f} "
          f"{r['x_alpha']:+9.5f} {r['y_alpha']:+9.5f} {r['time']:6.1f}s")

print(f"\n  Symmetric target:  beta_x = beta_y = {beta_sym:.4f} m, alpha_x = alpha_y = 0")
print(f"  Asymmetric target: beta_x = {beta_xm} m, alpha_x = {alpha_xm}, "
      f"beta_y = {beta_ym:.4f} m, alpha_y = {alpha_ym}")
