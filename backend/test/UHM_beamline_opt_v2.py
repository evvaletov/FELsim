# UH MkV FEL Beamline Optimization — 2 ps bunch, paper-aligned Twiss targets
#
# Revision of UHM_beamline_opt.py following Weinberg, Fisher & Li
# (arXiv:2510.14061v1).
#
# Beam parameters are identical to the original UHM_beamline_opt.py.  The only
# change is the undulator Twiss matching targets, now aligned with Table I:
#     Horizontal:  beta_x = 1.4 m,   alpha_x = 0.47  (radiation mode matching)
#     Vertical:    beta_y = 0.24 m,   alpha_y = 0     (natural undulator focusing)
#
# The original script targeted symmetric beta = 0.2418 m and alpha = 0 in both
# planes.  The horizontal matching to the radiation spot size (waist at undulator
# center) was not accounted for; this revision corrects that.
#
# Author: Eremey Valetov
# Date: 2026-02-06
# Reference: Weinberg, Fisher & Li, arXiv:2510.14061v1, Table I and §III

import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
import sympy as sp
import sympy.plotting as plot

from ebeam import beam
from beamline import lattice, beamline
from schematic import draw_beamline
from excelElements import ExcelElements
from beamOptimizer import beamOptimizer
from AlgebraicOptimization import AlgebraicOpti

# ── Beam parameters (identical to original UHM_beamline_opt.py) ──────────────
Energy = 40             # MeV
f = 2856e6              # Hz
bunch_spread = 2        # ps
energy_std_percent = 0.5  # %
h = 5e9                 # 1/s

epsilon_n = 8           # pi.mm.mrad
x_std = 0.8             # mm
y_std = 0.8             # mm
nb_particles = 1000
np.random.seed(42)

relat = lattice(1, fringeType=None)
relat.setE(E=Energy)
norm = relat.gamma * relat.beta
epsilon = epsilon_n / norm
print(f"gamma = {relat.gamma:.3f}")
print(f"beta  = {relat.beta:.6f}")
print(f"epsilon = {epsilon:.4f} pi.mm.mrad")
x_prime_std = epsilon / x_std
y_prime_std = epsilon / y_std

tof_std = bunch_spread * 1e-9 * f
energy_std = energy_std_percent * 10

ebeam = beam()
beam_dist = ebeam.gen_6d_gaussian(0, [x_std, x_prime_std, y_std, y_prime_std, tof_std, energy_std], nb_particles)
tof_dist = beam_dist[:, 4] / f
beam_dist[:, 5] += h * tof_dist

# ── Undulator matching targets ───────────────────────────────────────────────
# From Weinberg, Fisher & Li, arXiv:2510.14061v1, Table I.
K = 1.2
lambda_u = 2.3  # cm
N_u = 47

# Vertical: beam at waist, collimated through undulator (Eq. from [36])
beta_ym = relat.gamma * (lambda_u * 1e-2) / (2 * np.pi * K)
alpha_ym = 0.0
print(f"beta_y (matched) = {beta_ym:.4f} m")

# Horizontal: beam matched to radiation mode with waist at undulator center.
# Table I gives beta_x = 1.4 m, alpha_x = 0.47 rad at the undulator entrance.
beta_xm = 1.4   # m
alpha_xm = 0.47  # rad
print(f"beta_x (matched) = {beta_xm} m")
print(f"alpha_x (matched) = {alpha_xm} rad")

L_u = N_u * lambda_u * 1e-2
gamma_tw = (1 + alpha_xm**2) / beta_xm
z_waist = alpha_xm / gamma_tw
print(f"L_u = {L_u:.3f} m, L_u/2 = {L_u/2:.3f} m, z_waist = {z_waist:.3f} m")

# ── Load beamline ────────────────────────────────────────────────────────────
file_path = Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx'
excel = ExcelElements(file_path)
df = excel.get_dataframe()
beamlineUH = excel.create_beamline()
schem = draw_beamline()
line_UH = relat.changeBeamType("electron", Energy, beamlineUH)

print(f"Number of elements in beamline: {len(line_UH)}")
segments = 118
line = line_UH[:segments]
opti = beamOptimizer(line, beam_dist)

# ── Stage 1: First Quadrupole Doublet ────────────────────────────────────────

variables = {
    1: ["I", "current", lambda num: num],
    3: ["I2", "current", lambda num: num],
}
startPoint = {
    "I": {"bounds": (0, 10), "start": 1},
    "I2": {"bounds": (0, 10), "start": 1},
}
objectives = {
    8: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
        {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.0}],
    9: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
        {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]
}
result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                    plotBeam=False, printResults=True, plotProgress=False)

# ── Stage 2: First Chromaticity Quad ─────────────────────────────────────────

variables = {10: ["I", "current", lambda num: num]}
startPoint = {"I": {"bounds": (0, 10), "start": 1}}
objectives = {15: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}
result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                    plotBeam=False, printResults=True, plotProgress=False)

# ── Stage 3: Quadrupole Triplet 1 ───────────────────────────────────────────

variables = {
    16: ["I", "current", lambda num: num],
    18: ["I2", "current", lambda num: num],
    20: ["I3", "current", lambda num: num],
}
startPoint = {
    "I": {"bounds": (0, 10), "start": 2},
    "I2": {"bounds": (0, 10), "start": 5},
    "I3": {"bounds": (0, 10), "start": 3},
}
objectives = {
    25: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
         {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
    26: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
         {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]
}
result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                    plotBeam=False, printResults=True, plotProgress=False)

# ── Stage 4: Second Chromaticity Quad ────────────────────────────────────────

variables = {27: ["I", "current", lambda num: num]}
startPoint = {"I": {"bounds": (0, 10), "start": 1}}
objectives = {32: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}
result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                    plotBeam=False, printResults=True, plotProgress=False)

# ── Stage 5: Double Quadrupole Triplet ───────────────────────────────────────

variables = {
    37: ["I", "current", lambda num: num],
    35: ["I2", "current", lambda num: num],
    33: ["I3", "current", lambda num: num],
}
startPoint = {
    "I": {"bounds": (0, 10), "start": 2},
    "I2": {"bounds": (0, 10), "start": 2},
    "I3": {"bounds": (0, 10), "start": 2},
}
objectives = {
    37: [
        {"measure": ["x", "alpha"], "goal": 0, "weight": 1},
        {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
        {"measure": ["x", "envelope"], "goal": 2.0, "weight": 1},
        {"measure": ["y", "envelope"], "goal": 2.0, "weight": 1}
    ]
}
result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                    plotBeam=False, printResults=True, plotProgress=False)

line[43].current = line[33].current
line[41].current = line[35].current
line[39].current = line[37].current

# ── Stage 6: Third Chromaticity Quad ─────────────────────────────────────────

variables = {50: ["I", "current", lambda num: num]}
startPoint = {"I": {"bounds": (0, 10), "start": 1}}
objectives = {55: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}
result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                    plotBeam=False, printResults=True, plotProgress=False)

# ── Stage 7: IP Doublet (z = 7.11 m, element 59) ────────────────────────────

variables = {
    56: ["I", "current", lambda num: num],
    58: ["I2", "current", lambda num: num],
}
startPoint = {
    "I": {"bounds": (0, 10), "start": 2},
    "I2": {"bounds": (0, 10), "start": 2},
}
objectives = {
    59: [
        {"measure": ["x", "envelope"], "goal": 0.0, "weight": 1},
        {"measure": ["y", "envelope"], "goal": 0.0, "weight": 1}
    ]
}
result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                    plotBeam=False, printResults=True, plotProgress=False)

# ── Stage 8: Post-IP Doublet ─────────────────────────────────────────────────

variables = {
    61: ["I", "current", lambda num: num],
    63: ["I2", "current", lambda num: num],
}
startPoint = {
    "I": {"bounds": (0, 10), "start": 2},
    "I2": {"bounds": (0, 10), "start": 2},
}
objectives = {
    68: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
         {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
    69: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
         {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]
}
result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                    plotBeam=False, printResults=True, plotProgress=False)

# ── Stage 9: Fourth Chromaticity Quad ────────────────────────────────────────

variables = {70: ["I", "current", lambda num: num]}
startPoint = {"I": {"bounds": (0, 10), "start": 1}}
objectives = {75: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}
result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                    plotBeam=False, printResults=True, plotProgress=False)

# ── Stage 10: Quadrupole Triplet 2 ──────────────────────────────────────────

variables = {
    76: ["I", "current", lambda num: num],
    78: ["I2", "current", lambda num: num],
    80: ["I3", "current", lambda num: num],
}
startPoint = {
    "I": {"bounds": (0, 10), "start": 2},
    "I2": {"bounds": (0, 10), "start": 2},
    "I3": {"bounds": (0, 10), "start": 2},
}
objectives = {
    85: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
         {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
    86: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
         {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]
}
result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                    plotBeam=False, printResults=True, plotProgress=False)

# ── Stage 11: Chromaticity 5 + Final Triplet — Undulator Matching ────────────
# Joint optimization: chromaticity quad 5 (index 87) + final triplet (93, 95, 97).
# 4 variables for 4 asymmetric Twiss objectives + 1 secondary dispersion goal.
# The asymmetric targets (beta_x != beta_y, alpha_x != 0) make the 3-variable
# triplet-only problem overconstrained; adding chromaticity quad 5 resolves this.

variables = {
    87: ["Ic", "current", lambda num: num],
    93: ["I", "current", lambda num: num],
    95: ["I2", "current", lambda num: num],
    97: ["I3", "current", lambda num: num],
}
startPoint = {
    "Ic": {"bounds": (0, 10), "start": 4},
    "I": {"bounds": (0, 10), "start": 2},
    "I2": {"bounds": (0, 10), "start": 2},
    "I3": {"bounds": (0, 10), "start": 2},
}
objectives = {
    92: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 0.5}],
    117: [
        {"measure": ["x", "alpha"], "goal": alpha_xm, "weight": 1},
        {"measure": ["y", "alpha"], "goal": alpha_ym, "weight": 1},
        {"measure": ["x", "beta"], "goal": beta_xm, "weight": 1},
        {"measure": ["y", "beta"], "goal": beta_ym, "weight": 1}
    ]
}
result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                    plotBeam=False, printResults=True, plotProgress=False)

# ── Display ──────────────────────────────────────────────────────────────────
acceptance = {"shape": 'circle', "radius": 0.1, "origin": [0, 0]}
schem.plotBeamPositionTransform(beam_dist, line, 0.01, plot=True, showIndice=False,
                                defineLim=False, saveFig=7.11, shape=acceptance,
                                matchScaling=False, scatter=True)
