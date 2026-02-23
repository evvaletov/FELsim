# UH MkV FEL Beamline Optimization — 0.5 ps Parameter Sensitivity Study
#
# 1D sweeps of (energy_std_percent, h, epsilon_n) at fixed 0.5 ps bunch length,
# mapping where undulator Twiss matching succeeds or degrades.
#
# Baseline: arXiv:2510.14061v1 (σ_E=0.5%, h=5e9, ε_n=8 π·mm·mrad)
#
# Author: Eremey Valetov
# Date: 2026-02-11

import sys
import time
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ebeam import beam
from beamline import lattice, beamline
from excelElements import ExcelElements
from beamOptimizer import beamOptimizer

# ── Constants ─────────────────────────────────────────────────────────────────
ENERGY = 40             # MeV
RF_FREQ = 2856e6        # Hz
SEGMENTS = 118

# Undulator matching targets (arXiv:2510.14061v1, Table I)
K_UND = 1.2
LAMBDA_U = 2.3e-2       # m

# Baseline beam parameters (v2 study)
BASELINE = {
    'bunch_spread': 0.5,       # ps
    'energy_std_percent': 0.5, # %
    'h': 5e9,                  # 1/s
    'epsilon_n': 8,            # π·mm·mrad
    'x_std': 0.8,              # mm
    'y_std': 0.8,              # mm
}

# Quad element indices used across all 11 stages
QUAD_INDICES = [1, 3, 10, 16, 18, 20, 27, 33, 35, 37, 39, 41, 43,
                50, 56, 58, 61, 63, 70, 76, 78, 80, 87, 93, 95, 97]


def compute_twiss_targets():
    """Compute undulator Twiss targets from paper parameters."""
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    beta_ym = relat.gamma * LAMBDA_U / (2 * np.pi * K_UND)
    alpha_ym = 0.0
    beta_xm = 1.4
    alpha_xm = 0.47
    return beta_xm, alpha_xm, beta_ym, alpha_ym, relat


def run_optimization(bunch_spread=0.5, energy_std_percent=0.5, h=5e9,
                     epsilon_n=8, x_std=0.8, y_std=0.8,
                     nb_particles=500, seed=42,
                     chrom_upper_bound=10, n_restarts=1,
                     stage11_method='Nelder-Mead', stage11_kwargs=None):
    """Run 11-stage beamline optimization, return results dict.

    Parameters
    ----------
    chrom_upper_bound : float
        Upper bound for chromaticity quad currents (default 10 A).
        Applies to quads 10, 27, 50, 70, 87 (stages 2, 4, 6, 9, 11).
        Non-chromaticity quads always use 10 A.
    n_restarts : int
        Number of random restarts for Stage 11. Best result is kept.
        Ignored when stage11_method='glyfada'.
    stage11_method : str
        Optimization method for Stage 11: 'Nelder-Mead' or 'glyfada'.
    stage11_kwargs : dict or None
        Extra kwargs for glyfada (pop_size, max_gen, sigma, etc.).

    Returns dict with keys: mse, alpha_x, alpha_y, beta_x, beta_y,
    disp_resid, quad_currents (dict of index->current), time_s, nfev, converged.
    """
    beta_xm, alpha_xm, beta_ym, alpha_ym, relat = compute_twiss_targets()
    norm = relat.gamma * relat.beta
    epsilon = epsilon_n / norm
    qb = 10              # standard quad bound
    cb = chrom_upper_bound  # chromaticity quad bound

    np.random.seed(seed)
    x_prime_std = epsilon / x_std
    y_prime_std = epsilon / y_std
    tof_std = bunch_spread * 1e-9 * RF_FREQ
    energy_std = energy_std_percent * 10

    ebeam_gen = beam()
    beam_dist = ebeam_gen.gen_6d_gaussian(
        0, [x_std, x_prime_std, y_std, y_prime_std, tof_std, energy_std],
        nb_particles)
    tof_dist = beam_dist[:, 4] / RF_FREQ
    beam_dist[:, 5] += h * tof_dist

    # Load beamline
    file_path = (Path(__file__).resolve().parent.parent.parent
                 / 'beam_excel' / 'Beamline_elements.xlsx')
    excel = ExcelElements(file_path)
    beamlineUH = excel.create_beamline()
    line_UH = relat.changeBeamType("electron", ENERGY, beamlineUH)
    line = line_UH[:SEGMENTS]

    opti = beamOptimizer(line, beam_dist)
    t0 = time.perf_counter()

    # Stage 1: First Quadrupole Doublet
    variables = {
        1: ["I", "current", lambda num: num],
        3: ["I2", "current", lambda num: num],
    }
    startPoint = {
        "I": {"bounds": (0, qb), "start": 1},
        "I2": {"bounds": (0, qb), "start": 1},
    }
    objectives = {
        8: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
            {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.0}],
        9: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
            {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]
    }
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              plotBeam=False, printResults=False, plotProgress=False)

    # Stage 2: First Chromaticity Quad
    variables = {10: ["I", "current", lambda num: num]}
    startPoint = {"I": {"bounds": (0, cb), "start": 1}}
    objectives = {15: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              plotBeam=False, printResults=False, plotProgress=False)

    # Stage 3: Quadrupole Triplet 1
    variables = {
        16: ["I", "current", lambda num: num],
        18: ["I2", "current", lambda num: num],
        20: ["I3", "current", lambda num: num],
    }
    startPoint = {
        "I": {"bounds": (0, qb), "start": 2},
        "I2": {"bounds": (0, qb), "start": 5},
        "I3": {"bounds": (0, qb), "start": 3},
    }
    objectives = {
        25: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
             {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
        26: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
             {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]
    }
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              plotBeam=False, printResults=False, plotProgress=False)

    # Stage 4: Second Chromaticity Quad
    variables = {27: ["I", "current", lambda num: num]}
    startPoint = {"I": {"bounds": (0, cb), "start": 1}}
    objectives = {32: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              plotBeam=False, printResults=False, plotProgress=False)

    # Stage 5: Double Quadrupole Triplet
    variables = {
        37: ["I", "current", lambda num: num],
        35: ["I2", "current", lambda num: num],
        33: ["I3", "current", lambda num: num],
    }
    startPoint = {
        "I": {"bounds": (0, qb), "start": 2},
        "I2": {"bounds": (0, qb), "start": 2},
        "I3": {"bounds": (0, qb), "start": 2},
    }
    objectives = {
        37: [
            {"measure": ["x", "alpha"], "goal": 0, "weight": 1},
            {"measure": ["y", "alpha"], "goal": 0, "weight": 1},
            {"measure": ["x", "envelope"], "goal": 2.0, "weight": 1},
            {"measure": ["y", "envelope"], "goal": 2.0, "weight": 1}
        ]
    }
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              plotBeam=False, printResults=False, plotProgress=False)
    line[43].current = line[33].current
    line[41].current = line[35].current
    line[39].current = line[37].current

    # Stage 6: Third Chromaticity Quad
    variables = {50: ["I", "current", lambda num: num]}
    startPoint = {"I": {"bounds": (0, cb), "start": 1}}
    objectives = {55: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              plotBeam=False, printResults=False, plotProgress=False)

    # Stage 7: IP Doublet
    variables = {
        56: ["I", "current", lambda num: num],
        58: ["I2", "current", lambda num: num],
    }
    startPoint = {
        "I": {"bounds": (0, qb), "start": 2},
        "I2": {"bounds": (0, qb), "start": 2},
    }
    objectives = {
        59: [
            {"measure": ["x", "envelope"], "goal": 0.0, "weight": 1},
            {"measure": ["y", "envelope"], "goal": 0.0, "weight": 1}
        ]
    }
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              plotBeam=False, printResults=False, plotProgress=False)

    # Stage 8: Post-IP Doublet
    variables = {
        61: ["I", "current", lambda num: num],
        63: ["I2", "current", lambda num: num],
    }
    startPoint = {
        "I": {"bounds": (0, qb), "start": 2},
        "I2": {"bounds": (0, qb), "start": 2},
    }
    objectives = {
        68: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
             {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
        69: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
             {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]
    }
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              plotBeam=False, printResults=False, plotProgress=False)

    # Stage 9: Fourth Chromaticity Quad
    variables = {70: ["I", "current", lambda num: num]}
    startPoint = {"I": {"bounds": (0, cb), "start": 1}}
    objectives = {75: [{"measure": ["x", "dispersion"], "goal": 0, "weight": 1}]}
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              plotBeam=False, printResults=False, plotProgress=False)

    # Stage 10: Quadrupole Triplet 2
    variables = {
        76: ["I", "current", lambda num: num],
        78: ["I2", "current", lambda num: num],
        80: ["I3", "current", lambda num: num],
    }
    startPoint = {
        "I": {"bounds": (0, qb), "start": 2},
        "I2": {"bounds": (0, qb), "start": 2},
        "I3": {"bounds": (0, qb), "start": 2},
    }
    objectives = {
        85: [{"measure": ["x", "alpha"], "goal": 0, "weight": 1},
             {"measure": ["x", "beta"], "goal": 0.1, "weight": 0.5}],
        86: [{"measure": ["y", "alpha"], "goal": 0, "weight": 1},
             {"measure": ["y", "beta"], "goal": 0.1, "weight": 0.5}]
    }
    opti.calc("Nelder-Mead", variables, startPoint, objectives,
              plotBeam=False, printResults=False, plotProgress=False)

    # Stage 11: Chromaticity 5 + Final Triplet — Undulator Matching
    variables = {
        87: ["Ic", "current", lambda num: num],
        93: ["I", "current", lambda num: num],
        95: ["I2", "current", lambda num: num],
        97: ["I3", "current", lambda num: num],
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

    if stage11_method == 'glyfada':
        s11_kw = {'pop_size': 30, 'max_gen': 20, 'sigma': 0.05}
        if stage11_kwargs:
            s11_kw.update(stage11_kwargs)
        startPoint = {
            "Ic": {"bounds": (0, cb), "start": 4},
            "I": {"bounds": (0, qb), "start": 2},
            "I2": {"bounds": (0, qb), "start": 2},
            "I3": {"bounds": (0, qb), "start": 2},
        }
        result = opti.calc("glyfada", variables, startPoint, objectives,
                           plotBeam=False, printResults=False,
                           plotProgress=False, **s11_kw)
    else:
        # NM multi-start for Stage 11
        default_starts = [
            {"Ic": 4, "I": 2, "I2": 2, "I3": 2},
        ]
        rng = np.random.RandomState(seed + 999)
        for _ in range(n_restarts - 1):
            default_starts.append({
                "Ic": rng.uniform(0, cb),
                "I": rng.uniform(0, qb),
                "I2": rng.uniform(0, qb),
                "I3": rng.uniform(0, qb),
            })

        pre_s11_currents = {idx: line[idx].current for idx in [87, 93, 95, 97]}
        best_result = None

        for attempt, starts in enumerate(default_starts):
            for idx, cur in pre_s11_currents.items():
                line[idx].current = cur

            startPoint = {
                "Ic": {"bounds": (0, cb), "start": starts["Ic"]},
                "I": {"bounds": (0, qb), "start": starts["I"]},
                "I2": {"bounds": (0, qb), "start": starts["I2"]},
                "I3": {"bounds": (0, qb), "start": starts["I3"]},
            }
            result = opti.calc("Nelder-Mead", variables, startPoint, objectives,
                                plotBeam=False, printResults=False, plotProgress=False)
            if best_result is None or result.fun < best_result.fun:
                best_result = result
                best_s11_currents = {idx: line[idx].current for idx in [87, 93, 95, 97]}

        for idx, cur in best_s11_currents.items():
            line[idx].current = cur
        result = best_result

    elapsed = time.perf_counter() - t0

    # Extract Twiss at undulator entrance
    eb = beam()
    particles = beam_dist.copy()
    for elem in line:
        particles = np.array(elem.useMatrice(particles))
    _, _, twiss = eb.cal_twiss(particles)

    alpha_x = twiss.loc['x'][r"$\alpha$"]
    alpha_y = twiss.loc['y'][r"$\alpha$"]
    beta_x = twiss.loc['x'][r"$\beta$ (m)"]
    beta_y = twiss.loc['y'][r"$\beta$ (m)"]

    # Dispersion at element 92
    particles_92 = beam_dist.copy()
    for elem in line[:93]:
        particles_92 = np.array(elem.useMatrice(particles_92))
    _, _, twiss_92 = eb.cal_twiss(particles_92)
    disp_resid = twiss_92.loc['x'][r"$D$ (m)"]

    # Collect quad currents
    quad_currents = {idx: line[idx].current for idx in QUAD_INDICES}

    return {
        'mse': result.fun,
        'alpha_x': alpha_x,
        'alpha_y': alpha_y,
        'beta_x': beta_x,
        'beta_y': beta_y,
        'disp_resid': disp_resid,
        'quad_currents': quad_currents,
        'time_s': elapsed,
        'nfev': getattr(result, 'nfev', None),
        'converged': result.success,
    }


# ── CSV I/O ──────────────────────────────────────────────────────────────────

def csv_header():
    cols = ['param_value', 'mse', 'alpha_x', 'alpha_y', 'beta_x', 'beta_y',
            'disp_resid']
    cols += [f'quad_{idx}' for idx in QUAD_INDICES]
    cols += ['time_s']
    return cols


def result_to_row(param_value, res):
    row = [param_value, res['mse'], res['alpha_x'], res['alpha_y'],
           res['beta_x'], res['beta_y'], res['disp_resid']]
    row += [res['quad_currents'][idx] for idx in QUAD_INDICES]
    row += [res['time_s']]
    return row


def write_csv(filepath, header, rows):
    with open(filepath, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def read_csv(filepath):
    """Read scan CSV back into list of dicts for plotting."""
    rows = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: float(v) for k, v in row.items()})
    return rows


# ── Plotting ─────────────────────────────────────────────────────────────────

MSE_THRESHOLDS = {
    'Excellent': 1e-3,
    'Acceptable': 0.01,
    'Marginal': 0.1,
}

TWISS_TARGETS = None  # set at runtime


def plot_mse_vs_param(param_vals, mse_vals, xlabel, title, filepath):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(param_vals, mse_vals, 'ko-', markersize=5)
    colors = {'Excellent': 'green', 'Acceptable': 'orange', 'Marginal': 'red'}
    for label, thresh in MSE_THRESHOLDS.items():
        ax.axhline(thresh, color=colors[label], linestyle='--', alpha=0.7,
                    label=f'{label} ({thresh})')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('MSE')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(filepath, format='eps')
    plt.close(fig)
    print(f"  Saved {filepath}")


def plot_twiss_vs_param(param_vals, rows, xlabel, title_base, filepath):
    beta_xm, alpha_xm, beta_ym, alpha_ym, _ = compute_twiss_targets()
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    panels = [
        ('beta_x', r'$\beta_x$ (m)', beta_xm),
        ('beta_y', r'$\beta_y$ (m)', beta_ym),
        ('alpha_x', r'$\alpha_x$', alpha_xm),
        ('alpha_y', r'$\alpha_y$', alpha_ym),
    ]
    for ax, (key, ylabel, target) in zip(axes.flat, panels):
        vals = [r[key] for r in rows]
        ax.plot(param_vals, vals, 'ko-', markersize=5)
        ax.axhline(target, color='blue', linestyle='--', alpha=0.5,
                    label=f'Target = {target}')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(True, alpha=0.3)

    fig.suptitle(title_base, fontsize=14)
    fig.tight_layout()
    fig.savefig(filepath, format='eps')
    plt.close(fig)
    print(f"  Saved {filepath}")


def plot_currents_vs_param(param_vals, rows, xlabel, title_base, filepath,
                           key_quads=None):
    """Plot selected quad currents vs parameter."""
    if key_quads is None:
        key_quads = [87, 93, 95, 97, 1, 3, 56, 58]
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx in key_quads:
        col = f'quad_{idx}'
        vals = [r[col] for r in rows]
        ax.plot(param_vals, vals, 'o-', markersize=4, label=f'elem {idx}')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Current (A)')
    ax.set_title(title_base)
    ax.legend(ncol=2, fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(filepath, format='eps')
    plt.close(fig)
    print(f"  Saved {filepath}")


# ── Scan runners ─────────────────────────────────────────────────────────────

def run_scan(scan_name, param_name, param_values, outdir, nb_particles=500,
             seed=42, **overrides):
    """Run a 1D parameter scan and save CSV."""
    print(f"\n{'='*70}")
    print(f"Scan: {scan_name} — sweeping {param_name}")
    print(f"{'='*70}")

    header = csv_header()
    rows = []
    for i, val in enumerate(param_values):
        kwargs = dict(BASELINE)
        kwargs['nb_particles'] = nb_particles
        kwargs['seed'] = seed
        kwargs.update(overrides)
        kwargs[param_name] = val
        print(f"  [{i+1}/{len(param_values)}] {param_name} = {val}")
        try:
            res = run_optimization(**kwargs)
            row = result_to_row(val, res)
            quality = 'FAILED'
            for label, thresh in sorted(MSE_THRESHOLDS.items(),
                                        key=lambda x: x[1]):
                if res['mse'] < thresh:
                    quality = label.upper()
                    break
            print(f"    MSE = {res['mse']:.4e}  [{quality}]  "
                  f"β_x={res['beta_x']:.4f}  β_y={res['beta_y']:.4f}  "
                  f"α_x={res['alpha_x']:.4f}  α_y={res['alpha_y']:.4f}  "
                  f"({res['time_s']:.1f} s)")
        except Exception as e:
            print(f"    FAILED: {e}")
            row = [val] + [float('nan')] * (len(header) - 1)
        rows.append(row)

    csv_path = outdir / f'scan_{scan_name}.csv'
    write_csv(csv_path, header, rows)
    print(f"  Saved {csv_path}")
    return rows


def generate_plots(outdir):
    """Generate all plots from saved CSV data."""
    print(f"\n{'='*70}")
    print("Generating plots")
    print(f"{'='*70}")

    scans = {
        'energy_spread': {
            'xlabel': r'Energy spread $\sigma_E$ (%)',
            'title_mse': r'MSE vs Energy Spread ($\sigma_E$)',
            'title_twiss': r'Twiss Parameters vs Energy Spread ($\sigma_E$)',
            'title_curr': r'Quad Currents vs Energy Spread ($\sigma_E$)',
        },
        'chirp': {
            'xlabel': r'Chirp $h$ ($10^9$ /s)',
            'title_mse': r'MSE vs Chirp ($h$)',
            'title_twiss': r'Twiss Parameters vs Chirp ($h$)',
            'title_curr': r'Quad Currents vs Chirp ($h$)',
        },
        'emittance': {
            'xlabel': r'Normalized emittance $\varepsilon_n$ ($\pi \cdot$mm$\cdot$mrad)',
            'title_mse': r'MSE vs Normalized Emittance ($\varepsilon_n$)',
            'title_twiss': r'Twiss Parameters vs Normalized Emittance ($\varepsilon_n$)',
            'title_curr': r'Quad Currents vs Normalized Emittance ($\varepsilon_n$)',
        },
    }

    for name, info in scans.items():
        csv_path = outdir / f'scan_{name}.csv'
        if not csv_path.exists():
            print(f"  Skipping {name}: {csv_path} not found")
            continue

        rows = read_csv(csv_path)
        param_vals = [r['param_value'] for r in rows]
        mse_vals = [r['mse'] for r in rows]

        # For chirp scan, display in units of 1e9
        if name == 'chirp':
            param_vals_display = [v / 1e9 for v in param_vals]
        else:
            param_vals_display = param_vals

        plot_mse_vs_param(
            param_vals_display, mse_vals, info['xlabel'],
            info['title_mse'], outdir / f'mse_vs_{name}.eps')

        plot_twiss_vs_param(
            param_vals_display, rows, info['xlabel'],
            info['title_twiss'], outdir / f'twiss_vs_{name}.eps')

        if name == 'energy_spread':
            plot_currents_vs_param(
                param_vals_display, rows, info['xlabel'],
                info['title_curr'],
                outdir / f'currents_vs_{name}.eps')


# ── Main ─────────────────────────────────────────────────────────────────────

def run_w1(outdir, nb_particles=1000):
    """W1: Compare Table I optimizations with and without chirp."""
    print("\n" + "=" * 70)
    print("W1: Table I Optimizations — Chirp vs No-Chirp Comparison")
    print("=" * 70)

    configs = [
        ("0.5 ps, h=5e9", dict(bunch_spread=0.5, h=5e9)),
        ("0.5 ps, h=0",   dict(bunch_spread=0.5, h=0)),
        ("2 ps, h=5e9",   dict(bunch_spread=2.0, h=5e9)),
        ("2 ps, h=0",     dict(bunch_spread=2.0, h=0)),
    ]

    results = []
    for label, kwargs in configs:
        print(f"\n  Running: {label}  ({nb_particles} particles)")
        res = run_optimization(nb_particles=nb_particles, seed=42, **kwargs)
        res['label'] = label
        results.append(res)
        print(f"    MSE = {res['mse']:.4e}  "
              f"β_x={res['beta_x']:.4f}  β_y={res['beta_y']:.4f}  "
              f"α_x={res['alpha_x']:.4f}  α_y={res['alpha_y']:.6f}  "
              f"({res['time_s']:.1f} s)")

    # Print comparison table
    beta_xm, alpha_xm, beta_ym, alpha_ym, _ = compute_twiss_targets()
    print(f"\n{'─' * 90}")
    print(f"{'Configuration':<22} {'MSE':>10} {'β_x':>8} {'β_y':>8} "
          f"{'α_x':>8} {'α_y':>10} {'Disp':>8}")
    print(f"{'Target':.<22} {'':>10} {beta_xm:>8.4f} {beta_ym:>8.4f} "
          f"{alpha_xm:>8.4f} {alpha_ym:>10.6f} {'0':>8}")
    print(f"{'─' * 90}")
    for r in results:
        print(f"{r['label']:<22} {r['mse']:>10.4e} {r['beta_x']:>8.4f} "
              f"{r['beta_y']:>8.4f} {r['alpha_x']:>8.4f} "
              f"{r['alpha_y']:>10.6f} {r['disp_resid']:>8.4f}")
    print(f"{'─' * 90}")

    # Save CSV
    header = ['label', 'bunch_spread', 'h'] + csv_header()[1:]
    rows = []
    for (label, kwargs), res in zip(configs, results):
        row = [label, kwargs['bunch_spread'], kwargs['h'],
               res['mse'], res['alpha_x'], res['alpha_y'],
               res['beta_x'], res['beta_y'], res['disp_resid']]
        row += [res['quad_currents'][idx] for idx in QUAD_INDICES]
        row += [res['time_s']]
        rows.append(row)

    csv_path = outdir / 'w1_chirp_comparison.csv'
    write_csv(csv_path, header, rows)
    print(f"\nSaved {csv_path}")

    # Check if chirp has effect
    for label_with, label_without in [
        ("0.5 ps, h=5e9", "0.5 ps, h=0"),
        ("2 ps, h=5e9", "2 ps, h=0"),
    ]:
        with_chirp = [r for r in results if r['label'] == label_with][0]
        without_chirp = [r for r in results if r['label'] == label_without][0]
        bunch = label_with.split(",")[0]
        mse_diff = abs(with_chirp['mse'] - without_chirp['mse'])
        bx_diff = abs(with_chirp['beta_x'] - without_chirp['beta_x'])
        by_diff = abs(with_chirp['beta_y'] - without_chirp['beta_y'])
        print(f"\n  {bunch}: ΔMSE = {mse_diff:.4e}, "
              f"Δβ_x = {bx_diff:.6f}, Δβ_y = {by_diff:.6f}")
        if mse_diff < 1e-3 and bx_diff < 0.01 and by_diff < 0.01:
            print(f"    → Chirp has negligible effect on Twiss matching")
        else:
            print(f"    → Chirp affects Twiss matching — investigate")


def run_w2(outdir, nb_particles=500, chrom_upper_bound=15, n_restarts=5):
    """W2: Enhanced emittance scan with fallback multi-start.

    Strategy per point:
      1. Run with default 10 A bounds, no multi-start (fast — same as original)
      2. If MSE > 1e-3: also try 10 A + multi-start, and raised chrom. bounds +
         multi-start. Keep the best result across all attempts.
    """
    RETRY_THRESHOLD = 1e-3

    print("\n" + "=" * 70)
    print(f"W2: Emittance Scan — Fallback Multi-Start ({n_restarts} restarts, "
          f"{chrom_upper_bound} A chrom. bound)")
    print("=" * 70)

    emittance_values = np.concatenate([
        [1, 2, 3, 4, 5],
        np.arange(6, 21, 1.0),
    ])

    header = csv_header()
    rows = []
    for i, en in enumerate(emittance_values):
        print(f"  [{i+1}/{len(emittance_values)}] ε_n = {en}")

        # Strategy 1: default (fast)
        best_res = None
        best_label = "default"
        try:
            res = run_optimization(
                epsilon_n=en, nb_particles=nb_particles, seed=42)
            best_res = res
        except Exception:
            pass

        # If not excellent, try fallback strategies
        if best_res is None or best_res['mse'] > RETRY_THRESHOLD:
            fallbacks = [
                ("10A+ms", dict(chrom_upper_bound=10, n_restarts=n_restarts)),
                (f"{int(chrom_upper_bound)}A+ms",
                 dict(chrom_upper_bound=chrom_upper_bound,
                      n_restarts=n_restarts)),
            ]
            for label, kwargs in fallbacks:
                try:
                    res = run_optimization(
                        epsilon_n=en, nb_particles=nb_particles, seed=42,
                        **kwargs)
                    if best_res is None or res['mse'] < best_res['mse']:
                        best_res = res
                        best_label = label
                except Exception:
                    pass

        if best_res is not None:
            row = result_to_row(en, best_res)
            quality = 'FAILED'
            for qlabel, thresh in sorted(MSE_THRESHOLDS.items(),
                                          key=lambda x: x[1]):
                if best_res['mse'] < thresh:
                    quality = qlabel.upper()
                    break
            print(f"    MSE = {best_res['mse']:.4e}  [{quality}]  "
                  f"β_x={best_res['beta_x']:.4f}  β_y={best_res['beta_y']:.4f}  "
                  f"q10={best_res['quad_currents'][10]:.2f}  "
                  f"strat={best_label}  ({best_res['time_s']:.1f} s)")
        else:
            print(f"    FAILED: all strategies failed")
            row = [en] + [float('nan')] * (len(header) - 1)
        rows.append(row)

    csv_path = outdir / 'scan_emittance_w2.csv'
    write_csv(csv_path, header, rows)
    print(f"\nSaved {csv_path}")

    # Generate comparison plot with original scan
    orig_path = outdir / 'scan_emittance.csv'
    if orig_path.exists():
        orig_rows = read_csv(orig_path)
        orig_en = [r['param_value'] for r in orig_rows]
        orig_mse = [r['mse'] for r in orig_rows]
        new_en = [r[0] for r in rows]
        new_mse = [r[1] for r in rows]

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.semilogy(orig_en, orig_mse, 's--', color='gray', markersize=5,
                     label='Original (10 A bound, 1 start)', alpha=0.6)
        ax.semilogy(new_en, new_mse, 'ko-', markersize=5,
                     label=f'W2 (dual-strategy, {n_restarts} starts)')
        for label, thresh in MSE_THRESHOLDS.items():
            colors = {'Excellent': 'green', 'Acceptable': 'orange',
                      'Marginal': 'red'}
            ax.axhline(thresh, color=colors[label], linestyle='--',
                        alpha=0.7, label=f'{label} ({thresh})')
        ax.set_xlabel(r'Normalized emittance $\varepsilon_n$ ($\pi \cdot$mm$\cdot$mrad)')
        ax.set_ylabel('MSE')
        ax.set_title('Emittance Scan — Original vs Multi-Start + Raised Bounds')
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(outdir / 'mse_vs_emittance_w2_comparison.eps', format='eps')
        plt.close(fig)
        print(f"  Saved {outdir / 'mse_vs_emittance_w2_comparison.eps'}")


def run_w6(outdir, nb_particles=500, chrom_upper_bound=15, n_restarts=5,
           glyfada_kwargs=None):
    """W6: Glyfada vs Nelder-Mead benchmark on emittance scan.

    Stages 1–10 always use Nelder-Mead. Only Stage 11 is substituted.
    """
    print("\n" + "=" * 70)
    print("W6: Glyfada vs Nelder-Mead Benchmark (Stage 11)")
    print("=" * 70)

    emittance_points = [5, 8, 14]
    g_kw = {'pop_size': 30, 'max_gen': 20, 'sigma': 0.05}
    if glyfada_kwargs:
        g_kw.update(glyfada_kwargs)

    w6dir = outdir / 'W6'
    w6dir.mkdir(parents=True, exist_ok=True)

    # CSV setup
    header = ['epsilon_n', 'method', 'mse', 'alpha_x', 'alpha_y',
              'beta_x', 'beta_y', 'disp_resid',
              'quad_87', 'quad_93', 'quad_95', 'quad_97',
              'time_s', 'nfev']
    rows = []

    for i, en in enumerate(emittance_points):
        for method_label, method_name, extra_kw in [
            ('NM', 'Nelder-Mead', dict(chrom_upper_bound=chrom_upper_bound,
                                        n_restarts=n_restarts)),
            ('glyfada', 'glyfada', dict(chrom_upper_bound=chrom_upper_bound,
                                         stage11_kwargs=g_kw)),
        ]:
            print(f"  [{i+1}/{len(emittance_points)}] ε_n={en}, method={method_label}")
            try:
                res = run_optimization(
                    epsilon_n=en, nb_particles=nb_particles, seed=42,
                    stage11_method=method_name, **extra_kw)
                row = [en, method_label, res['mse'],
                       res['alpha_x'], res['alpha_y'],
                       res['beta_x'], res['beta_y'], res['disp_resid'],
                       res['quad_currents'][87], res['quad_currents'][93],
                       res['quad_currents'][95], res['quad_currents'][97],
                       res['time_s'], res['nfev']]
                quality = 'FAILED'
                for qlabel, thresh in sorted(MSE_THRESHOLDS.items(),
                                              key=lambda x: x[1]):
                    if res['mse'] < thresh:
                        quality = qlabel.upper()
                        break
                print(f"    MSE = {res['mse']:.4e}  [{quality}]  "
                      f"β_x={res['beta_x']:.4f}  β_y={res['beta_y']:.4f}  "
                      f"({res['time_s']:.1f} s, nfev={res['nfev']})")
            except Exception as e:
                print(f"    FAILED: {e}")
                row = [en, method_label] + [float('nan')] * (len(header) - 2)
            rows.append(row)

    csv_path = w6dir / 'benchmark_results.csv'
    write_csv(csv_path, header, rows)
    print(f"\nSaved {csv_path}")

    # Parse results for plotting
    nm_rows = [r for r in rows if r[1] == 'NM']
    gl_rows = [r for r in rows if r[1] == 'glyfada']
    eps_nm = [r[0] for r in nm_rows]
    eps_gl = [r[0] for r in gl_rows]
    mse_nm = [r[2] for r in nm_rows]
    mse_gl = [r[2] for r in gl_rows]
    time_nm = [r[12] for r in nm_rows]
    time_gl = [r[12] for r in gl_rows]

    # MSE bar chart
    x = np.arange(len(emittance_points))
    width = 0.35
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x - width/2, mse_nm, width, label='Nelder-Mead (5 restarts)', color='#4477AA')
    ax.bar(x + width/2, mse_gl, width, label='Glyfada', color='#EE6677')
    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(e)) for e in emittance_points])
    ax.set_xlabel(r'Normalised emittance $\varepsilon_n$ ($\pi\cdot$mm$\cdot$mrad)')
    ax.set_ylabel('MSE')
    ax.set_title('W6: Stage 11 MSE — Nelder-Mead vs Glyfada')
    for thresh_label, thresh in MSE_THRESHOLDS.items():
        colors = {'Excellent': 'green', 'Acceptable': 'orange', 'Marginal': 'red'}
        ax.axhline(thresh, color=colors[thresh_label], linestyle='--', alpha=0.6,
                    label=f'{thresh_label} ({thresh})')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')
    fig.tight_layout()
    fig.savefig(w6dir / 'mse_comparison.eps', format='eps')
    plt.close(fig)
    print(f"  Saved {w6dir / 'mse_comparison.eps'}")

    # Time bar chart
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x - width/2, time_nm, width, label='Nelder-Mead (5 restarts)', color='#4477AA')
    ax.bar(x + width/2, time_gl, width, label='Glyfada', color='#EE6677')
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(e)) for e in emittance_points])
    ax.set_xlabel(r'Normalised emittance $\varepsilon_n$ ($\pi\cdot$mm$\cdot$mrad)')
    ax.set_ylabel('Wall-clock time (s)')
    ax.set_title('W6: Stage 11 Wall-Clock Time — Nelder-Mead vs Glyfada')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')
    fig.tight_layout()
    fig.savefig(w6dir / 'time_comparison.eps', format='eps')
    plt.close(fig)
    print(f"  Saved {w6dir / 'time_comparison.eps'}")

    # Summary table
    beta_xm, alpha_xm, beta_ym, alpha_ym, _ = compute_twiss_targets()
    print(f"\n{'─' * 95}")
    print(f"{'ε_n':>5} {'Method':>10} {'MSE':>12} {'β_x':>8} {'β_y':>8} "
          f"{'α_x':>8} {'α_y':>10} {'Time(s)':>8} {'nfev':>8}")
    print(f"{'':>5} {'Target':>10} {'':>12} {beta_xm:>8.4f} {beta_ym:>8.4f} "
          f"{alpha_xm:>8.4f} {alpha_ym:>10.6f}")
    print(f"{'─' * 95}")
    for r in rows:
        print(f"{r[0]:>5.0f} {r[1]:>10} {r[2]:>12.4e} {r[5]:>8.4f} {r[6]:>8.4f} "
              f"{r[3]:>8.4f} {r[4]:>10.6f} {r[12]:>8.1f} {r[13]:>8}")
    print(f"{'─' * 95}")


def main():
    parser = argparse.ArgumentParser(
        description='0.5 ps parameter sensitivity study')
    parser.add_argument('--smoke', action='store_true',
                        help='Run baseline point only and verify against v2')
    parser.add_argument('--plots-only', action='store_true',
                        help='Regenerate plots from existing CSV data')
    parser.add_argument('--particles', type=int, default=None,
                        help='Particles per run (default: 500 for scans, 1000 for W1)')
    parser.add_argument('--scan', choices=['energy_spread', 'chirp', 'emittance'],
                        help='Run only one scan')
    parser.add_argument('--w1', action='store_true',
                        help='W1: chirp vs no-chirp comparison (1000 particles)')
    parser.add_argument('--w2', action='store_true',
                        help='W2: emittance scan with multi-start + raised bounds')
    parser.add_argument('--w6', action='store_true',
                        help='W6: Glyfada vs NM benchmark on emittance scan')
    parser.add_argument('--chrom-bound', type=float, default=15,
                        help='Upper chromaticity quad bound for W2/W6 (default: 15 A)')
    parser.add_argument('--multi-start', type=int, default=5,
                        help='Number of Stage-11 restarts for W2/W6 (default: 5)')
    args = parser.parse_args()

    outdir = Path(__file__).resolve().parent / 'results' / 'params_05ps'
    outdir.mkdir(parents=True, exist_ok=True)

    if args.w1:
        run_w1(outdir, nb_particles=args.particles or 1000)
        return

    if args.w2:
        run_w2(outdir, nb_particles=args.particles or 500,
               chrom_upper_bound=args.chrom_bound,
               n_restarts=args.multi_start)
        return

    if args.w6:
        run_w6(outdir, nb_particles=args.particles or 500,
               chrom_upper_bound=args.chrom_bound,
               n_restarts=args.multi_start)
        return

    if args.plots_only:
        generate_plots(outdir)
        return

    if args.smoke:
        print("Smoke test: running baseline point (σ_E=0.5%, h=5e9, ε_n=8)")
        print(f"  Particles: 1000 (matching v2 study)")
        res = run_optimization(nb_particles=1000, seed=42)

        beta_xm, alpha_xm, beta_ym, alpha_ym, _ = compute_twiss_targets()
        print(f"\n  MSE      = {res['mse']:.4e}")
        print(f"  β_x      = {res['beta_x']:.4f}  (target {beta_xm})")
        print(f"  β_y      = {res['beta_y']:.4f}  (target {beta_ym:.4f})")
        print(f"  α_x      = {res['alpha_x']:.4f}  (target {alpha_xm})")
        print(f"  α_y      = {res['alpha_y']:.6f}  (target {alpha_ym})")
        print(f"  Disp(92) = {res['disp_resid']:.4f}")
        print(f"  Time     = {res['time_s']:.1f} s")

        # Verification: MSE should be < 1e-3 (excellent)
        ok = True
        if res['mse'] > 1e-3:
            print(f"\n  FAIL: MSE {res['mse']:.4e} > 1e-3")
            ok = False
        if abs(res['beta_x'] - beta_xm) > 0.05:
            print(f"\n  FAIL: β_x off by {abs(res['beta_x'] - beta_xm):.4f}")
            ok = False
        if abs(res['beta_y'] - beta_ym) > 0.05:
            print(f"\n  FAIL: β_y off by {abs(res['beta_y'] - beta_ym):.4f}")
            ok = False
        if abs(res['alpha_x'] - alpha_xm) > 0.05:
            print(f"\n  FAIL: α_x off by {abs(res['alpha_x'] - alpha_xm):.4f}")
            ok = False

        if ok:
            print("\n  PASS: Baseline matches v2 study within tolerance")
        else:
            print("\n  SMOKE TEST FAILED")
            sys.exit(1)
        return

    # ── Full sweeps ──────────────────────────────────────────────────────────
    nb = args.particles or 500
    t_total = time.perf_counter()

    # Scan A: energy_std_percent
    if args.scan is None or args.scan == 'energy_spread':
        energy_spread_values = np.concatenate([
            np.linspace(0.1, 1.0, 7),    # fine near baseline
            np.linspace(1.5, 5.0, 8),     # coarser at high σ_E
        ])
        run_scan('energy_spread', 'energy_std_percent', energy_spread_values,
                 outdir, nb_particles=nb)

    # Scan B: h (chirp)
    if args.scan is None or args.scan == 'chirp':
        chirp_values = np.linspace(0, 40e9, 12)
        run_scan('chirp', 'h', chirp_values, outdir, nb_particles=nb)

    # Scan C: epsilon_n
    if args.scan is None or args.scan == 'emittance':
        emittance_values = np.linspace(1, 20, 10)
        run_scan('emittance', 'epsilon_n', emittance_values, outdir,
                 nb_particles=nb)

    elapsed_total = time.perf_counter() - t_total
    print(f"\nTotal sweep time: {elapsed_total:.0f} s "
          f"({elapsed_total/60:.1f} min)")

    # Generate plots
    generate_plots(outdir)

    print("\nDone. Results in:", outdir)


if __name__ == "__main__":
    main()
