"""
Verify achievability of 2 ps bunch Twiss parameters at MkIII undulator entrance.

Target parameters at undulator entrance:
    Energy: 40 MeV (γ ≈ 78.3)
    Energy spread: 0.5% rms
    Normalised emittance: εn,x = εn,y = 8 π·mm·mrad
    Bunch duration: 2 ps FWHM
    Peak current: 30 A
    βx = 1.4 m, αx = 0.47
    βy = 0.24 m, αy = 0

Author: Eremey Valetov
"""
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from pathlib import Path

from ebeam import beam
from beamline import lattice
from excelElements import ExcelElements
from beamOptimizer import beamOptimizer
from schematic import draw_beamline


# =============================================================================
# Physical parameters
# =============================================================================

ENERGY_MEV = 40.0
EPSILON_N = 8.0          # π·mm·mrad (normalised)
ENERGY_SPREAD_RMS = 0.5  # percent
BUNCH_DURATION_FWHM = 2  # ps
RF_FREQ = 2856e6         # Hz

# Target Twiss at undulator entrance
TARGET_BETA_X = 1.4      # m
TARGET_ALPHA_X = 0.47
TARGET_BETA_Y = 0.24     # m
TARGET_ALPHA_Y = 0.0

# Beamline file
EXCEL_PATH = Path(__file__).parent.parent.parent / "beam_excel" / "Beamline_elements_3.xlsx"


def compute_relativistic_params(energy_mev):
    """Compute γ, β, and normalisation factor."""
    relat = lattice(1, fringeType=None)
    relat.setE(E=energy_mev)
    return relat.gamma, relat.beta, relat.gamma * relat.beta


def generate_beam_from_twiss(alpha_x, beta_x, alpha_y, beta_y,
                              epsilon_geom, energy_spread_rms,
                              bunch_duration_ps, n_particles=1000):
    """
    Generate 6D particle distribution from specified Twiss parameters.

    Uses the covariance matrix approach rather than uncorrelated Gaussians,
    so that alpha (x-x' correlation) is properly represented.
    """
    ebeam_obj = beam()

    # Longitudinal parameters
    # FWHM = 2.355 * σ for Gaussian
    bunch_sigma_ps = bunch_duration_ps / 2.355
    tof_std = bunch_sigma_ps * 1e-12 * RF_FREQ  # in units of RF period fraction
    energy_std = energy_spread_rms * 10  # convert to 10^-3 ΔW/W

    # For longitudinal, assume uncorrelated (α_z = 0)
    # ε_z = σ_t * σ_δ (geometric)
    epsilon_z = tof_std * energy_std
    beta_z = tof_std**2 / epsilon_z if epsilon_z > 0 else 1.0

    twiss_params = {
        'x': {'alpha': alpha_x, 'beta': beta_x, 'epsilon': epsilon_geom, 'phi': 0.0},
        'y': {'alpha': alpha_y, 'beta': beta_y, 'epsilon': epsilon_geom, 'phi': 0.0},
        'z': {'alpha': 0.0, 'beta': beta_z, 'epsilon': epsilon_z, 'phi': 0.0},
    }

    return ebeam_obj.gen_6d_from_twiss(twiss_params, n_particles)


def load_beamline(excel_path, energy_mev):
    """Load beamline from Excel and configure for electron beam."""
    excel = ExcelElements(excel_path)
    beamline_raw = excel.create_beamline()

    relat = lattice(1, fringeType=None)
    return relat.changeBeamType("electron", energy_mev, beamline_raw)


def find_undulator_entrance_index(beamline_elements):
    """Find the element index corresponding to undulator entrance."""
    for i, elem in enumerate(beamline_elements):
        # The undulator is around z ≈ 12.39 m
        if hasattr(elem, 'z_end') and elem.z_end is not None:
            if elem.z_end > 12.3:
                return i
    # Fallback: use element 117 as in original script
    return min(117, len(beamline_elements) - 1)


def identify_all_quads(beamline_elements):
    """
    Identify all quadrupole indices in the beamline.
    Returns list of (index, element_type, current) tuples.
    """
    quads = []
    for i, elem in enumerate(beamline_elements):
        elem_type = type(elem).__name__
        if 'qp' in elem_type.lower() or 'quad' in elem_type.lower():
            current = getattr(elem, 'current', 0.0)
            quads.append((i, elem_type, current))
    return quads


def identify_final_matching_quads(beamline_elements, num_quads=4):
    """
    Identify quadrupole indices available for final undulator matching.
    Returns list of element indices.

    Parameters
    ----------
    beamline_elements : list
        Beamline element list
    num_quads : int
        Number of final quadrupoles to include (default 4)
    """
    all_quads = identify_all_quads(beamline_elements)
    quad_indices = [q[0] for q in all_quads]

    if len(quad_indices) >= num_quads:
        return quad_indices[-num_quads:]
    return quad_indices


def run_optimisation(beamline, beam_dist, quad_indices, target_index,
                     target_twiss, print_results=True, plot_progress=False,
                     method="Nelder-Mead", weights=None):
    """
    Run quadrupole optimisation to achieve target Twiss at undulator.

    Parameters
    ----------
    beamline : list
        Beamline element list
    beam_dist : ndarray
        Initial particle distribution
    quad_indices : list
        Indices of quadrupoles to optimise
    target_index : int
        Element index where target Twiss should be achieved
    target_twiss : dict
        Target Twiss parameters {'beta_x', 'alpha_x', 'beta_y', 'alpha_y'}
    method : str
        Optimisation method (Nelder-Mead, Powell, COBYLA, etc.)
    weights : dict, optional
        Weights for objectives {'alpha_x', 'alpha_y', 'beta_x', 'beta_y'}
    """
    if weights is None:
        weights = {'alpha_x': 1.0, 'alpha_y': 1.0, 'beta_x': 1.0, 'beta_y': 1.0}

    opti = beamOptimizer(beamline, beam_dist)

    # Set up optimisation variables (quadrupole currents)
    variables = {}
    start_point = {}
    var_names = [f"I{i}" for i in range(len(quad_indices))]

    for i, (var_name, quad_idx) in enumerate(zip(var_names, quad_indices)):
        variables[quad_idx] = [var_name, "current", lambda num: num]
        # Get current value as starting point
        current_val = getattr(beamline[quad_idx], 'current', 2.0)
        start_point[var_name] = {"bounds": (0.01, 10), "start": max(0.1, current_val)}

    # Set up objectives with weights
    # Normalise weights by target values to make objectives comparable
    w_alpha_x = weights['alpha_x']
    w_alpha_y = weights['alpha_y']
    w_beta_x = weights['beta_x'] / target_twiss['beta_x']**2 if target_twiss['beta_x'] != 0 else weights['beta_x']
    w_beta_y = weights['beta_y'] / target_twiss['beta_y']**2 if target_twiss['beta_y'] != 0 else weights['beta_y']

    objectives = {
        target_index: [
            {"measure": ["x", "alpha"], "goal": target_twiss['alpha_x'], "weight": w_alpha_x},
            {"measure": ["y", "alpha"], "goal": target_twiss['alpha_y'], "weight": w_alpha_y},
            {"measure": ["x", "beta"], "goal": target_twiss['beta_x'], "weight": w_beta_x},
            {"measure": ["y", "beta"], "goal": target_twiss['beta_y'], "weight": w_beta_y},
        ]
    }

    result = opti.calc(
        method,
        variables,
        start_point,
        objectives,
        plotBeam=False,
        printResults=print_results,
        plotProgress=plot_progress
    )

    return result, opti


# =============================================================================
# Sequential optimization (following UHM_beamline_opt.py structure)
# =============================================================================

# Define optimization steps: each step is a dict with:
#   'name': description
#   'quads': list of quad element indices
#   'targets': dict of target_index -> list of objectives
#   'start': starting current values (optional)
OPTIMIZATION_STEPS = [
    {
        'name': 'First Quadrupole Doublet',
        'quads': [1, 3],
        'targets': {
            8: [{'measure': ['x', 'alpha'], 'goal': 0, 'weight': 1}],
            9: [{'measure': ['y', 'alpha'], 'goal': 0, 'weight': 1},
                {'measure': ['y', 'beta'], 'goal': 0.1, 'weight': 0.5}]
        },
        'start': [1.0, 1.0]
    },
    {
        'name': 'First Chromaticity Quad',
        'quads': [10],
        'targets': {15: [{'measure': ['x', 'dispersion'], 'goal': 0, 'weight': 1}]},
        'start': [1.0]
    },
    {
        'name': 'Quadrupole Triplet (DC2)',
        'quads': [16, 18, 20],
        'targets': {
            25: [{'measure': ['x', 'alpha'], 'goal': 0, 'weight': 1},
                 {'measure': ['x', 'beta'], 'goal': 0.1, 'weight': 0.5}],
            26: [{'measure': ['y', 'alpha'], 'goal': 0, 'weight': 1},
                 {'measure': ['y', 'beta'], 'goal': 0.1, 'weight': 0.5}]
        },
        'start': [2.0, 5.0, 3.0]
    },
    {
        'name': 'Second Chromaticity Quad',
        'quads': [27],
        'targets': {32: [{'measure': ['x', 'dispersion'], 'goal': 0, 'weight': 1}]},
        'start': [1.0]
    },
    {
        'name': 'Double Quadrupole Triplet',
        'quads': [33, 35, 37],
        'targets': {
            37: [{'measure': ['x', 'alpha'], 'goal': 0, 'weight': 1},
                 {'measure': ['y', 'alpha'], 'goal': 0, 'weight': 1},
                 {'measure': ['x', 'envelope'], 'goal': 2.0, 'weight': 1},
                 {'measure': ['y', 'envelope'], 'goal': 2.0, 'weight': 1}]
        },
        'start': [2.0, 2.0, 2.0],
        'copy_to': [(37, 39), (35, 41), (33, 43)]  # Copy currents to symmetric quads
    },
    {
        'name': 'Third Chromaticity Quad',
        'quads': [50],
        'targets': {55: [{'measure': ['x', 'dispersion'], 'goal': 0, 'weight': 1}]},
        'start': [1.0]
    },
    {
        'name': 'IP Quadrupole Doublet',
        'quads': [56, 58],
        'targets': {
            59: [{'measure': ['x', 'envelope'], 'goal': 0.0, 'weight': 1},
                 {'measure': ['y', 'envelope'], 'goal': 0.0, 'weight': 1}]
        },
        'start': [2.0, 2.0]
    },
    {
        'name': 'Quadrupole Doublet (DC4)',
        'quads': [61, 63],
        'targets': {
            68: [{'measure': ['x', 'alpha'], 'goal': 0, 'weight': 1},
                 {'measure': ['x', 'beta'], 'goal': 0.1, 'weight': 0.5}],
            69: [{'measure': ['y', 'alpha'], 'goal': 0, 'weight': 1},
                 {'measure': ['y', 'beta'], 'goal': 0.1, 'weight': 0.5}]
        },
        'start': [2.0, 2.0]
    },
    {
        'name': 'Fourth Chromaticity Quad',
        'quads': [70],
        'targets': {75: [{'measure': ['x', 'dispersion'], 'goal': 0, 'weight': 1}]},
        'start': [1.0]
    },
    {
        'name': 'Quadrupole Doublet (DC5)',
        'quads': [76, 78],
        'targets': {
            85: [{'measure': ['x', 'alpha'], 'goal': 0, 'weight': 1},
                 {'measure': ['x', 'beta'], 'goal': 0.1, 'weight': 0.5}],
            86: [{'measure': ['y', 'alpha'], 'goal': 0, 'weight': 1},
                 {'measure': ['y', 'beta'], 'goal': 0.1, 'weight': 0.5}]
        },
        'start': [2.0, 2.0]
    },
    # Skip 5th chromaticity quad - include it in final matching for more DOF
    # Final step: include quads 80, 87, 93, 95, 97 for 5 DOF matching 4 constraints
    {
        'name': 'Final Undulator Matching (5 quads)',
        'quads': [80, 87, 93, 95, 97],
        'targets': None,  # Will be set dynamically based on user target
        'start': [2.0, 3.0, 2.0, 2.0, 2.0]
    },
]


def run_sequential_optimisation(beamline, beam_dist, target_twiss,
                                 undulator_idx=117, print_results=True,
                                 method='Nelder-Mead'):
    """
    Run sequential quadrupole optimisation through the beamline.

    This approach optimises each section separately, building up
    the optics progressively from linac exit to undulator entrance.

    Parameters
    ----------
    beamline : list
        Full beamline element list
    beam_dist : ndarray
        Initial particle distribution
    target_twiss : dict
        Target Twiss at undulator {'beta_x', 'alpha_x', 'beta_y', 'alpha_y'}
    undulator_idx : int
        Element index of undulator entrance
    print_results : bool
        Print results after each step
    method : str
        Optimisation method

    Returns
    -------
    dict
        Results including final quad currents and achieved Twiss
    """
    print("\n" + "=" * 60)
    print("Sequential Optimisation Mode")
    print("=" * 60)

    opti = beamOptimizer(beamline, beam_dist)
    results = []

    for step_num, step in enumerate(OPTIMIZATION_STEPS):
        print(f"\n--- Step {step_num + 1}: {step['name']} ---")

        # Set up variables
        variables = {}
        start_point = {}
        for i, quad_idx in enumerate(step['quads']):
            var_name = f"I{i}"
            variables[quad_idx] = [var_name, "current", lambda num: num]
            start_val = step['start'][i] if step.get('start') else 2.0
            start_point[var_name] = {"bounds": (0.01, 10), "start": start_val}

        # Set up objectives
        if step['targets'] is None:
            # Final step: use user-specified target Twiss
            objectives = {
                undulator_idx: [
                    {'measure': ['x', 'alpha'], 'goal': target_twiss['alpha_x'], 'weight': 1},
                    {'measure': ['y', 'alpha'], 'goal': target_twiss['alpha_y'], 'weight': 1},
                    {'measure': ['x', 'beta'], 'goal': target_twiss['beta_x'], 'weight': 1},
                    {'measure': ['y', 'beta'], 'goal': target_twiss['beta_y'], 'weight': 1},
                ]
            }
        else:
            objectives = step['targets']

        # Run optimisation for this step
        try:
            result = opti.calc(
                method,
                variables,
                start_point,
                objectives,
                plotBeam=False,
                printResults=print_results,
                plotProgress=False
            )
            results.append({'step': step['name'], 'success': True, 'result': result})

            # Copy currents to symmetric quads if specified
            if 'copy_to' in step:
                for src_idx, dst_idx in step['copy_to']:
                    beamline[dst_idx].current = beamline[src_idx].current
                    if print_results:
                        print(f"  Copied current from quad {src_idx} to {dst_idx}")

        except Exception as e:
            print(f"  Step failed: {e}")
            results.append({'step': step['name'], 'success': False, 'error': str(e)})

    # Collect final quad currents
    final_currents = {}
    for step in OPTIMIZATION_STEPS:
        for quad_idx in step['quads']:
            final_currents[quad_idx] = getattr(beamline[quad_idx], 'current', None)

    return {
        'results': results,
        'final_currents': final_currents,
        'beamline': beamline,
    }


def evaluate_final_twiss(beamline, beam_dist, target_index):
    """Propagate beam and compute Twiss at target location."""
    particles = beam_dist.copy()
    for i, elem in enumerate(beamline[:target_index + 1]):
        particles = np.array(elem.useMatrice(particles))

    ebeam_obj = beam()
    _, _, twiss_df = ebeam_obj.cal_twiss(particles)

    return {
        'beta_x': twiss_df.loc['x'][r"$\beta$ (m)"],
        'alpha_x': twiss_df.loc['x'][r"$\alpha$"],
        'beta_y': twiss_df.loc['y'][r"$\beta$ (m)"],
        'alpha_y': twiss_df.loc['y'][r"$\alpha$"],
        'epsilon_x': twiss_df.loc['x'][r"$\epsilon$ ($\pi$.mm.mrad)"],
        'epsilon_y': twiss_df.loc['y'][r"$\epsilon$ ($\pi$.mm.mrad)"],
    }


def evaluate_baseline(beamline, beam_dist, target_index):
    """
    Evaluate Twiss at target with current (unoptimised) quadrupole settings.
    Prints diagnostics about what the beamline currently produces.
    """
    print("\n" + "-" * 60)
    print("Baseline (unoptimised) Twiss at undulator entrance:")
    print("-" * 60)

    twiss = evaluate_final_twiss(beamline, beam_dist, target_index)
    print(f"  βx = {twiss['beta_x']:.4f} m")
    print(f"  αx = {twiss['alpha_x']:.4f}")
    print(f"  βy = {twiss['beta_y']:.4f} m")
    print(f"  αy = {twiss['alpha_y']:.4f}")
    print(f"  εx = {twiss['epsilon_x']:.4f} π·mm·mrad")
    print(f"  εy = {twiss['epsilon_y']:.4f} π·mm·mrad")

    return twiss


def print_comparison(target, achieved):
    """Print target vs achieved Twiss comparison."""
    print("\n" + "=" * 60)
    print("Twiss Parameter Comparison at Undulator Entrance")
    print("=" * 60)
    print(f"{'Parameter':<15} {'Target':>12} {'Achieved':>12} {'Diff':>12}")
    print("-" * 60)

    params = [
        ('βx [m]', 'beta_x'),
        ('αx', 'alpha_x'),
        ('βy [m]', 'beta_y'),
        ('αy', 'alpha_y'),
    ]

    for label, key in params:
        tgt = target[key]
        ach = achieved[key]
        diff = ach - tgt
        pct = 100 * diff / tgt if tgt != 0 else float('inf')
        print(f"{label:<15} {tgt:>12.4f} {ach:>12.4f} {diff:>+12.4f} ({pct:+.1f}%)")

    print("-" * 60)
    print(f"{'εx [π·mm·mrad]':<15} {'-':>12} {achieved['epsilon_x']:>12.4f}")
    print(f"{'εy [π·mm·mrad]':<15} {'-':>12} {achieved['epsilon_y']:>12.4f}")
    print("=" * 60)


def main(initial_beta_x=6.3, initial_alpha_x=0.0,
         initial_beta_y=6.3, initial_alpha_y=0.0,
         n_particles=2000, plot_result=True, num_quads=8,
         target_beta_x=None, target_beta_y=None,
         target_alpha_x=None, target_alpha_y=None,
         opt_method='Nelder-Mead', sequential=False):
    """
    Main verification routine.

    Parameters
    ----------
    initial_beta_x, initial_alpha_x : float
        Initial Twiss parameters in x-plane at linac exit
    initial_beta_y, initial_alpha_y : float
        Initial Twiss parameters in y-plane at linac exit
    n_particles : int
        Number of particles for simulation
    plot_result : bool
        Whether to plot the final beam evolution
    num_quads : int
        Number of final quadrupoles to include in optimization
    target_beta_x, target_alpha_x : float, optional
        Override target Twiss in x-plane
    target_beta_y, target_alpha_y : float, optional
        Override target Twiss in y-plane
    """
    # Use defaults from module constants if not overridden
    tgt_beta_x = target_beta_x if target_beta_x is not None else TARGET_BETA_X
    tgt_alpha_x = target_alpha_x if target_alpha_x is not None else TARGET_ALPHA_X
    tgt_beta_y = target_beta_y if target_beta_y is not None else TARGET_BETA_Y
    tgt_alpha_y = target_alpha_y if target_alpha_y is not None else TARGET_ALPHA_Y
    print("=" * 60)
    print("2 ps Bunch Undulator Matching Verification")
    print("=" * 60)

    # Compute relativistic parameters
    gamma, beta_rel, norm_factor = compute_relativistic_params(ENERGY_MEV)
    epsilon_geom = EPSILON_N / norm_factor

    print(f"\nBeam parameters:")
    print(f"  Energy: {ENERGY_MEV} MeV")
    print(f"  γ = {gamma:.2f}, β = {beta_rel:.6f}")
    print(f"  εn = {EPSILON_N} π·mm·mrad")
    print(f"  ε_geom = {epsilon_geom:.4f} π·mm·mrad")
    print(f"  Energy spread: {ENERGY_SPREAD_RMS}% rms")
    print(f"  Bunch duration: {BUNCH_DURATION_FWHM} ps FWHM")

    print(f"\nInitial Twiss (at linac exit):")
    print(f"  βx = {initial_beta_x:.2f} m, αx = {initial_alpha_x:.2f}")
    print(f"  βy = {initial_beta_y:.2f} m, αy = {initial_alpha_y:.2f}")

    print(f"\nTarget Twiss (at undulator entrance):")
    print(f"  βx = {tgt_beta_x} m, αx = {tgt_alpha_x}")
    print(f"  βy = {tgt_beta_y} m, αy = {tgt_alpha_y}")

    # Generate initial beam distribution
    print(f"\nGenerating {n_particles} particles...")
    beam_dist = generate_beam_from_twiss(
        initial_alpha_x, initial_beta_x,
        initial_alpha_y, initial_beta_y,
        epsilon_geom, ENERGY_SPREAD_RMS,
        BUNCH_DURATION_FWHM, n_particles
    )

    # Verify initial Twiss
    ebeam_obj = beam()
    _, _, initial_twiss = ebeam_obj.cal_twiss(beam_dist)
    print(f"\nVerified initial Twiss from generated distribution:")
    print(f"  βx = {initial_twiss.loc['x'][r'$\beta$ (m)']:.4f} m")
    print(f"  αx = {initial_twiss.loc['x'][r'$\alpha$']:.4f}")
    print(f"  βy = {initial_twiss.loc['y'][r'$\beta$ (m)']:.4f} m")
    print(f"  αy = {initial_twiss.loc['y'][r'$\alpha$']:.4f}")

    # Load beamline
    print(f"\nLoading beamline from {EXCEL_PATH}...")
    if not EXCEL_PATH.exists():
        print(f"ERROR: Excel file not found: {EXCEL_PATH}")
        return None

    beamline_elements = load_beamline(EXCEL_PATH, ENERGY_MEV)
    print(f"  Loaded {len(beamline_elements)} elements")

    # Find undulator entrance
    undulator_idx = find_undulator_entrance_index(beamline_elements)
    print(f"  Undulator entrance at element index {undulator_idx}")

    # Truncate beamline to undulator entrance
    beamline_truncated = beamline_elements[:undulator_idx + 1]

    # Identify matching quadrupoles
    quad_indices = identify_final_matching_quads(beamline_truncated, num_quads)
    print(f"  Using {len(quad_indices)} quadrupoles for matching: {quad_indices}")

    # Define target Twiss
    target_twiss = {
        'beta_x': tgt_beta_x,
        'alpha_x': tgt_alpha_x,
        'beta_y': tgt_beta_y,
        'alpha_y': tgt_alpha_y,
    }

    # Evaluate baseline (current beamline settings)
    baseline_twiss = evaluate_baseline(beamline_truncated, beam_dist, undulator_idx)

    # Run optimisation
    if sequential:
        # Sequential mode: optimise section by section
        seq_result = run_sequential_optimisation(
            beamline_truncated, beam_dist, target_twiss,
            undulator_idx=undulator_idx,
            print_results=True,
            method=opt_method
        )
        result = seq_result
    else:
        # Simultaneous mode: optimise all selected quads at once
        print("\n" + "-" * 60)
        print("Running quadrupole optimisation (simultaneous)...")
        print("-" * 60)

        result, opti = run_optimisation(
            beamline_truncated, beam_dist, quad_indices, undulator_idx,
            target_twiss, print_results=True, plot_progress=False,
            method=opt_method
        )

    # Evaluate final Twiss
    achieved_twiss = evaluate_final_twiss(beamline_truncated, beam_dist, undulator_idx)

    # Print comparison
    print_comparison(target_twiss, achieved_twiss)

    # Assess feasibility
    beta_x_err = abs(achieved_twiss['beta_x'] - tgt_beta_x) / tgt_beta_x if tgt_beta_x != 0 else 0
    beta_y_err = abs(achieved_twiss['beta_y'] - tgt_beta_y) / tgt_beta_y if tgt_beta_y != 0 else 0
    alpha_x_err = abs(achieved_twiss['alpha_x'] - tgt_alpha_x)
    alpha_y_err = abs(achieved_twiss['alpha_y'] - tgt_alpha_y)

    feasible = (beta_x_err < 0.1 and beta_y_err < 0.1 and
                alpha_x_err < 0.1 and alpha_y_err < 0.1)

    print(f"\nFeasibility assessment: {'ACHIEVABLE' if feasible else 'DIFFICULT'}")
    if not feasible:
        print("  Consider:")
        print("  - Different initial Twiss at linac exit")
        print("  - Additional quadrupoles in optimisation")
        print("  - Relaxing target tolerances")

    # Plot if requested
    if plot_result:
        print("\nPlotting beam evolution (close window to continue)...")
        schem = draw_beamline()
        acceptance = {"shape": 'circle', "radius": 10.0, "origin": [0, 0]}
        schem.plotBeamPositionTransform(
            beam_dist, beamline_truncated, 0.1,
            plot=True, showIndice=False, defineLim=False,
            shape=acceptance, matchScaling=False, scatter=True
        )

    return {
        'result': result,
        'target': target_twiss,
        'achieved': achieved_twiss,
        'feasible': feasible,
        'quad_indices': quad_indices,
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Verify 2 ps bunch Twiss parameters at MkIII undulator"
    )
    parser.add_argument('--beta-x0', type=float, default=6.3,
                        help='Initial βx at linac exit [m]')
    parser.add_argument('--alpha-x0', type=float, default=0.0,
                        help='Initial αx at linac exit')
    parser.add_argument('--beta-y0', type=float, default=6.3,
                        help='Initial βy at linac exit [m]')
    parser.add_argument('--alpha-y0', type=float, default=0.0,
                        help='Initial αy at linac exit')
    parser.add_argument('--particles', type=int, default=2000,
                        help='Number of particles')
    parser.add_argument('--num-quads', type=int, default=8,
                        help='Number of quadrupoles to optimise')
    parser.add_argument('--target-beta-x', type=float, default=None,
                        help='Target βx [m] (default: 1.4)')
    parser.add_argument('--target-beta-y', type=float, default=None,
                        help='Target βy [m] (default: 0.24)')
    parser.add_argument('--target-alpha-x', type=float, default=None,
                        help='Target αx (default: 0.47)')
    parser.add_argument('--target-alpha-y', type=float, default=None,
                        help='Target αy (default: 0.0)')
    parser.add_argument('--symmetric', action='store_true',
                        help='Use symmetric target (βx=βy=0.24m, α=0)')
    parser.add_argument('--method', type=str, default='Nelder-Mead',
                        choices=['Nelder-Mead', 'Powell', 'COBYLA', 'L-BFGS-B'],
                        help='Optimisation method')
    parser.add_argument('--sequential', action='store_true',
                        help='Use sequential optimisation (section by section)')
    parser.add_argument('--no-plot', action='store_true',
                        help='Disable plotting')

    args = parser.parse_args()

    # Handle target overrides
    if args.symmetric:
        # Use symmetric undulator matching (original script target)
        args.target_beta_x = 0.24
        args.target_beta_y = 0.24
        args.target_alpha_x = 0.0
        args.target_alpha_y = 0.0

    main(
        initial_beta_x=args.beta_x0,
        initial_alpha_x=args.alpha_x0,
        initial_beta_y=args.beta_y0,
        initial_alpha_y=args.alpha_y0,
        n_particles=args.particles,
        plot_result=not args.no_plot,
        num_quads=args.num_quads,
        target_beta_x=args.target_beta_x,
        target_beta_y=args.target_beta_y,
        target_alpha_x=args.target_alpha_x,
        target_alpha_y=args.target_alpha_y,
        opt_method=args.method,
        sequential=args.sequential
    )
