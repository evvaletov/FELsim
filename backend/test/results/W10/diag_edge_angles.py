#!/usr/bin/env python
"""Diagnostic: compare COSY σ_z with and without dipole edge angles."""
import sys, os, json, math, time, shutil
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
os.chdir(os.path.join(os.path.dirname(__file__), '..', '..'))

from pathlib import Path
from cosyAdapter import COSYAdapter
from cosyOptHelper import parse_beamline_felsim_indexed
from UHM_beamline_opt_05ps_params import run_optimization

EXCEL_PATH = '/home/evaletov/UH/fel-merge-workspace/FELsim/beam_excel/Beamline_elements.xlsx'
ENERGY = 40.0
SEGMENTS = 118
QUAD_HALF_APERTURE = 0.0135
DIPOLE_HALF_GAP = 0.007239
DIPOLE_HALF_WIDTH = 0.025
BETA_REL = 0.99992
C_LIGHT = 299792458.0
N_PARTICLES = 1000
RF_FREQ = 2.856e9
DIAG_DIR = Path(__file__).resolve().parent


def beam_sigma_z_ps(particles_cosy):
    return np.std(particles_cosy[:, 4]) / (BETA_REL * C_LIGHT) * 1e12


def apply_currents(beamline, currents):
    for idx, current in currents.items():
        if idx < len(beamline):
            beamline[idx]['current'] = current


def run_cosy(beamline, beam_cosy, label, test_dir):
    """Run COSY tracking with given beamline, saving all artifacts to test_dir."""
    os.makedirs(test_dir, exist_ok=True)

    sim = COSYAdapter(
        excel_path=EXCEL_PATH,
        mode='particle_tracking',
        fringe_field_order=0,
        quad_aperture=QUAD_HALF_APERTURE * 2,
        dipole_aperture=DIPOLE_HALF_GAP * 2,
        config={'simulation': {'dimensions': 3, 'KE': ENERGY}},
    )
    sim.enable_aperture_cuts(dipole_half_width=DIPOLE_HALF_WIDTH)
    native = sim.get_native_simulator()
    native.use_enge_coeffs = False
    native.beamline = beamline

    # Verify the beamline state
    dpw_angles = [(i, e.get('wedge_angle', 'N/A'))
                  for i, e in enumerate(native.beamline) if e['type'] == 'DPW']
    print(f"  DPW wedge angles in native.beamline: {dpw_angles[:4]}...")

    particles_felsim = native.transform_from_cosy_coordinates(beam_cosy)
    print(f"  Input: {beam_cosy.shape[0]} particles → {particles_felsim.shape[0]} after transform")

    # Run collect_evolution with separate output directory
    # Override the output_dir by calling native methods directly
    native.enable_particle_tracking(
        checkpoint_elements=list(range(1, len(native._detect_dipole_triplets()) + 1)))

    particles_cosy = native.transform_to_cosy_coordinates(particles_felsim)
    native.write_particle_file(particles_cosy, format='rray', output_dir=test_dir)

    result = native.run_simulation(output_dir=test_dir)
    if result.get('status') != 'success':
        print(f"  {label}: COSY FAILED: {result}")
        return None

    # Save input.fox for inspection
    shutil.copy(os.path.join(test_dir, 'input.fox'),
                DIAG_DIR / f'input_fox_{label.replace(" ", "_")}.fox')

    # Read final checkpoint
    grouped = native._detect_dipole_triplets()
    n_elem = len(grouped)
    checkpoints = native.read_checkpoints(
        list(range(1, n_elem + 1)),
        transform_to_felsim=True,
        validate=False,
        filter_invalid=True,
        output_dir=test_dir
    )

    if not checkpoints:
        print(f"  {label}: No checkpoints found")
        return None

    last_idx = max(checkpoints.keys())
    ps_final = checkpoints[last_idx]
    if ps_final is None or ps_final.shape[0] < 2:
        print(f"  {label}: Too few particles at final checkpoint")
        return None

    ps_cosy = native.transform_to_cosy_coordinates(ps_final)
    sz = beam_sigma_z_ps(ps_cosy)
    sy = np.std(ps_cosy[:, 2]) * 1e3
    sx = np.std(ps_cosy[:, 0]) * 1e3
    n_surv = ps_final.shape[0]
    print(f"  {label}: σ_z={sz:.2f} ps, σ_x={sx:.2f} mm, σ_y={sy:.2f} mm, "
          f"T={n_surv}/{beam_cosy.shape[0]} ({100*n_surv/beam_cosy.shape[0]:.1f}%)")
    return sz


# Get currents
cache_file = DIAG_DIR / 'diag_currents.json'
if cache_file.exists():
    print("Loading cached currents...")
    with open(cache_file) as f:
        currents = {int(k): v for k, v in json.load(f).items()}
else:
    print("Optimizing quad currents...")
    t0 = time.perf_counter()
    res = run_optimization(
        bunch_spread=2.0, energy_std_percent=0.5, h=0,
        epsilon_n=8, nb_particles=500, seed=42,
    )
    currents = res['quad_currents']
    print(f"  Done in {time.perf_counter()-t0:.1f}s, RMS={math.sqrt(res['mse']):.4e}")
    with open(cache_file, 'w') as f:
        json.dump(currents, f)

# Create beam in COSY coordinates (same method as W10)
from UHM_beamline_opt_05ps_params import compute_twiss_targets
_, _, _, _, relat = compute_twiss_targets()
norm = relat.gamma * relat.beta
epsilon_geom = 8 / norm  # pi.mm.mrad

sigma_x = 0.8e-3  # m
sigma_xp = epsilon_geom * 1e-6 / sigma_x  # rad
sigma_z = 2e-12 * BETA_REL * C_LIGHT  # m
sigma_delta = 0.005

rng = np.random.default_rng(42)
beam_cosy = np.zeros((N_PARTICLES, 6))
beam_cosy[:, 0] = rng.normal(0, sigma_x, N_PARTICLES)
beam_cosy[:, 1] = rng.normal(0, sigma_xp, N_PARTICLES)
beam_cosy[:, 2] = rng.normal(0, sigma_x, N_PARTICLES)
beam_cosy[:, 3] = rng.normal(0, sigma_xp, N_PARTICLES)
beam_cosy[:, 4] = rng.normal(0, sigma_z, N_PARTICLES)
beam_cosy[:, 5] = rng.normal(0, sigma_delta, N_PARTICLES)

print(f"Input beam: σ_x={np.std(beam_cosy[:,0])*1e3:.3f} mm, "
      f"σ_l={np.std(beam_cosy[:,4])*1e3:.4f} mm, "
      f"σ_δ={np.std(beam_cosy[:,5]):.5f}")

# Test 1: With wedge angles
print("\n=== Test 1: WITH dipole edge angles ===")
bl = parse_beamline_felsim_indexed(EXCEL_PATH)
apply_currents(bl, currents)
sz1 = run_cosy(bl[:SEGMENTS], beam_cosy, "with_edges",
               str(DIAG_DIR / 'test1_with_edges'))

# Test 2: Without wedge angles
print("\n=== Test 2: WITHOUT dipole edge angles ===")
bl = parse_beamline_felsim_indexed(EXCEL_PATH)
apply_currents(bl, currents)
for e in bl:
    if e['type'] == 'DPW':
        e['wedge_angle'] = 0.0
sz2 = run_cosy(bl[:SEGMENTS], beam_cosy, "no_edges",
               str(DIAG_DIR / 'test2_no_edges'))

print(f"\n{'='*60}")
print(f"  With edges:    σ_z = {sz1:.2f} ps" if sz1 else "  With edges: FAILED")
print(f"  Without edges: σ_z = {sz2:.2f} ps" if sz2 else "  Without edges: FAILED")
print(f"{'='*60}")
