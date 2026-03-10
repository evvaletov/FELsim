"""
Basic tests for RF-Track adapter.

Author: Eremey Valetov
"""
import numpy as np
import pytest
from pathlib import Path

# Check RF-Track availability
try:
    from rftrackAdapter import RFTrackAdapter, _RFTRACK_AVAILABLE
    if not _RFTRACK_AVAILABLE:
        pytest.skip("RF-Track not installed (pip install RF-Track)", allow_module_level=True)
except ImportError as e:
    pytest.skip(f"Import error: {e}", allow_module_level=True)

from simulatorBase import BeamlineElement, CoordinateSystem, SimulationMode
from simulatorFactory import SimulatorFactory
from physicalConstants import PhysicalConstants

BEAM_ENERGY = 45.0


def test_adapter_creation():
    """Test RFTrackAdapter instantiation."""
    print("\n=== Adapter Creation ===")

    sim = RFTrackAdapter(beam_energy=BEAM_ENERGY)
    print(f"Name: {sim.name}")
    print(f"Native coordinates: {sim.native_coordinates.value}")
    print(f"Mode: {sim.simulation_mode}")
    print(f"G_quad: {sim.G_quad:.4f} T/A/m")

    assert sim.name == "RF-Track"
    assert sim.native_coordinates == CoordinateSystem.RFTRACK
    assert sim.simulation_mode == SimulationMode.PARTICLE_TRACKING
    assert sim.G_quad == PhysicalConstants.G_quad_default

    sim_custom = RFTrackAdapter(beam_energy=30.0, G_quad=3.0)
    assert sim_custom.G_quad == 3.0
    assert sim_custom.beam_energy == 30.0

    print("Adapter creation: PASSED")


def test_factory_registration():
    """Test SimulatorFactory integration."""
    print("\n=== Factory Registration ===")

    available = SimulatorFactory.get_available_simulators()
    print(f"Available simulators: {available}")
    assert 'rftrack' in available

    sim = SimulatorFactory.create('rftrack', beam_energy=BEAM_ENERGY)
    assert sim.name == "RF-Track"

    info = SimulatorFactory.get_simulator_info('rftrack')
    print(f"Simulator info: {info}")
    assert info['class'] == 'RFTrackAdapter'

    print("Factory registration: PASSED")


def test_relativistic_params():
    """Test relativistic parameter calculations."""
    print("\n=== Relativistic Parameters ===")

    sim = RFTrackAdapter(beam_energy=BEAM_ENERGY)

    gamma_expected, beta_expected = PhysicalConstants.relativistic_parameters(
        BEAM_ENERGY, PhysicalConstants.E0_electron
    )

    print(f"γ: {sim._gamma:.4f} (expected: {gamma_expected:.4f})")
    print(f"β: {sim._beta:.6f} (expected: {beta_expected:.6f})")
    print(f"Pc: {sim._Pc:.3f} MeV/c")

    assert abs(sim._gamma - gamma_expected) < 1e-6
    assert abs(sim._beta - beta_expected) < 1e-6

    sim.set_beam_energy(100.0)
    gamma_100, _ = PhysicalConstants.relativistic_parameters(100.0, PhysicalConstants.E0_electron)
    assert abs(sim._gamma - gamma_100) < 1e-6

    print("Relativistic parameters: PASSED")


def test_current_to_k1():
    """Test quadrupole current to k1 conversion."""
    print("\n=== Current to k1 Conversion ===")

    sim = RFTrackAdapter(beam_energy=BEAM_ENERGY)

    k1_focus = sim._current_to_k1(1.0, 0.1, focusing=True)
    k1_defocus = sim._current_to_k1(1.0, 0.1, focusing=False)

    print(f"k1 (I=1A, focusing): {k1_focus:.4f} 1/m²")
    print(f"k1 (I=1A, defocusing): {k1_defocus:.4f} 1/m²")

    assert k1_focus > 0
    assert k1_defocus < 0
    assert abs(k1_focus) == abs(k1_defocus)

    # Verify against manual calculation
    mass_kg = sim.particle_mass * PhysicalConstants.MeV_to_J / PhysicalConstants.C**2
    k1_manual = abs(PhysicalConstants.Q * sim.G_quad * 1.0) / (
        mass_kg * PhysicalConstants.C * sim._beta * sim._gamma
    )
    assert abs(k1_focus - k1_manual) < 1e-10

    assert sim._current_to_k1(0.0, 0.1, focusing=True) == 0.0
    assert sim._current_to_k1(1.0, 0.0, focusing=True) == 0.0

    print("Current to k1 conversion: PASSED")


def test_coordinate_transforms():
    """Test coordinate transformations."""
    print("\n=== Coordinate Transformations ===")

    sim = RFTrackAdapter(beam_energy=BEAM_ENERGY)
    particles = np.array([
        [1.0, 0.1, 2.0, 0.2, 0.5, 0.3],
        [-0.5, -0.05, 1.0, -0.1, -0.2, 0.1],
    ])

    # FELsim -> RF-Track -> FELsim
    p_rft = sim.transform_coordinates(particles, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK)
    p_back = sim.transform_coordinates(p_rft, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM)

    error = np.max(np.abs(particles - p_back))
    print(f"FELsim round-trip error: {error:.2e}")
    assert error < 1e-10

    # COSY -> RF-Track -> COSY
    particles_cosy = np.array([[1e-3, 1e-4, 2e-3, 2e-4, 1e-3, 1e-4]])
    p_rft2 = sim.transform_coordinates(particles_cosy, CoordinateSystem.COSY, CoordinateSystem.RFTRACK)
    p_back2 = sim.transform_coordinates(p_rft2, CoordinateSystem.RFTRACK, CoordinateSystem.COSY)

    error2 = np.max(np.abs(particles_cosy - p_back2))
    print(f"COSY round-trip error: {error2:.2e}")
    assert error2 < 1e-10

    # Identity transform
    p_same = sim.transform_coordinates(particles, CoordinateSystem.FELSIM, CoordinateSystem.FELSIM)
    assert np.allclose(particles, p_same)

    print("Coordinate transformations: PASSED")


def test_beamline_setup():
    """Test beamline configuration."""
    print("\n=== Beamline Setup ===")

    sim = RFTrackAdapter(beam_energy=BEAM_ENERGY)

    elements = [
        BeamlineElement('DRIFT', 0.5),
        BeamlineElement('QUAD_F', 0.1, current=1.5),
        BeamlineElement('DRIFT', 0.3),
        BeamlineElement('QUAD_D', 0.1, current=1.5),
        BeamlineElement('DRIFT', 0.5),
    ]

    sim.set_beamline(elements)

    assert len(sim.beamline) == 5
    assert sim._lattice is not None
    assert sim._lattice.size() == 5

    total_length = sim._lattice.get_length()
    expected_length = sum(e.length for e in elements)
    print(f"Lattice length: {total_length:.4f} m (expected: {expected_length:.4f} m)")
    assert abs(total_length - expected_length) < 1e-6

    print("Beamline setup: PASSED")


def test_particle_generation():
    """Test particle distribution generation."""
    print("\n=== Particle Generation ===")

    sim = RFTrackAdapter(beam_energy=BEAM_ENERGY)

    # Gaussian
    p_gauss = sim.generate_particles(1000, distribution_type='gaussian',
                                      std_dev=[1.0, 0.1, 1.0, 0.1, 0.5, 0.2])
    assert p_gauss.shape == (1000, 6)
    print(f"Gaussian: shape={p_gauss.shape}, σx={np.std(p_gauss[:, 0]):.3f} mm")

    # Uniform
    p_uniform = sim.generate_particles(500, distribution_type='uniform',
                                        std_dev=[2.0, 0.2, 2.0, 0.2, 1.0, 0.5])
    assert p_uniform.shape == (500, 6)
    print(f"Uniform: shape={p_uniform.shape}, x_range=[{p_uniform[:, 0].min():.2f}, {p_uniform[:, 0].max():.2f}] mm")

    # Twiss-matched
    p_twiss = sim.generate_particles(500, distribution_type='twiss',
                                      twiss_x={'beta': 5.0, 'alpha': -1.0, 'emittance': 2.0},
                                      twiss_y={'beta': 8.0, 'alpha': 0.5, 'emittance': 2.0})
    assert p_twiss.shape == (500, 6)
    print(f"Twiss-matched: shape={p_twiss.shape}")

    print("Particle generation: PASSED")


def test_simulation():
    """Test particle tracking simulation."""
    print("\n=== Particle Tracking ===")

    sim = RFTrackAdapter(beam_energy=BEAM_ENERGY)

    elements = [
        BeamlineElement('DRIFT', 0.3),
        BeamlineElement('QUAD_F', 0.08, current=2.0),
        BeamlineElement('DRIFT', 0.2),
        BeamlineElement('QUAD_D', 0.08, current=2.0),
        BeamlineElement('DRIFT', 0.3),
    ]
    sim.set_beamline(elements)

    particles = sim.generate_particles(200, std_dev=[0.5, 0.05, 0.5, 0.05, 0.1, 0.1])

    result = sim.simulate(particles)

    assert result.success
    assert result.simulator_name == "RF-Track"
    assert result.final_particles is not None
    assert result.final_particles.shape[0] > 0

    print(f"Success: {result.success}")
    print(f"Particles: {result.metadata['num_good']} good, {result.metadata['num_lost']} lost")
    print(f"Lattice length: {result.metadata['lattice_length']:.4f} m")

    twiss = result.twiss_parameters_statistical['final']
    print(f"Final βx: {twiss['x']['beta']:.4f} m")
    print(f"Final βy: {twiss['y']['beta']:.4f} m")
    print(f"Final εx: {twiss['x']['emittance']:.4f} π·mm·mrad")

    assert 'x' in twiss and 'y' in twiss
    assert all(k in twiss['x'] for k in ['beta', 'alpha', 'gamma', 'emittance'])

    print("Particle tracking: PASSED")


def test_twiss_calculation():
    """Test Twiss parameter calculation from particles."""
    print("\n=== Twiss Calculation ===")

    sim = RFTrackAdapter(beam_energy=BEAM_ENERGY)

    # Generate Gaussian beam
    particles = sim.generate_particles(2000, std_dev=[1.0, 0.1, 1.0, 0.1, 0.5, 0.2])
    twiss = sim._calculate_twiss(particles)

    print(f"Calculated Twiss from Gaussian beam:")
    print(f"  βx={twiss['x']['beta']:.4f} m, αx={twiss['x']['alpha']:.4f}")
    print(f"  βy={twiss['y']['beta']:.4f} m, αy={twiss['y']['alpha']:.4f}")
    print(f"  εx={twiss['x']['emittance']:.4f}, εy={twiss['y']['emittance']:.4f}")

    # Basic sanity checks
    assert twiss['x']['beta'] > 0
    assert twiss['y']['beta'] > 0
    assert twiss['x']['emittance'] > 0
    assert twiss['y']['emittance'] > 0
    assert all(k in twiss['x'] for k in ['beta', 'alpha', 'gamma', 'emittance'])

    print("Twiss calculation: PASSED")


def test_g_quad_modification():
    """Test runtime modification of quadrupole gradient."""
    print("\n=== G_quad Modification ===")

    sim = RFTrackAdapter(beam_energy=BEAM_ENERGY)
    elements = [
        BeamlineElement('DRIFT', 0.2),
        BeamlineElement('QUAD_F', 0.1, current=1.0),
        BeamlineElement('DRIFT', 0.2),
    ]
    sim.set_beamline(elements)

    k1_default = sim._current_to_k1(1.0, 0.1, focusing=True)

    sim.set_quadrupole_gradient(5.0)
    assert sim.G_quad == 5.0

    k1_new = sim._current_to_k1(1.0, 0.1, focusing=True)
    ratio = k1_new / k1_default
    expected_ratio = 5.0 / PhysicalConstants.G_quad_default

    print(f"k1 ratio: {ratio:.4f} (expected: {expected_ratio:.4f})")
    assert abs(ratio - expected_ratio) < 1e-6

    print("G_quad modification: PASSED")


