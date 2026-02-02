"""
Test unified plotting with RF-Track adapter.

Author: Eremey Valetov
"""
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import traceback
from pathlib import Path

# Check RF-Track availability
try:
    from rftrackAdapter import RFTrackAdapter, _RFTRACK_AVAILABLE
    if not _RFTRACK_AVAILABLE:
        print("RF-Track not installed (pip install RF-Track)")
        sys.exit(1)
except ImportError as e:
    print(f"Import error: {e}")
    sys.exit(1)

from beamEvolution import BeamEvolution
from evolutionPlotter import EvolutionPlotter
from simulatorBase import CoordinateSystem, SimulationMode
from simulatorFactory import SimulatorFactory

# Import FELsim for comparison (optional)
try:
    from felsimAdapter import FELsimAdapter
    from ebeam import beam
    FELSIM_AVAILABLE = True
except ImportError:
    FELSIM_AVAILABLE = False

# Configuration
EXCEL_PATH = Path("../../beam_excel/Beamline_elements_3.xlsx")
BEAM_ENERGY = 45.0
MAX_ELEMENTS = 30  # Limit elements for faster testing


def create_test_particles(n=1000):
    """Generate test distribution in FELsim coordinates."""
    if FELSIM_AVAILABLE:
        ebeam = beam()
        std_dev = [1.0, 0.1, 1.0, 0.1, 1.0, 0.5]
        return ebeam.gen_6d_gaussian(0, std_dev, n)
    else:
        std_dev = np.array([1.0, 0.1, 1.0, 0.1, 1.0, 0.5])
        return np.random.randn(n, 6) * std_dev


def test_rftrack_adapter_creation():
    """Test RFTrackAdapter instantiation with Excel beamline."""
    print("\n=== RFTrackAdapter Creation ===")

    if not EXCEL_PATH.exists():
        print(f"Excel file not found: {EXCEL_PATH}")
        print("Creating adapter with manual beamline")
        sim = RFTrackAdapter(beam_energy=BEAM_ENERGY)
        _create_manual_beamline(sim)
    else:
        print(f"Loading beamline from: {EXCEL_PATH}")
        sim = RFTrackAdapter(
            excel_path=str(EXCEL_PATH),
            beam_energy=BEAM_ENERGY
        )

    print(f"Simulator: {sim.name}")
    print(f"Native coordinates: {sim.native_coordinates.value}")
    print(f"Mode: {sim.simulation_mode}")
    print(f"G_quad: {sim.G_quad:.4f} T/A/m")

    n_elements = len(sim.beamline)
    print(f"Beamline: {n_elements} elements")

    # Limit elements for faster testing
    if n_elements > MAX_ELEMENTS:
        sim.beamline = sim.beamline[:MAX_ELEMENTS]
        sim._build_lattice()
        print(f"Limited to {MAX_ELEMENTS} elements for testing")

    total_length = sim._lattice.get_length()
    print(f"Lattice length: {total_length:.4f} m")

    for i, elem in enumerate(sim.beamline[:5]):
        print(f"  [{i}] {elem.element_type}: L={elem.length:.4f} m")
    if len(sim.beamline) > 5:
        print(f"  ... and {len(sim.beamline) - 5} more elements")

    print("RFTrackAdapter created successfully")
    return sim


def _create_manual_beamline(sim):
    """Create a simple FODO beamline for testing."""
    from simulatorBase import BeamlineElement

    elements = []
    for _ in range(5):  # 5 FODO cells
        elements.extend([
            BeamlineElement('DRIFT', 0.2),
            BeamlineElement('QUAD_F', 0.08, current=2.0),
            BeamlineElement('DRIFT', 0.4),
            BeamlineElement('QUAD_D', 0.08, current=2.0),
            BeamlineElement('DRIFT', 0.2),
        ])
    sim.set_beamline(elements)


def test_particle_tracking(sim):
    """Test basic particle tracking."""
    print("\n=== Particle Tracking ===")

    particles = create_test_particles(500)
    print(f"Input particles: {particles.shape}")
    print(f"  x range: [{particles[:, 0].min():.4f}, {particles[:, 0].max():.4f}] mm")

    result = sim.simulate(particles)

    print(f"Simulation success: {result.success}")
    print(f"Particles: {result.metadata['num_good']} good, {result.metadata['num_lost']} lost")

    if result.final_particles is not None:
        final = result.final_particles
        print(f"Output particles: {final.shape}")
        print(f"  x range: [{final[:, 0].min():.4f}, {final[:, 0].max():.4f}] mm")

    return particles


def test_collect_evolution(sim, particles):
    """Test RFTrackAdapter.collect_evolution()."""
    print("\n=== Evolution Collection ===")

    try:
        n_elements = len(sim.beamline)
        checkpoint_list = list(range(min(n_elements, 20)))

        print(f"Running RF-Track simulation with {len(checkpoint_list)} checkpoints...")
        evolution = sim.collect_evolution(
            particles,
            checkpoint_elements=checkpoint_list
        )

        print(f"Evolution data:")
        print(f"  {len(evolution.s_positions)} s-positions")
        print(f"  {len(evolution.particles)} particle snapshots")
        print(f"  {len(evolution.twiss)} Twiss calculations")
        print(f"  {len(evolution.elements)} elements")
        print(f"  Total length: {evolution.total_length:.4f} m")

        # Verify data consistency
        assert len(evolution.s_positions) == len(evolution.particles), \
            f"Mismatch: {len(evolution.s_positions)} positions vs {len(evolution.particles)} particles"
        assert len(evolution.s_positions) == len(evolution.twiss), \
            f"Mismatch: {len(evolution.s_positions)} positions vs {len(evolution.twiss)} twiss"
        print("Data structures consistent")

        # Check Twiss fields
        if evolution.s_positions:
            s0 = evolution.s_positions[0]
            if s0 in evolution.twiss:
                twiss_x = evolution.twiss[s0].get('x', {})
                required_fields = ['beta', 'alpha', 'gamma', 'emittance']
                missing = [f for f in required_fields if f not in twiss_x]
                if missing:
                    print(f"Warning: missing Twiss fields: {missing}")
                else:
                    print("Twiss contains required fields")
                    print(f"  βx = {twiss_x['beta']:.4f} m")
                    print(f"  αx = {twiss_x['alpha']:.4f}")

        # Show s-position statistics
        if len(evolution.s_positions) > 1:
            spacings = np.diff(sorted(evolution.s_positions))
            print(f"  s-spacing: min={spacings.min():.4f}, max={spacings.max():.4f}, mean={spacings.mean():.4f} m")

        print("collect_evolution() completed successfully")
        return evolution

    except Exception as e:
        print(f"collect_evolution() failed: {e}")
        traceback.print_exc()
        return None


def test_twiss_dataframe(evolution):
    """Test BeamEvolution DataFrame export."""
    print("\n=== Twiss DataFrame Export ===")

    if evolution is None:
        print("Skipped: no evolution data")
        return None

    df = evolution.get_twiss_evolution()

    print(f"DataFrame shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    print(f"s range: [{df['s'].min():.4f}, {df['s'].max():.4f}] m")

    if 'beta_x' in df.columns and not df['beta_x'].isna().all():
        print(f"βx range: [{df['beta_x'].min():.4f}, {df['beta_x'].max():.4f}] m")
        print(f"Envelope_x range: [{df['envelope_x'].min():.4f}, {df['envelope_x'].max():.4f}] mm")

    return df


def test_plotter(evolution):
    """Test EvolutionPlotter with RF-Track data."""
    print("\n=== Evolution Plotting ===")

    if evolution is None:
        print("Skipped: no evolution data")
        return False

    if len(evolution.particles) == 0:
        print("Skipped: no particle data in evolution")
        return False

    plotter = EvolutionPlotter(axis_mode='local')

    print("Plotting RF-Track evolution (close window to continue)...")
    try:
        plotter.plot(
            evolution,
            show_phase_space=True,
            show_envelope=True,
            show_schematic=True,
            interactive=True,
            scatter=False,
            envelope_ylim=(0.0, 1.4)
        )
        print("Plot displayed successfully")
        return True
    except Exception as e:
        print(f"Plotting failed: {e}")
        traceback.print_exc()
        return False


def test_comparison_with_felsim(sim_rftrack, particles):
    """Compare RF-Track evolution with FELsim on same initial conditions."""
    print("\n=== Comparison with FELsim ===")

    if not FELSIM_AVAILABLE:
        print("Skipped: FELsim not available")
        return None

    try:
        print("Running RF-Track simulation...")
        evolution_rftrack = sim_rftrack.collect_evolution(
            particles.copy(),
            checkpoint_elements='all'
        )

        if evolution_rftrack is None or len(evolution_rftrack.twiss) == 0:
            print("RF-Track evolution returned no data")
            return None

        final_s = max(evolution_rftrack.s_positions)
        twiss_rftrack = evolution_rftrack.twiss.get(final_s, {})

        if 'x' in twiss_rftrack:
            print(f"\nRF-Track final Twiss (s={final_s:.4f} m):")
            print(f"  βx = {twiss_rftrack['x'].get('beta', 'N/A'):.4f} m")
            print(f"  βy = {twiss_rftrack['y'].get('beta', 'N/A'):.4f} m")
            print(f"  εx = {twiss_rftrack['x'].get('emittance', 'N/A'):.6f} π·mm·mrad")

        print("\nRF-Track simulation completed")
        return evolution_rftrack

    except Exception as e:
        print(f"Comparison failed: {e}")
        traceback.print_exc()
        return None


def run_all_tests():
    """Run complete RF-Track unified plotting test suite."""
    print("=" * 60)
    print("RF-Track Unified Plotting Test Suite")
    print("=" * 60)

    try:
        sim = test_rftrack_adapter_creation()
        particles = test_particle_tracking(sim)
        evolution = test_collect_evolution(sim, particles)
        test_twiss_dataframe(evolution)
        test_plotter(evolution)
        test_comparison_with_felsim(sim, particles)

        print("\n" + "=" * 60)
        print("RF-Track unified test suite completed")
        print("=" * 60)
        return True

    except Exception as e:
        print(f"\nTest suite failed: {e}")
        traceback.print_exc()
        return False


def visual_only():
    """Run only the visual plotting test."""
    print("\n=== RF-Track Visual Plot ===")

    if EXCEL_PATH.exists():
        sim = RFTrackAdapter(excel_path=str(EXCEL_PATH), beam_energy=BEAM_ENERGY)
        if len(sim.beamline) > MAX_ELEMENTS:
            sim.beamline = sim.beamline[:MAX_ELEMENTS]
            sim._build_lattice()
    else:
        sim = RFTrackAdapter(beam_energy=BEAM_ENERGY)
        _create_manual_beamline(sim)

    particles = create_test_particles(1000)

    print(f"Collecting evolution for {len(sim.beamline)} elements...")
    evolution = sim.collect_evolution(particles, checkpoint_elements='all')

    print("Displaying RF-Track plot...")
    plotter = EvolutionPlotter(axis_mode='local')
    plotter.plot(evolution, envelope_ylim=(0.0, 1.4))

    print("Visual test complete")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Test RF-Track unified plotting")
    parser.add_argument('--visual', action='store_true', help='Run visual plot only')
    parser.add_argument('--excel', type=str, default=str(EXCEL_PATH), help='Path to beamline Excel file')
    parser.add_argument('--max', type=int, default=MAX_ELEMENTS, help='Maximum number of elements')

    args = parser.parse_args()
    EXCEL_PATH = Path(args.excel)
    MAX_ELEMENTS = args.max

    if args.visual:
        visual_only()
    else:
        run_all_tests()
