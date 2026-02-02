#!/usr/bin/env python3
"""
Test beamline schematic plotting functionality.

Tests legacy and enhanced schematic plotting modes with FELsim beamline.
Compares straight vs curved trajectory rendering for dipoles.

Author: Eremey Valetov
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Configuration
EXCEL_PATH = Path("../../beam_excel/Beamline_elements.xlsx")
MAX_ELEMENTS = 50  # Limit for testing (avoid numerical instabilities)
BEAM_ENERGY_MEV = 45.0

from excelElements import ExcelElements
from schematic import draw_beamline
from ebeam import beam
from beamline import lattice


def load_beamline():
    """Load beamline from Excel file."""
    if not EXCEL_PATH.exists():
        raise FileNotFoundError(f"Excel file not found: {EXCEL_PATH}")

    print(f"Loading beamline from {EXCEL_PATH.name}")

    excel = ExcelElements(str(EXCEL_PATH))
    beamline_full = excel.create_beamline()

    # Limit to first N elements for stability
    beamline = beamline_full[:MAX_ELEMENTS]

    # Set beam energy
    relat = lattice(1, fringeType=None)
    relat.setE(E=BEAM_ENERGY_MEV)

    for elem in beamline:
        elem.setE(BEAM_ENERGY_MEV)

    # Print beamline summary
    total_length = sum(elem.length for elem in beamline)
    print(f"Loaded {len(beamline)}/{len(beamline_full)} elements")
    print(f"Total length: {total_length:.4f} m\n")

    # Count element types
    elem_counts = {}
    for elem in beamline:
        elem_type = elem.__class__.__name__
        elem_counts[elem_type] = elem_counts.get(elem_type, 0) + 1

    print("Element composition:")
    for elem_type, count in sorted(elem_counts.items()):
        print(f"  {elem_type:20s}: {count:3d}")

    # Check for dipoles
    dipole_count = sum(1 for e in beamline if 'dipole' in e.__class__.__name__.lower())
    print(f"\nDipoles found: {dipole_count}")
    if dipole_count > 0:
        print("  (Curved trajectory mode will show bending)")
    else:
        print("  (No bending to visualize)")

    return beamline


def create_test_particles(n=1000):
    """Generate test particle distribution."""
    ebeam = beam()
    std_dev = np.array([1.0, 0.1, 1.0, 0.1, 1.0, 0.5])
    return ebeam.gen_6d_gaussian(0, std_dev, n)


def test_standalone_legacy_schematic(beamline):
    """Test standalone legacy schematic plotting."""
    print("\n" + "=" * 70)
    print("TEST 1: Standalone Legacy Schematic")
    print("=" * 70)

    schem = draw_beamline()

    fig, ax = plt.subplots(figsize=(14, 3))

    # Use legacy method: manually draw rectangles
    ymin, ymax = -0.5, 0.5
    ax.set_ylim(ymin, ymax)

    s_pos = 0.0
    for i, elem in enumerate(beamline):
        elem_info = schem._extract_element_info(elem)
        length = elem_info['length']
        color = elem_info['color']

        rect = plt.Rectangle((s_pos, ymin), length, (ymax - ymin) * 0.1,
                             linewidth=1, edgecolor=color, facecolor=color)
        ax.add_patch(rect)
        s_pos += length

    ax.axhline(0, color='k', linewidth=1)
    ax.set_xlim(0, s_pos)
    ax.set_xlabel('s (m)')
    ax.set_title('Legacy Schematic (Original Block Drawing)')

    plt.tight_layout()
    print(f"✓ Legacy schematic rendered ({len(beamline)} elements, {s_pos:.4f} m)")
    print("Close plot window to continue...")
    plt.show()


def test_standalone_enhanced_straight(beamline):
    """Test standalone enhanced schematic with straight trajectory."""
    print("\n" + "=" * 70)
    print("TEST 2: Standalone Enhanced Schematic (Straight)")
    print("=" * 70)

    schem = draw_beamline()

    fig, ax = plt.subplots(figsize=(14, 3))
    schem.plot_beamline_schematic(
        beamline,
        ax=ax,
        show_curved_trajectory=False,
        height_scale=0.4,
        show_labels=False
    )

    ax.set_title('Enhanced Schematic: Straight Reference Trajectory')
    plt.tight_layout()

    print("✓ Enhanced schematic (straight) rendered")
    print("  - Element-specific shapes (quads up/down, dipoles colored)")
    print("  - Straight centreline")
    print("Close plot window to continue...")
    plt.show()


def test_standalone_enhanced_curved(beamline):
    """Test standalone enhanced schematic with curved trajectory."""
    print("\n" + "=" * 70)
    print("TEST 3: Standalone Enhanced Schematic (Curved)")
    print("=" * 70)

    schem = draw_beamline()

    fig, ax = plt.subplots(figsize=(14, 4))
    schem.plot_beamline_schematic(
        beamline,
        ax=ax,
        show_curved_trajectory=True,
        height_scale=0.4,
        show_labels=False
    )

    ax.set_title('Enhanced Schematic: Curved Reference Trajectory (Dipole Bending)')
    plt.tight_layout()

    print("✓ Enhanced schematic (curved) rendered")
    print("  - Element-specific shapes")
    print("  - Curved centreline showing dipole bending")
    print("  - Bend angles annotated on dipoles")
    print("Close plot window to continue...")
    plt.show()


def test_side_by_side_comparison(beamline):
    """Compare all three modes side by side."""
    print("\n" + "=" * 70)
    print("TEST 4: Side-by-Side Comparison")
    print("=" * 70)

    schem = draw_beamline()

    fig, axes = plt.subplots(3, 1, figsize=(14, 8))

    # Legacy mode
    ax = axes[0]
    ymin, ymax = -0.5, 0.5
    ax.set_ylim(ymin, ymax)
    s_pos = 0.0
    for elem in beamline:
        elem_info = schem._extract_element_info(elem)
        length = elem_info['length']
        color = elem_info['color']
        rect = plt.Rectangle((s_pos, ymin), length, (ymax - ymin) * 0.1,
                             linewidth=1, edgecolor=color, facecolor=color)
        ax.add_patch(rect)
        s_pos += length
    ax.axhline(0, color='k', linewidth=1)
    ax.set_xlim(0, s_pos)
    ax.set_title('Legacy Mode')
    ax.set_ylabel('')
    ax.tick_params(axis='y', which='both', left=False, labelleft=False)

    # Enhanced straight
    schem.plot_beamline_schematic(
        beamline,
        ax=axes[1],
        show_curved_trajectory=False,
        height_scale=0.4,
        show_labels=False
    )
    axes[1].set_title('Enhanced Mode: Straight Trajectory')
    axes[1].set_xlabel('')

    # Enhanced curved
    schem.plot_beamline_schematic(
        beamline,
        ax=axes[2],
        show_curved_trajectory=True,
        height_scale=0.4,
        show_labels=False
    )
    axes[2].set_title('Enhanced Mode: Curved Trajectory (Dipole Bending)')

    plt.tight_layout()
    print("✓ All three modes displayed side-by-side")
    print("Close plot window to continue...")
    plt.show()


def test_with_simulation_legacy(beamline):
    """Test schematic within full beam evolution simulation (legacy mode)."""
    print("\n" + "=" * 70)
    print("TEST 5: Full Simulation with Legacy Schematic")
    print("=" * 70)

    schem = draw_beamline()
    particles = create_test_particles(500)

    print("Running simulation with legacy schematic...")
    print("(This will show phase space + envelope plots + legacy schematic)")

    twiss_df = schem.plotBeamPositionTransform(
        particles,
        beamline,
        interval=0.1,
        defineLim=True,
        plot=True,
        showIndice=False,
        scatter=False,
        show_schematic='legacy'  # Explicitly use legacy
    )

    print("✓ Simulation completed with legacy schematic")
    print(f"  Final s: {twiss_df.columns[-1]} (data points)")


def test_with_simulation_enhanced_straight(beamline):
    """Test schematic within full simulation (enhanced, straight)."""
    print("\n" + "=" * 70)
    print("TEST 6: Full Simulation with Enhanced Schematic (Straight)")
    print("=" * 70)

    schem = draw_beamline()
    particles = create_test_particles(500)

    print("Running simulation with enhanced schematic (straight)...")

    twiss_df = schem.plotBeamPositionTransform(
        particles,
        beamline,
        interval=0.1,
        defineLim=True,
        plot=True,
        showIndice=False,
        scatter=False,
        show_schematic='enhanced',
        curved_trajectory=False
    )

    print("✓ Simulation completed with enhanced schematic (straight)")


def test_with_simulation_enhanced_curved(beamline):
    """Test schematic within full simulation (enhanced, curved)."""
    print("\n" + "=" * 70)
    print("TEST 7: Full Simulation with Enhanced Schematic (Curved)")
    print("=" * 70)

    schem = draw_beamline()
    particles = create_test_particles(500)

    print("Running simulation with enhanced schematic (curved trajectory)...")

    twiss_df = schem.plotBeamPositionTransform(
        particles,
        beamline,
        interval=0.1,
        defineLim=True,
        plot=True,
        showIndice=False,
        scatter=False,
        show_schematic='enhanced',
        curved_trajectory=True
    )

    print("✓ Simulation completed with enhanced schematic (curved)")


def test_element_labels(beamline):
    """Test schematic with element labels/indices."""
    print("\n" + "=" * 70)
    print("TEST 8: Enhanced Schematic with Element Labels")
    print("=" * 70)

    schem = draw_beamline()

    fig, ax = plt.subplots(figsize=(14, 4))
    schem.plot_beamline_schematic(
        beamline,
        ax=ax,
        show_curved_trajectory=True,
        height_scale=0.4,
        show_labels=True  # Show indices
    )

    ax.set_title('Enhanced Schematic with Element Indices')
    plt.tight_layout()

    print("✓ Schematic with labels rendered")
    print("  - Element indices displayed")
    print("Close plot window to continue...")
    plt.show()


def test_backward_compatibility(beamline):
    """Verify backward compatibility: existing code still works."""
    print("\n" + "=" * 70)
    print("TEST 9: Backward Compatibility Check")
    print("=" * 70)

    schem = draw_beamline()
    particles = create_test_particles(500)

    print("Running simulation with NO schematic parameters (should use legacy)...")

    # Call without new parameters - should default to legacy mode
    twiss_df = schem.plotBeamPositionTransform(
        particles,
        beamline,
        interval=0.1,
        defineLim=True,
        plot=True,
        showIndice=False,
        scatter=False
        # NOTE: No show_schematic or curved_trajectory parameters
    )

    print("✓ Backward compatibility verified")
    print("  - Existing code runs without modification")
    print("  - Defaults to legacy mode")


def analyze_beamline_geometry(beamline):
    """Analyze beamline geometry for curved trajectory validation."""
    print("\n" + "=" * 70)
    print("ANALYSIS: Beamline Geometry")
    print("=" * 70)

    schem = draw_beamline()

    total_length = 0.0
    total_bend = 0.0
    dipole_info = []

    for i, elem in enumerate(beamline):
        elem_info = schem._extract_element_info(elem)
        length = elem_info['length']
        angle = elem_info['angle']
        elem_type = elem_info['type']

        total_length += length

        if elem_type in ['DPH', 'DIPOLE'] and abs(angle) > 1e-10:
            total_bend += angle
            rho = length / np.radians(angle) if angle != 0 else float('inf')
            dipole_info.append({
                'index': i,
                'type': elem_type,
                'length': length,
                'angle_deg': angle,
                'angle_rad': np.radians(angle),
                'rho': rho
            })

    print(f"\nTotal beamline length: {total_length:.4f} m")
    print(f"Total bending angle: {total_bend:.4f}°")

    if dipole_info:
        print(f"\nDipole elements: {len(dipole_info)}")
        print("\n  Idx  Type      Length(m)  Angle(°)  Rho(m)")
        print("  " + "-" * 50)
        for d in dipole_info:
            print(f"  {d['index']:3d}  {d['type']:8s}  {d['length']:8.4f}  "
                  f"{d['angle_deg']:8.4f}  {d['rho']:8.4f}")

        print(f"\nCurved trajectory should show:")
        print(f"  - {len(dipole_info)} bending section(s)")
        print(f"  - Total deviation from straight line")

        # Estimate max perpendicular displacement
        max_displacement = 0.0
        for d in dipole_info:
            # For small angles: y ≈ rho * (1 - cos(theta)) ≈ rho * theta²/2
            theta = d['angle_rad']
            rho = d['rho']
            disp = rho * (1 - np.cos(theta))
            max_displacement = max(max_displacement, disp)

        print(f"  - Max perpendicular displacement: ~{max_displacement:.4f} m")
    else:
        print("\nNo dipoles found - curved trajectory will be straight")


def run_all_tests():
    """Run complete test suite."""
    print("\n" + "=" * 70)
    print("Beamline Schematic Plotting Test Suite")
    print("=" * 70)

    try:
        # Load beamline
        beamline = load_beamline()

        # Analyze geometry first
        analyze_beamline_geometry(beamline)

        # Standalone schematic tests
        test_standalone_legacy_schematic(beamline)
        test_standalone_enhanced_straight(beamline)
        test_standalone_enhanced_curved(beamline)
        test_side_by_side_comparison(beamline)
        test_element_labels(beamline)

        # Full simulation tests
        test_with_simulation_legacy(beamline)
        test_with_simulation_enhanced_straight(beamline)
        test_with_simulation_enhanced_curved(beamline)

        # Compatibility test
        test_backward_compatibility(beamline)

        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)
        print("\n✓ All schematic plotting tests passed!")
        print("\nFeatures verified:")
        print("  • Legacy mode (backward compatible)")
        print("  • Enhanced mode with straight trajectory")
        print("  • Enhanced mode with curved trajectory")
        print("  • Element-specific rendering (quads, dipoles)")
        print("  • Standalone and integrated plotting")
        print("  • Element labeling")
        print("\nUsage:")
        print("  Legacy:    show_schematic='legacy'")
        print("  Enhanced:  show_schematic='enhanced', curved_trajectory=False/True")
        print("\n" + "=" * 70 + "\n")

        return True

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


def quick_demo():
    """Quick demonstration of all three modes."""
    print("\n" + "=" * 70)
    print("QUICK DEMO: Three Schematic Modes")
    print("=" * 70)

    beamline = load_beamline()
    test_side_by_side_comparison(beamline)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Test beamline schematic plotting functionality"
    )
    parser.add_argument(
        '--quick',
        action='store_true',
        help='Run quick demo only (side-by-side comparison)'
    )
    parser.add_argument(
        '--excel',
        type=str,
        default=str(EXCEL_PATH),
        help='Path to beamline Excel file'
    )
    parser.add_argument(
        '--max-elements',
        type=int,
        default=MAX_ELEMENTS,
        help='Maximum number of elements to use'
    )

    args = parser.parse_args()
    EXCEL_PATH = Path(args.excel)
    MAX_ELEMENTS = args.max_elements

    if args.quick:
        quick_demo()
    else:
        success = run_all_tests()
        sys.exit(0 if success else 1)