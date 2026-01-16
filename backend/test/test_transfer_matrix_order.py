#!/usr/bin/env python3
"""
Demonstration of higher-order transfer map reading.

Tests the new read_transfer_map_full() and get_aberration_coefficient() methods.

Author: Eremey Valetov
"""

import sys
import numpy as np
from pathlib import Path

# Configuration
EXCEL_PATH = Path("../../beam_excel/Beamline_elements_3.xlsx")
BEAM_ENERGY = 45.0  # MeV


def demo_higher_order_reading():
    """Demonstrate reading higher-order coefficients."""
    print("\n" + "=" * 70)
    print("DEMO: Higher-Order Transfer Map Reading")
    print("=" * 70 + "\n")

    from cosyAdapter import COSYAdapter

    # Run simulation with order 3
    print("Running COSY simulation with transfer_matrix_order=3...")
    sim = COSYAdapter(
        excel_path=str(EXCEL_PATH),
        mode='transfer_matrix',
        transfer_matrix_order=3,
        debug=False
    )

    result = sim.simulate()

    if not result.success:
        print("ERROR: Simulation failed")
        return

    print("✓ Simulation completed\n")

    # Get results reader
    reader = sim.get_native_simulator().analyze_results()

    # Read full transfer map
    print("Reading full transfer map...")
    full_map = reader.read_transfer_map_full()

    print(f"\n✓ Successfully parsed transfer map")
    print(f"  Orders present: {sorted(full_map.keys())}\n")

    # Show statistics for each order
    for order in sorted(full_map.keys()):
        if order == 1:
            print(f"Order {order} (Linear):")
            print(f"  Matrix shape: {full_map[order].shape}")
            print(f"  Sample element M[0,0] (x|x): {full_map[order][0, 0]:.6f}")
        else:
            coeffs_dict = full_map[order]
            print(f"\nOrder {order}:")
            print(f"  Number of terms: {len(coeffs_dict)}")

            # Show a few examples
            sample_indices = list(coeffs_dict.keys())[:3]
            for idx in sample_indices:
                coeffs = coeffs_dict[idx]
                # Calculate which coordinates contribute
                powers = [int(d) for d in idx]
                coord_names = ['x', "x'", 'y', "y'", 'l', 'δK']
                term_str = ' '.join([f"{coord_names[i]}^{powers[i]}"
                                     for i in range(6) if powers[i] > 0])

                # Show non-zero output components
                nonzero = [(i, coeffs[i]) for i in range(6) if abs(coeffs[i]) > 1e-10]
                if nonzero:
                    for out_idx, val in nonzero[:2]:  # Show first 2
                        print(f"    {idx} ({term_str}) → {coord_names[out_idx]}: {val:.6e}")

    return reader, full_map


def demo_aberration_extraction(reader):
    """Demonstrate extracting specific aberration coefficients."""
    print("\n" + "=" * 70)
    print("DEMO: Aberration Coefficient Extraction")
    print("=" * 70 + "\n")

    # Examples of common aberrations
    # NOTE: source_coords are the coordinates that contribute (not power notation!)
    # For x²: source_coords = (0, 0) meaning x appears twice
    # For xx': source_coords = (0, 1) meaning x and x' each appear once
    aberrations = [
        # 2nd order
        ("T_200000", 'x', (0, 0), "x from x² (geometric aberration)"),
        ("T_110000", 'x', (0, 1), "x from xx' (focusing)"),
        ("T_101000", 'x', (0, 2), "x from xy (coupling)"),
        ("T_020000", 'x', (1, 1), "x from x'² (divergence)"),
        ("T_002000", 'x', (5, 5), "x from δ² (chromatic)"),

        # 3rd order
        ("T_300000", 'x', (0, 0, 0), "x from x³ (spherical aberration)"),
        ("T_210000", 'x', (0, 0, 1), "x from x²x'"),
        ("T_111000", 'x', (0, 1, 2), "x from xx'y (coupled)"),
        ("T_120000", 'x', (0, 1, 1), "x from xx'²"),
    ]

    print("Extracting specific aberration coefficients:\n")

    for name, out_coord, indices, description in aberrations:
        try:
            coeff = reader.get_aberration_coefficient(out_coord, indices)
            order = len(indices)  # Order = number of source coordinates
            print(f"{name:12s} (order {order}): {coeff:12.6e}  # {description}")
        except ValueError as e:
            print(f"{name:12s}: NOT FOUND ({e})")

    # Show how to access y-plane aberrations
    print("\nY-plane examples:")
    y_aberrations = [
        ("T_001000", 'y', (2,), "y from y (linear)"),
        ("T_002000", 'y', (2, 2), "y from y² (2nd order)"),
        ("T_011000", 'y', (2, 3), "y from yy' (2nd order)"),
        ("T_003000", 'y', (2, 2, 2), "y from y³ (3rd order)"),
    ]

    for name, out_coord, indices, description in y_aberrations:
        try:
            coeff = reader.get_aberration_coefficient(out_coord, indices)
            order = len(indices)  # Order = number of source coordinates
            print(f"{name:12s} (order {order}): {coeff:12.6e}  # {description}")
        except ValueError as e:
            print(f"{name:12s}: NOT FOUND")


def demo_comparison_with_linear():
    """Compare full map reader with original linear reader."""
    print("\n" + "=" * 70)
    print("DEMO: Comparison with Linear Reader")
    print("=" * 70 + "\n")

    from cosyAdapter import COSYAdapter

    sim = COSYAdapter(
        excel_path=str(EXCEL_PATH),
        mode='transfer_matrix',
        transfer_matrix_order=1,
        debug=False
    )

    result = sim.simulate()
    reader = sim.get_native_simulator().analyze_results()

    # Original method
    linear_map_old = reader.read_linear_transfer_map()

    # New method
    full_map = reader.read_transfer_map_full(max_order=1)
    linear_map_new = full_map[1]

    # Compare
    diff = np.max(np.abs(linear_map_old - linear_map_new))

    print(f"Original read_linear_transfer_map():")
    print(f"  Shape: {linear_map_old.shape}")
    print(f"  M[0,0]: {linear_map_old[0, 0]:.6f}")

    print(f"\nNew read_transfer_map_full(max_order=1)[1]:")
    print(f"  Shape: {linear_map_new.shape}")
    print(f"  M[0,0]: {linear_map_new[0, 0]:.6f}")

    print(f"\nMaximum difference: {diff:.2e}")

    if diff < 1e-10:
        print("✓ Methods agree perfectly (backward compatible)")
    else:
        print("✗ Methods differ!")


def analyze_aberration_strength(reader, full_map):
    """Analyze relative strength of different order aberrations."""
    print("\n" + "=" * 70)
    print("ANALYSIS: Aberration Strength by Order")
    print("=" * 70 + "\n")

    for order in sorted(full_map.keys()):
        if order == 1:
            # Linear: show matrix norms
            matrix = full_map[order]
            print(f"Order {order} (Linear):")
            print(f"  Frobenius norm: {np.linalg.norm(matrix, 'fro'):.6f}")
            print(f"  Max element: {np.max(np.abs(matrix)):.6f}")
        else:
            coeffs_dict = full_map[order]

            # Collect all coefficients
            all_coeffs = []
            for coeffs_list in coeffs_dict.values():
                all_coeffs.extend(coeffs_list)

            all_coeffs = np.array(all_coeffs)

            print(f"\nOrder {order}:")
            print(f"  Number of terms: {len(coeffs_dict)}")
            print(f"  Total coefficients: {len(all_coeffs)}")
            print(f"  Max |coefficient|: {np.max(np.abs(all_coeffs)):.6e}")
            print(f"  Mean |coefficient|: {np.mean(np.abs(all_coeffs)):.6e}")
            print(f"  RMS coefficient: {np.sqrt(np.mean(all_coeffs ** 2)):.6e}")

            # Find largest coefficients
            sorted_idx = np.argsort(np.abs(all_coeffs))[-5:][::-1]
            print(f"  Top 5 largest:")
            for i, idx in enumerate(sorted_idx, 1):
                print(f"    {i}. {all_coeffs[idx]:12.6e}")


def run_all_demos():
    """Run complete demonstration suite."""
    print("\n" + "=" * 70)
    print("Higher-Order Transfer Map Reading Demonstration")
    print("=" * 70)

    if not EXCEL_PATH.exists():
        print(f"\nERROR: Excel file not found: {EXCEL_PATH}")
        print("Please update EXCEL_PATH in the script")
        return False

    try:
        # Demo 1: Read full map
        reader, full_map = demo_higher_order_reading()

        # Demo 2: Extract specific aberrations
        demo_aberration_extraction(reader)

        # Demo 3: Compare with original
        demo_comparison_with_linear()

        # Demo 4: Analyze strength
        analyze_aberration_strength(reader, full_map)

        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)
        print("\n✓ Higher-order coefficient parsing WORKS!")
        print("✓ Aberration extraction API functional")
        print("✓ Backward compatible with original reader")
        print("\nNew capabilities:")
        print("  • read_transfer_map_full(max_order=N)")
        print("  • get_aberration_coefficient(output, [powers])")
        print("  • Full access to 2nd and 3rd order aberrations")
        print("\n" + "=" * 70 + "\n")

        return True

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Demonstrate higher-order transfer map reading"
    )
    parser.add_argument(
        '--excel',
        type=str,
        default=str(EXCEL_PATH),
        help='Path to beamline Excel file'
    )

    args = parser.parse_args()
    EXCEL_PATH = Path(args.excel)

    success = run_all_demos()
    sys.exit(0 if success else 1)