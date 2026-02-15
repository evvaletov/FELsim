"""Compare FELsim DPW-DPH-DPW sandwich with COSY consolidated DIL.

Tests two representative dipole types at FR 0 and FR 3 to quantify
the model difference and diagnose the Y-plane instability in cross-validation.

Author: Eremey Valetov
"""

import sys
from pathlib import Path
import numpy as np

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

from beamline import lattice, dipole, dipole_wedge
from cosyAdapter import COSYAdapter
from cosyOptHelper import parse_beamline_felsim_indexed

Energy = 40  # MeV


def felsim_sandwich_matrix(wedge_angle_entrance, dipole_angle, dipole_length,
                           wedge_angle_exit, wedge_length, pole_gap):
    """Compute FELsim's DPW * DPH * DPW product transfer matrix."""
    relat = lattice(1, fringeType=None)
    relat.setE(E=Energy)

    entrance = dipole_wedge(
        length=wedge_length, angle=wedge_angle_entrance,
        dipole_length=dipole_length, dipole_angle=dipole_angle,
        pole_gap=pole_gap, fringeType=None
    )
    entrance.setE(E=Energy)

    body = dipole(length=dipole_length, angle=dipole_angle, fringeType=None)
    body.setE(E=Energy)

    exit_w = dipole_wedge(
        length=wedge_length, angle=wedge_angle_exit,
        dipole_length=dipole_length, dipole_angle=dipole_angle,
        pole_gap=pole_gap, fringeType=None
    )
    exit_w.setE(E=Energy)

    M_ent = entrance._compute_numeric_matrix()
    M_body = body._compute_numeric_matrix()
    M_exit = exit_w._compute_numeric_matrix()

    return M_exit @ M_body @ M_ent, (M_ent, M_body, M_exit)


EXCEL_PATH = str(Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx')


def run_cosy_single_dipole(dipole_params, fr_order=0, use_enge=False):
    """Run COSY for a single consolidated dipole and extract the transfer matrix."""
    adapter = COSYAdapter(
        lattice_path=EXCEL_PATH, mode='transfer_matrix',
        config={'simulation': {'KE': Energy, 'order': 3, 'dimensions': 3}},
        fringe_field_order=fr_order,
        debug=False
    )
    sim = adapter.get_native_simulator()
    sim.use_enge_coeffs = use_enge

    # Build minimal beamline: just the consolidated dipole
    sim.beamline = [{
        'type': 'DIPOLE_CONSOLIDATED',
        'length': dipole_params['length'],
        'angle': dipole_params['angle'],
        'entrance_angle': dipole_params['e1'],
        'exit_angle': dipole_params['e2'],
        'pole_gap': dipole_params['pole_gap'],
        'entrance_enge_coeffs': dipole_params.get('entrance_enge'),
        'exit_enge_coeffs': dipole_params.get('exit_enge'),
        'original_elements': [{}, {}, {}],
    }]

    result = adapter.simulate()
    if not result.success:
        print(f"  COSY failed: {result.metadata}")
        return None

    reader = sim.analyze_results()
    R = reader.read_linear_transfer_map()
    return R


def compare_matrices(name, M_felsim, M_cosy, labels=('FELsim', 'COSY')):
    """Compare two 6x6 transfer matrices element by element."""
    print(f"\n  {name}: {labels[0]} vs {labels[1]}")
    print(f"  {'(i,j)':<8} {labels[0]:>12} {labels[1]:>12} {'Δ':>12} {'Δ (%)':>10}")
    print("  " + "-" * 58)

    max_diff = 0
    max_ij = (0, 0)
    plane_labels = {0: 'X', 1: "X'", 2: 'Y', 3: "Y'", 4: 'T', 5: 'δ'}

    for i in range(6):
        for j in range(6):
            f_val = M_felsim[i, j]
            c_val = M_cosy[i, j]
            diff = abs(f_val - c_val)

            if diff > max_diff:
                max_diff = diff
                max_ij = (i, j)

            if abs(f_val) > 1e-12 or abs(c_val) > 1e-12:
                if abs(f_val) > 1e-12:
                    pct = diff / abs(f_val) * 100
                else:
                    pct = float('inf')

                if diff > 1e-8:
                    print(f"  ({i+1},{j+1}){plane_labels.get(i,'?')}{plane_labels.get(j,'?')}"
                          f" {f_val:>12.6e} {c_val:>12.6e} {f_val-c_val:>+12.4e} {pct:>9.2f}%")

    print(f"  Max |Δ| = {max_diff:.4e} at ({max_ij[0]+1},{max_ij[1]+1})")
    return max_diff


def analyze_y_plane(M):
    """Check Y-plane stability of a 6x6 transfer matrix."""
    My = M[2:4, 2:4]
    half_trace = abs(np.trace(My) / 2)
    stable = half_trace <= 1.0
    return half_trace, stable


def test_dipole_type(name, dipole_params, wedge_params):
    """Run full comparison for one dipole type."""
    print(f"\n{'='*70}")
    print(f"Dipole: {name}")
    print(f"  Angle: {dipole_params['angle']}°, Length: {dipole_params['length']} m")
    print(f"  Entrance edge: {dipole_params['e1']}°, Exit edge: {dipole_params['e2']}°")
    print(f"  Pole gap: {dipole_params['pole_gap']} m, Wedge length: {wedge_params['length']} m")
    has_enge = dipole_params.get('entrance_enge') or dipole_params.get('exit_enge')
    if has_enge:
        print(f"  Enge coefficients: available")
    print(f"{'='*70}")

    # FELsim sandwich
    M_felsim, (M_ent, M_body, M_exit) = felsim_sandwich_matrix(
        wedge_angle_entrance=dipole_params['e1'],
        dipole_angle=dipole_params['angle'],
        dipole_length=dipole_params['length'],
        wedge_angle_exit=dipole_params['e2'],
        wedge_length=wedge_params['length'],
        pole_gap=dipole_params['pole_gap']
    )

    print(f"\n  FELsim sandwich matrix (4x4 transverse block):")
    for i in range(4):
        row = " ".join(f"{M_felsim[i,j]:+12.6e}" for j in range(4))
        print(f"    [{row}]")

    ht_f, stable_f = analyze_y_plane(M_felsim)
    print(f"  Y-plane |Tr/2| = {ht_f:.6f} ({'stable' if stable_f else 'UNSTABLE'})")

    # COSY FR 0
    print("\n  Running COSY FR 0 (no fringe)...")
    M_cosy_fr0 = run_cosy_single_dipole(dipole_params, fr_order=0, use_enge=False)
    if M_cosy_fr0 is not None:
        ht_c0, stable_c0 = analyze_y_plane(M_cosy_fr0)
        print(f"  COSY FR 0 Y-plane |Tr/2| = {ht_c0:.6f} ({'stable' if stable_c0 else 'UNSTABLE'})")
        compare_matrices("FR 0 vs FELsim", M_felsim, M_cosy_fr0, ('FELsim', 'COSY FR0'))

    # COSY FR 3
    print("\n  Running COSY FR 3 (third-order fringe)...")
    M_cosy_fr3 = run_cosy_single_dipole(dipole_params, fr_order=3, use_enge=False)
    if M_cosy_fr3 is not None:
        ht_c3, stable_c3 = analyze_y_plane(M_cosy_fr3)
        print(f"  COSY FR 3 Y-plane |Tr/2| = {ht_c3:.6f} ({'stable' if stable_c3 else 'UNSTABLE'})")
        compare_matrices("FR 3 vs FELsim", M_felsim, M_cosy_fr3, ('FELsim', 'COSY FR3'))

    if M_cosy_fr0 is not None and M_cosy_fr3 is not None:
        compare_matrices("FR 0 vs FR 3", M_cosy_fr0, M_cosy_fr3, ('COSY FR0', 'COSY FR3'))

    # COSY FR 3 + FC (Enge) if available
    if has_enge:
        print("\n  Running COSY FR 3 + FC (Enge coefficients)...")
        M_cosy_fc = run_cosy_single_dipole(dipole_params, fr_order=3, use_enge=True)
        if M_cosy_fc is not None:
            ht_fc, stable_fc = analyze_y_plane(M_cosy_fc)
            print(f"  COSY FR3+FC Y-plane |Tr/2| = {ht_fc:.6f} ({'stable' if stable_fc else 'UNSTABLE'})")
            compare_matrices("FR3+FC vs FELsim", M_felsim, M_cosy_fc, ('FELsim', 'COSY FR3+FC'))

    return {
        'felsim': M_felsim,
        'cosy_fr0': M_cosy_fr0,
        'cosy_fr3': M_cosy_fr3,
    }


def test_full_beamline_product():
    """Compute product of all 14 dipole sandwiches and compare total Y-plane effect."""
    print(f"\n{'='*70}")
    print("Full beamline: product of all 14 dipole sandwich matrices")
    print(f"{'='*70}")

    file_path = Path(__file__).resolve().parent.parent.parent / 'beam_excel' / 'Beamline_elements.xlsx'

    from excelElements import ExcelElements
    excel = ExcelElements(str(file_path))
    beamlineUH = excel.create_beamline()
    relat = lattice(1, fringeType=None)
    relat.setE(E=Energy)
    line = relat.changeBeamType("electron", Energy, beamlineUH)[:118]

    # Find all DPW-DPH-DPW triplets
    M_all_felsim = np.eye(6)
    triplet_count = 0
    i = 0
    while i < len(line):
        elem = line[i]
        if isinstance(elem, dipole_wedge) and i + 2 < len(line):
            if isinstance(line[i + 1], dipole) and isinstance(line[i + 2], dipole_wedge):
                M_ent = line[i]._compute_numeric_matrix()
                M_body = line[i + 1]._compute_numeric_matrix()
                M_exit = line[i + 2]._compute_numeric_matrix()
                M_triplet = M_exit @ M_body @ M_ent
                M_all_felsim = M_triplet @ M_all_felsim
                triplet_count += 1
                i += 3
                continue
        i += 1

    ht, stable = analyze_y_plane(M_all_felsim)
    print(f"\n  {triplet_count} dipole triplets found")
    print(f"  Product Y-plane |Tr/2| = {ht:.6f} ({'stable' if stable else 'UNSTABLE'})")
    print(f"  Y-plane submatrix:")
    My = M_all_felsim[2:4, 2:4]
    for i in range(2):
        print(f"    [{My[i,0]:+12.6e}  {My[i,1]:+12.6e}]")


if __name__ == "__main__":
    # Type 1: Transport dipole (1.5°, no Enge)
    # First triplet: entrance edge = 0°, exit edge = 1.5°
    transport_params = {
        'length': 0.0889, 'angle': 1.5,
        'e1': 0.0, 'e2': 1.5,
        'pole_gap': 0.014478,
    }
    transport_wedge = {'length': 0.01}

    # Type 2: Chicane dipole (4°, no Enge)
    chicane_params = {
        'length': 0.04064, 'angle': 4.0,
        'e1': 2.018, 'e2': 2.018,
        'pole_gap': 0.014478,
    }
    chicane_wedge = {'length': 0.01}

    # Type 3: FC dipole (11.25°, with Enge coefficients)
    fc_enge = [56.49, -50.79, 19.32, -3.621, 0.3315, -0.01193]
    fc_params = {
        'length': 0.037389, 'angle': 11.25,
        'e1': 0.0, 'e2': 11.25,
        'pole_gap': 0.012700,
        'entrance_enge': fc_enge,
        'exit_enge': None,
    }
    fc_wedge = {'length': 0.01}

    results = {}
    results['transport'] = test_dipole_type(
        "Transport (1.5°, e1=0°, e2=1.5°)", transport_params, transport_wedge)
    results['chicane'] = test_dipole_type(
        "Chicane (-4°, symmetric edges 2.018°)", chicane_params, chicane_wedge)
    results['fc1'] = test_dipole_type(
        "FC1 (11.25°, e1=0°, e2=11.25°, with Enge)", fc_params, fc_wedge)

    # Also test the 0° edge case (where FELsim gives nonzero kick but COSY gives zero)
    zero_edge_params = {
        'length': 0.0889, 'angle': 1.5,
        'e1': 1.5, 'e2': 0.0,  # reversed: entrance has angle, exit has 0
        'pole_gap': 0.014478,
    }
    results['zero_exit'] = test_dipole_type(
        "Transport (1.5°, e1=1.5°, e2=0° — zero exit edge)", zero_edge_params, transport_wedge)

    test_full_beamline_product()

    print("\n\nSummary")
    print("=" * 70)
    for name, res in results.items():
        print(f"\n{name}:")
        for key in ['felsim', 'cosy_fr0', 'cosy_fr3']:
            M = res.get(key)
            if M is not None:
                ht, stable = analyze_y_plane(M)
                My34 = M[3, 2]
                print(f"  {key:>12}: Y-kick M(4,3) = {My34:+.6e}, |Tr_y/2| = {ht:.6f}")
