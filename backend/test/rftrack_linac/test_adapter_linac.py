"""
Integration test: load the SLAC linac JSON via the RF-Track adapter and
verify the adapter produces results consistent with the standalone script.

Eremey Valetov, 2026-04-05
"""

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
BACKEND = REPO_ROOT / 'backend'
LATTICE_JSON = REPO_ROOT / 'var' / 'slac_linac.json'

sys.path.insert(0, str(BACKEND))

from rftrackAdapter import RFTrackAdapter

MC2_MEV = 0.510998950
K_INJECT = 1.0


def main():
    print(f"Loading linac lattice: {LATTICE_JSON}")
    adapter = RFTrackAdapter(lattice_path=str(LATTICE_JSON))

    # Verify the loaded element
    print(f"\nLoaded {len(adapter.beamline)} element(s):")
    for i, elem in enumerate(adapter.beamline):
        print(f"  [{i}] type={elem.element_type}, length={elem.length:.4f} m, "
              f"params={elem.parameters}")

    # Build the RF-Track lattice by converting the single RFC element
    import RF_Track as rft
    native = adapter._convert_element_to_native(adapter.beamline[0])
    print(f"\nNative element: {type(native).__name__}")
    print(f"  length = {native.get_length():.4f} m")
    print(f"  frequency = {native.get_frequency()/1e9:.4f} GHz")

    # Track a single electron at K=1 MeV
    P_in = np.sqrt((K_INJECT + MC2_MEV)**2 - MC2_MEV**2)
    bunch = rft.Bunch6d(MC2_MEV, 1.0, -1.0,
                        np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P_in]]))

    lat = rft.Lattice()
    lat.append(native)
    bunch_out = lat.track(bunch)
    M = bunch_out.get_phase_space('%x %xp %y %yp %t %Pc')

    if M.shape[0] == 0:
        print("PARTICLE LOST")
        return 1

    P_out = M[0, 5]
    K_out = np.sqrt(P_out**2 + MC2_MEV**2) - MC2_MEV
    print(f"\nTracking result:")
    print(f"  Input:  K = {K_INJECT:.3f} MeV, P = {P_in:.4f} MeV/c")
    print(f"  Output: K = {K_out:.3f} MeV, P = {P_out:.4f} MeV/c")

    # Compare to standalone reference (41.47 MeV at phid=0, autophased)
    expected = 41.47
    tol = 0.05
    delta = K_out - expected
    status = "PASS" if abs(delta) < tol else "FAIL"
    print(f"\nExpected (standalone): {expected} MeV")
    print(f"Delta:                 {delta:+.3f} MeV")
    print(f"Tolerance:             {tol} MeV")
    print(f"Status:                {status}")
    return 0 if abs(delta) < tol else 1


if __name__ == '__main__':
    sys.exit(main())
