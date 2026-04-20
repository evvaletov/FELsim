"""
Integration test: load the SLAC linac JSON via the RF-Track adapter and
verify the adapter produces results consistent with the standalone script.

Eremey Valetov, 2026-04-05
"""

import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
BACKEND = REPO_ROOT / 'backend'
LATTICE_JSON = REPO_ROOT / 'var' / 'slac_linac.json'

sys.path.insert(0, str(BACKEND))

MC2_MEV = 0.510998950
K_INJECT = 1.0
EXPECTED_K_OUT = 41.47  # MeV, from standalone script (phid=0, autophased)
TOLERANCE = 0.05        # MeV


def _track_linac():
    """Load SLAC linac JSON via adapter, track one electron, return K_out."""
    from rftrackAdapter import RFTrackAdapter
    import RF_Track as rft

    adapter = RFTrackAdapter(lattice_path=str(LATTICE_JSON))
    assert len(adapter.beamline) == 1, f"Expected 1 element, got {len(adapter.beamline)}"

    native = adapter._convert_element_to_native(adapter.beamline[0])
    P_in = np.sqrt((K_INJECT + MC2_MEV)**2 - MC2_MEV**2)
    bunch = rft.Bunch6d(MC2_MEV, 1.0, -1.0,
                        np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P_in]]))

    lat = rft.Lattice()
    lat.append(native)
    bunch_out = lat.track(bunch)
    M = bunch_out.get_phase_space('%x %xp %y %yp %t %Pc')
    assert M.shape[0] > 0, "Particle lost during tracking"

    P_out = M[0, 5]
    K_out = np.sqrt(P_out**2 + MC2_MEV**2) - MC2_MEV
    return K_out


def test_slac_linac_energy_gain():
    """Adapter-path energy gain matches standalone reference (41.47 MeV)."""
    K_out = _track_linac()
    assert abs(K_out - EXPECTED_K_OUT) < TOLERANCE, (
        f"K_out = {K_out:.4f} MeV, expected {EXPECTED_K_OUT} ± {TOLERANCE} MeV"
    )


def test_slac_linac_element_type():
    """Loaded element has correct type and length."""
    from rftrackAdapter import RFTrackAdapter
    adapter = RFTrackAdapter(lattice_path=str(LATTICE_JSON))
    elem = adapter.beamline[0]
    assert elem.element_type == 'RF_CAVITY'
    assert abs(elem.length - 3.048) < 1e-6


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
