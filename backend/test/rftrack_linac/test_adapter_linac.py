"""
Integration test: load the SLAC linac JSON via the RF-Track adapter and
verify the adapter produces results consistent with the standalone RF-Track
model in linac_standalone.py (same structure, autophased, phid = 0).

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

K_INJECT = 1.0
EXPECTED_K_OUT = 41.47  # MeV, linac_standalone.py (phid=0, autophased)
TOLERANCE = 0.05        # MeV


def _k_out(adapter, particles):
    fp = adapter.simulate(particles).final_particles
    assert fp.shape[0] == particles.shape[0], "Particle lost during tracking"
    return K_INJECT * (1.0 + fp[:, 5] * 1e-3)


def _track_linac():
    """Load SLAC linac JSON via adapter, track one electron, return K_out."""
    from rftrackAdapter import RFTrackAdapter

    adapter = RFTrackAdapter(lattice_path=str(LATTICE_JSON), beam_energy=K_INJECT,
                             aperture=None)
    assert len(adapter.beamline) == 1, f"Expected 1 element, got {len(adapter.beamline)}"
    return _k_out(adapter, np.zeros((1, 6)))[0]


def test_slac_linac_energy_gain():
    """Adapter-path energy gain matches the reference of linac_standalone.py
    (41.47 MeV)."""
    K_out = _track_linac()
    assert abs(K_out - EXPECTED_K_OUT) < TOLERANCE, (
        f"K_out = {K_out:.4f} MeV, expected {EXPECTED_K_OUT} ± {TOLERANCE} MeV"
    )


def test_slac_linac_phase_does_not_depend_on_the_beam():
    """The cavity is phased on the adapter's reference particle, so a second
    simulate(), element-by-element tracking, and a beam whose first particle
    is off momentum all give the on-axis particle the same energy."""
    from rftrackAdapter import RFTrackAdapter

    adapter = RFTrackAdapter(lattice_path=str(LATTICE_JSON), beam_energy=K_INJECT,
                             aperture=None)
    first = _k_out(adapter, np.zeros((1, 6)))[0]
    assert _k_out(adapter, np.zeros((1, 6)))[0] == pytest.approx(first, rel=1e-9)
    beam = np.zeros((2, 6))
    beam[0, 4:] = [30.0, 200.0]
    assert _k_out(adapter, beam)[1] == pytest.approx(first, rel=1e-9)
    evolution = adapter.collect_evolution(np.zeros((1, 6)))
    last = evolution.particles[max(evolution.s_positions)]
    assert K_INJECT * (1.0 + last[0, 5] * 1e-3) == pytest.approx(first, rel=1e-9)


def test_slac_linac_phased_again_for_a_new_beam_energy():
    """A new beam energy is a new reference particle, and the cavity is
    phased again on it: RF-Track keeps a cavity's phase once it is set
    unless it is unset first."""
    from rftrackAdapter import RFTrackAdapter

    adapter = RFTrackAdapter(lattice_path=str(LATTICE_JSON), beam_energy=K_INJECT,
                             aperture=None)
    adapter.simulate(np.zeros((1, 6)))
    adapter.set_beam_energy(2.0)
    fresh = RFTrackAdapter(lattice_path=str(LATTICE_JSON), beam_energy=2.0, aperture=None)
    np.testing.assert_allclose(adapter.simulate(np.zeros((1, 6))).final_particles,
                               fresh.simulate(np.zeros((1, 6))).final_particles, rtol=1e-12)


def test_slac_linac_element_type():
    """Loaded element has correct type and length."""
    from rftrackAdapter import RFTrackAdapter
    adapter = RFTrackAdapter(lattice_path=str(LATTICE_JSON))
    elem = adapter.beamline[0]
    assert elem.element_type == 'RF_CAVITY'
    assert abs(elem.length - 3.048) < 1e-6


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
