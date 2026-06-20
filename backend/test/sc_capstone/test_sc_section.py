"""Smoke tests for the section-aware DA-FMM tracker.

These guard the two invariants that make the SC capstone trustworthy:
  1. SC-off linear transport reproduces FELsim exactly (the tracker's maps are
     correct), and emittance is conserved to machine precision.
  2. SC-on increases the normalized emittance, and more charge => more growth.
"""
from __future__ import annotations

import numpy as np
import pytest

import common_distribution as cd
import sc_section as scs
from nosc_handoff import LATTICE, NODIP_RANGE, run_code

E = 45.0


@pytest.fixture(scope="module")
def beam():
    man = cd.BeamManifest(energy_mev=E, n_p=4000, seed=20260619,
                          sig_delta=0.0, eps_n_mm_mrad=8.0)
    return cd.make_felsim_distribution(man)


@pytest.fixture(scope="module")
def segs():
    return scs.extract_segments(LATTICE, *NODIP_RANGE, E)


def test_extract_segments_no_dipole(segs):
    kinds = {s.kind for s in segs}
    assert kinds <= {"drift", "quad"}
    assert sum(1 for s in segs if s.kind == "quad") == 6


def test_extract_segments_rejects_dipole():
    with pytest.raises(ValueError):
        scs.extract_segments(LATTICE, 0, 20, E)   # contains dipoles/wedges


def test_scoff_matches_felsim(beam, segs):
    da = scs.track_dafmm(beam, segs, E, q_nc=1.0, ds_target=0.01, spch_on=False)
    fel = run_code("felsim", beam, NODIP_RANGE, E, space_charge=False)
    assert da["sigx_exit_mm"] == pytest.approx(fel["sigx_mm"], rel=1e-4)
    assert da["sigy_exit_mm"] == pytest.approx(fel["sigy_mm"], rel=1e-4)
    # emittance conserved (no SC, linear maps, sigma_delta=0)
    assert da["epsnx_growth"] == pytest.approx(0.0, abs=1e-9)
    assert da["epsny_growth"] == pytest.approx(0.0, abs=1e-9)


def test_scon_grows_with_charge(beam, segs):
    lo = scs.track_dafmm(beam, segs, E, q_nc=0.3, ds_target=0.02, spch_on=True)
    hi = scs.track_dafmm(beam, segs, E, q_nc=3.0, ds_target=0.02, spch_on=True)
    assert hi["epsnx_growth"] > lo["epsnx_growth"] > 0.0


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
