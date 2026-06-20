"""Unit tests for the common SC distribution generator + manifest."""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

import common_distribution as cd
from common_distribution import (
    BeamManifest, build, make_felsim_distribution, array_sha256,
    write_manifest, regenerate_from_manifest, to_physical_si, gamma_beta,
)


def test_shape_and_finite():
    man = BeamManifest(n_p=2000)
    arr = make_felsim_distribution(man)
    assert arr.shape == (2000, 6)
    assert np.all(np.isfinite(arr))


def test_determinism_same_seed():
    a = make_felsim_distribution(BeamManifest(n_p=5000, seed=123))
    b = make_felsim_distribution(BeamManifest(n_p=5000, seed=123))
    assert array_sha256(a) == array_sha256(b)
    assert np.array_equal(a, b)


def test_seed_changes_array():
    a = make_felsim_distribution(BeamManifest(n_p=5000, seed=123))
    b = make_felsim_distribution(BeamManifest(n_p=5000, seed=124))
    assert array_sha256(a) != array_sha256(b)


def test_manifest_roundtrip_reproducible():
    man = BeamManifest(n_p=4000, seed=777, label="roundtrip")
    arr, man = build(man)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "m.json"
        write_manifest(man, p)
        arr2 = regenerate_from_manifest(p)
    assert np.array_equal(arr, arr2)


def test_emittance_recovered():
    """Generated beam should recover the requested normalized emittance and
    matched Twiss to sampling precision."""
    man = BeamManifest(n_p=200000, seed=42, eps_n_mm_mrad=8.0,
                       betx_m=3.0, alfx=0.0, bety_m=5.0, alfy=0.0)
    arr = make_felsim_distribution(man)
    g, b = gamma_beta(man.energy_mev)
    bg = b * g
    # geometric emittance in mm.mrad from second moments
    x, xp, y, yp = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
    eps_x = np.sqrt(np.mean(x**2) * np.mean(xp**2) - np.mean(x * xp)**2)
    eps_y = np.sqrt(np.mean(y**2) * np.mean(yp**2) - np.mean(y * yp)**2)
    eps_nx = bg * eps_x  # mm.mrad (x in mm, xp in mrad -> eps in mm.mrad)
    eps_ny = bg * eps_y
    assert eps_nx == pytest.approx(8.0, rel=0.02)
    assert eps_ny == pytest.approx(8.0, rel=0.02)
    # beta = <x^2>/eps_geom ; <x^2>[mm^2], eps_geom[mm.mrad] -> mm/mrad = m
    betx = np.mean(x**2) / eps_x
    bety = np.mean(y**2) / eps_y
    assert betx == pytest.approx(3.0, rel=0.02)
    assert bety == pytest.approx(5.0, rel=0.02)


def test_to_physical_si_consistency():
    man = BeamManifest(n_p=3000, seed=9)
    arr = make_felsim_distribution(man)
    si = to_physical_si(arr, man.energy_mev)
    # transverse: mm -> m
    assert np.allclose(si["x"], arr[:, 0] * 1e-3)
    assert np.allclose(si["xp"], arr[:, 1] * 1e-3)
    # longitudinal round-trips back to the same col4
    g, b = gamma_beta(man.energy_mev)
    v0 = b * cd.C
    col4_back = si["z"] / (v0 * cd.T_RF) * 1e3
    assert np.allclose(col4_back, arr[:, 4])
    # sig_z recovered
    assert si["z"].std() == pytest.approx(man.sig_z_m, rel=0.05)


def test_versions_collected():
    _, man = build(BeamManifest(n_p=500))
    assert "felsim_git" in man.versions
    assert "numpy" in man.versions
    assert man.array_sha256 is not None


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
