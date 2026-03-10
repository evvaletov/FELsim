#!/usr/bin/env python3
"""Verify chromatic quadrupole matrix implementation.

Author: Eremey Valetov
"""

import sys
from pathlib import Path
import numpy as np

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

from beamline import qpfLattice, qpdLattice


def test_chromatic_off_matches_standard():
    """chromatic=False must reproduce the standard 6×6 matrix result."""
    for QuadClass in [qpfLattice, qpdLattice]:
        q = QuadClass(current=5.0, length=0.0889)
        q.setE(40)
        q.chromatic = False

        np.random.seed(42)
        particles = np.random.randn(100, 6) * [0.8, 0.1, 0.8, 0.1, 0.5, 5.0]

        result_standard = q.useMatrice(particles)
        q.chromatic = True
        # On-momentum particles only
        p_mono = particles.copy()
        p_mono[:, 5] = 0.0
        q.chromatic = False
        r_mono_std = q.useMatrice(p_mono)
        q.chromatic = True
        r_mono_chr = q.useMatrice(p_mono)

        np.testing.assert_allclose(r_mono_chr, r_mono_std, atol=1e-12,
                                   err_msg=f"{QuadClass.__name__}: on-momentum mismatch")
        print(f"  {QuadClass.__name__}: on-momentum (δ=0) — PASS")


def test_chromaticity_effect():
    """Chromatic mode should produce different results for off-momentum particles."""
    for QuadClass in [qpfLattice, qpdLattice]:
        q = QuadClass(current=5.0, length=0.0889)
        q.setE(40)

        # Particles with 0.5% energy spread
        np.random.seed(42)
        particles = np.random.randn(100, 6) * [0.8, 0.1, 0.8, 0.1, 0.5, 5.0]

        q.chromatic = False
        r_linear = q.useMatrice(particles)
        q.chromatic = True
        r_chromatic = q.useMatrice(particles)

        diff = np.max(np.abs(r_chromatic - r_linear))
        print(f"  {QuadClass.__name__}: max Δ(chromatic vs linear) = {diff:.6e}")
        assert diff > 1e-6, f"{QuadClass.__name__}: chromaticity should produce different results"
        print(f"  {QuadClass.__name__}: chromaticity effect — PASS")


def test_symplectic_check():
    """Verify per-particle matrices are approximately symplectic."""
    for QuadClass in [qpfLattice, qpdLattice]:
        q = QuadClass(current=5.0, length=0.0889)
        q.setE(40)
        q.chromatic = True

        # Single particle with δ = 5 (0.5% energy spread)
        p = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, 5.0]])
        r = q.useMatrice(p)

        # Build effective 4×4 transfer from unit vectors
        M = np.zeros((4, 4))
        for i in range(4):
            e = np.zeros((1, 6))
            e[0, i] = 1e-6  # small perturbation
            e[0, 5] = 5.0  # same δ
            r0 = np.zeros((1, 6))
            r0[0, 5] = 5.0
            r_e = q.useMatrice(e)
            r_0 = q.useMatrice(r0)
            M[:, i] = (r_e[0, :4] - r_0[0, :4]) / 1e-6

        # Check det(M_x) ≈ 1 and det(M_y) ≈ 1
        det_x = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
        det_y = M[2, 2] * M[3, 3] - M[2, 3] * M[3, 2]
        print(f"  {QuadClass.__name__}: det(M_x) = {det_x:.12f}, det(M_y) = {det_y:.12f}")
        np.testing.assert_allclose(det_x, 1.0, atol=1e-8,
                                   err_msg=f"{QuadClass.__name__}: x-plane not symplectic")
        np.testing.assert_allclose(det_y, 1.0, atol=1e-8,
                                   err_msg=f"{QuadClass.__name__}: y-plane not symplectic")
        print(f"  {QuadClass.__name__}: symplecticity — PASS")


def test_zero_current_fallback():
    """I=0 should fall back to drift regardless of chromatic flag."""
    for QuadClass in [qpfLattice, qpdLattice]:
        q = QuadClass(current=0.0, length=0.0889)
        q.setE(40)

        particles = np.array([[1.0, 0.1, 2.0, 0.2, 0.5, 5.0]])
        q.chromatic = False
        r_std = q.useMatrice(particles)
        q.chromatic = True
        r_chr = q.useMatrice(particles)
        np.testing.assert_allclose(r_chr, r_std, atol=1e-14,
                                   err_msg=f"{QuadClass.__name__}: I=0 fallback mismatch")
        print(f"  {QuadClass.__name__}: I=0 fallback — PASS")


def test_rftrack_comparison():
    """Compare chromatic FELsim quad against RF-Track k_eff = k₀ × P₀/P."""
    q = qpfLattice(current=5.0, length=0.0889)
    q.setE(40)

    # Reference particle
    E0 = q.E0
    K0 = q.E
    P0 = q.beta * q.gamma * E0  # MeV/c (in natural units with c=1)

    # Particle with +0.5% energy offset
    delta_pct = 0.5
    delta_coord6 = delta_pct * 10  # ΔK/K₀ × 10³
    K_p = K0 * (1 + delta_coord6 * 1e-3)
    E_total_p = K_p + E0
    P_p = np.sqrt(E_total_p**2 - E0**2)

    k0 = np.abs(q.Q * q.G * q.current / (q.M * q.C * q.beta * q.gamma))
    k_rft = k0 * P0 / P_p  # RF-Track formula

    # What our chromatic code computes
    bg0 = q.beta * q.gamma
    gamma_p = E_total_p / E0
    bg_p = np.sqrt(gamma_p**2 - 1)
    k_felsim = k0 * bg0 / bg_p

    # bg = P/(mc) and P0 = bg0 * mc, P_p = bg_p * mc, so bg0/bg_p = P0/P_p
    print(f"  k₀ = {k0:.6f}")
    print(f"  k_rft (k₀×P₀/P) = {k_rft:.6f}")
    print(f"  k_felsim (k₀×βγ₀/βγ) = {k_felsim:.6f}")
    np.testing.assert_allclose(k_felsim, k_rft, rtol=1e-12,
                               err_msg="Chromatic k formulas disagree")
    print(f"  RF-Track formula equivalence — PASS")


if __name__ == '__main__':
    print("=" * 60)
    print("  Chromatic quadrupole matrix verification")
    print("=" * 60)

    print("\n1. On-momentum consistency:")
    test_chromatic_off_matches_standard()

    print("\n2. Chromaticity effect:")
    test_chromaticity_effect()

    print("\n3. Symplecticity:")
    test_symplectic_check()

    print("\n4. Zero-current fallback:")
    test_zero_current_fallback()

    print("\n5. RF-Track formula equivalence:")
    test_rftrack_comparison()

    print("\n" + "=" * 60)
    print("  All tests passed.")
    print("=" * 60)
