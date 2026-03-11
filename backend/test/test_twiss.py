"""
C4 V&V: Unit tests for Twiss parameter computation (ebeam.cal_twiss).

Tests against known analytical distributions where emittance, α, β are
predetermined by construction.

Author: Eremey Valetov
"""

import sys
from pathlib import Path

import numpy as np
import pytest

_BACKEND = Path(__file__).resolve().parent.parent
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))

from ebeam import beam


# ── Helpers ──────────────────────────────────────────────────────────────

def make_twiss_beam(beta, alpha, epsilon, n_particles=50000, seed=42):
    """Generate a 2D phase-space distribution (x, x') with known Twiss parameters.

    Uses the Courant-Snyder parameterization:
        x  = √(ε β) cos(φ) + ...
        x' = -√(ε/β) [α cos(φ) + sin(φ)] + ...
    Implemented via the σ-matrix factorization.
    """
    rng = np.random.default_rng(seed)
    gamma_tw = (1 + alpha**2) / beta
    # Sigma matrix: [[β, -α], [-α, γ]] × ε
    sigma = epsilon * np.array([[beta, -alpha], [-alpha, gamma_tw]])
    L = np.linalg.cholesky(sigma)
    uncorrelated = rng.standard_normal((n_particles, 2))
    return (L @ uncorrelated.T).T  # (N, 2)


def make_6d_beam(beta_x, alpha_x, eps_x, beta_y, alpha_y, eps_y,
                 sigma_z=1.0, sigma_delta=0.005, n_particles=50000, seed=42):
    """Generate a full 6D beam with known transverse Twiss and zero dispersion."""
    rng = np.random.default_rng(seed)
    xy_x = make_twiss_beam(beta_x, alpha_x, eps_x, n_particles, seed=seed)
    xy_y = make_twiss_beam(beta_y, alpha_y, eps_y, n_particles, seed=seed + 1)
    z = rng.normal(0, sigma_z, n_particles)
    delta = rng.normal(0, sigma_delta, n_particles)
    return np.column_stack([xy_x, xy_y, z, delta])


# ── Twiss computation tests ──────────────────────────────────────────────

class TestCalTwiss:
    @pytest.fixture
    def eb(self):
        return beam()

    def test_round_beam(self, eb):
        """Symmetric beam: β_x = β_y, α = 0, known ε."""
        beta, alpha, eps = 1.0, 0.0, 1e-6
        dist = make_6d_beam(beta, alpha, eps, beta, alpha, eps, n_particles=100000)
        _, _, twiss = eb.cal_twiss(dist)

        eps_x = twiss.loc['x'].iloc[0]
        beta_x = twiss.loc['x'].iloc[2]
        alpha_x = twiss.loc['x'].iloc[1]

        assert abs(eps_x - eps) / eps < 0.02  # 2% tolerance (statistical)
        assert abs(beta_x - beta) / beta < 0.02
        assert abs(alpha_x - alpha) < 0.02

    def test_asymmetric_twiss(self, eb):
        """Asymmetric beam with α ≠ 0."""
        beta_x, alpha_x, eps_x = 2.0, -0.5, 5e-7
        beta_y, alpha_y, eps_y = 0.5, 1.2, 1e-6
        dist = make_6d_beam(beta_x, alpha_x, eps_x, beta_y, alpha_y, eps_y,
                            n_particles=100000)
        _, _, twiss = eb.cal_twiss(dist)

        assert abs(twiss.loc['x'].iloc[0] - eps_x) / eps_x < 0.03
        assert abs(twiss.loc['x'].iloc[2] - beta_x) / beta_x < 0.03
        assert abs(twiss.loc['x'].iloc[1] - alpha_x) < 0.05

        assert abs(twiss.loc['y'].iloc[0] - eps_y) / eps_y < 0.03
        assert abs(twiss.loc['y'].iloc[2] - beta_y) / beta_y < 0.03
        assert abs(twiss.loc['y'].iloc[1] - alpha_y) < 0.05

    def test_gamma_relation(self, eb):
        """γ_Twiss = (1 + α²) / β."""
        beta_x, alpha_x, eps_x = 1.5, 0.8, 1e-6
        dist = make_6d_beam(beta_x, alpha_x, eps_x, 1.0, 0.0, 1e-6,
                            n_particles=100000)
        _, _, twiss = eb.cal_twiss(dist)

        beta_meas = twiss.loc['x'].iloc[2]
        alpha_meas = twiss.loc['x'].iloc[1]
        gamma_meas = twiss.loc['x'].iloc[3]
        gamma_expected = (1 + alpha_meas**2) / beta_meas
        assert abs(gamma_meas - gamma_expected) / gamma_expected < 1e-10

    def test_dispersion_extraction(self, eb):
        """Add known dispersion D to the beam; verify cal_twiss extracts it."""
        beta_x, alpha_x, eps_x = 1.0, 0.0, 1e-6
        sigma_delta = 0.005
        D_inject = 0.5  # m
        Dp_inject = 0.1

        dist = make_6d_beam(beta_x, alpha_x, eps_x, 1.0, 0.0, 1e-6,
                            sigma_delta=sigma_delta, n_particles=100000)
        # Add dispersion: x += D·δ, x' += D'·δ
        dist[:, 0] += D_inject * dist[:, 5]
        dist[:, 1] += Dp_inject * dist[:, 5]

        _, _, twiss = eb.cal_twiss(dist)
        D_meas = twiss.loc['x'].iloc[4]
        Dp_meas = twiss.loc['x'].iloc[5]

        assert abs(D_meas - D_inject) / D_inject < 0.05
        assert abs(Dp_meas - Dp_inject) / Dp_inject < 0.1

    def test_emittance_invariant_under_dispersion(self, eb):
        """Intrinsic emittance should not change when dispersion is added."""
        beta_x, alpha_x, eps_x = 1.0, 0.0, 1e-6
        sigma_delta = 0.005

        dist_no_D = make_6d_beam(beta_x, alpha_x, eps_x, 1.0, 0.0, 1e-6,
                                 sigma_delta=sigma_delta, n_particles=100000)
        _, _, twiss_no_D = eb.cal_twiss(dist_no_D)

        dist_with_D = dist_no_D.copy()
        dist_with_D[:, 0] += 0.5 * dist_with_D[:, 5]
        dist_with_D[:, 1] += 0.1 * dist_with_D[:, 5]
        _, _, twiss_with_D = eb.cal_twiss(dist_with_D)

        eps_no_D = twiss_no_D.loc['x'].iloc[0]
        eps_with_D = twiss_with_D.loc['x'].iloc[0]
        assert abs(eps_no_D - eps_with_D) / eps_no_D < 0.02

    def test_zero_energy_spread_no_dispersion(self, eb):
        """With σ_δ=0, dispersion should be zero (or handled gracefully)."""
        dist = make_6d_beam(1.0, 0.0, 1e-6, 1.0, 0.0, 1e-6,
                            sigma_delta=0.0, n_particles=10000)
        # This may produce NaN/Inf in dispersion calculation (division by σ_δ=0)
        _, _, twiss = eb.cal_twiss(dist)
        # The test verifies this doesn't crash; actual behavior may produce NaN
        # which is documented as an edge case
        assert twiss.shape == (3, 7)

    def test_two_particles_handles_gracefully(self, eb):
        """Two particles is the minimum for ddof=1 covariance."""
        dist = np.array([[0.001, 0.0001, 0.002, -0.0001, 0, 0.005],
                          [-0.001, -0.0001, -0.002, 0.0001, 0, -0.005]])
        _, _, twiss = eb.cal_twiss(dist)
        assert twiss.shape == (3, 7)

    def test_longitudinal_emittance(self, eb):
        """Longitudinal emittance should be computed from z-δ plane."""
        sigma_z, sigma_delta = 0.5, 0.01
        rng = np.random.default_rng(42)
        n = 100000
        dist = np.column_stack([
            rng.normal(0, 1e-3, n), rng.normal(0, 1e-4, n),
            rng.normal(0, 1e-3, n), rng.normal(0, 1e-4, n),
            rng.normal(0, sigma_z, n), rng.normal(0, sigma_delta, n)
        ])
        _, _, twiss = eb.cal_twiss(dist)
        eps_z = twiss.loc['z'].iloc[0]
        expected = sigma_z * sigma_delta  # uncorrelated → ε = σ_z × σ_δ
        assert abs(eps_z - expected) / expected < 0.02


# ── Twiss propagation through elements ───────────────────────────────────

class TestTwissPropagation:
    """Verify that tracking particles through an element and recomputing
    Twiss gives results consistent with the matrix transformation."""

    @pytest.fixture
    def eb(self):
        return beam()

    def test_drift_preserves_emittance(self, eb):
        dist = make_6d_beam(1.0, 0.0, 1e-6, 1.0, 0.0, 1e-6, n_particles=50000)
        _, _, twiss_before = eb.cal_twiss(dist)

        from beamline import driftLattice
        d = driftLattice(2.0)
        tracked = np.array(d.useMatrice(dist))
        _, _, twiss_after = eb.cal_twiss(tracked)

        eps_before = twiss_before.loc['x'].iloc[0]
        eps_after = twiss_after.loc['x'].iloc[0]
        assert abs(eps_before - eps_after) / eps_before < 0.01

    def test_quad_preserves_emittance(self, eb):
        dist = make_6d_beam(1.0, 0.0, 1e-6, 1.0, 0.0, 1e-6, n_particles=50000)
        _, _, twiss_before = eb.cal_twiss(dist)

        from beamline import qpfLattice
        q = qpfLattice(current=3.0, length=0.0889)
        tracked = np.array(q.useMatrice(dist))
        _, _, twiss_after = eb.cal_twiss(tracked)

        eps_before = twiss_before.loc['x'].iloc[0]
        eps_after = twiss_after.loc['x'].iloc[0]
        assert abs(eps_before - eps_after) / eps_before < 0.01

    def test_drift_beta_evolution(self, eb):
        """Through a drift, β(s) = β₀ - 2α₀s + γ₀s² (analytical formula)."""
        beta_0, alpha_0 = 1.0, 0.0
        eps = 1e-6
        gamma_0 = (1 + alpha_0**2) / beta_0
        s = 0.5  # drift length

        # Need nonzero sigma_delta to avoid NaN in dispersion extraction
        dist = make_6d_beam(beta_0, alpha_0, eps, 1.0, 0.0, 1e-6,
                            sigma_delta=0.001, n_particles=200000)

        from beamline import driftLattice
        d = driftLattice(s)
        tracked = np.array(d.useMatrice(dist))
        _, _, twiss = eb.cal_twiss(tracked)

        beta_expected = beta_0 - 2 * alpha_0 * s + gamma_0 * s**2
        beta_meas = twiss.loc['x'].iloc[2]
        assert abs(beta_meas - beta_expected) / beta_expected < 0.02
