"""C4: Cross-code benchmark — FELsim regression + RF-Track comparison.

Tier 1 (CI-friendly, no external deps):
  - Frozen Twiss evolution at key checkpoints (seed=42, 500 particles, 40 MeV)
  - Emittance conservation (y-plane; x-plane has apparent growth in
    dispersive regions due to x-δ coupling, which is expected physics)

Tier 2 (requires RF-Track):
  - FELsim vs RF-Track Twiss comparison in drift/quad-only regions
  - RMS envelope comparison over full beamline
  - y-emittance conservation in both codes

Author: Eremey Valetov
"""

import sys
from pathlib import Path

import numpy as np
import pytest

_BACKEND = Path(__file__).resolve().parent.parent
_PROJECT_ROOT = _BACKEND.parent
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))

JSON_PATH = _PROJECT_ROOT / "var" / "UH_FEL_beamline.json"

try:
    import RF_Track
    _RFTRACK_AVAILABLE = True
except ImportError:
    _RFTRACK_AVAILABLE = False


# ── Helpers ──────────────────────────────────────────────────────────────

def load_beamline(path):
    import latticeLoader
    return latticeLoader.create_beamline(str(path))


def compute_twiss(particles):
    """Compute Twiss from covariance matrix."""
    x, xp = particles[:, 0], particles[:, 1]
    y, yp = particles[:, 2], particles[:, 3]

    def _twiss_1d(u, up):
        u2 = np.var(u, ddof=1)
        up2 = np.var(up, ddof=1)
        uup = np.cov(u, up, ddof=1)[0, 1]
        emit = np.sqrt(max(u2 * up2 - uup**2, 0.0))
        if emit < 1e-30:
            return {'beta': 0, 'alpha': 0, 'emittance': 0}
        return {
            'beta': u2 / emit,
            'alpha': -uup / emit,
            'emittance': emit,
        }

    return {'x': _twiss_1d(x, xp), 'y': _twiss_1d(y, yp)}


def propagate_with_checkpoints(beamline, particles, checkpoints, energy_mev=40.0):
    """Track through beamline, recording Twiss at checkpoint indices."""
    for elem in beamline:
        elem.setE(energy_mev)

    results = {}
    state = particles.copy()
    s = 0.0
    for i, elem in enumerate(beamline):
        state = np.array(elem.useMatrice(state))
        s += elem.length
        if i in checkpoints:
            tw = compute_twiss(state)
            results[i] = {
                's': s,
                'beta_x': tw['x']['beta'], 'alpha_x': tw['x']['alpha'],
                'beta_y': tw['y']['beta'], 'alpha_y': tw['y']['alpha'],
                'emit_x': tw['x']['emittance'], 'emit_y': tw['y']['emittance'],
            }
    return results, state


# ── Frozen reference (seed=42, 500 particles, 40 MeV, compute_twiss) ────
#
# Note: x-emittance shows apparent growth through dispersive regions
# (chicane dipoles at elements 6-14, 23-31, 46-54, 66-74, 83-91, etc.)
# due to x-δ coupling. This is expected physics, not a bug.
# y-emittance is constant (no vertical dispersion in this beamline).

ENERGY = 40.0
N_PARTICLES = 500
SEED = 42

CHECKPOINT_INDICES = [0, 10, 20, 40, 60, 80, 87, 100, 120, 136]

FROZEN_TWISS = {
    0: {"s": 0.3587750000, "beta_x": 9.228324410989714e+00, "alpha_x": -3.334559842427499e-02, "beta_y": 9.685194665767733e+00, "alpha_y": -1.876938017195142e-02, "emit_x": 9.622477244496506e-08, "emit_y": 9.636959561757479e-08},
    10: {"s": 2.1219410000, "beta_x": 3.683009845309189e-01, "alpha_x": 2.397088135319600e+00, "beta_y": 1.639911982314170e-01, "alpha_y": -9.474945004189572e-01, "emit_x": 9.785260984839218e-08, "emit_y": 9.636959561757480e-08},
    20: {"s": 3.0884110000, "beta_x": 4.085196066620672e+00, "alpha_x": 8.736363985829305e+00, "beta_y": 5.081901308070959e+00, "alpha_y": 7.969942685472565e+00, "emit_x": 9.622285920787624e-08, "emit_y": 9.636959561757503e-08},
    40: {"s": 5.0748270000, "beta_x": 7.299971266737980e+00, "alpha_x": -1.542894356595095e+01, "beta_y": 2.506167176798369e+00, "alpha_y": 9.363019195163711e+00, "emit_x": 9.649344492814618e-08, "emit_y": 9.636959561757660e-08},
    60: {"s": 7.3256200000, "beta_x": 4.446696157707409e+00, "alpha_x": -1.523973851386648e+00, "beta_y": 9.662307763825356e-01, "alpha_y": 2.157144285064869e-01, "emit_x": 9.728646387461001e-08, "emit_y": 9.636959561757483e-08},
    80: {"s": 8.9467760000, "beta_x": 5.292482950906216e+01, "alpha_x": 1.113568110481251e+01, "beta_y": 8.710825063963560e+01, "alpha_y": -4.973337816985580e+01, "emit_x": 1.027987741515974e-07, "emit_y": 9.636959561769080e-08},
    87: {"s": 9.7177930000, "beta_x": 5.980608290395518e+00, "alpha_x": 5.401402290386918e+01, "beta_y": 3.031458560167554e+02, "alpha_y": -1.788557043126772e+03, "emit_x": 3.263157581891211e-07, "emit_y": 9.636959578302221e-08},
    100: {"s": 11.1467980000, "beta_x": 7.905022763155311e-04, "alpha_x": 5.222412006888344e-01, "beta_y": 1.962820319130078e+03, "alpha_y": 6.709290142948992e+02, "emit_x": 4.448129116433886e-07, "emit_y": 9.636959560706212e-08},
    120: {"s": 13.7414370000, "beta_x": 4.068014777593181e+04, "alpha_x": -1.595195742430243e+04, "beta_y": 1.825099409384378e+04, "alpha_y": -9.192392529920842e+03, "emit_x": 1.094451798409208e-07, "emit_y": 9.636959852971623e-08},
    136: {"s": 14.7601610000, "beta_x": 7.729867645823442e+04, "alpha_x": -2.198916541827052e+04, "beta_y": 1.018637470704576e+02, "alpha_y": -2.483880041623122e+03, "emit_x": 1.094451817757546e-07, "emit_y": 9.636959578302221e-08},
}

INIT_EMIT_X = 9.622477244496506e-08
INIT_EMIT_Y = 9.636959561757477e-08

# First dipole is at element 6 — elements 0–5 are dispersion-free
FIRST_DIPOLE_INDEX = 6


# ── Fixtures ─────────────────────────────────────────────────────────────

@pytest.fixture
def beamline():
    if not JSON_PATH.exists():
        pytest.skip("UH_FEL_beamline.json not found")
    return load_beamline(JSON_PATH)


@pytest.fixture
def particles():
    np.random.seed(SEED)
    return np.random.randn(N_PARTICLES, 6) * [1e-3, 1e-4, 1e-3, 1e-4, 1e-3, 0.005]


# ── Tier 1: FELsim regression (CI-friendly) ─────────────────────────────

class TestFELsimTwissRegression:
    """Frozen Twiss evolution: any code change that alters physics is caught."""

    def test_twiss_at_checkpoints(self, beamline, particles):
        """Twiss parameters at 10 key checkpoints must match frozen values."""
        results, _ = propagate_with_checkpoints(
            beamline, particles, set(CHECKPOINT_INDICES), ENERGY
        )

        for idx in CHECKPOINT_INDICES:
            ref = FROZEN_TWISS[idx]
            act = results[idx]

            np.testing.assert_allclose(
                act['s'], ref['s'], atol=1e-8,
                err_msg=f"s mismatch at element {idx}"
            )
            for param in ('beta_x', 'alpha_x', 'beta_y', 'alpha_y',
                          'emit_x', 'emit_y'):
                np.testing.assert_allclose(
                    act[param], ref[param], rtol=1e-10,
                    err_msg=f"{param} mismatch at element {idx}"
                )

    def test_y_emittance_conservation(self, beamline, particles):
        """y-emittance is conserved (no vertical dispersion in this beamline)."""
        results, _ = propagate_with_checkpoints(
            beamline, particles, set(CHECKPOINT_INDICES), ENERGY
        )
        for idx in CHECKPOINT_INDICES:
            np.testing.assert_allclose(
                results[idx]['emit_y'], INIT_EMIT_Y, rtol=1e-4,
                err_msg=f"y-emittance not conserved at element {idx}"
            )

    def test_x_emittance_pre_dipole(self, beamline, particles):
        """x-emittance is conserved before the first dipole (no dispersion)."""
        pre_dipole = [i for i in CHECKPOINT_INDICES if i < FIRST_DIPOLE_INDEX]
        if not pre_dipole:
            pytest.skip("No checkpoints before first dipole")

        results, _ = propagate_with_checkpoints(
            beamline, particles, set(pre_dipole), ENERGY
        )
        for idx in pre_dipole:
            np.testing.assert_allclose(
                results[idx]['emit_x'], INIT_EMIT_X, rtol=1e-4,
                err_msg=f"x-emittance not conserved at element {idx}"
            )

    def test_beamline_geometry(self, beamline):
        """Beamline length and element count must match known values."""
        total_length = sum(e.length for e in beamline)
        assert len(beamline) == 137
        np.testing.assert_allclose(total_length, 14.760161, atol=1e-5)

    def test_final_particles_deterministic(self, beamline, particles):
        """Two identical runs produce bitwise-identical results."""
        _, final1 = propagate_with_checkpoints(
            beamline, particles, set(), ENERGY
        )
        _, final2 = propagate_with_checkpoints(
            beamline, particles.copy(), set(), ENERGY
        )
        np.testing.assert_array_equal(final1, final2)


# ── Tier 2: FELsim vs RF-Track (requires RF-Track) ──────────────────────

@pytest.mark.skipif(not _RFTRACK_AVAILABLE, reason="RF-Track not installed")
class TestFELsimVsRFTrack:
    """FELsim vs RF-Track Twiss comparison.

    FELsim uses linear transfer matrices; RF-Track uses analytical
    sector-bend dipole model. In drift/quad sections the codes should
    agree closely; near dipoles, differences are expected due to
    edge-focusing and fringe-field models.
    """

    # Compare in the first 60 elements (before chicane C3 at elem ~66)
    # to avoid dipole-model differences contaminating nearby elements
    COMPARE_INDICES = sorted(set(range(0, 60, 5)) | {59})

    def test_drift_quad_agreement(self, beamline, particles):
        """In drift/quad-only regions, FELsim and RF-Track beta agree within 10%."""
        from beamline import dipole as dipole_cls, dipole_wedge as dpw_cls

        # Elements that are dipoles or within ±2 of a dipole
        dipole_neighborhood = set()
        for i, e in enumerate(beamline):
            if isinstance(e, (dipole_cls, dpw_cls)):
                for j in range(max(0, i - 2), min(len(beamline), i + 3)):
                    dipole_neighborhood.add(j)

        compare_dq = sorted(
            set(self.COMPARE_INDICES) - dipole_neighborhood
        )
        if len(compare_dq) < 3:
            pytest.skip("Too few drift/quad comparison points")

        # FELsim
        felsim_twiss, _ = propagate_with_checkpoints(
            beamline, particles, set(compare_dq), ENERGY
        )

        # RF-Track via collect_evolution
        from multiCodeSimulator import _felsim_to_generic
        from rftrackAdapter import RFTrackAdapter

        rt = RFTrackAdapter(beam_energy=ENERGY)
        generic_bl = [_felsim_to_generic(e) for e in beamline]
        rt.set_beamline(generic_bl)
        evolution = rt.collect_evolution(
            particles, checkpoint_elements=list(self.COMPARE_INDICES)
        )
        rt_df = evolution.get_twiss_evolution()

        for idx in compare_dq:
            f = felsim_twiss[idx]
            s_target = f['s']
            s_diff = np.abs(rt_df['s'].values - s_target)
            closest = s_diff.argmin()
            if s_diff[closest] > 0.01:
                continue

            rt_row = rt_df.iloc[closest]
            for param, label in [('beta_x', 'x'), ('beta_y', 'y')]:
                fval = f[param]
                rval = rt_row[param]
                if fval > 0.01 and rval > 0.01:
                    ratio = rval / fval
                    assert 0.9 < ratio < 1.1, (
                        f"elem {idx}: {param} FELsim={fval:.4g} "
                        f"RF-Track={rval:.4g} ratio={ratio:.4f}"
                    )

    def test_full_beamline_rms_envelope(self, beamline, particles):
        """Both codes produce finite output with similar RMS envelopes."""
        _, felsim_final = propagate_with_checkpoints(
            beamline, particles, set(), ENERGY
        )

        from multiCodeSimulator import _felsim_to_generic
        from rftrackAdapter import RFTrackAdapter

        rt = RFTrackAdapter(beam_energy=ENERGY)
        generic_bl = [_felsim_to_generic(e) for e in beamline]
        rt.set_beamline(generic_bl)
        rt_result = rt.simulate(particles=particles)
        assert rt_result.success

        assert np.all(np.isfinite(felsim_final))
        assert np.all(np.isfinite(rt_result.final_particles))

        # Transverse RMS within 2 orders of magnitude
        for col, name in [(0, 'x'), (2, 'y')]:
            f_rms = np.std(felsim_final[:, col])
            r_rms = np.std(rt_result.final_particles[:, col])
            if f_rms > 0 and r_rms > 0:
                ratio = f_rms / r_rms
                assert 0.01 < ratio < 100, (
                    f"{name}: FELsim RMS={f_rms:.4g} vs RF-Track RMS={r_rms:.4g}"
                )

    def test_y_emittance_both_codes(self, beamline, particles):
        """Both codes conserve y-emittance (no vertical dispersion)."""
        # FELsim at element 59 (before chicane C3)
        results, _ = propagate_with_checkpoints(
            beamline, particles, {59}, ENERGY
        )

        # RF-Track through same elements
        from multiCodeSimulator import _felsim_to_generic
        from rftrackAdapter import RFTrackAdapter

        rt = RFTrackAdapter(beam_energy=ENERGY)
        generic_bl = [_felsim_to_generic(e) for e in beamline[:60]]
        rt.set_beamline(generic_bl)
        rt_result = rt.simulate(particles=particles)
        assert rt_result.success

        rt_twiss = compute_twiss(rt_result.final_particles)

        # Both should conserve y-emittance
        np.testing.assert_allclose(
            results[59]['emit_y'], INIT_EMIT_Y, rtol=0.01,
            err_msg="FELsim y-emittance not conserved at elem 59"
        )
        np.testing.assert_allclose(
            rt_twiss['y']['emittance'], INIT_EMIT_Y, rtol=0.01,
            err_msg="RF-Track y-emittance not conserved at elem 59"
        )
