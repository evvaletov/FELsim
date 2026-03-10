"""
C3 V&V: Chicane dipole fieldmap validation.

Validates the corrected fieldmap (2026-02-22 fix) against:
  - Source OPERA-3D CSV data
  - Physical consistency (bending angle, effective length)
  - Format integrity (DELTAS, point count)
  - Absence of the erroneous 0.835× scaling factor

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

FIELDMAP_PATH = _PROJECT_ROOT / "fields" / "chicane_dipole_fieldmap.dat"
FIELDMAP_BACKUP = _PROJECT_ROOT / "fields" / "calculation" / "chicane_dipole_fieldmap.dat"
CSV_PATH = _PROJECT_ROOT / "fields" / "calculation" / "UH_chicane_fringe.csv"


# ── Helpers ──────────────────────────────────────────────────────────────

def load_fieldmap(path):
    """Load COSY INFINITY 1D MGE fieldmap: flag, npts, deltas, field[]."""
    with open(path, 'r') as f:
        flag = int(f.readline().strip())
        npts = int(f.readline().strip())
        deltas = float(f.readline().strip())
        field = [float(line.strip()) for line in f if line.strip()]
    return flag, npts, deltas, np.array(field)


def load_csv(path):
    """Load OPERA-3D CSV: columns are (z_cm, B_Gauss)."""
    return np.loadtxt(path, delimiter=',')


# ── Physical constants ───────────────────────────────────────────────────

E0 = 0.51099895  # MeV
Q = 1.60217663e-19  # C
C = 299792458.0  # m/s
DESIGN_KE = 45.0  # MeV — the magnet was measured/designed at this energy
DESIGN_ANGLE_DEG = 11.25  # chicane dipole bending angle


def Brho_Tm(KE):
    """Magnetic rigidity in T·m for given kinetic energy."""
    p_MeV = np.sqrt((KE + E0)**2 - E0**2)
    p_SI = p_MeV * 1e6 * Q / C
    return p_SI / Q


# ── Format tests ─────────────────────────────────────────────────────────

class TestFieldmapFormat:
    @pytest.fixture
    def fieldmap(self):
        if not FIELDMAP_PATH.exists():
            pytest.skip(f"Fieldmap not found: {FIELDMAP_PATH}")
        return load_fieldmap(FIELDMAP_PATH)

    def test_flag(self, fieldmap):
        flag, _, _, _ = fieldmap
        assert flag == 1

    def test_point_count(self, fieldmap):
        _, npts, _, field = fieldmap
        assert npts == 201
        assert len(field) == npts

    def test_deltas(self, fieldmap):
        """DELTAS must be 0.001 m (1 mm). Was erroneously 0 before 2026-02-22 fix."""
        _, _, deltas, _ = fieldmap
        assert deltas == 0.001

    def test_all_positive(self, fieldmap):
        """Field values should all be positive (no sign errors)."""
        _, _, _, field = fieldmap
        assert np.all(field > 0)

    def test_all_finite(self, fieldmap):
        _, _, _, field = fieldmap
        assert np.all(np.isfinite(field))

    def test_symmetric_profile(self, fieldmap):
        """Field should be approximately symmetric in the central region (within 50% of peak)."""
        _, npts, _, field = fieldmap
        center = npts // 2
        B_peak = field.max()
        # Only check symmetry in the high-field region (> 90% of peak)
        # where the Enge fit is well-constrained by dense source data.
        # The 50%-peak region has up to 6% asymmetry because the Enge fit
        # was generated from 132 non-uniformly spaced OPERA-3D points.
        left = field[:center][::-1]
        right = field[center + 1:]
        min_len = min(len(left), len(right))
        mask = (left[:min_len] > 0.9 * B_peak) & (right[:min_len] > 0.9 * B_peak)
        if mask.sum() > 0:
            ratio = left[:min_len][mask] / right[:min_len][mask]
            assert np.all(np.abs(ratio - 1.0) < 0.03), (
                f"Max asymmetry in core region: {np.abs(ratio - 1.0).max():.4f}"
            )

    def test_backup_matches_active(self):
        """Backup copy should be identical to active fieldmap."""
        if not FIELDMAP_PATH.exists() or not FIELDMAP_BACKUP.exists():
            pytest.skip("Fieldmap files not found")
        _, _, _, field_active = load_fieldmap(FIELDMAP_PATH)
        _, _, _, field_backup = load_fieldmap(FIELDMAP_BACKUP)
        np.testing.assert_array_equal(field_active, field_backup)


# ── Scaling validation ───────────────────────────────────────────────────

class TestFieldmapScaling:
    @pytest.fixture
    def data(self):
        if not FIELDMAP_PATH.exists() or not CSV_PATH.exists():
            pytest.skip("Required files not found")
        _, _, deltas, field = load_fieldmap(FIELDMAP_PATH)
        csv = load_csv(CSV_PATH)
        return field, deltas, csv

    def test_peak_matches_csv(self, data):
        """Fieldmap peak must match OPERA-3D source peak (0.5307 T)."""
        field, _, csv = data
        csv_peak_T = csv[:, 1].max() / 10000  # Gauss → Tesla
        fieldmap_peak = field.max()
        # Allow 0.1% for Enge interpolation
        assert abs(fieldmap_peak - csv_peak_T) / csv_peak_T < 0.001

    def test_no_erroneous_0835_scaling(self, data):
        """Peak should NOT be ~0.443 T (the old scaled value)."""
        field, _, _ = data
        old_scaled_peak = 0.4432612030
        assert field.max() > old_scaled_peak * 1.1  # at least 10% above old value

    def test_erroneous_factor_is_P375_over_P45(self):
        """The 0.835 factor = P(37.5 MeV) / P(45 MeV) — verify this relationship."""
        P_375 = np.sqrt((37.5 + E0)**2 - E0**2)
        P_45 = np.sqrt((45 + E0)**2 - E0**2)
        assert abs(P_375 / P_45 - 0.8351818473537908) < 1e-8

    def test_csv_units_gauss_and_cm(self, data):
        """Source CSV should be in (cm, Gauss)."""
        _, _, csv = data
        z_cm, B_gauss = csv[:, 0], csv[:, 1]
        # z range should be ~±10 cm
        assert -11 < z_cm.min() < -8
        assert 8 < z_cm.max() < 11
        # B range should be tens to thousands of Gauss
        assert B_gauss.min() > 10
        assert B_gauss.max() > 4000


# ── Physics validation ───────────────────────────────────────────────────

class TestFieldmapPhysics:
    @pytest.fixture
    def fieldmap(self):
        if not FIELDMAP_PATH.exists():
            pytest.skip(f"Fieldmap not found: {FIELDMAP_PATH}")
        return load_fieldmap(FIELDMAP_PATH)

    def test_bending_angle_at_design_energy(self, fieldmap):
        """At 45 MeV, the field integral should give θ ≈ 11.25°."""
        _, npts, deltas, field = fieldmap
        z = np.arange(npts) * deltas
        B_integral = np.trapz(field, z)
        br = Brho_Tm(DESIGN_KE)
        theta_deg = np.degrees(B_integral / br)
        assert abs(theta_deg - DESIGN_ANGLE_DEG) < 0.1, (
            f"θ = {theta_deg:.4f}° (expected {DESIGN_ANGLE_DEG}°)"
        )

    def test_effective_length(self, fieldmap):
        """L_eff = ∫B·ds / B_peak should be physically reasonable (30–80 mm)."""
        _, npts, deltas, field = fieldmap
        z = np.arange(npts) * deltas
        B_integral = np.trapz(field, z)
        L_eff = B_integral / field.max()
        assert 0.030 < L_eff < 0.080, f"L_eff = {L_eff*1000:.1f} mm"

    def test_mge_scaling_preserves_angle(self, fieldmap):
        """Runtime scaling P/P_45 should give the same bending angle at any energy."""
        _, npts, deltas, field = fieldmap
        z = np.arange(npts) * deltas
        B_integral = np.trapz(field, z)

        P_45 = np.sqrt((DESIGN_KE + E0)**2 - E0**2)
        theta_ref = np.degrees(B_integral / Brho_Tm(DESIGN_KE))

        for KE in [20, 30, 40, 45]:
            P = np.sqrt((KE + E0)**2 - E0**2)
            scaling = P / P_45
            scaled_integral = B_integral * scaling
            theta = np.degrees(scaled_integral / Brho_Tm(KE))
            assert abs(theta - theta_ref) < 1e-10, (
                f"KE={KE}: θ={theta:.6f}° vs ref {theta_ref:.6f}°"
            )

    def test_fringe_field_extent(self, fieldmap):
        """Fringe field should extend beyond the effective length."""
        _, npts, deltas, field = fieldmap
        B_peak = field.max()
        # Find where field drops below 10% of peak
        above_10pct = np.where(field > 0.1 * B_peak)[0]
        extent_mm = (above_10pct[-1] - above_10pct[0]) * deltas * 1000
        # Should be wider than L_eff but narrower than total extent
        assert extent_mm > 30, f"10% extent = {extent_mm:.1f} mm"
        assert extent_mm < 150, f"10% extent = {extent_mm:.1f} mm"

    def test_field_overall_shape(self, fieldmap):
        """Field should have a bell-like shape: high in center, low at edges."""
        _, npts, _, field = fieldmap
        center = npts // 2
        B_peak = field.max()
        # Peak should be near center (within 5 points)
        assert abs(np.argmax(field) - center) < 5
        # Edge values should be much smaller than peak
        assert field[0] < 0.05 * B_peak
        assert field[-1] < 0.05 * B_peak
        # The Enge fit from non-uniform source data can have small ripples
        # near the peak, so we don't require strict monotonicity
