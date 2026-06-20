"""
Common initial-distribution generator + reproducibility manifest for the
multi-code space-charge capstone (cosy-fmm DA-FMM / RF-Track / xsuite).

The whole point of the capstone is that every engine starts from the *same*
macroparticles, the same charge normalization, and the same recorded provenance,
so any downstream disagreement is physics (or a known code gap), never the seed.

This module produces one (N, 6) array in FELsim coordinates

    [x(mm), x'(mrad), y(mm), y'(mrad), DToF/T_RF*1e3, DK/K0*1e3]

(matched correlated Gaussian per transverse plane, uncorrelated Gaussian
longitudinally) plus a JSON manifest capturing every parameter, the code
versions, repo git hashes, and the SHA-256 of the generated array. The adapters
(FELsim / COSY / RF-Track / xsuite) each transform this FELsim array into their
own native coordinates internally, so feeding all of them the same array is
sufficient for an apples-to-apples comparison.

Author: Eremey Valetov
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from physicalConstants import PhysicalConstants

E0_E = PhysicalConstants.E0_electron           # MeV
C = float(PhysicalConstants.C)                 # m/s
F_RF = PhysicalConstants.f_RF_default          # Hz
T_RF = 1.0 / F_RF                              # s

FELSIM_REPO = Path(__file__).resolve().parents[3]          # .../FELsim
COSY_FMM_REPO = Path("/home/evaletov/COSY/cosy-fmm")


def gamma_beta(energy_mev: float) -> tuple[float, float]:
    gamma = 1.0 + energy_mev / E0_E
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    return gamma, beta


@dataclass
class BeamManifest:
    """Everything needed to regenerate the distribution bit-for-bit, plus the
    physical beam definition the capstone benchmarks."""
    # Beam
    energy_mev: float = 45.0
    q_nc: float = 1.0
    n_p: int = 6000
    seed: int = 20260619
    # Transverse Twiss (matched launch optics) + normalized emittance
    betx_m: float = 3.0
    alfx: float = 0.0
    bety_m: float = 3.0
    alfy: float = 0.0
    eps_n_mm_mrad: float = 8.0     # normalized geometric emittance, per plane
    # Longitudinal (uncorrelated Gaussian)
    sig_z_m: float = 6.0e-4        # rms bunch length (path length), ~2 ps at c
    sig_delta: float = 5.0e-3      # rms DK/K0
    # Provenance (filled by collect_versions / build)
    label: str = "capstone_default"
    versions: dict = field(default_factory=dict)
    array_sha256: Optional[str] = None

    def geom_eps_si(self) -> float:
        """Geometric emittance [m.rad] from the normalized value."""
        gamma, beta = gamma_beta(self.energy_mev)
        return self.eps_n_mm_mrad * 1e-6 / (beta * gamma)


def _git_hash(repo: Path) -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        h = out.stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(repo), "status", "--porcelain"],
            capture_output=True, text=True, timeout=10,
        ).stdout.strip()
        return f"{h}{'+dirty' if dirty else ''}" if h else "unknown"
    except Exception:
        return "unknown"


def collect_versions() -> dict:
    v = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "felsim_git": _git_hash(FELSIM_REPO),
        "cosy_fmm_git": _git_hash(COSY_FMM_REPO),
    }
    for mod in ("RF_Track", "xtrack", "xfields", "xpart"):
        try:
            m = __import__(mod)
            v[mod] = getattr(m, "__version__", "unknown")
        except Exception:
            v[mod] = "absent"
    return v


def array_sha256(arr: np.ndarray) -> str:
    """Stable content hash of the distribution (C-contiguous float64 bytes)."""
    a = np.ascontiguousarray(arr, dtype=np.float64)
    return hashlib.sha256(a.tobytes()).hexdigest()


def _matched_plane(rng, n_p, eps_si, beta_m, alpha):
    """Correlated Gaussian (x[m], x'[rad]) matched to (beta, alpha, eps)."""
    z1 = rng.standard_normal(n_p)
    z2 = rng.standard_normal(n_p)
    x = np.sqrt(eps_si * beta_m) * z1
    xp = np.sqrt(eps_si / beta_m) * (-alpha * z1 + z2)
    return x, xp


def make_felsim_distribution(man: BeamManifest) -> np.ndarray:
    """Generate the (N, 6) FELsim-coordinate macroparticle array.

    Sampling order is fixed (x-plane, y-plane, z, delta) so a given manifest
    always yields the identical array.
    """
    rng = np.random.default_rng(man.seed)
    eps_si = man.geom_eps_si()

    x_m, xp_r = _matched_plane(rng, man.n_p, eps_si, man.betx_m, man.alfx)
    y_m, yp_r = _matched_plane(rng, man.n_p, eps_si, man.bety_m, man.alfy)
    z_m = man.sig_z_m * rng.standard_normal(man.n_p)
    delta = man.sig_delta * rng.standard_normal(man.n_p)

    gamma, beta = gamma_beta(man.energy_mev)
    v0 = beta * C
    # FELsim col4 = DToF/T_RF*1e3, with xsuite zeta = -v0*DToF and zeta ~ -z.
    # So z[m] = v0*DToF  ->  DToF = z/v0  ->  col4 = z/(v0*T_RF)*1e3.
    col4 = z_m / (v0 * T_RF) * 1e3
    col5 = delta * 1e3  # DK/K0 * 1e3

    out = np.column_stack([
        x_m * 1e3,      # mm
        xp_r * 1e3,     # mrad
        y_m * 1e3,      # mm
        yp_r * 1e3,     # mrad
        col4,
        col5,
    ]).astype(np.float64)
    return out


def to_physical_si(felsim_arr: np.ndarray, energy_mev: float) -> dict:
    """Convert a FELsim array to SI physics coordinates for an external kernel
    (e.g. libspch_kick): x,y in m; x',y' in rad; z (lab path length) in m;
    delta = DK/K0. The transverse small-angle map matches the FODO tracker."""
    gamma, beta = gamma_beta(energy_mev)
    v0 = beta * C
    z_m = felsim_arr[:, 4] * 1e-3 * (v0 * T_RF)
    return {
        "x": felsim_arr[:, 0] * 1e-3,
        "xp": felsim_arr[:, 1] * 1e-3,
        "y": felsim_arr[:, 2] * 1e-3,
        "yp": felsim_arr[:, 3] * 1e-3,
        "z": z_m,
        "delta": felsim_arr[:, 5] * 1e-3,
    }


def build(man: BeamManifest) -> tuple[np.ndarray, BeamManifest]:
    """Generate the array and stamp the manifest with versions + array hash."""
    arr = make_felsim_distribution(man)
    man.versions = collect_versions()
    man.array_sha256 = array_sha256(arr)
    return arr, man


def write_manifest(man: BeamManifest, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(man), indent=2, sort_keys=True))


def load_manifest(path: Path) -> BeamManifest:
    d = json.loads(Path(path).read_text())
    return BeamManifest(**d)


def regenerate_from_manifest(path: Path) -> np.ndarray:
    """Reload a manifest and reproduce its array; assert the hash matches."""
    man = load_manifest(path)
    arr = make_felsim_distribution(man)
    h = array_sha256(arr)
    if man.array_sha256 is not None and h != man.array_sha256:
        raise ValueError(
            f"Reproduced array hash {h[:12]} != recorded {man.array_sha256[:12]}; "
            "manifest is not reproducible in this environment."
        )
    return arr


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Build a common SC distribution + manifest")
    ap.add_argument("--out", type=Path, default=Path("results/sc_capstone/dist_default"))
    ap.add_argument("--energy-mev", type=float, default=45.0)
    ap.add_argument("--q-nc", type=float, default=1.0)
    ap.add_argument("--n-p", type=int, default=6000)
    ap.add_argument("--seed", type=int, default=20260619)
    ap.add_argument("--eps-n", type=float, default=8.0, help="normalized eps [mm.mrad]")
    ap.add_argument("--betx", type=float, default=3.0)
    ap.add_argument("--bety", type=float, default=3.0)
    ap.add_argument("--label", type=str, default="capstone_default")
    args = ap.parse_args()

    man = BeamManifest(
        energy_mev=args.energy_mev, q_nc=args.q_nc, n_p=args.n_p, seed=args.seed,
        eps_n_mm_mrad=args.eps_n, betx_m=args.betx, bety_m=args.bety, label=args.label,
    )
    arr, man = build(man)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.save(out.with_suffix(".npy"), arr)
    write_manifest(man, out.with_suffix(".json"))

    g, b = gamma_beta(man.energy_mev)
    print(f"== common distribution: {man.label} ==")
    print(f"  N_p={man.n_p}  E={man.energy_mev} MeV (betagamma={b*g:.3f})  Q={man.q_nc} nC")
    print(f"  eps_n={man.eps_n_mm_mrad} mm.mrad  betx={man.betx_m} bety={man.bety_m} m")
    print(f"  sig_x={arr[:,0].std():.4f} mm  sig_xp={arr[:,1].std():.4f} mrad")
    print(f"  sig_y={arr[:,2].std():.4f} mm  sig_yp={arr[:,3].std():.4f} mrad")
    print(f"  array sha256={man.array_sha256[:16]}")
    print(f"  saved {out.with_suffix('.npy')} + {out.with_suffix('.json')}")
