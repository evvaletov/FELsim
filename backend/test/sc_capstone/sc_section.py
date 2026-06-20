"""
Section-aware DA-FMM space-charge tracker for a real FELsim drift+quad section.

Generalizes the cosy-fmm `spch_demo.py` FODO tracker (which is hardcoded to a
thin-lens FODO) to an arbitrary sequence of drifts and thick quadrupoles taken
from the FELsim transport line, so the cosy-fmm DA-FMM Coulomb treecode can be
exercised on the genuine no-dipole transport section [32, 46) rather than a toy
cell. The linear transport per slice uses FELsim's own quad strength
    k0 = |Q . G . I| / (M . C . beta . gamma)
(read straight off the FELsim element), sliced with the exact thick-quad
focusing/defocusing matrices, and the space-charge kick is the same
libspch_kick.so call (bare 1/r or Plummer-softened) the FODO demo uses.

Author: Eremey Valetov
"""
from __future__ import annotations

import ctypes
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import common_distribution as cd

QE = 1.602176634e-19
KE_COUL = 8.987551787368176e9
MEC2_J = cd.E0_E * 1e6 * QE

DAFMM_LIB = Path("/home/evaletov/COSY/cosy-fmm/demo/spch_demo/libspch_kick.so")
THETA_DEFAULT = 0.3


@dataclass
class Segment:
    kind: str       # 'drift' | 'quad'
    length: float   # m
    kx: float       # focusing strength in x [1/m^2], signed (>0 focusing)
    ky: float       # focusing strength in y [1/m^2], signed


def extract_segments(lattice_path: str, start: int, end: int,
                     energy_mev: float,
                     optics_energy_mev: float | None = None) -> list[Segment]:
    """Build the drift/quad segment list for elements [start, end) of the
    FELsim lattice. Raises if a non-drift/quad element (dipole, wedge, ...) is
    present — this tracker is for no-dipole sections.

    The quad strength k0 = |Q.G.I|/(M.C.beta.gamma) is energy-dependent, so a
    section tuned at 45 MeV is over-focused (and optically unstable) at, e.g.,
    1 MeV. Pass `optics_energy_mev` to compute k0 at that reference energy
    instead of the tracking energy — i.e. scale the currents to PRESERVE the
    optics across energies, so an energy comparison isolates the space-charge
    scaling (1/beta^2 gamma^3) from a trivial focusing-instability blow-up
    (this is the 're-match' Risk R1 mitigation, done analytically)."""
    import latticeLoader
    bl = latticeLoader.create_beamline(lattice_path)
    k_energy = optics_energy_mev if optics_energy_mev is not None else energy_mev
    for e in bl:
        e.setE(k_energy)
    segs = []
    for e in bl[start:end]:
        cls = type(e).__name__
        L = float(getattr(e, "length", 0.0) or 0.0)
        if cls == "driftLattice":
            segs.append(Segment("drift", L, 0.0, 0.0))
        elif cls in ("qpfLattice", "qpdLattice"):
            k0 = abs(e.Q * e.G * e.current / (e.M * e.C * e.beta * e.gamma))
            if cls == "qpfLattice":      # focus x, defocus y
                segs.append(Segment("quad", L, +k0, -k0))
            else:                         # qpd: defocus x, focus y
                segs.append(Segment("quad", L, -k0, +k0))
        else:
            raise ValueError(f"extract_segments: unsupported element {cls} in "
                             f"[{start},{end}); this tracker is no-dipole only")
    return segs


def _plane_map(x, xp, k, h):
    """Advance (x, x') a length h under constant focusing strength k.
    k>0 focusing, k<0 defocusing, k==0 drift. Vectorized over particles."""
    if k == 0.0:
        return x + xp * h, xp
    if k > 0.0:
        sk = np.sqrt(k)
        c, s = np.cos(sk * h), np.sin(sk * h)
        return c * x + (s / sk) * xp, -sk * s * x + c * xp
    sk = np.sqrt(-k)
    c, s = np.cosh(sk * h), np.sinh(sk * h)
    return c * x + (s / sk) * xp, sk * s * x + c * xp


def _slice_schedule(segs: list[Segment], ds_target: float):
    """Yield (h, kx, ky) micro-steps, subdividing each segment so step <= ds_target."""
    steps = []
    for sg in segs:
        if sg.length <= 0:
            continue
        n = max(1, int(np.ceil(sg.length / ds_target)))
        h = sg.length / n
        for _ in range(n):
            steps.append((h, sg.kx, sg.ky))
    return steps


def load_lib() -> ctypes.CDLL:
    if not DAFMM_LIB.exists():
        raise FileNotFoundError(f"{DAFMM_LIB} missing; build with `make libspch_kick.so`")
    lib = ctypes.CDLL(str(DAFMM_LIB))
    arr = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    lib.spch_kick_compute_eps.restype = None
    lib.spch_kick_compute_eps.argtypes = [
        ctypes.c_int, arr, arr, arr, arr,
        ctypes.c_double, ctypes.c_double, arr, arr, arr,
    ]
    return lib


def track_dafmm(felsim_arr: np.ndarray, segs: list[Segment], energy_mev: float,
                q_nc: float, ds_target: float = 0.02, theta: float = THETA_DEFAULT,
                spch_on: bool = True, softening_eps: float = 0.0) -> dict:
    """Track the common distribution through `segs` with leapfrog
    (half-map / SC kick / half-map). Returns initial/exit eps_n + growth."""
    si = cd.to_physical_si(felsim_arr, energy_mev)
    x, xp = si["x"].copy(), si["xp"].copy()
    y, yp = si["y"].copy(), si["yp"].copy()
    z = si["z"].copy()
    n_p = len(x)

    g, b = cd.gamma_beta(energy_mev)
    bg = b * g
    sc_scale = QE * KE_COUL / MEC2_J / (g * g - 1.0)
    qmacro = q_nc * 1e-9 / n_p

    lib = load_lib() if spch_on else None
    qarr = np.full(n_p, qmacro, dtype=np.float64)
    fx = np.empty(n_p); fy = np.empty(n_p); fz = np.empty(n_p)

    def epsn():
        dx, dxp = x - x.mean(), xp - xp.mean()
        dy, dyp = y - y.mean(), yp - yp.mean()
        ex = np.sqrt(max(np.mean(dx**2) * np.mean(dxp**2) - np.mean(dx*dxp)**2, 0.0))
        ey = np.sqrt(max(np.mean(dy**2) * np.mean(dyp**2) - np.mean(dy*dyp)**2, 0.0))
        return bg * ex * 1e6, bg * ey * 1e6  # m.rad -> mm.mrad

    enx0, eny0 = epsn()
    steps = _slice_schedule(segs, ds_target)
    t_sc = 0.0
    for h, kx, ky in steps:
        x, xp = _plane_map(x, xp, kx, h / 2)
        y, yp = _plane_map(y, yp, ky, h / 2)
        if spch_on:
            z_rest = g * z
            t0 = time.perf_counter()
            lib.spch_kick_compute_eps(
                n_p, np.ascontiguousarray(x), np.ascontiguousarray(y),
                np.ascontiguousarray(z_rest), qarr, theta, softening_eps,
                fx, fy, fz)
            t_sc += time.perf_counter() - t0
            xp = xp + sc_scale * fx * h
            yp = yp + sc_scale * fy * h
        x, xp = _plane_map(x, xp, kx, h / 2)
        y, yp = _plane_map(y, yp, ky, h / 2)

    enx, eny = epsn()
    return {
        "epsnx_init": enx0, "epsny_init": eny0,
        "epsnx_exit": enx, "epsny_exit": eny,
        "epsnx_growth": (enx - enx0) / enx0,
        "epsny_growth": (eny - eny0) / eny0,
        "sigx_exit_mm": float(x.std() * 1e3), "sigy_exit_mm": float(y.std() * 1e3),
        "n_steps": len(steps), "wall_sc": t_sc,
    }
