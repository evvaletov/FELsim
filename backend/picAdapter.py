"""PIC space-charge adapter for the FELsim multi-code framework.

A fourth space-charge engine alongside RF-Track and Xsuite. It tracks the
beamline with FELsim's own element maps and inserts space-charge kicks computed
by the cosy-pic mesh solver (pic-core-cli, via the cosy-pic ``tools/pic_sc.py``
relativistic driver) in a 2nd-order leapfrog (half-kick / map / kick).

cosy-pic owns the EXTERNAL transport space charge (gun -> undulator); Genesis4
owns the FEL interaction inside the undulator. cosy-pic models bunch self-field
space charge, NOT CSR -- use a CSR code at the chicane/compression waist if it
matters (cosy-pic git-bug 6b46f95).

Author: Eremey Valetov
"""

import os
import sys
from typing import Any, Dict, List, Optional

import numpy as np

from simulatorBase import (
    SimulatorBase, SimulationResult, BeamlineElement,
    CoordinateSystem, SimulationMode,
)
import latticeLoader
from ebeam import beam as ebeam_class

# cosy-pic Python SC engine (tools/pic_sc.py + particles_io.py). Optional dep:
# the adapter only registers when it imports cleanly and pic-core-cli is built.
_COSY_PIC = os.environ.get("COSY_PIC_DIR", os.path.expanduser("~/COSY/cosy-pic"))
_PIC_AVAILABLE = False
try:
    sys.path.insert(0, os.path.join(_COSY_PIC, "tools"))
    from particles_io import Grid          # noqa: E402
    from pic_sc import Beam, SpaceChargeDriver, gamma_beta, C_LIGHT  # noqa: E402
    _PIC_AVAILABLE = True
except Exception:  # pragma: no cover - missing optional dependency
    _PIC_AVAILABLE = False

_QE = 1.602176634e-19           # C
_ME_KG = 9.1093837015e-31       # kg
_MEC2_MEV = 0.51099895000       # MeV
_F_RF = 2.856e9                 # Hz (S-band)


class PicAdapter(SimulatorBase):
    """Mesh-based space-charge tracker (cosy-pic) for the multi-code framework."""

    def __init__(self, lattice_path: Optional[str] = None,
                 excel_path: Optional[str] = None,
                 beam_energy: float = 45.0,
                 space_charge: bool = True,
                 sc_mesh: tuple = (64, 64, 64),
                 sc_method: str = "pic3d",
                 bunch_charge_nc: float = 1.0,
                 sc_ds_max: float = 0.05,
                 cli_path: Optional[str] = None,
                 debug: bool = None, **kwargs):
        super().__init__(name="PIC", native_coordinates=CoordinateSystem.FELSIM)
        if not _PIC_AVAILABLE:
            raise ImportError(
                "cosy-pic SC engine unavailable; set COSY_PIC_DIR and build "
                "pic-core-cli (cmake --build build)")
        self.simulation_mode = SimulationMode.PARTICLE_TRACKING
        self.beam_energy = beam_energy
        self.space_charge_enabled = space_charge
        self.sc_mesh = tuple(sc_mesh)
        self.sc_method = sc_method
        self.bunch_charge_C = bunch_charge_nc * 1e-9
        # max SC integration step [m]: long drifts are split into sub-steps so
        # SC is applied along them, not only at element boundaries (A3).
        self.sc_ds_max = sc_ds_max
        self.debug = bool(debug)
        self._native_beamline: List[Any] = []
        self._ebeam = ebeam_class()
        self._cli_path = cli_path or os.environ.get("PIC_CORE_CLI")

        from beamline import driftLattice, qpfLattice, qpdLattice, dipole, dipole_wedge
        self._elem_map = {
            'DRIFT': driftLattice, 'QUAD_F': qpfLattice, 'QPF': qpfLattice,
            'QUAD_D': qpdLattice, 'QPD': qpdLattice, 'DIPOLE': dipole, 'DPH': dipole,
            'DIPOLE_WEDGE': dipole_wedge, 'DPW': dipole_wedge,
        }

        self._gamma0, self._beta0 = gamma_beta(beam_energy, _MEC2_MEV)
        self._t_rf = 1.0 / _F_RF
        # one driver, reused; its grid is rebuilt adaptively per SC station.
        self._drv = SpaceChargeDriver(
            self._placeholder_grid(), self._gamma0, self._beta0,
            cli_path=self._cli_path)

        path = lattice_path or excel_path
        if path:
            self._load_lattice(path)

    # ---- core interface ----------------------------------------------------

    def simulate(self, particles: Optional[np.ndarray] = None,
                 mode: Optional[SimulationMode] = None) -> SimulationResult:
        if mode and mode != SimulationMode.PARTICLE_TRACKING:
            raise NotImplementedError("PIC adapter only supports particle tracking")
        if particles is None:
            raise ValueError("particles required")
        if not self._native_beamline:
            raise ValueError("Beamline not set")

        P = np.asarray(particles, dtype=float).copy()
        sc = self.space_charge_enabled and self.bunch_charge_C > 0.0
        # With SC on, split long drifts so SC is integrated along them (A3).
        bl = self._expanded_beamline(self._native_beamline) if sc else self._native_beamline
        n = len(bl)
        n_oob = 0

        # 2nd-order leapfrog: half-kick, then (map, boundary-kick) per element.
        if sc and n > 0 and bl[0].length > 0:
            n_oob = max(n_oob, self._sc_kick(P, 0.5 * bl[0].length))
        for i in range(n):
            P = bl[i].useMatrice(P)
            if sc:
                ds = 0.5 * bl[i].length + (0.5 * bl[i + 1].length if i + 1 < n else 0.0)
                if ds > 0:
                    n_oob = max(n_oob, self._sc_kick(P, ds))

        _, _, twiss_df = self._ebeam.cal_twiss(P, ddof=1)
        twiss_dict = {axis: twiss_df.loc[axis].to_dict() for axis in twiss_df.index}
        return SimulationResult(
            simulator_name=self.name,
            success=True,
            twiss_parameters_statistical={'final': twiss_dict},
            final_particles=P,
            metadata={
                'num_particles': P.shape[0],
                'num_elements': n,
                'beam_energy_mev': self.beam_energy,
                'space_charge': sc,
                'sc_mesh': self.sc_mesh,
                'sc_ds_max': self.sc_ds_max,
                'num_sc_stations': n if sc else 0,
                'bunch_charge_nc': self.bunch_charge_C * 1e9,
                'max_particles_out_of_grid': n_oob,
            })

    # ---- space charge ------------------------------------------------------

    def set_space_charge(self, enabled: bool, mesh: Optional[tuple] = None,
                         method: Optional[str] = None,
                         bunch_charge_nc: Optional[float] = None):
        self.space_charge_enabled = enabled
        if mesh is not None:
            self.sc_mesh = tuple(mesh)
        if method is not None:
            self.sc_method = method
        if bunch_charge_nc is not None:
            self.bunch_charge_C = bunch_charge_nc * 1e-9

    def _sc_kick(self, P: np.ndarray, ds: float) -> int:
        """Apply one SC kick over path ``ds`` to the FELsim bunch P in place.
        Only the momentum-derived coords (x', y', dK/K0) change."""
        beam = self._felsim_to_beam(P)
        self._drv.grid = self._adaptive_grid(beam)
        n_oob = self._drv.sc_kick(beam, ds)
        self._beam_to_felsim(P, beam)
        return n_oob

    def _felsim_to_beam(self, P: np.ndarray) -> "Beam":
        ke0_j = self.beam_energy * 1e6 * _QE
        mec2_j = _MEC2_MEV * 1e6 * _QE
        ke_j = ke0_j * (1.0 + P[:, 5] * 1e-3)
        pmag = np.sqrt(ke_j * ke_j + 2.0 * ke_j * mec2_j) / C_LIGHT
        tx = P[:, 1] * 1e-3
        ty = P[:, 3] * 1e-3
        denom = np.sqrt(1.0 + tx * tx + ty * ty)
        npart = P.shape[0]
        return Beam(
            x=P[:, 0] * 1e-3, y=P[:, 2] * 1e-3,
            z=-self._beta0 * C_LIGHT * (P[:, 4] * 1e-3) * self._t_rf,
            px=pmag * tx / denom, py=pmag * ty / denom, pz=pmag / denom,
            q=np.full(npart, -_QE), m=np.full(npart, _ME_KG),
            w=np.full(npart, (self.bunch_charge_C / _QE) / npart),
        )

    def _beam_to_felsim(self, P: np.ndarray, beam: "Beam") -> None:
        # positions (x, y, dToF) are unchanged by a momentum kick; update angles + energy.
        mec2_j = _MEC2_MEV * 1e6 * _QE
        ke0_j = self.beam_energy * 1e6 * _QE
        P[:, 1] = (beam.px / beam.pz) * 1e3
        P[:, 3] = (beam.py / beam.pz) * 1e3
        pc = np.sqrt(beam.px**2 + beam.py**2 + beam.pz**2) * C_LIGHT
        ke_j = np.sqrt(pc * pc + mec2_j * mec2_j) - mec2_j
        P[:, 5] = (ke_j / ke0_j - 1.0) * 1e3

    def _adaptive_grid(self, beam: "Beam") -> "Grid":
        nx, ny, nz = self.sc_mesh
        sx = max(beam.x.std(), 1e-9)
        sy = max(beam.y.std(), 1e-9)
        sz = max(beam.z.std(), 1e-9) * self._gamma0   # rest-frame longitudinal
        return Grid(nx=nx, ny=ny, nz=nz,
                    dx=16 * sx / nx, dy=16 * sy / ny, dz=16 * sz / nz,
                    x0=-8 * sx, y0=-8 * sy, z0=-8 * sz,
                    bc_x=1, bc_y=1, bc_z=1)   # all Hockney open

    def _placeholder_grid(self) -> "Grid":
        return Grid(nx=2, ny=2, nz=2, dx=1.0, dy=1.0, dz=1.0,
                    x0=-1.0, y0=-1.0, z0=-1.0, bc_x=1, bc_y=1, bc_z=1)

    # ---- lattice + base interface -----------------------------------------

    def _load_lattice(self, lattice_path: str):
        native = latticeLoader.create_beamline(lattice_path)
        for elem in native:
            elem.setE(self.beam_energy)
        self._native_beamline = native

    def _expanded_beamline(self, beamline):
        """Split long drifts into <= sc_ds_max sub-drifts so the leapfrog SC
        integration has multiple stations along them (drifts compose exactly:
        driftLattice(L/n) applied n times == driftLattice(L)). Focusing/bending
        elements are short and kept whole (boundary SC kicks)."""
        from beamline import driftLattice
        if not (self.sc_ds_max and self.sc_ds_max > 0):
            return list(beamline)
        out = []
        for elem in beamline:
            L = elem.length
            if type(elem).__name__ == 'driftLattice' and L > self.sc_ds_max:
                n = int(np.ceil(L / self.sc_ds_max))
                for _ in range(n):
                    d = driftLattice(L / n)
                    d.setE(self.beam_energy)
                    out.append(d)
            else:
                out.append(elem)
        return out

    def set_beam_energy(self, energy_mev: float):
        self.beam_energy = energy_mev
        self._gamma0, self._beta0 = gamma_beta(energy_mev, _MEC2_MEV)
        self._drv.gamma = self._gamma0
        self._drv.beta = self._beta0
        for elem in self._native_beamline:
            elem.setE(energy_mev)

    def set_beamline(self, elements: List[Any]):
        if elements and isinstance(elements[0], BeamlineElement):
            self._native_beamline = [self._convert_element_to_native(e) for e in elements]
        else:
            self._native_beamline = list(elements)
        for elem in self._native_beamline:
            elem.setE(self.beam_energy)

    def get_native_beamline(self) -> List:
        return self._native_beamline

    def _convert_element_to_native(self, elem: BeamlineElement) -> Any:
        # Reuse FELsim's element mapping (same native lattice classes).
        from felsimAdapter import FELsimAdapter
        return FELsimAdapter._convert_element_to_native(self, elem)

    def transform_coordinates(self, particles: np.ndarray,
                              from_system: CoordinateSystem,
                              to_system: CoordinateSystem) -> np.ndarray:
        if from_system == to_system:
            return particles.copy()
        raise NotImplementedError(
            f"PIC adapter coordinate transform {from_system.value} -> {to_system.value}")

    def supports_mode(self, mode: SimulationMode) -> bool:
        return mode == SimulationMode.PARTICLE_TRACKING

    def get_capabilities(self) -> Dict:
        return {
            'name': self.name,
            'native_coordinates': self.native_coordinates.value,
            'space_charge': True,
            'sc_method': self.sc_method,
            'modes': [SimulationMode.PARTICLE_TRACKING.value],
        }
