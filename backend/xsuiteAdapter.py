"""
Xsuite adapter for the unified simulator interface.

Wraps xsuite (xtrack/xpart/xfields) behind SimulatorBase so it can be chained
in MultiCodeSimulator alongside FELsim, COSY, and RF-Track. Like the other
adapters, simulate() accepts and returns particles in FELsim coordinates and
converts to/from xsuite's native frame internally.

Supported elements: drift, quadrupole (thick, current→k1 via the same
calibration as FELsim/RF-Track), sector bend (basic, no edge/fringe model).
Optional space charge via xfields (frozen-Gaussian or 3D PIC), inserted
split-operator between short sub-slices of each element.

Dipole edges/fringe and RF cavities are not yet modelled (treated as drift
with a warning); the immediate use case is the focusing/transport line SC
benchmark, with and without bends.

Author: Eremey Valetov
"""

import logging
from typing import Any, Dict, List, Optional, Union

import numpy as np

from simulatorBase import (
    SimulatorBase, SimulationResult, BeamlineElement,
    CoordinateSystem, SimulationMode,
)
from physicalConstants import PhysicalConstants
from ebeam import beam as ebeam_class

logger = logging.getLogger(__name__)

try:
    import xtrack as xt
    import xpart as xp
    import xfields as xf
    _XSUITE_AVAILABLE = True
except ImportError:
    _XSUITE_AVAILABLE = False
    xt = xp = xf = None

# Horizontal-focusing sign for a QUAD_F in xsuite's k1 convention (electrons,
# q0=-1). Validated against FELsim on a FODO cell; flip if planes swap.
_FOCUS_SIGN = 1.0


class XsuiteAdapter(SimulatorBase):
    """Adapter exposing xsuite tracking through the SimulatorBase interface."""

    NATIVE_COORDINATES = CoordinateSystem.XSUITE

    def __init__(self,
                 lattice_path: Optional[str] = None,
                 excel_path: Optional[str] = None,
                 space_charge: bool = False,
                 sc_method: str = "frozen",
                 sc_mesh: tuple = (64, 64, 64),
                 n_slice_sc: int = 80,
                 bunch_charge_nc: float = 1.0,
                 beam_energy: float = 45.0,
                 particle_mass: Optional[float] = None,
                 particle_charge: float = -1.0,
                 debug: bool = None,
                 **kwargs):
        if not _XSUITE_AVAILABLE:
            raise ImportError(
                "xsuite not available. Install with `pip install xsuite`."
            )
        super().__init__(name="Xsuite", native_coordinates=CoordinateSystem.XSUITE)

        self.simulation_mode = SimulationMode.PARTICLE_TRACKING
        self.beam_energy = beam_energy
        self.particle_mass = particle_mass or PhysicalConstants.E0_electron
        self.particle_charge = particle_charge
        self.G_quad = PhysicalConstants.G_quad_default

        self.space_charge_enabled = space_charge
        self.sc_method = sc_method        # 'frozen' | 'pic3d'
        self.sc_mesh = tuple(sc_mesh)
        self.n_slice_sc = n_slice_sc
        self.bunch_charge_C = bunch_charge_nc * 1e-9

        self._ebeam = ebeam_class()
        self._update_energy()

        path = lattice_path or excel_path
        if path:
            import latticeLoader
            native = latticeLoader.create_beamline(path)
            for e in native:
                e.setE(self.beam_energy)
            self.set_beamline(native)

    # ------------------------------------------------------------------ energy
    def _update_energy(self):
        m = self.particle_mass
        self._gamma = 1.0 + self.beam_energy / m
        self._beta = np.sqrt(1.0 - 1.0 / self._gamma ** 2)
        self._betagamma = self._beta * self._gamma

    def set_beam_energy(self, energy_mev: float):
        super().set_beam_energy(energy_mev)
        self._update_energy()

    # --------------------------------------------------------------- beamline
    def set_beamline(self, elements: List[Union[BeamlineElement, Any]]):
        """Accept generic BeamlineElement objects or native FELsim elements."""
        if not elements:
            self.beamline = []
            return
        if isinstance(elements[0], BeamlineElement):
            self.beamline = list(elements)
        else:
            self.beamline = [self._convert_element_from_native(e) for e in elements]

    def _convert_element_from_native(self, elem) -> BeamlineElement:
        cls_name = type(elem).__name__
        type_map = {
            'driftLattice': 'DRIFT', 'qpfLattice': 'QUAD_F',
            'qpdLattice': 'QUAD_D', 'dipole': 'DIPOLE',
            'dipole_wedge': 'DIPOLE_WEDGE', 'rfCavityLattice': 'RF_CAVITY',
        }
        etype = type_map.get(cls_name, cls_name.upper())
        params = {}
        for attr in ('current', 'angle'):
            if hasattr(elem, attr):
                params[attr] = getattr(elem, attr)
        return BeamlineElement(etype, getattr(elem, 'length', 0.0), **params)

    def _convert_element_to_native(self, element: BeamlineElement) -> Any:
        # xsuite elements are built lazily in _build_line; nothing to do here.
        return element

    def _current_to_k1(self, current: float, focusing: bool) -> float:
        """k1 = |Q·G·I| / (m·c·β·γ), matching FELsim / RF-Track."""
        if current == 0:
            return 0.0
        mass_kg = self.particle_mass * PhysicalConstants.MeV_to_J / PhysicalConstants.C ** 2
        charge_C = abs(self.particle_charge) * PhysicalConstants.Q
        k1 = abs(charge_C * self.G_quad * current) / (
            mass_kg * PhysicalConstants.C * self._beta * self._gamma)
        return _FOCUS_SIGN * (k1 if focusing else -k1)

    # ----------------------------------------------------------- line builder
    def _xsuite_elements_for(self, elem: BeamlineElement, length: float):
        """Return the xsuite element(s) for one FELsim element of given length
        (length may be a sub-slice of the full element when space charge is on)."""
        etype = elem.element_type.upper()
        if etype == 'DRIFT' or length <= 0:
            return [xt.Drift(length=length)]
        if etype in ('QUAD_F', 'QPF', 'QUAD_D', 'QPD'):
            focusing = etype in ('QUAD_F', 'QPF')
            k1 = self._current_to_k1(elem.parameters.get('current', 0.0), focusing)
            return [xt.Quadrupole(length=length, k1=k1)]
        if etype in ('DIPOLE', 'DPH'):
            ang = float(elem.parameters.get('angle', 0.0))
            full = elem.length or length
            h = ang / full if full else 0.0
            return [xt.Bend(length=length, k0=h, h=h)]
        if etype in ('DIPOLE_WEDGE', 'DPW', 'RF_CAVITY'):
            logger.warning("Xsuite: %s not modelled; treated as drift", etype)
            return [xt.Drift(length=length)]
        logger.warning("Xsuite: unknown element %s; treated as drift", etype)
        return [xt.Drift(length=length)]

    def _make_sc_element(self, length, sig_x, sig_y, long_profile):
        if self.sc_method == "pic3d":
            nx, ny, nz = self.sc_mesh
            return xf.SpaceCharge3D(
                length=length, update_on_track=True,
                x_range=(-8 * sig_x, 8 * sig_x),
                y_range=(-8 * sig_y, 8 * sig_y),
                z_range=(-8 * long_profile.sigma_z, 8 * long_profile.sigma_z),
                nx=nx, ny=ny, nz=nz,
                solver="FFTSolver2p5D", gamma0=self._gamma)
        return xf.SpaceChargeBiGaussian(
            length=length, longitudinal_profile=long_profile,
            sigma_x=sig_x, sigma_y=sig_y, mean_x=0.0, mean_y=0.0,
            update_on_track=True)

    def _build_line(self, sc_on, sig_x=None, sig_y=None, sig_z=None, n_e=None):
        elements = []
        if not sc_on:
            for elem in self.beamline:
                elements += self._xsuite_elements_for(elem, elem.length)
        else:
            total_L = sum(e.length for e in self.beamline)
            ds = total_L / max(self.n_slice_sc, 1)
            long_profile = xf.LongitudinalProfileQGaussian(
                number_of_particles=n_e, sigma_z=sig_z, z0=0.0, q_parameter=1.0)
            for elem in self.beamline:
                L = elem.length
                if L <= 0:
                    elements += self._xsuite_elements_for(elem, 0.0)
                    continue
                n_sub = max(1, int(round(L / ds)))
                sub_L = L / n_sub
                for _ in range(n_sub):
                    elements += self._xsuite_elements_for(elem, sub_L)
                    elements.append(
                        self._make_sc_element(sub_L, sig_x, sig_y, long_profile))

        line = xt.Line(elements=elements)
        line.particle_ref = xp.Particles(
            mass0=self.particle_mass * 1e6, q0=self.particle_charge,
            kinetic_energy0=self.beam_energy * 1e6)
        line.build_tracker()
        return line

    # -------------------------------------------------------------- simulate
    def simulate(self, particles: Optional[np.ndarray] = None,
                 mode: Optional[SimulationMode] = None) -> SimulationResult:
        if mode and mode != SimulationMode.PARTICLE_TRACKING:
            raise NotImplementedError("Xsuite adapter only supports particle tracking")
        if particles is None:
            raise ValueError("particles required")
        if not self.beamline:
            raise ValueError("Beamline not set")
        self.validate_particles(particles)

        xs = self.transform_coordinates(
            particles, CoordinateSystem.FELSIM, CoordinateSystem.XSUITE)
        n_p = xs.shape[0]
        n_e = self.bunch_charge_C / PhysicalConstants.Q

        p = xp.Particles(
            mass0=self.particle_mass * 1e6, q0=self.particle_charge,
            kinetic_energy0=self.beam_energy * 1e6,
            x=xs[:, 0], px=xs[:, 1], y=xs[:, 2], py=xs[:, 3],
            zeta=xs[:, 4], delta=xs[:, 5],
            weight=(n_e / n_p if self.space_charge_enabled else 1.0))

        if self.space_charge_enabled:
            sig_x = float(np.std(xs[:, 0]))
            sig_y = float(np.std(xs[:, 2]))
            sig_z = float(np.std(xs[:, 4])) or 1e-6
            line = self._build_line(True, sig_x, sig_y, sig_z, n_e)
        else:
            line = self._build_line(False)

        line.track(p)

        alive = p.state > 0
        n_good = int(alive.sum())
        n_lost = n_p - n_good
        if n_good == 0:
            return SimulationResult(simulator_name=self.name, success=False,
                                    final_particles=np.empty((0, 6)),
                                    metadata={'num_lost': n_lost})

        xs_out = np.column_stack([
            np.asarray(p.x[alive]), np.asarray(p.px[alive]),
            np.asarray(p.y[alive]), np.asarray(p.py[alive]),
            np.asarray(p.zeta[alive]), np.asarray(p.delta[alive])])
        final = self.transform_coordinates(
            xs_out, CoordinateSystem.XSUITE, CoordinateSystem.FELSIM)

        twiss = self._calc_twiss(final)
        return SimulationResult(
            simulator_name=self.name, success=True,
            final_particles=final,
            twiss_parameters_statistical={'final': twiss},
            metadata={
                'num_particles': n_p, 'num_good': n_good, 'num_lost': n_lost,
                'beam_energy_mev': self.beam_energy,
                'space_charge': self.space_charge_enabled,
                'sc_method': self.sc_method if self.space_charge_enabled else None,
                'bunch_charge_nc': self.bunch_charge_C * 1e9,
            })

    def _calc_twiss(self, felsim_particles: np.ndarray) -> dict:
        if felsim_particles.shape[0] < 2:
            return {}
        _, _, twiss_df = self._ebeam.cal_twiss(felsim_particles, ddof=1)
        return {axis: twiss_df.loc[axis].to_dict() for axis in twiss_df.index}

    # ----------------------------------------------------------- coordinates
    def transform_coordinates(self, particles, from_system, to_system):
        from simulatorFactory import CoordinateTransformer
        return CoordinateTransformer.transform(
            particles, from_system, to_system, self.beam_energy,
            particle_mass_mev=self.particle_mass)

    # ------------------------------------------------------------- SC config
    def set_space_charge(self, enabled: bool, mesh: tuple = None,
                         method: str = None, bunch_charge_nc: float = None):
        self.space_charge_enabled = enabled
        if mesh is not None:
            self.sc_mesh = tuple(mesh)
        if method is not None:
            self.sc_method = method
        if bunch_charge_nc is not None:
            self.bunch_charge_C = bunch_charge_nc * 1e-9

    def supports_mode(self, mode: SimulationMode) -> bool:
        return mode == SimulationMode.PARTICLE_TRACKING
