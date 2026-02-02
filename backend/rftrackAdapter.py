"""
Adapter for RF-Track particle accelerator simulation code.

RF-Track is a tracking code developed at CERN for the design and optimisation
of particle accelerators. It solves fully relativistic equations of motion
and supports space charge, RF cavities, wakefields, and electromagnetic field maps.

See: https://pypi.org/project/RF-Track/

Author: Eremey Valetov
"""

import numpy as np
from typing import Dict, List, Optional, Any, Union
from simulatorBase import (
    SimulatorBase, SimulationResult, BeamlineElement,
    CoordinateSystem, SimulationMode
)
from beamEvolution import BeamEvolution, ElementInfo
from loggingConfig import get_logger_with_fallback
from physicalConstants import PhysicalConstants

# Attempt to import RF-Track
try:
    import RF_Track as rft
    _RFTRACK_AVAILABLE = True
except ImportError:
    _RFTRACK_AVAILABLE = False
    rft = None


class RFTrackAdapter(SimulatorBase):
    """
    Adapter providing unified interface to RF-Track simulator.

    RF-Track is a CERN-developed tracking code supporting:
    - Fully relativistic particle dynamics
    - Space charge effects (bunched and CW beams)
    - RF cavities and electromagnetic field maps (1D, 2D, 3D)
    - Synchrotron radiation and wakefields
    - Mixed particle species

    Coordinate System (RF-Track native):
        [x(m), x'(rad), y(m), y'(rad), t(s), δ]

    This adapter transforms to/from FELsim coordinates:
        [x(mm), x'(mrad), y(mm), y'(mrad), ΔToF/T(10^-3), δW/W(10^-3)]

    Examples
    --------
    Basic usage:
        >>> sim = RFTrackAdapter(beam_energy=45.0)
        >>> sim.set_beamline(elements)
        >>> particles = sim.generate_particles(1000)
        >>> result = sim.simulate(particles)

    With space charge:
        >>> sim = RFTrackAdapter(space_charge=True)
        >>> result = sim.simulate(particles)
    """

    # Class-level capabilities for factory introspection
    CAPABILITIES = {
        'particle_tracking': True,
        'transfer_matrix': False,
        'space_charge': True,
        'rf_cavities': True,
        'field_maps': True,
        'wakefields': True,
        'synchrotron_radiation': True,
    }

    NATIVE_COORDINATES = CoordinateSystem.RFTRACK

    # Default aperture for elements (5 cm) - prevents tracking hangs
    DEFAULT_APERTURE = 0.05

    def __init__(self,
                 excel_path: Optional[str] = None,
                 mode: str = 'particle_tracking',
                 space_charge: bool = False,
                 sc_mesh: Optional[tuple] = None,
                 beam_energy: float = 45.0,
                 particle_mass: Optional[float] = None,
                 particle_charge: float = -1.0,
                 aperture: float = 0.05,
                 G_quad: Optional[float] = None,
                 debug: bool = None):
        """
        Initialise RF-Track adapter.

        Parameters
        ----------
        excel_path : str, optional
            Path to Excel beamline definition file
        mode : str
            Simulation mode ('particle_tracking' only for RF-Track)
        space_charge : bool
            Enable space charge calculation
        sc_mesh : tuple, optional
            Space charge mesh size (nx, ny, nz). Default: (32, 32, 64)
        beam_energy : float
            Beam kinetic energy in MeV
        particle_mass : float, optional
            Particle rest mass in MeV/c². Default: electron mass
        particle_charge : float
            Particle charge in elementary charges. Default: -1 (electron)
        aperture : float
            Default aperture for elements in metres. Default: 0.05 (5 cm)
        G_quad : float, optional
            Quadrupole gradient calibration constant in T/A/m.
            Converts current to field gradient: gradient = G_quad × current.
            Default: 2.694 T/A/m (UH FEL quadrupoles)
        debug : bool, optional
            Enable debug logging
        """
        if not _RFTRACK_AVAILABLE:
            raise ImportError(
                "RF-Track is not installed. Install with: pip install RF-Track"
            )

        super().__init__(
            name="RF-Track",
            native_coordinates=CoordinateSystem.RFTRACK,
            debug=debug
        )

        self.logger, self.debug = get_logger_with_fallback(__name__, debug)

        # RF-Track only supports particle tracking
        if mode != 'particle_tracking':
            self.logger.warning(
                f"RF-Track only supports particle_tracking mode, ignoring '{mode}'"
            )
        self.simulation_mode = SimulationMode.PARTICLE_TRACKING

        # Beam parameters - use RF-Track's electron mass if not specified
        self.particle_mass = particle_mass if particle_mass is not None else rft.electronmass
        self.particle_charge = particle_charge
        self.beam_energy = beam_energy
        self._update_relativistic_params()

        # Default aperture for elements
        self.default_aperture = aperture

        # Quadrupole gradient calibration constant (T/A/m)
        # Default is for UH FEL quadrupoles (2.694 T/A/m)
        self.G_quad = G_quad if G_quad is not None else PhysicalConstants.G_quad_default

        # Space charge configuration
        self.space_charge_enabled = space_charge
        self.sc_mesh = sc_mesh or (32, 32, 64)
        self._space_charge_effect = None

        # RF-Track native objects
        self._lattice: Optional[rft.Lattice] = None
        self._bunch: Optional[rft.Bunch6d] = None
        self._native_elements: List[Any] = []

        # Element type mapping
        self._element_type_map = {
            'DRIFT': 'Drift',
            'QUAD_F': 'Quadrupole',
            'QPF': 'Quadrupole',
            'QUAD_D': 'Quadrupole',
            'QPD': 'Quadrupole',
            'DIPOLE': 'SBend',
            'DPH': 'SBend',
            'DIPOLE_WEDGE': 'SBend',
            'DPW': 'SBend',
            'SOLENOID': 'Solenoid',
            'RF_CAVITY': 'Cavity',
            'SEXTUPOLE': 'Sextupole',
        }

        # Load beamline if provided
        self.excel_path = excel_path
        if excel_path:
            self._load_from_excel(excel_path)

    def _update_relativistic_params(self):
        """Update relativistic parameters from beam energy."""
        self._gamma = 1 + self.beam_energy / self.particle_mass
        self._beta = np.sqrt(1 - 1/self._gamma**2)
        self._Pc = self._gamma * self._beta * self.particle_mass  # MeV/c

    def simulate(self,
                 particles: Optional[np.ndarray] = None,
                 mode: Optional[SimulationMode] = None) -> SimulationResult:
        """
        Run RF-Track particle tracking simulation.

        Parameters
        ----------
        particles : ndarray (N, 6)
            Initial distribution in FELsim coordinates:
            [x(mm), x'(mrad), y(mm), y'(mrad), ΔToF/T(10^-3), δW/W(10^-3)]
        mode : SimulationMode, optional
            Ignored (RF-Track only supports particle tracking)

        Returns
        -------
        SimulationResult
            Contains final particles, Twiss parameters, and metadata
        """
        if mode and mode != SimulationMode.PARTICLE_TRACKING:
            raise NotImplementedError(
                f"RF-Track only supports particle tracking, not {mode}"
            )

        if particles is None:
            raise ValueError("particles array required for simulation")

        self.validate_particles(particles)

        if self._lattice is None:
            raise ValueError(
                "Beamline not set. Call set_beamline() or provide excel_path"
            )

        particles_rftrack = self.transform_coordinates(
            particles, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK
        )

        self._bunch = rft.Bunch6d(
            self.particle_mass,
            self.particle_charge,
            self._Pc,
            particles_rftrack
        )

        self.logger.debug(
            f"Created bunch: {self._bunch.size()} particles, "
            f"Pc={self._Pc:.2f} MeV/c"
        )

        if self.space_charge_enabled and self._space_charge_effect is None:
            self._setup_space_charge()

        self.logger.info(
            f"Tracking {particles.shape[0]} particles through "
            f"{self._lattice.size()} elements (L={self._lattice.get_length():.3f} m)"
        )

        # Track returns the tracked bunch (original is not modified in place)
        tracked_bunch = self._lattice.track(self._bunch)

        final_rftrack = tracked_bunch.get_phase_space()
        final_particles = self.transform_coordinates(
            final_rftrack, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM
        )
        twiss = self._calculate_twiss(final_particles)

        n_good = tracked_bunch.get_ngood()
        n_lost = tracked_bunch.get_nlost()

        self.logger.info(
            f"Tracking complete: {n_good} good, {n_lost} lost particles"
        )

        return SimulationResult(
            simulator_name=self.name,
            success=True,
            final_particles=final_particles,
            twiss_parameters_statistical={'final': twiss},
            metadata={
                'num_particles': particles.shape[0],
                'num_good': n_good,
                'num_lost': n_lost,
                'beam_energy_mev': self.beam_energy,
                'momentum_mev_c': self._Pc,
                'space_charge': self.space_charge_enabled,
                'particle_mass_mev': self.particle_mass,
                'particle_charge': self.particle_charge,
                'lattice_length': self._lattice.get_length(),
            }
        )

    def optimize(self,
                 objectives: Dict,
                 variables: Dict,
                 initial_point: Dict,
                 method: Optional[str] = None,
                 **kwargs) -> SimulationResult:
        """
        Run optimization using RF-Track simulations.

        Uses scipy.optimize with RF-Track simulation as objective function.

        Parameters
        ----------
        objectives : dict
            Optimization objectives
        variables : dict
            Variable definitions
        initial_point : dict
            Initial values and bounds
        method : str, optional
            Optimization method (default: 'Nelder-Mead')
        **kwargs
            particles : ndarray (required)

        Returns
        -------
        SimulationResult
        """
        from scipy import optimize

        particles = kwargs.get('particles')
        if particles is None:
            raise ValueError("particles required for optimization")

        method = method or 'Nelder-Mead'

        def objective(x):
            for idx, (var_name, param_name, transform) in variables.items():
                value = transform(x[list(variables.keys()).index(idx)])
                self._modify_element(idx, **{param_name: value})
            self._build_lattice()
            result = self.simulate(particles)

            cost = 0.0
            twiss = result.twiss_parameters_statistical.get('final', {})
            for elem_idx, obj_list in objectives.items():
                if elem_idx == 'optimizer_settings':
                    continue
                for obj in obj_list:
                    axis, param = obj['measure']
                    goal = obj['goal']
                    weight = obj.get('weight', 1.0)
                    value = twiss.get(axis, {}).get(param, 0)
                    cost += weight * (value - goal)**2

            return cost

        x0 = [initial_point[v[0]]['start'] for v in variables.values()]
        bounds = [initial_point[v[0]].get('bounds') for v in variables.values()]

        result = optimize.minimize(
            objective, x0, method=method,
            bounds=bounds if any(bounds) else None
        )

        opt_vars = {
            var_name: transform(result.x[i])
            for i, (idx, (var_name, param_name, transform)) in enumerate(variables.items())
        }
        final = self.simulate(particles)

        return SimulationResult(
            simulator_name=self.name,
            success=result.success,
            twiss_parameters_statistical=final.twiss_parameters_statistical,
            final_particles=final.final_particles,
            optimization_variables=opt_vars,
            metadata={
                'method': method,
                'objective_value': result.fun,
                'num_iterations': getattr(result, 'nit', None),
                **final.metadata
            }
        )

    def _convert_element_to_native(self, element: BeamlineElement) -> Any:
        """
        Convert generic BeamlineElement to RF-Track element.

        Parameters
        ----------
        element : BeamlineElement
            Generic element representation

        Returns
        -------
        RF-Track element object (rft.Drift, rft.Quadrupole, etc.)
        """
        elem_type = element.element_type.upper()
        params = element.parameters
        length = element.length
        aperture = self.default_aperture

        if elem_type == 'DRIFT':
            elem = rft.Drift(length)

        elif elem_type in ['QUAD_F', 'QPF']:
            elem = rft.Quadrupole()
            elem.set_length(length)
            k1 = self._current_to_k1(params.get('current', 0.0), length, focusing=True)
            elem.set_strength(k1 * length)  # integrated strength

        elif elem_type in ['QUAD_D', 'QPD']:
            elem = rft.Quadrupole()
            elem.set_length(length)
            k1 = self._current_to_k1(params.get('current', 0.0), length, focusing=False)
            elem.set_strength(k1 * length)

        elif elem_type in ['DIPOLE', 'DPH', 'DIPOLE_WEDGE', 'DPW']:
            angle = params.get('angle', 0.0)
            if angle != 0 and length > 0:
                elem = rft.SBend()
                elem.set_length(length)
                # Use set_K0 (curvature = 1/rho = angle/length) instead of set_angle
                # for proper curvilinear coordinate tracking in RF-Track Lattice
                K0 = np.radians(angle) / length
                elem.set_K0(K0)
            else:
                elem = rft.Drift(length)

        elif elem_type == 'SOLENOID':
            elem = rft.Solenoid()
            elem.set_length(length)
            elem.set_Bz(params.get('field', 0.0))

        elif elem_type == 'SEXTUPOLE':
            elem = rft.Sextupole()
            elem.set_length(length)
            elem.set_strength(params.get('strength', 0.0))

        else:
            self.logger.warning(f"Unknown element type '{elem_type}', using drift")
            elem = rft.Drift(length)

        if hasattr(elem, 'set_aperture'):
            elem.set_aperture(aperture, aperture)

        if hasattr(elem, 'set_name') and 'name' in params:
            elem.set_name(params['name'])

        return elem

    def transform_coordinates(self,
                              particles: np.ndarray,
                              from_system: CoordinateSystem,
                              to_system: CoordinateSystem) -> np.ndarray:
        """
        Transform particle coordinates between systems.

        Coordinate systems:
        - FELSIM: [x(mm), x'(mrad), y(mm), y'(mrad), ΔToF/T(10^-3), δW/W(10^-3)]
        - RFTRACK: [x(mm), x'(mrad), y(mm), y'(mrad), t(mm/c), P(MeV/c)]

        Note: RF-Track Bunch6d uses mm, mrad, mm/c, and MeV/c units.
        Column 5 is momentum P, not energy E.

        Parameters
        ----------
        particles : ndarray (N, 6)
            Particle distribution
        from_system : CoordinateSystem
            Source coordinate system
        to_system : CoordinateSystem
            Target coordinate system

        Returns
        -------
        ndarray (N, 6)
            Transformed particle distribution
        """
        if from_system == to_system:
            return particles.copy()

        result = np.zeros_like(particles)

        if from_system == CoordinateSystem.FELSIM and to_system == CoordinateSystem.RFTRACK:
            # FELsim: [x(mm), x'(mrad), y(mm), y'(mrad), ΔToF/T×10³, δW/W×10³]
            # RF-Track: [x(mm), x'(mrad), y(mm), y'(mrad), t(mm/c), P(MeV/c)]
            # Transverse coordinates: same units, no conversion needed
            result[:, 0:4] = particles[:, 0:4]
            # Longitudinal time: ΔToF/T×10³ → t (mm/c)
            # For small deviations, approximate t ≈ 0 + relative offset
            # The FELsim coordinate is relative; RF-Track t is arrival time
            result[:, 4] = particles[:, 4]  # Keep relative timing
            # Energy → Momentum: δW/W×10³ → P = P_ref × (1 + δP/P)
            # For relativistic particles: δP/P ≈ δE/E
            result[:, 5] = self._Pc * (1.0 + particles[:, 5] * 1e-3)

        elif from_system == CoordinateSystem.RFTRACK and to_system == CoordinateSystem.FELSIM:
            # RF-Track: [x(mm), x'(mrad), y(mm), y'(mrad), t(mm/c), P(MeV/c)]
            # FELsim: [x(mm), x'(mrad), y(mm), y'(mrad), ΔToF/T×10³, δW/W×10³]
            # Transverse coordinates: same units
            result[:, 0:4] = particles[:, 0:4]
            # Time: keep as relative timing
            result[:, 4] = particles[:, 4]
            # Momentum → relative energy deviation
            # P → δW/W×10³ = (P/P_ref - 1) × 10³
            result[:, 5] = (particles[:, 5] / self._Pc - 1.0) * 1e3

        elif from_system == CoordinateSystem.COSY and to_system == CoordinateSystem.RFTRACK:
            # COSY uses m and rad; RF-Track Bunch6d uses mm and mrad
            result[:, 0:4] = particles[:, 0:4] * 1e3  # m/rad → mm/mrad
            result[:, 4] = particles[:, 4] / (self._beta * PhysicalConstants.C) * 1e3  # l(m) → t(mm/c)
            result[:, 5] = self._Pc * (1.0 + particles[:, 5])  # δ → P

        elif from_system == CoordinateSystem.RFTRACK and to_system == CoordinateSystem.COSY:
            # RF-Track uses mm/mrad; COSY uses m/rad
            result[:, 0:4] = particles[:, 0:4] * 1e-3  # mm/mrad → m/rad
            result[:, 4] = particles[:, 4] * (self._beta * PhysicalConstants.C) * 1e-3  # t(mm/c) → l(m)
            result[:, 5] = particles[:, 5] / self._Pc - 1.0  # P → δ

        else:
            raise NotImplementedError(
                f"Transformation {from_system.value} → {to_system.value} "
                "not implemented. Transform via FELSIM as intermediate."
            )

        return result

    def set_beamline(self, elements: List[Union[BeamlineElement, Any]]):
        """
        Set beamline from generic or native elements.

        Parameters
        ----------
        elements : list
            List of BeamlineElement or native RF-Track elements
        """
        super().set_beamline(elements)
        self._build_lattice()

    def _build_lattice(self):
        """Build RF-Track lattice from beamline elements."""
        self._lattice = rft.Lattice()
        self._native_elements = []

        for elem in self.beamline:
            native_elem = self._convert_element_to_native(elem)
            self._native_elements.append(native_elem)
            self._lattice.append(native_elem)

        # Set lattice aperture to match element apertures
        self._lattice.set_aperture(self.default_aperture, self.default_aperture)

        self.logger.debug(
            f"Built RF-Track lattice: {self._lattice.size()} elements, "
            f"L={self._lattice.get_length():.3f} m"
        )

    def _modify_element(self, index: int, **kwargs):
        """Modify element parameters in beamline."""
        if 0 <= index < len(self.beamline):
            for key, value in kwargs.items():
                self.beamline[index].parameters[key] = value

    def set_space_charge(self, enabled: bool, mesh: Optional[tuple] = None,
                         method: str = 'PIC'):
        """
        Configure space charge calculation.

        Parameters
        ----------
        enabled : bool
            Enable/disable space charge
        mesh : tuple, optional
            Mesh size (nx, ny, nz) for 3D solver
        method : str
            Space charge method: 'PIC', 'P2P', 'FreeSpace'
        """
        self.space_charge_enabled = enabled
        if mesh:
            self.sc_mesh = mesh
        self._sc_method = method

        if enabled:
            self._setup_space_charge()
        else:
            self._space_charge_effect = None

        self.logger.info(
            f"Space charge: {enabled}, method={method}, mesh={self.sc_mesh}"
        )

    def _setup_space_charge(self):
        """Attach space charge effect to lattice elements."""
        nx, ny, nz = self.sc_mesh
        method = getattr(self, '_sc_method', 'PIC')

        if method == 'P2P':
            self._space_charge_effect = rft.SpaceCharge_P2P()
        else:
            self._space_charge_effect = rft.SpaceCharge_PIC_FreeSpace(nx, ny, nz)

        if self._lattice is not None:
            for i in range(self._lattice.size()):
                elem = self._lattice[i]
                if hasattr(elem, 'add_collective_effect'):
                    elem.add_collective_effect(self._space_charge_effect)

    def collect_evolution(self,
                          particles: np.ndarray,
                          checkpoint_elements: Union[str, List[int]] = 'all') -> BeamEvolution:
        """
        Collect beam evolution data at element boundaries.

        Tracks element-by-element to capture phase space at each checkpoint.
        Uses RF-Track Lattice environment which tracks in curvilinear coordinates
        along the design orbit.

        Parameters
        ----------
        particles : ndarray (N, 6)
            Initial distribution in FELsim coordinates
        checkpoint_elements : str or list
            'all' or list of element indices for checkpoints

        Returns
        -------
        BeamEvolution
        """
        self.validate_particles(particles)

        if not self.beamline:
            raise ValueError("Beamline not set")

        evolution = BeamEvolution(
            simulator_name=self.name,
            num_particles=particles.shape[0],
            beam_energy=self.beam_energy
        )

        n_elements = len(self.beamline)
        if checkpoint_elements == 'all':
            checkpoint_set = set(range(n_elements))
        else:
            checkpoint_set = set(checkpoint_elements)

        # Initial state
        evolution.add_sample(0.0, particles.copy(), self._calculate_twiss(particles))

        # Convert to RF-Track coordinates
        particles_rftrack = self.transform_coordinates(
            particles, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK
        )

        # Track element by element
        s = 0.0
        for idx, elem in enumerate(self.beamline):
            # Build single-element lattice
            native_elem = self._convert_element_to_native(elem)
            single_lat = rft.Lattice()
            single_lat.append(native_elem)
            single_lat.set_aperture(self.default_aperture, self.default_aperture)

            # Create bunch for this segment
            bunch = rft.Bunch6d(
                self.particle_mass, self.particle_charge, self._Pc, particles_rftrack
            )

            # Track through element (returns tracked bunch)
            tracked_bunch = single_lat.track(bunch)

            # Update position
            s += elem.length

            # Get phase space from tracked bunch
            particles_rftrack = tracked_bunch.get_phase_space()

            # Record element info
            evolution.elements.append(ElementInfo(
                element_type=elem.element_type,
                s_start=s - elem.length,
                s_end=s,
                length=elem.length,
                color=self._get_element_color(elem.element_type),
                index=idx,
                parameters=elem.parameters
            ))

            # Checkpoint if requested
            if idx in checkpoint_set and particles_rftrack.size > 0:
                ps_felsim = self.transform_coordinates(
                    particles_rftrack, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM
                )
                evolution.add_sample(s, ps_felsim, self._calculate_twiss(ps_felsim))

        evolution.total_length = s
        return evolution

    def _load_from_excel(self, excel_path: str):
        """Load beamline from Excel specification."""
        from excelElements import ExcelElements

        try:
            excel = ExcelElements(excel_path)
            native_elements = excel.create_beamline()

            self.beamline = []
            for elem in native_elements:
                self.beamline.append(self._convert_element_from_native(elem))

            self._build_lattice()
            self.logger.info(f"Loaded {len(self.beamline)} elements from {excel_path}")

        except Exception as e:
            self.logger.error(f"Failed to load beamline from {excel_path}: {e}")
            raise

    def _convert_element_from_native(self, native_elem: Any) -> BeamlineElement:
        """Convert FELsim native element to generic BeamlineElement."""
        cls_name = type(native_elem).__name__

        type_map = {
            'driftLattice': 'DRIFT',
            'qpfLattice': 'QUAD_F',
            'qpdLattice': 'QUAD_D',
            'dipole': 'DIPOLE',
            'dipole_wedge': 'DIPOLE_WEDGE',
        }

        elem_type = type_map.get(cls_name, cls_name.upper())

        params = {}
        if hasattr(native_elem, 'current'):
            params['current'] = native_elem.current
        if hasattr(native_elem, 'angle'):
            params['angle'] = native_elem.angle
        if hasattr(native_elem, 'fringeType'):
            params['fringe_type'] = native_elem.fringeType

        return BeamlineElement(
            element_type=elem_type,
            length=native_elem.length,
            **params
        )

    def _calculate_twiss(self, particles: np.ndarray) -> Dict:
        """
        Calculate Twiss parameters from particle distribution (FELsim coords).

        FELsim coordinates: [x(mm), x'(mrad), y(mm), y'(mrad), ...]
        Returns beta in m, gamma in rad/m, emittance in π·mm·mrad.
        """
        if particles.shape[0] < 2:
            return {}

        twiss = {}
        for plane, (pos_idx, ang_idx) in [('x', (0, 1)), ('y', (2, 3))]:
            cov = np.cov(particles[:, pos_idx], particles[:, ang_idx], ddof=1)
            sig_x2, sig_xp2, sig_xxp = cov[0, 0], cov[1, 1], cov[0, 1]

            emit_sq = sig_x2 * sig_xp2 - sig_xxp**2
            emittance = np.sqrt(max(0, emit_sq))  # π·mm·mrad

            if emittance > 0:
                # beta = mm²/(mm·mrad) = mm/mrad = m
                beta = sig_x2 / emittance
                alpha = -sig_xxp / emittance
                # gamma = mrad²/(mm·mrad) = mrad/mm = rad/m
                gamma = sig_xp2 / emittance
            else:
                beta = alpha = gamma = 0.0

            twiss[plane] = {'beta': beta, 'alpha': alpha, 'gamma': gamma, 'emittance': emittance}

        return twiss

    def _current_to_k1(self, current: float, length: float, focusing: bool = True) -> float:
        """
        Convert quadrupole current to normalized gradient k1.

        Uses k = |Q·G·I| / (M·C·β·γ), consistent with FELsim's beamline.py.
        """
        if length <= 0 or current == 0:
            return 0.0

        mass_kg = self.particle_mass * PhysicalConstants.MeV_to_J / PhysicalConstants.C**2
        k1 = abs(PhysicalConstants.Q * self.G_quad * current) / (
            mass_kg * PhysicalConstants.C * self._beta * self._gamma
        )
        return k1 if focusing else -k1

    def _get_element_color(self, elem_type: str) -> str:
        """Map element type to display color."""
        colors = {
            'DRIFT': 'white',
            'QUAD_F': 'cornflowerblue',
            'QPF': 'cornflowerblue',
            'QUAD_D': 'lightcoral',
            'QPD': 'lightcoral',
            'DIPOLE': 'forestgreen',
            'DPH': 'forestgreen',
            'DIPOLE_WEDGE': 'lightgreen',
            'DPW': 'lightgreen',
            'SOLENOID': 'purple',
            'RF_CAVITY': 'gold',
            'SEXTUPOLE': 'orange',
        }
        return colors.get(elem_type.upper(), 'gray')

    def generate_particles(self,
                           num_particles: int = 1000,
                           distribution_type: str = 'gaussian',
                           **parameters) -> np.ndarray:
        """
        Generate initial particle distribution in FELsim coordinates.

        Parameters
        ----------
        num_particles : int
            Number of particles
        distribution_type : str
            'gaussian', 'uniform', 'waterbag', 'kv', or 'twiss'
        **parameters
            std_dev : list of 6 RMS values [mm, mrad, mm, mrad, -, -]
            twiss_x, twiss_y : dict with beta, alpha, emittance for 'twiss' type

        Returns
        -------
        ndarray (N, 6) in FELsim coordinates
        """
        std_dev = parameters.get('std_dev', [1.0, 0.1, 1.0, 0.1, 1.0, 0.1])
        mean = parameters.get('mean', 0.0)

        if distribution_type == 'gaussian':
            particles = np.random.randn(num_particles, 6) * std_dev + mean

        elif distribution_type == 'uniform':
            half_width = np.array(std_dev) * np.sqrt(3)
            particles = np.random.uniform(-half_width, half_width, (num_particles, 6))

        elif distribution_type == 'twiss':
            twiss_x = parameters.get('twiss_x', {'beta': 10, 'alpha': 0, 'emittance': 1})
            twiss_y = parameters.get('twiss_y', {'beta': 10, 'alpha': 0, 'emittance': 1})

            particles = np.zeros((num_particles, 6))
            for plane, twiss, idx in [('x', twiss_x, 0), ('y', twiss_y, 2)]:
                beta, alpha, emit = twiss['beta'], twiss['alpha'], twiss['emittance']
                u1, u2 = np.random.randn(num_particles), np.random.randn(num_particles)

                sigma_x = np.sqrt(emit * beta)
                sigma_xp = np.sqrt(emit / beta) if beta > 0 else 0

                particles[:, idx] = sigma_x * u1
                particles[:, idx+1] = sigma_xp * (-alpha * u1 / np.sqrt(beta) + u2) if beta > 0 else 0

            particles[:, 4] = np.random.randn(num_particles) * std_dev[4]
            particles[:, 5] = np.random.randn(num_particles) * std_dev[5]

        else:
            self.logger.warning(f"Distribution '{distribution_type}' not implemented, using Gaussian")
            particles = np.random.randn(num_particles, 6) * std_dev + mean

        return particles

    def set_beam_energy(self, energy_mev: float):
        """Set beam kinetic energy."""
        super().set_beam_energy(energy_mev)
        self.beam_energy = energy_mev
        self._update_relativistic_params()
        self.logger.debug(f"Energy: {energy_mev} MeV, γ={self._gamma:.2f}, Pc={self._Pc:.2f} MeV/c")

    def set_particle_type(self, mass_mev: float, charge: float):
        """Set particle species."""
        self.particle_mass = mass_mev
        self.particle_charge = charge
        self._update_relativistic_params()
        self.logger.info(f"Particle: m={mass_mev} MeV/c², q={charge}e")

    def set_quadrupole_gradient(self, G_quad: float):
        """Set quadrupole gradient calibration (T/A/m). Default: 2.694 for UH FEL."""
        self.G_quad = G_quad
        self.logger.info(f"Quadrupole gradient calibration: G = {G_quad:.4f} T/A/m")
        if self._lattice is not None and self.beamline:
            self._build_lattice()

    def supports_mode(self, mode: SimulationMode) -> bool:
        return mode == SimulationMode.PARTICLE_TRACKING

    def supports_optimization(self) -> bool:
        return True

    def get_capabilities(self) -> Dict[str, Any]:
        caps = super().get_capabilities()
        caps.update({
            'space_charge': self.space_charge_enabled,
            'rf_cavities': True,
            'field_maps': True,
            'particle_mass_mev': self.particle_mass,
            'particle_charge': self.particle_charge,
            'momentum_mev_c': self._Pc,
        })
        return caps

    def get_lattice(self) -> Any:
        """Get underlying RF-Track Lattice object."""
        return self._lattice

    def get_bunch(self) -> Any:
        """Get current RF-Track Bunch6d object."""
        return self._bunch


# Convenience function
def create_rftrack_simulator(**kwargs) -> RFTrackAdapter:
    """Create RF-Track simulator instance."""
    return RFTrackAdapter(**kwargs)
