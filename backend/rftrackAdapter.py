"""
Adapter for RF-Track particle accelerator simulation code.

RF-Track is a tracking code developed at CERN for the design and optimisation
of particle accelerators. It solves fully relativistic equations of motion
and supports space charge, RF cavities, wakefields, and electromagnetic field maps.

See: https://pypi.org/project/RF-Track/

Author: Eremey Valetov
"""

import math
import operator
import warnings

import numpy as np
from scipy.optimize import brentq
from typing import Dict, List, Optional, Any, Union
from simulatorBase import (
    SimulatorBase, SimulationResult, BeamlineElement,
    CoordinateSystem, SimulationMode, wedge_bend
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

    Coordinate System (RF-Track Bunch6d):
        [x(mm), x'(mrad), y(mm), y'(mrad), t(mm/c), P(MeV/c)]

    This adapter transforms to/from FELsim coordinates:
        [x(mm), x'(mrad), y(mm), y'(mrad), ΔToF/T_RF×10³, ΔK/K₀×10³]

    Bunches are built per particle as [X XP Y YP T P MASS Q N ID]
    (``_make_bunch``): Q = particle_charge, N = bunch charge / e / n for each
    of the n particles handed to a tracking call (n_initial for
    ``track_elements``), ID = the input row index (0-based). A zero-weight
    reference particle (on axis, nominal momentum, t = 0, ID -1) leads every
    bunch, or with space charge on is tracked in a bunch of its own: RF
    cavities are phased on it, and the longitudinal output of ``simulate``,
    ``collect_evolution`` and ``track_elements`` is the arrival time
    relative to it, so the reference stays at coord5 = 0 and an on-momentum
    particle keeps its coord5 through drifts and quadrupoles.
    ``transform_coordinates`` is a pure unit conversion.

    Quadrupole and sextupole strengths are normalised, as in FELsim and the
    Xsuite adapter, so QUAD_F focuses in x whatever the sign of
    particle_charge. Solenoids act through the charge and reverse with its
    sign; RF cavities are phased on the reference whatever
    the sign, and space charge is repulsive for either sign.

    Dipoles are RF-Track's own sector bends, tracked with the rest of the
    line in one Lattice: a DPH is SBend(L, angle, P/q), and a DPW a thin SBend that
    carries FELsim's edge followed by a TransferLine that carries the
    wedge's length (``_wedge_elements``). Every element of the line, the
    dipoles included, reproduces FELsim's transfer matrix; off momentum the
    dipole body focuses with K = h0 h, where FELsim's chromatic dipole has
    h^2. With physical apertures a dipole is cut as FELsim's
    ``apply_aperture`` cuts it, in y only at its own pole_gap/2
    (``_dipole_aperture``); every other aperture is RF-Track's circle
    (an ellipse when ax != ay).

    Examples
    --------
    Basic usage:
        >>> sim = RFTrackAdapter(beam_energy=45.0)
        >>> sim.set_beamline(elements)
        >>> particles = sim.generate_particles(1000)
        >>> result = sim.simulate(particles)

    With space charge (a bunch charge is required):
        >>> sim = RFTrackAdapter(space_charge=True, bunch_charge_nc=0.1)
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

    # Physical apertures for UH MkV FEL beamline [m]
    QUAD_HALF_APERTURE = 0.0135        # 27 mm bore / 2
    DIPOLE_HALF_GAP = None             # every dipole's, when set; else each its pole_gap / 2
    DIPOLE_HALF_WIDTH = None           # every dipole's, when set; else no cut in x, as FELsim
    WEDGE_POLE_GAP = 0.014478          # FELsim's dipole_wedge default pole gap
    BEAM_PIPE_RADIUS = 0.0127          # 1" beam pipe

    # Length [m] of the thin SBend that carries a DPW's edge: RF-Track's SBend
    # of zero length is the identity whatever its edge angles
    THIN_WEDGE_LENGTH = 1e-12

    # Longest space-charge step [m]: a thick element of length L takes
    # max(1, ceil(L / sc_ds_max)) steps
    SC_DS_MAX_DEFAULT = 0.05

    # Columns of every bunch the adapter builds and reads back
    BUNCH_COLUMNS = '%x %xp %y %yp %t %P %m %Q %N %id'
    REFERENCE_ID = -1

    # rft.cvar.SC_engine holds a raw pointer: keep the installed engine alive
    # after the adapter that installed it has gone
    _installed_sc_engine = None

    def __init__(self,
                 lattice_path: Optional[str] = None,
                 excel_path: Optional[str] = None,
                 mode: str = 'particle_tracking',
                 space_charge: bool = False,
                 sc_mesh: Optional[tuple] = None,
                 beam_energy: float = 45.0,
                 particle_mass: Optional[float] = None,
                 particle_charge: float = -1.0,
                 aperture: Optional[float] = 0.05,
                 G_quad: Optional[float] = None,
                 dipole_slices: Optional[int] = None,
                 rf_frequency: float = None,
                 debug: bool = None,
                 bunch_charge_nc: float = 0.0,
                 sc_ds_max: Optional[float] = SC_DS_MAX_DEFAULT):
        """
        Initialise RF-Track adapter.

        Parameters
        ----------
        lattice_path : str, optional
            Path to lattice file (Excel, JSON, or YAML)
        excel_path : str, optional
            Backward-compatible alias for lattice_path
        mode : str
            Simulation mode ('particle_tracking' only for RF-Track)
        space_charge : bool
            Enable space charge calculation (needs bunch_charge_nc > 0)
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
        dipole_slices : int, optional
            Deprecated and without effect (a FutureWarning when given):
            dipoles are RF-Track SBends. It set the analytical dipole
            correction that stood in for SBend, which a bunch built with its
            arguments in the wrong order had made look broken
            (docs/rftrack-bunch6d.md).
        rf_frequency : float, optional
            RF frequency in Hz for FELsim ↔ RF-Track coord5 conversion.
            Default: 2856 MHz (UH FEL S-band linac)
        debug : bool, optional
            Enable debug logging
        bunch_charge_nc : float
            Bunch charge magnitude in nC (the sign is particle_charge's),
            shared equally by the macro-particles handed to each tracking
            call. Default 0, with which space charge raises rather than
            running with no charge.
        sc_ds_max : float or None
            Longest space-charge step in metres: every element longer than
            the wedges' thin SBends takes max(1, ceil(L / sc_ds_max)) steps.
            Default 0.05 m: on a 0.5 ps, 60 pC bunch through the UH line it
            agrees with 0.0125 m to 0.22 % (eps_x, 1000 particles) and
            0.24 % (beta_x, 10000 particles) in the beam moments, where one
            step per element, the behaviour None restores, overstated the
            space-charge changes up to twofold. An integer ``sc_nsteps`` set on the
            adapter overrides both with that many steps per element.
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
        self.particle_charge = self._charged(particle_charge)
        self.beam_energy = beam_energy
        self._update_relativistic_params()

        # Default aperture for elements
        self.default_aperture = aperture

        # Quadrupole gradient calibration constant (T/A/m)
        # Default is for UH FEL quadrupoles (2.694 T/A/m)
        self.G_quad = G_quad if G_quad is not None else PhysicalConstants.G_quad_default

        self.dipole_slices = dipole_slices
        self._warned_energy_change = False
        self._warned_whole_bunch = False
        self._warned_wedges = set()

        # Space charge configuration
        self._set_bunch_charge(bunch_charge_nc)
        if space_charge:
            self._require_bunch_charge()
        self.space_charge_enabled = space_charge
        self.sc_mesh = sc_mesh or (32, 32, 64)
        self.sc_ds_max = sc_ds_max
        self.sc_nsteps = None  # SC steps per element; overrides sc_ds_max when set
        self._space_charge_effect = None

        # RF frequency for FELsim ↔ RF-Track coord5 conversion
        self._rf_frequency = rf_frequency if rf_frequency is not None else PhysicalConstants.f_RF_default
        self._t_scale = PhysicalConstants.C / self._rf_frequency  # c/f [m] ≈ 0.105

        # Physical aperture mode
        self._physical_apertures = False

        # RF-Track native objects
        self._lattice: Optional[rft.Lattice] = None
        self._bunch: Optional[rft.Bunch6d] = None
        self._lattice_index = np.empty(0, dtype=int)

        # NOTE: _element_type_map removed (was unused dead code).
        # Element dispatch is handled by _convert_element_to_native().

        # Load beamline if provided
        path = lattice_path or excel_path
        self.lattice_path = path
        self.excel_path = excel_path  # backward compat
        if path:
            self._load_lattice(path)

    @property
    def dipole_slices(self):
        """Deprecated and without effect; see the constructor."""
        return self._dipole_slices

    @dipole_slices.setter
    def dipole_slices(self, value):
        if value is not None:
            warnings.warn("RFTrackAdapter: dipole_slices is deprecated and has no effect; "
                          "dipoles are RF-Track SBends", FutureWarning, stacklevel=2)
        self._dipole_slices = value

    def _update_relativistic_params(self):
        """Update relativistic parameters from beam energy."""
        if self.beam_energy <= 0:
            raise ValueError(f"beam_energy must be positive, got {self.beam_energy}")
        self._gamma = 1 + self.beam_energy / self.particle_mass
        self._beta = np.sqrt(1 - 1/self._gamma**2)
        self._Pc = self._gamma * self._beta * self.particle_mass  # MeV/c

    def simulate(self,
                 particles: Optional[np.ndarray] = None,
                 mode: Optional[SimulationMode] = None,
                 n_initial: Optional[int] = None) -> SimulationResult:
        """
        Run RF-Track particle tracking simulation.

        Parameters
        ----------
        particles : ndarray (N, 6)
            Initial distribution in FELsim coordinates:
            [x(mm), x'(mrad), y(mm), y'(mrad), ΔToF/T_RF×10³, ΔK/K₀×10³]
        mode : SimulationMode, optional
            Ignored (RF-Track only supports particle tracking)
        n_initial : int, optional
            Number of macro-particles that share the bunch charge (default:
            the particles given), so that the survivors of an upstream
            section keep their charge; see ``track_elements``

        Returns
        -------
        SimulationResult
            Final particles (the survivors, coord5 relative to the reference
            particle; RF-Track moves the last survivor into a lost particle's
            row, so they are not in input order), Twiss parameters, and
            metadata, including 'survivor_ids' (the input row of each final
            particle), 'lost_particles' (see ``_track_rows``) and
            'lost_elements' (for each of its rows, the index in the beamline
            of the element that lost the particle; both None when nothing is
            lost)
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
        rows = self._bunch_rows(particles_rftrack, n_initial)
        self._bunch = self._make_bunch(rows[1:])

        self.logger.debug(
            f"Created bunch: {self._bunch.size()} particles, "
            f"Pc={self._Pc:.2f} MeV/c"
        )

        self.logger.info(
            f"Tracking {particles.shape[0]} particles through "
            f"{self._lattice.size()} elements (L={self._lattice.get_length():.3f} m)"
        )

        lost = []
        rows = self._track_rows(self._lattice, rows, self._lattice_index, lost)
        # RF-Track removes lost particles from the bunch; count survivors
        n_good = rows.shape[0] - 1
        lost = np.vstack(lost) if lost else None
        n_lost = particles.shape[0] - n_good

        if n_good > 0:
            final_particles = self.transform_coordinates(
                self._relative_to_reference(rows), CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM
            )
        else:
            final_particles = np.empty((0, 6))
        twiss = self._calculate_twiss(final_particles) if final_particles.shape[0] > 0 else {}

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
                'bunch_charge_nc': self.bunch_charge_C * 1e9,
                'sc_ds_max': self.sc_ds_max,
                'sc_nsteps': self.sc_nsteps,
                'particle_mass_mev': self.particle_mass,
                'particle_charge': self.particle_charge,
                'lattice_length': self._lattice.get_length(),
                'survivor_ids': rows[1:, 9].astype(int),
                'lost_particles': lost[:, :11] if lost is not None else None,
                'lost_elements': lost[:, 11].astype(int) if lost is not None else None,
            }
        )

    def _make_bunch(self, rows):
        """The one place the adapter builds an RF-Track bunch.

        ``rows`` is (n, 10), one macro-particle per row, in RF-Track's
        per-particle form [X XP Y YP T P MASS Q N ID] (reference manual 2.2):
        mm, mrad, mm, mrad, mm/c, MeV/c, MeV/c^2, charge in e, real particles
        per macro-particle, id. An empty array gives an empty bunch.
        """
        return rft.Bunch6d(np.asarray(rows, dtype=float).reshape(-1, 10))

    def _bunch_rows(self, ps, n_initial=None):
        """Rows for ``_make_bunch`` from an (n, 6) RF-Track phase space.

        Row 0 is the reference particle: on axis at the nominal momentum,
        t = 0, no weight (N = 0) and id REFERENCE_ID. Rows 1..n are the
        particles, with this adapter's mass and charge, the bunch charge
        shared equally among n_initial of them (default n; N = bunch
        charge / e / n_initial, so survivors keep their weight) and ids
        0..n-1 (the input row index). Weights and ids survive every rebuild
        of the bunch, and RF-Track keeps a surviving first particle first,
        so the reference stays in row 0.
        """
        ps = np.asarray(ps, dtype=float).reshape(-1, 6)
        n = ps.shape[0]
        n0 = n if n_initial is None else n_initial
        if n0 < n:
            raise ValueError(f"n_initial = {n_initial} is fewer than the {n} particles given")
        weight = self.bunch_charge_C / PhysicalConstants.Q / n0 if n0 else 0.0
        X = np.vstack([[0.0, 0.0, 0.0, 0.0, 0.0, self._Pc], ps])
        return np.column_stack([X, np.full(n + 1, self.particle_mass),
                                np.full(n + 1, float(self.particle_charge)),
                                np.r_[0.0, np.full(n, weight)],
                                np.r_[self.REFERENCE_ID, np.arange(n, dtype=float)]])

    def _new_lattice(self, elements, first=0):
        """RF-Track lattice of consecutive beamline elements, the first of
        them beamline[first], and the beamline index of the element each
        native element comes from (for the lost-particle table)."""
        lat, index = rft.Lattice(), []
        for i, elem in enumerate(elements, first):
            native = self._convert_element_to_native(elem, i)
            for ne in (native if isinstance(native, list) else [native]):
                lat.append(ne)
                index.append(i)
        return lat, np.asarray(index, dtype=int)

    def _track_rows(self, lat, rows, index=None, lost=None):
        """Track bunch rows (reference particle first) through an RF-Track
        lattice; the reference and the survivors' rows out.

        Particles lost in the lattice are appended to ``lost`` as RF-Track
        reports them, rows x xp y yp t P S m Q N id, except that S is the
        start in the adapter's line of the element in which the particle was
        lost, and followed by that element's beamline index, which also
        tells a zero-length element from the one after it (native element i
        comes from beamline[index[i]]; see ``_loss_elements``). The
        coordinates are where RF-Track caught the particle: at the entrance
        of a drift or sextupole, at the entrance or exit of a quadrupole,
        solenoid or sector bend; t runs on
        the clock the reference particle starts at 0. The reference particle
        is left out.

        Every lattice the adapter tracks goes through here. With space charge
        on, each of its elements gets its space-charge steps (``_sc_steps``)
        and this adapter's engine is installed just before tracking
        (RF-Track's engine is global); the
        reference is then tracked in a bunch of its own, since in the beam's
        bunch it would feel the beam's field and stretch the PIC mesh to
        reach it. RF elements are phased on the reference particle as it
        enters the lattice: RF-Track would otherwise autophase on the first
        track only, on whatever particle comes first, and keep that phase.
        Rows holding the reference alone track the reference alone.
        """
        if self.space_charge_enabled:
            self._require_bunch_charge()
            if self._space_charge_effect is None:
                self._setup_space_charge()
            self._enable_space_charge_on_lattice(lat)
            self._install_space_charge()
        rf_elements = lat.get_rf_elements()
        if len(rf_elements):
            for rf in rf_elements:
                rf.unset_t0()
            lat.autophase(self._make_bunch(rows[:1]))
        if self.space_charge_enabled:
            out, table = self._track_bunch(lat, rows[1:])
            ref, ref_table = self._track_bunch(lat, rows[:1])
            out, table = np.vstack([ref, out]), np.vstack([ref_table, table])
        else:
            out, table = self._track_bunch(lat, rows)
        is_ref = table[:, 10] == self.REFERENCE_ID
        if is_ref.any():
            out = np.vstack([self._reference_after_loss(lat, rows[0], table[is_ref]), out])
        table = table[~is_ref]
        if lost is not None and table.shape[0]:
            owner = index[self._loss_elements(lat, table)]
            table[:, 6] = self._element_starts()[owner]
            lost.append(np.column_stack([table, owner]))
        return out

    def _element_starts(self):
        """Where each beamline element starts [m]."""
        return np.cumsum([0.0] + [e.length for e in self.beamline])[:-1]

    def _track_bunch(self, lat, rows):
        """Rows of the survivors of one RF-Track track, and the lattice's
        lost-particle table (rows x xp y yp t P S m Q N id)."""
        if rows.shape[0] == 0:
            return rows, np.empty((0, 11))
        # RF-Track's "beam's first particle was lost" warning on a lossy
        # track is spurious here: the columns are absolute and the reference
        # stays in row 0.
        out = np.asarray(lat.track(self._make_bunch(rows)).get_phase_space(self.BUNCH_COLUMNS),
                         dtype=float).reshape(-1, 10)
        if out.shape[0] == rows.shape[0]:
            return out, np.empty((0, 11))
        return out, np.asarray(lat.get_lost_particles(), dtype=float).reshape(-1, 11)

    def _loss_elements(self, lat, table):
        """Index in ``lat`` of the element in which each row of its
        lost-particle table was lost.

        RF-Track reports a loss where it caught the particle: drifts and
        sextupoles check the whole path and report at their entrance,
        quadrupoles, solenoids and sector bends check their entrance and
        exit. So a row at a boundary belongs to
        the first of the elements ending there whose aperture the particle
        lies outside (one that checks its whole path would have caught it at
        its entrance), and otherwise to the element starting there.
        """
        lengths = np.array([lat[i].get_length() for i in range(lat.size())])
        ends = np.cumsum(lengths)
        tol = 1e-9      # [m] longer than a wedge's thin SBend, shorter than any other element
        owner = np.searchsorted(ends - lengths, table[:, 6] + tol, side='right') - 1
        owner = np.clip(owner, 0, lat.size() - 1)
        for r, (x, y, s) in enumerate(table[:, [0, 2, 6]]):
            first = owner[r]
            while first > 0 and abs(ends[first - 1] - s) <= tol:
                first -= 1
            owner[r] = next((k for k in range(first, owner[r])
                             if self._outside_aperture(lat[k], x, y)), owner[r])
        return owner

    @staticmethod
    def _outside_aperture(elem, x, y):
        """Whether (x, y) [mm] lies outside an RF-Track element's aperture."""
        shape = elem.get_aperture_shape()
        if shape == 'none':
            return False
        ax, ay = 1e3 * elem.get_aperture_x(), 1e3 * elem.get_aperture_y()
        if shape == 'rectangular':
            return abs(x) > ax or abs(y) > ay
        return (x / ax)**2 + (y / ay)**2 > 1.0

    @staticmethod
    def _relative_to_reference(rows):
        """The particles' (n, 6) RF-Track phase space from bunch rows, with t
        taken relative to the reference particle's arrival time."""
        ps = rows[1:, :6].copy()
        ps[:, 4] -= rows[0, 4]
        return ps

    def _reference_after_loss(self, lat, ref, hit):
        """Reference row at the exit of a lattice in which it was lost: from where it was
        lost (its row ``hit`` of the lost-particle table) it continues on
        axis, at the momentum it had there. Tracked with the beam, its row has
        gone to the last survivor, so the survivors leave input order."""
        t, P, S = hit[0, 4], hit[0, 5], hit[0, 6]
        beta = P / np.hypot(P, self.particle_mass)
        self.logger.warning(
            f"Reference particle lost at s = {S:.4f} m of a {lat.get_length():.4f} m "
            f"lattice; continuing it on axis to the exit"
        )
        out = ref.copy()
        out[:4] = 0.0
        out[4] = t + (lat.get_length() - S) * 1e3 / beta
        out[5] = P
        return out[None, :]

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

        opt_vars = {}
        for i, (idx, (var_name, param_name, transform)) in enumerate(variables.items()):
            value = transform(result.x[i])
            opt_vars[var_name] = value
            self._modify_element(idx, **{param_name: value})
        self._build_lattice()
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

    def _build_rf_cavity(self, length: float, params: dict) -> Any:
        """Build an RF-Track RF cavity from generic RFC parameters.

        Supports structure_type ∈ {'TW', 'SW', 'RFCA'}. For TW/SW the
        peak on-axis gradient is passed as a single Fourier coefficient
        (constant-gradient approximation); n_cells is derived from
        length and phase advance if not provided, assuming β_wave = 1.

        Expected params:
            frequency_hz         (required)
            phase_deg            (default 0.0 — on-crest under autophase)
            gradient_mv_per_m OR voltage_mv  (required)
            structure_type       (default 'TW')
            phase_advance_deg    (default 120.0)
            n_cells              (default: auto from length + phase advance)
        """
        freq = params.get('frequency_hz')
        if freq is None:
            raise ValueError("RF_CAVITY element missing 'frequency_hz'")
        phase_deg = float(params.get('phase_deg', 0.0))
        structure_type = str(params.get('structure_type', 'TW')).upper()
        phi_adv = float(params.get('phase_advance_deg', 120.0)) * np.pi / 180.0

        # Resolve gradient
        gradient_vpm = params.get('gradient_mv_per_m')
        voltage_mv = params.get('voltage_mv')
        if gradient_vpm is not None:
            e0_vpm = float(gradient_vpm) * 1e6
        elif voltage_mv is not None:
            if length <= 0:
                raise ValueError("RF_CAVITY: cannot derive gradient from voltage_mv with zero length")
            e0_vpm = float(voltage_mv) * 1e6 / length
        else:
            raise ValueError(
                "RF_CAVITY: provide 'gradient_mv_per_m' or 'voltage_mv'"
            )

        # Synchronous cell length assuming β_wave = 1
        c_m_s = 299792458.0
        l_cell_sync = c_m_s * phi_adv / (2.0 * np.pi * float(freq))

        if structure_type == 'TW':
            n_cells = params.get('n_cells')
            if n_cells is None:
                n_cells = max(1.0, length / l_cell_sync)
            else:
                # User-supplied n_cells: warn if it yields a length
                # significantly different from the element's length_m.
                actual_L = float(n_cells) * l_cell_sync
                if length > 0 and abs(actual_L - length) / length > 0.01:
                    self.logger.warning(
                        f"RFC TW: n_cells={n_cells} gives L={actual_L:.4f} m, "
                        f"but length_m={length} m (delta "
                        f"{(actual_L - length) * 1e3:+.1f} mm). Using "
                        f"n_cells-derived length."
                    )
            elem = rft.TW_Structure(
                float(e0_vpm), 0, float(freq), float(phi_adv), float(n_cells)
            )
            elem.set_phid(phase_deg)
        elif structure_type == 'SW':
            # SW_Structure(coefficients, freq, cell_length, n_cells).
            # Accept explicit cell_length_m or derive from phase advance
            # assuming β_wave = 1 (same convention as TW).
            cell_length = params.get('cell_length_m', l_cell_sync)
            n_cells = params.get('n_cells')
            if n_cells is None:
                n_cells = length / float(cell_length)
            elem = rft.SW_Structure(
                [float(e0_vpm)], float(freq),
                float(cell_length), float(n_cells)
            )
            elem.set_phid(phase_deg)
        else:  # RFCA — no direct RF-Track equivalent; approximate as TW
            self.logger.warning(
                "RFCA structure_type approximated as TW_Structure; "
                "true lumped-cavity behavior requires the elegant adapter."
            )
            n_cells = max(1.0, length / l_cell_sync)
            elem = rft.TW_Structure(
                float(e0_vpm), 0, float(freq), float(phi_adv), float(n_cells)
            )
            elem.set_phid(phase_deg)

        if 'name' in params:
            elem.set_name(str(params['name']))
        return elem

    def _convert_element_to_native(self, element: BeamlineElement,
                                   index: Optional[int] = None) -> Any:
        """
        Convert generic BeamlineElement to RF-Track element.

        Parameters
        ----------
        element : BeamlineElement
            Generic element representation
        index : int, optional
            Its index in the beamline: a DPW without dipole_length and
            dipole_angle takes its bend from the DPH beside it there
            (``wedge_bend``)

        Returns
        -------
        RF-Track element object (rft.Drift, rft.Quadrupole, etc.), or a list
        of them for a DPW
        """
        elem_type = element.element_type.upper()
        params = dict(element.parameters)
        length = element.length
        shape = 'circular'

        if elem_type == 'DRIFT':
            elem = rft.Drift(length)
            ap_x = ap_y = self.BEAM_PIPE_RADIUS if self._physical_apertures else self.default_aperture

        elif elem_type in ['QUAD_F', 'QPF']:
            elem = rft.Quadrupole()
            elem.set_length(length)
            k1 = self._current_to_k1(params.get('current', 0.0), length, focusing=True)
            # RF-Track's strength is S = (P/q) k1 L [MV/c/m] (manual 4.2.2);
            # set_K1L(P/q, k1 L) sets it for the reference momentum, and each
            # particle then sees k1 Pc/P (test/test_rftrack_quad_convention.py)
            elem.set_K1L(self._Pc / self.particle_charge, k1 * length)
            r = self.QUAD_HALF_APERTURE if self._physical_apertures else self.default_aperture
            ap_x = ap_y = r

        elif elem_type in ['QUAD_D', 'QPD']:
            elem = rft.Quadrupole()
            elem.set_length(length)
            k1 = self._current_to_k1(params.get('current', 0.0), length, focusing=False)
            elem.set_K1L(self._Pc / self.particle_charge, k1 * length)
            r = self.QUAD_HALF_APERTURE if self._physical_apertures else self.default_aperture
            ap_x = ap_y = r

        elif elem_type in ['DIPOLE_WEDGE', 'DPW']:
            wedge = self._wedge_elements(element, index)
            for elem in wedge:
                self._set_aperture(elem, *self._dipole_aperture(element))
                if 'name' in params:
                    elem.set_name(params['name'])
            return wedge

        elif elem_type in ['DIPOLE', 'DPH']:
            # FELsim's angle, signed as RF-Track's with P/q signed as the
            # charge: the reference particle follows the bend
            elem = rft.SBend(length, np.radians(params.get('angle', 0.0)),
                             self._Pc / self.particle_charge)
            ap_x, ap_y, shape = self._dipole_aperture(element)

        elif elem_type == 'SOLENOID':
            # Solenoid(L, B0 [T], R): in a Lattice, R only sets the aperture.
            # At B0 = 0 RF-Track's solenoid returns NaN, so that is a drift.
            ap_x = ap_y = self.default_aperture
            field = params.get('field', 0.0)
            if field == 0:
                elem = rft.Drift(length)
            else:
                elem = rft.Solenoid(length, field, ap_x if ap_x is not None else 1.0)
                if ap_x is None:
                    elem.set_aperture_shape('none')

        elif elem_type == 'SEXTUPOLE':
            # 'strength' is the integrated normalised k2 L [m^-2]
            elem = rft.Sextupole()
            elem.set_length(length)
            elem.set_K2L(self._Pc / self.particle_charge, params.get('strength', 0.0))
            ap_x = ap_y = self.default_aperture

        elif elem_type in ('RF_CAVITY', 'RFC'):
            elem = self._build_rf_cavity(length, params)
            ap_x = ap_y = self.default_aperture

        elif elem_type in ('ALPHA_MAGNET', 'AMG'):
            raise NotImplementedError(
                "RF-Track adapter has no alpha magnet element; track the "
                "section containing it with the COSY or FELsim backend"
            )

        else:
            self.logger.warning(f"Unknown element type '{elem_type}', using drift")
            elem = rft.Drift(length)
            ap_x = ap_y = self.BEAM_PIPE_RADIUS if self._physical_apertures else self.default_aperture

        if hasattr(elem, 'set_aperture'):
            self._set_aperture(elem, ap_x, ap_y, shape)

        if hasattr(elem, 'set_name') and 'name' in params:
            elem.set_name(params['name'])

        return elem

    @staticmethod
    def _set_aperture(elem, ap_x, ap_y, shape='circular'):
        """aperture=None leaves RF-Track's native 'none' shape: lossless."""
        if ap_x is not None and ap_y is not None:
            elem.set_aperture(ap_x, ap_y, shape)

    def _dipole_aperture(self, element):
        """Half-width and half-gap [m] and shape of a DPH's or DPW's aperture.

        With physical apertures, FELsim's: a dipole is cut in y only, at its
        pole_gap/2 (``beamline.lattice.apply_aperture`` tests |x| and |y|
        separately, and a dipole has no aperture_x), so a rectangle without
        bound in x. A DPH without a pole gap is not cut, and a DPW without
        one has dipole_wedge's default gap, WEDGE_POLE_GAP. DIPOLE_HALF_GAP
        and DIPOLE_HALF_WIDTH, when set, apply to every dipole. Otherwise the
        default aperture, the circle of every other element.
        """
        if not self._physical_apertures:
            return self.default_aperture, self.default_aperture, 'circular'
        gap = element.parameters.get('pole_gap')
        if gap is None and element.element_type.upper() in ('DIPOLE_WEDGE', 'DPW'):
            gap = self.WEDGE_POLE_GAP
        half_gap = self.DIPOLE_HALF_GAP if self.DIPOLE_HALF_GAP is not None else (
            gap / 2 if gap else None)
        if half_gap is None and self.DIPOLE_HALF_WIDTH is None:
            return None, None, 'rectangular'
        return (np.inf if self.DIPOLE_HALF_WIDTH is None else self.DIPOLE_HALF_WIDTH,
                np.inf if half_gap is None else half_gap, 'rectangular')

    def _wedge_elements(self, element, index=None):
        """The RF-Track elements of a DPW: a thin SBend that carries its
        edge, then a TransferLine over its length.

        FELsim's DPW is a thin edge, M21 = tan(eta)/R and M43 =
        -tan(eta - phi)/R with the triangle-model fringe angle phi, that
        takes the wedge's length only in s and in the time of flight (M12 =
        M34 = 0, M56 of a drift). The SBend, THIN_WEDGE_LENGTH long with
        curvature h = 1/|R| (FELsim's unsigned convention: the edge does not
        depend on the sign of the bend), has E1 = eta, E2 = 0, and hgap and
        fint set so that RF-Track's fringe terms on its two faces add up to
        phi (``_thin_wedge_fringe``); its edges scale as P0/P per particle,
        as FELsim's chromatic DPW. The TransferLine is the transverse
        identity (beta 1, alpha 0, no phase advance) and a drift in time.
        Without a bend (dipole_angle 0, R infinite) there is no edge.
        """
        p = element.parameters
        h = self._wedge_curvature(element, index)
        out = []
        if h:
            eta = np.radians(element.edge_angle)
            sb = rft.SBend(self.THIN_WEDGE_LENGTH, h * self.THIN_WEDGE_LENGTH,
                           self._Pc / self.particle_charge, eta, 0.0)
            k = self._thin_wedge_fringe(eta, h, element.length)
            if k:
                gap = p.get('pole_gap') or 0.0
                hgap = gap / 2 if gap > 0 else element.length
                sb.set_hgap(hgap)
                sb.set_fint(k / (2 * hgap * h))
            out.append(sb)
        if element.length > 0:
            twiss = np.array([[0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                              [element.length, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0]])
            out.append(rft.TransferLine(twiss, self._Pc))
        return out or [rft.Drift(0.0)]

    def _wedge_curvature(self, element, index=None):
        """1/|R| [1/m] of a DPW: from its own dipole_length and dipole_angle,
        as in FELsim (0 for a zero dipole_angle, an infinite R), or without
        them from the DPH beside it in the beamline (``wedge_bend``; the
        adapter warns once per wedge)."""
        elements, i = (self.beamline, index) if index is not None else ([element], 0)
        bend = wedge_bend(elements, i)
        if bend is None:
            raise ValueError(
                "RF-Track: a DIPOLE_WEDGE needs dipole_length and dipole_angle, "
                "or a DIPOLE next to it")
        length, angle, source = bend
        if source is not None and index not in self._warned_wedges:
            self._warned_wedges.add(index)
            self.logger.warning(
                f"RF-Track: DIPOLE_WEDGE {index} has no dipole_length and dipole_angle; "
                f"it takes those of the DIPOLE at {source} ({angle} deg over {length} m)")
        if not angle:
            return 0.0
        if not length > 0:
            raise ValueError(f"RF-Track: DIPOLE_WEDGE {index} bends by {angle} deg over "
                             f"dipole_length {length!r}; the length must be positive")
        return abs(np.radians(angle)) / length

    @staticmethod
    def _thin_wedge_fringe(eta, h, le):
        """k = 2 fint hgap h that gives a thin SBend with faces E1 = eta and
        E2 = 0 FELsim's vertical edge, M43 = -h tan(eta - phi), phi = (le/6)
        h a with a = (1 + sin^2 eta)/cos eta. RF-Track 2.5.5 gives each face
        of an SBend M43 = -h tan(E - psi), psi = k (1 + sin^2 E)/cos E, so k
        solves tan(eta - k a) - tan(k) = tan(eta - phi): about phi/(1 + a),
        where one face alone would take phi/a."""
        a = (1 + np.sin(eta)**2) / np.cos(eta)
        phi = le * h * a / 6
        if phi == 0:
            return 0.0
        target = np.tan(eta - phi)
        return brentq(lambda k: np.tan(eta - k * a) - np.tan(k) - target, 0.0, 4 * phi / (1 + a))

    def transform_coordinates(self,
                              particles: np.ndarray,
                              from_system: CoordinateSystem,
                              to_system: CoordinateSystem) -> np.ndarray:
        """
        Transform particle coordinates between systems.

        Coordinate systems:
        - FELSIM: [x(mm), x'(mrad), y(mm), y'(mrad), ΔToF/T_RF×10³, ΔK/K₀×10³]
        - RFTRACK: [x(mm), x'(mrad), y(mm), y'(mrad), ct(mm/c), P(MeV/c)]

        Coord5 conversion: ct_mm = coord5 × c/f_RF (the ×10³ and m→mm cancel).
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
        particles = np.asarray(particles, dtype=float)
        if from_system == to_system:
            return particles.copy()

        result = np.zeros_like(particles)

        if from_system == CoordinateSystem.FELSIM and to_system == CoordinateSystem.RFTRACK:
            # FELsim: [x(mm), x'(mrad), y(mm), y'(mrad), ΔToF/T×10³, ΔK/K₀×10³]
            # RF-Track: [x(mm), x'(mrad), y(mm), y'(mrad), t(mm/c), P(MeV/c)]
            # Transverse coordinates: same units, no conversion needed
            result[:, 0:4] = particles[:, 0:4]
            # Longitudinal: FELsim ΔToF/T_RF×10³ → RF-Track ct [mm]
            # ct_mm = coord5 × c/f_RF (the ×10³ and m→mm cancel)
            result[:, 4] = particles[:, 4] * self._t_scale
            # FELsim coord6 = ΔK/K₀ × 10³ → RF-Track P [MeV/c] (exact)
            K = self.beam_energy * (1.0 + particles[:, 5] * 1e-3)
            E = K + self.particle_mass
            result[:, 5] = np.sqrt(E**2 - self.particle_mass**2)

        elif from_system == CoordinateSystem.RFTRACK and to_system == CoordinateSystem.FELSIM:
            # RF-Track: [x(mm), x'(mrad), y(mm), y'(mrad), t(mm/c), P(MeV/c)]
            # FELsim: [x(mm), x'(mrad), y(mm), y'(mrad), ΔToF/T×10³, ΔK/K₀×10³]
            # Transverse coordinates: same units
            result[:, 0:4] = particles[:, 0:4]
            # Longitudinal: RF-Track ct [mm] → FELsim ΔToF/T_RF×10³
            result[:, 4] = particles[:, 4] / self._t_scale
            # RF-Track P [MeV/c] → FELsim coord6 = ΔK/K₀ × 10³ (exact)
            K = np.sqrt(particles[:, 5]**2 + self.particle_mass**2) - self.particle_mass
            result[:, 5] = (K / self.beam_energy - 1.0) * 1e3

        elif from_system == CoordinateSystem.COSY and to_system == CoordinateSystem.RFTRACK:
            # COSY uses m and rad; RF-Track Bunch6d uses mm and mrad
            result[:, 0:4] = particles[:, 0:4] * 1e3  # m/rad → mm/mrad
            # COSY l → RF-Track ct [mm]:
            # COSY l = -Δt · v₀ · γ/(1+γ), so Δt = -l·(1+γ)/(v₀·γ)
            # ct [mm] = c·Δt·10³ = -l·(1+γ)/(β·γ)·10³
            result[:, 4] = -particles[:, 4] * (1 + self._gamma) / (self._beta * self._gamma) * 1e3
            # COSY δ = ΔK/K₀ → RF-Track P [MeV/c] (exact)
            E0 = self.beam_energy + self.particle_mass  # total energy [MeV]
            E = E0 + self.beam_energy * particles[:, 5]  # E₀ + K₀δ
            result[:, 5] = np.sqrt(E**2 - self.particle_mass**2)  # P [MeV/c]

        elif from_system == CoordinateSystem.RFTRACK and to_system == CoordinateSystem.COSY:
            # RF-Track uses mm/mrad; COSY uses m/rad
            result[:, 0:4] = particles[:, 0:4] * 1e-3  # mm/mrad → m/rad
            # RF-Track ct [mm] → COSY l [m]:
            # l = -ct·β·γ/((1+γ)·10³)
            result[:, 4] = -particles[:, 4] * self._beta * self._gamma / ((1 + self._gamma) * 1e3)
            # RF-Track P [MeV/c] → COSY δ = ΔK/K₀ (exact)
            K = np.sqrt(particles[:, 5]**2 + self.particle_mass**2) - self.particle_mass
            result[:, 5] = K / self.beam_energy - 1.0

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
        self._warned_wedges = set()
        self._build_lattice()

    def track_elements(self, ps_rftrack, start_idx, end_idx, n_initial=None):
        """Track through beamline[start_idx:end_idx], in one RF-Track lattice.

        Input/output are RF-Track coordinates (N, 6) of the survivors, in
        RF-Track's order (see ``simulate``), with t relative to the reference
        particle, so the output of one range is the input of the next:
        track_elements(a, b) followed by track_elements(b, c, n_initial=n)
        equals track_elements(a, c) for n particles handed in at a, cuts
        inside a DPW-DPH-DPW triplet included. The reference enters the range
        as it would on the whole line (see ``_reference_at``), and n_initial
        (default: the particles handed in) is the number of macro-particles
        that share the bunch charge, which matters with space charge once
        particles have been lost: with space charge on and no n_initial for a
        range that starts after the first element, the adapter warns once
        that it takes the particles as the whole bunch.
        """
        if (self.space_charge_enabled and n_initial is None and start_idx > 0
                and not self._warned_whole_bunch):
            self._warned_whole_bunch = True
            self.logger.warning(
                "RF-Track track_elements with space charge: the particles handed in are taken "
                "as the whole bunch; pass n_initial when they are the survivors of an earlier range"
            )
        rows = self._bunch_rows(ps_rftrack, n_initial)
        rows[0] = self._reference_at(start_idx, end_idx)
        rows[1:, 4] += rows[0, 4]
        lat, index = self._new_lattice(self.beamline[start_idx:end_idx], start_idx)
        return self._relative_to_reference(self._track_rows(lat, rows, index))

    def _reference_at(self, start_idx, end_idx):
        """Row of the reference particle entering beamline[start_idx], for a
        range that ends at end_idx.

        RF cavities change the reference's momentum and are phased on its
        arrival time, and RF-Track's space-charge kicks change with the origin
        of time (at the level of the PIC noise). With space charge on, or when
        the line up to end_idx has an RF cavity, the reference is tracked to start_idx from the
        start of the line, so that it arrives with the momentum, orbit and
        time it has there on the whole line. Otherwise it would arrive on
        axis at the nominal momentum, and it starts at t = 0, which shifts
        every t alike.
        """
        rows = self._bunch_rows(np.empty((0, 6)))
        if self.space_charge_enabled or any(self._moves_reference(e)
                                            for e in self.beamline[:end_idx]):
            lat, index = self._new_lattice(self.beamline[:start_idx])
            rows = self._track_rows(lat, rows, index)
        return rows[0]

    def _moves_reference(self, elem):
        """Whether an element acts on a particle on the design orbit at the
        nominal momentum: an RF cavity."""
        return elem.element_type.upper() in ('RF_CAVITY', 'RFC')

    def _build_lattice(self):
        """Build RF-Track lattice from beamline elements."""
        self._lattice, self._lattice_index = self._new_lattice(self.beamline)
        self._warn_energy_change()

        self.logger.debug(
            f"Built RF-Track lattice: {self._lattice.size()} elements, "
            f"L={self._lattice.get_length():.3f} m"
        )
        if self.space_charge_enabled:
            if self._space_charge_effect is None:
                self._setup_space_charge()
            else:
                self._enable_space_charge_on_lattice(self._lattice)

    def _warn_energy_change(self):
        """Warn, once, of dipoles and sextupoles downstream of an RF cavity.

        Their strengths are angles and normalised strengths, set for the
        adapter's initial momentum, so the reference accelerated by the cavity
        is mis-bent. Quadrupoles and solenoids are set from their
        fields (G I, B0) and are right at any momentum.
        """
        if self._warned_energy_change:
            return
        types = [e.element_type.upper() for e in self.beamline]
        first_rf = next((i for i, t in enumerate(types) if t in ('RF_CAVITY', 'RFC')), None)
        if first_rf is None:
            return
        after = sorted({t for t in types[first_rf:]
                        if t in ('DIPOLE', 'DPH', 'DIPOLE_WEDGE', 'DPW', 'SEXTUPOLE')})
        if after:
            self._warned_energy_change = True
            self.logger.warning(
                f"RF-Track: {', '.join(after)} after an RF cavity is set for the "
                f"initial momentum {self._Pc:.4f} MeV/c, not for the accelerated reference"
            )

    def _modify_element(self, index: int, **kwargs):
        """Modify element parameters in beamline."""
        if 0 <= index < len(self.beamline):
            for key, value in kwargs.items():
                self.beamline[index].parameters[key] = value

    def set_space_charge(self, enabled: bool, mesh: Optional[tuple] = None,
                         method: str = 'PIC', bunch_charge_nc: Optional[float] = None):
        """
        Configure space charge calculation.

        Parameters
        ----------
        enabled : bool
            Enable/disable space charge; enabling needs a bunch charge
        mesh : tuple, optional
            Mesh size (nx, ny, nz) for 3D solver
        method : str
            Space charge method: 'PIC', 'P2P', 'FreeSpace'
        bunch_charge_nc : float, optional
            Bunch charge magnitude in nC; unchanged when None
        """
        if bunch_charge_nc is not None:
            self._set_bunch_charge(bunch_charge_nc)
        if enabled:
            self._require_bunch_charge()
        self.space_charge_enabled = enabled
        if mesh:
            self.sc_mesh = mesh
        self._sc_method = method

        if enabled:
            self._setup_space_charge()
        else:
            self._space_charge_effect = None
            if self._lattice is not None:
                for i in range(self._lattice.size()):
                    self._lattice[i].set_sc_nsteps(0)

        self.logger.info(
            f"Space charge: {enabled}, method={method}, mesh={self.sc_mesh}, "
            f"bunch charge={self.bunch_charge_C * 1e9:g} nC"
        )

    def _set_bunch_charge(self, bunch_charge_nc: float):
        try:
            q = float(bunch_charge_nc)
        except (TypeError, ValueError):
            q = np.nan
        if not (np.isfinite(q) and q >= 0):
            raise ValueError(
                f"bunch_charge_nc is a finite magnitude in nC (the sign comes "
                f"from particle_charge), got {bunch_charge_nc!r}"
            )
        self.bunch_charge_C = q * 1e-9

    @staticmethod
    def _charged(particle_charge: float) -> float:
        if particle_charge == 0:
            raise ValueError(
                "particle_charge must be non-zero: solenoids, quadrupoles "
                "(set with P/q) and space charge act through it"
            )
        return particle_charge

    def _require_bunch_charge(self):
        if not self.bunch_charge_C > 0:
            raise ValueError(
                "RF-Track space charge needs the bunch charge: pass "
                "bunch_charge_nc > 0 to the constructor or set_space_charge()"
            )

    def enable_physical_apertures(self, quad_half_aperture=None,
                                   dipole_half_gap=None, dipole_half_width=None,
                                   beam_pipe_radius=None):
        """Enable per-element physical apertures for particle tracking.

        Parameters
        ----------
        quad_half_aperture : float, optional
            Quadrupole half-aperture [m]. Default: 13.5 mm (27 mm bore/2).
        dipole_half_gap : float, optional
            Vertical half-gap [m] of every dipole. Default: each dipole's
            pole_gap/2, as in FELsim (see ``_dipole_aperture``).
        dipole_half_width : float, optional
            Horizontal half-width [m] of every dipole. Default: none, as in
            FELsim.
        beam_pipe_radius : float, optional
            Default beam pipe radius [m]. Default: 12.7 mm (1" pipe).
        """
        self._physical_apertures = True
        if quad_half_aperture is not None:
            self.QUAD_HALF_APERTURE = quad_half_aperture
        if dipole_half_gap is not None:
            self.DIPOLE_HALF_GAP = dipole_half_gap
        if dipole_half_width is not None:
            self.DIPOLE_HALF_WIDTH = dipole_half_width
        if beam_pipe_radius is not None:
            self.BEAM_PIPE_RADIUS = beam_pipe_radius

        self.logger.info(
            f"Physical apertures enabled: quad={self.QUAD_HALF_APERTURE*1e3:.1f} mm, "
            f"dipole half-gap={self._mm(self.DIPOLE_HALF_GAP, 'pole_gap/2')}, "
            f"dipole half-width={self._mm(self.DIPOLE_HALF_WIDTH, 'none')}, "
            f"pipe={self.BEAM_PIPE_RADIUS*1e3:.1f} mm"
        )

        if self._lattice is not None and self.beamline:
            self._build_lattice()

    @staticmethod
    def _mm(value, unset):
        return unset if value is None else f"{value * 1e3:.2f} mm"

    def disable_physical_apertures(self):
        """Disable per-element physical apertures, revert to default."""
        self._physical_apertures = False
        if self._lattice is not None and self.beamline:
            self._build_lattice()

    def _setup_space_charge(self):
        """Create this adapter's SC engine and switch on per-element kicks.

        RF-Track 2.5.x handles space charge through:
        1. A global SC engine set via ``rft.cvar.SC_engine``
        2. Per-element activation via ``elem.set_sc_nsteps(N)``

        The engine is global RF-Track state, so ``_track_rows`` installs this
        adapter's own engine before every track; adapters in one process can
        alternate. Parallel tracking needs separate processes.
        """
        nx, ny, nz = self.sc_mesh
        method = getattr(self, '_sc_method', 'PIC')

        if method == 'P2P':
            self._space_charge_effect = rft.SpaceCharge_P2P()
        else:
            self._space_charge_effect = rft.SpaceCharge_PIC_FreeSpace(nx, ny, nz)

        self._enable_space_charge_on_lattice(self._lattice)

    def _install_space_charge(self):
        """Make this adapter's engine RF-Track's global one. The class keeps a
        reference to it, because rft.cvar.SC_engine does not: a freed engine
        left installed segfaults the next track with space-charge steps."""
        rft.cvar.SC_engine = self._space_charge_effect
        RFTrackAdapter._installed_sc_engine = self._space_charge_effect

    @property
    def sc_ds_max(self) -> Optional[float]:
        return self._sc_ds_max

    @sc_ds_max.setter
    def sc_ds_max(self, value):
        if value is not None:
            try:
                value = np.nan if isinstance(value, (bool, np.bool_)) else float(value)
            except (TypeError, ValueError):
                value = np.nan
            if not (np.isfinite(value) and value > 0):
                raise ValueError(
                    "sc_ds_max is a positive step length in metres, or None "
                    "for one space-charge step per element")
        self._sc_ds_max = value

    @property
    def sc_nsteps(self) -> Optional[int]:
        return self._sc_nsteps

    @sc_nsteps.setter
    def sc_nsteps(self, value):
        if value is not None:
            try:
                if isinstance(value, (bool, np.bool_)):
                    raise TypeError
                value = operator.index(value)
            except TypeError:
                value = -1
            if value < 0:
                raise ValueError(
                    "sc_nsteps is a non-negative integer number of space-charge "
                    "steps per element, or None to follow sc_ds_max")
        self._sc_nsteps = value

    def _sc_steps(self, length: float) -> int:
        """Space-charge steps of a native element of this length: none on
        the thin SBends of the wedges (a PIC solve over 1e-12 m would be
        wasted; a zero-length element takes no kick either way), else
        sc_nsteps when set, else one per sc_ds_max, at least one."""
        if length <= self.THIN_WEDGE_LENGTH:
            return 0
        if self.sc_nsteps is not None:
            return self.sc_nsteps
        if self.sc_ds_max is None:
            return 1
        return max(1, math.ceil(length / self.sc_ds_max - 1e-9))

    def _enable_space_charge_on_lattice(self, lat):
        """Set every element's space-charge steps on a freshly built
        RF-Track lattice."""
        if not self.space_charge_enabled or lat is None:
            return
        for i in range(lat.size()):
            lat[i].set_sc_nsteps(self._sc_steps(lat[i].get_length()))

    def collect_evolution(self,
                          particles: np.ndarray,
                          checkpoint_elements: Union[str, List[int]] = 'all') -> BeamEvolution:
        """
        Collect beam evolution data at element boundaries.

        Tracks element-by-element to capture phase space at each checkpoint,
        each element in an RF-Track lattice of its own (curvilinear
        coordinates along the design orbit) built from the native elements
        ``simulate`` tracks. Samples are in FELsim coordinates with coord5
        relative to the reference particle, as in ``simulate``; once every
        particle is lost, the remaining elements report no survivors.

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

        rows = self._bunch_rows(self.transform_coordinates(
            particles, CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK
        ))

        # Track element by element; once every particle is lost the
        # remaining elements are recorded with no survivors
        s = 0.0
        n_initial = particles.shape[0]
        for idx, elem in enumerate(self.beamline):
            n_before = rows.shape[0] - 1

            if n_before > 0:
                rows = self._track_rows(self._new_lattice([elem], idx)[0], rows)

            # Update position
            s += elem.length

            n_good = rows.shape[0] - 1
            n_lost_elem = n_before - n_good

            if n_lost_elem > 0:
                self.logger.info(
                    f"Element {idx} ({elem.element_type}, s={s:.4f} m): "
                    f"lost {n_lost_elem}, remaining {n_good}/{n_initial} "
                    f"({n_good/n_initial:.1%})"
                )

            # Record element info with transmission data
            elem_params = dict(elem.parameters)
            elem_params['n_good'] = n_good
            elem_params['n_lost'] = n_lost_elem
            elem_params['transmission'] = n_good / n_initial if n_initial > 0 else 0

            evolution.elements.append(ElementInfo(
                element_type=elem.element_type,
                s_start=s - elem.length,
                s_end=s,
                length=elem.length,
                color=self._get_element_color(elem.element_type),
                index=idx,
                parameters=elem_params
            ))

            # Checkpoint if requested
            if idx in checkpoint_set and n_good > 0:
                ps_felsim = self.transform_coordinates(
                    self._relative_to_reference(rows), CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM
                )
                evolution.add_sample(s, ps_felsim, self._calculate_twiss(ps_felsim))

        evolution.total_length = s
        return evolution

    def _load_lattice(self, lattice_path: str):
        """Load beamline from Excel, JSON, or YAML lattice file."""
        import latticeLoader

        try:
            native_elements = latticeLoader.create_beamline(lattice_path)

            self.beamline = []
            for elem in native_elements:
                self.beamline.append(self._convert_element_from_native(elem))

            self._build_lattice()
            self.logger.info(f"Loaded {len(self.beamline)} elements from {lattice_path}")

        except Exception as e:
            self.logger.error(f"Failed to load beamline from {lattice_path}: {type(e).__name__}: {e}")
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
            'alphaMagnetLattice': 'ALPHA_MAGNET',
            'rfCavityLattice': 'RF_CAVITY',
        }

        elem_type = type_map.get(cls_name, cls_name.upper())

        params = {}
        if hasattr(native_elem, 'current'):
            params['current'] = native_elem.current
        if hasattr(native_elem, 'angle'):
            params['angle'] = native_elem.angle
        if hasattr(native_elem, 'pole_gap'):
            params['pole_gap'] = native_elem.pole_gap
        if hasattr(native_elem, 'dipole_length'):
            params['dipole_length'] = native_elem.dipole_length
        if hasattr(native_elem, 'dipole_angle'):
            params['dipole_angle'] = native_elem.dipole_angle
        if hasattr(native_elem, 'fringeType'):
            params['fringe_type'] = native_elem.fringeType
        if hasattr(native_elem, 'name') and native_elem.name:
            params['name'] = native_elem.name
        # RF cavity parameters
        for attr in ('frequency_hz', 'phase_deg', 'voltage_mv',
                     'gradient_mv_per_m', 'structure_type',
                     'phase_advance_deg', 'n_cells'):
            if hasattr(native_elem, attr):
                params[attr] = getattr(native_elem, attr)

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

            # Dispersion and dispersion-corrected emittance
            dispersion = dispersion_prime = 0.0
            disp_corrected = False
            if particles.shape[1] > 5:
                delta = particles[:, 5]
                var_delta = np.var(delta, ddof=1)
                if var_delta > 1e-30:
                    dispersion = np.cov(particles[:, pos_idx], delta, ddof=1)[0, 1] / var_delta
                    dispersion_prime = np.cov(particles[:, ang_idx], delta, ddof=1)[0, 1] / var_delta
                    sig_x2_corr = sig_x2 - dispersion**2 * var_delta
                    sig_xp2_corr = sig_xp2 - dispersion_prime**2 * var_delta
                    sig_xxp_corr = sig_xxp - dispersion * dispersion_prime * var_delta
                    emit_sq = max(0, sig_x2_corr * sig_xp2_corr - sig_xxp_corr**2)
                    disp_corrected = True
                else:
                    emit_sq = sig_x2 * sig_xp2 - sig_xxp**2
            else:
                emit_sq = sig_x2 * sig_xp2 - sig_xxp**2
            emittance = np.sqrt(max(0, emit_sq))  # π·mm·mrad

            if emittance > 0:
                # Use corrected moments when dispersion correction was applied
                s_x2 = sig_x2_corr if disp_corrected else sig_x2
                s_xp2 = sig_xp2_corr if disp_corrected else sig_xp2
                s_xxp = sig_xxp_corr if disp_corrected else sig_xxp
                beta = s_x2 / emittance
                alpha = -s_xxp / emittance
                gamma = s_xp2 / emittance
                if beta > 1e6:
                    self.logger.warning(f"Unphysical beta_{plane} = {beta:.1e} m — beam may be mismatched")
            else:
                beta = alpha = gamma = 0.0
                self.logger.warning(f"Zero emittance in {plane}-plane — degenerate beam distribution")

            twiss[plane] = {
                'beta': beta, 'alpha': alpha, 'gamma': gamma, 'emittance': emittance,
                'dispersion': dispersion, 'dispersion_prime': dispersion_prime,
            }

        return twiss

    def _current_to_k1(self, current: float, length: float, focusing: bool = True) -> float:
        """
        Convert quadrupole current to normalized gradient k1.

        Uses k = |Q·G·I| / (M·C·β·γ), consistent with FELsim's beamline.py.
        """
        if length <= 0 or current == 0:
            return 0.0

        mass_kg = self.particle_mass * PhysicalConstants.MeV_to_J / PhysicalConstants.C**2
        charge_C = abs(self.particle_charge) * PhysicalConstants.Q
        k1 = abs(charge_C * self.G_quad * current) / (
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
                particles[:, idx+1] = sigma_xp * (-alpha * u1 + u2) if beta > 0 else 0

            particles[:, 4] = np.random.randn(num_particles) * std_dev[4]
            particles[:, 5] = np.random.randn(num_particles) * std_dev[5]

        else:
            self.logger.warning(f"Distribution '{distribution_type}' not implemented, using Gaussian")
            particles = np.random.randn(num_particles, 6) * std_dev + mean

        return particles

    def set_beam_energy(self, energy_mev: float):
        """Set beam kinetic energy. A built lattice is rebuilt: dipoles,
        wedges and sextupoles are set for the reference momentum."""
        super().set_beam_energy(energy_mev)
        self.beam_energy = energy_mev
        self._update_relativistic_params()
        self.logger.debug(f"Energy: {energy_mev} MeV, γ={self._gamma:.2f}, Pc={self._Pc:.2f} MeV/c")
        if self._lattice is not None and self.beamline:
            self._build_lattice()

    def set_particle_type(self, mass_mev: float, charge: float):
        """Set particle species."""
        self.particle_charge = self._charged(charge)
        self.particle_mass = mass_mev
        self._update_relativistic_params()
        self.logger.info(f"Particle: m={mass_mev} MeV/c², q={charge}e")
        if self._lattice is not None and self.beamline:
            self._build_lattice()

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
            'bunch_charge_nc': self.bunch_charge_C * 1e9,
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
        """RF-Track bunch of the particles last handed to ``simulate``, as built
        (without the reference particle)."""
        return self._bunch


# Convenience function
def create_rftrack_simulator(**kwargs) -> RFTrackAdapter:
    """Create RF-Track simulator instance."""
    return RFTrackAdapter(**kwargs)
