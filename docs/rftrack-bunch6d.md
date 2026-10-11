# RF-Track adapter: Bunch6d argument order

## The bug

RF-Track's bunch constructor is `Bunch6d(mass, population, charge, X)`.
`backend/rftrackAdapter.py` called it as `Bunch6d(mass, particle_charge, Pc, X)`
since its first version, so every macro-particle carried a charge of +Pc
(about 45.5 e at 45 MeV) and the bunch had no population: RF-Track 2.5.5
stores N = -1/n per macro-particle and `get_population()` returns 0.

Consequences through the adapter:

- Space charge did nothing, since the bunch carried no charge.
- RF cavities gave the wrong energy gain: 58.6 MeV instead of 41.5 MeV out
  of the SLAC linac in `var/slac_linac.json` for a 1 MeV electron.
- Elements whose strength RF-Track scales with P/q looked wrong, which led
  to workarounds: dipoles tracked as drifts with an analytical sector-bend
  correction (the "broken SBend"), and quadrupole strengths that were right
  only because the charge happened to be Pc. With a correctly built bunch,
  RF-Track's `SBend` and quadrupole conventions are exact.
- Transverse transport and particle losses on a line of drifts, quadrupoles
  and dipoles without space charge were correct, because the workarounds
  compensated.

## The fix

The adapter no longer calls the positional constructor. Every bunch is built
from per-particle rows

    [X XP Y YP T P MASS Q N ID]

(`_bunch_rows`, `_make_bunch`), so the question of which positional argument
is the population and which the charge does not arise, and the code works the
same way on RF-Track 2.5.x and 2.7.x.

- `MASS` is the particle mass and `Q` is `particle_charge` (default -1).
- `N` is the number of real particles a macro-particle stands for:
  `N = Q_bunch / e / n_initial`, with the bunch charge given as
  `bunch_charge_nc` (a magnitude in nC; the sign comes from
  `particle_charge`). Space charge without a bunch charge raises.
- `ID` is the row index of the particle in the input array.
- A zero-weight reference particle (on axis, nominal momentum, t = 0)
  leads every bunch. RF cavities are phased on it, and the longitudinal
  output is the arrival time relative to it. With space charge on it is
  tracked in a bunch of its own, outside the beam's mesh.

Quadrupoles and sextupoles are set as (P/q) k1 L, solenoids as
`Solenoid(L, B0, R)`, and dipoles are RF-Track's own `SBend(L, angle, P/q)`
with the wedge edges as thin `SBend`s, all in one `Lattice`.
`dipole_slices` is deprecated and has no effect. With space charge on, each
element gets max(1, ceil(L / sc_ds_max)) space-charge steps
(`sc_ds_max`, default 0.05 m).

### Weights through losses

`n_initial` is the number of particles in the bunch before any losses.
`simulate()` and `track_elements()` take it, and `MultiCodeSimulator` passes
the size of the beam that entered its first section to every section whose
adapter accepts it. A particle that survives a lossy element, a
`track_elements` range or an earlier code section therefore keeps the charge
`Q_bunch / n_initial`, and the bunch charge falls with the losses instead of
being redistributed over the survivors. With space charge on, calling
`track_elements` on a range that does not start at the first element without
`n_initial` warns once that the given particles are taken as the whole bunch.

### Known limit

Between `track_elements` ranges and between multi-code sections the beam is
passed as an (n, 6) array, so the bunch is rebuilt with equal weights
`Q_bunch / e / n_initial`. A beam whose macro-particles carry unequal
weights, such as a gun dump with per-particle charges, is not supported; the
weights would be replaced by equal shares at the first boundary.

## SLAC linac energy test

`backend/test/rftrack_linac/test_adapter_linac.py` loads `var/slac_linac.json`
(one 3.048 m S-band travelling-wave structure) through the adapter and tracks
a single electron injected on axis at 1 MeV kinetic energy, with no space
charge and no aperture. It checks the energy only: the kinetic energy out is
41.47 MeV within 0.05 MeV, the value of the standalone RF-Track script with
the cavity autophased. It also checks that the phase does not depend on the
beam, and that a new beam energy rephases the cavity.

    cd backend
    python -m pytest test/rftrack_linac/test_adapter_linac.py -v

The other adapter tests:

    python -m pytest test/test_rftrack.py test/test_rftrack_bunch.py \
        test/test_rftrack_elements.py test/test_rftrack_quad_convention.py \
        test/test_rftrack_dipoles.py test/test_rftrack_space_charge.py \
        test/test_rftrack_unified.py test/test_multicode.py -v

The COSY-backed tests in `test_multicode.py` run COSY INFINITY. The COSY
simulator looks for the `cosy` binary and `cosy.fox` in `backend/`, the
working directory, `fields/`, `/usr/local/bin` and `~/COSY/10.2/UNIX`, copies
them into `results/` and runs them there, so these tests fail on a machine
without COSY INFINITY. `test_rftrack_unified.py` and the FELsim sections of `test_multicode.py`
build an `ebeam`, which calls `plt.cm.get_cmap`; that was removed in
matplotlib 3.9, so these tests need matplotlib 3.8 or older.

## RF-Track versions

The tests above pass on RF-Track 2.5.5 and 2.7.2. On 2.7.2 the SLAC energy
agrees with 2.5.5 to 1e-14 relative, and transport without space charge on
the UH line is identical except for coord5 at 7e-11. Space charge differs
slightly between the versions: the drift known-answer slope by 5e-6
relative at 1 nC, and up to 2 to 3 percent of the space-charge kick for a
2 nC bunch on a 32^3 mesh.

On 2.7.2 a space-charge result is also not always bit-reproducible between
processes when RF-Track runs on many threads: for 2 nC on a 32^3 mesh two
processes differed by up to 8.8e-4 in x' against a 3e-2 kick. Within one
process it is reproducible. The engine-isolation tests in
`test_rftrack_space_charge.py` run each case in its own process and compare
to 1e-12, so each takes its reference from the first track in the same
process, before any other engine exists, rather than from a separate
process.
