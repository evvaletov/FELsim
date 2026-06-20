#!/usr/bin/env python3
"""Smoke test for the cosy-pic PIC space-charge adapter in the FELsim
multi-code framework. Tracks a Gaussian bunch through the UH FEL beamline with
space charge OFF and ON via the factory-created 'pic' adapter and checks:
  - the adapter is selectable via SimulatorFactory ('pic');
  - SC-off and SC-on runs both succeed and keep the bunch inside the grid;
  - SC has a measurable effect (final distribution / emittance changes).
Run: PIC_CORE_CLI=<...> COSY_PIC_DIR=<...> python3 test_pic_adapter.py
"""
import sys
import numpy as np

from simulatorFactory import SimulatorFactory
from ebeam import beam as ebeam_class

LATTICE = "../var/UH_FEL_beamline.json"
N = 600


def emit(twiss, axis):
    # geometric emittance from the statistical twiss dict (key is LaTeX-named)
    d = twiss['final'][axis]
    for k in d:
        if 'epsilon' in k.lower() or 'mm.mrad' in k:
            return float(d[k])
    return float('nan')


def main():
    b = ebeam_class()
    mean = [0, 0, 0, 0, 0, 0]
    std = [1.0, 0.1, 1.0, 0.1, 1.0, 1.0]   # mm, mrad, mm, mrad, dToF*1e3, dK*1e3
    particles = b.gen_6d_gaussian(mean, std, N)

    # higher charge so the 45 MeV SC effect is clearly measurable in a smoke test
    sim = SimulatorFactory.create("pic", lattice_path=LATTICE,
                                  beam_energy=45.0, bunch_charge_nc=5.0,
                                  sc_mesh=(32, 32, 32))

    sim.set_space_charge(False)
    r_off = sim.simulate(particles=particles)
    sim.set_space_charge(True, bunch_charge_nc=5.0)
    r_on = sim.simulate(particles=particles)

    if not (r_off.success and r_on.success):
        sys.exit("FAIL: a run did not succeed")
    oob = r_on.metadata.get('max_particles_out_of_grid', -1)
    print(f"sc_on oob={oob}  elements={r_on.metadata['num_elements']}")
    if oob != 0:
        sys.exit(f"FAIL: {oob} particles left the rest-frame grid")

    ex_off, ey_off = emit(r_off.twiss_parameters_statistical, 'x'), emit(r_off.twiss_parameters_statistical, 'y')
    ex_on, ey_on = emit(r_on.twiss_parameters_statistical, 'x'), emit(r_on.twiss_parameters_statistical, 'y')
    print(f"eps_x [pi.mm.mrad]: off={ex_off:.4e} on={ex_on:.4e} ({100*(ex_on/ex_off-1):+.2f}%)")
    print(f"eps_y [pi.mm.mrad]: off={ey_off:.4e} on={ey_on:.4e} ({100*(ey_on/ey_off-1):+.2f}%)")

    # SC must have a measurable effect on the final distribution
    dpos = np.abs(r_on.final_particles - r_off.final_particles).max()
    print(f"max |final_on - final_off| = {dpos:.4e}")
    if dpos <= 0.0:
        sys.exit("FAIL: space charge had no effect (SC-on == SC-off)")

    print("test_pic_adapter PASS")


if __name__ == "__main__":
    main()
