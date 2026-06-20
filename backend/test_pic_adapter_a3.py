#!/usr/bin/env python3
"""A3 test for PicAdapter: intra-drift SC splitting, integration convergence,
and longitudinal space charge.
  1. drift-split is transport-preserving (expanded beamline == native, SC off);
  2. finer sc_ds_max -> more SC stations AND a converging result;
  3. longitudinal SC changes the energy spread (the Ez kick is applied).
Run: PIC_CORE_CLI=<...> COSY_PIC_DIR=<...> python3 test_pic_adapter_a3.py
"""
import sys
import numpy as np

from simulatorFactory import SimulatorFactory
from ebeam import beam as ebeam_class

LATTICE = "../var/UH_FEL_beamline.json"
N = 300


def emit(twiss, axis):
    d = twiss['final'][axis]
    for k in d:
        if 'epsilon' in k.lower() or 'mm.mrad' in k:
            return float(d[k])
    return float('nan')


def make_bunch():
    b = ebeam_class()
    return b.gen_6d_gaussian([0, 0, 0, 0, 0, 0], [1.0, 0.1, 1.0, 0.1, 1.0, 1.0], N)


def main():
    P0 = make_bunch()

    # 1. drift-split transport identity (no SC): expanded beamline tracks
    #    bit-identically to the native one.
    sim = SimulatorFactory.create("pic", lattice_path=LATTICE, beam_energy=45.0,
                                  sc_ds_max=0.05)
    native = sim._native_beamline
    expanded = sim._expanded_beamline(native)
    print(f"native elements={len(native)}  expanded(split drifts)={len(expanded)}")
    a = P0.copy()
    for e in native:
        a = e.useMatrice(a)
    b = P0.copy()
    for e in expanded:
        b = e.useMatrice(b)
    dmax = float(np.abs(a - b).max())
    print(f"drift-split transport identity: max|native-expanded| = {dmax:.2e}")
    if dmax > 1e-9:
        sys.exit("FAIL: drift splitting changed the transport")
    if len(expanded) <= len(native):
        sys.exit("FAIL: no drifts were split (expected expansion)")

    # 2. SC-integration convergence with resolution.
    def run(ds):
        s = SimulatorFactory.create("pic", lattice_path=LATTICE, beam_energy=45.0,
                                    bunch_charge_nc=5.0, sc_mesh=(32, 32, 32),
                                    sc_ds_max=ds)
        s.set_space_charge(True, bunch_charge_nc=5.0)
        r = s.simulate(particles=P0)
        return emit(r.twiss_parameters_statistical, 'x'), r.metadata['num_sc_stations'], r
    m_c, n_c, _ = run(0.20)
    m_m, n_m, _ = run(0.05)
    m_f, n_f, r_f = run(0.0125)
    print(f"stations: coarse={n_c} mid={n_m} fine={n_f}")
    print(f"eps_x:    coarse={m_c:.4e} mid={m_m:.4e} fine={m_f:.4e}")
    if not (n_c < n_m < n_f):
        sys.exit("FAIL: finer sc_ds_max did not add SC stations")
    if r_f.metadata.get('max_particles_out_of_grid', -1) != 0:
        sys.exit("FAIL: bunch left the grid at fine resolution")
    # converging: the step coarse->mid->fine shrinks
    if abs(m_f - m_m) > abs(m_m - m_c):
        sys.exit(f"FAIL: not converging (|fine-mid|={abs(m_f-m_m):.2e} "
                 f">= |mid-coarse|={abs(m_m-m_c):.2e})")

    # 3. longitudinal SC effect: energy spread (coord 5) changes with SC.
    s0 = SimulatorFactory.create("pic", lattice_path=LATTICE, beam_energy=45.0,
                                 bunch_charge_nc=5.0, sc_mesh=(32, 32, 32))
    s0.set_space_charge(False)
    dE_off = float(s0.simulate(particles=P0).final_particles[:, 5].std())
    s0.set_space_charge(True, bunch_charge_nc=5.0)
    dE_on = float(s0.simulate(particles=P0).final_particles[:, 5].std())
    print(f"energy spread (coord5 std): off={dE_off:.4e} on={dE_on:.4e}")
    if abs(dE_on - dE_off) <= 0.0:
        sys.exit("FAIL: longitudinal SC had no effect on energy spread")

    print("test_pic_adapter_a3 PASS")


if __name__ == "__main__":
    main()
