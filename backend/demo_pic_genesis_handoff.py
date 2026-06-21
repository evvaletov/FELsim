#!/usr/bin/env python3
"""A5 demo: cosy-pic -> Genesis4 particle-distribution handoff for the UH FEL.

cosy-pic owns space charge in the EXTERNAL transport (gun -> undulator entrance);
Genesis4 owns space charge INSIDE the undulator. The handoff plane is the
undulator entrance. This demo:

  1. tracks the UH FEL line with the cosy-pic SC engine (PicAdapter) UP TO the
     undulator entrance (the transport region only -- SC applied once, here);
  2. emits the SC-evolved distribution at that plane as ASTRA + openPMD files
     (the formats Genesis4 reads via &importdistribution / &importbeam);
  3. asserts the NO-DOUBLE-COUNT invariant:
        (a) transport SC actually moved the beam (handoff emittance with SC on
            != SC off) -- so it is applied, exactly once, before the plane;
        (b) the written file round-trips to the SAME distribution Genesis4 will
            read -- no SC re-applied and nothing lost in the handoff;
        (c) cosy-pic applies NO SC at/after the undulator entrance (the tracked
            sub-line ends at the plane) -- Genesis4 takes over there.

Run: PIC_CORE_CLI=<...> COSY_PIC_DIR=<...> \
     MPLBACKEND=Agg PYTHONPATH=$(pwd) python3 demo_pic_genesis_handoff.py
"""
import os
import sys

import numpy as np

from simulatorFactory import SimulatorFactory
from ebeam import beam as ebeam_class

sys.path.insert(0, os.path.join(
    os.environ.get('COSY_PIC_DIR', os.path.expanduser('~/COSY/cosy-pic')), 'tools'))
from genesis_handoff import (write_astra, read_astra, write_openpmd,  # noqa: E402
                             transverse_emittance)

LATTICE = "../var/UH_FEL_beamline.json"
OUTDIR = "test/results/pic_genesis_handoff"
ENERGY = 45.0
CHARGE_NC = 0.1
N = 2000
SEED = 42


def undulator_entrance_index(native):
    for i, e in enumerate(native):
        nm = getattr(e, 'name', '') or ''
        if 'UND' in type(e).__name__.upper() or 'UND' in str(nm).upper():
            return i
    raise RuntimeError("no undulator found in lattice")


def beam_emit(beam):
    return (transverse_emittance(beam.x, beam.px / beam.pz),
            transverse_emittance(beam.y, beam.py / beam.pz))


def si_moments(x, y, z, px, py, pz, q, w):
    """Order-independent SI moments -- robust to SC distortion (no pz division),
    so they prove the written distribution is read back verbatim."""
    cov = lambda a, b: float(np.cov(a, b)[0, 1])  # noqa: E731
    epsx = np.sqrt(max(0.0, cov(x, x) * cov(px, px) - cov(x, px)**2))
    epsy = np.sqrt(max(0.0, cov(y, y) * cov(py, py) - cov(y, py)**2))
    return dict(n=len(x), qtot=float((q * w).sum()),
                cx=float(x.mean()), cz=float(z.mean()),
                sx=float(x.std()), sz=float(z.std()),
                spx=float(px.std()), epsx_si=float(epsx), epsy_si=float(epsy))


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    b = ebeam_class()
    np.random.seed(SEED)
    P0 = b.gen_6d_gaussian([0, 0, 0, 0, 0, 0], [1.0, 0.1, 1.0, 0.1, 1.0, 1.0], N)

    sim = SimulatorFactory.create("pic", lattice_path=LATTICE, beam_energy=ENERGY,
                                  bunch_charge_nc=CHARGE_NC, sc_mesh=(32, 32, 32),
                                  sc_ds_max=0.1)
    native = sim.get_native_beamline()
    iu = undulator_entrance_index(native)
    transport = native[:iu]               # gun -> undulator entrance (exclusive)
    print(f"lattice: {len(native)} elements; undulator entrance at index {iu} "
          f"({getattr(native[iu], 'name', '')}); transport region = {len(transport)} elements")

    # ---- track the transport with SC on, then SC off (same sub-line) --------
    sim.set_beamline(list(transport))
    sim.set_space_charge(True, bunch_charge_nc=CHARGE_NC)
    r_on = sim.simulate(particles=P0.copy())
    P_on = r_on.final_particles

    sim.set_beamline(list(transport))
    sim.set_space_charge(False)
    r_off = sim.simulate(particles=P0.copy())
    P_off = r_off.final_particles

    beam_on = sim._felsim_to_beam(P_on)
    beam_off = sim._felsim_to_beam(P_off)
    ex_on, ey_on = beam_emit(beam_on)
    ex_off, ey_off = beam_emit(beam_off)
    print(f"handoff-plane eps (geom):  SC off  eps_x={ex_off:.4e} eps_y={ey_off:.4e}")
    print(f"                           SC on   eps_x={ex_on:.4e} eps_y={ey_on:.4e}")

    # ---- emit the SC-evolved distribution (Genesis4 input) ------------------
    astra = os.path.join(OUTDIR, "undulator_entrance.astra")
    opmd = os.path.join(OUTDIR, "undulator_entrance.h5")
    write_astra(astra, beam_on.x, beam_on.y, beam_on.z,
                beam_on.px, beam_on.py, beam_on.pz, beam_on.q, beam_on.w)
    have_opmd = write_openpmd(opmd, beam_on.x, beam_on.y, beam_on.z,
                              beam_on.px, beam_on.py, beam_on.pz, beam_on.q, beam_on.w)
    print(f"wrote {astra}" + (f" and {opmd}" if have_opmd else " (openPMD skipped: no h5py)"))

    # ---- NO-DOUBLE-COUNT assertions -----------------------------------------
    rc = 0
    # (a) transport SC actually moved the beam (applied once, before the plane)
    if abs(ex_on / ex_off - 1.0) < 1e-6 and abs(ey_on / ey_off - 1.0) < 1e-6:
        print("FAIL (a): transport SC had no effect at the handoff plane")
        rc = 1
    else:
        print(f"PASS (a): transport SC moved eps_x by {100*(ex_on/ex_off-1):+.2f}% "
              f"(applied exactly once, before the undulator)")

    # (b) the written file is EXACTLY what Genesis4 will read (round-trip).
    # Compare SI moments (robust to SC distortion, no pz division).
    m_w = si_moments(beam_on.x, beam_on.y, beam_on.z,
                     beam_on.px, beam_on.py, beam_on.pz, beam_on.q, beam_on.w)
    m_r = si_moments(*read_astra(astra))
    worst = 0.0
    for k in ('n', 'qtot', 'cx', 'cz', 'sx', 'sz', 'spx', 'epsx_si', 'epsy_si'):
        ref = abs(m_w[k]) if abs(m_w[k]) > 0 else 1.0
        worst = max(worst, abs(m_r[k] - m_w[k]) / ref)
    if worst < 1e-9:
        print(f"PASS (b): handoff round-trips as a set (n={m_r['n']}, all SI moments "
              f"-- charge, sizes, emittances -- preserved to {worst:.1e}); Genesis4 "
              f"receives the transport-SC-evolved distribution")
    else:
        print(f"FAIL (b): handoff not faithful (worst SI-moment rel error {worst:.1e})")
        rc = 1

    # (c) cosy-pic applied NO SC at/after the undulator entrance: the tracked
    # sub-line is native[:iu], so the undulator is excluded by construction, and
    # SC stations ran along the transport region (disjoint from Genesis4's region).
    und_in_transport = any('UND' in (type(e).__name__ + str(getattr(e, 'name', ''))).upper()
                           for e in transport)
    stations = r_on.metadata.get('num_sc_stations', 0)
    if (not und_in_transport) and stations > 0:
        print(f"PASS (c): cosy-pic SC ran on the transport region only "
              f"({stations} SC stations, undulator excluded); Genesis4 owns SC inside "
              f"the undulator (elements {iu}..{len(native)-1}) -- disjoint, no double count")
    else:
        print(f"FAIL (c): SC region not disjoint from the undulator "
              f"(und_in_transport={und_in_transport}, stations={stations})")
        rc = 1

    print("\ndemo_pic_genesis_handoff:", "PASS" if rc == 0 else "FAIL")
    return rc


if __name__ == "__main__":
    sys.exit(main())
