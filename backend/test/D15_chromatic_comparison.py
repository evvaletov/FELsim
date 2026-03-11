#!/usr/bin/env python3
"""D15: Three-way comparison — linear FELsim vs chromatic FELsim vs RF-Track.

Tests whether chromatic quadrupole matrices close the gap between FELsim
and RF-Track for beams with nonzero energy spread.

Author: Eremey Valetov
"""

import sys
from pathlib import Path
import numpy as np

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

from ebeam import beam
from beamline import lattice, qpfLattice, qpdLattice
from excelElements import ExcelElements
from simulatorBase import CoordinateSystem
from rftrackAdapter import RFTrackAdapter

import json

EXCEL_PATH = backend_dir.parent / 'beam_excel' / 'Beamline_elements.xlsx'
CURRENTS_PATH = backend_dir / 'test' / 'results' / 'felsim_nm_warm.json'
ENERGY = 40
EPSILON_N = 8
X_STD = 0.8
Y_STD = 0.8
FREQ = 2856e6
BUNCH_SPREAD = 2
NB_PARTICLES = 500
SEGMENTS = 118
SEED = 42


def load_currents():
    if CURRENTS_PATH.exists():
        with open(CURRENTS_PATH) as f:
            data = json.load(f)
        return {int(k): v for k, v in data['currents'].items()}
    return {}


def build_felsim_line(chromatic=False, currents=None):
    relat = lattice(1, fringeType=None)
    relat.setE(E=ENERGY)
    excel = ExcelElements(str(EXCEL_PATH))
    beamlineUH = excel.create_beamline()
    line = relat.changeBeamType("electron", ENERGY, beamlineUH)[:SEGMENTS]
    if currents:
        for idx, current in currents.items():
            if idx < len(line):
                line[idx].current = abs(current)
    if chromatic:
        for elem in line:
            elem.chromatic = True
    return line, relat


def build_rftrack_adapter(currents=None):
    sim = RFTrackAdapter(
        lattice_path=str(EXCEL_PATH),
        beam_energy=ENERGY,
        space_charge=False,
        aperture=0.5,
    )
    sim.beamline = sim.beamline[:SEGMENTS]
    if currents:
        for idx, current in currents.items():
            if idx < len(sim.beamline):
                sim.beamline[idx].parameters['current'] = abs(current)
    sim._build_lattice()
    return sim


def generate_beam(relat, energy_spread_pct):
    np.random.seed(SEED)
    norm = relat.gamma * relat.beta
    epsilon = EPSILON_N / norm
    sigma_p = energy_spread_pct * 10  # coord6 = ΔK/K₀ × 10³
    beam_dist = beam().gen_6d_gaussian(
        0, [X_STD, epsilon / X_STD, Y_STD, epsilon / Y_STD,
            BUNCH_SPREAD * 1e-9 * FREQ, sigma_p],
        NB_PARTICLES)
    return beam_dist


def felsim_track(line, beam_dist):
    ps = beam_dist.copy()
    for elem in line:
        ps = np.array(elem.useMatrice(ps))
    return ps


def rftrack_track(sim, beam_dist):
    ps_rft_in = sim.transform_coordinates(
        beam_dist.copy(), CoordinateSystem.FELSIM, CoordinateSystem.RFTRACK)
    ps_rft_out = sim.track_elements(ps_rft_in, 0, len(sim.beamline))
    return sim.transform_coordinates(
        ps_rft_out, CoordinateSystem.RFTRACK, CoordinateSystem.FELSIM)


def compute_twiss(ps):
    eb = beam()
    _, _, tw = eb.cal_twiss(ps, ddof=1)
    return tw


def main():
    print("=" * 90)
    print("  D15: Three-way comparison — linear FELsim vs chromatic FELsim vs RF-Track")
    print("=" * 90)
    print()

    currents = load_currents()
    print(f"  Currents: {len(currents)} quadrupoles loaded from {CURRENTS_PATH.name}")
    print()

    line_lin, relat = build_felsim_line(chromatic=False, currents=currents)
    line_chr, _ = build_felsim_line(chromatic=True, currents=currents)
    sim = build_rftrack_adapter(currents=currents)

    spreads = [0.0, 0.1, 0.25, 0.5, 1.0]

    # Part 1: Beta functions
    print("  Beta functions at undulator entrance (element 117):")
    print()
    print(f"  {'σ_p(%)':>8}"
          f" {'βx_lin':>9} {'βx_chr':>9} {'βx_rft':>9} {'Δβx_lin%':>9} {'Δβx_chr%':>9}"
          f" {'βy_lin':>9} {'βy_chr':>9} {'βy_rft':>9} {'Δβy_lin%':>9} {'Δβy_chr%':>9}")
    print(f"  {'-' * 98}")

    for sp in spreads:
        beam_dist = generate_beam(relat, energy_spread_pct=sp)
        if sp == 0.0:
            beam_dist[:, 5] = 0.0
            beam_dist[:, 4] = 0.0

        ps_lin = felsim_track(line_lin, beam_dist)
        ps_chr = felsim_track(line_chr, beam_dist)
        ps_rft = rftrack_track(sim, beam_dist)

        tw_lin = compute_twiss(ps_lin)
        tw_chr = compute_twiss(ps_chr)
        tw_rft = compute_twiss(ps_rft)

        bx_l = tw_lin.loc['x'][r"$\beta$ (m)"]
        bx_c = tw_chr.loc['x'][r"$\beta$ (m)"]
        bx_r = tw_rft.loc['x'][r"$\beta$ (m)"]
        by_l = tw_lin.loc['y'][r"$\beta$ (m)"]
        by_c = tw_chr.loc['y'][r"$\beta$ (m)"]
        by_r = tw_rft.loc['y'][r"$\beta$ (m)"]

        # Relative error vs RF-Track (the "truth" for this comparison)
        dbx_l = (bx_l - bx_r) / bx_r * 100 if bx_r > 1e-15 else 0
        dbx_c = (bx_c - bx_r) / bx_r * 100 if bx_r > 1e-15 else 0
        dby_l = (by_l - by_r) / by_r * 100 if by_r > 1e-15 else 0
        dby_c = (by_c - by_r) / by_r * 100 if by_r > 1e-15 else 0

        print(f"  {sp:8.2f}"
              f" {bx_l:9.4f} {bx_c:9.4f} {bx_r:9.4f} {dbx_l:9.2f} {dbx_c:9.2f}"
              f" {by_l:9.4f} {by_c:9.4f} {by_r:9.4f} {dby_l:9.2f} {dby_c:9.2f}")

    # Part 2: Alpha functions
    print()
    print(f"  {'σ_p(%)':>8}"
          f" {'αx_lin':>9} {'αx_chr':>9} {'αx_rft':>9} {'Δαx_lin':>9} {'Δαx_chr':>9}"
          f" {'αy_lin':>9} {'αy_chr':>9} {'αy_rft':>9} {'Δαy_lin':>9} {'Δαy_chr':>9}")
    print(f"  {'-' * 98}")

    for sp in spreads:
        beam_dist = generate_beam(relat, energy_spread_pct=sp)
        if sp == 0.0:
            beam_dist[:, 5] = 0.0
            beam_dist[:, 4] = 0.0

        ps_lin = felsim_track(line_lin, beam_dist)
        ps_chr = felsim_track(line_chr, beam_dist)
        ps_rft = rftrack_track(sim, beam_dist)

        tw_lin = compute_twiss(ps_lin)
        tw_chr = compute_twiss(ps_chr)
        tw_rft = compute_twiss(ps_rft)

        ax_l = tw_lin.loc['x'][r"$\alpha$"]
        ax_c = tw_chr.loc['x'][r"$\alpha$"]
        ax_r = tw_rft.loc['x'][r"$\alpha$"]
        ay_l = tw_lin.loc['y'][r"$\alpha$"]
        ay_c = tw_chr.loc['y'][r"$\alpha$"]
        ay_r = tw_rft.loc['y'][r"$\alpha$"]

        print(f"  {sp:8.2f}"
              f" {ax_l:9.4f} {ax_c:9.4f} {ax_r:9.4f} {ax_l - ax_r:9.4f} {ax_c - ax_r:9.4f}"
              f" {ay_l:9.4f} {ay_c:9.4f} {ay_r:9.4f} {ay_l - ay_r:9.4f} {ay_c - ay_r:9.4f}")

    # Part 3: Emittance
    print()
    print(f"  {'σ_p(%)':>8}"
          f" {'εx_lin':>9} {'εx_chr':>9} {'εx_rft':>9} {'Δεx_l%':>9} {'Δεx_c%':>9}"
          f" {'εy_lin':>9} {'εy_chr':>9} {'εy_rft':>9} {'Δεy_l%':>9} {'Δεy_c%':>9}")
    print(f"  {'-' * 98}")

    for sp in spreads:
        beam_dist = generate_beam(relat, energy_spread_pct=sp)
        if sp == 0.0:
            beam_dist[:, 5] = 0.0
            beam_dist[:, 4] = 0.0

        ps_lin = felsim_track(line_lin, beam_dist)
        ps_chr = felsim_track(line_chr, beam_dist)
        ps_rft = rftrack_track(sim, beam_dist)

        tw_lin = compute_twiss(ps_lin)
        tw_chr = compute_twiss(ps_chr)
        tw_rft = compute_twiss(ps_rft)

        ex_l = tw_lin.loc['x'][r"$\epsilon$ ($\pi$.mm.mrad)"]
        ex_c = tw_chr.loc['x'][r"$\epsilon$ ($\pi$.mm.mrad)"]
        ex_r = tw_rft.loc['x'][r"$\epsilon$ ($\pi$.mm.mrad)"]
        ey_l = tw_lin.loc['y'][r"$\epsilon$ ($\pi$.mm.mrad)"]
        ey_c = tw_chr.loc['y'][r"$\epsilon$ ($\pi$.mm.mrad)"]
        ey_r = tw_rft.loc['y'][r"$\epsilon$ ($\pi$.mm.mrad)"]

        dex_l = (ex_l - ex_r) / ex_r * 100 if ex_r > 1e-15 else 0
        dex_c = (ex_c - ex_r) / ex_r * 100 if ex_r > 1e-15 else 0
        dey_l = (ey_l - ey_r) / ey_r * 100 if ey_r > 1e-15 else 0
        dey_c = (ey_c - ey_r) / ey_r * 100 if ey_r > 1e-15 else 0

        print(f"  {sp:8.2f}"
              f" {ex_l:9.4f} {ex_c:9.4f} {ex_r:9.4f} {dex_l:9.2f} {dex_c:9.2f}"
              f" {ey_l:9.4f} {ey_c:9.4f} {ey_r:9.4f} {dey_l:9.2f} {dey_c:9.2f}")

    # Part 4: Beam sizes
    print()
    print(f"  {'σ_p(%)':>8}"
          f" {'σx_lin':>9} {'σx_chr':>9} {'σx_rft':>9} {'Δσx_l%':>9} {'Δσx_c%':>9}"
          f" {'σy_lin':>9} {'σy_chr':>9} {'σy_rft':>9} {'Δσy_l%':>9} {'Δσy_c%':>9}")
    print(f"  {'-' * 98}")

    for sp in spreads:
        beam_dist = generate_beam(relat, energy_spread_pct=sp)
        if sp == 0.0:
            beam_dist[:, 5] = 0.0
            beam_dist[:, 4] = 0.0

        ps_lin = felsim_track(line_lin, beam_dist)
        ps_chr = felsim_track(line_chr, beam_dist)
        ps_rft = rftrack_track(sim, beam_dist)

        sx_l = np.std(ps_lin[:, 0], ddof=1)
        sx_c = np.std(ps_chr[:, 0], ddof=1)
        sx_r = np.std(ps_rft[:, 0], ddof=1)
        sy_l = np.std(ps_lin[:, 2], ddof=1)
        sy_c = np.std(ps_chr[:, 2], ddof=1)
        sy_r = np.std(ps_rft[:, 2], ddof=1)

        dsx_l = (sx_l - sx_r) / sx_r * 100 if sx_r > 1e-15 else 0
        dsx_c = (sx_c - sx_r) / sx_r * 100 if sx_r > 1e-15 else 0
        dsy_l = (sy_l - sy_r) / sy_r * 100 if sy_r > 1e-15 else 0
        dsy_c = (sy_c - sy_r) / sy_r * 100 if sy_r > 1e-15 else 0

        print(f"  {sp:8.2f}"
              f" {sx_l:9.4f} {sx_c:9.4f} {sx_r:9.4f} {dsx_l:9.2f} {dsx_c:9.2f}"
              f" {sy_l:9.4f} {sy_c:9.4f} {sy_r:9.4f} {dsy_l:9.2f} {dsy_c:9.2f}")

    print()


if __name__ == "__main__":
    main()
    print("Done.")
