"""
Standalone RF-Track reference for test_adapter_linac.py.

Builds the SLAC 3-m S-band travelling-wave structure of var/slac_linac.json
directly in RF-Track, without the adapter, and tracks one electron injected
on axis at 1 MeV kinetic energy:

- TW_Structure with a single Fourier coefficient, the peak on-axis field
  E0 = gradient_mv_per_m (13.3 MV/m), harmonic 0, f = 2856 MHz, phase
  advance phase_advance_deg (2 pi / 3) per cell;
- cell length c * phase_advance / (2 pi f) (synchronous with beta = 1), and
  a fractional number of cells that gives the element length (3.048 m);
- phid = phase_deg (0). RF-Track already phases the structure on the first
  particle, so phid = 0 is the crest; the autophase call is kept as a guard
  and does not change the result.

It prints the kinetic energy out, 41.47 MeV.

    cd backend
    python test/rftrack_linac/linac_standalone.py
"""

import json
from pathlib import Path

import numpy as np
import RF_Track as rft

MC2_MEV = 0.51099895
CLIGHT = 299792458.0
LATTICE_JSON = Path(__file__).resolve().parents[3] / 'var' / 'slac_linac.json'


def build_structure(elem):
    p = elem['parameters']
    freq = p['frequency_hz']
    phase_adv = np.deg2rad(p['phase_advance_deg'])
    l_cell = CLIGHT * phase_adv / (2.0 * np.pi * freq)
    tw = rft.TW_Structure(p['gradient_mv_per_m'] * 1e6, 0, freq, phase_adv,
                          elem['length_m'] / l_cell)
    tw.set_phid(p['phase_deg'])
    return tw


def electron(k_mev):
    """One electron on axis, as a per-particle row [X XP Y YP T P MASS Q N ID]."""
    pc = np.sqrt((k_mev + MC2_MEV)**2 - MC2_MEV**2)
    return rft.Bunch6d(np.array([[0.0, 0.0, 0.0, 0.0, 0.0, pc, MC2_MEV, -1.0, 1.0, 0.0]]))


def main():
    line = json.loads(LATTICE_JSON.read_text())['beamline']
    k_in = line['beam_parameters']['particle']['kinetic_energy_mev']
    tw = build_structure(line['elements'][0])
    lat = rft.Lattice()
    lat.append(tw)
    lat.autophase(electron(k_in))
    out = lat.track(electron(k_in)).get_phase_space('%Pc')
    if np.size(out) == 0:
        raise SystemExit('electron lost')
    k_out = np.hypot(float(np.ravel(out)[0]), MC2_MEV) - MC2_MEV
    print(f'structure length {tw.get_length():.4f} m')
    print(f'K_in {k_in:.3f} MeV, K_out {k_out:.4f} MeV')


if __name__ == '__main__':
    main()
