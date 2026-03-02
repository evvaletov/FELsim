#!/usr/bin/env python3
"""
Minimal Working Example: RF-Track v2.5.5 SBend tracking bug.

Author: Eremey Valetov
Date: 2026-03-02

This script demonstrates that RF-Track's SBend element produces incorrect
tracking results:
  - Without P/Q reference: acts as identity (no field effect)
  - With P/Q reference: produces catastrophically wrong output (~1000mm
    displacements) because the transfer matrix treats Bunch6d's 6th
    coordinate (absolute P in MeV/c) as momentum deviation δ = ΔP/P₀

Expected behavior: a sector bend with angle θ and bending radius ρ should
apply the matrix:
    x_out  = cos(θ) x + ρ sin(θ) x'
    x'_out = -sin(θ)/ρ x + cos(θ) x'

Environment: RF-Track 2.5.5, Python 3.12, Linux x86_64 (RHEL 9)
"""

import RF_Track as rft
import numpy as np
import sys

print(f"Python: {sys.version.split()[0]}")
print()

# Parameters
me = 0.511       # MeV/c^2 (electron mass)
E0 = 40.0        # MeV kinetic energy
P0 = np.sqrt((E0 + me)**2 - me**2)  # MeV/c
P_Q = P0                              # MV/c (charge magnitude = 1)

L = 0.20         # m (dipole length)
angle_deg = 22.5
angle_rad = np.radians(angle_deg)
K0 = angle_rad / L
rho = 1.0 / K0
Brho = P_Q / 299.792458  # T*m

print(f"Electron: E0={E0} MeV, P0={P0:.4f} MeV/c")
print(f"Dipole:   L={L} m, angle={angle_deg}°, K0={K0:.4f} 1/m, rho={rho:.4f} m")
print()

# Tracking helper
ps_offaxis = np.array([[1.0, 0.0, 1.0, 0.0, 0.0, P0]])

def track(elem, particles):
    lat = rft.Lattice()
    lat.append(elem)
    b = rft.Bunch6d(me, -1.0, P0, particles.copy())
    tracked = lat.track(b)
    r = np.array(tracked.get_phase_space())
    if r.size == 0:
        return None  # all particles lost
    return r[0]

def report(label, elem):
    r = track(elem, ps_offaxis)
    if r is None:
        print(f"  {label:38s}  ALL PARTICLES LOST")
    else:
        print(f"  {label:38s}  x={r[0]:12.4f} mm  x'={r[1]:12.4f} mrad")

# Expected
cos_th = np.cos(angle_rad)
sin_th = np.sin(angle_rad)
print(f"Expected sector-bend output for x=1mm off-axis particle:")
print(f"  x  = cos(θ)×1mm = {cos_th:.6f} mm")
print(f"  x' = -sin(θ)/ρ×1mm = {-sin_th/rho:.6f} mrad")
print()

# Reference: Drift and Quadrupole (both track correctly)
d = rft.Drift(L)
r = track(d, ps_offaxis)
print(f"  {'Drift(L) — reference':38s}  x={r[0]:12.4f} mm  x'={r[1]:12.4f} mrad")

q = rft.Quadrupole()
q.set_length(L)
q.set_strength(K0**2 * L)
q.set_aperture(100, 100)
r = track(q, ps_offaxis)
print(f"  {'Quad(K1=K0²) — reference':38s}  x={r[0]:12.4f} mm  x'={r[1]:12.4f} mrad")
print()

# --- Test 1: SBend without P/Q ---
print("=" * 72)
print("Test 1: SBend without P/Q → no field effect")
print("=" * 72)

sb = rft.SBend()
sb.set_length(L)
sb.set_K0(K0)
sb.set_aperture(100, 100)
report("SBend set_K0 only", sb)
print("  → Identity: set_K0 alone sets neither field nor curved reference.")
print()

# --- Test 2: SBend with P/Q → catastrophically wrong ---
print("=" * 72)
print("Test 2: SBend with P/Q → catastrophically wrong output")
print("=" * 72)

sb = rft.SBend()
sb.set_length(L)
sb.set_K0(K0)
sb.set_P_over_Q(P_Q)
sb.set_aperture(100, 100)
report("SBend K0 + P/Q", sb)

sb = rft.SBend()
sb.set_length(L)
sb.set_h(K0)
sb.set_K0(K0)
sb.set_P_over_Q(P_Q)
sb.set_aperture(100, 100)
report("SBend h + K0 + P/Q", sb)

sb = rft.SBend(L, angle_rad, P_Q)
sb.set_aperture(100, 100)
report("SBend(L, angle, P_Q)", sb)

rb = rft.RBend(L, angle_rad, P_Q)
rb.set_aperture(100, 100)
report("RBend(L, angle, P_Q)", rb)

print("  → ~1000mm displacements for a 1mm off-axis particle.")
print()

# --- Test 3: On-axis reference particle ---
print("=" * 72)
print("Test 3: On-axis reference particle through SBend constructor")
print("=" * 72)
ps_onaxis = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P0]])
sb = rft.SBend(L, angle_rad, P_Q)
sb.set_aperture(100, 100)
r = track(sb, ps_onaxis)
if r is not None:
    print(f"  x = {r[0]:.4f} mm  (expected 0)")
    print(f"  The on-axis reference particle should traverse the sector bend")
    print(f"  along the curved reference orbit and exit at x=0.")
    print(f"  The large displacement suggests the P={P0:.2f} MeV/c in the 6th")
    print(f"  Bunch6d coordinate is interpreted as δ = ΔP/P₀ ≈ {P0:.1f},")
    print(f"  i.e., a {P0*100:.0f}% momentum deviation. Dispersion terms then")
    print(f"  produce the huge position offset.")
print()

# --- Test 4: set_E1/set_E2 has no effect ---
print("=" * 72)
print("Test 4: set_E1/set_E2 has no effect (setter-only SBend)")
print("=" * 72)

sb_no_e = rft.SBend()
sb_no_e.set_length(L)
sb_no_e.set_K0(K0)
sb_no_e.set_aperture(100, 100)

sb_with_e = rft.SBend()
sb_with_e.set_length(L)
sb_with_e.set_K0(K0)
sb_with_e.set_E1(angle_rad)
sb_with_e.set_E2(angle_rad)
sb_with_e.set_hgap(0.01)
sb_with_e.set_fint(0.5)
sb_with_e.set_aperture(100, 100)

r_no = track(sb_no_e, ps_offaxis)
r_with = track(sb_with_e, ps_offaxis)
print(f"  Without E1/E2:  x={r_no[0]:10.6f}  x'={r_no[1]:10.6f}")
print(f"  With E1/E2:     x={r_with[0]:10.6f}  x'={r_with[1]:10.6f}")
print(f"  Difference:     Δx={abs(r_with[0]-r_no[0]):.1e}")
print("  → E1, E2, hgap, fint have zero effect on tracking.")
print()

# --- Conclusion ---
print("=" * 72)
print("SUMMARY")
print("=" * 72)
print("RF-Track v2.5.5 SBend:")
print("  1. Without P/Q: acts as identity (no field, no bending)")
print("  2. With P/Q: applies field but interprets absolute P as δ,")
print("     producing ~1000mm displacements for typical electron beams")
print("  3. set_E1/set_E2 have no effect on tracking")
print("  4. RBend has the same issues")
print()
print("For comparison, Quadrupole and Drift track correctly.")
print()
print("The root cause appears to be that the SBend transfer matrix reads")
print("the 6th Bunch6d coordinate (absolute P in MeV/c) as momentum")
print("deviation δ = ΔP/P₀. For a 40 MeV electron (P₀ ≈ 40.5 MeV/c),")
print("the SBend treats this as δ ≈ 40.5 (a 4050% momentum deviation),")
print("which multiplied by the dispersion produces the ~900mm offset.")
print()
print("Suggested fix: in SBend tracking, convert the 6th coordinate from")
print("absolute P to δ = (P - P_ref) / P_ref before applying the transfer")
print("matrix, where P_ref is the beam's reference momentum from Bunch6d.")
