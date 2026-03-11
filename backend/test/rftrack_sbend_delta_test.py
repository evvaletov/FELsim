#!/usr/bin/env python3
"""Phase 1.1: Systematic investigation of RF-Track SBend behaviour.

Key insight from manual (§4.2.3): P_Q is magnetic rigidity P/q [MV/c].
For electrons with q = -1e, P_Q = P₀/(-1) = -P₀.
We have been passing P_Q = +P₀ — wrong sign!

Tests:
  A. P_Q = +P₀ (our current usage — known broken)
  B. P_Q = -P₀ (correct for electrons)
  C. E1/E2 via constructor (not setters)
  D. Full DPW-DPH-DPW triplet
  E. nsteps variation

Author: Eremey Valetov
"""

import numpy as np
import sys

import RF_Track as rft

me = 0.511          # MeV/c^2
E0 = 40.0           # MeV kinetic energy
P0 = np.sqrt((E0 + me)**2 - me**2)

L = 0.20            # m
angle_deg = 22.5
angle_rad = np.radians(angle_deg)
K0 = angle_rad / L
rho = 1.0 / K0

cos_th = np.cos(angle_rad)
sin_th = np.sin(angle_rad)

print(f"E0 = {E0} MeV, P0 = {P0:.4f} MeV/c")
print(f"L = {L} m, angle = {angle_deg}°, ρ = {rho:.6f} m")
print(f"P/q for electron (q=-1): {-P0:.4f} MV/c")
print()

# Sector bend analytical expectations (on-axis, on-momentum → x=0)
# Off-axis x=1mm: x_out = cos(θ), x'_out = -sin(θ)/ρ
# Vertical: drift (y_out = y_in, y'_out = y'_in)

def track(elem, particles, P_ref=P0, charge=-1.0):
    lat = rft.Lattice()
    lat.append(elem)
    b = rft.Bunch6d(me, charge, P_ref, particles.copy())
    tracked = lat.track(b)
    ps = np.array(tracked.get_phase_space())
    if ps.size == 0:
        return None
    return ps[0]


def report(label, r, x_exp=None, xp_exp=None):
    if r is None:
        print(f"  {label:40s}  LOST")
        return
    line = f"  {label:40s}  x={r[0]:10.4f}  x'={r[1]:10.4f}"
    if x_exp is not None:
        line += f"  (exp x={x_exp:.4f})"
    if xp_exp is not None:
        line += f"  (exp x'={xp_exp:.4f})"
    print(line)


# Test particles (Bunch6d format: x[mm], x'[mrad], y[mm], y'[mrad], t[mm/c], P[MeV/c])
ps_onaxis = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P0]])
ps_offx = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, P0]])
ps_offy = np.array([[0.0, 0.0, 1.0, 0.0, 0.0, P0]])
ps_offxp = np.array([[0.0, 1.0, 0.0, 0.0, 0.0, P0]])

# ── Test A: P_Q = +P₀ (current usage — should be broken) ─────────────────
print("=" * 78)
print("Test A: P_Q = +P₀ (current, WRONG for electrons)")
print("=" * 78)

sb_pos = rft.SBend(L, angle_rad, +P0)
sb_pos.set_aperture(10, 10)

report("On-axis",  track(sb_pos, ps_onaxis), 0.0, 0.0)
report("x=1mm",    track(sb_pos, ps_offx),   cos_th, -sin_th/rho)
print()

# ── Test B: P_Q = -P₀ (correct for electrons, q=-1) ──────────────────────
print("=" * 78)
print("Test B: P_Q = -P₀ (CORRECT for electrons)")
print("=" * 78)

sb_neg = rft.SBend(L, angle_rad, -P0)
sb_neg.set_aperture(10, 10)

r_on = track(sb_neg, ps_onaxis)
report("On-axis", r_on, 0.0, 0.0)

r_offx = track(sb_neg, ps_offx)
report("x=1mm", r_offx, cos_th, -sin_th/rho)

r_offy = track(sb_neg, ps_offy)
if r_offy is not None:
    print(f"  {'y=1mm':40s}  y={r_offy[2]:10.4f}  y'={r_offy[3]:10.4f}"
          f"  (exp y=1.0000, drift)")

r_offxp = track(sb_neg, ps_offxp)
report("x'=1mrad", r_offxp, rho*sin_th, cos_th)

# Off-momentum test
delta_p = 0.005  # 0.5%
ps_offmom = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P0 * (1 + delta_p)]])
r_offmom = track(sb_neg, ps_offmom)
# Dispersion: x = ρ(1-cos(θ))×δ in metres → ×1000 for mm
x_disp_exp = rho * (1 - cos_th) * delta_p * 1000  # mm
if r_offmom is not None:
    print(f"  {'δ=+0.5%':40s}  x={r_offmom[0]:10.4f}  "
          f"(exp x={x_disp_exp:.4f} mm dispersion)")
    print(f"  {'':40s}  P={r_offmom[5]:10.4f}  (input={P0*(1+delta_p):.4f})")

# Negative delta
ps_offmom_neg = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P0 * (1 - delta_p)]])
r_offmom_neg = track(sb_neg, ps_offmom_neg)
x_disp_neg = rho * (1 - cos_th) * (-delta_p) * 1000
if r_offmom_neg is not None:
    print(f"  {'δ=-0.5%':40s}  x={r_offmom_neg[0]:10.4f}  "
          f"(exp x={x_disp_neg:.4f})")
print()

# Check accuracy
if r_on is not None and r_offx is not None:
    err_on = abs(r_on[0])
    err_off_x = abs(r_offx[0] - cos_th)
    err_off_xp = abs(r_offx[1] - (-sin_th/rho))
    print(f"  Accuracy: on-axis |x| = {err_on:.2e} mm")
    print(f"            x=1mm  |Δx| = {err_off_x:.2e},  |Δx'| = {err_off_xp:.2e}")
    if err_on < 0.01 and err_off_x < 0.01:
        print("  *** P_Q = -P₀ WORKS for sector bend body ***")
    else:
        print("  *** Still broken ***")
print()

# ── Test C: E1/E2 via constructor with P_Q = -P₀ ─────────────────────────
print("=" * 78)
print("Test C: Edge angles via constructor (P_Q = -P₀)")
print("=" * 78)

# Sector bend (E1=E2=0) vs rectangular bend (E1=E2=angle/2)
sb_sector = rft.SBend(L, angle_rad, -P0, 0, 0)
sb_sector.set_aperture(10, 10)

sb_rect = rft.SBend(L, angle_rad, -P0, angle_rad/2, angle_rad/2)
sb_rect.set_aperture(10, 10)

# Edge angles like our DPW: E1 = wedge_angle at entrance, E2 = wedge at exit
# For symmetric dipole with equal wedges:
wedge = np.radians(10.0)  # 10° wedge
sb_wedge = rft.SBend(L, angle_rad, -P0, wedge, wedge)
sb_wedge.set_aperture(10, 10)

r_sec = track(sb_sector, ps_offx)
r_rect = track(sb_rect, ps_offx)
r_wedge = track(sb_wedge, ps_offx)

report("Sector (E1=E2=0)", r_sec, cos_th, -sin_th/rho)
if r_rect is not None:
    report("Rect (E1=E2=θ/2)", r_rect)
if r_wedge is not None:
    report("Wedge (E1=E2=10°)", r_wedge)

# Compare y-plane (edge effects should change y focusing)
r_sec_y = track(sb_sector, ps_offy)
r_rect_y = track(sb_rect, ps_offy)
r_wedge_y = track(sb_wedge, ps_offy)

if all(r is not None for r in [r_sec_y, r_rect_y, r_wedge_y]):
    print(f"  Y-plane:")
    print(f"    Sector:  y'={r_sec_y[3]:10.6f}")
    print(f"    Rect:    y'={r_rect_y[3]:10.6f}")
    print(f"    Wedge:   y'={r_wedge_y[3]:10.6f}")
    has_edge_effect = abs(r_rect_y[3] - r_sec_y[3]) > 1e-6
    print(f"    Edge effect: {'YES' if has_edge_effect else 'NO'} "
          f"(Δy' = {r_rect_y[3] - r_sec_y[3]:.6f})")

# X-plane edge effects
if all(r is not None for r in [r_sec, r_rect]):
    dx = r_rect[0] - r_sec[0]
    dxp = r_rect[1] - r_sec[1]
    print(f"  X-plane edge Δ: Δx={dx:.6f}, Δx'={dxp:.6f}")
print()

# ── Test D: E1/E2 via setters (in case constructor differs) ───────────────
print("=" * 78)
print("Test D: Edge angles via setters (P_Q = -P₀)")
print("=" * 78)

sb_set_edge = rft.SBend(L, angle_rad, -P0)
sb_set_edge.set_E1(angle_rad/2)
sb_set_edge.set_E2(angle_rad/2)
sb_set_edge.set_aperture(10, 10)

r_set = track(sb_set_edge, ps_offx)
report("Setter E1=E2=θ/2", r_set)

if r_rect is not None and r_set is not None:
    match = abs(r_set[0] - r_rect[0]) < 1e-6
    print(f"  Constructor vs setter: {'MATCH' if match else 'DIFFER'} "
          f"(Δx = {r_set[0] - r_rect[0]:.2e})")
print()

# ── Test E: hgap/fint fringe field parameters ─────────────────────────────
print("=" * 78)
print("Test E: Fringe field parameters (hgap, fint)")
print("=" * 78)

sb_fringe = rft.SBend(L, angle_rad, -P0, angle_rad/2, angle_rad/2)
sb_fringe.set_hgap(0.00724)  # UH MkV half-gap
sb_fringe.set_fint(0.5)
sb_fringe.set_aperture(10, 10)

r_fringe = track(sb_fringe, ps_offx)
report("With hgap/fint", r_fringe)
if r_rect is not None and r_fringe is not None:
    dy_p = abs(r_fringe[3]) - abs(r_rect[3]) if r_fringe is not None else 0
    print(f"  Fringe effect on x': Δx' = {r_fringe[1] - r_rect[1]:.6f} mrad")
print()

# ── Test F: DPW-DPH-DPW triplet as single SBend ──────────────────────────
print("=" * 78)
print("Test F: DPW-DPH-DPW triplet")
print("=" * 78)
# UH MkV typical: entrance wedge, sector body, exit wedge
# This is what we'd consolidate each triplet into
e1_deg = 11.25  # example wedge
e2_deg = 11.25
e1_rad = np.radians(e1_deg)
e2_rad = np.radians(e2_deg)

sb_triplet = rft.SBend(L, angle_rad, -P0, e1_rad, e2_rad)
sb_triplet.set_aperture(10, 10)

r_triplet = track(sb_triplet, ps_offx)
report(f"SBend(E1={e1_deg}°,E2={e2_deg}°)", r_triplet)

# Compare with manual thin-lens edge kick + drift (current workaround)
# Edge kick: K1L = -K0 * tan(e_angle)
K1L_entrance = -K0 * np.tan(e1_rad)
K1L_exit = -K0 * np.tan(e2_rad)
print(f"  Thin-lens K1L: entrance={K1L_entrance:.6f}, exit={K1L_exit:.6f}")
print()

# ── Test G: Negative angle (chicane dipoles bend opposite) ────────────────
print("=" * 78)
print("Test G: Negative bending angle (chicane return dipoles)")
print("=" * 78)

sb_neg_angle = rft.SBend(L, -angle_rad, -P0)
sb_neg_angle.set_aperture(10, 10)

r_neg = track(sb_neg_angle, ps_offx)
cos_neg = np.cos(-angle_rad)  # = cos(θ) (even function)
sin_neg = np.sin(-angle_rad)  # = -sin(θ)
rho_neg = 1.0 / (-angle_rad / L)  # negative rho
report("Neg angle, on-momentum", r_neg, cos_neg, -sin_neg/rho_neg)

r_neg_on = track(sb_neg_angle, ps_onaxis)
report("Neg angle, on-axis", r_neg_on, 0.0, 0.0)
print()

# ── Test H: set_nsteps variation ──────────────────────────────────────────
print("=" * 78)
print("Test H: nsteps variation (P_Q = -P₀)")
print("=" * 78)

for ns in [1, 10, 100]:
    sb_ns = rft.SBend(L, angle_rad, -P0)
    sb_ns.set_nsteps(ns)
    sb_ns.set_aperture(10, 10)
    r_ns = track(sb_ns, ps_offx)
    report(f"nsteps={ns}", r_ns, cos_th, -sin_th/rho)
print()

# ══════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════
print("=" * 78)
print("SUMMARY")
print("=" * 78)

sb_final = rft.SBend(L, angle_rad, -P0)
sb_final.set_aperture(10, 10)

r_final_on = track(sb_final, ps_onaxis)
r_final_off = track(sb_final, ps_offx)

if r_final_on is not None and r_final_off is not None:
    err_on = abs(r_final_on[0])
    err_x = abs(r_final_off[0] - cos_th)
    err_xp = abs(r_final_off[1] - (-sin_th/rho))
    body_ok = err_on < 0.01 and err_x < 0.01 and err_xp < 0.1

    sb_e = rft.SBend(L, angle_rad, -P0, angle_rad/2, angle_rad/2)
    sb_e.set_aperture(10, 10)
    r_e = track(sb_e, ps_offx)
    r_s = track(sb_final, ps_offx)
    edge_ok = r_e is not None and r_s is not None and abs(r_e[0] - r_s[0]) > 1e-6

    print(f"  Body focusing:  {'OK' if body_ok else 'BROKEN'}")
    print(f"  Edge kicks:     {'OK' if edge_ok else 'BROKEN'}")
    print(f"  Dispersion:     check Test B δ results above")

    if body_ok and edge_ok:
        print()
        print("  *** FIX: Use P_Q = -P₀ (negative rigidity for electrons) ***")
        print("  *** E1/E2 via constructor provide edge kicks ***")
        print("  *** DPW-DPH-DPW triplets can be consolidated into single SBend ***")
        print("  → Proceed to Phase 1.2")
    elif body_ok:
        print()
        print("  *** PARTIAL FIX: P_Q = -P₀ fixes body, but E1/E2 broken ***")
        print("  → Keep thin-quad DPW workaround, but use P_Q = -P₀ for DPH")
        print("  → Proceed to Phase 1.2 with hybrid approach")
    else:
        print()
        print("  *** P_Q = -P₀ does NOT fix SBend ***")
        print("  → Hypothesis rejected, keep current workaround (MSE ≈ 4)")
else:
    print("  Particles lost even with P_Q = -P₀")
    print("  → Hypothesis rejected")
