#!/usr/bin/env python3
"""Deep diagnostic of RF-Track SBend. Safe version (no Volume/Multipole)."""

import numpy as np
import RF_Track as rft

me = 0.511
E0 = 40.0
P0 = np.sqrt((E0 + me)**2 - me**2)

L = 0.20
angle_deg = 22.5
angle_rad = np.radians(angle_deg)
K0 = angle_rad / L
rho = 1.0 / K0
cos_th = np.cos(angle_rad)
sin_th = np.sin(angle_rad)

print(f"P0 = {P0:.6f} MeV/c, L = {L} m, angle = {angle_rad:.6f} rad")
print(f"K0 = {K0:.6f} 1/m, ρ = {rho:.6f} m")
print(f"Expected Bρ = P0/c = {P0/299.792458:.6f} T·m")
print(f"Expected B = {P0/(299.792458*rho):.6f} T")
print()


def track_lat(elem, ps, P_ref=P0, charge=-1.0):
    lat = rft.Lattice()
    lat.append(elem)
    b = rft.Bunch6d(me, charge, P_ref, ps.copy())
    tracked = lat.track(b)
    out = np.array(tracked.get_phase_space())
    if out.size == 0:
        return None
    return out[0]


# ── Query internal SBend parameters ──────────────────────────────────────
print("=" * 78)
print("Internal SBend parameters")
print("=" * 78)

for pq_label, pq_val in [("P_Q = +P0", +P0), ("P_Q = -P0", -P0)]:
    sb = rft.SBend(L, angle_rad, pq_val)
    print(f"\n  {pq_label}:")
    print(f"    get_angle()   = {sb.get_angle():.6f} rad")
    print(f"    get_h()       = {sb.get_h():.6f} 1/m  (exp: {K0:.6f})")
    try:
        print(f"    get_Brho()    = {sb.get_Brho():.6f} T·m  (exp: {P0/299.792458:.6f})")
    except Exception as e:
        print(f"    get_Brho()    = {e}")
    print(f"    get_Bfield()  = {sb.get_Bfield():.6f} T  (exp: {P0/(299.792458*rho):.6f})")
    try:
        print(f"    get_K1()      = {sb.get_K1():.6f}")
    except Exception as e:
        print(f"    get_K1()      = {e}")

# ── Check δ via get_phase_space ──────────────────────────────────────────
print("\n" + "=" * 78)
print("δ computation check (before and after tracking)")
print("=" * 78)

ps_on = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P0]])
ps_off = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, P0]])

for pq_label, pq_val in [("P_Q = +P0", +P0), ("P_Q = -P0", -P0)]:
    print(f"\n  {pq_label}:")
    sb = rft.SBend(L, angle_rad, pq_val)
    sb.set_aperture(10, 10)

    b = rft.Bunch6d(me, -1.0, P0, ps_on.copy())
    d_before = np.array(b.get_phase_space('%d'))
    P_before = np.array(b.get_phase_space('%P'))
    print(f"    Before: P={P_before[0,0]:.4f}, δ={d_before[0,0]:.4f} ‰")

    lat = rft.Lattice()
    lat.append(sb)
    tracked = lat.track(b)
    out = np.array(tracked.get_phase_space())
    d_after = np.array(tracked.get_phase_space('%d'))
    P_after = np.array(tracked.get_phase_space('%P'))
    if out.size > 0:
        print(f"    After:  x={out[0,0]:.4f}, P={P_after[0,0]:.4f}, δ={d_after[0,0]:.4f} ‰")

# ── P_Q sweep ────────────────────────────────────────────────────────────
print("\n" + "=" * 78)
print("P_Q sweep: seeking correct body tracking")
print("=" * 78)

ps_offx = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, P0]])

for pq_mult in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0,
                 -0.001, -0.01, -0.1, -0.5, -1.0, -2.0, -5.0, -10.0]:
    pq = P0 * pq_mult
    sb = rft.SBend(L, angle_rad, pq)
    sb.set_aperture(100, 100)
    r = track_lat(sb, ps_offx)
    if r is not None:
        err = abs(r[0] - cos_th)
        marker = " <<<" if err < 0.1 else ""
        print(f"  P_Q/P0={pq_mult:8.3f}  x={r[0]:12.4f}  "
              f"x'={r[1]:12.4f}  |err_x|={err:.2e}{marker}")
    else:
        print(f"  P_Q/P0={pq_mult:8.3f}  LOST")

# ── Setter-only construction ─────────────────────────────────────────────
print("\n" + "=" * 78)
print("Setter-only SBend construction")
print("=" * 78)

B_expected = P0 / (299.792458 * rho)

setups = [
    ("set_K0 + set_Bfield",
     lambda sb: (sb.set_K0(K0), sb.set_Bfield(B_expected))),
    ("set_h + set_Bfield",
     lambda sb: (sb.set_h(K0), sb.set_Bfield(B_expected))),
    ("set_K0 + set_P_over_Q(+P0)",
     lambda sb: (sb.set_K0(K0), sb.set_P_over_Q(P0))),
    ("set_K0 + set_P_over_Q(-P0)",
     lambda sb: (sb.set_K0(K0), sb.set_P_over_Q(-P0))),
]

for label, setup_fn in setups:
    sb = rft.SBend()
    sb.set_length(L)
    setup_fn(sb)
    sb.set_aperture(10, 10)
    r = track_lat(sb, ps_offx)
    if r is not None:
        print(f"  {label:40s}  x={r[0]:12.4f}  x'={r[1]:12.4f}")
    else:
        print(f"  {label:40s}  LOST")

# ── Proton test ──────────────────────────────────────────────────────────
print("\n" + "=" * 78)
print("Proton test (positive charge)")
print("=" * 78)

mp = 938.272
E0_p = 200
P0_p = np.sqrt((E0_p + mp)**2 - mp**2)

# For proton: P/q = P0_p / (+1) = P0_p
sb_p = rft.SBend(L, angle_rad, P0_p)
sb_p.set_aperture(100, 100)
print(f"  Proton: E0={E0_p} MeV, P0={P0_p:.4f} MeV/c")

ps_p_on = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P0_p]])
ps_p_off = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, P0_p]])

r_on = track_lat(sb_p, ps_p_on, P_ref=P0_p, charge=+1.0)
r_off = track_lat(sb_p, ps_p_off, P_ref=P0_p, charge=+1.0)
if r_on is not None:
    print(f"  On-axis:  x={r_on[0]:12.4f} (exp 0)")
if r_off is not None:
    print(f"  x=1mm:    x={r_off[0]:12.4f} (exp {cos_th:.6f})")

# ── nsteps test ──────────────────────────────────────────────────────────
print("\n" + "=" * 78)
print("nsteps variation")
print("=" * 78)

for ns in [0, 1, 2, 5, 10, 50, 100]:
    sb = rft.SBend(L, angle_rad, P0)
    sb.set_nsteps(ns)
    sb.set_aperture(100, 100)
    r = track_lat(sb, ps_offx)
    if r is not None:
        print(f"  nsteps={ns:3d}  x={r[0]:12.4f}  x'={r[1]:12.4f}")
    else:
        print(f"  nsteps={ns:3d}  LOST")

# ── Angle sweep (is the error proportional to angle?) ────────────────────
print("\n" + "=" * 78)
print("Angle sweep (P_Q = +P0)")
print("=" * 78)

for angle_d in [1, 5, 10, 22.5, 45, 90]:
    ang = np.radians(angle_d)
    sb = rft.SBend(L, ang, P0)
    sb.set_aperture(1000, 1000)
    r_on_a = track_lat(sb, np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P0]]))
    r_off_a = track_lat(sb, np.array([[1.0, 0.0, 0.0, 0.0, 0.0, P0]]))
    if r_on_a is not None and r_off_a is not None:
        c_exp = np.cos(ang)
        body_resp = r_off_a[0] - r_on_a[0]
        print(f"  θ={angle_d:5.1f}°  on-axis x={r_on_a[0]:12.4f}  "
              f"body Δx={body_resp:8.4f} (exp cos(θ)={c_exp:.4f})  "
              f"ratio={body_resp/c_exp:.4f}")
    else:
        print(f"  θ={angle_d:5.1f}°  LOST")
