#!/usr/bin/env python3
"""SBend workaround: sliced dipole using Corrector + Drift.

Since SBend body tracking is fundamentally broken in RF-Track 2.5.5,
we build dipoles from working elements:
  - Corrector: thin dipole kick (T·mm integrated field)
  - Drift: free-space propagation
  - Quadrupole: edge kicks (thin-lens)

A sector bend is split into N slices, each: half-drift → kick → half-drift.
This is a symplectic split-operator method that reproduces body focusing,
dispersion, and path-length effects to O(L/N)² accuracy.

Author: Eremey Valetov
"""

import numpy as np
import RF_Track as rft

me = 0.511
E0 = 40.0
P0 = np.sqrt((E0 + me)**2 - me**2)
Brho = P0 / 299.792458  # T·m

L = 0.20
angle_deg = 22.5
angle_rad = np.radians(angle_deg)
K0 = angle_rad / L
rho = 1.0 / K0
cos_th = np.cos(angle_rad)
sin_th = np.sin(angle_rad)

print(f"P0 = {P0:.6f} MeV/c, Bρ = {Brho:.6f} T·m")
print(f"L = {L} m, θ = {angle_deg}°, ρ = {rho:.6f} m")
print(f"Expected on-axis: x=0, x'=0")
print(f"Expected x=1mm: x={cos_th:.6f}, x'={-sin_th/rho:.6f}")
print()


def track(lat, ps, P_ref=P0, charge=-1.0):
    b = rft.Bunch6d(me, charge, P_ref, ps.copy())
    tracked = lat.track(b)
    out = np.array(tracked.get_phase_space())
    if out.size == 0:
        return None
    return out[0]


# ── Test Corrector element ───────────────────────────────────────────────
print("=" * 78)
print("Test 0: Corrector element basics")
print("=" * 78)

# Thin corrector with full bend angle
# Kick: Δx' = -Kx / Brho_particle  [sign depends on charge convention]
# For electron: Kx in T·mm, Brho = P/(|q|c) in T·m = P/299.792 in T·m
# Integrated field B·L = ρ × B × θ/B = ρ × θ × B... no, simpler:
# B·L = Brho × θ = (P₀/c) × θ = 0.1351 × 0.3927 = 0.05306 T·m = 53.06 T·mm

BdL = Brho * angle_rad * 1000  # T·mm (physical integrated field)
# RF-Track Corrector normalisation: pass BdL/P₀, not raw BdL
# (same pattern as Quadrupole set_strength taking k1*L, not P/q*k1*L)
Kx = BdL / P0
print(f"  Integrated dipole field: {BdL:.4f} T·mm, Kx = BdL/P0 = {Kx:.6f}")

# Test: thin corrector should give Δx' = angle_rad
corr = rft.Corrector(0, Kx, 0)
corr.set_aperture(1000, 1000)
lat_c = rft.Lattice()
lat_c.append(corr)

ps_test = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P0]])
r = track(lat_c, ps_test)
if r is not None:
    kick_mrad = r[1]  # x' after thin kick
    kick_rad = kick_mrad / 1000
    print(f"  Thin corrector kick: {kick_mrad:.6f} mrad = {kick_rad:.6f} rad "
          f"(expected ≈ {angle_rad:.6f} rad)")
    if abs(kick_rad) > 0:
        print(f"  Sign: {'positive' if kick_rad > 0 else 'negative'}")
        # Check momentum dependence: off-momentum particle
        ps_offp = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P0 * 1.005]])
        r_off = track(lat_c, ps_offp)
        if r_off is not None:
            kick_off = r_off[1] / 1000
            ratio = kick_off / kick_rad
            print(f"  P+0.5% kick ratio: {ratio:.6f} (expected {1/1.005:.6f} ≈ chromatic)")
else:
    print("  Corrector: LOST")

# Try with negative Kx for opposite sign
corr_neg = rft.Corrector(0, -Kx, 0)
corr_neg.set_aperture(1000, 1000)
lat_cn = rft.Lattice()
lat_cn.append(corr_neg)
r_neg = track(lat_cn, ps_test)
if r_neg is not None:
    print(f"  Negative Kx kick: {r_neg[1]:.6f} mrad")
print()


# ── Sliced sector bend ───────────────────────────────────────────────────
def build_sliced_dipole(L, angle_rad, N, aperture=1000):
    """Build a sector bend from N Corrector+Drift slices (split-operator).

    Each slice: Drift(L/2N) → Corrector(θ/N) → Drift(L/2N)

    Returns RF-Track Lattice.
    """
    ds = L / N       # drift per slice
    dtheta = angle_rad / N
    dKx = Brho * dtheta * 1000 / P0  # normalised kick per slice

    lat = rft.Lattice()
    for i in range(N):
        # Half drift
        d1 = rft.Drift(ds / 2)
        d1.set_aperture(aperture / 1000, aperture / 1000)
        lat.append(d1)

        # Thin kick
        c = rft.Corrector(0, dKx, 0)
        c.set_aperture(aperture / 1000, aperture / 1000)
        lat.append(c)

        # Half drift
        d2 = rft.Drift(ds / 2)
        d2.set_aperture(aperture / 1000, aperture / 1000)
        lat.append(d2)

    return lat


print("=" * 78)
print("Test 1: Sliced sector bend vs analytical")
print("=" * 78)

ps_onaxis = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P0]])
ps_offx = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, P0]])

for N in [1, 5, 10, 20, 50, 100]:
    lat = build_sliced_dipole(L, angle_rad, N)
    r_on = track(lat, ps_onaxis)
    r_off = track(lat, ps_offx)
    if r_on is not None and r_off is not None:
        body_x = r_off[0] - r_on[0]
        body_xp = r_off[1] - r_on[1]
        err_on = abs(r_on[0])
        err_body = abs(body_x - cos_th)
        err_xp = abs(body_xp - (-sin_th / rho))
        print(f"  N={N:3d}: on-axis x={r_on[0]:10.4f}  "
              f"body Δx={body_x:10.6f} (exp {cos_th:.6f}, err={err_body:.2e})  "
              f"Δx'={body_xp:10.6f} (exp {-sin_th/rho:.6f})")
    else:
        print(f"  N={N:3d}: LOST")
print()

# ── Test 2: Dispersion check ─────────────────────────────────────────────
print("=" * 78)
print("Test 2: Dispersion (off-momentum particle)")
print("=" * 78)

delta_p = 0.005
ps_offmom = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P0 * (1 + delta_p)]])
x_disp_exp = rho * (1 - cos_th) * delta_p * 1000  # mm

for N in [10, 50, 100]:
    lat = build_sliced_dipole(L, angle_rad, N)
    r = track(lat, ps_offmom)
    r_ref = track(lat, ps_onaxis)
    if r is not None and r_ref is not None:
        x_disp = r[0] - r_ref[0]  # subtract geometric offset
        err_disp = abs(x_disp - x_disp_exp)
        print(f"  N={N:3d}: x_disp={x_disp:10.6f} mm  "
              f"(exp {x_disp_exp:.6f}, err={err_disp:.4f})")
print()

# ── Test 3: Edge kicks ───────────────────────────────────────────────────
print("=" * 78)
print("Test 3: DPW-DPH-DPW triplet (edges + sliced body)")
print("=" * 78)

wedge_angle = np.radians(11.25)  # half the bend angle
K0_local = angle_rad / L

# Edge kick: thin Quadrupole with K1L = -K0 * tan(wedge_angle)
K1L_edge = -K0_local * np.tan(wedge_angle)

def build_dpw_dph_dpw(L, angle_rad, wedge_e1, wedge_e2, N_body, aperture=1000):
    """Build DPW-DPH-DPW triplet: edge_kick + sliced_body + edge_kick."""
    lat = rft.Lattice()
    ap = aperture / 1000

    # Entrance edge kick
    K1L_e1 = -K0_local * np.tan(wedge_e1)
    q_e1 = rft.Quadrupole()
    q_e1.set_length(1e-10)
    q_e1.set_strength(K1L_e1)
    q_e1.set_aperture(ap, ap)
    lat.append(q_e1)

    # Sliced body
    ds = L / N_body
    dtheta = angle_rad / N_body
    dKx = Brho * dtheta * 1000 / P0  # normalised

    for i in range(N_body):
        d1 = rft.Drift(ds / 2)
        d1.set_aperture(ap, ap)
        lat.append(d1)
        c = rft.Corrector(0, dKx, 0)
        c.set_aperture(ap, ap)
        lat.append(c)
        d2 = rft.Drift(ds / 2)
        d2.set_aperture(ap, ap)
        lat.append(d2)

    # Exit edge kick
    K1L_e2 = -K0_local * np.tan(wedge_e2)
    q_e2 = rft.Quadrupole()
    q_e2.set_length(1e-10)
    q_e2.set_strength(K1L_e2)
    q_e2.set_aperture(ap, ap)
    lat.append(q_e2)

    return lat


ps_offy = np.array([[0.0, 0.0, 1.0, 0.0, 0.0, P0]])
lat_triplet = build_dpw_dph_dpw(L, angle_rad, wedge_angle, wedge_angle, 20)
r_triplet = track(lat_triplet, ps_offx)
r_triplet_y = track(lat_triplet, ps_offy)
if r_triplet is not None:
    print(f"  Triplet (N=20): x={r_triplet[0]:.6f}  x'={r_triplet[1]:.6f}")
if r_triplet_y is not None:
    print(f"  Y-plane:        y={r_triplet_y[2]:.6f}  y'={r_triplet_y[3]:.6f}")
    # Expected: edge kicks focus vertically: Δy' ≈ 2×tan(wedge)/ρ×y
    expected_yp = -2 * np.tan(wedge_angle) / rho * 1.0  # vertical defocusing from sector edges
    # Actually for sector bend edges, vertical kick is: Δy' = -tan(e)/ρ × y (defocusing)
    # Two edges: Δy'_total = -2tan(e)/ρ × y
    print(f"  Expected y' ≈ {expected_yp:.6f} (+ body drift effect)")
print()

# ── Test 4: Convergence study ────────────────────────────────────────────
print("=" * 78)
print("Test 4: Convergence with number of slices")
print("=" * 78)

# Compare sliced dipole transfer matrix with analytical sector bend
# Use multi-particle tracking to extract effective R-matrix

def extract_R_matrix(lat_builder, N):
    """Extract 4×4 transverse R-matrix using basis vectors."""
    lat = lat_builder(N)
    eps = 0.01  # mm or mrad perturbation
    ref = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, P0]])
    r_ref = track(lat, ref)
    if r_ref is None:
        return None

    R = np.zeros((4, 4))
    for j in range(4):
        ps_plus = ref.copy()
        ps_plus[0, j] = +eps
        ps_minus = ref.copy()
        ps_minus[0, j] = -eps
        r_plus = track(lat, ps_plus)
        r_minus = track(lat, ps_minus)
        if r_plus is None or r_minus is None:
            return None
        for i in range(4):
            R[i, j] = (r_plus[i] - r_minus[i]) / (2 * eps)

    return R, r_ref


# Analytical sector bend 4×4 matrix
R_analytical = np.array([
    [cos_th,           rho * sin_th,  0, 0],
    [-sin_th / rho,    cos_th,        0, 0],
    [0,                0,             1, L],
    [0,                0,             0, 1],
])

print("  Analytical R-matrix:")
for i in range(4):
    print(f"    [{R_analytical[i,0]:10.6f} {R_analytical[i,1]:10.6f} "
          f"{R_analytical[i,2]:10.6f} {R_analytical[i,3]:10.6f}]")
print()

for N in [5, 10, 20, 50, 100]:
    result = extract_R_matrix(lambda n: build_sliced_dipole(L, angle_rad, n), N)
    if result is not None:
        R, ref = result
        err = np.max(np.abs(R - R_analytical))
        det = np.linalg.det(R[:2, :2])
        print(f"  N={N:3d}: max|ΔR|={err:.2e}  det(x-block)={det:.8f}  "
              f"R11={R[0,0]:.6f} R12={R[0,1]:.6f}")
    else:
        print(f"  N={N:3d}: extraction failed")
print()

# ══════════════════════════════════════════════════════════════════════════
print("=" * 78)
print("VERDICT")
print("=" * 78)

result = extract_R_matrix(lambda n: build_sliced_dipole(L, angle_rad, n), 20)
if result is not None:
    R, ref = result
    err = np.max(np.abs(R - R_analytical))
    offset = abs(ref[0])  # should be ~0 for on-axis
    print(f"  N=20 slices: max|ΔR| = {err:.2e}, on-axis offset = {offset:.4f} mm")
    if err < 0.01:
        print("  *** SLICED DIPOLE CONVERGES TO ANALYTICAL MATRIX ***")
        print("  → Proceed to Phase 1.2: implement sliced dipoles in rftrackAdapter")
    else:
        print(f"  R-matrix error still significant ({err:.2e})")
        if err < 0.1:
            print("  → Usable with more slices. Proceed with caution.")
        else:
            print("  → Need alternative approach")
