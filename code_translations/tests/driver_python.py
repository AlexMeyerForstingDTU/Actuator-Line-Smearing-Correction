"""Test driver for the Python translation.

Usage: driver_python.py <case>
Runs one test case in the current folder and writes results.txt.
The correction itself writes the diagnostics files test.nwRR.blBB.
The test cases and inputs repeat test_common.f90.
"""
import os
import sys

import numpy as np

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from al_nw_smearing_correction import (init_nw_smearing_correction, nw_smearing_correction,
                                       smearing_factor, nw_hr_components, nw_phi_components,
                                       nw_xy_components, nw_circulation, nw_calc_relax,
                                       nw_interp1d, nw2vstn, vxyz_vstn)

NS, NV, LEN2D = 63, 64, 181
D2R, R2D = np.pi / 180, 180 / np.pi
WIND, A_IND, RHO = 8.0, 0.3, 1.225
OMEGA = 9.2 * 2 * np.pi / 60

# Chord and twist of the blade at the table radii
RTAB = [2.8667, 5.6, 8.3333, 11.75, 15.85, 19.95, 24.05, 28.15, 32.25, 36.35, 40.45, 44.55,
        48.65, 52.75, 56.1667, 58.9, 61.6333]
TTAB = [13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188,
        3.125, 2.319, 1.526, 0.863, 0.370, 0.106]
CTAB = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010,
        2.764, 2.518, 2.313, 2.086, 1.419]


def get_case(name):
    """Parameters of a named test case: frozen inflow marched in time."""
    c = {"eps": 10.0,            # Smearing length [m]
         "cone": 0.0,            # Cone [deg]
         "dt": 0.1,              # Time step [s]
         "n_steps": 300,
         "l1point": False,       # Ang Li correction
         "bent": False,          # Pre-bent, swept blade with out-of-plane offsets
         "root_at_zero": False,  # First trailing point at r = 0
         "max_resid": 1e-6,      # Convergence limit [m/s]
         "pitch_step_at": 0,     # From this step on: pitch = 2 deg (0 = never)
         "n_rotors": 1,
         "n_blades": 1,
         "pitch_22": 0.0}        # Pitch of rotor 2, blade 2 [deg]
    if name == "eps10_dt001":
        c["dt"] = 0.01
        c["n_steps"] = 3000
    elif name == "eps10":
        pass
    elif name == "cone25":
        c["cone"] = 2.5
    elif name == "eps25":
        c["eps"] = 2.5
    elif name == "nocorr":
        c["eps"] = 0.01
    elif name == "bent":
        c["cone"] = 2.5
        c["bent"] = True
        c["l1point"] = True
    elif name == "pitchstep":
        c["n_steps"] = 400
        c["pitch_step_at"] = 301
    elif name == "noconv":
        c["n_steps"] = 2
        c["max_resid"] = 0.0
    elif name == "multi":
        c["n_steps"] = 100
        c["n_rotors"] = 2
        c["n_blades"] = 2
        c["pitch_22"] = 1.0
    elif name == "rootzero":
        c["n_steps"] = 1
        c["root_at_zero"] = True
    else:
        sys.exit("Error: unknown test case.")
    return c


def table_interp(x, xt, yt):
    """Piecewise linear table lookup, constant beyond the ends."""
    if x <= xt[0]:
        return yt[0]
    if x >= xt[-1]:
        return yt[-1]
    k = 1
    while x > xt[k]:
        k = k + 1
    return yt[k - 1] + (x - xt[k - 1]) / (xt[k] - xt[k - 1]) * (yt[k] - yt[k - 1])


def build_blade(c):
    """Geometry, polars and frozen cylindrical inflow of the test blade (a dict)."""
    b = {}
    # 1. Trailing points (uniform) and control points in between
    r0 = 0.0 if c["root_at_zero"] else 1.5
    b["rv"] = np.array([r0 + (63 - r0) * k / (NV - 1) for k in range(NV)])
    b["rc"] = 0.5 * (b["rv"][0:NS] + b["rv"][1:NV])
    b["chord"] = np.array([table_interp(r, RTAB, CTAB) for r in b["rc"]])
    b["twist"] = np.array([table_interp(r, RTAB, TTAB) for r in b["rc"]])

    # 2. Blade shape: straight, or pre-bent and swept with out-of-plane offsets
    b["fbend"] = np.zeros(NS)
    b["sweep"] = np.zeros(NS)
    b["dsweep"] = np.zeros(NS)
    b["cyl_cp"] = np.zeros((NS, 3))
    b["cyl_vtp"] = np.zeros((NV, 3))
    b["cyl_cp"][:, 0] = b["rc"]
    b["cyl_vtp"][:, 0] = b["rv"]
    if c["bent"]:
        t = b["rc"] / 63
        b["fbend"] = 4 * t * t
        b["sweep"] = 3 * t * t * t
        b["dsweep"] = 0.5 * b["sweep"]
        b["cyl_cp"][:, 2] = -1.5 * t * t
        t = b["rv"] / 63
        b["cyl_vtp"][:, 2] = -1.5 * t * t

    # 3. Polars: cylinder inboard of 10 m, smooth thin-airfoil-like outboard
    b["alpha_2d"] = np.zeros((LEN2D, NS))
    b["cl_2d"] = np.zeros((LEN2D, NS))
    b["cd_2d"] = np.zeros((LEN2D, NS))
    for l in range(NS):
        b["alpha_2d"][:, l] = -90 + np.arange(LEN2D)
        if b["rc"][l] < 10:
            b["cl_2d"][:, l] = 0.0
            b["cd_2d"][:, l] = 0.4
        else:
            b["cl_2d"][:, l] = 1.4 * np.tanh(2 * np.pi * (b["alpha_2d"][:, l] + 3) * D2R / 1.4)
            b["cd_2d"][:, l] = 0.01

    # 4. Frozen "CFD" inflow in the cylindrical frame: uniform a = 0.3
    b["vcyl"] = np.zeros((NS, 3))
    b["vcyl"][:, 2] = WIND * (1 - A_IND)
    return b


def cfd_vstn(b, cone, pitch):
    """CFD velocities in section axes for a given cone and pitch."""
    return vxyz_vstn(b["vcyl"], 0, 0, 0, cone, pitch, b["sweep"], b["fbend"])


def power_thrust(b, Fstn, factor):
    """Rotor power and thrust from the section forces of one blade."""
    P = 0.0
    T = 0.0
    for l in range(NS):
        P = P + factor * Fstn[l, 1] * OMEGA * b["rc"][l] * (b["rv"][l + 1] - b["rv"][l])
        T = T + factor * Fstn[l, 2] * (b["rv"][l + 1] - b["rv"][l])
    return P, T


def write_values(out, name, values):
    """Write one named result line: name v1 v2 ..."""
    out.write(name)
    for value in np.ravel(values):
        out.write(" %.16e" % value)
    out.write("\n")


def run_units(out):
    """Unit tests of the individual routines."""
    # 1. Smearing factor, with and without an out-of-plane offset
    write_values(out, "smearing_factor", [smearing_factor(60.0, 0.1, 0.0, 0.0, 10.0),
                                          smearing_factor(60.0, 0.1, 2.0, 0.3, 10.0),
                                          smearing_factor(10.0, -0.5, 1.0, 0.5, 2.5)])

    # 2. Geometry constants, r_vtp = 10 m: phi, phi_s, dx, dy, a(1:4) per h/r
    values = []
    for h_r in [0.1, 0.3, 0.5, 0.7, 0.8, -0.5, -3.0]:
        phi, phi_s, dx, dy, a = nw_hr_components(10 * h_r, h_r)
        values = values + [phi, phi_s, dx, dy, a[0], a[1], a[2], a[3]]
    write_values(out, "hr_components", values)

    # 3. Phi* for four (h/r, helix angle) pairs
    values = []
    for h_r, angle in [(0.1, 5.0), (0.3, 30.0), (-0.5, 10.0), (-3.0, 60.0)]:
        phi, phi_s, dx, dy, a = nw_hr_components(10 * h_r, h_r)
        values.append(nw_phi_components(h_r, a, angle * D2R, phi, phi_s))
    write_values(out, "phi_components", values)

    # 4. Indicial functions: one step, then the steady state of 2000 steps
    X, Y = nw_xy_components(0.09, 0.0123, -0.0008125, 0.05, 2.0, 0.3, -0.02)
    write_values(out, "xy_single", [X, Y])
    X = 0.0
    Y = 0.0
    for n in range(2000):
        X, Y = nw_xy_components(0.09, 0.0123, -0.0008125, 0.05, 2.0, X, Y)
    write_values(out, "xy_steady", [X, Y])

    # 5. Interpolation, inside and beyond both ends
    x_tab = np.array([0.0, 1.0, 2.0, 3.0])
    y_tab = np.array([0.0, 10.0, 15.0, 30.0])
    values = []
    for x in [-1, 0, 0.5, 2, 2.5, 3, 4]:
        values.append(nw_interp1d(x_tab, y_tab, x))
    write_values(out, "interp1d", values)

    # 6. nw2vstn called repeatedly on a coned blade (note Sect. 4.1 and 6.1)
    values = []
    for n in range(15):
        v1 = nw2vstn(np.array([0.0]), np.array([-1.0]), 2.5, 0.0, np.array([0.0]), np.array([0.0]))
        values.append(v1[0, 2])
    write_values(out, "nw2vstn_repeat", values)
    write_values(out, "nw2vstn_repeat_last", v1[0, :])

    # 7. nw2vstn and the full rotation chain with all angles non-zero
    v2 = nw2vstn(np.array([1.5, -0.3]), np.array([-0.7, 2.0]), 2.5, 3.0,
                 np.array([4.0, -2.0]), np.array([-2.0, 5.0]))
    write_values(out, "nw2vstn_general", v2)
    v = vxyz_vstn(np.array([[1.0, -2.0, 0.5]]), 10.0, 5.0, 30.0, 2.5, 3.0, 4.0, -2.0)
    write_values(out, "vxyz_vstn", v)

    # 8. Circulation of the test blade without correction (cone 0, pitch 0)
    b = build_blade(get_case("eps10"))
    un = np.full(NS, WIND * (1 - A_IND))
    ut = -(OMEGA * b["rc"])
    vrel = np.sqrt(un**2 + ut**2)
    alpha = np.arctan2(un, -ut) * R2D - b["twist"]
    G_hist = np.zeros(NS)
    dyn_hist = np.zeros((NS, 3))
    G_steady, G, dG, dyn = nw_circulation(0.1, vrel, alpha, G_hist, dyn_hist, b["alpha_2d"],
                                          b["cl_2d"], b["chord"])
    write_values(out, "circulation_gamma_steady", G_steady)
    write_values(out, "circulation_gamma", G)
    write_values(out, "circulation_dgamma", dG)
    # Steady state of the time filter: Gamma -> Gamma_steady when sum(A_j) = 1
    for n in range(300):
        G_steady, G, dG, dyn = nw_circulation(1.0, vrel, alpha, G_hist, dyn_hist, b["alpha_2d"],
                                              b["cl_2d"], b["chord"])
        G_hist = G_steady
        dyn_hist = dyn
    ratio = 0.0
    for l in range(NS):
        if abs(G_steady[l]) > 1e-3:
            ratio = max(ratio, abs(G[l] / G_steady[l] - 1))
    write_values(out, "circulation_ratio", [ratio])

    # 9. Relaxation factor of the test blade at two time steps
    phi_pair = np.zeros((NV, NS))
    dx_pair = np.zeros((NV, NS))
    dy_pair = np.zeros((NV, NS))
    for k in range(NV):
        for l in range(NS):
            h = b["rv"][k] - b["rc"][l]
            phi_pair[k, l], phi_s, dx_pair[k, l], dy_pair[k, l], a = \
                nw_hr_components(h, h / b["rv"][k])
    vrel_vtp = np.zeros(NV)
    vrel_vtp[0] = vrel[0]
    vrel_vtp[NV - 1] = vrel[NS - 1]
    for k in range(1, NV - 1):
        vrel_vtp[k] = nw_interp1d(b["rc"], vrel, b["rv"][k])
    values = []
    for dt in [0.01, 0.1]:
        dbeta = (vrel_vtp * dt) / b["rv"]
        values.append(nw_calc_relax(vrel, alpha, un, -ut, b["chord"], dt, dx_pair, dy_pair, dbeta,
                                    phi_pair, 0.2))
    write_values(out, "calc_relax", values)


def run_case(c, out):
    """March one test case in time and write its results."""
    # 1. Inputs and the correction state
    b = build_blade(c)
    state = init_nw_smearing_correction(c["n_rotors"], c["n_blades"], NS, "test")
    state.nw_max_resid = c["max_resid"]

    # 2. March in time; F[i][j] holds the forces of rotor i + 1, blade j + 1
    F = [[None, None], [None, None]]
    P_series = []
    T_series = []
    for it in range(1, c["n_steps"] + 1):
        time = it * c["dt"]
        pitch = 0.0
        if c["pitch_step_at"] > 0 and it >= c["pitch_step_at"]:
            pitch = 2.0
        for i in range(1, c["n_rotors"] + 1):
            for j in range(1, c["n_blades"] + 1):
                p = pitch
                if i == 2 and j == 2:
                    p = c["pitch_22"]
                vstn = cfd_vstn(b, c["cone"], p)
                F[i - 1][j - 1] = nw_smearing_correction(
                    state, i, j, time, c["dt"], RHO, b["cyl_cp"], b["cyl_vtp"], vstn, b["vcyl"],
                    c["eps"], OMEGA, b["chord"], b["twist"], p, b["cl_2d"], b["cd_2d"],
                    b["alpha_2d"], c["cone"], b["fbend"], b["sweep"], b["dsweep"], c["l1point"])
        if c["pitch_step_at"] > 0 and it > c["n_steps"] - 110:
            P, T = power_thrust(b, F[0][0], 3)
            P_series.append(P)
            T_series.append(T)

    # 3. Results
    if c["root_at_zero"]:
        write_values(out, "kaya_after", [1.0 if state.kaya else 0.0])
    elif c["n_rotors"] > 1:
        P_blade = []
        T_blade = []
        for i in range(2):
            for j in range(2):
                P, T = power_thrust(b, F[i][j], 1)
                P_blade.append(P)
                T_blade.append(T)
        write_values(out, "P_blade", P_blade)
        write_values(out, "T_blade", T_blade)
    else:
        P, T = power_thrust(b, F[0][0], 3)
        write_values(out, "P", [P])
        write_values(out, "T", [T])
        write_values(out, "Fstn2", F[0][0][:, 1])
        write_values(out, "Fstn3", F[0][0][:, 2])
        write_values(out, "vn_history", state.blade[0][0].vn_history)
        write_values(out, "vt_history", state.blade[0][0].vt_history)
        if c["pitch_step_at"] > 0:
            write_values(out, "P_series", P_series)
            write_values(out, "T_series", T_series)

    # Close the diagnostics files
    for i in range(c["n_rotors"]):
        for j in range(c["n_blades"]):
            state.blade[i][j].fid.close()


if __name__ == "__main__":
    case_name = sys.argv[1]
    with open("results.txt", "w") as results:
        if case_name == "units":
            run_units(results)
        else:
            run_case(get_case(case_name), results)
