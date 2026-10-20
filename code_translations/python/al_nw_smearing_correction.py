"""Fast smearing correction for the actuator line (Python version).

Translated from al_nw_smearing_correction.f. Method: Meyer Forsting et al.
(2019, 2020), with the near-wake model of Pirrung et al. (2016, 2017a, 2017b)
and the 2020 corrigenda. Equation numbers refer to these papers.

Usage:
    state = init_nw_smearing_correction(n_rotors, n_blades, n_sections, project)
    then every time step, for each rotor i and blade j (both from 1):
    Fstn = nw_smearing_correction(state, i, j, time, dt, density, ...)

Array shapes: sections (n_s), trailing points (n_v = n_s + 1), vortex-section
pairs (n_v, n_s), vectors (n, 3), polars (len_2d, n_s).
"""
from time import process_time

import numpy as np

D2R = np.pi / 180
R2D = 180 / np.pi

# Fit constants of the mixing factor k_Phi, Pirrung et al. (2017b). Row i gives a_hr,i.
# P: h/r > 0, Eq. (14). N: h/r < 0, Eq. (13).
NW_P = np.array([
    [-1.64636754184988, 8.14821474772595, -12.1784861715161, 5.02652654794305, 21.7713111885100],
    [-0.499006672527007, 6.08465268675599, -15.1712005446736, 14.8254135311781, -2.42318681495122],
    [3.90835606394205, -18.7662330478906, 39.1243282362892, -29.4870086106427, -60.2947309276912],
    [-1.60622843516794, 7.42952711674470, -15.8594777074799, 11.6870163960092, -195.060865844227]])
NW_N = np.array([
    [1.01933206359188, -0.135668151345576, 0.395524705910666, 0.0801822997046969, 44.8347503471217],
    [12.9874512870379, 49.9999999999916, 0.00235345218088083, 11.3116074964299, 3935.34323353354],
    [-0.690159136848837, 101.238775997879, -0.00154467589177718, 3.99520376143011, 0.394541495293350],
    [-0.269253409815015, 49.9999999999843, -0.00247975579295973, 0.403642518874381, 1.16610462658971]])


class BladeState:
    """History and geometry constants of one rotor blade."""

    def __init__(self, n_s, n_v):
        self.refresh = True                         # Recompute the geometry constants
        self.cp_loop_ind = np.zeros((n_v, 2), int)  # First and last CP felt by each vortex
        self.cp_loop_ind[:, 1] = n_s - 1
        self.W_history = np.zeros(n_s)
        self.vn_history = np.zeros(n_s)
        self.vt_history = np.zeros(n_s)
        self.Gamma_history = np.zeros(n_s)
        self.dyn_comp_history = np.zeros((n_s, 3))
        self.X_history = np.zeros((n_v, n_s))
        self.Y_history = np.zeros((n_v, n_s))
        self.phi = np.zeros((n_v, n_s))             # Root-corrected Phi_C
        self.phi_s = np.zeros((n_v, n_s))
        self.dx = np.zeros((n_v, n_s))
        self.dy = np.zeros((n_v, n_s))
        self.a = np.zeros((n_v, n_s, 4))            # Fit coefficients a_hr,i
        self.smear_fac = np.zeros((n_v, n_s))
        self.fid = None                             # Diagnostics file


class NWState:
    """Settings and the state of all rotor blades."""

    def __init__(self, n_rotors, n_blades, n_sections):
        self.kaya = True              # Switch: correction on/off
        self.nw_max_resid = 1e-6      # Convergence limit of the correction velocity [m/s]
        self.x_aoa = 0.5              # Chordwise position of the AoA in the CFD [-]
        self.n_s = n_sections
        self.n_v = n_sections + 1
        self.blade = [[BladeState(self.n_s, self.n_v) for j in range(n_blades)]
                      for i in range(n_rotors)]


def init_nw_smearing_correction(n_rotors, n_blades, n_sections, project=""):
    """Create the state of all blades; open the diagnostics files <project>.nwRR.blBB."""
    state = NWState(n_rotors, n_blades, n_sections)
    if project:
        for i in range(n_rotors):
            for j in range(n_blades):
                state.blade[i][j].fid = open(f"{project}.nw{i + 1:02d}.bl{j + 1:02d}", "w")
    return state


def nw_smearing_correction(state, i, j, time, dt, density, cyl_cp, cyl_vtp, vstn_cfd, vcyl_cfd,
                           eps_cp, omega, chord, twist, pitch, cl_2d, cd_2d, alpha_2d,
                           cone, fbend, sweep, dsweep, l1point):
    """Corrected blade forces of rotor i, blade j for one time step.

    The smeared AL forces leave part of the induction out of the CFD velocities;
    the near-wake model supplies that part and corrects the angle of attack.
    Inputs: positions (r, theta, z) in the cylindrical rotor frame; CFD velocities
    at the control points (CP) in section axes (vstn_cfd) and in the cylindrical
    frame (vcyl_cfd); angles in degrees; omega in rad/s.
    Returns Fstn (n_s, 3): spanwise, tangential and normal force per unit length,
    or None when the correction is switched off.
    """
    if not state.kaya:
        return None
    n_s, n_v = state.n_s, state.n_v
    r_cp = cyl_cp[:, 0]
    r_vtp = cyl_vtp[:, 0]
    # Important: r = 0 at the first trailing point is a singularity
    if r_vtp[0] < 0.1:
        print("!!!! WARNING smearing correction switched off !!!!")
        print("First vortex trailing point is zero, please modify")
        state.kaya = False
        return None
    cpu_t1 = process_time()
    b = state.blade[i - 1][j - 1]

    # 1. Preamble: start from the previous time step
    relax_safety = 0.2
    recompute_relax = True
    W = b.W_history.copy()
    vn = b.vn_history.copy()
    vt = b.vt_history.copy()
    scone = np.sin((cone + fbend) * D2R)
    v_rot = -(omega * r_cp)                    # Rotational speed in section axes
    vn_rot = -v_rot * np.sin(-pitch * D2R)
    vt_rot = v_rot * np.cos(-pitch * D2R) * np.cos(dsweep * D2R)

    # 2. Geometry constants: first call, or after a refresh
    if b.refresh:
        core_cutoff = 1.82 * eps_cp            # 99th percentile of the Gaussian
        b.cp_loop_ind[:, 0] = 0
        b.cp_loop_ind[:, 1] = n_s - 1
        for k in range(n_v):
            for l in range(n_s):
                h = r_vtp[k] - r_cp[l]
                h_r = h / r_vtp[k]
                b.phi[k, l], b.phi_s[k, l], b.dx[k, l], b.dy[k, l], b.a[k, l, :] = \
                    nw_hr_components(h, h_r)
                b.smear_fac[k, l] = smearing_factor(r_vtp[k], h_r, 0, 0, eps_cp)
                # CPs felt by vortex k: inboard and outboard ones within the cut-off
                if h > 0:
                    if abs(h) > core_cutoff:
                        b.cp_loop_ind[k, 0] = l + 1
                elif abs(h) < core_cutoff:
                    b.cp_loop_ind[k, 1] = l
        b.refresh = False

    # 3. Helix angles from the CFD only (MF2019 step 2), fixed during the sub-iteration
    phi_cp_cfd = np.arctan2(vcyl_cfd[:, 2], -(v_rot + vcyl_cfd[:, 1]))
    for l in range(n_s):
        if phi_cp_cfd[l] > np.pi / 2:
            phi_cp_cfd[l] = np.pi - phi_cp_cfd[l]
    phi_vtp = np.zeros(n_v)
    phi_vtp[0] = phi_cp_cfd[0]
    phi_vtp[n_v - 1] = phi_cp_cfd[n_s - 1]
    for k in range(1, n_v - 1):
        phi_vtp[k] = nw_interp1d(r_cp, phi_cp_cfd, r_vtp[k])
    phi_star = np.zeros((n_v, n_s))
    for k in range(n_v):
        l1 = b.cp_loop_ind[k, 0]                # CPs l1 ... l2-1 are felt by vortex k
        l2 = b.cp_loop_ind[k, 1] + 1
        h_r = (r_vtp[k] - r_cp[l1:l2]) / r_vtp[k]
        a = b.a[k, l1:l2, :].T                  # Rows a[0] ... a[3]: a_hr,1 ... a_hr,4
        phi_star[k, l1:l2] = nw_phi_components(h_r, a, phi_vtp[k], b.phi[k, l1:l2],
                                               b.phi_s[k, l1:l2])
        delta = -(cyl_vtp[k, 2] - cyl_cp[l1:l2, 2])     # Normal to the disc, +ve upstream
        b.smear_fac[k, l1:l2] = smearing_factor(r_vtp[k], h_r, delta, phi_vtp[k], eps_cp)

    # 4. Sub-iterate until the correction velocities converge
    sub_iter = 0
    converged = False
    while not converged:
        # Velocities and angles with the current correction
        vstn = nw2vstn(vt, vn, cone, pitch, sweep, fbend)
        un = vn_rot + vstn_cfd[:, 2] + vstn[:, 2]      # Normal velocity
        ut = vt_rot + vstn_cfd[:, 1] + vstn[:, 1]      # Tangential velocity (< 0)
        phi_cp = np.arctan2(un, -ut)
        alpha_cp = phi_cp * R2D - twist
        vrel_cp = np.sqrt(un**2 + ut**2)
        vrel_vtp = np.zeros(n_v)
        vrel_vtp[0] = vrel_cp[0]
        vrel_vtp[n_v - 1] = vrel_cp[n_s - 1]
        for k in range(1, n_v - 1):
            vrel_vtp[k] = nw_interp1d(r_cp, vrel_cp, r_vtp[k])
        Gamma_steady, Gamma, dGamma, dyn_comp = nw_circulation(
            dt, vrel_cp, alpha_cp, b.Gamma_history, b.dyn_comp_history, alpha_2d, cl_2d, chord)
        dbeta = (vrel_vtp * dt) / r_vtp               # Advance in beta*, 2017b Eq. (3)

        # Induction of the vorticity shed this time step, smeared part only
        W_iter = np.zeros(n_s)
        Wn_iter = np.zeros(n_s)
        Wt_iter = np.zeros(n_s)
        X = np.zeros((n_v, n_s))
        Y = np.zeros((n_v, n_s))
        for k in range(n_v):
            l1 = b.cp_loop_ind[k, 0]
            l2 = b.cp_loop_ind[k, 1] + 1
            X[k, l1:l2], Y[k, l1:l2] = nw_xy_components(
                phi_star[k, l1:l2], b.dx[k, l1:l2], b.dy[k, l1:l2], dbeta[k], dGamma[k],
                b.X_history[k, l1:l2], b.Y_history[k, l1:l2])
            W_k = b.smear_fac[k, l1:l2] * (X[k, l1:l2] + Y[k, l1:l2])
            W_iter[l1:l2] += W_k
            # Normal and tangential parts with each vortex's helix angle (MF2019 Eq. 19)
            Wn_iter[l1:l2] += W_k * np.cos(phi_vtp[k])
            Wt_iter[l1:l2] += W_k * np.sin(phi_vtp[k])

        if recompute_relax:
            relax = nw_calc_relax(vrel_cp, alpha_cp, un, -ut, chord, dt, b.dx, b.dy, dbeta, b.phi,
                                  relax_safety)
            recompute_relax = False

        # Relax and measure the change
        W_iter = relax * W + (1 - relax) * W_iter
        section_res = np.sqrt((W_iter - W)**2)
        phi_cp_cyl = np.arctan2(vcyl_cfd[:, 2] + vn, -(v_rot + vcyl_cfd[:, 1] + vt))
        # Tip vortex within 45 deg of the rotor plane: axial residual, else tangential
        if phi_vtp[n_v - 1] < np.pi / 4:
            resmax = np.max(section_res * np.cos(phi_cp_cyl))
        else:
            resmax = np.max(section_res * np.sin(phi_cp_cyl))

        sub_iter = sub_iter + 1
        converged = resmax < state.nw_max_resid or sub_iter == 500
        if sub_iter == 100:
            recompute_relax = True
        if sub_iter == 500:
            if relax_safety < 0.8:
                relax_safety = relax_safety + 0.1
                relax = min(relax + 0.1 * (1 - relax), 0.99)
                recompute_relax = True
            # No convergence: keep the near-wake induction of the last time step
            W = b.W_history.copy()
            vn = b.vn_history.copy()
            vt = b.vt_history.copy()
            dyn_comp = b.dyn_comp_history.copy()
            Gamma_steady = b.Gamma_history.copy()
            X = b.X_history.copy()
            Y = b.Y_history.copy()
        else:
            W = W_iter
            vn = relax * vn + (1 - relax) * Wn_iter
            vt = relax * vt + (1 - relax) * Wt_iter

    # 5. Blade forces with the corrected inflow
    vstn = nw2vstn(vt, vn, cone, pitch, sweep, fbend)
    un = vn_rot + vstn_cfd[:, 2] + vstn[:, 2]
    ut = vt_rot + vstn_cfd[:, 1] + vstn[:, 1]
    phi_cp = np.arctan2(un, -ut)
    alpha_cp = phi_cp * R2D - twist
    alpha_cp_cfd = phi_cp_cfd * R2D - twist
    vrel_cp = np.sqrt(un**2 + ut**2)
    # Ang Li correction for coned or bent blades (Li et al., 2022): lift and drag
    # at the 3/4-chord AoA, force direction at the 1/4-chord inflow angle
    if l1point:
        thetadot = -omega * scone * R2D
        alpha_cp_34 = alpha_cp + (0.75 - state.x_aoa) * chord / vrel_cp * thetadot
        alpha_cp_14 = alpha_cp + (0.25 - state.x_aoa) * chord / vrel_cp * thetadot
    else:
        alpha_cp_34 = alpha_cp
        alpha_cp_14 = alpha_cp
    phi_cp_14 = (alpha_cp_14 + twist) * D2R
    cl = np.zeros(n_s)
    cd = np.zeros(n_s)
    for l in range(n_s):
        cl[l] = nw_interp1d(alpha_2d[:, l], cl_2d[:, l], alpha_cp_34[l])
        cd[l] = nw_interp1d(alpha_2d[:, l], cd_2d[:, l], alpha_cp_34[l])
    lift = 0.5 * density * vrel_cp**2 * chord * cl
    drag = 0.5 * density * vrel_cp**2 * chord * cd
    Fstn = np.zeros((n_s, 3))
    Fstn[:, 1] = lift * np.sin(phi_cp_14) - drag * np.cos(phi_cp_14)
    Fstn[:, 2] = lift * np.cos(phi_cp_14) + drag * np.sin(phi_cp_14)

    # 6. Diagnostics line and history for the next time step
    if b.fid is not None:
        b.fid.write("%15.7E%5d%15.7E" % (time, sub_iter, resmax))
        for value in np.concatenate((alpha_cp_cfd, alpha_cp, phi_vtp, vt, vn, dGamma)):
            b.fid.write("%15.7E" % value)
        b.fid.write("%15.7E\n" % (process_time() - cpu_t1))
    b.W_history = W
    b.dyn_comp_history = dyn_comp
    b.Gamma_history = Gamma_steady
    b.vn_history = vn
    b.vt_history = vt
    b.X_history = X
    b.Y_history = Y
    return Fstn


def smearing_factor(r, h_r, delta, phi, epsilon):
    """Smearing factor of a vortex at its trailing point (beta = 0).

    Meyer Forsting et al. (2019) Eq. (21)-(22) with beta = 0, as in Meyer
    Forsting et al. (2020) Sect. 2.3. delta: distance normal to the rotor disc
    (+ve upstream), phi: helix angle [rad]. Works on scalars and arrays.
    """
    dist_orth_sq = (r * h_r)**2 + (delta * np.cos(phi))**2
    return np.exp(-dist_orth_sq / epsilon**2)


def nw_hr_components(h, h_r):
    """Constants of one vortex-section pair, h = r_vtp - r_cp and h_r = h/r_vtp.

    Returns phi (root-corrected Phi_C), phi_s (straight vortex), the amplitudes
    dx and dy (without Gamma and Phi) and the fit coefficients a (4 values).
    Pirrung et al. (2016) and (2017b).
    """
    # Phi of Wang and Coton, 2016 Eq. (4)
    if h_r > 0:
        phi = np.pi / 4 * abs((1 + 0.5 * h_r) * np.log(1 - h_r))
    else:
        phi = np.log(1 - h_r) / (1.5 + np.log(1 - 0.5 * h_r))

    # Amplitudes, 2016 Eq. (25). Important: D_Y is negative.
    dw = 1 / (4 * np.pi * abs(h) * h_r)
    dx = 1.359 * dw
    dy = -0.359 / 4 * dw

    # Root correction, 2016 Eq. (27)-(29). Important: applied once, through Phi (Phi_C, Eq. 31).
    root_corr = ((dx * (1 - np.exp(-np.pi / (2 * phi))) + dy * (1 - np.exp(-2 * np.pi / phi)))
                 / (dx + dy))
    phi = root_corr * phi

    # Phi of a straight vortex, 2017b Eq. (7)
    phi_s = 0.788 * abs(h_r)

    # Fit coefficients of the mixing factor, 2017b Eq. (10) and (12)
    N, P = NW_N, NW_P
    if h_r < 0:
        a = (N[:, 0] + N[:, 1] * np.exp(N[:, 2] * h_r) + N[:, 3] * np.exp(N[:, 4] * h_r)
             - N[:, 3] - N[:, 1])
    else:
        a = P[:, 0] + P[:, 1] * h_r + P[:, 2] * h_r**2 + P[:, 3] * h_r**3
    return phi, phi_s, dx, dy, a


def nw_phi_components(h_r, a, phi, phi_c, phi_s):
    """Phi* for the helix angle phi [rad]: weighted mix of straight and circular Phi.

    Pirrung et al. (2017b) Eq. (8)-(11), with the root-corrected Phi_C of the
    corrigendum. a[0:4] are the fit coefficients. Works on scalars and arrays.
    """
    k_neg = (a[0] + a[1] * np.exp(a[2] * (np.pi / 2 - phi)) + a[3] * np.exp(-8 * (np.pi / 2 - phi))
             - a[1] - a[3])
    k_pos = a[0] + a[1] * phi + a[2] * phi**2 + a[3] * phi**3
    k_phi = np.where(h_r < 0, k_neg, k_pos)     # Eq. (11) for h/r < 0, Eq. (9) otherwise
    # The fit is only valid below 89 degrees: keep the weight in [0, 1]
    k_phi = np.clip(k_phi, 0, 1)
    return k_phi * phi_s + (1 - k_phi) * phi_c


def nw_xy_components(phi_star, dx, dy, dbeta, dGamma, X_history, Y_history):
    """Advance the slow (X) and fast (Y) indicial functions by one time step.

    Pirrung et al. (2016) Eq. (25)-(26) with Phi* of (2017b): new element plus
    the decayed history. Works on scalars and arrays.
    """
    X_exp = np.exp(-dbeta / phi_star)
    Y_exp = np.exp(-4 * dbeta / phi_star)
    X = dGamma * dx * phi_star * (1 - X_exp) + X_history * X_exp
    Y = dGamma * dy * phi_star * (1 - Y_exp) + Y_history * Y_exp
    return X, Y


def nw_circulation(dt, vrel, alpha, Gamma_history, dyn_comp_history, alpha_2d, cl_2d, chord):
    """Bound circulation, its time-filtered value and the trailed vortex strengths.

    alpha in degrees. Returns Gamma_steady (n_s), Gamma (n_s), dGamma (n_s + 1)
    and the filter components dyn_comp (n_s, 3).
    """
    # Quasi-steady circulation from the static polar
    n_s = len(vrel)
    cl = np.zeros(n_s)
    for l in range(n_s):
        cl[l] = nw_interp1d(alpha_2d[:, l], cl_2d[:, l], alpha[l])
    Gamma_steady = 0.5 * vrel * chord * cl

    # Time filter: Madsen and Gaunaa (2004), Pirrung et al. (2017a) Eq. (15)-(19).
    # The constants are A_j/2 and b_j; A3 = 0.2625 (paper: 0.2656) so that the A_j sum to 1.
    t = dt / (chord / (2 * vrel))
    dyn_comp = np.zeros((n_s, 3))
    dyn_comp[:, 0] = (dyn_comp_history[:, 0] * np.exp(-0.3064 * t)
                      + 0.27735 * (Gamma_steady + Gamma_history) * (1 - np.exp(-0.3064 * t)))
    dyn_comp[:, 1] = (dyn_comp_history[:, 1] * np.exp(-0.0439 * t)
                      + 0.0914 * (Gamma_steady + Gamma_history) * (1 - np.exp(-0.0439 * t)))
    dyn_comp[:, 2] = (dyn_comp_history[:, 2] * np.exp(-3.227 * t)
                      + 0.13125 * (Gamma_steady + Gamma_history) * (1 - np.exp(-3.227 * t)))
    Gamma = dyn_comp[:, 0] + dyn_comp[:, 1] + dyn_comp[:, 2]

    # Trailed strengths: the end vortices take all Gamma of the end sections
    dGamma = np.zeros(n_s + 1)
    dGamma[0] = Gamma[0]
    dGamma[1:n_s] = Gamma[1:n_s] - Gamma[0:n_s - 1]
    dGamma[n_s] = -Gamma[n_s - 1]
    return Gamma_steady, Gamma, dGamma, dyn_comp


def nw_calc_relax(vrel, alpha, vn, vt, chord, dt, dx, dy, dbeta, phi, safety):
    """Relaxation factor of the sub-iteration, Pirrung et al. (2017a) Eq. (27)-(31).

    Uses the corrected Eq. (28). alpha in degrees; vn and vt are the normal and
    tangential speeds (vt > 0); phi is the root-corrected Phi_C.
    """
    relax = 0
    for l in range(len(vrel)):
        # Gamma_dyn/Gamma_QS: Eq. (15) averaged over the time step
        tau = (2 * vrel[l]) / chord[l]
        d = (1 - 0.5547 / (tau * dt * 0.3064) * (1 - np.exp(-tau * dt * 0.3064))
             - 0.1828 / (tau * dt * 0.0439) * (1 - np.exp(-tau * dt * 0.0439))
             - 0.2625 / (tau * dt * 3.227) * (1 - np.exp(-tau * dt * 3.227)))
        # B_v of the inner (l) and outer (l + 1) vortex of section l, Eq. (28)
        A1 = (dx[l, l] * phi[l, l] * (1 - np.exp(-dbeta[l] / phi[l, l]))
              + dy[l, l] * phi[l, l] * (1 - np.exp(-4 * dbeta[l] / phi[l, l])))
        A2 = (dx[l + 1, l] * phi[l + 1, l] * (1 - np.exp(-dbeta[l + 1] / phi[l + 1, l]))
              + dy[l + 1, l] * phi[l + 1, l] * (1 - np.exp(-4 * dbeta[l + 1] / phi[l + 1, l])))
        # Gradient, Eq. (27): alpha converted to radians
        gradw = (d * np.pi * chord[l] * (A1 - A2)
                 * ((alpha[l] * np.pi / 180 * vn[l]) / vrel[l]
                    + vrel[l] / (abs(vt[l]) * ((vn[l] / abs(vt[l]))**2 + 1))))
        relax = max(relax, -(1 + gradw) / (1 - gradw))      # Eq. (31)
    relax = relax + safety * (1 - relax)
    return min(relax, 0.99)


def nw_interp1d(x, y, xin):
    """Linear interpolation in the ascending table (x, y), linear extrapolation outside."""
    if xin < x[0]:
        return y[0] + (xin - x[0]) / (x[1] - x[0]) * (y[1] - y[0])
    if xin > x[-1]:
        return y[-1] + (xin - x[-1]) / (x[-1] - x[-2]) * (y[-1] - y[-2])
    return np.interp(xin, x, y)


def nw2vstn(vt, vn, cone, pitch, sweep, fbend):
    """Correction velocities (cylindrical vt, vn) in section axes, shape (n, 3)."""
    v = np.column_stack((np.zeros(len(vt)), vt, vn))
    return vxyz_vstn(v, 0, 0, 0, cone, pitch, sweep, fbend)


def vxyz_vstn(v, yaw, tilt, azimut, cone, pitch, sweep, fbend):
    """Rotate vectors (rows of v) from the Cartesian or cylindrical frame to section axes.

    Order: yaw, tilt, azimuth, cone, pitch, sweep, bend. Angles in degrees;
    sweep and fbend can be given per row.
    """
    v = rot_y(v, yaw)
    v = rot_x(v, tilt)
    v = rot_z(v, -azimut)
    v = rot_y(v, cone)
    v = rot_x(v, -pitch)
    v = rot_z(v, -sweep)
    return rot_y(v, fbend)


def rot_x(v, angle):
    """Rotation about the first axis (tilt, pitch)."""
    c = np.cos(angle * np.pi / 180)
    s = np.sin(angle * np.pi / 180)
    return np.column_stack((v[:, 0], c * v[:, 1] + s * v[:, 2], -s * v[:, 1] + c * v[:, 2]))


def rot_y(v, angle):
    """Rotation about the second axis (yaw, cone, bend)."""
    c = np.cos(angle * np.pi / 180)
    s = np.sin(angle * np.pi / 180)
    return np.column_stack((c * v[:, 0] - s * v[:, 2], v[:, 1], s * v[:, 0] + c * v[:, 2]))


def rot_z(v, angle):
    """Rotation about the third axis (azimuth, sweep)."""
    c = np.cos(angle * np.pi / 180)
    s = np.sin(angle * np.pi / 180)
    return np.column_stack((c * v[:, 0] - s * v[:, 1], s * v[:, 0] + c * v[:, 1], v[:, 2]))
