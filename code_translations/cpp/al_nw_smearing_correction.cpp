// Fast smearing correction for the actuator line (C++ version).
// See al_nw_smearing_correction.hpp for the usage.
#include "al_nw_smearing_correction.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>

static const double pi = 4.0 * atan(1.0);
static const double d2r = pi / 180.0;
static const double r2d = 180.0 / pi;

// Fit constants of the mixing factor k_Phi, Pirrung et al. (2017b). Row i gives a_hr,i.
// P: h/r > 0, Eq. (14). N: h/r < 0, Eq. (13).
static const double nw_P[4][5] = {
    {-1.64636754184988, 8.14821474772595, -12.1784861715161, 5.02652654794305, 21.7713111885100},
    {-0.499006672527007, 6.08465268675599, -15.1712005446736, 14.8254135311781, -2.42318681495122},
    {3.90835606394205, -18.7662330478906, 39.1243282362892, -29.4870086106427, -60.2947309276912},
    {-1.60622843516794, 7.42952711674470, -15.8594777074799, 11.6870163960092, -195.060865844227}};
static const double nw_N[4][5] = {
    {1.01933206359188, -0.135668151345576, 0.395524705910666, 0.0801822997046969, 44.8347503471217},
    {12.9874512870379, 49.9999999999916, 0.00235345218088083, 11.3116074964299, 3935.34323353354},
    {-0.690159136848837, 101.238775997879, -0.00154467589177718, 3.99520376143011, 0.394541495293350},
    {-0.269253409815015, 49.9999999999843, -0.00247975579295973, 0.403642518874381, 1.16610462658971}};

// Rotation about the first axis (tilt, pitch), angle in degrees
static Vec3 rot_x(const Vec3& a, double angle) {
    double c = cos(angle * pi / 180.0);
    double s = sin(angle * pi / 180.0);
    Vec3 v = {a[0], c * a[1] + s * a[2], -s * a[1] + c * a[2]};
    return v;
}

// Rotation about the second axis (yaw, cone, bend), angle in degrees
static Vec3 rot_y(const Vec3& a, double angle) {
    double c = cos(angle * pi / 180.0);
    double s = sin(angle * pi / 180.0);
    Vec3 v = {c * a[0] - s * a[2], a[1], s * a[0] + c * a[2]};
    return v;
}

// Rotation about the third axis (azimuth, sweep), angle in degrees
static Vec3 rot_z(const Vec3& a, double angle) {
    double c = cos(angle * pi / 180.0);
    double s = sin(angle * pi / 180.0);
    Vec3 v = {c * a[0] - s * a[1], s * a[0] + c * a[1], a[2]};
    return v;
}

NWState init_nw_smearing_correction(int n_rotors, int n_blades, int n_sections,
                                    const std::string& project) {
    NWState state;
    state.kaya = true;
    state.nw_max_resid = 1e-6;
    state.x_aoa = 0.5;
    state.n_s = n_sections;
    state.n_v = n_sections + 1;
    state.n_blades = n_blades;
    int n_s = state.n_s;
    int n_v = state.n_v;
    for (int i = 1; i <= n_rotors; i++) {
        for (int j = 1; j <= n_blades; j++) {
            BladeState b;
            b.refresh = true;
            b.cp_first.assign(n_v, 0);
            b.cp_last.assign(n_v, n_s - 1);
            b.W_history.assign(n_s, 0.0);
            b.vn_history.assign(n_s, 0.0);
            b.vt_history.assign(n_s, 0.0);
            b.Gamma_history.assign(n_s, 0.0);
            b.dyn_comp_history.assign(n_s, Vec3{0.0, 0.0, 0.0});
            b.X_history.assign(n_v, Vec(n_s, 0.0));
            b.Y_history.assign(n_v, Vec(n_s, 0.0));
            b.phi.assign(n_v, Vec(n_s, 0.0));
            b.phi_s.assign(n_v, Vec(n_s, 0.0));
            b.dx.assign(n_v, Vec(n_s, 0.0));
            b.dy.assign(n_v, Vec(n_s, 0.0));
            b.a.assign(n_v, std::vector<Coef4>(n_s));
            b.smear_fac.assign(n_v, Vec(n_s, 0.0));
            b.fid = NULL;
            if (!project.empty()) {
                char name[256];
                snprintf(name, sizeof(name), "%s.nw%02d.bl%02d", project.c_str(), i, j);
                b.fid = fopen(name, "w");
            }
            state.blade.push_back(b);
        }
    }
    return state;
}

void nw_smearing_correction(NWState& state, int i, int j, double time, double dt,
                            double density, const std::vector<Vec3>& cyl_cp,
                            const std::vector<Vec3>& cyl_vtp, const std::vector<Vec3>& vstn_cfd,
                            const std::vector<Vec3>& vcyl_cfd, double eps_cp, double omega,
                            const Vec& chord, const Vec& twist, double pitch,
                            const Grid& cl_2d, const Grid& cd_2d, const Grid& alpha_2d,
                            double cone, const Vec& fbend, const Vec& sweep, const Vec& dsweep,
                            bool l1point, std::vector<Vec3>& Fstn) {
    if (!state.kaya) return;
    int n_s = state.n_s;
    int n_v = state.n_v;
    Vec r_cp(n_s), r_vtp(n_v);
    for (int l = 0; l < n_s; l++) r_cp[l] = cyl_cp[l][0];
    for (int k = 0; k < n_v; k++) r_vtp[k] = cyl_vtp[k][0];
    // Important: r = 0 at the first trailing point is a singularity
    if (r_vtp[0] < 0.1) {
        printf("!!!! WARNING smearing correction switched off !!!!\n");
        printf("First vortex trailing point is zero, please modify\n");
        state.kaya = false;
        return;
    }
    clock_t cpu_t1 = clock();
    BladeState& b = state.blade[(i - 1) * state.n_blades + (j - 1)];

    // 1. Preamble: start from the previous time step
    double relax_safety = 0.2;
    double relax = 0.0;
    bool recompute_relax = true;
    Vec W = b.W_history;
    Vec vn = b.vn_history;
    Vec vt = b.vt_history;
    Vec scone(n_s), v_rot(n_s), vn_rot(n_s), vt_rot(n_s);
    for (int l = 0; l < n_s; l++) {
        scone[l] = sin((cone + fbend[l]) * d2r);
        v_rot[l] = -(omega * r_cp[l]);  // Rotational speed in section axes
        vn_rot[l] = -v_rot[l] * sin(-pitch * d2r);
        vt_rot[l] = v_rot[l] * cos(-pitch * d2r) * cos(dsweep[l] * d2r);
    }

    // 2. Geometry constants: first call, or after a refresh
    if (b.refresh) {
        double core_cutoff = 1.82 * eps_cp;  // 99th percentile of the Gaussian
        for (int k = 0; k < n_v; k++) {
            b.cp_first[k] = 0;
            b.cp_last[k] = n_s - 1;
            for (int l = 0; l < n_s; l++) {
                double h = r_vtp[k] - r_cp[l];
                double h_r = h / r_vtp[k];
                nw_hr_components(h, h_r, b.phi[k][l], b.phi_s[k][l], b.dx[k][l], b.dy[k][l],
                                 b.a[k][l]);
                b.smear_fac[k][l] = smearing_factor(r_vtp[k], h_r, 0.0, 0.0, eps_cp);
                // CPs felt by vortex k: inboard and outboard ones within the cut-off
                if (h > 0.0) {
                    if (fabs(h) > core_cutoff) b.cp_first[k] = l + 1;
                } else {
                    if (fabs(h) < core_cutoff) b.cp_last[k] = l;
                }
            }
        }
        b.refresh = false;
    }

    // 3. Helix angles from the CFD only (MF2019 step 2), fixed during the sub-iteration
    Vec phi_cp_cfd(n_s), phi_vtp(n_v);
    for (int l = 0; l < n_s; l++) {
        phi_cp_cfd[l] = atan2(vcyl_cfd[l][2], -(v_rot[l] + vcyl_cfd[l][1]));
        if (phi_cp_cfd[l] > pi / 2.0) phi_cp_cfd[l] = pi - phi_cp_cfd[l];
    }
    phi_vtp[0] = phi_cp_cfd[0];
    phi_vtp[n_v - 1] = phi_cp_cfd[n_s - 1];
    for (int k = 1; k < n_v - 1; k++) phi_vtp[k] = nw_interp1d(r_cp, phi_cp_cfd, r_vtp[k]);
    Grid phi_star(n_v, Vec(n_s, 0.0));
    for (int k = 0; k < n_v; k++) {
        for (int l = b.cp_first[k]; l <= b.cp_last[k]; l++) {
            double h = r_vtp[k] - r_cp[l];
            double h_r = h / r_vtp[k];
            phi_star[k][l] = nw_phi_components(h_r, b.a[k][l], phi_vtp[k], b.phi[k][l],
                                               b.phi_s[k][l]);
            double delta = -(cyl_vtp[k][2] - cyl_cp[l][2]);  // Normal to the disc, +ve upstream
            b.smear_fac[k][l] = smearing_factor(r_vtp[k], h_r, delta, phi_vtp[k], eps_cp);
        }
    }

    // 4. Sub-iterate until the correction velocities converge
    int sub_iter = 0;
    bool converged = false;
    double resmax = 0.0;
    Vec un(n_s), ut(n_s), phi_cp(n_s), alpha_cp(n_s), vrel_cp(n_s), vt_pos(n_s);
    Vec vrel_vtp(n_v), dbeta(n_v), Gamma_steady, Gamma, dGamma;
    std::vector<Vec3> vstn, dyn_comp;
    Grid X(n_v, Vec(n_s, 0.0)), Y(n_v, Vec(n_s, 0.0));
    while (!converged) {
        // Velocities and angles with the current correction
        vstn = nw2vstn(vt, vn, cone, pitch, sweep, fbend);
        for (int l = 0; l < n_s; l++) {
            un[l] = vn_rot[l] + vstn_cfd[l][2] + vstn[l][2];  // Normal velocity
            ut[l] = vt_rot[l] + vstn_cfd[l][1] + vstn[l][1];  // Tangential velocity (< 0)
            phi_cp[l] = atan2(un[l], -ut[l]);
            alpha_cp[l] = phi_cp[l] * r2d - twist[l];
            vrel_cp[l] = sqrt(un[l] * un[l] + ut[l] * ut[l]);
        }
        vrel_vtp[0] = vrel_cp[0];
        vrel_vtp[n_v - 1] = vrel_cp[n_s - 1];
        for (int k = 1; k < n_v - 1; k++) vrel_vtp[k] = nw_interp1d(r_cp, vrel_cp, r_vtp[k]);
        nw_circulation(dt, vrel_cp, alpha_cp, b.Gamma_history, b.dyn_comp_history, alpha_2d, cl_2d,
                       chord, Gamma_steady, Gamma, dGamma, dyn_comp);
        for (int k = 0; k < n_v; k++) dbeta[k] = (vrel_vtp[k] * dt) / r_vtp[k];  // 2017b Eq. (3)

        // Induction of the vorticity shed this time step, smeared part only
        Vec W_iter(n_s, 0.0), Wn_iter(n_s, 0.0), Wt_iter(n_s, 0.0);
        for (int k = 0; k < n_v; k++) {
            for (int l = 0; l < n_s; l++) {
                X[k][l] = 0.0;
                Y[k][l] = 0.0;
            }
            for (int l = b.cp_first[k]; l <= b.cp_last[k]; l++) {
                nw_xy_components(phi_star[k][l], b.dx[k][l], b.dy[k][l], dbeta[k], dGamma[k],
                                 b.X_history[k][l], b.Y_history[k][l], X[k][l], Y[k][l]);
                double W_k = b.smear_fac[k][l] * (X[k][l] + Y[k][l]);
                W_iter[l] = W_iter[l] + W_k;
                // Normal and tangential parts with each vortex's helix angle (MF2019 Eq. 19)
                Wn_iter[l] = Wn_iter[l] + W_k * cos(phi_vtp[k]);
                Wt_iter[l] = Wt_iter[l] + W_k * sin(phi_vtp[k]);
            }
        }

        if (recompute_relax) {
            for (int l = 0; l < n_s; l++) vt_pos[l] = -ut[l];
            relax = nw_calc_relax(vrel_cp, alpha_cp, un, vt_pos, chord, dt, b.dx, b.dy, dbeta,
                                  b.phi, relax_safety);
            recompute_relax = false;
        }

        // Relax and measure the change
        for (int l = 0; l < n_s; l++) {
            W_iter[l] = relax * W[l] + (1.0 - relax) * W_iter[l];
            double section_res = sqrt((W_iter[l] - W[l]) * (W_iter[l] - W[l]));
            double phi_cp_cyl = atan2(vcyl_cfd[l][2] + vn[l], -(v_rot[l] + vcyl_cfd[l][1] + vt[l]));
            // Tip vortex within 45 deg of the rotor plane: axial residual, else tangential
            double res;
            if (phi_vtp[n_v - 1] < pi / 4.0)
                res = section_res * cos(phi_cp_cyl);
            else
                res = section_res * sin(phi_cp_cyl);
            if (l == 0 || res > resmax) resmax = res;
        }

        sub_iter = sub_iter + 1;
        converged = resmax < state.nw_max_resid || sub_iter == 500;
        if (sub_iter == 100) recompute_relax = true;
        if (sub_iter == 500) {
            if (relax_safety < 0.8) {
                relax_safety = relax_safety + 0.1;
                relax = std::min(relax + 0.1 * (1.0 - relax), 0.99);
                recompute_relax = true;
            }
            // No convergence: keep the near-wake induction of the last time step
            W = b.W_history;
            vn = b.vn_history;
            vt = b.vt_history;
            dyn_comp = b.dyn_comp_history;
            Gamma_steady = b.Gamma_history;
            X = b.X_history;
            Y = b.Y_history;
        } else {
            W = W_iter;
            for (int l = 0; l < n_s; l++) {
                vn[l] = relax * vn[l] + (1.0 - relax) * Wn_iter[l];
                vt[l] = relax * vt[l] + (1.0 - relax) * Wt_iter[l];
            }
        }
    }

    // 5. Blade forces with the corrected inflow
    vstn = nw2vstn(vt, vn, cone, pitch, sweep, fbend);
    Vec alpha_cp_cfd(n_s);
    Fstn.assign(n_s, Vec3{0.0, 0.0, 0.0});
    for (int l = 0; l < n_s; l++) {
        un[l] = vn_rot[l] + vstn_cfd[l][2] + vstn[l][2];
        ut[l] = vt_rot[l] + vstn_cfd[l][1] + vstn[l][1];
        phi_cp[l] = atan2(un[l], -ut[l]);
        alpha_cp[l] = phi_cp[l] * r2d - twist[l];
        alpha_cp_cfd[l] = phi_cp_cfd[l] * r2d - twist[l];
        vrel_cp[l] = sqrt(un[l] * un[l] + ut[l] * ut[l]);
        // Ang Li correction for coned or bent blades (Li et al., 2022): lift and drag
        // at the 3/4-chord AoA, force direction at the 1/4-chord inflow angle
        double alpha_cp_34 = alpha_cp[l];
        double alpha_cp_14 = alpha_cp[l];
        if (l1point) {
            double thetadot = -omega * scone[l] * r2d;
            alpha_cp_34 = alpha_cp[l] + (0.75 - state.x_aoa) * chord[l] / vrel_cp[l] * thetadot;
            alpha_cp_14 = alpha_cp[l] + (0.25 - state.x_aoa) * chord[l] / vrel_cp[l] * thetadot;
        }
        double phi_cp_14 = (alpha_cp_14 + twist[l]) * d2r;
        double cl = nw_interp1d(alpha_2d[l], cl_2d[l], alpha_cp_34);
        double cd = nw_interp1d(alpha_2d[l], cd_2d[l], alpha_cp_34);
        double lift = 0.5 * density * (vrel_cp[l] * vrel_cp[l]) * chord[l] * cl;
        double drag = 0.5 * density * (vrel_cp[l] * vrel_cp[l]) * chord[l] * cd;
        Fstn[l][1] = lift * sin(phi_cp_14) - drag * cos(phi_cp_14);
        Fstn[l][2] = lift * cos(phi_cp_14) + drag * sin(phi_cp_14);
    }

    // 6. Diagnostics line and history for the next time step
    if (b.fid != NULL) {
        fprintf(b.fid, "%15.7E%5d%15.7E", time, sub_iter, resmax);
        for (int l = 0; l < n_s; l++) fprintf(b.fid, "%15.7E", alpha_cp_cfd[l]);
        for (int l = 0; l < n_s; l++) fprintf(b.fid, "%15.7E", alpha_cp[l]);
        for (int k = 0; k < n_v; k++) fprintf(b.fid, "%15.7E", phi_vtp[k]);
        for (int l = 0; l < n_s; l++) fprintf(b.fid, "%15.7E", vt[l]);
        for (int l = 0; l < n_s; l++) fprintf(b.fid, "%15.7E", vn[l]);
        for (int k = 0; k < n_v; k++) fprintf(b.fid, "%15.7E", dGamma[k]);
        fprintf(b.fid, "%15.7E\n", double(clock() - cpu_t1) / CLOCKS_PER_SEC);
    }
    b.W_history = W;
    b.dyn_comp_history = dyn_comp;
    b.Gamma_history = Gamma_steady;
    b.vn_history = vn;
    b.vt_history = vt;
    b.X_history = X;
    b.Y_history = Y;
}

// Smearing factor of a vortex at its trailing point (beta = 0): Meyer Forsting et al.
// (2019) Eq. (21)-(22) with beta = 0, as in Meyer Forsting et al. (2020) Sect. 2.3.
// delta: distance normal to the rotor disc (+ve upstream), phi: helix angle [rad].
double smearing_factor(double r, double h_r, double delta, double phi, double epsilon) {
    double rh = r * h_r;
    double dc = delta * cos(phi);
    return exp(-(rh * rh + dc * dc) / (epsilon * epsilon));
}

// Constants of one vortex-section pair, h = r_vtp - r_cp and h_r = h/r_vtp:
// phi (root-corrected Phi_C), phi_s (straight vortex), the amplitudes dx and dy
// (without Gamma and Phi) and the fit coefficients a. Pirrung et al. (2016), (2017b).
void nw_hr_components(double h, double h_r, double& phi, double& phi_s, double& dx, double& dy,
                      Coef4& a) {
    // Phi of Wang and Coton, 2016 Eq. (4)
    if (h_r > 0.0)
        phi = pi / 4.0 * fabs((1.0 + 0.5 * h_r) * log(1.0 - h_r));
    else
        phi = log(1.0 - h_r) / (1.5 + log(1.0 - 0.5 * h_r));

    // Amplitudes, 2016 Eq. (25). Important: D_Y is negative.
    double dw = 1.0 / (4.0 * pi * fabs(h) * h_r);
    dx = 1.359 * dw;
    dy = -0.359 / 4.0 * dw;

    // Root correction, 2016 Eq. (27)-(29). Important: applied once, through Phi (Phi_C, Eq. 31).
    double root_corr = (dx * (1.0 - exp(-pi / (2.0 * phi))) + dy * (1.0 - exp(-2.0 * pi / phi)))
                       / (dx + dy);
    phi = root_corr * phi;

    // Phi of a straight vortex, 2017b Eq. (7)
    phi_s = 0.788 * fabs(h_r);

    // Fit coefficients of the mixing factor, 2017b Eq. (10) and (12)
    for (int n = 0; n < 4; n++) {
        if (h_r < 0.0)
            a[n] = nw_N[n][0] + nw_N[n][1] * exp(nw_N[n][2] * h_r)
                   + nw_N[n][3] * exp(nw_N[n][4] * h_r) - nw_N[n][3] - nw_N[n][1];
        else
            a[n] = nw_P[n][0] + nw_P[n][1] * h_r + nw_P[n][2] * (h_r * h_r)
                   + nw_P[n][3] * (h_r * h_r * h_r);
    }
}

// Phi* for the helix angle phi [rad]: weighted mix of straight and circular Phi.
// Pirrung et al. (2017b) Eq. (8)-(11), with the root-corrected Phi_C of the corrigendum.
double nw_phi_components(double h_r, const Coef4& a, double phi, double phi_c, double phi_s) {
    double k_phi;
    if (h_r < 0.0)
        k_phi = a[0] + a[1] * exp(a[2] * (pi / 2.0 - phi)) + a[3] * exp(-8.0 * (pi / 2.0 - phi))
                - a[1] - a[3];
    else
        k_phi = a[0] + a[1] * phi + a[2] * (phi * phi) + a[3] * (phi * phi * phi);
    // The fit is only valid below 89 degrees: keep the weight in [0, 1]
    k_phi = std::max(std::min(k_phi, 1.0), 0.0);
    return k_phi * phi_s + (1.0 - k_phi) * phi_c;
}

// Advance the slow (X) and fast (Y) indicial functions by one time step: new element
// plus the decayed history. Pirrung et al. (2016) Eq. (25)-(26) with Phi* of (2017b).
void nw_xy_components(double phi_star, double dx, double dy, double dbeta, double dGamma,
                      double X_history, double Y_history, double& X, double& Y) {
    double X_exp = exp(-dbeta / phi_star);
    double Y_exp = exp(-4.0 * dbeta / phi_star);
    X = dGamma * dx * phi_star * (1.0 - X_exp) + X_history * X_exp;
    Y = dGamma * dy * phi_star * (1.0 - Y_exp) + Y_history * Y_exp;
}

// Bound circulation, its time-filtered value and the trailed strengths (n_s + 1)
void nw_circulation(double dt, const Vec& vrel, const Vec& alpha, const Vec& Gamma_history,
                    const std::vector<Vec3>& dyn_comp_history, const Grid& alpha_2d,
                    const Grid& cl_2d, const Vec& chord, Vec& Gamma_steady, Vec& Gamma,
                    Vec& dGamma, std::vector<Vec3>& dyn_comp) {
    int n = vrel.size();
    Gamma_steady.assign(n, 0.0);
    Gamma.assign(n, 0.0);
    dGamma.assign(n + 1, 0.0);
    dyn_comp.assign(n, Vec3{0.0, 0.0, 0.0});
    for (int l = 0; l < n; l++) {
        // Quasi-steady circulation from the static polar
        double cl = nw_interp1d(alpha_2d[l], cl_2d[l], alpha[l]);
        Gamma_steady[l] = 0.5 * vrel[l] * chord[l] * cl;
        // Time filter: Madsen and Gaunaa (2004), Pirrung et al. (2017a) Eq. (15)-(19).
        // The constants are A_j/2 and b_j; A3 = 0.2625 (paper: 0.2656) so that the A_j sum to 1.
        double t = dt / (chord[l] / (2.0 * vrel[l]));
        double G_sum = Gamma_steady[l] + Gamma_history[l];
        dyn_comp[l][0] = dyn_comp_history[l][0] * exp(-0.3064 * t)
                         + 0.27735 * G_sum * (1.0 - exp(-0.3064 * t));
        dyn_comp[l][1] = dyn_comp_history[l][1] * exp(-0.0439 * t)
                         + 0.0914 * G_sum * (1.0 - exp(-0.0439 * t));
        dyn_comp[l][2] = dyn_comp_history[l][2] * exp(-3.227 * t)
                         + 0.13125 * G_sum * (1.0 - exp(-3.227 * t));
        Gamma[l] = dyn_comp[l][0] + dyn_comp[l][1] + dyn_comp[l][2];
    }
    // Trailed strengths: the end vortices take all Gamma of the end sections
    dGamma[0] = Gamma[0];
    for (int l = 1; l < n; l++) dGamma[l] = Gamma[l] - Gamma[l - 1];
    dGamma[n] = -Gamma[n - 1];
}

// Relaxation factor of the sub-iteration, Pirrung et al. (2017a) Eq. (27)-(31) with the
// corrected Eq. (28). alpha in degrees; vn and vt are the normal and tangential speeds (vt > 0).
double nw_calc_relax(const Vec& vrel, const Vec& alpha, const Vec& vn, const Vec& vt,
                     const Vec& chord, double dt, const Grid& dx, const Grid& dy,
                     const Vec& dbeta, const Grid& phi, double safety) {
    double relax = 0.0;
    for (int l = 0; l < int(vrel.size()); l++) {
        // Gamma_dyn/Gamma_QS: Eq. (15) averaged over the time step
        double tau = (2.0 * vrel[l]) / chord[l];
        double d = 1 - 0.5547 / (tau * dt * 0.3064) * (1 - exp(-tau * dt * 0.3064))
                     - 0.1828 / (tau * dt * 0.0439) * (1 - exp(-tau * dt * 0.0439))
                     - 0.2625 / (tau * dt * 3.227) * (1 - exp(-tau * dt * 3.227));
        // B_v of the inner (l) and outer (l + 1) vortex of section l, Eq. (28)
        double A1 = dx[l][l] * phi[l][l] * (1.0 - exp(-dbeta[l] / phi[l][l]))
                    + dy[l][l] * phi[l][l] * (1.0 - exp(-4.0 * dbeta[l] / phi[l][l]));
        double A2 = dx[l+1][l] * phi[l+1][l] * (1.0 - exp(-dbeta[l+1] / phi[l+1][l]))
                    + dy[l+1][l] * phi[l+1][l] * (1.0 - exp(-4.0 * dbeta[l+1] / phi[l+1][l]));
        // Gradient, Eq. (27): alpha converted to radians
        double q = vn[l] / fabs(vt[l]);
        double gradw = d * pi * chord[l] * (A1 - A2)
                       * ((alpha[l] * pi / 180.0 * vn[l]) / vrel[l]
                          + vrel[l] / (fabs(vt[l]) * (q * q + 1.0)));
        relax = std::max(relax, -(1.0 + gradw) / (1.0 - gradw));  // Eq. (31)
    }
    relax = relax + safety * (1.0 - relax);
    return std::min(relax, 0.99);
}

// Linear interpolation in the ascending table (x, y), linear extrapolation outside
double nw_interp1d(const Vec& x, const Vec& y, double xin) {
    int n = x.size();
    if (xin < x[0])
        return y[0] + (xin - x[0]) / (x[1] - x[0]) * (y[1] - y[0]);
    if (xin > x[n - 1])
        return y[n - 1] + (xin - x[n - 1]) / (x[n - 1] - x[n - 2]) * (y[n - 1] - y[n - 2]);
    int i = 1;
    while (x[i] < xin) i++;
    return y[i - 1] + (xin - x[i - 1]) / (x[i] - x[i - 1]) * (y[i] - y[i - 1]);
}

// Rotate a vector from the Cartesian or cylindrical frame to section axes.
// Order: yaw, tilt, azimuth, cone, pitch, sweep, bend.
Vec3 vxyz_vstn(const Vec3& a, double yaw, double tilt, double azimut, double cone, double pitch,
               double sweep, double fbend) {
    Vec3 v = rot_y(a, yaw);
    v = rot_x(v, tilt);
    v = rot_z(v, -azimut);
    v = rot_y(v, cone);
    v = rot_x(v, -pitch);
    v = rot_z(v, -sweep);
    v = rot_y(v, fbend);
    return v;
}

// Correction velocities (cylindrical vt, vn) in section axes
std::vector<Vec3> nw2vstn(const Vec& vt, const Vec& vn, double cone, double pitch,
                          const Vec& sweep, const Vec& fbend) {
    std::vector<Vec3> vstn(vt.size());
    for (int l = 0; l < int(vt.size()); l++) {
        Vec3 v = {0.0, vt[l], vn[l]};
        vstn[l] = vxyz_vstn(v, 0.0, 0.0, 0.0, cone, pitch, sweep[l], fbend[l]);
    }
    return vstn;
}
