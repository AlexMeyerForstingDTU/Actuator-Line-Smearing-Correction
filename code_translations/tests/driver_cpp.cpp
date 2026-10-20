// Test driver for the C++ translation.
//
// Usage: driver_cpp <case>
// Runs one test case in the current folder and writes results.txt.
// The correction itself writes the diagnostics files test.nwRR.blBB.
// The test cases and inputs repeat test_common.f90.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "al_nw_smearing_correction.hpp"

static const int NS = 63, NV = NS + 1, LEN2D = 181;
static const double PI = 4.0 * atan(1.0), D2R = PI / 180.0, R2D = 180.0 / PI;
static const double WIND = 8.0, A_IND = 0.3, RHO = 1.225;
static const double OMEGA = 9.2 * 2.0 * PI / 60.0;

// Chord and twist of the blade at the table radii
static const int NTAB = 17;
static const double RTAB[NTAB] = {2.8667, 5.6, 8.3333, 11.75, 15.85, 19.95, 24.05, 28.15, 32.25,
                                  36.35, 40.45, 44.55, 48.65, 52.75, 56.1667, 58.9, 61.6333};
static const double TTAB[NTAB] = {13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
                                  6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106};
static const double CTAB[NTAB] = {3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
                                  3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419};

// One test case: frozen inflow marched in time
struct TestCase {
    double eps, cone, dt;  // Smearing length [m], cone [deg], time step [s]
    int n_steps;
    bool l1point;          // Ang Li correction
    bool bent;             // Pre-bent, swept blade with out-of-plane offsets
    bool root_at_zero;     // First trailing point at r = 0
    double max_resid;      // Convergence limit [m/s]
    int pitch_step_at;     // From this step on: pitch = 2 deg (0 = never)
    int n_rotors, n_blades;
    double pitch_22;       // Pitch of rotor 2, blade 2 [deg]
};

static TestCase get_case(const std::string& name) {
    TestCase c = {10.0, 0.0, 0.1, 300, false, false, false, 1e-6, 0, 1, 1, 0.0};
    if (name == "eps10_dt001") {
        c.dt = 0.01;
        c.n_steps = 3000;
    } else if (name == "eps10") {
    } else if (name == "cone25") {
        c.cone = 2.5;
    } else if (name == "eps25") {
        c.eps = 2.5;
    } else if (name == "nocorr") {
        c.eps = 0.01;
    } else if (name == "bent") {
        c.cone = 2.5;
        c.bent = true;
        c.l1point = true;
    } else if (name == "pitchstep") {
        c.n_steps = 400;
        c.pitch_step_at = 301;
    } else if (name == "noconv") {
        c.n_steps = 2;
        c.max_resid = 0.0;
    } else if (name == "multi") {
        c.n_steps = 100;
        c.n_rotors = 2;
        c.n_blades = 2;
        c.pitch_22 = 1.0;
    } else if (name == "rootzero") {
        c.n_steps = 1;
        c.root_at_zero = true;
    } else {
        printf("Error: unknown test case.\n");
        exit(1);
    }
    return c;
}

// Piecewise linear table lookup, constant beyond the ends
static double table_interp(double x, const double* xt, const double* yt, int n) {
    if (x <= xt[0]) return yt[0];
    if (x >= xt[n - 1]) return yt[n - 1];
    int k = 1;
    while (x > xt[k]) k++;
    return yt[k - 1] + (x - xt[k - 1]) / (xt[k] - xt[k - 1]) * (yt[k] - yt[k - 1]);
}

// Geometry, polars and frozen cylindrical inflow of the test blade
struct TestBlade {
    Vec rv, rc, chord, twist, fbend, sweep, dsweep;
    std::vector<Vec3> cyl_cp, cyl_vtp, vcyl;
    Grid alpha_2d, cl_2d, cd_2d;  // [l][m]
};

static TestBlade build_blade(const TestCase& c) {
    TestBlade b;
    // 1. Trailing points (uniform) and control points in between
    double r0 = c.root_at_zero ? 0.0 : 1.5;
    b.rv.assign(NV, 0.0);
    for (int k = 0; k < NV; k++) b.rv[k] = r0 + (63.0 - r0) * double(k) / double(NV - 1);
    b.rc.assign(NS, 0.0);
    b.chord.assign(NS, 0.0);
    b.twist.assign(NS, 0.0);
    for (int l = 0; l < NS; l++) {
        b.rc[l] = 0.5 * (b.rv[l] + b.rv[l + 1]);
        b.chord[l] = table_interp(b.rc[l], RTAB, CTAB, NTAB);
        b.twist[l] = table_interp(b.rc[l], RTAB, TTAB, NTAB);
    }
    // 2. Blade shape: straight, or pre-bent and swept with out-of-plane offsets
    b.fbend.assign(NS, 0.0);
    b.sweep.assign(NS, 0.0);
    b.dsweep.assign(NS, 0.0);
    b.cyl_cp.assign(NS, Vec3{0.0, 0.0, 0.0});
    b.cyl_vtp.assign(NV, Vec3{0.0, 0.0, 0.0});
    for (int l = 0; l < NS; l++) b.cyl_cp[l][0] = b.rc[l];
    for (int k = 0; k < NV; k++) b.cyl_vtp[k][0] = b.rv[k];
    if (c.bent) {
        for (int l = 0; l < NS; l++) {
            double t = b.rc[l] / 63.0;
            b.fbend[l] = 4.0 * t * t;
            b.sweep[l] = 3.0 * t * t * t;
            b.dsweep[l] = 0.5 * b.sweep[l];
            b.cyl_cp[l][2] = -1.5 * t * t;
        }
        for (int k = 0; k < NV; k++) {
            double t = b.rv[k] / 63.0;
            b.cyl_vtp[k][2] = -1.5 * t * t;
        }
    }
    // 3. Polars: cylinder inboard of 10 m, smooth thin-airfoil-like outboard
    b.alpha_2d.assign(NS, Vec(LEN2D, 0.0));
    b.cl_2d.assign(NS, Vec(LEN2D, 0.0));
    b.cd_2d.assign(NS, Vec(LEN2D, 0.0));
    for (int l = 0; l < NS; l++) {
        for (int m = 0; m < LEN2D; m++) {
            b.alpha_2d[l][m] = -90.0 + double(m);
            if (b.rc[l] < 10.0) {
                b.cl_2d[l][m] = 0.0;
                b.cd_2d[l][m] = 0.4;
            } else {
                b.cl_2d[l][m] = 1.4 * tanh(2.0 * PI * (b.alpha_2d[l][m] + 3.0) * D2R / 1.4);
                b.cd_2d[l][m] = 0.01;
            }
        }
    }
    // 4. Frozen "CFD" inflow in the cylindrical frame: uniform a = 0.3
    b.vcyl.assign(NS, Vec3{0.0, 0.0, WIND * (1.0 - A_IND)});
    return b;
}

// CFD velocities in section axes for a given cone and pitch
static std::vector<Vec3> cfd_vstn(const TestBlade& b, double cone, double pitch) {
    std::vector<Vec3> vstn(NS);
    for (int l = 0; l < NS; l++)
        vstn[l] = vxyz_vstn(b.vcyl[l], 0.0, 0.0, 0.0, cone, pitch, b.sweep[l], b.fbend[l]);
    return vstn;
}

// Rotor power and thrust from the section forces of one blade
static void power_thrust(const TestBlade& b, const std::vector<Vec3>& Fstn, double factor,
                         double& P, double& T) {
    P = 0.0;
    T = 0.0;
    for (int l = 0; l < NS; l++) {
        P = P + factor * Fstn[l][1] * OMEGA * b.rc[l] * (b.rv[l + 1] - b.rv[l]);
        T = T + factor * Fstn[l][2] * (b.rv[l + 1] - b.rv[l]);
    }
}

// Write one named result line: name v1 v2 ...
static void write_values(FILE* out, const char* name, const Vec& values) {
    fprintf(out, "%s", name);
    for (int n = 0; n < int(values.size()); n++) fprintf(out, " %.16e", values[n]);
    fprintf(out, "\n");
}

// Unit tests of the individual routines
static void run_units(FILE* out) {
    double phi, phi_s, dx, dy, X, Y;
    Coef4 a;
    Vec values;

    // 1. Smearing factor, with and without an out-of-plane offset
    values.push_back(smearing_factor(60.0, 0.1, 0.0, 0.0, 10.0));
    values.push_back(smearing_factor(60.0, 0.1, 2.0, 0.3, 10.0));
    values.push_back(smearing_factor(10.0, -0.5, 1.0, 0.5, 2.5));
    write_values(out, "smearing_factor", values);

    // 2. Geometry constants, r_vtp = 10 m: phi, phi_s, dx, dy, a(1:4) per h/r
    const double hr_list[7] = {0.1, 0.3, 0.5, 0.7, 0.8, -0.5, -3.0};
    values.clear();
    for (int n = 0; n < 7; n++) {
        nw_hr_components(10.0 * hr_list[n], hr_list[n], phi, phi_s, dx, dy, a);
        values.push_back(phi);
        values.push_back(phi_s);
        values.push_back(dx);
        values.push_back(dy);
        for (int m = 0; m < 4; m++) values.push_back(a[m]);
    }
    write_values(out, "hr_components", values);

    // 3. Phi* for four (h/r, helix angle) pairs
    const double phi_hr[4] = {0.1, 0.3, -0.5, -3.0};
    const double phi_deg[4] = {5.0, 30.0, 10.0, 60.0};
    values.clear();
    for (int n = 0; n < 4; n++) {
        nw_hr_components(10.0 * phi_hr[n], phi_hr[n], phi, phi_s, dx, dy, a);
        values.push_back(nw_phi_components(phi_hr[n], a, phi_deg[n] * D2R, phi, phi_s));
    }
    write_values(out, "phi_components", values);

    // 4. Indicial functions: one step, then the steady state of 2000 steps
    nw_xy_components(0.09, 0.0123, -0.0008125, 0.05, 2.0, 0.3, -0.02, X, Y);
    write_values(out, "xy_single", Vec{X, Y});
    X = 0.0;
    Y = 0.0;
    for (int n = 0; n < 2000; n++) {
        double X_hist = X, Y_hist = Y;
        nw_xy_components(0.09, 0.0123, -0.0008125, 0.05, 2.0, X_hist, Y_hist, X, Y);
    }
    write_values(out, "xy_steady", Vec{X, Y});

    // 5. Interpolation, inside and beyond both ends
    const Vec x_tab = {0.0, 1.0, 2.0, 3.0}, y_tab = {0.0, 10.0, 15.0, 30.0};
    const double x_in[7] = {-1.0, 0.0, 0.5, 2.0, 2.5, 3.0, 4.0};
    values.clear();
    for (int n = 0; n < 7; n++) values.push_back(nw_interp1d(x_tab, y_tab, x_in[n]));
    write_values(out, "interp1d", values);

    // 6. nw2vstn called repeatedly on a coned blade (note Sect. 4.1 and 6.1)
    std::vector<Vec3> v1;
    values.clear();
    for (int n = 0; n < 15; n++) {
        v1 = nw2vstn(Vec{0.0}, Vec{-1.0}, 2.5, 0.0, Vec{0.0}, Vec{0.0});
        values.push_back(v1[0][2]);
    }
    write_values(out, "nw2vstn_repeat", values);
    write_values(out, "nw2vstn_repeat_last", Vec{v1[0][0], v1[0][1], v1[0][2]});

    // 7. nw2vstn and the full rotation chain with all angles non-zero
    std::vector<Vec3> v2 = nw2vstn(Vec{1.5, -0.3}, Vec{-0.7, 2.0}, 2.5, 3.0, Vec{4.0, -2.0},
                                   Vec{-2.0, 5.0});
    write_values(out, "nw2vstn_general",
                 Vec{v2[0][0], v2[0][1], v2[0][2], v2[1][0], v2[1][1], v2[1][2]});
    Vec3 v = vxyz_vstn(Vec3{1.0, -2.0, 0.5}, 10.0, 5.0, 30.0, 2.5, 3.0, 4.0, -2.0);
    write_values(out, "vxyz_vstn", Vec{v[0], v[1], v[2]});

    // 8. Circulation of the test blade without correction (cone 0, pitch 0)
    TestBlade b = build_blade(get_case("eps10"));
    Vec un(NS), ut(NS), vrel(NS), alpha(NS), vt_pos(NS);
    for (int l = 0; l < NS; l++) {
        un[l] = WIND * (1.0 - A_IND);
        ut[l] = -(OMEGA * b.rc[l]);
        vrel[l] = sqrt(un[l] * un[l] + ut[l] * ut[l]);
        alpha[l] = atan2(un[l], -ut[l]) * R2D - b.twist[l];
        vt_pos[l] = -ut[l];
    }
    Vec G_steady, G, dG, G_hist(NS, 0.0);
    std::vector<Vec3> dyn, dyn_hist(NS, Vec3{0.0, 0.0, 0.0});
    nw_circulation(0.1, vrel, alpha, G_hist, dyn_hist, b.alpha_2d, b.cl_2d, b.chord, G_steady, G,
                   dG, dyn);
    write_values(out, "circulation_gamma_steady", G_steady);
    write_values(out, "circulation_gamma", G);
    write_values(out, "circulation_dgamma", dG);
    // Steady state of the time filter: Gamma -> Gamma_steady when sum(A_j) = 1
    for (int n = 0; n < 300; n++) {
        nw_circulation(1.0, vrel, alpha, G_hist, dyn_hist, b.alpha_2d, b.cl_2d, b.chord, G_steady,
                       G, dG, dyn);
        G_hist = G_steady;
        dyn_hist = dyn;
    }
    double ratio = 0.0;
    for (int l = 0; l < NS; l++)
        if (fabs(G_steady[l]) > 1e-3) ratio = std::max(ratio, fabs(G[l] / G_steady[l] - 1.0));
    write_values(out, "circulation_ratio", Vec{ratio});

    // 9. Relaxation factor of the test blade at two time steps
    Grid phi_pair(NV, Vec(NS)), dx_pair(NV, Vec(NS)), dy_pair(NV, Vec(NS));
    for (int k = 0; k < NV; k++)
        for (int l = 0; l < NS; l++)
            nw_hr_components(b.rv[k] - b.rc[l], (b.rv[k] - b.rc[l]) / b.rv[k], phi_pair[k][l],
                             phi_s, dx_pair[k][l], dy_pair[k][l], a);
    Vec vrel_vtp(NV), dbeta(NV);
    vrel_vtp[0] = vrel[0];
    vrel_vtp[NV - 1] = vrel[NS - 1];
    for (int k = 1; k < NV - 1; k++) vrel_vtp[k] = nw_interp1d(b.rc, vrel, b.rv[k]);
    const double relax_dt[2] = {0.01, 0.1};
    values.clear();
    for (int n = 0; n < 2; n++) {
        for (int k = 0; k < NV; k++) dbeta[k] = (vrel_vtp[k] * relax_dt[n]) / b.rv[k];
        values.push_back(nw_calc_relax(vrel, alpha, un, vt_pos, b.chord, relax_dt[n], dx_pair,
                                       dy_pair, dbeta, phi_pair, 0.2));
    }
    write_values(out, "calc_relax", values);
}

// March one test case in time and write its results
static void run_case(const TestCase& c, FILE* out) {
    // 1. Inputs and the correction state
    TestBlade b = build_blade(c);
    NWState state = init_nw_smearing_correction(c.n_rotors, c.n_blades, NS, "test");
    state.nw_max_resid = c.max_resid;

    // 2. March in time; F[i][j] holds the forces of rotor i + 1, blade j + 1
    std::vector<std::vector<std::vector<Vec3>>> F(2, std::vector<std::vector<Vec3>>(2));
    Vec P_series, T_series;
    for (int it = 1; it <= c.n_steps; it++) {
        double time = it * c.dt;
        double pitch = 0.0;
        if (c.pitch_step_at > 0 && it >= c.pitch_step_at) pitch = 2.0;
        for (int i = 1; i <= c.n_rotors; i++) {
            for (int j = 1; j <= c.n_blades; j++) {
                double p = pitch;
                if (i == 2 && j == 2) p = c.pitch_22;
                std::vector<Vec3> vstn = cfd_vstn(b, c.cone, p);
                nw_smearing_correction(state, i, j, time, c.dt, RHO, b.cyl_cp, b.cyl_vtp, vstn,
                                       b.vcyl, c.eps, OMEGA, b.chord, b.twist, p, b.cl_2d,
                                       b.cd_2d, b.alpha_2d, c.cone, b.fbend, b.sweep, b.dsweep,
                                       c.l1point, F[i - 1][j - 1]);
            }
        }
        if (c.pitch_step_at > 0 && it > c.n_steps - 110) {
            double P, T;
            power_thrust(b, F[0][0], 3.0, P, T);
            P_series.push_back(P);
            T_series.push_back(T);
        }
    }

    // 3. Results
    if (c.root_at_zero) {
        write_values(out, "kaya_after", Vec{state.kaya ? 1.0 : 0.0});
    } else if (c.n_rotors > 1) {
        Vec P_blade, T_blade;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                double P, T;
                power_thrust(b, F[i][j], 1.0, P, T);
                P_blade.push_back(P);
                T_blade.push_back(T);
            }
        }
        write_values(out, "P_blade", P_blade);
        write_values(out, "T_blade", T_blade);
    } else {
        double P, T;
        power_thrust(b, F[0][0], 3.0, P, T);
        write_values(out, "P", Vec{P});
        write_values(out, "T", Vec{T});
        Vec F2(NS), F3(NS);
        for (int l = 0; l < NS; l++) {
            F2[l] = F[0][0][l][1];
            F3[l] = F[0][0][l][2];
        }
        write_values(out, "Fstn2", F2);
        write_values(out, "Fstn3", F3);
        write_values(out, "vn_history", state.blade[0].vn_history);
        write_values(out, "vt_history", state.blade[0].vt_history);
        if (c.pitch_step_at > 0) {
            write_values(out, "P_series", P_series);
            write_values(out, "T_series", T_series);
        }
    }
}

int main(int argc, char** argv) {
    std::string name = argc > 1 ? argv[1] : "";
    FILE* out = fopen("results.txt", "w");
    if (name == "units")
        run_units(out);
    else
        run_case(get_case(name), out);
    fclose(out);
    return 0;
}
