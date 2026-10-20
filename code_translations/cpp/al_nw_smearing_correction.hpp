// Fast smearing correction for the actuator line (C++ version).
//
// Translated from al_nw_smearing_correction_v3.f. Method: Meyer Forsting et al.
// (2019, 2020), with the near-wake model of Pirrung et al. (2016, 2017a, 2017b)
// and the 2020 corrigenda. Equation numbers refer to these papers.
//
// Usage:
//   NWState state = init_nw_smearing_correction(n_rotors, n_blades, n_sections, project);
//   then every time step, for each rotor i and blade j (both from 1):
//   nw_smearing_correction(state, i, j, time, dt, density, ..., Fstn);
//
// Indexing (from 0): sections [l], trailing points [k], vortex-section pairs [k][l],
// vectors [l][0..2], polars [l][m] (one table per section).
#pragma once

#include <array>
#include <cstdio>
#include <string>
#include <vector>

typedef std::vector<double> Vec;
typedef std::vector<Vec> Grid;          // 2D array: pairs [k][l] or polars [l][m]
typedef std::array<double, 3> Vec3;
typedef std::array<double, 4> Coef4;

// History and geometry constants of one rotor blade
struct BladeState {
    bool refresh;                       // Recompute the geometry constants
    std::vector<int> cp_first, cp_last; // First and last CP felt by each vortex
    Vec W_history, vn_history, vt_history, Gamma_history;
    std::vector<Vec3> dyn_comp_history;
    Grid X_history, Y_history;
    Grid phi, phi_s, dx, dy;            // phi = root-corrected Phi_C
    std::vector<std::vector<Coef4>> a;  // Fit coefficients a_hr,i
    Grid smear_fac;
    FILE* fid;                          // Diagnostics file, NULL = none
};

// Settings and the state of all rotor blades
struct NWState {
    bool kaya;                          // Switch: correction on/off
    double nw_max_resid;                // Convergence limit of the correction velocity [m/s]
    double x_aoa;                       // Chordwise position of the AoA in the CFD [-]
    int n_s, n_v, n_blades;             // Sections (CPs), trailing points, blades per rotor
    std::vector<BladeState> blade;      // Blade (i, j) is blade[(i - 1)*n_blades + j - 1]
};

// Create the state; open the diagnostics files <project>.nwRR.blBB ("" = none)
NWState init_nw_smearing_correction(int n_rotors, int n_blades, int n_sections,
                                    const std::string& project);

// Corrected forces Fstn[l] = (spanwise, tangential, normal) of rotor i, blade j.
// Positions (r, theta, z) in the cylindrical rotor frame; CFD velocities at the CPs
// in section axes (vstn_cfd) and in the cylindrical frame (vcyl_cfd); angles in
// degrees, omega in rad/s. Fstn is not changed when the correction is switched off.
void nw_smearing_correction(NWState& state, int i, int j, double time, double dt,
                            double density, const std::vector<Vec3>& cyl_cp,
                            const std::vector<Vec3>& cyl_vtp, const std::vector<Vec3>& vstn_cfd,
                            const std::vector<Vec3>& vcyl_cfd, double eps_cp, double omega,
                            const Vec& chord, const Vec& twist, double pitch,
                            const Grid& cl_2d, const Grid& cd_2d, const Grid& alpha_2d,
                            double cone, const Vec& fbend, const Vec& sweep, const Vec& dsweep,
                            bool l1point, std::vector<Vec3>& Fstn);

// Near-wake model components
double smearing_factor(double r, double h_r, double delta, double phi, double epsilon);
void nw_hr_components(double h, double h_r, double& phi, double& phi_s, double& dx, double& dy,
                      Coef4& a);
double nw_phi_components(double h_r, const Coef4& a, double phi, double phi_c, double phi_s);
void nw_xy_components(double phi_star, double dx, double dy, double dbeta, double dGamma,
                      double X_history, double Y_history, double& X, double& Y);
void nw_circulation(double dt, const Vec& vrel, const Vec& alpha, const Vec& Gamma_history,
                    const std::vector<Vec3>& dyn_comp_history, const Grid& alpha_2d,
                    const Grid& cl_2d, const Vec& chord, Vec& Gamma_steady, Vec& Gamma,
                    Vec& dGamma, std::vector<Vec3>& dyn_comp);
double nw_calc_relax(const Vec& vrel, const Vec& alpha, const Vec& vn, const Vec& vt,
                     const Vec& chord, double dt, const Grid& dx, const Grid& dy,
                     const Vec& dbeta, const Grid& phi, double safety);

// Interpolation and coordinate transformations (angles in degrees)
double nw_interp1d(const Vec& x, const Vec& y, double xin);
Vec3 vxyz_vstn(const Vec3& a, double yaw, double tilt, double azimut, double cone, double pitch,
               double sweep, double fbend);
std::vector<Vec3> nw2vstn(const Vec& vt, const Vec& vn, double cone, double pitch,
                          const Vec& sweep, const Vec& fbend);
