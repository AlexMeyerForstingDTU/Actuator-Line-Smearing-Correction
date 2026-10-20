function [Fstn, state] = nw_smearing_correction(state, i, j, time, dt, density, cyl_cp, cyl_vtp, ...
    vstn_cfd, vcyl_cfd, eps_cp, omega, chord, twist, pitch, cl_2d, cd_2d, alpha_2d, ...
    cone, fbend, sweep, dsweep, l1point)
% Corrected blade forces of rotor i, blade j for one time step.
%   The smeared actuator-line forces leave part of the induction out of the CFD velocities;
%   the near-wake model supplies that part and corrects the angle of attack.
%   Method: Meyer Forsting et al. (2019, 2020), near-wake model of Pirrung et al. (2016,
%   2017a, 2017b) with the 2020 corrigenda. Translated from al_nw_smearing_correction.f.
%
%   time, dt, density   time [s], time step [s], air density [kg/m^3]
%   cyl_cp, cyl_vtp     CP (n_s x 3) and trailing-point (n_v x 3) positions (r, theta, z)
%   vstn_cfd, vcyl_cfd  CFD velocities at the CPs in section axes (s, t, n) and cylindrical axes
%   eps_cp, omega       smearing length [m], rotor speed [rad/s]
%   chord, twist        chord [m] and twist [deg] at the CPs (n_s x 1)
%   pitch, cone         blade pitch (+ve to feather) and cone (+ve upstream) [deg]
%   cl_2d, cd_2d        polars (len_2d x n_s) at the angles alpha_2d [deg]
%   fbend, sweep, dsweep local bend, sweep and delta sweep [deg] (n_s x 1)
%   l1point             Ang Li correction on/off
%   Fstn                forces per unit length at the CPs: spanwise, tangential, normal (n_s x 3)

Fstn = [];
if ~state.kaya
    return
end
n_s = state.n_s;
n_v = state.n_v;
d2r = pi/180;
r2d = 180/pi;
r_cp = cyl_cp(:, 1);
r_vtp = cyl_vtp(:, 1);
% Important: r = 0 at the first trailing point is a singularity
if r_vtp(1) < 0.1
    fprintf('!!!! WARNING smearing correction switched off !!!!\n');
    fprintf('First vortex trailing point is zero, please modify\n');
    state.kaya = false;
    return
end
cpu_t1 = cputime;
b = state.blade(i, j);

%% 1. Preamble: start from the previous time step
relax_safety = 0.2;
recompute_relax = true;
W = b.W_history;
vn = b.vn_history;
vt = b.vt_history;
scone = sin((cone + fbend)*d2r);
v_rot = -(omega*r_cp);                      % rotational speed in section axes
vn_rot = -v_rot*sin(-pitch*d2r);
vt_rot = v_rot*cos(-pitch*d2r).*cos(dsweep*d2r);

%% 2. Geometry constants: first call, or after a refresh
if b.refresh
    core_cutoff = 1.82*eps_cp;              % 99th percentile of the Gaussian
    b.cp_loop_ind(:, 1) = 1;
    b.cp_loop_ind(:, 2) = n_s;
    for k = 1:n_v
        for l = 1:n_s
            h = r_vtp(k) - r_cp(l);
            h_r = h/r_vtp(k);
            [b.phi(k, l), b.phi_s(k, l), b.dx(k, l), b.dy(k, l), a] = nw_hr_components(h, h_r);
            b.a(k, l, :) = a;
            b.smear_fac(k, l) = smearing_factor(r_vtp(k), h_r, 0, 0, eps_cp);
            % CPs felt by vortex k: inboard and outboard ones within the cut-off
            if h > 0
                if abs(h) > core_cutoff
                    b.cp_loop_ind(k, 1) = l + 1;
                end
            elseif abs(h) < core_cutoff
                b.cp_loop_ind(k, 2) = l;
            end
        end
    end
    b.refresh = false;
end

%% 3. Helix angles from the CFD only (MF2019 step 2), fixed during the sub-iteration
phi_cp_cfd = atan2(vcyl_cfd(:, 3), -(v_rot + vcyl_cfd(:, 2)));
for l = 1:n_s
    if phi_cp_cfd(l) > pi/2
        phi_cp_cfd(l) = pi - phi_cp_cfd(l);
    end
end
phi_vtp = [phi_cp_cfd(1); interp1(r_cp, phi_cp_cfd, r_vtp(2:n_v-1), 'linear', 'extrap'); ...
    phi_cp_cfd(n_s)];
phi_star = zeros(n_v, n_s);
for k = 1:n_v
    for l = b.cp_loop_ind(k, 1):b.cp_loop_ind(k, 2)
        h = r_vtp(k) - r_cp(l);
        h_r = h/r_vtp(k);
        phi_star(k, l) = nw_phi_components(h_r, b.a(k, l, :), phi_vtp(k), b.phi(k, l), ...
            b.phi_s(k, l));
        delta = -(cyl_vtp(k, 3) - cyl_cp(l, 3));     % normal to the disc, +ve upstream
        b.smear_fac(k, l) = smearing_factor(r_vtp(k), h_r, delta, phi_vtp(k), eps_cp);
    end
end

%% 4. Sub-iterate until the correction velocities converge
sub_iter = 0;
converged = false;
while ~converged
    % Velocities and angles with the current correction
    vstn = nw2vstn(vt, vn, cone, pitch, sweep, fbend);
    un = vn_rot + vstn_cfd(:, 3) + vstn(:, 3);       % normal velocity
    ut = vt_rot + vstn_cfd(:, 2) + vstn(:, 2);       % tangential velocity (< 0)
    phi_cp = atan2(un, -ut);
    alpha_cp = phi_cp*r2d - twist;
    vrel_cp = sqrt(un.^2 + ut.^2);
    vrel_vtp = [vrel_cp(1); interp1(r_cp, vrel_cp, r_vtp(2:n_v-1), 'linear', 'extrap'); ...
        vrel_cp(n_s)];
    [Gamma_steady, ~, dGamma, dyn_comp] = nw_circulation(dt, vrel_cp, alpha_cp, b.Gamma_history, ...
        b.dyn_comp_history, alpha_2d, cl_2d, chord);
    dbeta = (vrel_vtp*dt)./r_vtp;                  % advance in beta*, 2017b Eq. (3)

    % Induction of the vorticity shed this time step, smeared part only
    W_iter = zeros(n_s, 1);
    Wn_iter = zeros(n_s, 1);
    Wt_iter = zeros(n_s, 1);
    X = zeros(n_v, n_s);
    Y = zeros(n_v, n_s);
    for k = 1:n_v
        l = b.cp_loop_ind(k, 1):b.cp_loop_ind(k, 2);      % CPs felt by vortex k
        [X(k, l), Y(k, l)] = nw_xy_components(phi_star(k, l), b.dx(k, l), b.dy(k, l), dbeta(k), ...
            dGamma(k), b.X_history(k, l), b.Y_history(k, l));
        W_k = b.smear_fac(k, l).*(X(k, l) + Y(k, l));
        W_iter(l) = W_iter(l) + W_k(:);
        % Normal and tangential parts with each vortex's helix angle (MF2019 Eq. 19)
        Wn_iter(l) = Wn_iter(l) + W_k(:)*cos(phi_vtp(k));
        Wt_iter(l) = Wt_iter(l) + W_k(:)*sin(phi_vtp(k));
    end

    if recompute_relax
        relax = nw_calc_relax(vrel_cp, alpha_cp, un, -ut, chord, dt, b.dx, b.dy, dbeta, b.phi, ...
            relax_safety);
        recompute_relax = false;
    end

    % Relax and measure the change
    W_iter = relax*W + (1 - relax)*W_iter;
    section_res = sqrt((W_iter - W).^2);
    phi_cp_cyl = atan2(vcyl_cfd(:, 3) + vn, -(v_rot + vcyl_cfd(:, 2) + vt));
    % Tip vortex within 45 deg of the rotor plane: axial residual, else tangential
    if phi_vtp(n_v) < pi/4
        resmax = max(section_res.*cos(phi_cp_cyl));
    else
        resmax = max(section_res.*sin(phi_cp_cyl));
    end

    sub_iter = sub_iter + 1;
    converged = resmax < state.nw_max_resid || sub_iter == 500;
    if sub_iter == 100
        recompute_relax = true;
    end
    if sub_iter == 500
        if relax_safety < 0.8
            relax_safety = relax_safety + 0.1;
            relax = min(relax + 0.1*(1 - relax), 0.99);
            recompute_relax = true;
        end
        % No convergence: keep the near-wake induction of the last time step
        W = b.W_history;
        vn = b.vn_history;
        vt = b.vt_history;
        dyn_comp = b.dyn_comp_history;
        Gamma_steady = b.Gamma_history;
        X = b.X_history;
        Y = b.Y_history;
    else
        W = W_iter;
        vn = relax*vn + (1 - relax)*Wn_iter;
        vt = relax*vt + (1 - relax)*Wt_iter;
    end
end

%% 5. Blade forces with the corrected inflow
vstn = nw2vstn(vt, vn, cone, pitch, sweep, fbend);
un = vn_rot + vstn_cfd(:, 3) + vstn(:, 3);
ut = vt_rot + vstn_cfd(:, 2) + vstn(:, 2);
phi_cp = atan2(un, -ut);
alpha_cp = phi_cp*r2d - twist;
alpha_cp_cfd = phi_cp_cfd*r2d - twist;
vrel_cp = sqrt(un.^2 + ut.^2);
% Ang Li correction for coned or bent blades (Li et al., 2022): lift and drag at the
% 3/4-chord AoA, force direction at the 1/4-chord inflow angle
if l1point
    thetadot = -omega*scone*r2d;
    alpha_cp_34 = alpha_cp + (0.75 - state.x_aoa)*chord./vrel_cp.*thetadot;
    alpha_cp_14 = alpha_cp + (0.25 - state.x_aoa)*chord./vrel_cp.*thetadot;
else
    alpha_cp_34 = alpha_cp;
    alpha_cp_14 = alpha_cp;
end
phi_cp_14 = (alpha_cp_14 + twist)*d2r;
cl = zeros(n_s, 1);
cd = zeros(n_s, 1);
for l = 1:n_s
    cl(l) = interp1(alpha_2d(:, l), cl_2d(:, l), alpha_cp_34(l), 'linear', 'extrap');
    cd(l) = interp1(alpha_2d(:, l), cd_2d(:, l), alpha_cp_34(l), 'linear', 'extrap');
end
lift = 0.5*density*vrel_cp.^2.*chord.*cl;
drag = 0.5*density*vrel_cp.^2.*chord.*cd;
Fstn = zeros(n_s, 3);
Fstn(:, 2) = lift.*sin(phi_cp_14) - drag.*cos(phi_cp_14);
Fstn(:, 3) = lift.*cos(phi_cp_14) + drag.*sin(phi_cp_14);

%% 6. Diagnostics line and history for the next time step
if b.fid > 0
    fprintf(b.fid, '%15.7E%5d%15.7E', time, sub_iter, resmax);
    fprintf(b.fid, '%15.7E', [alpha_cp_cfd; alpha_cp; phi_vtp; vt; vn; dGamma]);
    fprintf(b.fid, '%15.7E\n', cputime - cpu_t1);
end
b.W_history = W;
b.dyn_comp_history = dyn_comp;
b.Gamma_history = Gamma_steady;
b.vn_history = vn;
b.vt_history = vt;
b.X_history = X;
b.Y_history = Y;
state.blade(i, j) = b;
end
