function driver_matlab(case_name)
% Test driver for the MATLAB translation.
%   Runs one test case in the current folder and writes results.txt. The correction itself
%   writes the diagnostics files test.nwRR.blBB. The test cases and inputs repeat
%   test_common.f90. Needs the folder translations/matlab on the MATLAB path.
out = fopen('results.txt', 'w');
if strcmp(case_name, 'units')
    run_units(out);
else
    run_case(get_case(case_name), out);
end
fclose(out);
end

function p = test_parameters()
% Constants of the test blade and the inflow
p.ns = 63;                  % sections (control points)
p.nv = 64;                  % trailing points
p.len2d = 181;              % points per polar
p.wind = 8;                 % wind speed [m/s]
p.a_ind = 0.3;              % induction of the frozen "CFD" inflow [-]
p.rho = 1.225;              % air density [kg/m^3]
p.omega = 9.2*2*pi/60;      % rotor speed [rad/s]
end

function c = get_case(name)
% Parameters of a named test case: frozen inflow marched in time
c.eps = 10;                 % smearing length [m]
c.cone = 0;                 % cone [deg]
c.dt = 0.1;                 % time step [s]
c.n_steps = 300;
c.l1point = false;          % Ang Li correction
c.bent = false;             % pre-bent, swept blade with out-of-plane offsets
c.root_at_zero = false;     % first trailing point at r = 0
c.max_resid = 1e-6;         % convergence limit [m/s]
c.pitch_step_at = 0;        % from this step on: pitch = 2 deg (0 = never)
c.n_rotors = 1;
c.n_blades = 1;
c.pitch_22 = 0;             % pitch of rotor 2, blade 2 [deg]
switch name
    case 'eps10_dt001'
        c.dt = 0.01;
        c.n_steps = 3000;
    case 'eps10'
    case 'cone25'
        c.cone = 2.5;
    case 'eps25'
        c.eps = 2.5;
    case 'nocorr'
        c.eps = 0.01;
    case 'bent'
        c.cone = 2.5;
        c.bent = true;
        c.l1point = true;
    case 'pitchstep'
        c.n_steps = 400;
        c.pitch_step_at = 301;
    case 'noconv'
        c.n_steps = 2;
        c.max_resid = 0;
    case 'multi'
        c.n_steps = 100;
        c.n_rotors = 2;
        c.n_blades = 2;
        c.pitch_22 = 1;
    case 'rootzero'
        c.n_steps = 1;
        c.root_at_zero = true;
    otherwise
        error('Unknown test case: %s', name);
end
end

function y = table_interp(x, xt, yt)
% Piecewise linear table lookup, constant beyond the ends
if x <= xt(1)
    y = yt(1);
elseif x >= xt(end)
    y = yt(end);
else
    k = 2;
    while x > xt(k)
        k = k + 1;
    end
    y = yt(k-1) + (x - xt(k-1))/(xt(k) - xt(k-1))*(yt(k) - yt(k-1));
end
end

function b = build_blade(c, p)
% Geometry, polars and frozen cylindrical inflow of the test blade
rtab = [2.8667, 5.6, 8.3333, 11.75, 15.85, 19.95, 24.05, 28.15, 32.25, 36.35, 40.45, 44.55, ...
    48.65, 52.75, 56.1667, 58.9, 61.6333];
ttab = [13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, ...
    3.125, 2.319, 1.526, 0.863, 0.370, 0.106];
ctab = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, ...
    2.764, 2.518, 2.313, 2.086, 1.419];
ns = p.ns;
nv = p.nv;

%% 1. Trailing points (uniform) and control points in between
r0 = 1.5;
if c.root_at_zero
    r0 = 0;
end
b.rv = r0 + (63 - r0)*(0:nv-1)'/(nv - 1);
b.rc = 0.5*(b.rv(1:ns) + b.rv(2:nv));
b.chord = zeros(ns, 1);
b.twist = zeros(ns, 1);
for l = 1:ns
    b.chord(l) = table_interp(b.rc(l), rtab, ctab);
    b.twist(l) = table_interp(b.rc(l), rtab, ttab);
end

%% 2. Blade shape: straight, or pre-bent and swept with out-of-plane offsets
b.fbend = zeros(ns, 1);
b.sweep = zeros(ns, 1);
b.dsweep = zeros(ns, 1);
b.cyl_cp = [b.rc, zeros(ns, 2)];
b.cyl_vtp = [b.rv, zeros(nv, 2)];
if c.bent
    t = b.rc/63;
    b.fbend = 4*t.*t;
    b.sweep = 3*t.*t.*t;
    b.dsweep = 0.5*b.sweep;
    b.cyl_cp(:, 3) = -1.5*t.*t;
    t = b.rv/63;
    b.cyl_vtp(:, 3) = -1.5*t.*t;
end

%% 3. Polars: cylinder inboard of 10 m, smooth thin-airfoil-like outboard
b.alpha_2d = repmat((-90:90)', 1, ns);
b.cl_2d = zeros(p.len2d, ns);
b.cd_2d = zeros(p.len2d, ns);
for l = 1:ns
    if b.rc(l) < 10
        b.cd_2d(:, l) = 0.4;
    else
        b.cl_2d(:, l) = 1.4*tanh(2*pi*(b.alpha_2d(:, l) + 3)*(pi/180)/1.4);
        b.cd_2d(:, l) = 0.01;
    end
end

%% 4. Frozen "CFD" inflow in the cylindrical frame: uniform a = 0.3
b.vcyl = [zeros(ns, 2), p.wind*(1 - p.a_ind)*ones(ns, 1)];
end

function [P, T] = power_thrust(b, Fstn, factor, omega)
% Rotor power and thrust from the section forces of one blade
P = 0;
T = 0;
for l = 1:size(Fstn, 1)
    P = P + factor*Fstn(l, 2)*omega*b.rc(l)*(b.rv(l+1) - b.rv(l));
    T = T + factor*Fstn(l, 3)*(b.rv(l+1) - b.rv(l));
end
end

function write_values(out, name, values)
% Write one named result line: name v1 v2 ...
fprintf(out, '%s', name);
fprintf(out, ' %.16e', values);
fprintf(out, '\n');
end

function run_units(out)
% Unit tests of the individual routines
p = test_parameters();
d2r = pi/180;
r2d = 180/pi;

%% 1. Smearing factor, with and without an out-of-plane offset
write_values(out, 'smearing_factor', [smearing_factor(60, 0.1, 0, 0, 10), ...
    smearing_factor(60, 0.1, 2, 0.3, 10), smearing_factor(10, -0.5, 1, 0.5, 2.5)]);

%% 2. Geometry constants, r_vtp = 10 m: phi, phi_s, dx, dy, a(1:4) per h/r
hr_list = [0.1, 0.3, 0.5, 0.7, 0.8, -0.5, -3];
values = zeros(1, 56);
for n = 1:7
    [phi, phi_s, dx, dy, a] = nw_hr_components(10*hr_list(n), hr_list(n));
    values(8*n-7:8*n) = [phi, phi_s, dx, dy, a'];
end
write_values(out, 'hr_components', values);

%% 3. Phi* for four (h/r, helix angle) pairs
phi_hr = [0.1, 0.3, -0.5, -3];
phi_deg = [5, 30, 10, 60];
values = zeros(1, 4);
for n = 1:4
    [phi, phi_s, ~, ~, a] = nw_hr_components(10*phi_hr(n), phi_hr(n));
    values(n) = nw_phi_components(phi_hr(n), a, phi_deg(n)*d2r, phi, phi_s);
end
write_values(out, 'phi_components', values);

%% 4. Indicial functions: one step, then the steady state of 2000 steps
[X, Y] = nw_xy_components(0.09, 0.0123, -0.0008125, 0.05, 2, 0.3, -0.02);
write_values(out, 'xy_single', [X, Y]);
X = 0;
Y = 0;
for n = 1:2000
    [X, Y] = nw_xy_components(0.09, 0.0123, -0.0008125, 0.05, 2, X, Y);
end
write_values(out, 'xy_steady', [X, Y]);

%% 5. Interpolation, inside and beyond both ends
write_values(out, 'interp1d', interp1([0, 1, 2, 3], [0, 10, 15, 30], [-1, 0, 0.5, 2, 2.5, 3, 4], ...
    'linear', 'extrap'));

%% 6. nw2vstn called repeatedly on a coned blade (note Sect. 4.1 and 6.1)
values = zeros(1, 15);
for n = 1:15
    v1 = nw2vstn(0, -1, 2.5, 0, 0, 0);
    values(n) = v1(1, 3);
end
write_values(out, 'nw2vstn_repeat', values);
write_values(out, 'nw2vstn_repeat_last', v1(1, :));

%% 7. nw2vstn and the full rotation chain with all angles non-zero
v2 = nw2vstn([1.5; -0.3], [-0.7; 2], 2.5, 3, [4; -2], [-2; 5]);
write_values(out, 'nw2vstn_general', [v2(1, :), v2(2, :)]);
write_values(out, 'vxyz_vstn', vxyz_vstn([1, -2, 0.5], 10, 5, 30, 2.5, 3, 4, -2));

%% 8. Circulation of the test blade without correction (cone 0, pitch 0)
b = build_blade(get_case('eps10'), p);
un = p.wind*(1 - p.a_ind)*ones(p.ns, 1);
ut = -(p.omega*b.rc);
vrel = sqrt(un.^2 + ut.^2);
alpha = atan2(un, -ut)*r2d - b.twist;
[G_steady, G, dG] = nw_circulation(0.1, vrel, alpha, zeros(p.ns, 1), zeros(p.ns, 3), ...
    b.alpha_2d, b.cl_2d, b.chord);
write_values(out, 'circulation_gamma_steady', G_steady);
write_values(out, 'circulation_gamma', G);
write_values(out, 'circulation_dgamma', dG);
% Steady state of the time filter: Gamma -> Gamma_steady when sum(A_j) = 1
G_hist = zeros(p.ns, 1);
dyn_hist = zeros(p.ns, 3);
for n = 1:300
    [G_steady, G, ~, dyn] = nw_circulation(1, vrel, alpha, G_hist, dyn_hist, b.alpha_2d, ...
        b.cl_2d, b.chord);
    G_hist = G_steady;
    dyn_hist = dyn;
end
lifting = abs(G_steady) > 1e-3;
write_values(out, 'circulation_ratio', max(abs(G(lifting)./G_steady(lifting) - 1)));

%% 9. Relaxation factor of the test blade at two time steps
phi_pair = zeros(p.nv, p.ns);
dx_pair = zeros(p.nv, p.ns);
dy_pair = zeros(p.nv, p.ns);
for k = 1:p.nv
    for l = 1:p.ns
        h = b.rv(k) - b.rc(l);
        [phi_pair(k, l), ~, dx_pair(k, l), dy_pair(k, l)] = nw_hr_components(h, h/b.rv(k));
    end
end
vrel_vtp = [vrel(1); interp1(b.rc, vrel, b.rv(2:p.nv-1), 'linear', 'extrap'); vrel(p.ns)];
relax_dt = [0.01, 0.1];
values = zeros(1, 2);
for n = 1:2
    dbeta = (vrel_vtp*relax_dt(n))./b.rv;
    values(n) = nw_calc_relax(vrel, alpha, un, -ut, b.chord, relax_dt(n), dx_pair, dy_pair, ...
        dbeta, phi_pair, 0.2);
end
write_values(out, 'calc_relax', values);
end

function run_case(c, out)
% March one test case in time and write its results
p = test_parameters();

%% 1. Inputs and the correction state
b = build_blade(c, p);
state = init_nw_smearing_correction(c.n_rotors, c.n_blades, p.ns, 'test');
state.nw_max_resid = c.max_resid;

%% 2. March in time; F{i, j} holds the forces of rotor i, blade j
F = cell(2, 2);
P_series = zeros(1, 110);
T_series = zeros(1, 110);
n_series = 0;
for it = 1:c.n_steps
    time = it*c.dt;
    pitch = 0;
    if c.pitch_step_at > 0 && it >= c.pitch_step_at
        pitch = 2;
    end
    for i = 1:c.n_rotors
        for j = 1:c.n_blades
            pitch_ij = pitch;
            if i == 2 && j == 2
                pitch_ij = c.pitch_22;
            end
            vstn = vxyz_vstn(b.vcyl, 0, 0, 0, c.cone, pitch_ij, b.sweep, b.fbend);
            [F{i, j}, state] = nw_smearing_correction(state, i, j, time, c.dt, p.rho, b.cyl_cp, ...
                b.cyl_vtp, vstn, b.vcyl, c.eps, p.omega, b.chord, b.twist, pitch_ij, b.cl_2d, ...
                b.cd_2d, b.alpha_2d, c.cone, b.fbend, b.sweep, b.dsweep, c.l1point);
        end
    end
    if c.pitch_step_at > 0 && it > c.n_steps - 110
        n_series = n_series + 1;
        [P_series(n_series), T_series(n_series)] = power_thrust(b, F{1, 1}, 3, p.omega);
    end
end

%% 3. Results
if c.root_at_zero
    write_values(out, 'kaya_after', double(state.kaya));
elseif c.n_rotors > 1
    P_blade = zeros(1, 4);
    T_blade = zeros(1, 4);
    for i = 1:2
        for j = 1:2
            [P_blade(2*i+j-2), T_blade(2*i+j-2)] = power_thrust(b, F{i, j}, 1, p.omega);
        end
    end
    write_values(out, 'P_blade', P_blade);
    write_values(out, 'T_blade', T_blade);
else
    [P, T] = power_thrust(b, F{1, 1}, 3, p.omega);
    write_values(out, 'P', P);
    write_values(out, 'T', T);
    write_values(out, 'Fstn2', F{1, 1}(:, 2));
    write_values(out, 'Fstn3', F{1, 1}(:, 3));
    write_values(out, 'vn_history', state.blade(1, 1).vn_history);
    write_values(out, 'vt_history', state.blade(1, 1).vt_history);
    if c.pitch_step_at > 0
        write_values(out, 'P_series', P_series);
        write_values(out, 'T_series', T_series);
    end
end

%% 4. Close the diagnostics files
for i = 1:c.n_rotors
    for j = 1:c.n_blades
        fclose(state.blade(i, j).fid);
    end
end
end
