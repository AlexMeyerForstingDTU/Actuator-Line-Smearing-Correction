function state = init_nw_smearing_correction(n_rotors, n_blades, n_sections, project)
% Create the state of the smearing correction for n_rotors x n_blades actuator lines.
%   n_sections: number of control points (CP) per blade; trailing points: n_sections + 1.
%   project: prefix of the diagnostics files <project>.nwRR.blBB, '' for no output.
%   Translated from al_nw_smearing_correction.f (init_nw_smearing_correction).

%% Settings
state.kaya = true;            % switch: correction on/off
state.nw_max_resid = 1e-6;    % convergence limit of the correction velocity [m/s]
state.x_aoa = 0.5;            % chordwise position of the AoA in the CFD [-]
state.n_s = n_sections;
state.n_v = n_sections + 1;

%% History and geometry constants of one blade
n_s = state.n_s;
n_v = state.n_v;
blade.refresh = true;                      % recompute the geometry constants
blade.cp_loop_ind = [ones(n_v, 1), n_s*ones(n_v, 1)];   % first and last CP felt by each vortex
blade.W_history = zeros(n_s, 1);
blade.vn_history = zeros(n_s, 1);
blade.vt_history = zeros(n_s, 1);
blade.Gamma_history = zeros(n_s, 1);
blade.dyn_comp_history = zeros(n_s, 3);
blade.X_history = zeros(n_v, n_s);
blade.Y_history = zeros(n_v, n_s);
blade.phi = zeros(n_v, n_s);               % root-corrected Phi_C
blade.phi_s = zeros(n_v, n_s);
blade.dx = zeros(n_v, n_s);
blade.dy = zeros(n_v, n_s);
blade.a = zeros(n_v, n_s, 4);              % fit coefficients a_hr,i
blade.smear_fac = zeros(n_v, n_s);
blade.fid = -1;                            % diagnostics file, -1 = none

%% All blades, with their diagnostics files
for i = 1:n_rotors
    for j = 1:n_blades
        state.blade(i, j) = blade;
        if ~isempty(project)
            state.blade(i, j).fid = fopen(sprintf('%s.nw%02d.bl%02d', project, i, j), 'w');
        end
    end
end
end
