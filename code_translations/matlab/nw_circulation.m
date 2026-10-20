function [Gamma_steady, Gamma, dGamma, dyn_comp] = nw_circulation(dt, vrel, alpha, ...
    Gamma_history, dyn_comp_history, alpha_2d, cl_2d, chord)
% Bound circulation, its time-filtered value and the trailed vortex strengths.
%   vrel, alpha [deg], Gamma_history, chord: n_s x 1. dyn_comp_history: n_s x 3.
%   alpha_2d, cl_2d: polars (len_2d x n_s). dGamma: n_s + 1 trailed strengths.

%% Quasi-steady circulation from the static polar
n_s = length(vrel);
cl = zeros(n_s, 1);
for l = 1:n_s
    cl(l) = interp1(alpha_2d(:, l), cl_2d(:, l), alpha(l), 'linear', 'extrap');
end
Gamma_steady = 0.5*vrel.*chord.*cl;

%% Time filter
% Madsen and Gaunaa (2004), Pirrung et al. (2017a) Eq. (15)-(19). The constants are A_j/2 and
% b_j; A3 = 0.2625 (paper: 0.2656) so that the A_j sum to 1.
t = dt./(chord./(2*vrel));
dyn_comp = zeros(n_s, 3);
dyn_comp(:, 1) = dyn_comp_history(:, 1).*exp(-0.3064*t) ...
    + 0.27735*(Gamma_steady + Gamma_history).*(1 - exp(-0.3064*t));
dyn_comp(:, 2) = dyn_comp_history(:, 2).*exp(-0.0439*t) ...
    + 0.0914*(Gamma_steady + Gamma_history).*(1 - exp(-0.0439*t));
dyn_comp(:, 3) = dyn_comp_history(:, 3).*exp(-3.227*t) ...
    + 0.13125*(Gamma_steady + Gamma_history).*(1 - exp(-3.227*t));
Gamma = dyn_comp(:, 1) + dyn_comp(:, 2) + dyn_comp(:, 3);

%% Trailed strengths: the end vortices take all Gamma of the end sections
dGamma = [Gamma(1); diff(Gamma); -Gamma(n_s)];
end
