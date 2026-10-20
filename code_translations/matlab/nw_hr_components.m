function [phi, phi_s, dx, dy, a] = nw_hr_components(h, h_r)
% Constants of one vortex-section pair, h = r_vtp - r_cp and h_r = h/r_vtp.
%   phi: root-corrected Phi_C, phi_s: Phi of a straight vortex, dx and dy: amplitudes without
%   Gamma and Phi, a: fit coefficients a_hr,i of the mixing factor (4 x 1).
%   Pirrung et al. (2016) and (2017b).

% Fit constants of the mixing factor, Pirrung et al. (2017b). Row i gives a_hr,i.
% P: h/r > 0, Eq. (14). N: h/r < 0, Eq. (13).
P = [-1.64636754184988, 8.14821474772595, -12.1784861715161, 5.02652654794305, 21.7713111885100
    -0.499006672527007, 6.08465268675599, -15.1712005446736, 14.8254135311781, -2.42318681495122
    3.90835606394205, -18.7662330478906, 39.1243282362892, -29.4870086106427, -60.2947309276912
    -1.60622843516794, 7.42952711674470, -15.8594777074799, 11.6870163960092, -195.060865844227];
N = [1.01933206359188, -0.135668151345576, 0.395524705910666, 0.0801822997046969, 44.8347503471217
    12.9874512870379, 49.9999999999916, 0.00235345218088083, 11.3116074964299, 3935.34323353354
    -0.690159136848837, 101.238775997879, -0.00154467589177718, 3.99520376143011, 0.394541495293350
    -0.269253409815015, 49.9999999999843, -0.00247975579295973, 0.403642518874381, 1.16610462658971];

% Phi of Wang and Coton, 2016 Eq. (4)
if h_r > 0
    phi = pi/4*abs((1 + 0.5*h_r)*log(1 - h_r));
else
    phi = log(1 - h_r)/(1.5 + log(1 - 0.5*h_r));
end

% Amplitudes, 2016 Eq. (25). Important: D_Y is negative.
dw = 1/(4*pi*abs(h)*h_r);
dx = 1.359*dw;
dy = -0.359/4*dw;

% Root correction, 2016 Eq. (27)-(29). Important: applied once, through Phi (Phi_C, Eq. 31).
root_corr = (dx*(1 - exp(-pi/(2*phi))) + dy*(1 - exp(-2*pi/phi)))/(dx + dy);
phi = root_corr*phi;

% Phi of a straight vortex, 2017b Eq. (7)
phi_s = 0.788*abs(h_r);

% Fit coefficients of the mixing factor, 2017b Eq. (10) and (12)
if h_r < 0
    a = N(:, 1) + N(:, 2).*exp(N(:, 3)*h_r) + N(:, 4).*exp(N(:, 5)*h_r) - N(:, 4) - N(:, 2);
else
    a = P(:, 1) + P(:, 2)*h_r + P(:, 3)*h_r^2 + P(:, 4)*h_r^3;
end
end
