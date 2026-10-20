function phi_star = nw_phi_components(h_r, a, phi, phi_c, phi_s)
% Phi* for the helix angle phi [rad]: weighted mix of the straight and the circular Phi.
%   Pirrung et al. (2017b) Eq. (8)-(11), with the root-corrected Phi_C of the corrigendum.
%   a: fit coefficients a_hr,i (4 values).
if h_r < 0
    k_phi = a(1) + a(2)*exp(a(3)*(pi/2 - phi)) + a(4)*exp(-8*(pi/2 - phi)) - a(2) - a(4);
else
    k_phi = a(1) + a(2)*phi + a(3)*phi^2 + a(4)*phi^3;
end
% The fit is only valid below 89 degrees: keep the weight in [0, 1]
k_phi = max(min(k_phi, 1), 0);
phi_star = k_phi*phi_s + (1 - k_phi)*phi_c;
end
