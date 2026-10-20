function [X, Y] = nw_xy_components(phi_star, dx, dy, dbeta, dGamma, X_history, Y_history)
% Advance the slow (X) and fast (Y) indicial functions by one time step.
%   Pirrung et al. (2016) Eq. (25)-(26) with Phi* of (2017b): new element plus the decayed
%   history. Works element-wise on arrays; dbeta and dGamma are scalars.
X_exp = exp(-dbeta./phi_star);
Y_exp = exp(-4*dbeta./phi_star);
X = dGamma*dx.*phi_star.*(1 - X_exp) + X_history.*X_exp;
Y = dGamma*dy.*phi_star.*(1 - Y_exp) + Y_history.*Y_exp;
end
