function relax = nw_calc_relax(vrel, alpha, vn, vt, chord, dt, dx, dy, dbeta, phi, safety)
% Relaxation factor of the sub-iteration, Pirrung et al. (2017a) Eq. (27)-(31).
%   Uses the corrected Eq. (28). alpha [deg]; vn and vt are the normal and tangential speeds
%   (vt > 0); dx, dy and phi (root-corrected Phi_C) are n_v x n_s; safety: safety margin.
relax = 0;
for l = 1:length(vrel)
    % Gamma_dyn/Gamma_QS: Eq. (15) averaged over the time step
    tau = (2*vrel(l))/chord(l);
    d = 1 - 0.5547/(tau*dt*0.3064)*(1 - exp(-tau*dt*0.3064)) ...
        - 0.1828/(tau*dt*0.0439)*(1 - exp(-tau*dt*0.0439)) ...
        - 0.2625/(tau*dt*3.227)*(1 - exp(-tau*dt*3.227));
    % B_v of the inner (l) and outer (l + 1) vortex of section l, Eq. (28)
    A1 = dx(l, l)*phi(l, l)*(1 - exp(-dbeta(l)/phi(l, l))) ...
        + dy(l, l)*phi(l, l)*(1 - exp(-4*dbeta(l)/phi(l, l)));
    A2 = dx(l+1, l)*phi(l+1, l)*(1 - exp(-dbeta(l+1)/phi(l+1, l))) ...
        + dy(l+1, l)*phi(l+1, l)*(1 - exp(-4*dbeta(l+1)/phi(l+1, l)));
    % Gradient, Eq. (27): alpha converted to radians
    gradw = d*pi*chord(l)*(A1 - A2)*((alpha(l)*pi/180*vn(l))/vrel(l) ...
        + vrel(l)/(abs(vt(l))*((vn(l)/abs(vt(l)))^2 + 1)));
    relax = max(relax, -(1 + gradw)/(1 - gradw));     % Eq. (31)
end
relax = relax + safety*(1 - relax);
relax = min(relax, 0.99);
end
