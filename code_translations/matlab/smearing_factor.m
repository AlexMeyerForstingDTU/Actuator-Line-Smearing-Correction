function smear_fac = smearing_factor(r, h_r, delta, phi, epsilon)
% Smearing factor of a vortex at its trailing point (beta = 0).
%   Meyer Forsting et al. (2019) Eq. (21)-(22) with beta = 0, as in Meyer Forsting et al.
%   (2020) Sect. 2.3. r: trailing-point radius [m], h_r = h/r, delta: distance normal to the
%   rotor disc (+ve upstream) [m], phi: helix angle [rad], epsilon: smearing length [m].
dist_orth_sq = (r*h_r).^2 + (delta.*cos(phi)).^2;
smear_fac = exp(-dist_orth_sq/epsilon^2);
end
