function vstn = nw2vstn(vt, vn, cone, pitch, sweep, fbend)
% Correction velocities (cylindrical vt, vn; n x 1) in section axes (n x 3).
vstn = vxyz_vstn([zeros(size(vt)), vt, vn], 0, 0, 0, cone, pitch, sweep, fbend);
end
