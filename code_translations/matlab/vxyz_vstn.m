function v = vxyz_vstn(v, yaw, tilt, azimut, cone, pitch, sweep, fbend)
% Rotate vectors (rows of v) from the Cartesian or cylindrical frame to section axes.
%   Order: yaw, tilt, azimuth, cone, pitch, sweep, bend. Angles in degrees; sweep and fbend
%   can be given per row (n x 1).
v = rot_y(v, yaw);
v = rot_x(v, tilt);
v = rot_z(v, -azimut);
v = rot_y(v, cone);
v = rot_x(v, -pitch);
v = rot_z(v, -sweep);
v = rot_y(v, fbend);
end

function v = rot_x(v, angle)
% Rotation about the first axis (tilt, pitch)
c = cos(angle*pi/180);
s = sin(angle*pi/180);
v = [v(:, 1), c.*v(:, 2) + s.*v(:, 3), -s.*v(:, 2) + c.*v(:, 3)];
end

function v = rot_y(v, angle)
% Rotation about the second axis (yaw, cone, bend)
c = cos(angle*pi/180);
s = sin(angle*pi/180);
v = [c.*v(:, 1) - s.*v(:, 3), v(:, 2), s.*v(:, 1) + c.*v(:, 3)];
end

function v = rot_z(v, angle)
% Rotation about the third axis (azimuth, sweep)
c = cos(angle*pi/180);
s = sin(angle*pi/180);
v = [c.*v(:, 1) - s.*v(:, 2), s.*v(:, 1) + c.*v(:, 2), v(:, 3)];
end
