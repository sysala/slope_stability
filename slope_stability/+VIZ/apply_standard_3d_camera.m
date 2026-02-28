function apply_standard_3d_camera(ax, verts, opts)
%APPLY_STANDARD_3D_CAMERA Unified 3D camera setup for VIZ plots.
%
%   VIZ.apply_standard_3d_camera(ax, verts)
%   VIZ.apply_standard_3d_camera(ax, verts, opts)
%
% Inputs:
%   ax    - target axes handle
%   verts - N x 3 array of plotted 3D vertices
%   opts  - optional struct:
%       projection      ('orthographic' default)
%       view_dir        ([-1.2, 1, -2] default)
%       cam_dist_factor (2.3 default)
%       pad_factor      (0.05 default)
%       ydir            ('normal' default)

if nargin < 3
    opts = struct();
end
if size(verts, 2) ~= 3
    error('verts must be N x 3.');
end

projection = local_get_opt(opts, 'projection', 'orthographic');
view_dir = local_get_opt(opts, 'view_dir', [-1.2, 1, -2]);
cam_dist_factor = local_get_opt(opts, 'cam_dist_factor', 2.3);
pad_factor = local_get_opt(opts, 'pad_factor', 0.05);
ydir_mode = local_get_opt(opts, 'ydir', 'normal');

mins = min(verts, [], 1);
maxs = max(verts, [], 1);
span = max(maxs - mins);
if span <= 0
    span = 1;
end

pad = pad_factor * span;
xlim(ax, [mins(1) - pad, maxs(1) + pad]);
ylim(ax, [mins(2) - pad, maxs(2) + pad]);
zlim(ax, [mins(3) - pad, maxs(3) + pad]);

set(ax, 'Projection', projection);
set(ax, 'Ydir', ydir_mode);
axis(ax, 'equal');
axis(ax, 'tight');

center = (mins + maxs) / 2;
if norm(view_dir) <= eps
    view_dir = [-1.2, 1, -2];
end
view_dir = view_dir / norm(view_dir);
cam_pos = center + cam_dist_factor * span * view_dir;
cam_tgt = center;
set(ax, 'CameraPosition', cam_pos);
set(ax, 'CameraTarget', cam_tgt);

% Enforce x-horizontal and y-vertical image orientation.
d = cam_tgt - cam_pos;
nd = norm(d);
if nd <= eps
    set(ax, 'CameraUpVector', [0 1 0]);
    set(ax, 'Xdir', 'normal');
    box(ax, 'off');
    return;
end
d = d / nd;

x_world = [1 0 0];
y_world = [0 1 0];

x_proj = x_world - dot(x_world, d) * d;
nx = norm(x_proj);
if nx <= 1e-10
    set(ax, 'CameraUpVector', [0 1 0]);
    set(ax, 'Xdir', 'normal');
    box(ax, 'off');
    return;
end
x_proj = x_proj / nx;

up = cross(x_proj, d);
nu = norm(up);
if nu <= 1e-10
    set(ax, 'CameraUpVector', [0 1 0]);
    set(ax, 'Xdir', 'normal');
    box(ax, 'off');
    return;
end
up = up / nu;

y_proj = y_world - dot(y_world, d) * d;
if dot(y_proj, up) < 0
    up = -up;
end
set(ax, 'CameraUpVector', up);

right = cross(d, up);
nr = norm(right);
if nr > 1e-10
    right = right / nr;
    if dot(x_proj, right) < 0
        set(ax, 'Xdir', 'reverse');
    else
        set(ax, 'Xdir', 'normal');
    end
else
    set(ax, 'Xdir', 'normal');
end

box(ax, 'off');

end

function val = local_get_opt(opts, name, default_val)
if isfield(opts, name) && ~isempty(opts.(name))
    val = opts.(name);
else
    val = default_val;
end
end
