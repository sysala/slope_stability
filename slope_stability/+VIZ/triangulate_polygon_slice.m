function [TR2, P3] = triangulate_polygon_slice(poly3, h, opts)
%TRIANGULATE_POLYGON_SLICE  Constrained Delaunay mesh of a planar 3D polygon slice.
%
% [TR2, P3] = triangulate_polygon_slice(poly3, h, opts)
% - poly3: (3xN) CLOSED polygon (first==last), one coord constant
% - h    : target size
% - opts.grid_factor (default 1/sqrt(2)), opts.plot2d=false, opts.plot3d=false
% - returns: TR2 (2D triangulation), P3 (3xNpts 3D node coords with fixed coord)

if nargin < 3, opts = struct(); end
grid_factor = get_opt(opts,'grid_factor', 1/sqrt(2));
plot2d      = get_opt(opts,'plot2d',      false);
plot3d      = get_opt(opts,'plot3d',      false);
snap_tol    = get_opt(opts,'snap_tol',    []);

% --- ensure closed, then drop the duplicate last vertex
if any(poly3(:,1) ~= poly3(:,end)), error('poly3 must be closed (first==last).'); end
poly3u = poly3(:,1:end-1);

% --- detect fixed axis and project to 2D
s = [std(poly3u(1,:)), std(poly3u(2,:)), std(poly3u(3,:))];
[~, fixed_ax] = min(s);
free_axes = setdiff(1:3, fixed_ax);
c_fixed = poly3u(fixed_ax,1);
Pbd = poly3u(free_axes,:).';        % Nb x 2

% --- tolerances
bbmin = min(Pbd,[],1); bbmax = max(Pbd,[],1);
bbox_diag = max(1, norm(bbmax - bbmin));
if isempty(snap_tol), snap_tol = 1e-12 * bbox_diag; end

% --- clean consecutive duplicates
Pbd = dedup_consecutive(Pbd, snap_tol);
Nb = size(Pbd,1);  if Nb < 3, error('Polygon needs at least 3 distinct vertices.'); end

% --- densify boundary to spacing ~ h (no duplicate first vertex at end)
[BP, C] = densify_boundary_no_dup(Pbd, h);

% --- interior sampling ~ h/sqrt(2)
hx = grid_factor * h; hy = grid_factor * h;
xv = bbmin(1):hx:bbmax(1);
yv = bbmin(2):hy:bbmax(2);
[Xg, Yg] = meshgrid(xv, yv);

% polyshape & inside mask (VECTOR inputs to support older MATLAB)
ps = polyshape(Pbd(:,1), Pbd(:,2), 'Simplify', true);
inside_vec = isinterior(ps, Xg(:), Yg(:));
Gin = [Xg(inside_vec), Yg(inside_vec)];

% fallback: if no interior grid point landed inside, add centroid
if isempty(Gin)
    Gin = mean(Pbd,1);
end

% --- constrained Delaunay
AllPts = [BP; Gin];               % first |BP| are boundary points
dt = delaunayTriangulation(AllPts, C);

% --- keep triangles inside polygon (incenters test)
tri = triangulation(dt.ConnectivityList, dt.Points);
ic = incenter(tri);
keep = isinterior(ps, ic(:,1), ic(:,2));
T = dt.ConnectivityList(keep,:);
TR2 = triangulation(T, dt.Points);

% --- lift to 3D
Npts = size(TR2.Points,1);
P3 = zeros(3, Npts);
P3(free_axes,:) = TR2.Points.';   % fill free axes
P3(fixed_ax,:)  = c_fixed;        % constant coord

% --- optional plots
if plot2d
    figure; triplot(TR2); axis equal; title('2D constrained Delaunay');
end

% ===== helpers =====
function val = get_opt(s, name, def)
    if isfield(s,name) && ~isempty(s.(name)), val = s.(name); else, val = def; end
end

function Q = dedup_consecutive(P, tol)
    keep = true(size(P,1),1);
    for i = 2:size(P,1)
        if norm(P(i,:) - P(i-1,:)) <= tol, keep(i) = false; end
    end
    Q = P(keep,:);
    if size(Q,1) > 1 && norm(Q(end,:) - Q(1,:)) <= tol
        Q(end,:) = []; % keep only one of the wrap endpoints
    end
end

function [BP, C] = densify_boundary_no_dup(P, hseg)
    % P: Nb x 2, open polygon (first!=last)
    Nb = size(P,1);
    BP = P(1,:);          % start with the first point
    C  = [];              % constraints (edges between rows of BP)
    last_idx = 1;
    for i = 1:Nb
        j = i+1; if j > Nb, j = 1; end
        v1 = P(i,:); v2 = P(j,:);
        L  = norm(v2 - v1);
        nseg = max(1, ceil(L / max(hseg, eps)));
        for k = 1:nseg
            t = k / nseg;
            newp = (1-t)*v1 + t*v2;
            if (i == Nb) && (k == nseg)
                % close to the very first vertex, do NOT add a duplicate point
                new_idx = 1;
            else
                BP = [BP; newp]; %#ok<AGROW>
                new_idx = size(BP,1);
            end
            C = [C; last_idx, new_idx]; %#ok<AGROW>
            last_idx = new_idx;
        end
    end
end

end
