function out = plot_slice(coord, edges, y0, varargin)
%PLOT_Y_SLICE_FROM_EDGES Slice a 3D edge-boundary by the plane y = y0 and plot in x–z.
%
% out = plot_y_slice_from_edges(coord, edges, y0, 'Name', value, ...)
%
% Inputs
%   coord : 3 x N double  (rows are x; y; z)
%   edges : M x 2 int     (1-based indices into columns of coord)
%   y0    : scalar double (target y-plane)
%
% Name–Value options
%   'TolY'  : tolerance in y to treat as on-plane (default: 1e-9 * range(y), min 1e-12)
%   'TolXZ' : tolerance to de-duplicate x–z points (default: 1e-9 * max(range(x),range(z)))
%   'Axes'  : target axes handle (default: new figure)
%   'Plot'  : logical, whether to plot (default: true)
%
% Outputs (struct)
%   out.points_xz         : K x 2 intersection points [x z]
%   out.points_edge_ids   : K x 1 indices of edges that produced the points
%   out.segments_xz       : S x 4 segments lying in the plane [x1 z1 x2 z2]
%   out.segments_edge_ids : S x 1 indices of edges that lie in the plane
%   out.y0                : the slice value
%   out.toly, out.tolxz   : tolerances used
%
% Example:
%   out = plot_y_slice_from_edges(coord, edges_on_edge, 0.25);

% --- Parse & sanity -------------------------------------------------------
p = inputParser;
p.addParameter('TolY',  [], @(v) isempty(v) || isscalar(v));
p.addParameter('TolXZ', [], @(v) isempty(v) || isscalar(v));
p.addParameter('Axes',  [], @(h) isempty(h) || isgraphics(h,'axes'));
p.addParameter('Plot',  true, @(b) islogical(b) && isscalar(b));
p.parse(varargin{:});
TolY  = p.Results.TolY;
TolXZ = p.Results.TolXZ;
ax    = p.Results.Axes;
doplot = p.Results.Plot;

x = coord(1,:).';  % N x 1
y = coord(2,:).';
z = coord(3,:).';

if isempty(TolY)
    ry = max(y) - min(y);
    TolY = max(1e-9*max(ry,1), 1e-12);
end
if isempty(TolXZ)
    rx = max(x) - min(x);
    rz = max(z) - min(z);
    TolXZ = max(1e-9*max([rx, rz, 1]), 1e-12);
end

e1 = edges(:,1);  % M x 1
e2 = edges(:,2);

x1 = x(e1);  y1 = y(e1);  z1 = z(e1);
x2 = x(e2);  y2 = y(e2);  z2 = z(e2);

dy1 = y1 - y0;
dy2 = y2 - y0;

on1   = abs(dy1) <= TolY;
on2   = abs(dy2) <= TolY;
copl  = on1 & on2;                            % whole edge lies in plane
touch = xor(on1, on2);                        % one endpoint on plane
cross = (dy1 .* dy2) < 0;                     % strict sign change (not touching)

% --- Segments that lie entirely in the plane -----------------------------
seg_ids = find(copl);
segments_xz = [x1(copl) z1(copl) x2(copl) z2(copl)];

% --- Points from crossings and touches -----------------------------------
pt_edge_mask = (touch | cross);
ids = find(pt_edge_mask);
den = (y2 - y1);              % safe: for coplanar edges we didn't enter here

% param t for intersection p = p1 + t*(p2-p1)
t = (y0 - y1(pt_edge_mask)) ./ den(pt_edge_mask);

% For 'touch' cases where an endpoint is exactly on plane, t should be 0 or 1.
% Numerical safety: clamp t to [0,1].
t = max(0, min(1, t));

xi = x1(pt_edge_mask) + t .* (x2(pt_edge_mask) - x1(pt_edge_mask));
zi = z1(pt_edge_mask) + t .* (z2(pt_edge_mask) - z1(pt_edge_mask));
points_xz = [xi zi];

% De-duplicate intersection points within TolXZ (optional but cleaner)
if ~isempty(points_xz)
    [~, ia] = uniquetol(points_xz, TolXZ, 'ByRows', true);
    points_xz = points_xz(ia, :);
    ids       = ids(ia);
end

% --- Prepare outputs ------------------------------------------------------
out = struct();
out.points_xz = points_xz;
out.points_edge_ids = ids(:);
out.segments_xz = segments_xz;
out.segments_edge_ids = seg_ids(:);
out.y0 = y0;
out.toly = TolY;
out.tolxz = TolXZ;

% --- Plot -----------------------------------------------------------------
if doplot
    if isempty(ax) || ~isvalid(ax)
        figure('Name', sprintf('Slice y = %g', y0));
        ax = axes; %#ok<LAXES>
    end
    hold(ax, 'on');
    if ~isempty(segments_xz)
        % Plot segments that lie in the plane
        for k = 1:size(segments_xz,1)
            plot(ax, segments_xz(k,[1 3]), segments_xz(k,[2 4]), '-', 'LineWidth', 1.0);
        end
    end
    if ~isempty(points_xz)
        plot(ax, points_xz(:,1), points_xz(:,2), '.', 'MarkerSize', 12);
    end
    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'x');
    ylabel(ax, 'z');
    title(ax, sprintf('Intersection with plane y = %g', y0));
    legend(ax, {'segments on plane','edge crossings'}, 'Location','bestoutside');
    hold(ax, 'off');
end
end
