function [poly3, loops3] = slice_by_plane(faces, edges_merged, coord, plane_id, plane_val, opts)
%SLICE_DOMAIN_POLYGON  Intersect a 3D polyhedral boundary with an axis-aligned plane.
%
% [poly3, loops3] = slice_domain_polygon(faces, edges_merged, coord, plane_id, plane_val, opts)
%
% Inputs
%   faces         : struct from group_faces_by_edge_pairs (needs faces.edges{f})
%   edges_merged  : (E x 2) int, boundary edges (vertex indices, 1-based)
%   coord         : (3 x N) double, vertex coordinates
%   plane_id      : 1 -> x = plane_val (YZ plane)
%                   2 -> y = plane_val (XZ plane)
%                   3 -> z = plane_val (XY plane)
%   plane_val     : scalar value of the fixed coordinate
%   opts (optional):
%       .plane_tol   : tolerance for plane hit (default: 1e-9 * bbox_diag)
%       .snap_tol    : tolerance to glue endpoints into a loop (default: 1e-8 * bbox_diag)
%       .plot        : true/false to plot polygon (default: false)
%
% Outputs
%   poly3   : (3 x M) the selected polygon loop in 3D (closed: first==last)
%   loops3  : 1xL cell of all loops found (each 3 x Mi), before selection
%
% Algorithm (as requested):
%   - For each edge, test intersection with the plane (ignore edges fully in plane)
%   - For each face, if exactly two of its edges intersect → one segment between the two points
%   - Gather all such segments and connect them into closed loops
%   - Return the loop with the largest perimeter as poly3 (also return all loops)

if nargin < 6, opts = struct(); end

% tolerances scaled to geometry size
bbmin = min(coord,[],2); bbmax = max(coord,[],2);
bbox_diag = max(1, norm(bbmax - bbmin));
plane_tol = get_opt(opts,'plane_tol', 1e-9*bbox_diag);
snap_tol  = get_opt(opts,'snap_tol',  1e-8*bbox_diag);
do_plot   = get_opt(opts,'plot',     false);

% Axis triplet for projection (use the two free axes for 2D work)
switch plane_id
    case 1, axis_fixed = 1; proj = [2,3];  % x = c → work in (y,z)
    case 2, axis_fixed = 2; proj = [1,3];  % y = c → work in (x,z)
    case 3, axis_fixed = 3; proj = [1,2];  % z = c → work in (x,y)
    otherwise, error('plane_id must be 1, 2, or 3.');
end

% ---------- 1) Edge-plane intersections (compute once for all edges) ----------
E = size(edges_merged,1);
edge_hit  = false(E,1);
edge_pint = nan(3,E);     % intersection point if hit
edge_inpl = false(E,1);   % edge lies fully in plane → ignore

for e = 1:E
    a = edges_merged(e,1); b = edges_merged(e,2);
    p = coord(:,a); q = coord(:,b);
    [hit, P, in_plane] = edge_axis_plane_intersection(p, q, axis_fixed, plane_val, plane_tol);
    edge_hit(e)  = hit;
    edge_pint(:,e) = P;
    edge_inpl(e) = in_plane;
end

% ---------- 2) Per-face: keep those with exactly two intersecting edges ----------
segP = [];  % 3 x (2*M segments) stacked as [P1 P2 P1 P2 ...]
M = 0;

for f = 1:numel(faces.edges)
    eids = faces.edges{f}(:);
    if isempty(eids), continue; end
    % consider only edges that intersect but are not fully in plane
    mask = edge_hit(eids) & ~edge_inpl(eids);
    if nnz(mask) == 2
        ids2 = eids(mask);
        P1 = edge_pint(:, ids2(1));
        P2 = edge_pint(:, ids2(2));
        if norm(P1 - P2) > snap_tol
            segP = [segP, P1, P2]; %#ok<AGROW>
            M = M + 1;
        end
    end
end

if M == 0
    poly3 = zeros(3,0);
    loops3 = {};
    if do_plot, warning('No intersection segments found.'); end
    return;
end

% reshape segments: 3 x 2 x M
segP = reshape(segP, [3, 2, M]);

% ---------- 3) Connect segments into polygon loop(s) in 2D (projected) ----------
% snap points to merge duplicates
pts2 = squeeze(segP(proj,:,:));   % 2 x 2 x M
pts2 = reshape(pts2, 2, 2*M);     % 2 x (2M)
pts3 = reshape(segP, 3, 2*M);     % 3 x (2M)

% snapping by rounding to grid
grid2 = round(pts2 / snap_tol) * snap_tol;

% unique nodes in 2D
[uniq2, ~, ic] = unique(grid2.', 'rows', 'stable'); % K x 2
ic = reshape(ic, 2, M); % node ids per segment endpoints

K = size(uniq2,1);
% adjacency list of nodes via segments
adj = cell(K,1);
seg_at_node = cell(K,1);
for s = 1:M
    i = ic(1,s); j = ic(2,s);
    if i == j, continue; end
    adj{i}(end+1) = j; seg_at_node{i}(end+1) = s; %#ok<AGROW>
    adj{j}(end+1) = i; seg_at_node{j}(end+1) = s; %#ok<AGROW>
end

% walk loops
usedSeg = false(1,M);
loops3 = {};
for s = 1:M
    if usedSeg(s), continue; end
    % start loop with segment s
    i = ic(1,s); j = ic(2,s);
    loop_nodes = [i, j];
    usedSeg(s) = true;

    % extend forward from j
    cur = j;
    prev = i;
    while true
        [nextNode, segIdx] = next_unvisited(cur, prev, adj, seg_at_node, usedSeg);
        if isempty(nextNode), break; end
        loop_nodes(end+1) = nextNode; %#ok<AGROW>
        usedSeg(segIdx) = true;
        prev = cur; cur = nextNode;
        if cur == loop_nodes(1) % closed
            break;
        end
    end

    % if not closed, try extend backward from i
    if loop_nodes(end) ~= loop_nodes(1)
        cur = i; prev = j;
        while true
            [nextNode, segIdx] = next_unvisited(cur, prev, adj, seg_at_node, usedSeg);
            if isempty(nextNode), break; end
            loop_nodes = [nextNode, loop_nodes]; %#ok<AGROW>
            usedSeg(segIdx) = true;
            prev = cur; cur = nextNode;
            if cur == loop_nodes(end) % closed (head met tail)
                break;
            end
        end
    end

    % form 3D loop
    loop2 = uniq2(loop_nodes, :).';             % 2 x L
    L = size(loop2,2);
    loop3 = zeros(3, L);
    loop3(proj,:) = loop2;
    loop3(axis_fixed,:) = plane_val;
    % close if needed
    if any(loop3(:,1) ~= loop3(:,end)), loop3(:,end+1) = loop3(:,1); end
    loops3{end+1} = loop3; %#ok<AGROW>
end

% ---------- 4) Select the "main" polygon (largest perimeter) ----------
if isempty(loops3)
    poly3 = zeros(3,0);
else
    perims = cellfun(@(P) sum(vecnorm(diff(P,1,2),2,1)), loops3);
    [~, imax] = max(perims);
    poly3 = loops3{imax};
end

% ---------- optional plot ----------
if do_plot && ~isempty(poly3)
    hold_state = ishold; hold on;
    plot3(poly3(1,:), poly3(2,:), poly3(3,:), '-', 'LineWidth', 2);
    axis equal; grid on;
    if ~hold_state, hold off; end
end

% ---------------- helpers ----------------
function val = get_opt(s, name, def)
    if isfield(s,name) && ~isempty(s.(name)), val = s.(name); else, val = def; end
end

function [hit, P, in_plane] = edge_axis_plane_intersection(p, q, ax, c, tol)
    dp = p(ax) - c; dq = q(ax) - c;
    % edge fully in plane?
    if abs(dp) <= tol && abs(dq) <= tol
        hit = false; P = [NaN;NaN;NaN]; in_plane = true; return;
    end
    in_plane = false;
    % straddle or endpoint on plane
    if (dp <= 0 && dq >= 0) || (dp >= 0 && dq <= 0)
        denom = (dp - dq);
        if abs(denom) < eps
            hit = false; P = [NaN;NaN;NaN]; return;
        end
        t = dp / (dp - dq);  % p + t*(q-p)
        if t < -1e-12 || t > 1+1e-12
            hit = false; P = [NaN;NaN;NaN]; return;
        end
        P = p + t*(q - p);
        P(ax) = c; % clamp exact
        hit = true; return;
    else
        hit = false; P = [NaN;NaN;NaN];
    end
end

function [nxt, segIdx] = next_unvisited(cur, prev, adj, seg_at_node, usedSeg)
    % pick a neighbor via an unused segment, not going back to 'prev' unless forced
    nbrs = adj{cur}; sids = seg_at_node{cur};
    % prefer not to go back to prev
    order = [find(nbrs ~= prev), find(nbrs == prev)];
    nxt = [];
    segIdx = [];
    for k = order(:)'
        if ~usedSeg(sids(k))
            nxt = nbrs(k);
            segIdx = sids(k);
            return;
        end
    end
end

end
