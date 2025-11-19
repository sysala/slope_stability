function out = plot_deviatoric_slice(coord, elem, norm_E, transformed_points, axDim, cutCoord, varargin)
%PLOT_DEVIATORIC_SLICE  Exact 2-D slice of a 3-D tet mesh at x/y/z = const, with nodal 3-D interpolation.
%
% Usage (with nodal values you already have):
%   out = plot_deviatoric_slice(coord, elem, norm_E, transformed_points, 2, 30, ...
%                               'NodeValues', node_E, 'GridSpacing', 0.5);
%
% Usage (without nodal values: will derive node_E from GP data once):
%   out = plot_deviatoric_slice(coord, elem, norm_E, transformed_points, 2, 30, ...
%                               'GridSpacing', 0.5);
%
% Inputs
%   coord               3 x Np   node coordinates (P2 allowed; geometry uses elem(1:4,:))
%   elem                10 x Nt  tet connectivity (rows 1..4 are vertex ids of the corners)
%   norm_E              1 x Ni   value at integration points (GPs)
%   transformed_points  3 x Ni   GP positions in global coords
%   axDim               {1|2|3}  plane normal axis (1:x, 2:y, 3:z)
%   cutCoord            scalar   plane coordinate (e.g. y=30 -> axDim=2)
%
% Name-Value options
%   'NodeValues'        Np x 1 or 1 x Np vector of nodal values to use (preferred)
%   'GridSpacing'       target spacing for interior sampling (default: longSide/300)
%   'Tol'               geometric snap tolerance (default: 1e-8 * mesh scale)
%
% Returns
%   out.Ploops    cell array of loops (Mx2) describing slice boundary (no polyshape used)
%   out.V, out.T  triangulation vertices (2D) and faces
%   out.C         interpolated values at V
%   out.node_E    nodal values actually used (either provided or derived from GP)
%   out.edges     constraint edges on V

% ---------- sanity & defaults ----------
coord = double(coord); elem = double(elem);
norm_E = double(norm_E(:)); transformed_points = double(transformed_points);
assert(size(coord,1)==3, 'coord must be 3 x Np');
assert(size(elem,1)>=4, 'elem must have at least 4 rows (tet corners in 1..4)');
assert(ismember(axDim,[1 2 3]), 'axDim must be 1, 2, or 3');
Ni = size(transformed_points,2); assert(numel(norm_E)==Ni, 'GP sizes mismatch');

Np = size(coord,2);
ext = max(coord,[],2) - min(coord,[],2);
scale = max(ext); if scale==0, scale=1; end
bbox    = [min(coord,[],2) max(coord,[],2)];
keep2D = [1 2; 1 3; 2 3]; keep = keep2D(axDim,:);
labs = {{'y','z'},{'x','z'},{'x','y'}}; lab = labs{axDim};
planeName = sprintf('%s = %.6g', 'xyz',(axDim), cutCoord);

% defaults
GridSpacing = [];
Tol = 1e-8*scale;
NodeValues = [];

% parse name-values
if ~isempty(varargin)
    for k=1:2:numel(varargin)
        switch lower(varargin{k})
            case 'gridspacing', GridSpacing = varargin{k+1};
            case 'tol',         Tol = varargin{k+1};
            case 'nodevalues',  NodeValues = varargin{k+1};
            otherwise, error('Unknown option: %s', varargin{k});
        end
    end
end

% sensible default spacing
if isempty(GridSpacing)
    longSide = max(bbox(keep,2) - bbox(keep,1));
    GridSpacing = max(longSide/300, eps);
end
snapTol = max(Tol, 1e-12*scale);

% ---------- 1) exact slice boundary: surface tris ∩ plane ----------
% surface faces
tetV = elem(1:4,:).';
F = [tetV(:,[1 2 3]); tetV(:,[1 2 4]); tetV(:,[1 3 4]); tetV(:,[2 3 4])];
Fsort = sort(F,2); [~,~,ic] = unique(Fsort,'rows'); counts = accumarray(ic,1);
surfFaces = F(counts(ic)==1,:);

% intersect each surface triangle with the plane, collect 2D segments
segs = zeros(0,4);
for i = 1:size(surfFaces,1)
    tri = coord(:, surfFaces(i,:));     % 3x3
    seg3 = local_tri_plane_intersection(tri, axDim, cutCoord, snapTol);
    if ~isempty(seg3)
        p2 = seg3(keep,:).';            % 2x2 rows: [x y]
        if norm(p2(1,:)-p2(2,:)) > 100*eps(scale)
            segs(end+1,:) = [p2(1,:) p2(2,:)]; %#ok<AGROW>
        end
    end
end
if isempty(segs)
    warning('Plane does not cut the mesh surface. Nothing to plot.');
    out = local_empty_out(Ni);
    return
end

% snap endpoints and make unique vertex list
E = [segs(:,1:2); segs(:,3:4)];
Esn = round(E/snapTol)*snapTol;
[U, ~, im] = unique(Esn, 'rows', 'stable');   % unique snapped vertices
edges = [im(1:size(segs,1)) im(size(segs,1)+1:end)];  % edges as indices into U
% remove degenerate/duplicate edges
edges = edges(edges(:,1)~=edges(:,2),:);
edges = unique(sort(edges,2), 'rows');

% split edges at any coincident vertices lying on them
edges = local_split_edges(U, edges, snapTol);

% stitch into closed loops (degree-2 graph), pure combinatorial
loopsIdx = local_build_loops_degree2(U, edges);
% to coordinates
Ploops = cellfun(@(idx) U(idx,:), loopsIdx, 'UniformOutput', false);
% remove tiny loops
Ploops = Ploops( cellfun(@(L) (size(L,1)>=3) && local_polyarea(L)>(snapTol^2), Ploops) );
if isempty(Ploops)
    warning('Failed to assemble a valid slice polygon from intersections.');
    out = local_empty_out(Ni); return
end

% ---------- 2) constrained triangulation (no polyshape) ----------
% (a) resample each loop to ~GridSpacing for better triangle quality
PloopsRS = cell(size(Ploops));
for i=1:numel(Ploops)
    PloopsRS{i} = local_resample_loop(Ploops{i}, GridSpacing);
end

% (b) collect boundary vertices & build constraint edges
bndPts = []; loopSizes = zeros(numel(PloopsRS),1);
for i=1:numel(PloopsRS)
    L = PloopsRS{i};
    if isempty(L), continue; end
    if norm(L(1,:)-L(end,:)) < 10*snapTol, L = L(1:end-1,:); end
    if size(L,1) < 3, continue; end
    bndPts = [bndPts; L]; %#ok<AGROW>
    loopSizes(i) = size(L,1);
end
if isempty(bndPts)
    out = local_empty_out(Ni); return
end

% (c) interior grid points and point-in-polygon via even-odd rule on loops
xmin = min(bndPts(:,1)); xmax = max(bndPts(:,1));
ymin = min(bndPts(:,2)); ymax = max(bndPts(:,2));
xv = xmin:GridSpacing:xmax; if numel(xv)<2, xv = linspace(xmin,xmax,25); end
yv = ymin:GridSpacing:ymax; if numel(yv)<2, yv = linspace(ymin,ymax,25); end
[Xg,Yg] = meshgrid(xv,yv);
cand = [Xg(:) Yg(:)];
in = local_points_in_loops(cand, PloopsRS);   % even-odd rule across all loops
pin = cand(in,:);

% (d) assemble vertices, snap-unique, and remap constraint edges
Vall = [bndPts; pin];
Vsn = round(Vall/snapTol)*snapTol;
[V, ~, imap] = unique(Vsn, 'rows', 'stable');

Cedges = zeros(0,2);
pos = 0;
for i=1:numel(PloopsRS)
    m = loopSizes(i);
    if m<3, continue; end
    idx = pos + (1:m);
    ui = imap(idx).';
    ui = local_unique_consecutive(ui);
    if numel(ui)>=3
        Cedges = [Cedges; [ui(:) [ui(2:end) ui(1)].']]; %#ok<AGROW>
    end
    pos = pos + m;
end
% remove duplicate edges
if ~isempty(Cedges)
    Es = sort(Cedges,2);
    [~,ia] = unique(Es,'rows','stable');
    Cedges = Cedges(ia,:);
end

% (e) constrained Delaunay + keep only triangles whose centroids are inside (even-odd)
dt = delaunayTriangulation(V, Cedges);
T = dt.ConnectivityList;
cent = (V(T(:,1),:) + V(T(:,2),:) + V(T(:,3),:))/3;
inside = local_points_in_loops(cent, PloopsRS);
T = T(inside,:);
% prune degenerate slivers
if ~isempty(T)
    A = local_tri_area_2d(V,T);
    T = T(A > (GridSpacing^2)*1e-4, :);
end
if isempty(T)
    warning('Triangulation produced no interior faces.'); 
    out = struct('Ploops',{Ploops},'V',V,'T',[],'C',[],'node_E',[],'edges',Cedges);
    ax=gca; hold(ax,'on'); local_plot_loops(Ploops, lab, planeName); return
end

% ---------- 3) 3-D interpolation from nodal values ----------
if ~isempty(NodeValues)
    node_E = double(NodeValues(:));
    assert(numel(node_E)==Np, 'NodeValues length must match number of nodes (size(coord,2)).');
else
    % Derive nodal values once from GP data in 3-D
    F3 = scatteredInterpolant( ...
        transformed_points(1,:).', transformed_points(2,:).', transformed_points(3,:).', ...
        norm_E, 'natural', 'nearest');
    node_E = F3(coord(1,:).', coord(2,:).', coord(3,:).');
end

% Build 3-D coordinates for each 2-D vertex on the slice plane
V3 = zeros(size(V,1),3);
switch axDim
    case 1, V3(:,1) = cutCoord; V3(:,2:3) = V;          % [x, y, z] with x = const
    case 2, V3(:,2) = cutCoord; V3(:,[1 3]) = V;        % y = const
    case 3, V3(:,3) = cutCoord; V3(:,1:2) = V;          % z = const
end

% Interpolate from nodal values to slice vertices (3-D scatteredInterpolant over nodes)
Fnode = scatteredInterpolant(coord(1,:).', coord(2,:).', coord(3,:).', node_E, 'natural', 'nearest');
Cvals = Fnode(V3(:,1), V3(:,2), V3(:,3));

% ---------- 4) plot ----------
ax = gca; hold(ax,'on');
patch('Faces',T,'Vertices',V,'FaceVertexCData',Cvals, ...
      'FaceColor','interp','EdgeColor','none','Parent',ax);
local_plot_loops(Ploops, lab, planeName);
view(ax,2); axis(ax,'equal'); axis(ax,'tight'); colormap(ax,parula); colorbar(ax);
xlabel(ax, lab{1}); ylabel(ax, lab{2});
title(ax, sprintf('Value on slice (%s) — nodal 3D interpolation', planeName));

% output
out = struct('Ploops',{Ploops}, 'V',V, 'T',T, 'C',Cvals, 'node_E',node_E, 'edges',Cedges);
end

% ======================= helpers =======================

function seg = local_tri_plane_intersection(tri, axDim, c, tol)
% tri: 3x3 (xyz rows), returns 3x2 endpoints in 3D or []
d = tri(axDim,:) - c; d(abs(d)<tol)=0;
if all(d>0)||all(d<0)||all(d==0), seg=[]; return; end
edges = [1 2; 2 3; 3 1];
pts = [];
for e=1:3
    i = edges(e,1); j = edges(e,2);
    di=d(i); dj=d(j); pi=tri(:,i); pj=tri(:,j);
    if di==0 && dj==0
        continue
    elseif di==0
        pts(:,end+1)=pi; %#ok<AGROW>
    elseif dj==0
        pts(:,end+1)=pj; %#ok<AGROW>
    elseif di*dj<0
        t=di/(di-dj); pts(:,end+1)=pi+t*(pj-pi); %#ok<AGROW>
    end
end
if isempty(pts), seg=[]; return; end
pts2 = uniquetol(pts.', tol, 'ByRows', true).';
switch size(pts2,2)
    case 0, seg=[];
    case 1, seg=[];
    case 2, seg=pts2;
    otherwise
        % pick farthest pair (small k, avoid pdist)
        k=size(pts2,2); best=[1 2]; bestd=-inf;
        for a=1:k-1, for b=a+1:k
            d2=sum((pts2(:,a)-pts2(:,b)).^2);
            if d2>bestd, bestd=d2; best=[a b]; end
        end, end
        seg=pts2(:,best);
end
end

function edges2 = local_split_edges(U, edges, tol)
% Split any edge [a b] at snapped vertices that lie on its segment.
% Ensures no "coincident point splits" occur later.
nU = size(U,1);
edges2 = zeros(0,2);
for k=1:size(edges,1)
    a = edges(k,1); b = edges(k,2);
    pa = U(a,:); pb = U(b,:);
    e = pb - pa; len2 = sum(e.^2);
    if len2==0, continue; end
    % bounding box filter
    mn = min(pa,pb)-tol; mx = max(pa,pb)+tol;
    cand = find(U(:,1)>=mn(1) & U(:,1)<=mx(1) & U(:,2)>=mn(2) & U(:,2)<=mx(2));
    ts = [0; 1]; idxs = [a; b];
    for j = cand.'
        if j==a || j==b, continue; end
        v = U(j,:) - pa;
        % collinearity by cross product magnitude
        cross = abs(e(1)*v(2) - e(2)*v(1));
        if cross <= tol*max(1,sqrt(len2))
            t = (e(1)*v(1) + e(2)*v(2))/len2;
            if t>0 && t<1
                ts(end+1,1)=t; idxs(end+1,1)=j; %#ok<AGROW>
            end
        end
    end
    % sort and create sub-edges
    [ts,ord] = sort(ts); idxs = idxs(ord);
    for s=1:numel(ts)-1
        i1 = idxs(s); i2 = idxs(s+1);
        if i1~=i2
            edges2(end+1,:) = [i1 i2]; %#ok<AGROW>
        end
    end
end
% cleanup duplicates and orientation
edges2 = unique(sort(edges2,2),'rows');
end

function loops = local_build_loops_degree2(U, edges)
% Build cycles in a graph where every vertex has degree 0 or 2 (after splits).
adj = cell(size(U,1),1);
for k=1:size(edges,1)
    a=edges(k,1); b=edges(k,2);
    if a==b, continue; end
    adj{a}(end+1)=b; %#ok<AGROW>
    adj{b}(end+1)=a; %#ok<AGROW>
end
used = false(size(edges,1),1);
edgeAt = containers.Map('KeyType','char','ValueType','int32');
for k=1:size(edges,1)
    key = sprintf('%d-%d', min(edges(k,:)), max(edges(k,:)));
    edgeAt(key) = k;
end
loops = {};
visitedV = false(size(U,1),1);
for v=1:size(U,1)
    if visitedV(v), continue; end
    if numel(adj{v})~=2, visitedV(v)=true; continue; end
    % start a loop at v
    path = v;
    prev = -1; curr = v;
    while true
        nbrs = adj{curr};
        % choose next neighbor not equal to prev
        nxt = nbrs(1); if nxt==prev && numel(nbrs)>=2, nxt = nbrs(2); end
        path(end+1) = nxt; %#ok<AGROW>
        prev = curr; curr = nxt;
        if curr==v, break; end
        if numel(adj{curr})~=2, break; end
        if numel(path) > 1e6, break; end % safety
    end
    % close & clean
    if numel(path)>=4 && path(end)==path(1)
        path(end) = [];
        path = local_unique_consecutive(path);
        if numel(path)>=3
            loops{end+1} = path; %#ok<AGROW>
            visitedV(path) = true;
        end
    end
end
end

function area = local_polyarea(L)
% signed area (2D)
x=L(:,1); y=L(:,2);
area = 0.5 * sum(x.*y([2:end 1]) - y.*x([2:end 1]));
area = abs(area);
end

function Ld = local_resample_loop(L, step)
if isempty(L), Ld=L; return; end
if norm(L(1,:)-L(end,:)) < eps, L = L(1:end-1,:); end
n = size(L,1); Ld = zeros(0,2);
for i=1:n
    p=L(i,:); q=L(mod(i,n)+1,:);
    e = q-p; len = norm(e);
    m = max(1, ceil(len/max(step,eps)));
    if m==1, Ld=[Ld; p]; %#ok<AGROW>
    else
        t = linspace(0,1,m+1).'; t = t(1:end-1);
        Ld=[Ld; p + t.*e]; %#ok<AGROW>
    end
end
end

function inside = local_points_in_loops(P, loops)
% Even-odd rule: inside if contained in an odd number of loops.
inside = false(size(P,1),1);
for i=1:numel(loops)
    L = loops{i};
    if isempty(L), continue; end
    [IN, ON] = inpolygon(P(:,1), P(:,2), L(:,1), L(:,2));
    inside = xor(inside, IN | ON);
end
end

function A = local_tri_area_2d(V,T)
p1=V(T(:,1),:); p2=V(T(:,2),:); p3=V(T(:,3),:);
A = 0.5*abs( (p2(:,1)-p1(:,1)).*(p3(:,2)-p1(:,2)) - (p2(:,2)-p1(:,2)).*(p3(:,1)-p1(:,1)) );
end

function idx2 = local_unique_consecutive(idx)
if isempty(idx), idx2=idx; return; end
mask = [true; diff(idx(:))~=0];
idx2 = idx(mask);
end

function local_plot_loops(Ploops, lab, planeName)
ax = gca;
for i=1:numel(Ploops)
    L = Ploops{i};
    plot(ax, L([1:end 1],1), L([1:end 1],2), 'k-', 'LineWidth', 1.2); hold(ax,'on');
end
xlabel(ax, lab{1}); ylabel(ax, lab{2});
title(ax, sprintf('Slice outline (%s)', planeName));
end

function out = local_empty_out(Ni)
out = struct('Ploops',{{}},'V',[],'T',[],'C',[],'node_E',[],'edges',[],'usedPointsMask',false(1,Ni));
end
