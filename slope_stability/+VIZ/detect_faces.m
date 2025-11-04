function faces = detect_faces(edges_merged, coord, opts)
%GROUP_FACES_BY_EDGE_PAIRS  Identify faces as groups of boundary edges
% using edge-pair plane normals and simple BFS-style grouping.
%
% faces = group_faces_by_edge_pairs(edges_merged, coord, opts)
%
% Inputs
%   edges_merged : (E x 2) int   boundary edges (vertex indices, 1-based)
%   coord        : (3 x N) double vertex coordinates
%   opts (optional):
%       .normal_cos_tol   : accept if abs(dot(n_i, n_seed)) >= 1 - tol (def: 1e-4)
%       .sin_colinear_tol : skip pairs if ||cross(t1,t2)|| <= tol (def: 1e-8)
%
% Output struct
%   faces.edges   : 1xF cell, each is vector of edge indices (into edges_merged)
%   faces.pairs   : 1xF cell, each is vector of pair indices used to form the face
%   faces.normals : (3 x F) representative unit normal per face
%
% Algorithm (as requested):
%   1) list all pairs of edges that share a vertex; compute unit plane normal
%   2) pick first pair as a seed; grow a group by adding pairs that:
%        - share at least one edge with the group, and
%        - have nearly parallel normals (abs dot ~ 1 up to threshold)
%      repeat until no change → one face
%   3) repeat with remaining pairs until none left

if nargin < 3, opts = struct(); end
normal_cos_tol   = get_opt(opts, 'normal_cos_tol',   1e-4);
sin_colinear_tol = get_opt(opts, 'sin_colinear_tol', 1e-8);

E = size(edges_merged,1);
V = size(coord,2);
I = edges_merged(:,1);
J = edges_merged(:,2);

% incident edges for each vertex
incEdges = accumarray([I;J], [(1:E)'; (1:E)'], [V,1], @(x){x}, {[]});

% ---- build all edge pairs that share a vertex + their unit normals ----
pairs_e1 = []; pairs_e2 = []; pairs_v = []; pair_normals = [];

for v = 1:V
    eList = incEdges{v};
    m = numel(eList);
    if m < 2, continue; end
    for a = 1:m-1
        eA = eList(a);
        for b = a+1:m
            eB = eList(b);
            % directions of edges from the common vertex v
            tA = edge_dir_from_vertex(eA, v, I, J, coord);
            tB = edge_dir_from_vertex(eB, v, I, J, coord);
            n  = cross(tA, tB);
            s  = norm(n);
            if s <= sin_colinear_tol
                continue; % edges are (near) colinear → no plane
            end
            n = n / s; % unit normal (sign arbitrary; we use abs dot later)
            pairs_e1(end+1,1) = eA; %#ok<AGROW>
            pairs_e2(end+1,1) = eB; %#ok<AGROW>
            pairs_v(end+1,1)  = v;  %#ok<AGROW>
            pair_normals(:,end+1) = n; %#ok<AGROW>
        end
    end
end

P = numel(pairs_e1);
assigned = false(P,1);

faces.edges   = {};
faces.pairs   = {};
faces.normals = [];

% ---- group pairs into faces ----
while true
    seed = find(~assigned, 1, 'first');
    if isempty(seed), break; end

    n_seed = pair_normals(:,seed);
    group_pairs = seed;
    assigned(seed) = true;

    group_edge_mask = false(E,1);
    group_edge_mask(pairs_e1(seed)) = true;
    group_edge_mask(pairs_e2(seed)) = true;

    changed = true;
    while changed
        changed = false;
        for p = 1:P
            if assigned(p), continue; end
            % must share at least one edge with the current group
            if ~(group_edge_mask(pairs_e1(p)) || group_edge_mask(pairs_e2(p)))
                continue;
            end
            % normals must be (nearly) parallel
            if abs(dot(pair_normals(:,p), n_seed)) >= 1 - normal_cos_tol
                assigned(p) = true;
                group_pairs(end+1) = p; %#ok<AGROW>
                group_edge_mask(pairs_e1(p)) = true;
                group_edge_mask(pairs_e2(p)) = true;
                changed = true;
            end
        end
    end

    % edges belonging to this face = all unique edges present in its pairs
    face_edges = unique([pairs_e1(group_pairs); pairs_e2(group_pairs)])';
    faces.edges{end+1}   = face_edges; %#ok<AGROW>
    faces.pairs{end+1}   = group_pairs; %#ok<AGROW>
    faces.normals(:,end+1) = n_seed;   %#ok<AGROW>
end

% ------------------- helpers -------------------
function val = get_opt(s, name, def)
    if isfield(s, name) && ~isempty(s.(name)), val = s.(name); else, val = def; end
end

function t = edge_dir_from_vertex(eidx, vtx, I, J, coord)
    a = I(eidx); b = J(eidx);
    if vtx == a
        t = coord(:,b) - coord(:,a);
    elseif vtx == b
        t = coord(:,a) - coord(:,b);
    else
        error('Internal: edge %d does not touch vertex %d.', eidx, vtx);
    end
    n = norm(t);
    if n == 0, error('Zero-length edge encountered.'); end
    t = t / n;
end
end
