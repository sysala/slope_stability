function [mergedEdges, out] = find_outside_edges_p2_tet(coord, elem, plane_tol_deg, colinear_tol_deg, do_plot)
% FIND_OUTSIDE_EDGES_P2_TET
% Robustly find & plot outside (crease) edges of a 3-D P2 tetra mesh.
%
% INPUTS
%   coord : (3 x Np) node coordinates (includes midside nodes)
%   elem  : (10 x Nt) connectivity of P2 tets; ONLY elem(1:4,:) used here
%
% OPTIONAL
%   plane_tol_deg     : faces with normal angle <= this are considered coplanar (default 3)
%   colinear_tol_deg  : merge segments whose directions differ by <= this (default 5)
%   do_plot           : true/false to plot merged edges (default true)
%
% OUTPUTS
%   mergedEdges : (K x 2) vertex indices of merged outside edges (start,end)
%   out         : struct with diagnostics (see fields below)
%
% Steps:
%   1) Build all faces from the first 4 vertices of each P2 tetra.
%   2) Boundary faces = faces with incidence to exactly 1 tetra (via incidence matrix).
%   3) Compute **unit outward normals** (vectorized) for boundary faces; drop degenerate tiny faces.
%   4) Build boundary edge→faces adjacency.
%   5) Build a graph connecting **coplanar** neighboring boundary faces (angle ≤ plane_tol_deg).
%      Coplanar clusters are connected components of this graph.
%   6) Mark crease edges: boundary edges whose two incident faces belong to **different clusters**
%      (or edges with only one incident face).
%   7) Merge nearly colinear crease segments into long outside edges.
%
% Michal-ready: save as `find_outside_edges_p2_tet.m` and call directly.

    if nargin < 3 || isempty(plane_tol_deg),    plane_tol_deg    = 3;   end
    if nargin < 4 || isempty(colinear_tol_deg), colinear_tol_deg = 5;   end
    if nargin < 5 || isempty(do_plot),          do_plot          = true; end

    % ---------- basics ----------
    V = elem(1:4, :).';      % Nt x 4
    Nt = size(V,1);
    Np = size(coord,2); %#ok<NASGU>

    % ---------- all oriented faces for each tet ----------
    f1 = V(:, [2 3 4]);  opp1 = V(:,1);
    f2 = V(:, [1 4 3]);  opp2 = V(:,2);
    f3 = V(:, [1 2 4]);  opp3 = V(:,3);
    f4 = V(:, [1 3 2]);  opp4 = V(:,4);

    faces_oriented = [f1; f2; f3; f4];       % (4Nt x 3)
    opp_node       = [opp1; opp2; opp3; opp4];
    face2tet       = repmat((1:Nt).', 4, 1); % (4Nt x 1) maps row->tet

    faces_sorted = sort(faces_oriented,2);
    [faces_unique, ia, ic] = unique(faces_sorted, 'rows', 'stable'); % F unique faces
    F = size(faces_unique,1);

    % ---------- face–tet incidence & boundary faces ----------
    S_face_tet = sparse(ic, face2tet, 1, F, Nt);
    deg_face   = full(S_face_tet * ones(Nt,1));     % times a face appears among tets
    bnd_uid    = find(deg_face == 1);               % unique ids of boundary faces

    raw_bnd_idx     = ia(bnd_uid);                  % raw row index for an oriented representative
    bnd_faces_orient= faces_oriented(raw_bnd_idx,:);% oriented triplets (for normals)
    Fb = size(bnd_faces_orient,1);

    % ---------- vectorized outward unit normals ----------
    A   = coord(:, bnd_faces_orient(:,1)).';
    B   = coord(:, bnd_faces_orient(:,2)).';
    C   = coord(:, bnd_faces_orient(:,3)).';
    P0  = coord(:, opp_node(raw_bnd_idx)).';
    cent= (A + B + C)/3;

    N   = cross(B - A, C - A, 2);           % raw normals
    toOpp = P0 - cent;
    s   = sign(sum(N .* toOpp, 2)); s(s==0)=1;
    N   = -N .* s;                          % flip to point outward w.r.t the tetra
    Atri= 0.5*vecnorm(N,2,2);               % triangle areas
    Nlen= vecnorm(N,2,2);  nz = Nlen>0;
    Nunit = zeros(Fb,3); Nunit(nz,:) = N(nz,:)./Nlen(nz);

    % drop degenerate tiny faces (numerical slivers)
    bbox = [min(coord,[],2) max(coord,[],2)];
    diag2 = sum((bbox(:,2)-bbox(:,1)).^2);
    area_tol = max(1e-14*diag2, eps);
    keep = Atri > area_tol;
    bnd_faces_orient = bnd_faces_orient(keep,:);
    bnd_uid          = bnd_uid(keep);
    Nunit            = Nunit(keep,:);
    Fb = size(bnd_faces_orient,1);
    cent = cent(keep,:);

    % ---------- build boundary edges & edge→faces ----------
    Eraw = [ bnd_faces_orient(:,[1 2]);
             bnd_faces_orient(:,[2 3]);
             bnd_faces_orient(:,[3 1]) ];
    Eraw_sorted = sort(Eraw,2);
    face_ids_rep = repelem((1:Fb).', 3, 1);

    [edges_unique, iae, ice] = unique(Eraw_sorted, 'rows', 'stable'); %#ok<ASGLU>
    Me = size(edges_unique,1);
    % Edge->faces incidence on boundary-only edges
    edge2faces = accumarray(ice, face_ids_rep, [Me,1], @(x){x});

    % ---------- build coplanar-face adjacency graph ----------
    % Two boundary faces are considered "coplanar-neighbors" if:
    %   (i) they share a boundary edge; AND
    %  (ii) angle between their unit normals <= plane_tol_deg
    cos_plane = cosd(plane_tol_deg);
    ii = []; jj = [];
    for e = 1:Me
        fids = edge2faces{e};
        if numel(fids) ~= 2, continue; end
        f1 = fids(1); f2 = fids(2);
        c  = dot(Nunit(f1,:), Nunit(f2,:));
        c  = max(-1,min(1,c));
        if c >= cos_plane
            ii(end+1) = f1; %#ok<AGROW>
            jj(end+1) = f2; %#ok<AGROW>
        end
    end
    Apl = sparse([ii jj], [jj ii], 1, Fb, Fb);  % symmetric
    G   = graph(Apl);
    comp= conncomp(G);                          % coplanar clusters (1..K)

    % ---------- mark crease edges (different clusters or single-face edges) ----------
    sharp_mask = false(Me,1);
    for e = 1:Me
        fids = edge2faces{e};
        if isempty(fids)
            continue
        elseif numel(fids)==1
            sharp_mask(e) = true; % silhouette / open boundary
        else
            sharp_mask(e) = (comp(fids(1)) ~= comp(fids(2)));
        end
    end
    sharp_edges = edges_unique(sharp_mask,:);
    
    % ---------- merge colinear crease segments ----------
    mergedEdges = merge_colinear(sharp_edges, coord, colinear_tol_deg);

    % ---------- outputs ----------
    out.faces_unique          = faces_unique;
    out.boundary_face_ids     = bnd_uid;
    out.boundary_faces_sorted = sort(bnd_faces_orient,2);
    out.boundary_normals_unit = Nunit;
    out.boundary_centroids    = cent;
    out.edge2faces            = edge2faces;
    out.sharp_edges_unmerged  = sharp_edges;
    out.params = struct('plane_tol_deg', plane_tol_deg, ...
                        'colinear_tol_deg', colinear_tol_deg, ...
                        'area_tol', area_tol);

    % ---------- plot ----------
    if do_plot
        figure('Color','w'); hold on; axis equal vis3d
        xlabel('x'); ylabel('y'); zlabel('z');
        title(sprintf('Outside edges (plane tol %g°, colinear tol %g°)', ...
                      plane_tol_deg, colinear_tol_deg));
        for k = 1:size(mergedEdges,1)
            i = mergedEdges(k,1); j = mergedEdges(k,2);
            p = coord(:,[i j]);
            plot3(p(1,:),p(2,:),p(3,:),'-','LineWidth',2);
        end
        grid on; view(3);
    end
end

% ======================= helpers =======================
function merged = merge_colinear(edgePairs, coord, colinear_tol_deg)
    % Merge chains of nearly colinear segments into long edges.
    if isempty(edgePairs)
        merged = zeros(0,2); return
    end

    % Undirected edge graph on vertices
    nNodes = size(coord,2);
    M = size(edgePairs,1);
    deg = accumarray([edgePairs(:,1); edgePairs(:,2)], 1, [nNodes,1]);
    incident = accumarray([edgePairs(:,1); edgePairs(:,2)], ...
                          [1:M  1:M].', [nNodes,1], @(x){x});

    used = false(M,1);
    merged = zeros(0,2);

    cos_tol = cosd(colinear_tol_deg);

    % function to get directed unit vector from a->b
    function u = dir_ab(a,b)
        v = (coord(:,b) - coord(:,a)).'; L = norm(v);
        if L==0, u = [1 0 0]; else, u = v/L; end
    end

    % grow a path from (v, prev_dir), preferring colinear continuation
    function [vend, used_here] = grow(v, prev_dir)
        used_here = [];
        while true
            candE = incident{v};
            candE = candE(~used(candE));
            if isempty(candE), break; end

            % among candidate edges touching v, pick the most colinear
            bestE = 0; bestCos = -Inf; v_next = v; dir_next = prev_dir;
            for ee = candE(:).'
                a = edgePairs(ee,1); b = edgePairs(ee,2);
                if v == a
                    u = dir_ab(a,b); v2 = b;
                elseif v == b
                    u = dir_ab(b,a); v2 = a;
                else
                    continue
                end
                c = dot(u, prev_dir);
                if c > bestCos
                    bestCos = c; bestE = ee; v_next = v2; dir_next = u;
                end
            end

            if bestCos < cos_tol
                break
            end

            used(bestE) = true;
            used_here(end+1) = bestE; %#ok<AGROW>
            v = v_next; prev_dir = dir_next;
        end
        vend = v;
    end

    % Start at junctions and endpoints (deg ~= 2)
    endpoints = find(deg~=2);
    started = false(M,1);

    % Walk from endpoints first (good for polyline chains)
    for v0 = endpoints(:).'
        candE = incident{v0};
        candE = candE(~used(candE));
        for e0 = candE(:).'
            if used(e0), continue; end
            a = edgePairs(e0,1); b = edgePairs(e0,2);
            if v0 == a, u0 = dir_ab(a,b); s = a; t = b;
            else,        u0 = dir_ab(b,a); s = b; t = a;
            end
            used(e0) = true; started(e0) = true;
            [vL, usedL] = grow(s, -u0);
            [vR, usedR] = grow(t,  u0);
            used([usedL usedR]) = true;
            if vL ~= vR
                merged(end+1,:) = sort([vL vR]); %#ok<AGROW>
            else
                merged(end+1,:) = sort([s t]);   %#ok<AGROW>
            end
        end
    end

    % Handle residual loops (all deg==2), start anywhere
    leftover = find(~used);
    while ~isempty(leftover)
        e0 = leftover(1);
        a = edgePairs(e0,1); b = edgePairs(e0,2);
        u0 = dir_ab(a,b);
        used(e0) = true;
        [vL, usedL] = grow(a, -u0);
        [vR, usedR] = grow(b,  u0);
        used([usedL usedR]) = true;
        if vL ~= vR
            merged(end+1,:) = sort([vL vR]); %#ok<AGROW>
        else
            merged(end+1,:) = sort([a b]);   %#ok<AGROW>
        end
        leftover = find(~used);
    end

    % unique (unordered)
    merged = unique(sort(merged,2),'rows');
end
