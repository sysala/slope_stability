function edges_out = merge_collinear_edges(edges_in, coord, angle_tol_deg, min_edge_len, max_passes)
%MERGE_COLLINEAR_EDGES Merge consecutive edges that are collinear at a shared node.
%
% edges_in       : [M x 2] int, undirected edges (node indices).
% coord          : [3 x N] double, node coordinates.
% angle_tol_deg  : (opt) scalar, tolerance around 180° at the shared node. Default 5.
% min_edge_len   : (opt) scalar, ignore edges shorter than this (units of coord). Default auto.
% max_passes     : (opt) int, max global passes. Default 5*M.
%
% Returns:
%   edges_out    : [K x 2] unique, undirected merged edges.
%
% Notes:
% - Merges only at nodes of degree 2 (avoids T-junctions).
% - Tests collinearity using vectors from the shared node v to its neighbors:
%       t1 = (a - v), t2 = (b - v)  -> must be nearly opposite (dot ~ -1).
% - Robust in 3D (dot + cross checks).
%
% Example:
%   E2 = merge_collinear_edges(E, coord, 3, []);  % 3° tolerance

    if nargin < 3 || isempty(angle_tol_deg), angle_tol_deg = 5; end
    if nargin < 4, min_edge_len = []; end

    % Clean & normalize input edge list
    E = sort(edges_in, 2);
    E = unique(E, 'rows');

    N = size(coord, 2);
    assert(size(coord,1) == 3, 'coord must be 3xN');

    % Auto min length (based on bbox) if not provided
    if isempty(min_edge_len)
        bbox_diag = norm(max(coord,[],2) - min(coord,[],2));
        min_edge_len = max(1e-12, 1e-12 * bbox_diag);
    end

    % Precompute tolerances
    cos_tol = cosd(angle_tol_deg);   % want dot(t1, t2) ~ -1 -> dot(t1, t2) <= -cos_tol
    sin_tol = sind(angle_tol_deg);   % small cross magnitude

    if nargin < 5 || isempty(max_passes), max_passes = 5 * size(E,1); end
    pass = 0;

    while pass < max_passes
        pass = pass + 1;
        merged_any = false;

        % Degree per node
        deg = accumarray(E(:), 1, [N, 1], @sum, 0);

        % Candidates: degree exactly 2
        cand = find(deg == 2);
        if isempty(cand), break; end

        % We'll mark edges to drop via a logical mask
        drop = false(size(E,1),1);
        adds = [];  % new merged edges to append

        % Build quick incidence map: edge i connects E(i,1) -- E(i,2)
        % For lookups, we’ll just scan; M ~ few 100, OK. (Vectorized incidence is overkill here.)
        for idx = 1:numel(cand)
            v = cand(idx);

            % Find the two incident edges (current E and not already dropped)
            inc_mask = ~drop & (E(:,1) == v | E(:,2) == v);
            if nnz(inc_mask) ~= 2
                continue; % degree may have changed this pass
            end
            inc_edges = E(inc_mask, :);

            % Neighbors a, b (the other endpoint in each edge)
            ab = [ ...
                inc_edges(1,1) + inc_edges(1,2) - v; ...
                inc_edges(2,1) + inc_edges(2,2) - v ];
            a = ab(1); b = ab(2);
            if a == v || b == v || a == b, continue; end

            Va = coord(:,a) - coord(:,v);  na = norm(Va);
            Vb = coord(:,b) - coord(:,v);  nb = norm(Vb);

            if na < min_edge_len || nb < min_edge_len
                continue; % too small to trust
            end

            t1 = Va / na;  % tangent from v to a
            t2 = Vb / nb;  % tangent from v to b

            % Collinearity at v: t1 and t2 should be opposite (angle ~ 180°).
            c = cross(t1, t2);
            dotv = dot(t1, t2);
            is_straight = (norm(c) <= sin_tol) && (dotv <= -cos_tol);

            if ~is_straight
                continue;
            end

            % Merge: remove (a–v) and (v–b), add (a–b)
            rows = find(inc_mask);
            drop(rows) = true;

            new_edge = sort([a, b]);
            adds = [adds; new_edge]; %#ok<AGROW>
            merged_any = true;
        end

        if ~merged_any
            break;
        end

        % Commit: drop old, add new, then uniquify
        E = [E(~drop, :); adds];
        E = sort(E, 2);
        E = unique(E, 'rows');
    end

    if pass == max_passes
        warning('merge_collinear_edges:MaxPasses','Reached max passes; stopping.');
    end

    edges_out = E;
end
