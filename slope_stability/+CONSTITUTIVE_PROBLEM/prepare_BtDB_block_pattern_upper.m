function P = prepare_BtDB_block_pattern_upper(B, n_strain, n_int, WEIGHT)
%--------------------------------------------------------------------------
% prepare_BtDB_block_pattern_upper  Precompute the upper-triangle assembly
% pattern for the product  K = B' * D * B  where D is block-diagonal with
% n_strain x n_strain blocks and B is fixed.
%
% This is a generalisation of the diagonal-D prebuilt pattern described in
% SPARSERSB_PREBUILT_GUIDE.md.  The block-diagonal D has n_int blocks of
% size n_strain x n_strain.  Each block is characterised by n_strain^2
% values stored column-major in a vector d of length n_strain^2 * n_int,
% laid out as  d = DS(:)  where DS is  (n_strain^2) x n_int.
%
% Usage (one-time precompute):
%   P = prepare_BtDB_block_pattern_upper(B, n_strain, n_int, WEIGHT);
%
% Hot-path assembly:
%   d = r * DS_elast(:) + (1-r) * DS(:);
%   vals_pair = P.Wraw .* d(P.Rraw);
%   vals_u    = accumarray(P.g, vals_pair, [numel(P.Iu), 1]);
%   A         = sparsersb(P.Iu, P.Ju, vals_u, P.n, P.n, "unique", "sym");
%
% For the full (unsymmetric) sparse matrix:
%   vals_full = [vals_u; vals_u(P.off)];
%   A_sparse  = sparse(P.Ifull, P.Jfull, vals_full, P.n, P.n);
%
% INPUTS:
%   B       - Sparse strain-displacement matrix, size (n_strain*n_int) x n.
%   n_strain- Number of strain components (6 for 3D, 3 for 2D).
%   n_int   - Number of integration points.
%   WEIGHT  - Quadrature weights, size (1, n_int).
%
% OUTPUT:
%   P  - Struct with fields:
%          n      - Number of columns of B (DOF count).
%          Iu, Ju - Unique upper-triangle (i,j) pairs.
%          g      - Assignment vector  raw-pair -> unique-pair index.
%          Rraw   - Index into d = DS(:)  for each raw pair.
%          Wraw   - Precomputed weight = B_val_a * B_val_b * WEIGHT(q).
%          off    - Indices of off-diagonal unique pairs (Iu ~= Ju).
%          Ifull, Jfull - Full-matrix (upper+lower) row/col indices.
%--------------------------------------------------------------------------

n = size(B, 2);  % number of DOFs (columns)

% ---- Extract B in sorted-by-row COO format ----
[iB_raw, jB_raw, vB_raw] = find(B);
[iB_raw, perm] = sort(iB_raw);
jB_raw = jB_raw(perm);
vB_raw = vB_raw(perm);

m = n_strain * n_int;
row_counts = accumarray(iB_raw, 1, [m, 1]);
row_ptr    = [0; cumsum(row_counts)];       % 0-based pointers

% ---- For each strain-component sub-matrix B_a, build fixed-width arrays
%      B_a = B(a : n_strain : end, :)  is  n_int x n  sparse.
%      Each row has at most max_nnz_a non-zeros. ----

% Determine max nnz per row in each sub-matrix
max_nnz = zeros(n_strain, 1);
for a = 1:n_strain
    rows_a = a : n_strain : m;
    max_nnz(a) = max(row_counts(rows_a));
end

% Build dense (n_int x max_nnz_a) col-index and value arrays per component
cols_cell = cell(n_strain, 1);
vals_cell = cell(n_strain, 1);
for a = 1:n_strain
    mxa = max_nnz(a);
    c_mat = zeros(n_int, mxa);
    v_mat = zeros(n_int, mxa);
    for q = 1:n_int
        r = (q - 1) * n_strain + a;
        k = row_counts(r);
        if k == 0, continue; end
        idx = row_ptr(r) + (1:k);
        c_mat(q, 1:k) = jB_raw(idx);
        v_mat(q, 1:k) = vB_raw(idx);
    end
    cols_cell{a} = c_mat;
    vals_cell{a} = v_mat;
end

% ---- Generate raw pairs, vectorised over integration points ----
% Pre-estimate total raw pairs
total_est = 0;
for a = 1:n_strain
    for b = 1:n_strain
        total_est = total_est + max_nnz(a) * max_nnz(b) * n_int;
    end
end

Iraw = zeros(total_est, 1);
Jraw = zeros(total_est, 1);
Rraw = zeros(total_est, 1);
Wraw = zeros(total_est, 1);
pos  = 0;

W_col = WEIGHT(:);             % column vector (n_int x 1)
qq    = (0 : n_int - 1)';     % 0-based integration-point indices

for a = 1:n_strain
    ca = cols_cell{a};   % n_int x max_nnz(a)
    va = vals_cell{a};
    mxa = max_nnz(a);

    for b = 1:n_strain
        cb = cols_cell{b};
        vb = vals_cell{b};
        mxb = max_nnz(b);

        sub = (b - 1) * n_strain + a;   % linear index within DS column (1-based)

        for ia = 1:mxa
            ci = ca(:, ia);   % n_int x 1  column indices (0 where padded)
            vi = va(:, ia);

            for ib = 1:mxb
                cj = cb(:, ib);
                vj = vb(:, ib);

                valid = (ci > 0) & (cj > 0);
                nv = nnz(valid);
                if nv == 0, continue; end

                rr = pos + (1:nv);
                Iraw(rr) = ci(valid);
                Jraw(rr) = cj(valid);
                Wraw(rr) = vi(valid) .* vj(valid) .* W_col(valid);
                Rraw(rr) = sub + n_strain^2 * qq(valid);  % 1-based linear index into DS(:)
                pos = pos + nv;
            end
        end
    end
end

Iraw = Iraw(1:pos);
Jraw = Jraw(1:pos);
Rraw = Rraw(1:pos);
Wraw = Wraw(1:pos);

% ---- Fold to upper triangle (swap so that I <= J) ----
swap = Iraw > Jraw;
tmp       = Iraw(swap);
Iraw(swap) = Jraw(swap);
Jraw(swap) = tmp;

% ---- Reduce to unique (I, J) pairs ----
[IJu, ~, g] = unique([Iraw, Jraw], 'rows');
Iu = IJu(:, 1);
Ju = IJu(:, 2);

off = find(Iu ~= Ju);

% ---- Pack output struct ----
P = struct();
P.n     = n;
P.Iu    = Iu;
P.Ju    = Ju;
P.g     = g;
P.Rraw  = Rraw;
P.Wraw  = Wraw;
P.off   = off;
P.Ifull = [Iu; Ju(off)];
P.Jfull = [Ju; Iu(off)];

end
