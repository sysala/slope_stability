function results = run_benchmarks(varargin)
%RUN_BENCHMARKS Benchmark FEM-like kernels in MATLAB or Octave.
%
% Usage:
%   run_benchmarks()
%   run_benchmarks('profile', 'small')
%   run_benchmarks('profile', 'medium', 'repeats', 3, 'out_file', 'results.json')
%
% Kernels:
%   1) K_tangent-like assembly: K = B' * D_p * B
%   2) Large sparse matrix-vector multiply: y = A * x
%   3) Dense column dot-product matrix: G = X' * X

opts = struct();
opts.profile = 'medium';
opts.repeats = 3;
opts.warmup = 1;
opts.seed = 42;
opts.out_file = '';
opts.sparse_backend = 'native';

opts = parse_name_value(opts, varargin{:});
cfg = get_profile(opts.profile);
set_seed(opts.seed);

[engine_name, engine_version] = detect_engine();
timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');

backend = opts.sparse_backend;
if strcmp(backend, 'sparsersb')
    if exist('OCTAVE_VERSION', 'builtin') == 0
        error('run_benchmarks:BackendUnsupported', ...
            'sparse_backend="sparsersb" is only supported in Octave.');
    end
    try
        pkg('load', 'sparsersb');
    catch err
        error('run_benchmarks:BackendLoadError', ...
            'Failed to load sparsersb package: %s', err.message);
    end
end

fprintf('\n');
fprintf('=== slope_stability kernel benchmark ===\n');
fprintf('Engine     : %s\n', engine_name);
fprintf('Version    : %s\n', engine_version);
fprintf('Timestamp  : %s\n', timestamp);
fprintf('Profile    : %s\n', cfg.profile_name);
fprintf('Backend    : %s\n', backend);
fprintf('Repeats    : %d (warmup: %d)\n', opts.repeats, opts.warmup);

% --- Case 1: K_tangent = B' * D_p * B -----------------------------------
fprintf('\n[1/3] Building tangent-assembly test matrices...\n');
[B, D_p, tangent_meta] = build_tangent_case(cfg);
fprintf('      B: %d x %d (nnz=%d), D_p: %d x %d (nnz=%d)\n', ...
    tangent_meta.B_rows, tangent_meta.B_cols, tangent_meta.B_nnz, ...
    tangent_meta.D_rows, tangent_meta.D_cols, tangent_meta.D_nnz);

B_backend = apply_sparse_backend(B, backend);
D_backend = apply_sparse_backend(D_p, backend);
tangent_fun = @() tangent_kernel(B_backend, D_backend, cfg.symmetrize_tangent);
tangent_stats = time_kernel(tangent_fun, opts.repeats, opts.warmup);
tangent_stats.kernel = 'B'' * D_p * B';
tangent_stats.meta = tangent_meta;

clear B D_p B_backend D_backend;

% --- Case 2: Sparse A*x ---------------------------------------------------
fprintf('[2/3] Building sparse matvec test matrices...\n');
[A, x, matvec_meta] = build_sparse_matvec_case(cfg);
fprintf('      A: %d x %d (nnz=%d), x: %d x 1\n', ...
    matvec_meta.A_rows, matvec_meta.A_cols, matvec_meta.A_nnz, matvec_meta.x_len);

A_backend = apply_sparse_backend(A, backend);
matvec_fun = @() sparse_matvec_kernel(A_backend, x);
matvec_stats = time_kernel(matvec_fun, opts.repeats, opts.warmup);
matvec_stats.kernel = 'A * x (sparse)';
matvec_stats.meta = matvec_meta;

clear A A_backend x;

% --- Case 3: Dense X'*X ---------------------------------------------------
fprintf('[3/3] Building dense Gram test matrix...\n');
[X, dense_meta] = build_dense_gram_case(cfg);
fprintf('      X: %d x %d\n', dense_meta.X_rows, dense_meta.X_cols);

dense_fun = @() dense_gram_kernel(X);
dense_stats = time_kernel(dense_fun, opts.repeats, opts.warmup);
dense_stats.kernel = 'X'' * X (dense Gram)';
dense_stats.meta = dense_meta;

clear X;

cases = {tangent_stats, matvec_stats, dense_stats};

results = struct();
results.engine = engine_name;
results.version = engine_version;
results.timestamp = timestamp;
results.profile = cfg.profile_name;
results.sparse_backend = backend;
results.seed = opts.seed;
results.repeats = opts.repeats;
results.warmup = opts.warmup;
results.cases = cases;

fprintf('\nSummary (seconds):\n');
fprintf('%-24s %10s %10s %10s\n', 'Kernel', 'median', 'min', 'max');
fprintf('%-24s %10.6f %10.6f %10.6f\n', cases{1}.kernel, cases{1}.median_s, cases{1}.min_s, cases{1}.max_s);
fprintf('%-24s %10.6f %10.6f %10.6f\n', cases{2}.kernel, cases{2}.median_s, cases{2}.min_s, cases{2}.max_s);
fprintf('%-24s %10.6f %10.6f %10.6f\n', cases{3}.kernel, cases{3}.median_s, cases{3}.min_s, cases{3}.max_s);

if ~isempty(opts.out_file)
    write_results(opts.out_file, results);
    fprintf('\nSaved results to: %s\n', opts.out_file);
end
fprintf('\n');

end

function opts = parse_name_value(opts, varargin)
if mod(numel(varargin), 2) ~= 0
    error('run_benchmarks:InvalidInput', 'Expected name/value pairs.');
end
for i = 1:2:numel(varargin)
    key = lower(as_char(varargin{i}));
    val = varargin{i + 1};
    if ~isfield(opts, key)
        error('run_benchmarks:UnknownOption', 'Unknown option: %s', key);
    end
    opts.(key) = val;
end
opts.profile = lower(as_char(opts.profile));
opts.out_file = as_char(opts.out_file);
opts.sparse_backend = lower(as_char(opts.sparse_backend));
if ~strcmp(opts.sparse_backend, 'native') && ~strcmp(opts.sparse_backend, 'sparsersb')
    error('run_benchmarks:BadBackend', ...
        'Unknown sparse_backend "%s". Use native or sparsersb.', opts.sparse_backend);
end
end

function s = as_char(x)
if ischar(x)
    s = x;
elseif isa(x, 'string')
    s = char(x);
else
    error('run_benchmarks:TypeError', 'Expected char/string input.');
end
end

function set_seed(seed)
if exist('rng', 'file') == 2 || exist('rng', 'builtin') == 5
    rng(seed, 'twister');
else
    rand('seed', seed);
    randn('seed', seed);
end
end

function [name, verstr] = detect_engine()
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    name = 'Octave';
    verstr = OCTAVE_VERSION();
else
    name = 'MATLAB';
    verstr = version();
end
end

function cfg = get_profile(profile_name)
cfg = struct();
cfg.profile_name = profile_name;
cfg.symmetrize_tangent = true;

switch lower(profile_name)
    case 'small'
        cfg.n_nodes = 4000;
        cfg.n_elem = 6000;
        cfg.nodes_per_elem = 10;
        cfg.local_window = 28;

        cfg.n_matvec = 200000;
        cfg.diag_count = 9;
        cfg.diag_span = 120;
        cfg.diag_boost = 6.0;

        cfg.gram_rows = 30000;
        cfg.gram_cols = 64;

    case 'medium'
        cfg.n_nodes = 8000;
        cfg.n_elem = 12000;
        cfg.nodes_per_elem = 10;
        cfg.local_window = 36;

        cfg.n_matvec = 400000;
        cfg.diag_count = 11;
        cfg.diag_span = 180;
        cfg.diag_boost = 8.0;

        cfg.gram_rows = 60000;
        cfg.gram_cols = 96;

    case 'large'
        cfg.n_nodes = 12000;
        cfg.n_elem = 18000;
        cfg.nodes_per_elem = 10;
        cfg.local_window = 48;

        cfg.n_matvec = 700000;
        cfg.diag_count = 13;
        cfg.diag_span = 220;
        cfg.diag_boost = 8.0;

        cfg.gram_rows = 90000;
        cfg.gram_cols = 128;

    otherwise
        error('run_benchmarks:BadProfile', ...
            'Unknown profile "%s". Use small, medium, or large.', profile_name);
end
end

function [B, D_p, meta] = build_tangent_case(cfg)
dim = 3;
n_strain = 6;
ndof_per_elem = dim * cfg.nodes_per_elem;
n_rows_local = n_strain;
nz_local = n_rows_local * ndof_per_elem;

elem_nodes = generate_local_connectivity(cfg.n_nodes, cfg.n_elem, cfg.nodes_per_elem, cfg.local_window);

elem_dof = zeros(ndof_per_elem, cfg.n_elem);
cursor = 1;
for p = 1:cfg.nodes_per_elem
    base = dim * elem_nodes(p, :) - (dim - 1);
    elem_dof(cursor, :) = base;
    elem_dof(cursor + 1, :) = base + 1;
    elem_dof(cursor + 2, :) = base + 2;
    cursor = cursor + 3;
end

local_rows = kron((1:n_rows_local)', ones(1, ndof_per_elem));
local_cols = repmat(1:ndof_per_elem, n_rows_local, 1);

iB = bsxfun(@plus, local_rows(:), n_rows_local * (0:cfg.n_elem - 1));
jB = elem_dof(local_cols(:), :);
vB = randn(nz_local, cfg.n_elem);

B = sparse(iB(:), jB(:), vB(:), n_rows_local * cfg.n_elem, dim * cfg.n_nodes);

D0 = [ ...
    5.0, 1.2, 1.1, 0.1, 0.0, 0.0; ...
    1.2, 4.7, 1.0, 0.0, 0.1, 0.0; ...
    1.1, 1.0, 4.9, 0.0, 0.0, 0.1; ...
    0.1, 0.0, 0.0, 2.2, 0.1, 0.0; ...
    0.0, 0.1, 0.0, 0.1, 2.0, 0.1; ...
    0.0, 0.0, 0.1, 0.0, 0.1, 2.1];

[ii, jj, vv] = find(sparse(D0));
iD = bsxfun(@plus, ii, n_rows_local * (0:cfg.n_elem - 1));
jD = bsxfun(@plus, jj, n_rows_local * (0:cfg.n_elem - 1));
weights = 0.85 + 0.30 * rand(1, cfg.n_elem);
vD = vv * weights;

D_p = sparse(iD(:), jD(:), vD(:), n_rows_local * cfg.n_elem, n_rows_local * cfg.n_elem);

meta = struct();
meta.B_rows = size(B, 1);
meta.B_cols = size(B, 2);
meta.B_nnz = nnz(B);
meta.D_rows = size(D_p, 1);
meta.D_cols = size(D_p, 2);
meta.D_nnz = nnz(D_p);
meta.n_nodes = cfg.n_nodes;
meta.n_elem = cfg.n_elem;
meta.nodes_per_elem = cfg.nodes_per_elem;
end

function elem = generate_local_connectivity(n_nodes, n_elem, nodes_per_elem, local_window)
if local_window < nodes_per_elem
    error('run_benchmarks:BadConnectivity', 'local_window must be >= nodes_per_elem.');
end
if n_nodes < local_window
    error('run_benchmarks:BadConnectivity', 'n_nodes must be >= local_window.');
end

elem = zeros(nodes_per_elem, n_elem);
max_start = n_nodes - local_window + 1;
for e = 1:n_elem
    start_idx = randi(max_start);
    local_pick = randperm(local_window, nodes_per_elem);
    elem(:, e) = start_idx - 1 + local_pick(:);
end
end

function [A, x, meta] = build_sparse_matvec_case(cfg)
offsets = unique(round(linspace(-cfg.diag_span, cfg.diag_span, cfg.diag_count)));
if ~any(offsets == 0)
    offsets = unique([offsets, 0]);
end

vals = randn(cfg.n_matvec, numel(offsets));
zero_col = find(offsets == 0, 1);
vals(:, zero_col) = vals(:, zero_col) + cfg.diag_boost;

A = spdiags(vals, offsets, cfg.n_matvec, cfg.n_matvec);
A = A + spdiags(0.01 * rand(cfg.n_matvec, 1), 0, cfg.n_matvec, cfg.n_matvec);
x = randn(cfg.n_matvec, 1);

meta = struct();
meta.A_rows = size(A, 1);
meta.A_cols = size(A, 2);
meta.A_nnz = nnz(A);
meta.x_len = numel(x);
meta.diag_count = numel(offsets);
meta.diag_span = cfg.diag_span;
end

function [X, meta] = build_dense_gram_case(cfg)
X = randn(cfg.gram_rows, cfg.gram_cols);

meta = struct();
meta.X_rows = size(X, 1);
meta.X_cols = size(X, 2);
meta.X_bytes = numel(X) * 8;
end

function checksum = tangent_kernel(B, D_p, do_symmetrize)
K = B' * D_p * B;
if do_symmetrize
    K = 0.5 * (K + K');
end
k = min(64, size(K, 1));
checksum = full(sum(diag(K(1:k, 1:k)))) + double(nnz(K));
end

function checksum = sparse_matvec_kernel(A, x)
y = A * x;
k = min(2048, numel(y));
checksum = full(sum(y(1:k)));
end

function checksum = dense_gram_kernel(X)
G = X' * X;
k = min(64, size(G, 1));
checksum = full(sum(diag(G(1:k, 1:k))));
end

function M = apply_sparse_backend(M, backend)
if strcmp(backend, 'native')
    return;
end

if strcmp(backend, 'sparsersb')
    M = sparsersb(M);
    return;
end

error('run_benchmarks:BadBackend', ...
    'Unknown sparse backend "%s".', backend);
end

function stats = time_kernel(fun, repeats, warmup)
for i = 1:warmup
    checksum = fun();
end

times = zeros(repeats, 1);
for i = 1:repeats
    t0 = tic();
    checksum = fun();
    times(i) = toc(t0);
end

stats = struct();
stats.times_s = times(:)';
stats.median_s = median(times);
stats.mean_s = mean(times);
stats.min_s = min(times);
stats.max_s = max(times);
stats.checksum = double(checksum);
end

function write_results(out_file, results)
[out_dir, ~, ~] = fileparts(out_file);
if ~isempty(out_dir) && exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

json_ok = (exist('jsonencode', 'builtin') ~= 0) || (exist('jsonencode', 'file') == 2);
if json_ok
    try
        txt = jsonencode(results);
    catch
        txt = encode_json_value(results);
    end
else
    txt = encode_json_value(results);
end

fid = fopen(out_file, 'w');
if fid < 0
    error('run_benchmarks:IOError', 'Cannot open file: %s', out_file);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', txt);
end

function txt = encode_json_value(v)
if isstruct(v)
    if numel(v) > 1
        parts = cell(1, numel(v));
        for i = 1:numel(v)
            parts{i} = encode_json_value(v(i));
        end
        txt = ['[', strjoin(parts, ','), ']'];
        return;
    end
    f = fieldnames(v);
    parts = cell(1, numel(f));
    for i = 1:numel(f)
        key = f{i};
        parts{i} = ['"', escape_json_string(key), '":', encode_json_value(v.(key))];
    end
    txt = ['{', strjoin(parts, ','), '}'];
    return;
end

if iscell(v)
    parts = cell(1, numel(v));
    for i = 1:numel(v)
        parts{i} = encode_json_value(v{i});
    end
    txt = ['[', strjoin(parts, ','), ']'];
    return;
end

if ischar(v)
    txt = ['"', escape_json_string(v), '"'];
    return;
end

if isa(v, 'string')
    if isscalar(v)
        txt = ['"', escape_json_string(char(v)), '"'];
    else
        c = cellstr(v);
        parts = cell(1, numel(c));
        for i = 1:numel(c)
            parts{i} = ['"', escape_json_string(c{i}), '"'];
        end
        txt = ['[', strjoin(parts, ','), ']'];
    end
    return;
end

if islogical(v)
    if isscalar(v)
        if v
            txt = 'true';
        else
            txt = 'false';
        end
    else
        parts = cell(1, numel(v));
        for i = 1:numel(v)
            if v(i)
                parts{i} = 'true';
            else
                parts{i} = 'false';
            end
        end
        txt = ['[', strjoin(parts, ','), ']'];
    end
    return;
end

if isnumeric(v)
    if isempty(v)
        txt = '[]';
        return;
    end
    if isscalar(v)
        txt = sprintf('%.17g', v);
    else
        vals = v(:);
        parts = cell(1, numel(vals));
        for i = 1:numel(vals)
            parts{i} = sprintf('%.17g', vals(i));
        end
        txt = ['[', strjoin(parts, ','), ']'];
    end
    return;
end

txt = 'null';
end

function s = escape_json_string(s)
s = strrep(s, '\', '\\');
s = strrep(s, '"', '\"');
s = strrep(s, sprintf('\n'), '\n');
s = strrep(s, sprintf('\r'), '\r');
s = strrep(s, sprintf('\t'), '\t');
end
