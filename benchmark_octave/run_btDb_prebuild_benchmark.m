function results = run_btDb_prebuild_benchmark(varargin)
%RUN_BTDB_PREBUILD_BENCHMARK Benchmark naive vs prebuilt Bt*diag(d)*B assembly.
%
% Problem:
%   A(d) = B' * diag(d) * B
% for many vectors d with fixed sparse matrix B.
%
% Compared kernels:
%   1) Naive:
%      - native backend:   A = B' * spdiags(d,0,m,m) * B
%      - sparsersb backend: A = Br' * Dr * Br
%   2) Prebuilt:
%      - precompute pair pattern from B rows once
%      - for each d evaluate values via accumarray
%      - construct result matrix from prebuilt (I,J,values)
%
% Usage examples:
%   run_btDb_prebuild_benchmark()
%   run_btDb_prebuild_benchmark('profile','medium','backend','native')
%   run_btDb_prebuild_benchmark('profile','medium','backend','sparsersb')
%   run_btDb_prebuild_benchmark(...,'out_file','results.json')

opts = struct();
opts.profile = 'medium';
opts.backend = 'native';  % native | sparsersb
opts.repeats = 3;
opts.warmup = 1;
opts.n_d = 12;            % number of d vectors per timed call
opts.seed = 42;
opts.out_file = '';

opts = parse_name_value(opts, varargin{:});
cfg = get_profile(opts.profile);
set_seed(opts.seed);

backend = lower(opts.backend);
if ~strcmp(backend, 'native') && ~strcmp(backend, 'sparsersb')
    error('run_btDb_prebuild_benchmark:BadBackend', ...
        'Unknown backend "%s". Use native or sparsersb.', backend);
end

if strcmp(backend, 'sparsersb')
    if exist('OCTAVE_VERSION', 'builtin') == 0
        error('run_btDb_prebuild_benchmark:OctaveOnly', ...
            'backend="sparsersb" is available only in Octave.');
    end
    try
        pkg('load', 'sparsersb');
    catch err
        error('run_btDb_prebuild_benchmark:PkgLoadFailed', ...
            'Failed to load sparsersb: %s', err.message);
    end
end

[engine_name, engine_version] = detect_engine();
timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');

fprintf('\n');
fprintf('=== Bt*diag(d)*B prebuild benchmark ===\n');
fprintf('Engine      : %s\n', engine_name);
fprintf('Version     : %s\n', engine_version);
fprintf('Timestamp   : %s\n', timestamp);
fprintf('Profile     : %s\n', cfg.profile_name);
fprintf('Backend     : %s\n', backend);
fprintf('d vectors   : %d per call\n', opts.n_d);
fprintf('Repeats     : %d (warmup: %d)\n', opts.repeats, opts.warmup);

[B, B_meta] = build_fem_like_B(cfg);
m = size(B, 1);
n = size(B, 2);
d_batch = 0.5 + rand(m, opts.n_d);  % positive diagonal entries

fprintf('\nMatrix B:\n');
fprintf('  size: %d x %d, nnz=%d, avg nnz/row=%.2f\n', ...
    m, n, nnz(B), nnz(B) / m);

fprintf('Precomputing row-pair pattern...\n');
t0 = tic();
P = precompute_BtDB_pattern_upper(B);
t_precompute = toc(t0);
fprintf('  upper pair entries: %d, unique (I,J): %d, precompute: %.6f s\n', ...
    P.total_pairs, numel(P.Iu), t_precompute);

if strcmp(backend, 'sparsersb')
    Br = sparsersb(B);
    diag_i = (1:m)';
else
    Br = [];
    diag_i = [];
end

naive_fun = @() kernel_naive(B, Br, diag_i, d_batch, backend);
prebuilt_fun = @() kernel_prebuilt(P, d_batch, backend);

naive_stats = time_kernel(naive_fun, opts.repeats, opts.warmup);
naive_stats.kernel = 'naive B''*diag(d)*B';
naive_stats.meta = B_meta;

prebuilt_stats = time_kernel(prebuilt_fun, opts.repeats, opts.warmup);
prebuilt_stats.kernel = 'prebuilt pairs + accumarray';
prebuilt_stats.meta = B_meta;

cases = {naive_stats, prebuilt_stats};

results = struct();
results.engine = engine_name;
results.version = engine_version;
results.timestamp = timestamp;
results.profile = cfg.profile_name;
results.backend = backend;
results.seed = opts.seed;
results.repeats = opts.repeats;
results.warmup = opts.warmup;
results.n_d = opts.n_d;
results.precompute_s = t_precompute;
results.B = B_meta;
results.precompute_meta = struct('total_pairs', P.total_pairs, 'unique_pairs', numel(P.Iu));
results.cases = cases;

fprintf('\nSummary (seconds):\n');
fprintf('%-34s %10s %10s %10s\n', 'Kernel', 'median', 'min', 'max');
fprintf('%-34s %10.6f %10.6f %10.6f\n', naive_stats.kernel, naive_stats.median_s, naive_stats.min_s, naive_stats.max_s);
fprintf('%-34s %10.6f %10.6f %10.6f\n', prebuilt_stats.kernel, prebuilt_stats.median_s, prebuilt_stats.min_s, prebuilt_stats.max_s);

speedup = naive_stats.median_s / prebuilt_stats.median_s;
fprintf('\nSpeedup (naive / prebuilt): %.3fx\n', speedup);

if ~isempty(opts.out_file)
    write_results(opts.out_file, results);
    fprintf('Saved results to: %s\n', opts.out_file);
end
fprintf('\n');

end

function [B, meta] = build_fem_like_B(cfg)
m = cfg.m_rows;
n = cfg.n_cols;
k = cfg.row_nnz;
w = cfg.local_window;

if w < k
    error('run_btDb_prebuild_benchmark:BadWindow', ...
        'local_window must be >= row_nnz.');
end
if n < w
    error('run_btDb_prebuild_benchmark:BadSize', ...
        'n_cols must be >= local_window.');
end

i = repelem((1:m)', k, 1);
j = zeros(m * k, 1);
v = randn(m * k, 1);

max_start = n - w + 1;
for r = 1:m
    start_idx = randi(max_start);
    cols = start_idx - 1 + randperm(w, k);
    idx = (r - 1) * k + (1:k);
    j(idx) = cols(:);
end

B = sparse(i, j, v, m, n);

meta = struct();
meta.m_rows = m;
meta.n_cols = n;
meta.nnz = nnz(B);
meta.row_nnz = k;
meta.local_window = w;
end

function P = precompute_BtDB_pattern_upper(B)
[iB, jB, vB] = find(B);
m = size(B, 1);
n = size(B, 2);

[iB, p] = sort(iB);
jB = jB(p);
vB = vB(p);

row_counts = accumarray(iB, 1, [m, 1]);
row_ptr = [1; cumsum(row_counts) + 1];

pair_counts = row_counts .* (row_counts + 1) / 2;
total_pairs = sum(pair_counts);

Iraw = zeros(total_pairs, 1);
Jraw = zeros(total_pairs, 1);
Rraw = zeros(total_pairs, 1);
Wraw = zeros(total_pairs, 1);

pos = 1;
for r = 1:m
    k = row_counts(r);
    if k == 0
        continue;
    end

    idx = row_ptr(r):(row_ptr(r + 1) - 1);
    cols = jB(idx);
    vals = vB(idx);

    [CJ, CI] = meshgrid(cols, cols);
    [VJ, VI] = meshgrid(vals, vals);
    mask = triu(true(k));

    t = nnz(mask);
    rr = pos:(pos + t - 1);
    Iraw(rr) = CI(mask);
    Jraw(rr) = CJ(mask);
    Rraw(rr) = r;
    Wraw(rr) = VI(mask) .* VJ(mask);
    pos = pos + t;
end

if pos <= total_pairs
    keep = 1:(pos - 1);
    Iraw = Iraw(keep);
    Jraw = Jraw(keep);
    Rraw = Rraw(keep);
    Wraw = Wraw(keep);
end

[IJu, ~, g] = unique([Iraw, Jraw], 'rows');
Iu = IJu(:, 1);
Ju = IJu(:, 2);
off = find(Iu ~= Ju);

P = struct();
P.m = m;
P.n = n;
P.total_pairs = numel(Wraw);
P.Iu = Iu;
P.Ju = Ju;
P.g = g;
P.Rraw = Rraw;
P.Wraw = Wraw;
P.off = off;
P.Ifull = [Iu; Ju(off)];
P.Jfull = [Ju; Iu(off)];
end

function checksum = kernel_naive(B, Br, diag_i, d_batch, backend)
m = size(B, 1);
acc = 0.0;

for k = 1:size(d_batch, 2)
    d = d_batch(:, k);
    if strcmp(backend, 'native')
        D = spdiags(d, 0, m, m);
        A = B' * D * B;
    else
        D = sparsersb(diag_i, diag_i, d, m, m);
        A = Br' * D * Br;
    end
    acc = acc + double(nnz(A));
end

checksum = acc;
end

function checksum = kernel_prebuilt(P, d_batch, backend)
acc = 0.0;

for k = 1:size(d_batch, 2)
    d = d_batch(:, k);
    vals_pair = P.Wraw .* d(P.Rraw);
    vals_u = accumarray(P.g, vals_pair, [numel(P.Iu), 1]);

    if strcmp(backend, 'native')
        A = sparse(P.Iu, P.Ju, vals_u, P.n, P.n);
        A = A + triu(A, 1)';
    else
        vals_full = [vals_u; vals_u(P.off)];
        A = sparsersb(P.Ifull, P.Jfull, vals_full, P.n, P.n);
    end
    acc = acc + double(nnz(A));
end

checksum = acc;
end

function stats = time_kernel(fun, repeats, warmup)
for i = 1:warmup
    checksum = fun(); %#ok<NASGU>
end

times = zeros(repeats, 1);
for i = 1:repeats
    t0 = tic();
    checksum = fun(); %#ok<NASGU>
    times(i) = toc(t0);
end

stats = struct();
stats.times_s = times(:)';
stats.median_s = median(times);
stats.mean_s = mean(times);
stats.min_s = min(times);
stats.max_s = max(times);
stats.checksum = checksum;
end

function opts = parse_name_value(opts, varargin)
if mod(numel(varargin), 2) ~= 0
    error('run_btDb_prebuild_benchmark:InvalidInput', ...
        'Expected name/value pairs.');
end

for k = 1:2:numel(varargin)
    key = as_char(varargin{k});
    if ~isfield(opts, key)
        error('run_btDb_prebuild_benchmark:UnknownOption', ...
            'Unknown option: %s', key);
    end
    opts.(key) = varargin{k + 1};
end

opts.profile = lower(as_char(opts.profile));
opts.backend = lower(as_char(opts.backend));
opts.repeats = double(opts.repeats);
opts.warmup = double(opts.warmup);
opts.n_d = double(opts.n_d);
opts.seed = double(opts.seed);
opts.out_file = as_char(opts.out_file);
end

function cfg = get_profile(profile_name)
cfg = struct();
cfg.profile_name = profile_name;

switch lower(profile_name)
    case 'small'
        cfg.m_rows = 8000;
        cfg.n_cols = 3000;
        cfg.row_nnz = 10;
        cfg.local_window = 40;
    case 'medium'
        cfg.m_rows = 20000;
        cfg.n_cols = 7000;
        cfg.row_nnz = 12;
        cfg.local_window = 56;
    case 'large'
        cfg.m_rows = 35000;
        cfg.n_cols = 12000;
        cfg.row_nnz = 12;
        cfg.local_window = 64;
    otherwise
        error('run_btDb_prebuild_benchmark:BadProfile', ...
            'Unknown profile "%s". Use small, medium, or large.', profile_name);
end
end

function s = as_char(x)
if ischar(x)
    s = x;
elseif isstring(x) && isscalar(x)
    s = char(x);
else
    error('run_btDb_prebuild_benchmark:TypeError', ...
        'Expected char/string input.');
end
end

function set_seed(seed)
if exist('rng', 'file') == 2 || exist('rng', 'builtin') ~= 0
    rng(seed, 'twister');
else
    rand('seed', seed);
    randn('seed', seed);
end
end

function [name, verstr] = detect_engine()
if exist('OCTAVE_VERSION', 'builtin')
    name = 'Octave';
    verstr = OCTAVE_VERSION();
else
    name = 'MATLAB';
    verstr = version();
end
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
    error('run_btDb_prebuild_benchmark:IOError', ...
        'Cannot open file: %s', out_file);
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, txt, 'char');
end

function txt = encode_json_value(v)
if isnumeric(v)
    if isempty(v)
        txt = '[]';
    elseif isscalar(v)
        if isfinite(v)
            txt = sprintf('%.17g', v);
        elseif isnan(v)
            txt = 'null';
        elseif v > 0
            txt = '1e999';
        else
            txt = '-1e999';
        end
    else
        parts = cell(1, numel(v));
        for i = 1:numel(v)
            parts{i} = encode_json_value(v(i));
        end
        txt = ['[', strjoin(parts, ','), ']'];
    end
elseif islogical(v)
    if isscalar(v)
        if v
            txt = 'true';
        else
            txt = 'false';
        end
    else
        parts = cell(1, numel(v));
        for i = 1:numel(v)
            parts{i} = encode_json_value(v(i));
        end
        txt = ['[', strjoin(parts, ','), ']'];
    end
elseif isstruct(v)
    if numel(v) ~= 1
        parts = cell(1, numel(v));
        for i = 1:numel(v)
            parts{i} = encode_json_value(v(i));
        end
        txt = ['[', strjoin(parts, ','), ']'];
        return;
    end
    keys = fieldnames(v);
    parts = cell(1, numel(keys));
    for i = 1:numel(keys)
        key = keys{i};
        parts{i} = ['"', escape_json_string(key), '":', encode_json_value(v.(key))];
    end
    txt = ['{', strjoin(parts, ','), '}'];
elseif iscell(v)
    parts = cell(1, numel(v));
    for i = 1:numel(v)
        parts{i} = encode_json_value(v{i});
    end
    txt = ['[', strjoin(parts, ','), ']'];
elseif ischar(v)
    txt = ['"', escape_json_string(v), '"'];
elseif isstring(v)
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
else
    error('Unsupported type for JSON encoding.');
end
end

function s = escape_json_string(s)
s = strrep(s, '\', '\\');
s = strrep(s, '"', '\"');
s = strrep(s, sprintf('\b'), '\b');
s = strrep(s, sprintf('\f'), '\f');
s = strrep(s, sprintf('\n'), '\n');
s = strrep(s, sprintf('\r'), '\r');
s = strrep(s, sprintf('\t'), '\t');
end
