function results = run_sparsersb_inputpath_benchmark(varargin)
%RUN_SPARSERSB_INPUTPATH_BENCHMARK Compare native sparse vs sparsersb paths.
%
% This benchmark answers:
%   - Is direct sparsersb(I,J,V,...) faster than sparsersb(sparse(...))?
%   - How does either compare to calling native sparse matrix matvec?
%
% Timed kernels use repeated SpMV to mimic iterative solver usage.
%
% Usage:
%   run_sparsersb_inputpath_benchmark()
%   run_sparsersb_inputpath_benchmark('profile','medium','matvec_repeats',30,'out_file','results.json')

opts = struct();
opts.profile = 'medium';
opts.repeats = 3;
opts.warmup = 1;
opts.matvec_repeats = 30;
opts.seed = 42;
opts.out_file = '';

opts = parse_name_value(opts, varargin{:});
cfg = get_profile(opts.profile);
set_seed(opts.seed);

if exist('OCTAVE_VERSION', 'builtin') == 0
    error('run_sparsersb_inputpath_benchmark:OctaveOnly', ...
        'This benchmark requires Octave with the sparsersb package.');
end

try
    pkg('load', 'sparsersb');
catch err
    error('run_sparsersb_inputpath_benchmark:PkgLoadFailed', ...
        'Failed to load sparsersb package: %s', err.message);
end

[engine_name, engine_version] = detect_engine();
timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');

fprintf('\n');
fprintf('=== sparsersb input-path benchmark ===\n');
fprintf('Engine         : %s\n', engine_name);
fprintf('Version        : %s\n', engine_version);
fprintf('Timestamp      : %s\n', timestamp);
fprintf('Profile        : %s\n', cfg.profile_name);
fprintf('Matvec repeats : %d\n', opts.matvec_repeats);
fprintf('Repeats        : %d (warmup: %d)\n', opts.repeats, opts.warmup);

[iA, jA, vA, n, X, meta] = build_spmv_case_from_coo(cfg, opts.matvec_repeats);
A_sparse = sparse(iA, jA, vA, n, n);
A_rsb_sparse = sparsersb(A_sparse);
A_rsb_coo = sparsersb(iA, jA, vA, n, n);

fprintf('\nMatrix:\n');
fprintf('  size=%d x %d, nnz=%d\n', n, n, meta.A_nnz);

cases = {};

% 1) Native sparse, prebuilt.
k1 = @() spmv_repeated(A_sparse, X);
s = time_kernel(k1, opts.repeats, opts.warmup);
s.kernel = '30*(A*x), sparse prebuilt';
s.meta = meta;
cases{end + 1} = s;

% 2) Native sparse, build + calls.
k2 = @() spmv_repeated(sparse(iA, jA, vA, n, n), X);
s = time_kernel(k2, opts.repeats, opts.warmup);
s.kernel = 'build sparse + 30*(A*x)';
s.meta = meta;
cases{end + 1} = s;

% 3) RSB from sparse, prebuilt.
k3 = @() spmv_repeated(A_rsb_sparse, X);
s = time_kernel(k3, opts.repeats, opts.warmup);
s.kernel = '30*(A_rsb*x), rsb prebuilt from sparse';
s.meta = meta;
cases{end + 1} = s;

% 4) RSB from sparse, conversion + calls (with prebuilt sparse).
k4 = @() spmv_repeated(sparsersb(A_sparse), X);
s = time_kernel(k4, opts.repeats, opts.warmup);
s.kernel = 'build rsb from sparse + 30*(A_rsb*x)';
s.meta = meta;
cases{end + 1} = s;

% 5) Full path via sparse then rsb + calls.
k5 = @() spmv_repeated(sparsersb(sparse(iA, jA, vA, n, n)), X);
s = time_kernel(k5, opts.repeats, opts.warmup);
s.kernel = 'build sparse + build rsb + 30*(A_rsb*x)';
s.meta = meta;
cases{end + 1} = s;

% 6) RSB from COO, prebuilt.
k6 = @() spmv_repeated(A_rsb_coo, X);
s = time_kernel(k6, opts.repeats, opts.warmup);
s.kernel = '30*(A_rsb*x), rsb prebuilt from coo';
s.meta = meta;
cases{end + 1} = s;

% 7) RSB from COO, direct build + calls.
k7 = @() spmv_repeated(sparsersb(iA, jA, vA, n, n), X);
s = time_kernel(k7, opts.repeats, opts.warmup);
s.kernel = 'build rsb from coo + 30*(A_rsb*x)';
s.meta = meta;
cases{end + 1} = s;

results = struct();
results.engine = engine_name;
results.version = engine_version;
results.timestamp = timestamp;
results.profile = cfg.profile_name;
results.seed = opts.seed;
results.repeats = opts.repeats;
results.warmup = opts.warmup;
results.matvec_repeats = opts.matvec_repeats;
results.cases = cases;

fprintf('\nSummary (seconds):\n');
fprintf('%-44s %10s %10s %10s\n', 'Kernel', 'median', 'min', 'max');
for i = 1:numel(cases)
    c = cases{i};
    fprintf('%-44s %10.6f %10.6f %10.6f\n', c.kernel, c.median_s, c.min_s, c.max_s);
end

if ~isempty(opts.out_file)
    write_results(opts.out_file, results);
    fprintf('\nSaved results to: %s\n', opts.out_file);
end
fprintf('\n');

end

function [iA, jA, vA, n, X, meta] = build_spmv_case_from_coo(cfg, matvec_repeats)
n = cfg.n_matvec;
offsets = unique(round(linspace(-cfg.diag_span, cfg.diag_span, cfg.diag_count)));
if ~any(offsets == 0)
    offsets = unique([offsets, 0]);
end

vals = randn(n, numel(offsets));
zero_col = find(offsets == 0, 1);
vals(:, zero_col) = vals(:, zero_col) + cfg.diag_boost;

A = spdiags(vals, offsets, n, n);
A = A + spdiags(0.01 * rand(n, 1), 0, n, n);
[iA, jA, vA] = find(A);
X = randn(n, matvec_repeats);

meta = struct();
meta.A_rows = n;
meta.A_cols = n;
meta.A_nnz = nnz(A);
meta.diag_count = numel(offsets);
meta.diag_span = cfg.diag_span;
meta.matvec_repeats = matvec_repeats;
end

function checksum = spmv_repeated(A, X)
acc = 0.0;
for i = 1:size(X, 2)
    y = A * X(:, i);
    k = min(2048, numel(y));
    acc = acc + sum(y(1:k)) * (1 + 1e-3 * i);
end
checksum = full(acc);
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
    error('run_sparsersb_inputpath_benchmark:InvalidInput', 'Expected name/value pairs.');
end

for i = 1:2:numel(varargin)
    key = as_char(varargin{i});
    if ~isfield(opts, key)
        error('run_sparsersb_inputpath_benchmark:UnknownOption', 'Unknown option: %s', key);
    end
    opts.(key) = varargin{i + 1};
end

opts.profile = lower(as_char(opts.profile));
opts.repeats = double(opts.repeats);
opts.warmup = double(opts.warmup);
opts.matvec_repeats = double(opts.matvec_repeats);
opts.seed = double(opts.seed);
opts.out_file = as_char(opts.out_file);
end

function s = as_char(x)
if ischar(x)
    s = x;
elseif isstring(x) && isscalar(x)
    s = char(x);
else
    error('run_sparsersb_inputpath_benchmark:TypeError', 'Expected char/string input.');
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

function cfg = get_profile(profile_name)
cfg = struct();
cfg.profile_name = profile_name;

switch lower(profile_name)
    case 'small'
        cfg.n_matvec = 200000;
        cfg.diag_count = 9;
        cfg.diag_span = 120;
        cfg.diag_boost = 6.0;
    case 'medium'
        cfg.n_matvec = 400000;
        cfg.diag_count = 11;
        cfg.diag_span = 180;
        cfg.diag_boost = 8.0;
    case 'large'
        cfg.n_matvec = 700000;
        cfg.diag_count = 13;
        cfg.diag_span = 220;
        cfg.diag_boost = 8.0;
    otherwise
        error('run_sparsersb_inputpath_benchmark:BadProfile', ...
            'Unknown profile "%s". Use small, medium, or large.', profile_name);
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
    error('run_sparsersb_inputpath_benchmark:IOError', 'Cannot open file: %s', out_file);
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
