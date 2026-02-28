function results = run_btDb_prebuilt_line_breakdown(varargin)
%RUN_BTDB_PREBUILT_LINE_BREAKDOWN Time each line of prebuilt A assembly.
%
% Timed lines (constructor='full'):
%   vals_pair = P.Wraw .* d(P.Rraw);
%   vals_u = accumarray(P.g, vals_pair, [numel(P.Iu), 1]);
%   vals_full = [vals_u; vals_u(P.off)];
%   A = sparsersb(P.Ifull, P.Jfull, vals_full, P.n, P.n);
%
% Timed lines (constructor='sym_unique'):
%   vals_pair = P.Wraw .* d(P.Rraw);
%   vals_u = accumarray(P.g, vals_pair, [numel(P.Iu), 1]);
%   A = sparsersb(P.Iu, P.Ju, vals_u, P.n, P.n, 'unique', 'sym');
%
% Usage:
%   run_btDb_prebuilt_line_breakdown()
%   run_btDb_prebuilt_line_breakdown('profile','medium','repeats',5,'n_d',12)
%   run_btDb_prebuilt_line_breakdown(...,'out_file','results.json')

opts = struct();
opts.profile = 'medium';
opts.repeats = 5;
opts.warmup = 1;
opts.n_d = 12;
opts.seed = 42;
opts.out_file = '';
opts.constructor = 'full';  % full | sym_unique

opts = parse_name_value(opts, varargin{:});
cfg = get_profile(opts.profile);
ctor_mode = lower(as_char(opts.constructor));
if ~strcmp(ctor_mode, 'full') && ~strcmp(ctor_mode, 'sym_unique')
    error('run_btDb_prebuilt_line_breakdown:BadConstructor', ...
        'Unknown constructor "%s". Use full or sym_unique.', ctor_mode);
end
set_seed(opts.seed);

if exist('OCTAVE_VERSION', 'builtin') == 0
    error('run_btDb_prebuilt_line_breakdown:OctaveOnly', ...
        'This benchmark uses sparsersb and must run in Octave.');
end

try
    pkg('load', 'sparsersb');
catch err
    error('run_btDb_prebuilt_line_breakdown:PkgLoadFailed', ...
        'Failed to load sparsersb: %s', err.message);
end

[engine_name, engine_version] = detect_engine();
timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');

fprintf('\n');
fprintf('=== Bt*diag(d)*B prebuilt line breakdown ===\n');
fprintf('Engine    : %s\n', engine_name);
fprintf('Version   : %s\n', engine_version);
fprintf('Timestamp : %s\n', timestamp);
fprintf('Profile   : %s\n', cfg.profile_name);
fprintf('Ctor mode : %s\n', ctor_mode);
fprintf('d vectors : %d per repeat\n', opts.n_d);
fprintf('Repeats   : %d (warmup: %d)\n', opts.repeats, opts.warmup);

[B, B_meta] = build_fem_like_B(cfg);
m = size(B, 1);
d_batch = 0.5 + rand(m, opts.n_d);

fprintf('\nMatrix B:\n');
fprintf('  size: %d x %d, nnz=%d\n', size(B, 1), size(B, 2), nnz(B));

t0 = tic();
P = precompute_BtDB_pattern_upper(B);
t_precompute = toc(t0);
fprintf('Precompute pattern: %.6f s (pairs=%d, unique=%d)\n', ...
    t_precompute, P.total_pairs, numel(P.Iu));

[line_totals, checksums, line_names] = time_lines(P, d_batch, opts.repeats, opts.warmup, ctor_mode); %#ok<ASGLU>

median_totals = median(line_totals, 1);
mean_totals = mean(line_totals, 1);
total_median = sum(median_totals);
median_per_d = median_totals / opts.n_d;
percent = 100 * median_totals / total_median;

fprintf('\nLine timing summary (per repeat, %d d-vectors):\n', opts.n_d);
fprintf('%-46s %12s %10s\n', 'Line', 'median_s', 'share_%');
for i = 1:numel(line_names)
    fprintf('%-46s %12.6f %10.2f\n', line_names{i}, median_totals(i), percent(i));
end
fprintf('%-46s %12.6f %10.2f\n', 'TOTAL', total_median, 100.0);

fprintf('\nLine timing summary (per single d-vector):\n');
fprintf('%-46s %12s\n', 'Line', 'median_s_per_d');
for i = 1:numel(line_names)
    fprintf('%-46s %12.6f\n', line_names{i}, median_per_d(i));
end

results = struct();
results.engine = engine_name;
results.version = engine_version;
results.timestamp = timestamp;
results.profile = cfg.profile_name;
results.seed = opts.seed;
results.repeats = opts.repeats;
results.warmup = opts.warmup;
results.n_d = opts.n_d;
results.constructor = ctor_mode;
results.precompute_s = t_precompute;
results.B = B_meta;
results.line_names = line_names;
results.line_totals_s = line_totals;
results.median_totals_s = median_totals;
results.mean_totals_s = mean_totals;
results.median_per_d_s = median_per_d;
results.share_percent = percent;
results.total_median_s = total_median;

if ~isempty(opts.out_file)
    write_results(opts.out_file, results);
    fprintf('\nSaved results to: %s\n', opts.out_file);
end
fprintf('\n');

end

function [line_totals, checksums, line_names] = time_lines(P, d_batch, repeats, warmup, ctor_mode)
if strcmp(ctor_mode, 'full')
    line_names = { ...
        'vals_pair = P.Wraw .* d(P.Rraw)', ...
        'vals_u = accumarray(P.g, vals_pair, ...)', ...
        'vals_full = [vals_u; vals_u(P.off)]', ...
        'A = sparsersb(P.Ifull, P.Jfull, vals_full, ...)'};
else
    line_names = { ...
        'vals_pair = P.Wraw .* d(P.Rraw)', ...
        'vals_u = accumarray(P.g, vals_pair, ...)', ...
        'A = sparsersb(P.Iu, P.Ju, vals_u, ..., ''unique'', ''sym'')'};
end

for i = 1:warmup
    [~, checksum] = run_lines_once(P, d_batch, ctor_mode); %#ok<NASGU>
end

line_totals = zeros(repeats, numel(line_names));
checksums = zeros(repeats, 1);
for i = 1:repeats
    [line_totals(i, :), checksums(i)] = run_lines_once(P, d_batch, ctor_mode);
end
end

function [line_sums, checksum] = run_lines_once(P, d_batch, ctor_mode)
if strcmp(ctor_mode, 'full')
    line_sums = zeros(1, 4);
else
    line_sums = zeros(1, 3);
end
checksum = 0.0;

for k = 1:size(d_batch, 2)
    d = d_batch(:, k);

    t = tic();
    vals_pair = P.Wraw .* d(P.Rraw);
    line_sums(1) = line_sums(1) + toc(t);

    t = tic();
    vals_u = accumarray(P.g, vals_pair, [numel(P.Iu), 1]);
    line_sums(2) = line_sums(2) + toc(t);

    if strcmp(ctor_mode, 'full')
        t = tic();
        vals_full = [vals_u; vals_u(P.off)];
        line_sums(3) = line_sums(3) + toc(t);

        t = tic();
        A = sparsersb(P.Ifull, P.Jfull, vals_full, P.n, P.n);
        line_sums(4) = line_sums(4) + toc(t);
    else
        t = tic();
        A = sparsersb(P.Iu, P.Ju, vals_u, P.n, P.n, 'unique', 'sym');
        line_sums(3) = line_sums(3) + toc(t);
    end

    checksum = checksum + double(nnz(A));
end
end

function [B, meta] = build_fem_like_B(cfg)
m = cfg.m_rows;
n = cfg.n_cols;
k = cfg.row_nnz;
w = cfg.local_window;

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

function opts = parse_name_value(opts, varargin)
if mod(numel(varargin), 2) ~= 0
    error('run_btDb_prebuilt_line_breakdown:InvalidInput', ...
        'Expected name/value pairs.');
end

for k = 1:2:numel(varargin)
    key = as_char(varargin{k});
    if ~isfield(opts, key)
        error('run_btDb_prebuilt_line_breakdown:UnknownOption', ...
            'Unknown option: %s', key);
    end
    opts.(key) = varargin{k + 1};
end

opts.profile = lower(as_char(opts.profile));
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
        error('run_btDb_prebuilt_line_breakdown:BadProfile', ...
            'Unknown profile "%s". Use small, medium, or large.', profile_name);
end
end

function s = as_char(x)
if ischar(x)
    s = x;
elseif isstring(x) && isscalar(x)
    s = char(x);
else
    error('run_btDb_prebuilt_line_breakdown:TypeError', ...
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
    error('run_btDb_prebuilt_line_breakdown:IOError', ...
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
