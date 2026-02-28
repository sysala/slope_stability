function results = run_comsol_matfree_benchmark(varargin)
%RUN_COMSOL_MATFREE_BENCHMARK Benchmark explicit vs mat-free tangent matvecs.
%
% This benchmark uses matrices assembled from the
% slope_stability_3D_hetero_seepage_SSR_comsol setup:
%   - finite-element B matrix from COMSOL mesh path
%   - constitutive D_p assembled from DS at U = 0, lambda = 1
%
% Timed kernels:
%   1) Build K = Bq' * D_p * Bq and run 30 matvecs (K*x)
%   2) Run 30 matvecs (K*x) with prebuilt K
%   3) Run 30 mat-free products Bq' * (D_p * (Bq*x))
%
% Usage:
%   run_comsol_matfree_benchmark()
%   run_comsol_matfree_benchmark('element_limit', 8000, 'out_file', 'result.json')

opts = struct();
opts.mesh_file = fullfile('slope_stability', 'meshes', 'comsol_mesh.h5');
opts.element_limit = 8000;
opts.matvec_repeats = 30;
opts.repeats = 3;
opts.warmup = 1;
opts.seed = 42;
opts.out_file = '';

opts = parse_name_value(opts, varargin{:});
set_seed(opts.seed);

[engine_name, engine_version] = detect_engine();
timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');

repo_root = fileparts(fileparts(mfilename('fullpath')));
slope_root = fullfile(repo_root, 'slope_stability');
addpath(slope_root);

mesh_path = as_char(opts.mesh_file);
if ~is_absolute_path(mesh_path)
    mesh_path = fullfile(repo_root, mesh_path);
end
if exist(mesh_path, 'file') ~= 2
    error('run_comsol_matfree_benchmark:MissingMesh', ...
        'Mesh file not found: %s', mesh_path);
end

fprintf('\n');
fprintf('=== COMSOL mat-free benchmark ===\n');
fprintf('Engine         : %s\n', engine_name);
fprintf('Version        : %s\n', engine_version);
fprintf('Timestamp      : %s\n', timestamp);
fprintf('Mesh file      : %s\n', mesh_path);
fprintf('Element limit  : %d\n', opts.element_limit);
fprintf('Matvec repeats : %d\n', opts.matvec_repeats);
fprintf('Repeats        : %d (warmup: %d)\n', opts.repeats, opts.warmup);

[Bq, D_p, matvec_inputs, meta] = build_comsol_case(mesh_path, opts.element_limit, opts.matvec_repeats);

fprintf('\nCase matrices:\n');
fprintf('  Bq: %d x %d (nnz=%d)\n', size(Bq, 1), size(Bq, 2), nnz(Bq));
fprintf('  D_p: %d x %d (nnz=%d)\n', size(D_p, 1), size(D_p, 2), nnz(D_p));

explicit_total_fun = @() kernel_explicit_total(Bq, D_p, matvec_inputs);
explicit_total = time_kernel(explicit_total_fun, opts.repeats, opts.warmup);
explicit_total.kernel = 'build K + 30*(K*x)';
explicit_total.meta = meta;

K_cached = Bq' * D_p * Bq;
explicit_only_fun = @() kernel_explicit_only(K_cached, matvec_inputs);
explicit_only = time_kernel(explicit_only_fun, opts.repeats, opts.warmup);
explicit_only.kernel = '30*(K*x), K prebuilt';
explicit_only.meta = meta;

matfree_fun = @() kernel_matfree(Bq, D_p, matvec_inputs);
matfree = time_kernel(matfree_fun, opts.repeats, opts.warmup);
matfree.kernel = '30*(B''*(D_p*(B*x)))';
matfree.meta = meta;

cases = {explicit_total, explicit_only, matfree};

results = struct();
results.engine = engine_name;
results.version = engine_version;
results.timestamp = timestamp;
results.seed = opts.seed;
results.repeats = opts.repeats;
results.warmup = opts.warmup;
results.matvec_repeats = opts.matvec_repeats;
results.mesh_file = mesh_path;
results.element_limit = opts.element_limit;
results.cases = cases;
results.meta = meta;

fprintf('\nSummary (seconds):\n');
fprintf('%-28s %10s %10s %10s\n', 'Kernel', 'median', 'min', 'max');
fprintf('%-28s %10.6f %10.6f %10.6f\n', cases{1}.kernel, cases{1}.median_s, cases{1}.min_s, cases{1}.max_s);
fprintf('%-28s %10.6f %10.6f %10.6f\n', cases{2}.kernel, cases{2}.median_s, cases{2}.min_s, cases{2}.max_s);
fprintf('%-28s %10.6f %10.6f %10.6f\n', cases{3}.kernel, cases{3}.median_s, cases{3}.min_s, cases{3}.max_s);

if ~isempty(opts.out_file)
    write_results(opts.out_file, results);
    fprintf('\nSaved results to: %s\n', opts.out_file);
end
fprintf('\n');

end

function [Bq, D_p, matvec_inputs, meta] = build_comsol_case(mesh_path, element_limit, matvec_repeats)
elem_type = 'P2';
[Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);
[~, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);

[coord, elem, ~, Q, material, ~] = MESH.load_mesh_P2(mesh_path, 1);
n_elem_total = size(elem, 2);
if element_limit > 0 && element_limit < n_elem_total
    keep = 1:element_limit;
    elem = elem(:, keep);
    material = material(keep);
end

n_nodes = size(coord, 2);
n_elem = size(elem, 2);
n_q = numel(WF);
n_int = n_elem * n_q;
n_strain = 6;

mat_props = [ ...
    15, 30, 0, 10000, 0.33, 19, 19; ...
    15, 38, 0, 50000, 0.30, 22, 22; ...
    10, 35, 0, 50000, 0.30, 21, 21; ...
    18, 32, 0, 20000, 0.33, 20, 20];

fields = {'c0', 'phi', 'psi', 'young', 'poisson', 'gamma_sat', 'gamma_unsat'};
materials = cellfun(@(x) cell2struct(num2cell(x), fields, 2), ...
    num2cell(mat_props, 2), 'UniformOutput', false);
saturation = true(1, n_int);
[c0, phi, psi, shear, bulk, lame, ~] = ASSEMBLY.heterogenous_materials(material, saturation, n_q, materials);

[~, B, WEIGHT] = ASSEMBLY.elastic_stiffness_matrix_3D(elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

builder = CONSTITUTIVE_PROBLEM.CONSTITUTIVE( ...
    B, c0, phi, psi, 'B', shear, bulk, lame, WEIGHT, n_strain, n_int, 3);
builder.reduction(1.0);
U0 = zeros(3, n_nodes);
builder.constitutive_problem_stress_tangent(U0);

vD = builder.vD_pre .* builder.DS;
D_p = sparse(builder.iD(:), builder.jD(:), vD(:), n_strain * n_int, n_strain * n_int);

q_idx = find(Q(:));
Bq = B(:, q_idx);
matvec_inputs = randn(numel(q_idx), matvec_repeats);

meta = struct();
meta.n_nodes = n_nodes;
meta.n_elem_total = n_elem_total;
meta.n_elem_used = n_elem;
meta.n_q = n_q;
meta.n_int = n_int;
meta.n_unknown = numel(q_idx);
meta.B_rows = size(Bq, 1);
meta.B_cols = size(Bq, 2);
meta.B_nnz = nnz(Bq);
meta.D_rows = size(D_p, 1);
meta.D_cols = size(D_p, 2);
meta.D_nnz = nnz(D_p);
meta.matvec_repeats = matvec_repeats;
end

function checksum = kernel_explicit_total(Bq, D_p, matvec_inputs)
K = Bq' * D_p * Bq;
checksum = kernel_explicit_only(K, matvec_inputs) + 1e-12 * double(nnz(K));
end

function checksum = kernel_explicit_only(K, matvec_inputs)
acc = 0.0;
n_cols = size(matvec_inputs, 2);
for i = 1:n_cols
    y = K * matvec_inputs(:, i);
    k = min(64, numel(y));
    acc = acc + double(sum(y(1:k))) * (1.0 + 0.001 * i);
end
checksum = acc;
end

function checksum = kernel_matfree(Bq, D_p, matvec_inputs)
acc = 0.0;
n_cols = size(matvec_inputs, 2);
for i = 1:n_cols
    y = Bq' * (D_p * (Bq * matvec_inputs(:, i)));
    k = min(64, numel(y));
    acc = acc + double(sum(y(1:k))) * (1.0 + 0.001 * i);
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
    error('run_comsol_matfree_benchmark:InvalidInput', 'Expected name/value pairs.');
end

for k = 1:2:numel(varargin)
    key = as_char(varargin{k});
    if ~isfield(opts, key)
        error('run_comsol_matfree_benchmark:UnknownOption', 'Unknown option: %s', key);
    end
    opts.(key) = varargin{k + 1};
end

opts.mesh_file = as_char(opts.mesh_file);
opts.element_limit = double(opts.element_limit);
opts.matvec_repeats = double(opts.matvec_repeats);
opts.repeats = double(opts.repeats);
opts.warmup = double(opts.warmup);
opts.seed = double(opts.seed);
opts.out_file = as_char(opts.out_file);
end

function s = as_char(x)
if ischar(x)
    s = x;
elseif isstring(x) && isscalar(x)
    s = char(x);
else
    error('run_comsol_matfree_benchmark:TypeError', 'Expected char/string input.');
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

function tf = is_absolute_path(p)
if isempty(p)
    tf = false;
    return;
end

if ispc
    tf = (numel(p) >= 2 && p(2) == ':') || (numel(p) >= 2 && p(1) == '\' && p(2) == '\');
else
    tf = (p(1) == '/');
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
    error('run_comsol_matfree_benchmark:IOError', 'Cannot open file: %s', out_file);
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
