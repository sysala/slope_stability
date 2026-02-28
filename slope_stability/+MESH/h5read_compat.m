function data = h5read_compat(file_path, dataset_path)
%H5READ_COMPAT Read an HDF5 dataset in MATLAB and Octave.
%
% Uses h5read when available (MATLAB, some Octave builds). For Octave
% builds without h5read, falls back to load('-hdf5', ...).

if nargin < 2
    error('Usage: data = MESH.h5read_compat(file_path, dataset_path)');
end

file_path = char(file_path);
dataset_path = char(dataset_path);

if exist('h5read', 'file') == 2 || exist('h5read', 'builtin') == 5
    data = h5read(file_path, dataset_path);
    return;
end

all_data = load('-hdf5', file_path);

dataset_key = regexprep(dataset_path, '^/+', '');
if isempty(dataset_key)
    error('Dataset path must be non-empty.');
end

if isfield(all_data, dataset_key)
    data = all_data.(dataset_key);
    return;
end

dataset_key_flat = strrep(dataset_key, '/', '_');
if isfield(all_data, dataset_key_flat)
    data = all_data.(dataset_key_flat);
    return;
end

parts = strsplit(dataset_key, '/');
dataset_leaf = parts{end};
if isfield(all_data, dataset_leaf)
    data = all_data.(dataset_leaf);
    return;
end

available_fields = fieldnames(all_data);
error(['Dataset "%s" not found in "%s". Available fields from ', ...
       'load(''-hdf5''): %s'], dataset_path, file_path, ...
      strjoin(available_fields, ', '));

end
