function build_hypre_boomeramg_mex(hypre_root)
%BUILD_HYPRE_BOOMERAMG_MEX Build hypre_boomeramg_mex against local HYPRE.
%
%   LINEAR_SOLVERS.build_hypre_boomeramg_mex()
%   LINEAR_SOLVERS.build_hypre_boomeramg_mex(hypre_root)
%
% Default hypre_root:
%   third_party/hypre-openmp   (relative to repository root)
%
% Requires:
%   - mex configured
%   - libHYPRE available in hypre_root/lib
%   - HYPRE headers in hypre_root/include

this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(this_dir);
repo_dir = fileparts(root_dir);

if nargin < 1 || isempty(hypre_root)
    hypre_root = fullfile(repo_dir, 'third_party', 'hypre-openmp');
end
hypre_root = char(hypre_root);

inc_dir = fullfile(hypre_root, 'include');
lib_dir = fullfile(hypre_root, 'lib');
src = fullfile(root_dir, 'mex', 'hypre_boomeramg_mex.cpp');

if ~isfolder(inc_dir)
    error('Missing include directory: %s', inc_dir);
end
if ~isfolder(lib_dir)
    error('Missing library directory: %s', lib_dir);
end
if ~isfile(src)
    error('Missing source file: %s', src);
end

fprintf('Building hypre_boomeramg_mex\n');
fprintf('  source: %s\n', src);
fprintf('  include: %s\n', inc_dir);
fprintf('  lib: %s\n', lib_dir);

out_dir = root_dir;
mex_file = fullfile(out_dir, ['hypre_boomeramg_mex.', mexext]);
if isfile(mex_file)
    delete(mex_file);
end

try
    mex( ...
        '-R2018a', ...
        sprintf('CXXFLAGS=$CXXFLAGS -O3 -std=c++17 -fopenmp'), ...
        sprintf('LDFLAGS=$LDFLAGS -fopenmp -Wl,-rpath,%s', lib_dir), ...
        ['-I', inc_dir], ...
        ['-L', lib_dir], ...
        '-lHYPRE', ...
        src, ...
        '-outdir', out_dir);
catch ME
    % Some MATLAB setups incorrectly raise ENOTMEX right after successful link.
    if contains(ME.message, 'is not a MEX file') && isfile(mex_file)
        warning('build_hypre_boomeramg_mex:postcheck', ...
            ['mex post-build check reported ENOTMEX, but output exists. ', ...
             'Trying runtime load check...']);
        addpath(out_dir);
        clear mex;
        hypre_boomeramg_mex('info'); %#ok<NASGU>
        clear mex;
    else
        rethrow(ME);
    end
end

fprintf('Built: %s\n', fullfile(out_dir, ['hypre_boomeramg_mex.', mexext]));

end
