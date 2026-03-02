function tf = is_hypre_mex_available()
%IS_HYPRE_MEX_AVAILABLE True if hypre_boomeramg_mex binary exists on disk.
%
%   tf = LINEAR_SOLVERS.is_hypre_mex_available()
%
% Octave does not resolve exist() for mex files inside +package folders,
% so we check the file path directly.

this_dir = fileparts(mfilename('fullpath'));
tf = (exist(fullfile(this_dir, ['hypre_boomeramg_mex.' mexext]), 'file') == 3);
end
