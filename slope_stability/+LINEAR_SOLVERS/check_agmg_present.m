function is_agmg_present = check_agmg_present(agmg_folder)
%CHECK_AGMG_PRESENT Adds AGMG to path and checks for agmg.m presence.
%
%   is_agmg_present = CHECK_AGMG_PRESENT()
%   Adds the folder 'agmg' to the MATLAB path and checks whether
%   the file 'agmg.m' exists in that folder.
%
%   Output:
%     is_agmg_present - 1 if agmg.m exists in the agmg folder, otherwise 0.

    addpath(agmg_folder);  % Add AGMG to the path

    % Check if agmg.m exists in the specified folder
    if exist(fullfile(agmg_folder, 'agmg.m'), 'file') == 2
        is_agmg_present = 1;
        % fprintf('AGMG found in folder "%s".\n', agmg_folder);
    else
        is_agmg_present = 0;
        % fprintf('AGMG not found in folder "%s".\n', agmg_folder);
    end
end
