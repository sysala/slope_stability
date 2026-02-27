% Get all .fig files in the current folder
figFiles = dir('*.fig');

% Loop over each .fig file
for i = 1:length(figFiles)
    % Open the figure
    figName = figFiles(i).name;
    openfig(figName);
    
    % Get the filename without extension
    [~, name, ~] = fileparts(figName);
    
    % Export as PNG with 600 dpi
    pngName = [name '.png'];
    exportgraphics(gcf, pngName, 'Resolution', 600);
    
    % Close the figure
    close(gcf);
    
    fprintf('Exported: %s\n', pngName);
end