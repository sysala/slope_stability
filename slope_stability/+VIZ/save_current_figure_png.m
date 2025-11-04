function save_current_figure_png(filename, fig_width, fig_height)
%SAVE_CURRENT_FIGURE_PNG Save current figure to PNG with specified size and 600 dpi
%
%   save_current_figure_png(filename, fig_width, fig_height)
%
%   Inputs:
%       filename   - output filename (string, e.g. 'figure.png')
%       fig_width  - width in centimeters
%       fig_height - height in centimeters
%
%   Example:
%       save_current_figure_png('out.png', 12, 8);

    if nargin < 3
        error('Usage: save_current_figure_png(filename, fig_width, fig_height)');
    end

    % Get current figure
    fig = gcf;

    % Set figure size in centimeters
    set(fig, 'Units', 'centimeters');
    set(fig, 'Position', [2 2 fig_width fig_height]); % position on screen
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [fig_width fig_height]);
    set(fig, 'PaperPositionMode', 'manual');
    set(fig, 'PaperPosition', [0 0 fig_width fig_height]);

    % Save at 600 dpi
    print(fig, filename, '-dpng', '-r600');
end
