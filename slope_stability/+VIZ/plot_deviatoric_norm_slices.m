function [figs, slice_info] = plot_deviatoric_norm_slices(B, U, elem, coord, Xi, surf, plane_vals_cell, h, clim)
%PLOT_DEVIATORIC_NORM_SLICES  Plot |dev(E)| on slices for given plane values.
%
% Inputs
%   B, U                 -- standard FE strain-displacement matrix and DOFs
%   elem, coord, Xi      -- mesh (P2 tetra), node coords (3 x N), IPs rule
%   surf                 -- surface/skin (as used in your VIZ.* utilities)
%   plane_vals_cell      -- 1x3 cell: {x_vals, y_vals, z_vals}
%                           e.g., {[], [35 40], []} to slice at y=35 and y=40
%   h                    -- target edge size for 2D triangulation of slice
%   clim                 -- [cmin cmax] or [] for auto
%
% Outputs
%   figs       -- array of figure handles for each produced slice
%   slice_info -- struct array with meta (plane_id, plane_val, TR2, P3, etc.)
%
% Notes
%   • For a constant-x slice, the plot axes are (y, z).
%   • For a constant-y slice, the plot axes are (x, z).
%   • For a constant-z slice, the plot axes are (x, y).
%
%   Title is set like "x = value", "y = value", or "z = value".
%   All text uses LaTeX with a default fontsize of 15.
warnState = warning('off','all'); 
    if nargin < 9, clim = []; end
    if isempty(h), h = 1; end

    % --- style defaults
    fs = 15;
    set(groot, 'defaultAxesFontSize', fs);
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');

    % --- volumetric / deviatoric projector
    iota = [1; 1; 1; 0; 0; 0];
    VOL  = iota * iota.';
    DEV  = diag([1, 1, 1, 1/2, 1/2, 1/2]) - VOL / 3;

    % --- strain at all IPs, deviatoric norm
    E = B * U(:);
    E = reshape(E, 6, []);                 % columns = integration points
    dev_E = DEV * E;
    norm_E = sqrt(max(0, sum(E .* dev_E, 1))).';   % as column vector

    % --- integration point coordinates (3 x nIP)
    transformed_points = ASSEMBLY.integration_points(elem, coord, Xi);

    % --- boundary processing (do once)
    edges_merged = VIZ.compute_bounding_edges(surf, coord, 0);
    faces        = VIZ.detect_faces(edges_merged, coord);

    % --- collect all requested (plane_id, plane_val) pairs
    req = [];
    for pid = 1:3
        vals = plane_vals_cell{pid};
        if ~isempty(vals)
            req = [req; [pid * ones(numel(vals),1), vals(:)]]; %#ok<AGROW>
        end
    end

    figs = gobjects(0);
    slice_info = struct('plane_id',{},'plane_val',{},'TR2',{},'P3',{},'poly3',{},'vals',{});

    % --- iterate over all requested slices
    for k = 1:size(req,1)
        plane_id = req(k,1);
        plane_val = req(k,2);

        % slice, triangulate
        [poly3, ~] = VIZ.slice_by_plane(faces, edges_merged, coord, plane_id, plane_val);
        [TR2, P3]  = VIZ.triangulate_polygon_slice(poly3, h, struct('plot2d', false));

        % axis labels per plane:
        switch plane_id
            case 1 % x = const  -> axes (y,z)
                ttl = sprintf('$x = %.6g$', plane_val);
                xl = '$y$'; yl = '$z$';
            case 2 % y = const  -> axes (x,z)
                ttl = sprintf('$y = %.6g$', plane_val);
                xl = '$x$'; yl = '$z$';
            case 3 % z = const  -> axes (x,y)
                ttl = sprintf('$z = %.6g$', plane_val);
                xl = '$x$'; yl = '$y$';
            otherwise
                error('Invalid plane_id: %d', plane_id);
        end

        opts = struct('title', ttl, ...
                      'xlabel', xl, 'ylabel', yl, ...
                      'clim', clim, 'rotate', false, 'fontsize', fs, ...
                      'method', 'linear', 'colormap', 'parula');

        % plot
        fig = figure('Name', sprintf('Slice %d @ %.6g', plane_id, plane_val));
        [vals, hC, hB, F] = VIZ.interpolate_and_plot_on_slice( ...
            TR2, P3, transformed_points, norm_E', poly3, opts);

        % enforce clim & colormap uniformly (in case opts not applied inside)
        if ~isempty(clim) && isfloat(clim) && numel(clim) == 2
            caxis(clim);
        end
        colormap(parula);
        % if ~isempty(hB) && isgraphics(hB), set(hB, 'Interpret', 'latex'); end

        xlabel(xl, 'Interpreter', 'latex', 'FontSize', fs);
        ylabel(yl, 'Interpreter', 'latex', 'FontSize', fs);
        title(ttl, 'Interpreter', 'latex', 'FontSize', fs);

        % store
        figs(end+1,1) = fig; %#ok<AGROW>
        slice_info(end+1) = struct('plane_id', plane_id, ...
                                   'plane_val', plane_val, ...
                                   'TR2', TR2, 'P3', P3, ...
                                   'poly3', poly3, 'vals', vals); %#ok<AGROW>
        hold off
        drawnow
        pause(0.5)
    end
    warnState = warning('on','all'); 
end
