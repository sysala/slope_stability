function [figs, slice_info] = plot_deviatoric_norm_slices(B, U, elem, coord, Xi, surf, plane_vals_cell, h, clim)
%PLOT_DEVIATORIC_NORM_SLICES  Plot |dev(E)| on slices for given plane values.
%
% Inputs
%   B, U                 -- standard FE strain-displacement matrix and DOFs
%   elem, coord          -- mesh (P2 tetra), node coords (3 x N), IPs rule
%   Xi                   -- local coordinates of quadrature points, size: (3, n_q)
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
    local_set_default('defaultAxesFontSize', fs);
    local_set_default('defaultTextInterpreter', 'latex');
    local_set_default('defaultAxesTickLabelInterpreter', 'latex');
    local_set_default('defaultLegendInterpreter', 'latex');

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

    figs = [];
    slice_info = struct('plane_id',{},'plane_val',{},'TR2',{},'P3',{},'poly3',{},'vals',{});
    cmap_name = local_colormap_name();

    % Octave fallback: avoid expensive 3D scattered interpolation and
    % triangulation classes not available on many Octave builds.
    if exist('OCTAVE_VERSION', 'builtin') ~= 0
        for k = 1:size(req,1)
            plane_id = req(k,1);
            plane_val = req(k,2);

            [poly3, ~] = VIZ.slice_by_plane(faces, edges_merged, coord, plane_id, plane_val);
            if isempty(poly3)
                continue;
            end
            [TR2, P3] = VIZ.triangulate_polygon_slice(poly3, h, struct('plot2d', false)); %#ok<NASGU>
            [ttl, xl, yl, free_axes] = local_plane_labels(plane_id, plane_val);

            % Sample integration points close to the requested plane.
            d = abs(transformed_points(plane_id,:) - plane_val);
            slice_tol = max(1e-6, 0.75 * h);
            idx = d <= slice_tol;
            if ~any(idx)
                [~, ord] = sort(d, 'ascend');
                take = min(10000, numel(ord));
                idx = false(size(d));
                idx(ord(1:take)) = true;
            end

            pts2 = transformed_points(free_axes, idx).';
            vals = norm_E(idx);
            if size(pts2,1) > 30000
                rng(0);
                take = randperm(size(pts2,1), 30000);
                pts2 = pts2(take, :);
                vals = vals(take);
            end

            fig = figure('Name', sprintf('Slice %d @ %.6g (Octave)', plane_id, plane_val));
            CL = TR2.ConnectivityList;
            V2 = TR2.Points;

            % Interpolate sampled values onto triangulation nodes.
            vals_on_nodes = local_interp2d(pts2, vals, V2);

            patch('Faces', CL, 'Vertices', V2, ...
                'FaceVertexCData', double(vals_on_nodes), ...
                'FaceColor', 'interp', 'EdgeColor', 'none');
            hold on;
            poly2 = poly3(free_axes,:).';
            plot(poly2(:,1), poly2(:,2), 'k-', 'LineWidth', 1.2);
            hold off;

            axis equal;
            axis tight;
            grid on;
            box on;
            if ~isempty(clim) && isfloat(clim) && numel(clim) == 2
                caxis(clim);
            end
            colormap(cmap_name);
            cb = colorbar; %#ok<NASGU>
            xlabel(xl, 'Interpreter', 'latex', 'FontSize', fs);
            ylabel(yl, 'Interpreter', 'latex', 'FontSize', fs);
            title(ttl, 'Interpreter', 'latex', 'FontSize', fs);

            figs(end+1,1) = fig; %#ok<AGROW>
            slice_info(end+1) = struct('plane_id', plane_id, ...
                                       'plane_val', plane_val, ...
                                       'TR2', TR2, 'P3', P3, ...
                                       'poly3', poly3, 'vals', vals_on_nodes); %#ok<AGROW>
            drawnow;
        end
        warning(warnState);
        return;
    end

    % --- iterate over all requested slices
    for k = 1:size(req,1)
        plane_id = req(k,1);
        plane_val = req(k,2);

        % slice, triangulate
        [poly3, ~] = VIZ.slice_by_plane(faces, edges_merged, coord, plane_id, plane_val);
        [TR2, P3]  = VIZ.triangulate_polygon_slice(poly3, h, struct('plot2d', false));

        [ttl, xl, yl] = local_plane_labels(plane_id, plane_val);

        opts = struct('title', ttl, ...
                      'xlabel', xl, 'ylabel', yl, ...
                      'clim', clim, 'rotate', false, 'fontsize', fs, ...
                      'method', 'linear', 'colormap', cmap_name);

        % plot
        fig = figure('Name', sprintf('Slice %d @ %.6g', plane_id, plane_val));
        [vals, hC, hB, F] = VIZ.interpolate_and_plot_on_slice( ...
            TR2, P3, transformed_points, norm_E', poly3, opts);

        % enforce clim & colormap uniformly (in case opts not applied inside)
        if ~isempty(clim) && isfloat(clim) && numel(clim) == 2
            caxis(clim);
        end
        colormap(cmap_name);
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

function local_set_default(prop_name, prop_value)
try
    set(0, prop_name, prop_value);
catch
    % Skip unsupported graphics defaults on this runtime (e.g., Octave).
end
end

function [ttl, xl, yl, free_axes] = local_plane_labels(plane_id, plane_val)
switch plane_id
    case 1 % x = const  -> axes (y,z)
        ttl = sprintf('$x = %.6g$', plane_val);
        xl = '$y$'; yl = '$z$';
        free_axes = [2 3];
    case 2 % y = const  -> axes (x,z)
        ttl = sprintf('$y = %.6g$', plane_val);
        xl = '$x$'; yl = '$z$';
        free_axes = [1 3];
    case 3 % z = const  -> axes (x,y)
        ttl = sprintf('$z = %.6g$', plane_val);
        xl = '$x$'; yl = '$y$';
        free_axes = [1 2];
    otherwise
        error('Invalid plane_id: %d', plane_id);
end
end

function cmap_name = local_colormap_name()
if exist('parula', 'file') == 2 || exist('parula', 'builtin') == 5
    cmap_name = 'parula';
else
    cmap_name = 'jet';
end
end

function vals_q = local_interp2d(pts2, vals, q2)
% Robust 2D interpolation with nearest fallback.
if size(pts2,1) < 3
    vals_q = repmat(mean(vals), size(q2,1), 1);
    return;
end
vals_q = griddata(pts2(:,1), pts2(:,2), vals, q2(:,1), q2(:,2), 'linear');
bad = isnan(vals_q);
if any(bad)
    vals_q(bad) = griddata(pts2(:,1), pts2(:,2), vals, q2(bad,1), q2(bad,2), 'nearest');
end
if any(isnan(vals_q))
    vals_q(isnan(vals_q)) = mean(vals);
end
end
