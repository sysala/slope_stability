function [vals_on_nodes, hPatch, hBnd, F_out] = ...
    interpolate_and_plot_on_slice(TR2, P3, transformed_points, norm_E, poly3, opts)
% INTERPOLATE_AND_PLOT_ON_SLICE_2D
% Map scattered 3D data to a 2D slice mesh and plot in 2D (no grid shown).
%
% [vals_on_nodes, hPatch, hBnd, F_out] = interpolate_and_plot_on_slice_2d( ...
%     TR2, P3, transformed_points, norm_E, poly3, opts)
%
% Required:
%   TR2                : triangulation (2D) of slice (Points Nx2, ConnectivityList Mx3)
%   P3                 : (3 x N) 3D coordinates of TR2 nodes (same order as TR2.Points)
%   transformed_points : (3 x M) scattered 3D samples
%   norm_E             : (1 x M) values at transformed_points
%   poly3              : (3 x K) closed slice boundary (first == last)
%
% opts (all optional):
%   .fast_method   : 'scattered' (default) | 'idw'
%   .method        : 'natural' | 'linear' | 'nearest'  (for scattered; default 'natural')
%   .F             : prebuilt scatteredInterpolant to reuse
%   .max_points    : subsample scattered data before building interpolant ([])
%   .single_precision : false (default) | true
%   .title         : char, LaTeX title (default '')
%   .xlabel        : char, LaTeX x-label (default '')
%   .ylabel        : char, LaTeX y-label (default '')
%   .clim          : [] or [vmin vmax] (default [])
%   .rotate        : false (default) | true   % swap axes for display
%   .fontsize      : 15 (default)
%   .colormap      : 'parula' (default) or other
%   .idw_k         : 16  (neighbors for IDW)
%   .idw_power     : 2   (IDW power)
%
% Outputs:
%   vals_on_nodes : (N x 1) interpolated values at TR2 nodes
%   hPatch        : handle to colored patch (no edges)
%   hBnd          : handle to boundary line (red)
%   F_out         : interpolant (reuse in next calls via opts.F)

if nargin < 6, opts = struct(); end
fast_method = get_opt(opts,'fast_method','scattered');
method      = get_opt(opts,'method','natural');
F_out       = [];
max_points  = get_opt(opts,'max_points',[]);
use_single  = get_opt(opts,'single_precision',false);

ttl         = get_opt(opts,'title','');
xl          = get_opt(opts,'xlabel','');
yl          = get_opt(opts,'ylabel','');
climv       = get_opt(opts,'clim',[]);
rotateFlag  = get_opt(opts,'rotate',false);
fs          = get_opt(opts,'fontsize',15);
cm          = get_opt(opts,'colormap','parula');

idw_k       = get_opt(opts,'idw_k',16);
idw_p       = get_opt(opts,'idw_power',2);

% --- checks
Nnodes = size(P3,2);
if size(TR2.Points,1) ~= Nnodes
    error('P3 has %d nodes, TR2.Points has %d rows. They must match.', Nnodes, size(TR2.Points,1));
end
if size(transformed_points,1) ~= 3
    error('transformed_points must be 3 x M.');
end
if size(norm_E,1) ~= 1 || size(norm_E,2) ~= size(transformed_points,2)
    error('norm_E must be 1 x M, matching transformed_points.');
end

% --- targets (slice nodes)
Xn = P3(1,:).'; Yn = P3(2,:).'; Zn = P3(3,:).';
if use_single, Xn=single(Xn); Yn=single(Yn); Zn=single(Zn); end

% --- interpolate ONCE (with 'nearest' extrapolation) or IDW
switch lower(fast_method)
    case 'scattered'
        if isfield(opts,'F') && isa(opts.F,'scatteredInterpolant')
            F = opts.F;
        else
            M = size(transformed_points,2);
            idx = 1:M;
            if ~isempty(max_points) && max_points < M
                rng(0); idx = randperm(M, max_points);
            end
            Xs = transformed_points(1,idx).';
            Ys = transformed_points(2,idx).';
            Zs = transformed_points(3,idx).';
            Vs = norm_E(1,idx).';
            if use_single, Xs=single(Xs); Ys=single(Ys); Zs=single(Zs); Vs=single(Vs); end
            F = scatteredInterpolant(Xs, Ys, Zs, Vs, method, 'nearest');
        end
        vals_on_nodes = F(Xn, Yn, Zn);
        F_out = F;

    case 'idw'
        M = size(transformed_points,2);
        idx = 1:M;
        if ~isempty(max_points) && max_points < M
            rng(0); idx = randperm(M, max_points);
        end
        Xs = transformed_points(1,idx).';
        Ys = transformed_points(2,idx).';
        Zs = transformed_points(3,idx).';
        Vs = norm_E(1,idx).';
        if use_single, Xs=single(Xs); Ys=single(Ys); Zs=single(Zs); Vs=single(Vs); end
        try
            ns = createns([Xs Ys Zs],'NSMethod','kdtree'); %#ok<CREATENS>
            [nbrIdx, d] = knnsearch(ns,[Xn Yn Zn],'K',idw_k); %#ok<KNNSearch>
        catch
            nbrIdx = dsearchn([Xs Ys Zs],[Xn Yn Zn]); %#ok<DSEARCHN>
            d = ones(size(nbrIdx)); idw_k = 1;
        end
        w = 1 ./ max(d, eps).^idw_p; w = w ./ sum(w,2);
        vals_on_nodes = sum(w .* Vs(nbrIdx), 2);
    otherwise
        error('Unknown fast_method: %s', fast_method);
end

% --- 2D vertices (optionally swap x/y for display)
V2 = TR2.Points;            % Nx2
if rotateFlag, V2 = V2(:,[2 1]); end

% --- boundary from poly3 (project the two free axes)
poly3u = poly3;                                 % closed expected
if any(poly3u(:,1) ~= poly3u(:,end)), poly3u(:,end+1) = poly3u(:,1); end
srow = [std(poly3u(1,:)), std(poly3u(2,:)), std(poly3u(3,:))];
[~, fixed_ax] = min(srow);
free_axes = setdiff(1:3,fixed_ax);
poly2 = poly3u(free_axes,:).';
if rotateFlag, poly2 = poly2(:,[2 1]); end

% --- draw (no grid)
CL = TR2.ConnectivityList;
hPatch = patch('Faces', CL, 'Vertices', V2, ...
    'FaceVertexCData', double(vals_on_nodes), ...
    'FaceColor','interp', 'EdgeColor','none', 'FaceAlpha', 1.0);
hold on;
hBnd = plot(poly2(:,1), poly2(:,2), 'k-', 'LineWidth', 1);
hold off;

axis equal; axis tight; box on; grid on;
colormap(cm);
cb = colorbar; set(cb,'TickLabelInterpreter','latex','FontSize',fs);
if ~isempty(climv), caxis(climv); end

ax = gca;
set(ax,'TickLabelInterpreter','latex','FontSize',fs);
xlabel(ax, xl, 'Interpreter','latex','FontSize',fs);
ylabel(ax, yl, 'Interpreter','latex','FontSize',fs);
if ~isempty(ttl), title(ax, ttl, 'Interpreter','latex','FontSize',fs); end
end

function val = get_opt(s, name, def)
if isfield(s,name) && ~isempty(s.(name)), val = s.(name); else, val = def; end
end
