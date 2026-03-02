%% Debug-only visualization runner for slope_stability_3D_hetero_seepage_SSR_comsol
% Loads a saved snapshot of postprocessing inputs and runs visualization only.

snapshot_file = fullfile('tmp', 'slope_stability_3D_hetero_seepage_SSR_comsol_vis_inputs.mat');
if ~isfile(snapshot_file)
    error('Snapshot file not found: %s', snapshot_file);
end

out_dir = fullfile('tmp', 'vis_debug_outputs');
if ~isfolder(out_dir)
    mkdir(out_dir);
end

close all;
load(snapshot_file);
fprintf('Loaded visualization snapshot: %s\n', snapshot_file);

if direct_on && ~isempty(U2)
    fprintf('\nDirect continuation postprocessing\n');
    VIZ.draw_mesh_3D(coord, surf);
    VIZ.draw_quantity_3D_old(coord, surf, zeros(size(coord)), pw, x1, x2, x3, y1, y2, z);
    VIZ.plot_displacements_3D(U2, coord, elem, surf, 0.05 * max(abs(coord(:))) / max(abs(U2(:))));
    VIZ.plot_deviatoric_strain_3D(U2, coord, elem, surf, B);
    figure; hold on; box on; grid on;
    plot(omega_hist2, lambda_hist2, '-o');
    title('Direct continuation method', 'Interpreter', 'latex');
    xlabel('variable - $\xi$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end

if indirect_on && ~isempty(U3)
    fprintf('\nIndirect continuation postprocessing\n');
    VIZ.draw_mesh_3D(coord, surf);
    drawnow;
    pause(0.2);

    VIZ.plot_pore_pressure_3D(pw, coord, surf);
    drawnow;
    pause(0.2);

    VIZ.plot_displacements_3D(U3, coord, elem, surf, 0.05 * max(abs(coord(:))) / max(abs(U3(:))));
    drawnow;
    pause(0.2);

    VIZ.plot_deviatoric_strain_3D(U3, coord, elem, surf, B);
    cl = caxis(gca);
    clim = [cl(1), 0.25 * cl(2)];
    caxis(clim);
    drawnow;
    pause(0.2);

    plane_vals = {[], [35], [1e-16, 21.6506]};
    [figs, info] = VIZ.plot_deviatoric_norm_slices(B, U3, elem, coord, Xi, surf, plane_vals, 1); %#ok<NASGU,ASGLU>

    figure; hold on; box on; grid on;
    plot(omega_hist3, lambda_hist3, '-o');
    title('Indirect continuation method', 'Interpreter', 'latex');
    xlabel('control variable - $\omega$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end

drawnow;
fig_handles = sort(get(0, 'children'));
saved_count = 0;
for k = 1:numel(fig_handles)
    f = fig_handles(k);
    if ~(ishghandle(f) && strcmp(get(f, 'type'), 'figure'))
        continue;
    end
    out_png = fullfile(out_dir, sprintf('fig_%02d.png', k));
    try
        print(f, out_png, '-dpng', '-r200');
        saved_count = saved_count + 1;
    catch
        % Skip invalid/stale graphics handles.
    end
end
fprintf('\nVisualization debug completed. Saved %d figure(s) to: %s\n', saved_count, out_dir);
