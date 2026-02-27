filename = "convergence_hydro";
set(gcf, 'Position', [100, 100, 1050, 600]);
exportgraphics(gcf, filename + ".pdf", "ContentType", "vector", "BackgroundColor", "none");
savefig(gcf, filename + ".fig");
