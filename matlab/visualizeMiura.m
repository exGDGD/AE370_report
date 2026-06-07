function figureHandle = visualizeMiura(model, outputPath)
%VISUALIZEMIURA Plot isometric, top, and side views of a Miura model.

arguments
    model struct
    outputPath string = ""
end

figureHandle = figure("Color", "white", "Name", "Finite-thickness Miura-ori");
layout = tiledlayout(1, 3, "TileSpacing", "compact", "Padding", "compact");
title(layout, sprintf("Miura-ori: active area %.4f, bbox %.3f x %.3f x %.3f", ...
    model.activeArea, model.bboxSize(1), model.bboxSize(2), model.bboxSize(3)));

views = {[35, 26], [0, 90], [0, 0]};
titles = {"Isometric folded view", "Top view", "Side view"};
colors = [
    0.41, 0.65, 0.83
    0.94, 0.64, 0.37
    0.53, 0.74, 0.47
    0.77, 0.59, 0.83
];

for viewIndex = 1:3
    nexttile;
    hold on;
    for facetIndex = 1:numel(model.facets)
        facet = model.facets(facetIndex);
        if size(facet.activeVertices, 1) < 3
            continue
        end
        row = facet.index(1);
        column = facet.index(2);
        color = colors(mod(row + column - 2, size(colors, 1)) + 1, :);
        patch("Vertices", facet.activeVertices, "Faces", 1:size(facet.activeVertices, 1), ...
            "FaceColor", color, "FaceAlpha", 0.82, "EdgeColor", [0.22, 0.25, 0.30]);
    end
    mesh = model.mesh;
    for row = 1:size(mesh, 1)
        plot3(squeeze(mesh(row, :, 1)), squeeze(mesh(row, :, 2)), squeeze(mesh(row, :, 3)), ...
            "Color", [0.16, 0.20, 0.27], "LineWidth", 0.7);
    end
    for column = 1:size(mesh, 2)
        plot3(squeeze(mesh(:, column, 1)), squeeze(mesh(:, column, 2)), squeeze(mesh(:, column, 3)), ...
            "Color", [0.16, 0.20, 0.27], "LineWidth", 0.7);
    end
    axis equal;
    grid on;
    xlabel("x"); ylabel("y"); zlabel("z");
    title(titles{viewIndex});
    view(views{viewIndex});
end

if strlength(outputPath) > 0
    exportgraphics(figureHandle, outputPath, "Resolution", 180);
end
end
