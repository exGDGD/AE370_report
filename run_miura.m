%% Final Miura report search over the justified integer range
%
% Workflow:
%   1. Coarse mixed-integer search over rows, cols
%   2. Refine the top coarse layouts.
%   3. Generate convergence and local landscape figures for the final best.

projectDirectory = fileparts(mfilename("fullpath"));
addpath(fullfile(projectDirectory, "matlab"));
cd(projectDirectory);

reportDirectory = fullfile(projectDirectory, "result");
if ~exist(reportDirectory, "dir")
    mkdir(reportDirectory);
end

% Set true when search settings or integer ranges changed and old MAT files
% should not be reused.
forceRerun = true;

%% 1. Coarse integer search
coarseOptions.box = [5.2, 3.0, 0.8];
coarseOptions.thickness = 0.06;
coarseOptions.clearance = 0;
coarseOptions.rowRange = [6, 6];
coarseOptions.colRange = [4, 4];
coarseOptions.gammaRangeDegrees = [35, 89.9];
coarseOptions.edgeTiltRange = [0.35, 0.999];
coarseOptions.aspectRatioRange = [0.35, 3.00];
coarseOptions.bisectionIterations = 5;
coarseOptions.pathSamples = 3;
coarseOptions.maxFunctionEvaluations = 80;
coarseOptions.maxIterations = 40;
coarseOptions.optimizerTolerance = 8e-5;
coarseOptions.scaleTolerance = 8e-4;
coarseOptions.optimizer = "fminsearch";
coarseOptions.recordHistory = false;
coarseOptions.outputDirectory = fullfile(reportDirectory, ...
    sprintf("coarse_%dto%d_%dto%d", coarseOptions.rowRange(1), coarseOptions.rowRange(2), ...
    coarseOptions.colRange(1), coarseOptions.colRange(2)));
coarseOptions.showFigure = false;

coarseResultPath = fullfile(coarseOptions.outputDirectory, "best_result.mat");
if exist(coarseResultPath, "file") && ~forceRerun
    fprintf("\nReusing coarse search.\n");
    loaded = load(coarseResultPath, "best", "layouts", "options");
    coarseResults.best = loaded.best;
    coarseResults.layouts = loaded.layouts;
    coarseResults.options = loaded.options;
else
    fprintf("\nRunning coarse 4..10 mixed-integer search...\n");
    coarseResults = optimizeMiura(coarseOptions);
end

coarseTable = readtable(fullfile(coarseOptions.outputDirectory, "rankings.csv"));
topK = min(6, height(coarseTable));
topLayouts = unique(coarseTable{1:topK, ["rows", "cols"]}, "rows", "stable");

%% 2. Refine top coarse layouts
refinedResults = struct([]);
for index = 1:size(topLayouts, 1)
    rows = topLayouts(index, 1);
    cols = topLayouts(index, 2);

    refineOptions = coarseOptions;
    refineOptions.rowRange = [rows, rows];
    refineOptions.colRange = [cols, cols];
    refineOptions.bisectionIterations = 9;
    refineOptions.pathSamples = 7;
    refineOptions.maxFunctionEvaluations = 220;
    refineOptions.maxIterations = 80;
    refineOptions.optimizerTolerance = 1e-5;
    refineOptions.scaleTolerance = 1e-4;
    refineOptions.recordHistory = true;
    refineOptions.outputDirectory = fullfile(reportDirectory, sprintf("refined_%dx%d", rows, cols));
    refineOptions.showFigure = false;

    resultPath = fullfile(refineOptions.outputDirectory, "best_result.mat");
    if exist(resultPath, "file") && ~forceRerun
        fprintf("\nReusing refined layout %d x %d.\n", rows, cols);
        loaded = load(resultPath, "best", "layouts", "options", "optimizationHistories");
        runResult.best = loaded.best;
        runResult.layouts = loaded.layouts;
        runResult.options = loaded.options;
        runResult.optimizationHistories = loaded.optimizationHistories;
    else
        fprintf("\nRefining layout %d x %d...\n", rows, cols);
        runResult = optimizeMiura(refineOptions);
    end

    refinedResults(index).layout = [rows, cols]; %#ok<SAGROW>
    refinedResults(index).best = runResult.best; %#ok<SAGROW>
end

densityValues = arrayfun(@(item) item.best.areaDensity, refinedResults);
[~, order] = sort(densityValues, "descend");
refinedResults = refinedResults(order);

summaryTable = table();
for index = 1:numel(refinedResults)
    summaryTable = [summaryTable; struct2table(refinedResults(index).best)]; %#ok<AGROW>
end
writetable(summaryTable, fullfile(reportDirectory, "final_refined_rankings.csv"));
save(fullfile(reportDirectory, "final_refined_results.mat"), ...
    "coarseResults", "refinedResults", "summaryTable");

finalBest = refinedResults(1).best;

%% 3. Local objective landscape around the final best
landscapeOptions = coarseOptions;
landscapeOptions.bisectionIterations = 8;
landscapeOptions.scaleTolerance = 2e-4;

gammaValues = linspace(max(35, finalBest.gammaDegrees - 7), ...
    min(78, finalBest.gammaDegrees + 7), 31);
tiltValues = linspace(max(0.35, finalBest.normalizedEdgeTilt - 0.06), ...
    min(0.985, finalBest.normalizedEdgeTilt + 0.02), 31);
density = nan(numel(tiltValues), numel(gammaValues));

fprintf("\nSampling local landscape around final best %d x %d...\n", ...
    finalBest.rows, finalBest.cols);
for tiltIndex = 1:numel(tiltValues)
    for gammaIndex = 1:numel(gammaValues)
        candidate = evaluateMiuraDesign(finalBest.rows, finalBest.cols, ...
            gammaValues(gammaIndex), tiltValues(tiltIndex), ...
            finalBest.aspectRatio, landscapeOptions);
        if candidate.feasible
            density(tiltIndex, gammaIndex) = candidate.areaDensity;
        end
    end
end

gammaVariableNames = matlab.lang.makeValidName(compose("gamma_%0.3f", gammaValues));
landscapeTable = array2table(density, "VariableNames", gammaVariableNames);
landscapeTable.edgeTilt = tiltValues(:);
landscapeTable = movevars(landscapeTable, "edgeTilt", "Before", 1);
writetable(landscapeTable, fullfile(reportDirectory, "final_best_local_landscape.csv"));

figureHandle = figure("Color", "white", "Name", "Final best local landscape", "Visible", "off");
contourf(gammaValues, tiltValues, density, 18, "LineColor", "none");
hold on;
plot(finalBest.gammaDegrees, finalBest.normalizedEdgeTilt, "wo", ...
    "MarkerFaceColor", "k", "MarkerSize", 6);
colorbar;
grid on;
xlabel("gamma (deg)");
ylabel("normalized edge tilt");
title(sprintf("Local objective landscape near final best %d x %d", ...
    finalBest.rows, finalBest.cols));
exportgraphics(figureHandle, fullfile(reportDirectory, "final_best_local_landscape.png"), ...
    "Resolution", 220);
close(figureHandle);

fprintf("\nFinal report search complete.\n");
fprintf("Best refined layout: %d x %d\n", finalBest.rows, finalBest.cols);
fprintf("Best active area / volume: %.8f\n", finalBest.areaDensity);
fprintf("outputs written to:\n%s\n", reportDirectory);
