function results = optimizeMiura(options)
%OPTIMIZEMIURA Optimize active Miura area inside a launcher bounding box.
% Integer row and column counts are enumerated. For each layout, a multi-start
% continuous optimization solves for gamma, normalized edge tilt, and a:b.
% Scale is eliminated by a one-dimensional bisection against the bounding box.

if nargin < 1
    options = struct();
end
options = setDefault(options, "box", [5.2, 3.0, 0.8]);
options = setDefault(options, "thickness", 0.06);
options = setDefault(options, "clearance", 0);
options = setDefault(options, "rowRange", [2, 5]);
options = setDefault(options, "colRange", [2, 7]);
options = setDefault(options, "bisectionIterations", 10);
options = setDefault(options, "pathSamples", 11);
options = setDefault(options, "outputDirectory", "optimization_matlab");
options = setDefault(options, "showFigure", true);
options = setDefault(options, "maxFunctionEvaluations", 220);
options = setDefault(options, "maxIterations", 80);
options = setDefault(options, "optimizerTolerance", 1e-5);
options = setDefault(options, "scaleTolerance", 1e-4);
options = setDefault(options, "optimizer", "fminsearch");
options = setDefault(options, "recordHistory", true);
options = setDefault(options, "gammaRangeDegrees", [35, 89.9]);
options = setDefault(options, "edgeTiltRange", [0.35, 0.999]);
options = setDefault(options, "aspectRatioRange", [0.45, 2.20]);
options = setDefault(options, "skipFinalCollisionChecks", false);
options.optimizer = lower(string(options.optimizer));
mustBePositive(options.box); mustBePositive(options.thickness);
mustBeNonnegative(options.clearance);
mustBeInteger(options.rowRange); mustBePositive(options.rowRange);
mustBeInteger(options.colRange); mustBePositive(options.colRange);
if ~any(options.optimizer == ["auto", "sqp", "fminsearch"])
    error("optimizeMiura:InvalidOptimizer", ...
        "options.optimizer must be 'auto', 'sqp', or 'fminsearch'.");
end

if options.rowRange(1) > options.rowRange(2) || options.colRange(1) > options.colRange(2)
    error("optimizeMiura:InvalidRange", "Layout ranges must be increasing.");
end

function options = setDefault(options, name, value)
%SETDEFAULT Add an option value only when the caller did not provide it.
if ~isfield(options, name)
    options.(name) = value;
end
end

% Enumerate integer row/column layouts. Each layout has its own continuous
% optimization problem for gamma, normalized edge tilt, and aspect ratio.
layouts = struct([]);
optimizationHistories = struct([]);
layoutIndex = 0;
for rows = options.rowRange(1):options.rowRange(2)
    for cols = options.colRange(1):options.colRange(2)
        layoutIndex = layoutIndex + 1;
        [candidate, layoutHistories] = optimizeLayout(rows, cols, options);
        if isempty(layouts)
            layouts = candidate;
        else
            layouts(layoutIndex) = candidate; %#ok<AGROW>
        end
        optimizationHistories = appendHistories(optimizationHistories, layoutHistories);
        fprintf("layout %d x %d: density %.6f\n", rows, cols, layouts(layoutIndex).areaDensity);
    end
end

% Sort all integer layouts by objective value and then run the expensive path
% collision gate in that order.
[~, order] = sort([layouts.areaDensity], "descend");
layouts = layouts(order);
for index = 1:numel(layouts)
    if options.skipFinalCollisionChecks
        layouts(index).pathCollisionFree = layouts(index).feasible;
    else
        layouts(index).pathCollisionFree = verifyPath(layouts(index), options);
    end
end

valid = find([layouts.pathCollisionFree], 1, "first");
if isempty(valid)
    error("optimizeMiura:NoFeasibleCandidate", "No path-valid candidate was found.");
end
best = layouts(valid);

% Rebuild the final model with expensive verification enabled so the reported
% residuals and collision count are from the exact final geometry.
model = makeModel(best, options);
best.collisionCount = model.collisionCount;
best.rigidFoldResidual = model.rigidFoldResidual;
best.hingeRotationResidual = model.hingeRotationResidual;
layouts(valid) = best;

if ~exist(options.outputDirectory, "dir")
    mkdir(options.outputDirectory);
end
writetable(struct2table(layouts), fullfile(options.outputDirectory, "rankings.csv"));
historyTable = historiesToTable(optimizationHistories);
if ~isempty(historyTable)
    writetable(historyTable, fullfile(options.outputDirectory, "objective_history.csv"));
end
save(fullfile(options.outputDirectory, "best_result.mat"), "best", "layouts", "options", "optimizationHistories");
writeSummary(best, options, fullfile(options.outputDirectory, "best_summary.txt"));
if options.showFigure
    visualizeMiura(model, fullfile(options.outputDirectory, "best_miura.png"));
end
plotRankings(layouts, fullfile(options.outputDirectory, "search_overview.png"));
plotConvergence(optimizationHistories, best, fullfile(options.outputDirectory, "objective_convergence.png"));

fprintf("\nBest verified candidate\n");
fprintf("layout:                  %d x %d\n", best.rows, best.cols);
fprintf("gamma:                   %.6f deg\n", best.gammaDegrees);
fprintf("normalized edge tilt:    %.8f\n", best.normalizedEdgeTilt);
fprintf("aspect ratio a:b:        %.8f\n", best.aspectRatio);
fprintf("a, b:                    %.8f, %.8f\n", best.a, best.b);
fprintf("setback q:               %.8f\n", best.setback);
fprintf("active area:             %.8f\n", best.activeArea);
fprintf("active area / volume:    %.8f\n", best.areaDensity);
fprintf("bbox:                    %.5f x %.5f x %.5f\n", best.bboxSize);
fprintf("hinge dihedral range:    %.4f to %.4f deg\n", best.hingeDihedralRangeDegrees);
fprintf("rigid-fold residual:     %.3e\n", best.rigidFoldResidual);
fprintf("Rodrigues residual:      %.3e\n", best.hingeRotationResidual);
fprintf("non-adjacent collisions: %d\n", best.collisionCount);

results.best = best;
results.layouts = layouts;
results.options = options;
results.optimizationHistories = optimizationHistories;
end

function [candidate, histories] = optimizeLayout(rows, cols, options)
%OPTIMIZELAYOUT Solve the continuous subproblem for one integer layout.
% Variables are [gamma, normalizedEdgeTilt, log(a/b)]. The log aspect ratio
% makes stretching and shrinking symmetric in the optimizer.
lower = [deg2rad(options.gammaRangeDegrees(1)), options.edgeTiltRange(1), log(options.aspectRatioRange(1))];
upper = [deg2rad(options.gammaRangeDegrees(2)), options.edgeTiltRange(2), log(options.aspectRatioRange(2))];

% Multi-start points reduce the chance of getting stuck in a poor local
% optimum. The last two starts intentionally bias toward deeply folded states.
starts = [
    deg2rad(45), 0.65, log(0.70)
    deg2rad(55), 0.82, log(1.00)
    deg2rad(65), 0.93, log(1.35)
    deg2rad(75), 0.965, log(0.75)
    deg2rad(85), 0.985, log(1.20)
    deg2rad(58), 0.975, log(0.55)
    deg2rad(60), 0.990, log(0.85)
    deg2rad(65), 0.995, log(1.00)
    deg2rad(72), 0.990, log(1.55)
    deg2rad(88), 0.995, log(0.90)
];
if isfield(options, "customStarts")
    starts = options.customStarts;
end

bestScore = inf;
bestPoint = starts(1, :);
histories = emptyHistory();
hasFmincon = exist("fmincon", "file") == 2;
if options.optimizer == "sqp" && ~hasFmincon
    error("optimizeMiura:MissingOptimizationToolbox", ...
        "SQP requires fmincon from MATLAB Optimization Toolbox. Set options.optimizer='fminsearch' or install Optimization Toolbox.");
end
useSqp = (options.optimizer == "sqp") || (options.optimizer == "auto" && hasFmincon);
for startIndex = 1:size(starts, 1)
    start = starts(startIndex, :);
    history = emptyHistory();
    history.rows = rows;
    history.cols = cols;
    history.startIndex = startIndex;
    if useSqp
        % Use constrained SQP when Optimization Toolbox is available.
        history.algorithm = "sqp";
        solverOptions = optimoptions("fmincon", "Display", "none", ...
            "Algorithm", "sqp", "MaxFunctionEvaluations", options.maxFunctionEvaluations, ...
            "MaxIterations", options.maxIterations, ...
            "StepTolerance", options.optimizerTolerance, ...
            "FunctionTolerance", options.optimizerTolerance, ...
            "OutputFcn", @recordIteration);
        point = fmincon(@objective, start, [], [], [], [], lower, upper, [], solverOptions);
    else
        % Base MATLAB fallback: map unconstrained variables into the bounded
        % design space and call Nelder-Mead through fminsearch.
        history.algorithm = "fminsearch";
        unconstrainedStart = inverseLogistic(start, lower, upper);
        solverOptions = optimset("Display", "off", ...
            "MaxFunEvals", options.maxFunctionEvaluations, ...
            "MaxIter", options.maxIterations, ...
            "TolX", options.optimizerTolerance, ...
            "TolFun", options.optimizerTolerance, ...
            "OutputFcn", @recordIteration);
        unconstrainedPoint = fminsearch(@(value) objective(logistic(value, lower, upper)), ...
            unconstrainedStart, solverOptions);
        point = logistic(unconstrainedPoint, lower, upper);
    end
    score = objective(point);
    history.finalGammaDegrees = rad2deg(point(1));
    history.finalNormalizedEdgeTilt = point(2);
    history.finalAspectRatio = exp(point(3));
    history.finalAreaDensity = -score;
    history.finalScore = score;
    histories = appendHistories(histories, history);
    if score < bestScore
        bestScore = score;
        bestPoint = point;
    end
end

candidate = evaluatePoint(rows, cols, bestPoint, options);

    function score = objective(point)
        % fminsearch/fmincon minimize, so negate the active-area density.
        evaluated = evaluatePoint(rows, cols, point, options);
        if ~evaluated.feasible
            score = 1e6;
        else
            score = -evaluated.areaDensity;
        end
    end

    function stop = recordIteration(~, optimValues, state)
        %RECORDITERATION Collect objective values for convergence plots.
        stop = false;
        if ~options.recordHistory || strcmp(state, "done") || ~isfield(optimValues, "fval")
            return
        end
        if ~isempty(history.score) && history.score(end) == optimValues.fval
            if ~isfield(optimValues, "funccount") || history.functionCount(end) == optimValues.funccount
                return
            end
        end
        history.iteration(end + 1, 1) = optimValues.iteration;
        if isfield(optimValues, "funccount")
            history.functionCount(end + 1, 1) = optimValues.funccount;
        else
            history.functionCount(end + 1, 1) = numel(history.iteration);
        end
        history.score(end + 1, 1) = optimValues.fval;
        history.areaDensity(end + 1, 1) = -optimValues.fval;
    end
end

function candidate = evaluatePoint(rows, cols, point, options)
%EVALUATEPOINT Evaluate one continuous design point for a fixed layout.
% Repeated calls use a fast model: no SAT and no residual checks.
gamma = point(1);
normalizedEdgeTilt = point(2);
aspectRatio = exp(point(3));
scale = maximizeScale(rows, cols, gamma, normalizedEdgeTilt, aspectRatio, options);

candidate = blankCandidate(rows, cols);
if ~isfinite(scale)
    return
end
parameters = makeParameters(rows, cols, scale * aspectRatio, scale, gamma, normalizedEdgeTilt, options, false, false);
model = buildMiuraModel(parameters);

% A candidate is feasible during optimization if it fits the launcher and every
% facet retains nonzero active area after setback clipping.
candidate.feasible = all(model.bboxSize <= options.box + 1e-9) && all([model.facets.activeArea] > 0);
candidate.rows = rows;
candidate.cols = cols;
candidate.gammaDegrees = rad2deg(gamma);
candidate.normalizedEdgeTilt = normalizedEdgeTilt;
candidate.aspectRatio = aspectRatio;
candidate.scale = scale;
candidate.a = scale * aspectRatio;
candidate.b = scale;
candidate.setback = model.setback;
candidate.activeArea = model.activeArea;
candidate.idealArea = model.idealArea;
candidate.activeAreaRatio = model.activeAreaRatio;
candidate.areaDensity = model.activeArea / prod(options.box);
candidate.bboxSize = model.bboxSize;
candidate.minimumOpenAngleDegrees = rad2deg(model.minimumOpenAngle);
candidate.hingeDihedralRangeDegrees = rad2deg(model.hingeDihedralRange);
candidate.rigidFoldResidual = model.rigidFoldResidual;
candidate.hingeRotationResidual = model.hingeRotationResidual;
candidate.collisionCount = model.collisionCount;
candidate.pathCollisionFree = false;
end

function scale = maximizeScale(rows, cols, gamma, normalizedEdgeTilt, aspectRatio, options)
%MAXIMIZESCALE Find the largest scale that fits the launcher box.
% This removes scale from the optimizer variables and turns it into a
% one-dimensional bisection problem.
probe = makeParameters(rows, cols, aspectRatio, 1, gamma, normalizedEdgeTilt, options, false, false);
probeModel = buildMiuraModel(probe);

% Minimum scale must leave a nonzero active polygon after hinge setbacks.
minimumScale = max([
    2.01 * probeModel.setback / sin(gamma)
    2.01 * probeModel.setback / (aspectRatio * sin(gamma))
    1e-6
]);
if ~fitsAtScale(minimumScale)
    scale = nan;
    return
end

low = minimumScale;
high = max(2 * low, max(options.box));

% Expand high until it is outside the bounding box.
while high < 100 * max(options.box) && fitsAtScale(high)
    low = high;
    high = 2 * high;
end
if high >= 100 * max(options.box) && fitsAtScale(high)
    scale = nan;
    return
end

% Bisection stops either at the hard iteration limit or when the interval is
% smaller than the relative scale tolerance.
for iteration = 1:options.bisectionIterations
    middle = (low + high) / 2;
    if fitsAtScale(middle)
        low = middle;
    else
        high = middle;
    end
    if (high - low) <= options.scaleTolerance * max(1, abs(low))
        break
    end
end
scale = low;

    function fits = fitsAtScale(value)
        % Scale feasibility only needs active area and bounding box, not SAT or
        % residual verification.
        parameters = makeParameters(rows, cols, value * aspectRatio, value, gamma, normalizedEdgeTilt, options, false, false);
        model = buildMiuraModel(parameters);
        fits = all(model.bboxSize <= options.box + 1e-9) && all([model.facets.activeArea] > 0);
    end
end

function isValid = verifyPath(candidate, options)
%VERIFYPATH Run final SAT collision checks along the folding path.
% This is deliberately outside the repeated objective calls for speed.
isValid = true;
if ~candidate.feasible
    isValid = false;
    return
end
for normalizedEdgeTilt = linspace(0, candidate.normalizedEdgeTilt, options.pathSamples)
    parameters = makeParameters(candidate.rows, candidate.cols, candidate.a, candidate.b, ...
        deg2rad(candidate.gammaDegrees), normalizedEdgeTilt, options, true, false);
    model = buildMiuraModel(parameters);
    if any([model.facets.activeArea] <= 0) || model.collisionCount > 0
        isValid = false;
        return
    end
end
% Square-cut setbacks prevent adjacent hinge interference analytically.
% Non-adjacent SAT checks are intentionally reserved for this final path gate.
end

function parameters = makeParameters(rows, cols, a, b, gamma, normalizedEdgeTilt, options, checkCollisions, checkResiduals)
%MAKEPARAMETERS Pack scalar options into the struct expected by buildMiuraModel.
parameters.rows = rows;
parameters.cols = cols;
parameters.a = a;
parameters.b = b;
parameters.gamma = gamma;
parameters.normalizedEdgeTilt = normalizedEdgeTilt;
parameters.thickness = options.thickness;
parameters.clearance = options.clearance;
parameters.checkCollisions = checkCollisions;
parameters.checkResiduals = checkResiduals;
end

function model = makeModel(candidate, options)
%MAKEMODEL Rebuild a candidate with all expensive final checks enabled.
checkCollisions = ~options.skipFinalCollisionChecks;
checkResiduals = ~options.skipFinalCollisionChecks;
model = buildMiuraModel(makeParameters(candidate.rows, candidate.cols, candidate.a, candidate.b, ...
    deg2rad(candidate.gammaDegrees), candidate.normalizedEdgeTilt, options, checkCollisions, checkResiduals));
end

function candidate = blankCandidate(rows, cols)
%BLANKCANDIDATE Return a failed candidate with all expected fields present.
candidate = struct("feasible", false, "rows", rows, "cols", cols, "gammaDegrees", nan, ...
    "normalizedEdgeTilt", nan, "aspectRatio", nan, "scale", nan, "a", nan, "b", nan, ...
    "setback", nan, "activeArea", -inf, "idealArea", nan, "activeAreaRatio", nan, ...
    "areaDensity", -inf, "bboxSize", [nan, nan, nan], "minimumOpenAngleDegrees", nan, ...
    "hingeDihedralRangeDegrees", [nan, nan], "rigidFoldResidual", nan, ...
    "hingeRotationResidual", nan, "collisionCount", -1, "pathCollisionFree", false);
end

function constrained = logistic(unconstrained, lower, upper)
%LOGISTIC Map unconstrained variables into finite lower/upper bounds.
constrained = lower + (upper - lower) ./ (1 + exp(-unconstrained));
end

function unconstrained = inverseLogistic(constrained, lower, upper)
%INVERSELOGISTIC Map bounded variables back to unconstrained space.
ratio = (constrained - lower) ./ (upper - constrained);
unconstrained = log(ratio);
end

function histories = emptyHistory()
%EMPTYHISTORY Return the scalar struct used for one optimizer start.
histories = struct("rows", double.empty(0, 1), "cols", double.empty(0, 1), ...
    "startIndex", double.empty(0, 1), "algorithm", strings(0, 1), ...
    "iteration", [], "functionCount", [], "score", [], "areaDensity", [], ...
    "finalGammaDegrees", nan, "finalNormalizedEdgeTilt", nan, ...
    "finalAspectRatio", nan, "finalAreaDensity", nan, "finalScore", nan);
end

function histories = appendHistories(histories, newHistories)
%APPENDHISTORIES Append scalar or vector history structs while preserving empty initialization.
if isempty(newHistories) || isempty([newHistories.rows])
    return
end
if isempty(histories) || isempty([histories.rows])
    histories = newHistories;
else
    histories = [histories, newHistories]; %#ok<AGROW>
end
end

function historyTable = historiesToTable(histories)
%HISTORIESTOTABLE Flatten nested optimizer histories into one CSV-friendly table.
historyTable = table();
if isempty(histories) || isempty([histories.rows])
    return
end
for historyIndex = 1:numel(histories)
    history = histories(historyIndex);
    if isempty(history.iteration)
        continue
    end
    count = numel(history.iteration);
    partial = table( ...
        repmat(history.rows, count, 1), ...
        repmat(history.cols, count, 1), ...
        repmat(history.startIndex, count, 1), ...
        repmat(string(history.algorithm), count, 1), ...
        history.iteration(:), ...
        history.functionCount(:), ...
        history.score(:), ...
        history.areaDensity(:), ...
        'VariableNames', {'rows', 'cols', 'startIndex', 'algorithm', ...
        'iteration', 'functionCount', 'score', 'areaDensity'});
    historyTable = [historyTable; partial]; %#ok<AGROW>
end
end

function writeSummary(best, options, outputPath)
%WRITESUMMARY Save a compact text summary for the best candidate.
file = fopen(outputPath, "w");
cleanup = onCleanup(@() fclose(file));
fprintf(file, "layout: %d x %d\n", best.rows, best.cols);
fprintf(file, "launcher box: %.8f %.8f %.8f\n", options.box);
fprintf(file, "thickness: %.8f\n", options.thickness);
fprintf(file, "gamma degrees: %.10f\n", best.gammaDegrees);
fprintf(file, "normalized edge tilt: %.10f\n", best.normalizedEdgeTilt);
fprintf(file, "aspect ratio a:b: %.10f\n", best.aspectRatio);
fprintf(file, "a b: %.10f %.10f\n", best.a, best.b);
fprintf(file, "setback q: %.10f\n", best.setback);
fprintf(file, "active area: %.10f\n", best.activeArea);
fprintf(file, "active area ratio: %.10f\n", best.activeAreaRatio);
fprintf(file, "active area / launcher volume: %.10f\n", best.areaDensity);
fprintf(file, "bbox: %.10f %.10f %.10f\n", best.bboxSize);
fprintf(file, "hinge dihedral range degrees: %.10f %.10f\n", best.hingeDihedralRangeDegrees);
fprintf(file, "rigid-fold residual: %.3e\n", best.rigidFoldResidual);
fprintf(file, "Rodrigues residual: %.3e\n", best.hingeRotationResidual);
fprintf(file, "non-adjacent collisions: %d\n", best.collisionCount);
fprintf(file, "optimizer: %s\n", string(options.optimizer));
end

function plotRankings(layouts, outputPath)
%PLOTRANKINGS Save a simple scatter plot of layout size versus objective.
figureHandle = figure("Color", "white", "Name", "Miura layout optimization", "Visible", "off");
facetCounts = [layouts.rows] .* [layouts.cols];
scatter(facetCounts, [layouts.areaDensity], 45, [0.33, 0.62, 0.82], "filled");
grid on;
xlabel("Number of facets");
ylabel("Active area / launcher volume");
title("Miura layout optimization");
exportgraphics(figureHandle, outputPath, "Resolution", 180);
end

function plotConvergence(histories, best, outputPath)
%PLOTCONVERGENCE Save objective-density traces for all starts of the best layout.
if isempty(histories) || isempty([histories.rows])
    return
end
isBestLayout = [histories.rows] == best.rows & [histories.cols] == best.cols;
selected = histories(isBestLayout);
if isempty(selected)
    return
end

figureHandle = figure("Color", "white", "Name", "Objective convergence", "Visible", "off");
hold on;
for index = 1:numel(selected)
    history = selected(index);
    if isempty(history.areaDensity)
        continue
    end
    plot(history.functionCount, history.areaDensity, "LineWidth", 1.2, ...
        "DisplayName", sprintf("start %d", history.startIndex));
end
yline(best.areaDensity, "k--", "best verified", "LineWidth", 1.0, ...
    "LabelHorizontalAlignment", "left");
grid on;
xlabel("Function evaluations");
ylabel("Active area / launcher volume");
title(sprintf("Convergence history for %d x %d layout", best.rows, best.cols));
legend("Location", "best");
exportgraphics(figureHandle, outputPath, "Resolution", 180);
close(figureHandle);
end
