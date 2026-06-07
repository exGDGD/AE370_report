function candidate = evaluateMiuraDesign(rows, cols, gammaDegrees, normalizedEdgeTilt, aspectRatio, options)
%EVALUATEMIURADESIGN Evaluate one fixed Miura design point.
% This helper mirrors the fast objective evaluation used by optimizeMiura:
% scale is chosen by bounding-box bisection, while SAT collision checks and
% residual checks are skipped unless the caller enables them after the fact.

if nargin < 6
    options = struct();
end
options = setDefault(options, "box", [5.2, 3.0, 0.8]);
options = setDefault(options, "thickness", 0.06);
options = setDefault(options, "clearance", 0);
options = setDefault(options, "bisectionIterations", 10);
options = setDefault(options, "scaleTolerance", 1e-4);

gamma = deg2rad(gammaDegrees);
scale = maximizeScale(rows, cols, gamma, normalizedEdgeTilt, aspectRatio, options);

candidate = blankCandidate(rows, cols);
if ~isfinite(scale)
    return
end

parameters = makeParameters(rows, cols, scale * aspectRatio, scale, gamma, normalizedEdgeTilt, options);
model = buildMiuraModel(parameters);

candidate.feasible = all(model.bboxSize <= options.box + 1e-9) && all([model.facets.activeArea] > 0);
candidate.rows = rows;
candidate.cols = cols;
candidate.gammaDegrees = gammaDegrees;
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
candidate.collisionCount = model.collisionCount;

    function options = setDefault(options, name, value)
        if ~isfield(options, name)
            options.(name) = value;
        end
    end

    function scale = maximizeScale(rows, cols, gamma, normalizedEdgeTilt, aspectRatio, options)
        probe = makeParameters(rows, cols, aspectRatio, 1, gamma, normalizedEdgeTilt, options);
        probeModel = buildMiuraModel(probe);
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
        while high < 100 * max(options.box) && fitsAtScale(high)
            low = high;
            high = 2 * high;
        end
        if high >= 100 * max(options.box) && fitsAtScale(high)
            scale = nan;
            return
        end

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
            parameters = makeParameters(rows, cols, value * aspectRatio, value, gamma, normalizedEdgeTilt, options);
            model = buildMiuraModel(parameters);
            fits = all(model.bboxSize <= options.box + 1e-9) && all([model.facets.activeArea] > 0);
        end
    end

    function parameters = makeParameters(rows, cols, a, b, gamma, normalizedEdgeTilt, options)
        parameters.rows = rows;
        parameters.cols = cols;
        parameters.a = a;
        parameters.b = b;
        parameters.gamma = gamma;
        parameters.normalizedEdgeTilt = normalizedEdgeTilt;
        parameters.thickness = options.thickness;
        parameters.clearance = options.clearance;
        parameters.checkCollisions = false;
        parameters.checkResiduals = false;
    end

    function candidate = blankCandidate(rows, cols)
        candidate = struct("feasible", false, "rows", rows, "cols", cols, ...
            "gammaDegrees", nan, "normalizedEdgeTilt", nan, "aspectRatio", nan, ...
            "scale", nan, "a", nan, "b", nan, "setback", nan, ...
            "activeArea", -inf, "idealArea", nan, "activeAreaRatio", nan, ...
            "areaDensity", -inf, "bboxSize", [nan, nan, nan], ...
            "minimumOpenAngleDegrees", nan, "hingeDihedralRangeDegrees", [nan, nan], ...
            "collisionCount", -1);
    end
end
