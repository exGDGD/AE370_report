function model = buildMiuraModel(parameters)
%BUILDMIURAMODEL Build inset active panels and finite-thickness prisms.
% The output struct contains geometry, active area, bounding box, optional
% residual checks, and optional non-adjacent collision count.

% Optional fields:
% clearance: extra hinge setback added to the theoretical minimum.
% checkCollisions: true enables non-adjacent SAT collision checks.
% checkResiduals: true enables rigid-fold and Rodrigues residual checks.
if ~isfield(parameters, "clearance")
    parameters.clearance = 0;
end
if ~isfield(parameters, "checkCollisions")
    parameters.checkCollisions = true;
end
if ~isfield(parameters, "checkResiduals")
    parameters.checkResiduals = true;
end

% Validate the required geometric parameters.
mustBeInteger(parameters.rows); mustBePositive(parameters.rows);
mustBeInteger(parameters.cols); mustBePositive(parameters.cols);
mustBePositive(parameters.a); mustBePositive(parameters.b);
mustBePositive(parameters.gamma); mustBePositive(parameters.thickness);
mustBeNonnegative(parameters.clearance);

% Generate the rigid Miura vertex grid.
% normalizedEdgeTilt: measure of how much the facets are fold
mesh = rigidMiuraVertices(parameters.rows, parameters.cols, parameters.a, ...
    parameters.b, parameters.gamma, parameters.normalizedEdgeTilt);

% First build zero-setback facets to measure adjacent opening angles.
baseFacets = createFacets(mesh, parameters.rows, parameters.cols, 0, parameters.thickness);
openAngles = adjacentOpenAngles(baseFacets, parameters.rows, parameters.cols);
minimumOpenAngle = min(openAngles);
setback = parameters.thickness / (2 * tan(minimumOpenAngle / 2)) + parameters.clearance;

% Rebuild facets with the required hinge setback applied.
facets = createFacets(mesh, parameters.rows, parameters.cols, setback, parameters.thickness);

% Bounding box is computed from all finite-thickness prism vertices.
allPrismVertices = vertcat(facets.prism);
bboxMin = min(allPrismVertices, [], 1);
bboxMax = max(allPrismVertices, [], 1);
dihedralAngles = pi - openAngles;

model.parameters = parameters;
model.mesh = mesh;
model.facets = facets;
model.setback = setback;
model.minimumOpenAngle = minimumOpenAngle;
model.hingeDihedralRange = [min(dihedralAngles), max(dihedralAngles)];
model.activeArea = sum([facets.activeArea]);
model.idealArea = parameters.rows * parameters.cols * parameters.a * parameters.b * sin(parameters.gamma);
model.activeAreaRatio = model.activeArea / model.idealArea;
model.bboxMin = bboxMin;
model.bboxMax = bboxMax;
model.bboxSize = bboxMax - bboxMin;
if parameters.checkResiduals
    model.rigidFoldResidual = rigidFoldResidual(mesh, parameters.a, parameters.b, parameters.gamma);
    model.hingeRotationResidual = hingeRotationResidual(baseFacets, parameters.rows, parameters.cols);
else
    model.rigidFoldResidual = nan;
    model.hingeRotationResidual = nan;
end
if parameters.checkCollisions
    model.collisionCount = nonAdjacentCollisionCount(facets);
else
    model.collisionCount = nan;
end
end

function facets = createFacets(mesh, rows, cols, setback, thickness)
%CREATEFACETS Convert a vertex grid into facet structs.
% Each facet stores original vertices, inset active vertices, active area,
% normal vector, and the finite-thickness prism.
facets(rows * cols) = struct("index", [], "vertices", [], "normal", [], ...
    "activeVertices", [], "activeArea", 0, "prism", []);
counter = 0;
for row = 1:rows
    for column = 1:cols
        counter = counter + 1;
        vertices = [
            squeeze(mesh(row, column, :)).'
            squeeze(mesh(row, column + 1, :)).'
            squeeze(mesh(row + 1, column + 1, :)).'
            squeeze(mesh(row + 1, column, :)).'
        ];
        edgeSetbacks = [
            setback * (row > 1)
            setback * (column < cols)
            setback * (row < rows)
            setback * (column > 1)
        ];
        activeVertices = insetConvexPolygon(vertices, edgeSetbacks);
        normal = polygonNormal(vertices);
        facets(counter).index = [row, column];
        facets(counter).vertices = vertices;
        facets(counter).normal = normal;
        facets(counter).activeVertices = activeVertices;
        facets(counter).activeArea = polygonArea3d(activeVertices);
        facets(counter).prism = makePrism(activeVertices, thickness);
    end
end
end

function inset = insetConvexPolygon(vertices, setbacks)
%INSETCONVEXPOLYGON Offset each polygon edge inward by its setback.
% The 3D polygon is projected to a local 2D basis, clipped by edge
% half-planes, and then mapped back to 3D.
origin = vertices(1, :);
basisX = unitVector(vertices(2, :) - origin);
normal = polygonNormal(vertices);
basisY = unitVector(cross(normal, basisX));
points = [(vertices - origin) * basisX.', (vertices - origin) * basisY.'];
if signedArea2d(points) < 0
    points = flipud(points);
    setbacks = flipud(setbacks);
end

clipped = points;
for edgeIndex = 1:size(points, 1)
    following = mod(edgeIndex, size(points, 1)) + 1;
    edge = points(following, :) - points(edgeIndex, :);
    inward = unitVector([-edge(2), edge(1)]);
    threshold = dot(inward, points(edgeIndex, :)) + setbacks(edgeIndex);
    clipped = clipHalfPlane(clipped, inward, threshold);
    if size(clipped, 1) < 3
        inset = zeros(0, 3);
        return
    end
end
inset = origin + clipped(:, 1) * basisX + clipped(:, 2) * basisY;
end

function clipped = clipHalfPlane(points, normal, threshold)
%CLIPHALFPLANE Clip a convex 2D polygon to dot(normal,x) >= threshold.
clipped = zeros(0, 2);
for index = 1:size(points, 1)
    following = mod(index, size(points, 1)) + 1;
    currentPoint = points(index, :);
    nextPoint = points(following, :);
    currentValue = dot(normal, currentPoint) - threshold;
    nextValue = dot(normal, nextPoint) - threshold;
    currentInside = currentValue >= -1e-10;
    nextInside = nextValue >= -1e-10;
    if currentInside
        clipped(end + 1, :) = currentPoint; %#ok<AGROW>
    end
    if currentInside ~= nextInside
        weight = currentValue / (currentValue - nextValue);
        clipped(end + 1, :) = currentPoint + weight * (nextPoint - currentPoint); %#ok<AGROW>
    end
end
end

function prism = makePrism(vertices, thickness)
%MAKEPRISM Thicken a polygon by extruding half-thickness along its normal.
if size(vertices, 1) < 3
    prism = zeros(0, 3);
    return
end
offset = polygonNormal(vertices) * thickness / 2;
prism = [vertices + offset; vertices - offset];
end

function angles = adjacentOpenAngles(facets, rows, cols)
%ADJACENTOPENANGLES Compute material opening angles for shared-hinge pairs.
angles = zeros((rows - 1) * cols + rows * (cols - 1), 1);
counter = 0;
for row = 1:rows
    for column = 1:cols
        first = facets(facetNumber(row, column, cols));
        if row < rows
            counter = counter + 1;
            angles(counter) = openAngle(first, facets(facetNumber(row + 1, column, cols)));
        end
        if column < cols
            counter = counter + 1;
            angles(counter) = openAngle(first, facets(facetNumber(row, column + 1, cols)));
        end
    end
end
end

function angle = openAngle(first, second)
%OPENANGLE Return pi for coplanar panels and smaller values as they close.
normalAngle = acos(max(-1, min(1, dot(first.normal, second.normal))));
angle = pi - normalAngle;
end

function residual = rigidFoldResidual(mesh, a, b, gamma)
%RIGIDFOLDRESIDUAL Check edge lengths, sector angle, and cell closure.
residual = 0;
targetDot = abs(a * b * cos(gamma));
for row = 1:(size(mesh, 1) - 1)
    for column = 1:(size(mesh, 2) - 1)
        first = squeeze(mesh(row, column + 1, :) - mesh(row, column, :)).';
        second = squeeze(mesh(row + 1, column, :) - mesh(row, column, :)).';
        closure = squeeze(mesh(row, column, :) + mesh(row, column + 1, :) - ...
            mesh(row, column, :) + mesh(row + 1, column, :) - mesh(row, column, :) - ...
            mesh(row + 1, column + 1, :)).';
        residual = max([residual, abs(norm(first) - a), abs(norm(second) - b), ...
            abs(abs(dot(first, second)) - targetDot), norm(closure)]);
    end
end
end

function residual = hingeRotationResidual(facets, rows, cols)
%HINGEROTATIONRESIDUAL Check that adjacent normals match by hinge rotation.
residual = 0;
for row = 1:rows
    for column = 1:cols
        first = facets(facetNumber(row, column, cols));
        if row < rows
            residual = max(residual, pairResidual(first, facets(facetNumber(row + 1, column, cols))));
        end
        if column < cols
            residual = max(residual, pairResidual(first, facets(facetNumber(row, column + 1, cols))));
        end
    end
end
end

function residual = pairResidual(first, second)
%PAIRRESIDUAL Rodrigues-rotate one facet normal about the shared hinge axis.
shared = zeros(0, 3);
for firstIndex = 1:4
    for secondIndex = 1:4
        if norm(first.vertices(firstIndex, :) - second.vertices(secondIndex, :)) < 1e-10
            shared(end + 1, :) = first.vertices(firstIndex, :); %#ok<AGROW>
        end
    end
end
if size(shared, 1) ~= 2
    residual = inf;
    return
end
axis = unitVector(shared(2, :) - shared(1, :));
angle = acos(max(-1, min(1, dot(first.normal, second.normal))));
positive = rodriguesRotate(first.normal, axis, angle);
negative = rodriguesRotate(first.normal, axis, -angle);
residual = min(norm(positive - second.normal), norm(negative - second.normal));
end

function count = nonAdjacentCollisionCount(facets)
%NONADJACENTCOLLISIONCOUNT Count prism intersections for non-neighbor facets.
count = 0;
for firstIndex = 1:numel(facets)
    first = facets(firstIndex);
    for secondIndex = (firstIndex + 1):numel(facets)
        second = facets(secondIndex);
        if sum(abs(first.index - second.index)) <= 1
            continue
        end
        firstMin = min(first.prism, [], 1); firstMax = max(first.prism, [], 1);
        secondMin = min(second.prism, [], 1); secondMax = max(second.prism, [], 1);
        if any(firstMax <= secondMin + 1e-10) || any(secondMax <= firstMin + 1e-10)
            continue
        end
        if prismsOverlap(first.prism, second.prism)
            count = count + 1;
        end
    end
end
end

function overlaps = prismsOverlap(first, second)
%PRISMSOVERLAP Separating Axis Theorem test for two convex prisms.
% Mere point, edge, or face contact is allowed; positive overlap on all
% tested axes is counted as intersection.
firstDirections = prismDirections(first);
secondDirections = prismDirections(second);
axes = zeros(0, 3);
directions = [firstDirections; secondDirections];
for firstIndex = 1:size(directions, 1)
    for secondIndex = (firstIndex + 1):size(directions, 1)
        axis = cross(directions(firstIndex, :), directions(secondIndex, :));
        if norm(axis) > 1e-10
            axes(end + 1, :) = unitVector(axis); %#ok<AGROW>
        end
    end
end
overlaps = true;
for axisIndex = 1:size(axes, 1)
    axis = axes(axisIndex, :).';
    firstProjection = first * axis;
    secondProjection = second * axis;
    overlap = min(max(firstProjection), max(secondProjection)) - ...
        max(min(firstProjection), min(secondProjection));
    if overlap <= 1e-10
        overlaps = false;
        return
    end
end
end

function directions = prismDirections(prism)
%PRISMDIRECTIONS Return edge and thickness directions used to form SAT axes.
topCount = size(prism, 1) / 2;
directions = zeros(0, 3);
for index = 1:topCount
    following = mod(index, topCount) + 1;
    edge = prism(following, :) - prism(index, :);
    if norm(edge) > 1e-10
        directions(end + 1, :) = unitVector(edge); %#ok<AGROW>
    end
end
thicknessDirection = prism(topCount + 1, :) - prism(1, :);
if norm(thicknessDirection) > 1e-10
    directions(end + 1, :) = unitVector(thicknessDirection);
end
end

function rotated = rodriguesRotate(vector, axis, angle)
%RODRIGUESROTATE Rotate a vector about an axis by Rodrigues' formula.
axis = unitVector(axis);
rotated = vector * cos(angle) + cross(axis, vector) * sin(angle) + ...
    axis * dot(axis, vector) * (1 - cos(angle));
end

function area = polygonArea3d(vertices)
%POLYGONAREA3D Compute area of a planar 3D polygon.
if size(vertices, 1) < 3
    area = 0;
    return
end
normal = polygonNormal(vertices);
area = 0;
for index = 1:size(vertices, 1)
    following = mod(index, size(vertices, 1)) + 1;
    area = area + dot(cross(vertices(index, :), vertices(following, :)), normal);
end
area = abs(area) / 2;
end

function normal = polygonNormal(vertices)
%POLYGONNORMAL Unit normal of a planar polygon.
normal = unitVector(cross(vertices(2, :) - vertices(1, :), vertices(end, :) - vertices(1, :)));
end

function area = signedArea2d(points)
%SIGNEDAREA2D Signed doubled area for orientation of a 2D polygon.
following = [2:size(points, 1), 1];
area = sum(points(:, 1) .* points(following, 2) - points(following, 1) .* points(:, 2));
end

function number = facetNumber(row, column, cols)
%FACETNUMBER Convert row and column indices to a linear facet index.
number = (row - 1) * cols + column;
end

function normalized = unitVector(vector)
%UNITVECTOR Normalize a vector and fail loudly for degenerate input.
lengthValue = norm(vector);
if lengthValue < 1e-12
    error("buildMiuraModel:ZeroVector", "Cannot normalize a zero-length vector.");
end
normalized = vector / lengthValue;
end
