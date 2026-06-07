function vertices = rigidMiuraVertices(rows, cols, a, b, gamma, normalizedEdgeTilt)
%RIGIDMIURAVERTICES Generate a one-DOF rigid-folded Miura mesh.
% gamma is the fixed acute facet angle in radians.
% normalizedEdgeTilt = beta / gamma, where beta is the out-of-plane tilt of
% the a-directed edge. It is not a hinge dihedral angle.

arguments
    rows (1,1) double {mustBeInteger, mustBePositive}
    cols (1,1) double {mustBeInteger, mustBePositive}
    a (1,1) double {mustBePositive}
    b (1,1) double {mustBePositive}
    gamma (1,1) double {mustBePositive}
    normalizedEdgeTilt (1,1) double {mustBeGreaterThanOrEqual(normalizedEdgeTilt, 0), mustBeLessThan(normalizedEdgeTilt, 1)}
end

if gamma >= pi / 2
    error("rigidMiuraVertices:InvalidGamma", "Require 0 < gamma < pi/2.");
end

beta = normalizedEdgeTilt * gamma;
rowX = cos(gamma) / cos(beta);
rowY = sqrt(max(0, 1 - rowX^2));

columnEdges = zeros(cols, 3);
for column = 1:cols
    columnEdges(column, :) = [a * cos(beta), 0, (-1)^(column - 1) * a * sin(beta)];
end

rowEdges = zeros(rows, 3);
for row = 1:rows
    rowEdges(row, :) = [(-1)^(row - 1) * b * rowX, b * rowY, 0];
end

vertices = zeros(rows + 1, cols + 1, 3);
for column = 2:(cols + 1)
    vertices(1, column, :) = vertices(1, column - 1, :) + reshape(columnEdges(column - 1, :), 1, 1, 3);
end
for row = 2:(rows + 1)
    vertices(row, 1, :) = vertices(row - 1, 1, :) + reshape(rowEdges(row - 1, :), 1, 1, 3);
    for column = 2:(cols + 1)
        vertices(row, column, :) = vertices(row, column - 1, :) + reshape(columnEdges(column - 1, :), 1, 1, 3);
    end
end
end
