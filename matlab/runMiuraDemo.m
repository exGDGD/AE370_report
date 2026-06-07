function model = runMiuraDemo()
%RUNMIURADEMO Generate and visualize a default finite-thickness Miura model.

parameters.rows = 4;
parameters.cols = 6;
parameters.a = 1.0;
parameters.b = 1.0;
parameters.gamma = deg2rad(60);
parameters.normalizedEdgeTilt = 0.80;
parameters.thickness = 0.06;
parameters.clearance = 0;

model = buildMiuraModel(parameters);
visualizeMiura(model, "miura_demo.png");

fprintf("facets:                  %d\n", parameters.rows * parameters.cols);
fprintf("setback q:               %.8f\n", model.setback);
fprintf("active area:             %.8f\n", model.activeArea);
fprintf("active area ratio:       %.4f%%\n", 100 * model.activeAreaRatio);
fprintf("bbox:                    %.5f x %.5f x %.5f\n", model.bboxSize);
fprintf("hinge dihedral range:    %.4f to %.4f deg\n", rad2deg(model.hingeDihedralRange));
fprintf("rigid-fold residual:     %.3e\n", model.rigidFoldResidual);
fprintf("Rodrigues residual:      %.3e\n", model.hingeRotationResidual);
fprintf("non-adjacent collisions: %d\n", model.collisionCount);
end
