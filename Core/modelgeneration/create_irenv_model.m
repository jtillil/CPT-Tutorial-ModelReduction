addpath(genpath('modelfiles'))

model = load('modelBC_full.mat').model;
model = rmfield(model, 'obs');

if ~isfield(model, 'obs')
    [model.ir, model.contr, model.obs] = compute_ir_indices(model, false);
end

model = model2irenv(model);

save("modelfiles/modelBC_full.mat", "model")