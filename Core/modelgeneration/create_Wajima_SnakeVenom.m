% load blood coagulation model
% model

% change relevant properties
model.I = config2I(model.I, repmat("dyn", [1, model.I.nstates]), model.L);

model.scenario = 'in_vivo_snake_venom';

model.I.input = model.I.AVenom;
model.I.output = model.I.Fg;

multiple.multiple = 0;
model.multiple = multiple;

dose_snake_venom = 0.0015; % [mg]
SF_mg_to_nmol = 1e-3 / 2e5 * 1e9;
u_ref = SF_mg_to_nmol * dose_snake_venom;
model.X0(model.I.AVenom) = 0 + u_ref;

model.par(model.I.v10) = 25000.0;
model.par(model.I.k10) = 1800.0;

model.par(model.I.v14) = 21000.0;
model.par(model.I.k14) = 30000.0;

model.t_ref = 0.0:0.1:40.0; % 40h
[model.t_ref, model.X_ref] = simModel(model.t_ref, model.X0, model.par, model.I, model.param, model.multiple, model.odefun, model.jacfun);

% save model
save("modelfiles/modelBCSnake_minimal.mat", "model")

% plot
plot(model.t_ref, model.X_ref(:, model.I.output));
