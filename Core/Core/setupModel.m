function setupModel(name, namesimple, scenario, config)

%% setup
addpath(genpath('../Core'))

%% model 
% model.name = 'Hornberg2005EGFRsignalling';
% model.scenario = '-';

% model.name = 'SimpleParallelPathways';
% model.scenario = 'with_crosstalk';

% model.name = 'MMEnzymeKinetics';
% model.scenario = 'nothing';

% model.scenario = 'basic_irek';
% model.scenario = 'basic_irek_enzyme_excess';

model.name = name;
model.scenario = scenario;
model.savenameroot = ['./results/' namesimple];
if exist([model.name '_model_set_up_details'], "file") == 2
    model = feval([model.name '_model_set_up_details'],model);
end
model.ode_is_matlabfun = 0;
model = set_up_the_model(model);

%% calculate indices

if config.ir
    [ir, contr, obs] = compute_ir_indices(model, false);
    model.ir = ir;
    model.contr = contr;
    model.obs = obs;
end

if config.non_ir
    model.env = compute_non_ir_indices(model, "env", false);
    model.pneg = compute_non_ir_indices(model, "pneg", false);
    model.cneg = compute_non_ir_indices(model, "cneg", false);
    model.pss = compute_non_ir_indices(model, "pss", false);
end

if config.analyze
    % model = compute_and_analyse_indices(model, 'compute');
    model.threshold = 0.1;
    model.analysis = analyse_all_indices(model);
end

%% add relevant components
model = model2minimal(model);
model.ode_is_matlabfun = 1;

%% save model

save(['../Core/modelfiles/model' namesimple '_' scenario '_full.mat'], 'model');

end
