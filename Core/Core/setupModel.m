function setupModel(name, namesimple, scenario, config)

%% setup
addpath(genpath('../../Core'))

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
model = feval([model.name '_model_set_up_details'],model);
model = set_up_the_model(model);

%% calculate indices

if config.ir
    [ir, contr, obs] = compute_ir_indices(model, false);
    model.ir = ir;
    model.contr = contr;
    model.obs = obs;
end

if config.non_ir
    model = compute_non_ir_indices(model, ["env" "pneg" "cneg" "pss"], false);
end

if config.analyze
    % model = compute_and_analyse_indices(model, 'compute');
    model.threshold = 0.1;
    model = analyse_all_indices(model);
end

%% save model

save(['../modelfiles/model' namesimple '.mat'], 'model');

end
