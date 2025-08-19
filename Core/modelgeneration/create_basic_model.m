%% paths

addpath(genpath('modelspecification'))
addpath(genpath('modelfiles'))
addpath(genpath('Indices'))

%% model setup

% model.name = 'Hornberg2005EGFRsignalling';
% model.name = 'Wajim12009BloodCoagulation';
model.name = 'Gulati2014BClumped';
% model.name = 'SimpleParallelPathways';
% model.name = 'MMEnzymeKinetics';

% model.scenario = 'Cpss_EpCenv';
% model.scenario = 'with_crosstalk';
% model.scenario = 'in_vivo';

model.name_short = 'Gulati2014BClumped';
model.scenario_short = model.scenario;

% model = feval([model.name '_model_set_up_details'],model);
model = set_up_the_model(model);

%% generate all model components

model = compute_and_analyse_indices(model, 'compute');
fprintf('\nFinished indices. Compute irenv states.')
model = model2irenv(model);
fprintf('\nFinished irenv states. Compute adjacency matrix.')
% model = model2adjacency(model);
fprintf('\nFinished adjacency matrix.')

%% save model

save(['modelfiles/model' model.name_short '_' model.scenario_short '_full.mat'], 'model')
