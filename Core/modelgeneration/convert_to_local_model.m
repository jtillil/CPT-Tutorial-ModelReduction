%% paths

% addpath(genpath('modelspecification'))
% addpath(genpath('modelfiles'))
% addpath(genpath('Indices'))

%% model setup

localmodel.name = 'Hornberg2005EGFRsignalling';
% model.name = 'SimpleParallelPathways';
localmodel.scenario = '-';

localmodel = feval([localmodel.name '_model_set_up_details'],localmodel);
localmodel = set_up_the_model(localmodel);

%% convert additional model components

load("modelEGFR_full.mat")

localmodel.analysis = model.analysis;
localmodel.irenv_arith = model.irenv_arith;
localmodel.irenv_geom = model.irenv_geom;

%% save model

save("../modelfiles/modelEGFR_full_local.mat", 'model')
