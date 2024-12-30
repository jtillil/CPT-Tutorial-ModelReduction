%% Setup
addpath(genpath("../../CPT-Tutorial-ModelReduction"))

config.ir       = true;
config.non_ir   = true;
config.analyze  = true;

%% parallel pathways
name = 'SimpleParallelPathways';
namesimple = 'SPP';
scenario = 'no_crosstalk';

setupModel(name, namesimple, scenario, config)

scenario = 'with_crosstalk';

setupModel(name, namesimple, scenario, config)

%% enzyme kinetics
name = 'MMEnzymeKinetics';
namesimple = 'MMEK';
scenario = 'Cpss_Eenv';

setupModel(name, namesimple, scenario, config)

scenario = 'Cpss_EpCenv';

setupModel(name, namesimple, scenario, config)

% %% blood coagulation
% name = 'SimpleParallelPathways';
% namesimple = 'SPP';
% scenario = 'no_crosstalk';
% 
% setupModel(name, namesimple, scenario, config)
% 
% scenario = 'with_crosstalk';
% 
% setupModel(name, namesimple, scenario, config)