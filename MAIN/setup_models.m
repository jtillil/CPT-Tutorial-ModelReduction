%% Setup
clear; clc;
addpath(genpath("../../CPT-Tutorial-ModelReduction"))

config.ir       = true;
config.non_ir   = true;
config.analyze  = true;

%% parallel pathways
% name = 'SimpleParallelPathways';
% namesimple = 'SPP';
% scenario = 'no_crosstalk';
% 
% setupModel(name, namesimple, scenario, config)
% 
% scenario = 'with_crosstalk';
% 
% setupModel(name, namesimple, scenario, config)

%% enzyme kinetics
% name = 'MMEnzymeKinetics';
% namesimple = 'MMEK';
% scenario = 'Cpss_Eenv';
% 
% setupModel(name, namesimple, scenario, config)
% 
% scenario = 'Cpss_EpCenv';
% 
% setupModel(name, namesimple, scenario, config)

%% blood coagulation
% name = 'Wajima2009BloodCoagulation';
% namesimple = 'BC';
% scenario = 'in_vivo_snakevenom_40h';
% 
% setupModel(name, namesimple, scenario, config)

%% blood coagulation: snake venom 40h (Wajima 2009)
% model = struct;
% model.name = 'Wajima2009BloodCoagulation';
% model.namesimple = 'BC';
% model.scenario = 'in_vivo_snakevenom_40h';
% model = Wajima2009BloodCoagulation_model_set_up_details(model);
% model = set_up_the_model(model);

model = load("Wajima2009BloodCoagulation_in_vivo_snakevenom_40h_ir_index.mat").model;
model = model2minimal(model);
model.ir = load("Wajima2009BloodCoagulation_in_vivo_snakevenom_40h_ir_index.mat").indx;
model.contr = load("Wajima2009BloodCoagulation_in_vivo_snakevenom_40h_contr_index.mat").indx;
model.obs = load("Wajima2009BloodCoagulation_in_vivo_snakevenom_40h_obs_index.mat").indx;
model.pneg = load("Wajima2009BloodCoagulation_in_vivo_snakevenom_40h_pneg_index.mat").indx;
model.cneg = load("Wajima2009BloodCoagulation_in_vivo_snakevenom_40h_cneg_index.mat").indx;
model.env = load("Wajima2009BloodCoagulation_in_vivo_snakevenom_40h_env_index.mat").indx;
model.pss = load("Wajima2009BloodCoagulation_in_vivo_snakevenom_40h_pss_index.mat").indx;
model.threshold = 0.1;
model.analysis = analyse_all_indices(model);

save("../Core/modelfiles/modelBC_SV40_from_JKn_2024_t0_1.mat", "model")

model.threshold = 0.05;
model.analysis = analyse_all_indices(model);

save("../Core/modelfiles/modelBC_SV40_from_JKn_2024_t0_05.mat", "model")

model.threshold = 0.01;
model.analysis = analyse_all_indices(model);

save("../Core/modelfiles/modelBC_SV40_from_JKn_2024_t0_01.mat", "model")

%% blood coagulation
% model = struct;
% model.name = 'Wajima2009BloodCoagulation';
% model.namesimple = 'BC';
% model.scenario = 'SV';
% 
% model.I = Wajima2009BloodCoagulation_indexing();
% model.I = config2I(model.I, repmat("dyn", [1, model.I.nstates]), 0);
% model.X0 = Wajima2009BloodCoagulation_initialvalues(model);
% model.par = Wajima2009BloodCoagulation_parameters(model);
% model.ode = @(t, X) Wajima2009BloodCoagulation_ode(t,X,par,model);
% 
% % add relevant components
% model = create_full_model(model);
% model = model2minimal(model);
% 
% % calculate indices
% if config.ir
%     [ir, contr, obs] = compute_ir_indices(model, false);
%     model.ir = ir;
%     model.contr = contr;
%     model.obs = obs;
% end
% if config.non_ir
%     model.env = compute_non_ir_indices(model, "env", false);
%     model.pneg = compute_non_ir_indices(model, "pneg", false);
%     model.cneg = compute_non_ir_indices(model, "cneg", false);
%     model.pss = compute_non_ir_indices(model, "pss", false);
% end
% if config.analyze
%     % model = compute_and_analyse_indices(model, 'compute');
%     model.threshold = 0.1;
%     model.analysis = analyse_all_indices(model);
% end
% 
% % add relevant components
% model = model2minimal(model);
% 
% % save model
% save(['../Core/modelfiles/model' namesimple '_' scenario '_full.mat'], 'model');

%% blood coagulation: snake venom 40h (Gulati 2014)
% model = struct;
% model.name = 'Gulati2014BloodCoagulation';
% model.namesimple = 'BC_Gulati2014';
% model.scenario = 'in_vivo';
% 
% model.I = Gulati2014BloodCoagulation_indexing();
% model.I = config2I(model.I, repmat("dyn", [1, model.I.nstates]), 0);
% model.X0 = Gulati2014BloodCoagulation_initialvalues(model);
% model.par = Gulati2014BloodCoagulation_parameters(model);
% model.ode = @(t, X) Gulati2014BloodCoagulation_ode(t,X,par,model);
% 
% % add relevant components
% model = create_full_model(model);
% model = model2minimal(model);
% 
% % calculate indices
% if config.ir
%     [ir, contr, obs] = compute_ir_indices(model, false);
%     model.ir = ir;
%     model.contr = contr;
%     model.obs = obs;
% end
% if config.non_ir
%     model.env = compute_non_ir_indices(model, "env", false);
%     model.pneg = compute_non_ir_indices(model, "pneg", false);
%     model.cneg = compute_non_ir_indices(model, "cneg", false);
%     model.pss = compute_non_ir_indices(model, "pss", false);
% end
% if config.analyze
%     % model = compute_and_analyse_indices(model, 'compute');
%     model.threshold = 0.1;
%     model.analysis = analyse_all_indices(model);
% end
% 
% % add relevant components
% model = model2minimal(model);
% 
% % save model
% save(['../Core/modelfiles/model' namesimple '_' scenario '_full.mat'], 'model');
