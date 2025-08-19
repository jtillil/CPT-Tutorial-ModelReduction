%% Setup
clear; clc;
addpath(genpath("../../CPT-Tutorial-ModelReduction"))

config.ir       = true;
config.non_ir   = true;
config.analyze  = true;

%% parallel pathways
name = 'SimpleParallelPathways';
namesimple = 'SPP';
% scenario = 'no_crosstalk';
% 
% setupModel(name, namesimple, scenario, config)
% 
% scenario = 'with_crosstalk';
% 
% setupModel(name, namesimple, scenario, config)

scenario = 'C_B_neglectable';

setupModel(name, namesimple, scenario, config)

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

%% add param struct to 

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")

model.I = config2I(model.I, repmat("dyn", [1, model.I.nstates]), []);
model = model2minimal(model);

save("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat", "model")

%% ir indices

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")

[model.ir, model.contr, model.obs, model.t_ind]  = compute_ir_indices_matlabfun(model);

model.threshold = 0.1;
model.analysis = analyse_all_indices(model);

f = figure;
plot(model.t_ind, model.ir.index(:, model.analysis.ir.I_sorted_max_nindex_above_threshold), 'LineWidth',1.5)
xlim([0, 0.15])
legend(model.I.nmstatelegend{model.analysis.ir.I_sorted_max_nindex_above_threshold})

f = figure;
I = model.I;
plot(model.t_ind, model.ir.nindex(:, [I.AVenom, I.Fg, I.CVenom, I.V, I.II, I.Va, I.P, I.Xa_Va, I.IIa_Tmod, I.APC_PS, I.IIa]), 'LineWidth',1.5)
xlim([0, 0.15])

% save("../Core/modelfiles/modelBC_SV40_from_JKn_2024_irmatlabfun.mat", "model")

%% ir indices

load("../Core/modelfiles/modelBC_SV1_from_JKn_2024.mat")

[model.ir, model.contr, model.obs, model.t_ind]  = compute_ir_indices_matlabfun(model);

save("../Core/modelfiles/modelBC_SV1_from_JKn_2024_irmatlabfun.mat", "model")

%% non_ir indices

% load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")

model.env  = compute_non_ir_indices_matlabfun(model,'env');
model.pneg = compute_non_ir_indices_matlabfun(model,'pneg');
model.cneg = compute_non_ir_indices_matlabfun(model,'cneg');
model.pss  = compute_non_ir_indices_matlabfun(model,'pss');

save("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat", "model")

%% 1h model

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")

config = repmat("dyn", [1 model.I.nstates]);
model.t_ref = [0 1];
[model.t_ref, model.X_ref] = simModel_simple(model, config);

save("../Core/modelfiles/modelBC_SV1_from_JKn_2024.mat", "model")

%% make variability

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
load("modelBC_population_from_UFa_2023.mat")

variability.X_ref_pop =  cell([1 1001]);
variability.X0_pop(2:1001, :) = variability.X0_pop;
variability.X0_pop(1, :) = model.X0';
variability.par_pop(2:1001, :) = variability.par_pop;
variability.par_pop(1, :) = model.par';
X0_pop(:, model.I.AVenom) = 0.0075;

tic
for i = 1:1001
    config = repmat("dyn", [1, model.I.nstates]);
    model.X0 = variability.X0_pop(i, :)';
    model.par = variability.par_pop(i, :)';
    [~, Xout] = simModel_simple(model, config);
    variability.X_ref_pop{i} = Xout;
end
toc

save("../Core/modelfiles/modelBC_population_from_UFa_2023.mat", "variability")