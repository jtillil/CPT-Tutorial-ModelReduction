%% Setup
clear; clc;
addpath(genpath("../."))

%% derive 8-state model - reduce via ir-indices of 1h scenario for 40h model (using pneg, env)

model_BC_SV1 = load("modelBC_SV1_from_JKn_2024.mat").model;
model = load("modelBC_SV40_from_JKn_2024.mat").model;

redmodel = mor_sequential_JKn_2018(model, flip(model_BC_SV1.analysis.ir.I_sorted_max_nindex), 0.1);

save("./results/modelBC_SV40_from_JKn_2024_reduced_8_state.mat", "redmodel")

%% derive 13-state model - reduce via ir-indices of 40h scenario for 40h model (using pneg, env)

model = load("modelBC_SV40_from_JKn_2024.mat").model;

redmodel = mor_sequential_JKn_2018(model, flip(model.analysis.ir.I_sorted_max_nindex), 0.1);

save("./results/modelBC_SV40_from_JKn_2024_reduced_13_state.mat", "redmodel")

%% ir-indices of 8-state model

% load model
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")

[ir, contr, obs, t_ind] = compute_ir_indices_matlabfun(redmodel);

save("./results/modelBC_SV40_from_JKn_2024_indices_8_state_full.mat", "ir", "contr", "obs", "t_ind")

%% ir-indices of 13-state model

% load model
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")

[ir, contr, obs, t_ind] = compute_ir_indices_matlabfun(redmodel);

save("./results/modelBC_SV40_from_JKn_2024_indices_13_state_full.mat", "ir", "contr", "obs", "t_ind")
