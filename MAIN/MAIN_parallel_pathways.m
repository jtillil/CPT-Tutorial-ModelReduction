%% Setup
addpath(genpath("../../../CPT-Tutorial-ModelReduction"))
reduced_errors = struct;

%% Scenario 1 - no crosstalk
% load model file
load("modelSPP_no_crosstalk_full.mat")

% lumping
[lump_matrices,inv_lump_matrices,errors,out_states] = lumping_Aarons(model);

reduced_errors.lumping_no_crosstalk(1,1) = calculate_lumping_error(model, lump_matrices{1});

% S env

% S 

%% Scenario 2 - with crosstalk
% load model file
load("modelSPP_with_crosstalk_full.mat")

%% Analyze both scenarios


%% Save results

save("./results/errors_SPP.mat", "reduced_errors")