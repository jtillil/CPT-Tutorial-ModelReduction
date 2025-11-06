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

%% generate virtual population

% setup
seed = 1234;                % RNG seed
rng(seed, "twister");       % RNG
% Npop = 1000;              % number of virtual individuals
Npop = 10;
CV = 0.4;                   % coefficient of variation of the parameters

% load model
load("modelBC_SV40_from_JKn_2024.mat")
par = model.par;
X0 = model.X0;

virtual_pop_par = zeros(Npop+1, length(par));
for i = 1:length(par)
    if par(i) ~= 0
        m = par(i);
        v = (CV*par(i))^2;
        mu = log(m^2 / sqrt(v + m^2));
        sigma = sqrt(log(v / (m^2) + 1));
        % disp(m)
        % disp(sqrt(v))
        % disp(sigma)

        virtual_pop_par(2:end, i) = lognrnd(mu, sigma, [Npop 1]);
        % lognpdf_par(:, i) = lognpdf(virtual_pop_par(:, i), mu, sigma);
    else
        % lognpdf_par(:, i) = 1;
    end
end

% steady-state initials and reference solution for virtual pop
X0_no_input = X0;
X0_no_input(model.I.input) = 0;
virtual_pop_par_initials_ss = zeros(Npop+1, length(X0));
X_ref_var_ss = cell([1 Npop+1]);
for npop = 2:(Npop + 1)
    [~, curr_X_ref, ~, ~] = simModel([0 2000], X0_no_input, virtual_pop_par(npop, :)', model.I, model.param, model.multiple, model.odefun, model.jacfun);
    virtual_pop_par_initials_ss(npop, :) = curr_X_ref(end, :);
    virtual_pop_par_initials_ss(npop, model.I.input) = X0(model.I.input);
    virtual_pop_par_initials_ss(npop, (virtual_pop_par_initials_ss(npop, :) < 1e-10)) = 0;
    [~, X_ref_var_ss{npop}, ~, ~] = simModel(model.t_ref, virtual_pop_par_initials_ss(npop, :)', virtual_pop_par(npop, :)', model.I, model.param, model.multiple, model.odefun, model.jacfun);
end

X_ref_var_ss{1} = model.X_ref;
virtual_pop_par_initials_ss(1, :) = X0';
virtual_pop_par(1, :) = par';

variability.X0_pop = virtual_pop_par_initials_ss;       % steady-state initials for every set of parameters
variability.par_pop = virtual_pop_par;                  % parameter sets of virtual population
variability.X_ref_pop = X_ref_var_ss;                   % ODE system trajectories with virtual population parameters and steady-state initials
variability.variability = 1;                            % wheter to reduce the model with respect to the virtual population

save(['../Core/modelfiles/modelBC_SV40_population_CV' char(num2str(100*CV)) '.mat'], "variability")

%% derive 25-state model - reduce model for virtual population

model_BC_SV1 = load("modelBC_SV1_from_JKn_2024.mat").model;
model = load("modelBC_SV40_from_JKn_2024.mat").model;
load("modelBC_SV40_population_CV40.mat")

redmodel_seq1 = mor_sequential_JKn_2018(model, flip(model_BC_SV1.analysis.ir.I_sorted_max_nindex), 0.1);

redmodel_backwards = mor_backwards_UFa_2023(redmodel_seq1, flip(model.analysis.ir.I_sorted_max_nindex), 0.1, variability, 95);

redmodel = mor_sequential_JKn_2018(redmodel_backwards, flip(model.analysis.ir.I_sorted_max_nindex), 0.1, variability, 95);

save("./results/modelBC_SV40_redvar_CV40_popprct95.mat", "redmodel")