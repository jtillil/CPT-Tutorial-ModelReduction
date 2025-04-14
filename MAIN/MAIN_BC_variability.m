clear; clc;
addpath(genpath("../../CPT-Tutorial-ModelReduction"))

% load model
load("modelBC_SV40_from_JKn_2024.mat")
config = repmat("dyn", [1, model.I.nstates]);
model.I = config2I(model.I, config, []);
par = model.par;
X0 = model.X0;

% generate virtual pop
seed = 1234;
rng(seed, "twister");
Npop = 100;

virtual_pop_par = zeros(Npop, length(par));
% lognpdf_par = zeros(Npop, length(par));
for i = 1:length(par)
    if par(i) ~= 0
        m = par(i);
        v = (0.2*par(i))^2;
        mu = log(m^2 / sqrt(v + m^2));
        sigma = sqrt(log(v / (m^2) + 1));

        virtual_pop_par(:, i) = lognrnd(mu, sigma, [Npop 1]);
        % lognpdf_par(:, i) = lognpdf(virtual_pop_par(:, i), mu, sigma);
    else
        % lognpdf_par(:, i) = 1;
    end
end

% virtual_pop_X0 = zeros(Npop, length(X0));
% % lognpdf_X0 = zeros(Npop, length(X0));
% for i = 1:length(X0)
%     if X0(i) ~= 0
%         m = X0(i);
%         v = (0.4*X0(i))^2;
%         mu = log(m^2 / sqrt(v + m^2));
%         sigma = sqrt(log(v / (m^2) + 1));
% 
%         virtual_pop_X0(:, i) = lognrnd(mu, sigma, [Npop 1]);
%         % lognpdf_X0(:, i) = lognpdf(virtual_pop_X0(:, i), mu, sigma);
%     else
%         % lognpdf_X0(:, i) = 1;
%     end
% end
% virtual_pop_X0(:, model.I.input) = X0(model.I.input);
% % lognpdf_X0(:, model.I.input) = 1;

% steady-state initials and reference solution for virtual pop
X0_no_input = X0;
X0_no_input(model.I.input) = 0;
virtual_pop_par_initials_ss = zeros(Npop, length(X0));
X_ref_var_ss = cell([1 Npop]);
for npop = 1:Npop
    [~, curr_X_ref, ~, ~] = simModel([0 2000], X0_no_input, virtual_pop_par(npop, :)', model.I, model.param, model.multiple, model.odefun, model.jacfun);
    virtual_pop_par_initials_ss(npop, :) = curr_X_ref(end, :);
    virtual_pop_par_initials_ss(npop, model.I.input) = X0(model.I.input);
    virtual_pop_par_initials_ss(npop, (virtual_pop_par_initials_ss(npop, :) < 1e-10)) = 0;
    [~, X_ref_var_ss{npop}, ~, ~] = simModel(model.t_ref, virtual_pop_par_initials_ss(npop, :)', virtual_pop_par(npop, :)', model.I, model.param, model.multiple, model.odefun, model.jacfun);
end

% virtual_pop_X0(:, model.I.input) = 0;
% virtual_pop_X0_ss = zeros(Npop, length(X0));
% X_ref_var_ss = cell([1 Npop]);
% for npop = 1:Npop
%     [~, curr_X_ref, ~, ~] = simModel([0 2000], virtual_pop_X0(npop, :)', model.par, model.I, model.param, model.multiple, model.odefun, model.jacfun);
%     virtual_pop_X0_ss(npop, :) = curr_X_ref(end, :);
%     virtual_pop_X0_ss(npop, model.I.input) = X0(model.I.input);
%     virtual_pop_X0_ss(npop, (virtual_pop_X0_ss(npop, :) < 1e-10)) = 0;
%     [~, X_ref_var_ss{npop}, ~, ~] = simModel(model.t_ref, virtual_pop_X0_ss(npop, :)', model.par, model.I, model.param, model.multiple, model.odefun, model.jacfun);
% end

% read out distribution of steady-state initials

% choose 5 most extreme steady-state individuals
% lognpdf_par_indv = prod(lognpdf_par, 2);
% [vals, idx] = sort(lognpdf_par_indv);
% extreme_vals_par_5 = vals(1:5);
% extreme_indv_par_5 = idx(1:5);

% lognpdf_X0_indv = prod(lognpdf_X0, 2);
% [vals, idx] = sort(lognpdf_X0_indv);
% extreme_vals_X0_5 = vals(1:5);
% extreme_indv_X0_5 = idx(1:5);

% compute indices for 5 extreme individuals
% indv_nr = 2; % can be 1, 2, 3, 4, 5
% indv_id = extreme_indv_X0_5(indv_nr);
% model.X0 = virtual_pop_X0(indv_nr, :)';
% [ir, contr, obs] = compute_ir_indices_matlabfun(model);
% save(['results/irindices_BCSV40_from_JKn_full_extreme' char(string(indv_nr)) '_seed' char(string(seed)) '.mat'], "ir", "contr", "obs")

% reduce model for ref + 5 extreme indvs
% TODO after indices done

% reduced model solutions for virtual pop after steady-state
load("BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv.mat")
redconfig_8_state = redmodel.redobj.redconfig;
errors = zeros(Npop, 1);
for i = 1:Npop
    [obj, ~, ~, ~, ~] = objfun(model.t_ref, X_ref_var_ss{i}, virtual_pop_par_initials_ss(i, :)', virtual_pop_par(i, :)', model.I, [], model.param, model.multiple, model.odefun, model.jacfun, redconfig_8_state, "MRSE");
    errors(i) = obj.errout;
end
mean(errors)
sum(errors < 0.1)