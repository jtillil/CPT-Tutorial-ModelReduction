clear; clc;
addpath(genpath("../../CPT-Tutorial-ModelReduction"))

% load model
load("modelBC_SV40_from_JKn_2024.mat")
config = repmat("dyn", [1, model.I.nstates]);
model.I = config2I(model.I, config, []);
model.L = [];
par = model.par;
X0 = model.X0;

% generate virtual pop
seed = 1234;
rng(seed, "twister");
Npop = 100;

virtual_pop_par = zeros(Npop+1, length(par));
% lognpdf_par = zeros(Npop, length(par));
for i = 1:length(par)
    if par(i) ~= 0
        m = par(i);
        v = (0.2*par(i))^2;
        mu = log(m^2 / sqrt(v + m^2));
        sigma = sqrt(log(v / (m^2) + 1));

        virtual_pop_par(2:end, i) = lognrnd(mu, sigma, [Npop 1]);
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
% load("BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv.mat")
% redconfig_8_state = redmodel.redobj.redconfig;
% errors = zeros(Npop, 1);
% for i = 1:Npop
%     [obj, ~, ~, ~, ~] = objfun(model.t_ref, X_ref_var_ss{i}, virtual_pop_par_initials_ss(i, :)', virtual_pop_par(i, :)', model.I, [], model.param, model.multiple, model.odefun, model.jacfun, redconfig_8_state, "MRSE");
%     errors(i) = obj.errout;
% end
% mean(errors)
% sum(errors < 0.1)

%%%% DO NEW NORMAL REDUCTION TO HAVE SEQUENCE FOR BACKTRACKING
% morexh_model(model, ...
%     'BC_SV40_from_JKn_2024', ...
%     'from_start', ...
%     0.1, ...
%     inf, ...
%     100, ...
%     120, ...
%     'linear', ...
%     'MRSE', ...
%     0, ...
%     0, ...
%     0, ...
%     ["dyn", "pneg", "env"], ...
%     0, ...
%     0, ...
%     0, ...
%     0, ...
%     [], ...
%     [], ...
%     [], ...
%     0, ...
%     [], [], 0)

%%%% DO NEW GREEDY REDUCTION WITH VARIABILITY
morexh_model(model, ...
    'BC_SV40_from_JKn_2024', ...
    'from_start', ...
    0.1, ...
    inf, ...
    100, ...
    120, ...
    'linear', ...
    'MRSE', ...
    0, ...
    0, ...
    0, ...
    ["dyn", "pneg", "env"], ...
    1, ...
    90, ...
    0, ...
    0, ...
    virtual_pop_par_initials_ss, ...
    virtual_pop_par, ...
    X_ref_var_ss, ...
    0, ...
    [], [], 0)


%% BACKTRACK
% indexing
I = model.I;
model.L = [];

load("BC_SV40_from_JKn_2024_exh_t120_MRSE_0.1_Inf_linear_dynpnegenv.mat")

% initiate first config
% exhaustive_mor.configs = repmat("dyn", 1, I.nstates);
% exhaustive_mor.configs(logical(model.state_unimportant)) = "pneg";
mor_options = redmodel.mor_options;
mor_options.err_out = 0.1;
mor_options.err_int = inf;

mor_options.variability = 1;
mor_options.var_obj_prctile = 95;
mor_options.virtual_pop_par = virtual_pop_par;
mor_options.virtual_pop_X0 = virtual_pop_par_initials_ss;
mor_options.X_ref_var = X_ref_var_ss;

validindices = find(redmodel.exhaustive_mor.objvals(:, 2) < mor_options.err_out & redmodel.exhaustive_mor.objvals(:, 3) < mor_options.err_int);
best_idx = validindices(end);
lastconfig = redmodel.exhaustive_mor.configs(best_idx, :);
exhaustive_mor.configs(1, :) = lastconfig;

% initiate exhaustive_mor
if mor_options.variability
    exhaustive_mor.objvals = [I.nstates 0 0 0 0 0];
else
    exhaustive_mor.objvals = [I.nstates 0 0 0];
end
exhaustive_mor.criterion = 0;
exhaustive_mor.statefromtoreduced = {0 "dyn" "dyn" 0};
exhaustive_mor.time = 0;

% make backwards mor
backwards_mor = exhaustive_mor;
backwards_mor.morexh_iteration = best_idx;

[backwards_mor, ~] = mor_exh_backwards(model, redmodel, backwards_mor, best_idx, mor_options, []);

%%%% DO NEW GREEDY REDUCTION WITH VARIABILITY

%%%% DO NEW NORMAL REDUCTION TO HAVE SEQUENCE FOR BACKTRACKING

%% errors of backtracked model

load("BC_SV40_from_JKn_2024_exh_t120_MRSE_0.1_Inf_linear_dynpnegenv_bintermediate.mat")
redconfig_24_state_backwards = backwards_mor.configs(end, :);
errors = zeros(Npop, 1);
for i = 1:Npop
    [obj, ~, ~, ~, ~] = objfun(model.t_ref, X_ref_var_ss{i}, virtual_pop_par_initials_ss(i, :)', virtual_pop_par(i, :)', model.I, [], model.param, model.multiple, model.odefun, model.jacfun, redconfig_24_state_backwards, "MRSE");
    errors(i) = obj.errout;
end
max(errors)
mean(errors)
sum(errors < 0.1)

%% remaining states

config = backwards_mor.configs(end, :);

model.I.nmstate(config == "dyn")