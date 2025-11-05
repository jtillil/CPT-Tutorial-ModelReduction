%% Setup
clear; clc;
addpath(genpath("../../../CPT-Tutorial-ModelReduction"))
reduced_errors = struct;

load("modelMMEK_Cpss_Eenv_full.mat")
AUC_Cpss_Eenv = trapz(model.t_ref, model.X_ref, 1);

load("modelMMEK_Cpss_EpCenv_full.mat")
AUC_Cpss_EpCenv = trapz(model.t_ref, model.X_ref, 1);

%% Scenario 1
% load model file
load("modelMMEK_Cpss_Eenv_full.mat")
model.multiple.multiple = 0;
config = repmat("dyn", [1 model.I.nstates]);
model.I = config2I(model.I, config, model.L);

% lumping Cpss Eenv 2 (S), 3 (E)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(2, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_S_E(1,1) = calculate_lumping_error(model, opt);

% lumping Cpss EpCenv 1 (A), 2 (S)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(1, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_A_S(1,1) = calculate_lumping_error(model, opt);

% scaled lumping Cpss EpCenv 1 (A), 2 (S)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(1, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
opt.invlumpmat = opt.lumpmat';
for col = 1:size(opt.invlumpmat, 2)
    opt.invlumpmat(:, col) = (AUC_Cpss_Eenv.*opt.invlumpmat(:, col)') ./ (AUC_Cpss_Eenv*opt.invlumpmat(:, col));
end
reduced_errors.scaled_lumping_Cpss_Eenv_A_S(1,1) = calculate_lumping_error(model, opt);

% scaled lumping Cpss EpCenv 1 (A), 2 (S)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(1, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
opt.invlumpmat = opt.lumpmat';
for col = 1:size(opt.invlumpmat, 2)
    opt.invlumpmat(:, col) = (AUC_Cpss_EpCenv.*opt.invlumpmat(:, col)') ./ (AUC_Cpss_EpCenv*opt.invlumpmat(:, col));
end
reduced_errors.scaled_lumping_Cpss_EpCenv_A_S(1,1) = calculate_lumping_error(model, opt);

% lumping indices 3 (E), 4 (C)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(2, model.I.S) = 1;
opt.lumpmat(3, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_E_C(1,1) = calculate_lumping_error(model, opt);

% A env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.A) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.A_env(1,1) = obj.errout;

% S env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.S_env(1,1) = obj.errout;

% S pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.S_pss(1,1) = obj.errout;

% C env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_env(1,1) = obj.errout;

% C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_pss(1,1) = obj.errout;

% E env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.E) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.E_env(1,1) = obj.errout;

% E pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.E) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.E_pss(1,1) = obj.errout;

% C pss E env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
config(model.I.E) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_pss_E_env(1,1) = obj.errout;

%% Scenario 1 - C pss E env ir-indices
load("modelMMEK_Cpss_Eenv_full.mat")
model.multiple.multiple = 0;
config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.C) = "pss";
% config(model.I.E) = "env";
model.I = config2I(model.I, config, []);
[irred.ir, irred.contr, irred.obs, irred.t_ir] = compute_ir_indices_matlabfun(model);

size = 8;
lw = 1;
lwt = 0.5;

% nir-indices
figure
% plot(model.t_ref, irred.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
plot(irred.t_ir, irred.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-1 31])
ylim([-0.01 1])
legend('A', 'S', 'E', 'C', 'P', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_Eenv_ir.pdf")

%%%% new
load("modelMMEK_Cpss_Eenv_full.mat")
model.multiple.multiple = 0;
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
% config(model.I.E) = "env";
model.I = config2I(model.I, config, []);
[irred.ir, irred.contr, irred.obs, irred.t_ir] = compute_ir_indices_matlabfun(model);

size = 8;
lw = 1;
lwt = 0.5;

% nir-indices
figure
% plot(model.t_ref, irred.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
plot(irred.t_ir, irred.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-1 31])
ylim([-0.01 1])
legend('A', 'S', 'E', 'C', 'P', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_Eenv_ir_C_pss.pdf")

%% Scenario 2 - C pss E pCenv ir-indices

%%%% full model
load("modelMMEK_Cpss_EpCenv_full.mat")
model.multiple.multiple = 0;
config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.C) = "pss";
% config(model.I.A) = "env";
model.I = config2I(model.I, config, []);
[irred.ir, irred.contr, irred.obs, irred.t_ir] = compute_ir_indices_matlabfun(model);

size = 8;
lw = 1;
lwt = 0.5;

% nir-indices
figure
% plot(model.t_ref, irred.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
plot(irred.t_ir, irred.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-1 31])
ylim([-0.01 1])
legend('A', 'S', 'E', 'C', 'P', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_EpCenv_ir.pdf")

%%%% C in pss
load("modelMMEK_Cpss_EpCenv_full.mat")
model.multiple.multiple = 0;
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
% config(model.I.A) = "env";
model.I = config2I(model.I, config, []);
[irred.ir, irred.contr, irred.obs, irred.t_ir] = compute_ir_indices_matlabfun(model);

size = 8;
lw = 1;
lwt = 0.5;

% nir-indices
figure
% plot(model.t_ref, irred.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
plot(irred.t_ir, irred.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-1 31])
ylim([-0.01 1])
legend('A', 'S', 'E', 'C', 'P', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_EpCenv_ir_C_pss.pdf")

%% Scenario 2
% load model file
load("modelMMEK_Cpss_EpCenv_full.mat")
model.multiple.multiple = 0;
config = repmat("dyn", [1 model.I.nstates]);
model.I = config2I(model.I, config, model.L);

% lumping Cpss Eenv 2 (S), 3 (E)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(2, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_S_E(1,2) = calculate_lumping_error(model, opt);

% lumping Cpss EpCenv 1 (A), 2 (S)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(1, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_A_S(1,2) = calculate_lumping_error(model, opt);

% scaled lumping Cpss EpCenv 1 (A), 2 (S)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(1, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
opt.invlumpmat = opt.lumpmat';
for col = 1:size(opt.invlumpmat, 2)
    opt.invlumpmat(:, col) = (AUC_Cpss_Eenv.*opt.invlumpmat(:, col)') ./ (AUC_Cpss_Eenv*opt.invlumpmat(:, col));
end
reduced_errors.scaled_lumping_Cpss_Eenv_A_S(1,2) = calculate_lumping_error(model, opt);

% scaled lumping Cpss EpCenv 1 (A), 2 (S)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(1, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
opt.invlumpmat = opt.lumpmat';
for col = 1:size(opt.invlumpmat, 2)
    opt.invlumpmat(:, col) = (AUC_Cpss_EpCenv.*opt.invlumpmat(:, col)') ./ (AUC_Cpss_EpCenv*opt.invlumpmat(:, col));
end
reduced_errors.scaled_lumping_Cpss_EpCenv_A_S(1,2) = calculate_lumping_error(model, opt);

% lumping indices 3 (E), 4 (C)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(2, model.I.S) = 1;
opt.lumpmat(3, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_E_C(1,2) = calculate_lumping_error(model, opt);

% A env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.A) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.A_env(1,2) = obj.errout;

% S env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.S_env(1,2) = obj.errout;

% S pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.S_pss(1,2) = obj.errout;

% C env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_env(1,2) = obj.errout;

% C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_pss(1,2) = obj.errout;

% E env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.E) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.E_env(1,2) = obj.errout;

% E pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.E) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.E_pss(1,2) = obj.errout;

% C pss E env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
config(model.I.E) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_pss_E_env(1,2) = obj.errout;

%% Scenario 1 with E in conservationlaw E+C
% load model file
load("modelMMEK_Cpss_Eenv_full.mat")
model.multiple.multiple = 0;

% add conservation law
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.E) = "con1";
model.I = config2I(model.I, config, model.L);

% lumping Cpss Eenv 2 (S), 3 (E)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(2, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_S_E(1,3) = calculate_lumping_error(model, opt);

% lumping Cpss EpCenv 1 (A), 2 (S)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(1, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_A_S(1,3) = calculate_lumping_error(model, opt);

% scaled lumping Cpss EpCenv 1 (A), 2 (S)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(1, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
opt.invlumpmat = opt.lumpmat';
for col = 1:size(opt.invlumpmat, 2)
    opt.invlumpmat(:, col) = (AUC_Cpss_Eenv.*opt.invlumpmat(:, col)') ./ (AUC_Cpss_Eenv*opt.invlumpmat(:, col));
end
reduced_errors.scaled_lumping_Cpss_Eenv_A_S(1,3) = calculate_lumping_error(model, opt);

% scaled lumping Cpss EpCenv 1 (A), 2 (S)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(1, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
opt.invlumpmat = opt.lumpmat';
for col = 1:size(opt.invlumpmat, 2)
    opt.invlumpmat(:, col) = (AUC_Cpss_EpCenv.*opt.invlumpmat(:, col)') ./ (AUC_Cpss_EpCenv*opt.invlumpmat(:, col));
end
reduced_errors.scaled_lumping_Cpss_EpCenv_A_S(1,3) = calculate_lumping_error(model, opt);

% lumping indices 3 (E), 4 (C)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(2, model.I.S) = 1;
opt.lumpmat(3, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_E_C(1,3) = calculate_lumping_error(model, opt);

% A env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.A) = "env";
config(model.I.E) = "con1";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.A_env(1,3) = obj.errout;

% S env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "env";
config(model.I.E) = "con1";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.S_env(1,3) = obj.errout;

% S pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "pss";
config(model.I.E) = "con1";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.S_pss(1,3) = obj.errout;

% C env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "env";
config(model.I.E) = "con1";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_env(1,3) = obj.errout;

% C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
config(model.I.E) = "con1";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_pss(1,3) = obj.errout;

% E env
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.E) = "env";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.E_env(1,3) = NaN;

% E pss
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.E) = "pss";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.E_pss(1,3) = NaN;

% C pss E env
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.C) = "pss";
% config(model.I.E) = "env";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_pss_E_env(1,3) = NaN;

%% Scenario 2 with E in conservationlaw E+C
% load model file
load("modelMMEK_Cpss_EpCenv_full.mat")
model.multiple.multiple = 0;

% add conservation law
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.E) = "con1";
model.I = config2I(model.I, config, model.L);

% lumping Cpss Eenv 2 (S), 3 (E)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(2, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_S_E(1,4) = calculate_lumping_error(model, opt);

% lumping Cpss EpCenv 1 (A), 2 (S)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(1, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_A_S(1,4) = calculate_lumping_error(model, opt);

% scaled lumping Cpss EpCenv 1 (A), 2 (S)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(1, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
opt.invlumpmat = opt.lumpmat';
for col = 1:size(opt.invlumpmat, 2)
    opt.invlumpmat(:, col) = (AUC_Cpss_Eenv.*opt.invlumpmat(:, col)') ./ (AUC_Cpss_Eenv*opt.invlumpmat(:, col));
end
reduced_errors.scaled_lumping_Cpss_Eenv_A_S(1,4) = calculate_lumping_error(model, opt);

% scaled lumping Cpss EpCenv 1 (A), 2 (S)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(1, model.I.S) = 1;
opt.lumpmat(2, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
opt.invlumpmat = opt.lumpmat';
for col = 1:size(opt.invlumpmat, 2)
    opt.invlumpmat(:, col) = (AUC_Cpss_EpCenv.*opt.invlumpmat(:, col)') ./ (AUC_Cpss_EpCenv*opt.invlumpmat(:, col));
end
reduced_errors.scaled_lumping_Cpss_EpCenv_A_S(1,4) = calculate_lumping_error(model, opt);

% lumping indices 3 (E), 4 (C)
opt.lumpmat = zeros(4, 5);
opt.lumpmat(1, model.I.A) = 1;
opt.lumpmat(2, model.I.S) = 1;
opt.lumpmat(3, model.I.E) = 1;
opt.lumpmat(3, model.I.C) = 1;
opt.lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_E_C(1,4) = calculate_lumping_error(model, opt);

% A env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.A) = "env";
config(model.I.E) = "con1";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.A_env(1,4) = obj.errout;

% S env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "env";
config(model.I.E) = "con1";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.S_env(1,4) = obj.errout;

% S pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "pss";
config(model.I.E) = "con1";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.S_pss(1,4) = obj.errout;

% C env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "env";
config(model.I.E) = "con1";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_env(1,4) = obj.errout;

% C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
config(model.I.E) = "con1";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_pss(1,4) = obj.errout;

% E env
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.E) = "env";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.E_env(1,4) = NaN;

% E pss
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.E) = "pss";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.E_pss(1,4) = NaN;

% C pss E env
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.C) = "pss";
% config(model.I.E) = "env";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_pss_E_env(1,4) = NaN;

%% Save results
save("./results/errors_MMEK.mat", "reduced_errors")
