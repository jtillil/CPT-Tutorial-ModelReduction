%% Setup
clear; clc;
addpath(genpath("../../CPT-Tutorial-ModelReduction"))
reduced_errors = struct;

%% Scenario 1 - no crosstalk
% load model file
load("modelSPP_no_crosstalk_full.mat")
model.multiple.multiple = 0;
config = repmat("dyn", [1 model.I.nstates]);
model.I = config2I(model.I, config, model.L);

% lumping merge both pathways 3 (B), 4 (C)
options.lumpmat = zeros(4, 5);
options.lumpmat(1, model.I.A) = 1;
options.lumpmat(2, model.I.S) = 1;
options.lumpmat(3, model.I.B) = 1;
options.lumpmat(3, model.I.C) = 1;
options.lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_B_C(1,1) = calculate_lumping_error(model, options);

% lumping no crosstalk 4 (C), 5 (D)
options.lumpmat = zeros(4, 5);
options.lumpmat(1, model.I.A) = 1;
options.lumpmat(2, model.I.S) = 1;
options.lumpmat(3, model.I.B) = 1;
options.lumpmat(4, model.I.C) = 1;
options.lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_C_D(1,1) = calculate_lumping_error(model, options);

% lumping with crosstalk 2 (S), 3 (B)
options.lumpmat = zeros(4, 5);
options.lumpmat(1, model.I.A) = 1;
options.lumpmat(2, model.I.S) = 1;
options.lumpmat(2, model.I.B) = 1;
options.lumpmat(3, model.I.C) = 1;
options.lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_S_B(1,1) = calculate_lumping_error(model, options);

% lumping try S C 2 (S), 4 (C)
options.lumpmat = zeros(4, 5);
options.lumpmat(1, model.I.A) = 1;
options.lumpmat(2, model.I.S) = 1;
options.lumpmat(3, model.I.B) = 1;
options.lumpmat(2, model.I.C) = 1;
options.lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_S_C(1,1) = calculate_lumping_error(model, options);

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

% C cneg
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "cneg";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_cneg(1,1) = obj.errout;

% C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_pss(1,1) = obj.errout;

% B cneg
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "cneg";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.B_cneg(1,1) = obj.errout;

% B pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.B_pss(1,1) = obj.errout;

% B ssenv
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "ssenv";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.B_ssenv(1,1) = obj.errout;

% B cneg C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "cneg";
config(model.I.C) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.B_cneg_C_pss(1,1) = obj.errout;

% B pss C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "pss";
config(model.I.C) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.B_pss_C_pss(1,1) = obj.errout;

% B cneg C cneg
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "cneg";
config(model.I.C) = "cneg";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.B_cneg_C_cneg(1,1) = obj.errout;

% S env C cneg
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "env";
config(model.I.C) = "cneg";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.S_env_C_cneg(1,1) = obj.errout;

%% Scenario 1 - no crosstalk C pss ir-indices
% load("modelSPP_no_crosstalk_full.mat")
% model.multiple.multiple = 0;
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.C) = "pss";
% % config(model.I.E) = "env";
% model.I = config2I(model.I, config, []);
% [irred.ir, irred.contr, irred.obs, irred.t_ir] = compute_ir_indices_matlabfun(model);
% 
% size = 8;
% lw = 1;
% lwt = 0.5;
% 
% % nir-indices
% figure
% % plot(model.t_ref, irred.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
% plot(irred.t_ir, irred.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
% yline(0.1, 'k--', 'LineWidth', lwt)
% xlim([-0.002 0.052])
% ylim([-0.01 1])
% legend('A', 'S', 'B', 'C', 'D', 'threshold', 'Location','northeast')
% xlabel("t [min]")
% ylabel("nir-index")
% 
% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]
% 
% exportgraphics(gcf, "./figures/SPP_no_crosstalk_ir_C_pss.pdf")

%% Scenario 2 - with crosstalk
% load model file
load("modelSPP_with_crosstalk_full.mat")
model.multiple.multiple = 0;
config = repmat("dyn", [1 model.I.nstates]);
model.I = config2I(model.I, config, model.L);

% lumping merge both pathways 3 (B), 4 (C)
options.lumpmat = zeros(4, 5);
options.lumpmat(1, model.I.A) = 1;
options.lumpmat(2, model.I.S) = 1;
options.lumpmat(3, model.I.B) = 1;
options.lumpmat(3, model.I.C) = 1;
options.lumpmat(4, model.I.D) = 1;
% reduced_errors.lumping_B_C(1,2) = calculate_lumping_error(model, options);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    lumping_error(model, options);

% lumping no crosstalk 4 (C), 5 (D)
options.lumpmat = zeros(4, 5);
options.lumpmat(1, model.I.A) = 1;
options.lumpmat(2, model.I.S) = 1;
options.lumpmat(3, model.I.B) = 1;
options.lumpmat(4, model.I.C) = 1;
options.lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_C_D(1,2) = calculate_lumping_error(model, options);

% lumping with crosstalk 2 (S), 3 (B)
options.lumpmat = zeros(4, 5);
options.lumpmat(1, model.I.A) = 1;
options.lumpmat(2, model.I.S) = 1;
options.lumpmat(2, model.I.B) = 1;
options.lumpmat(3, model.I.C) = 1;
options.lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_S_B(1,2) = calculate_lumping_error(model, options);

% lumping try S C 2 (S), 4 (C)
options.lumpmat = zeros(4, 5);
options.lumpmat(1, model.I.A) = 1;
options.lumpmat(2, model.I.S) = 1;
options.lumpmat(3, model.I.B) = 1;
options.lumpmat(2, model.I.C) = 1;
options.lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_S_C(1,2) = calculate_lumping_error(model, options);

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

% C cneg
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "cneg";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_cneg(1,2) = obj.errout;

% C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.C_pss(1,2) = obj.errout;

% B cneg
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "cneg";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.B_cneg(1,2) = obj.errout;

% B pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.B_pss(1,2) = obj.errout;

% B ssenv
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "ssenv";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.B_ssenv(1,2) = obj.errout;

% B cneg C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "cneg";
config(model.I.C) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.B_cneg_C_pss(1,2) = obj.errout;

% B pss C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "pss";
config(model.I.C) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.B_pss_C_pss(1,2) = obj.errout;

% B cneg C cneg
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "cneg";
config(model.I.C) = "cneg";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.B_cneg_C_cneg(1,2) = obj.errout;

% S env C cneg
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "env";
config(model.I.C) = "cneg";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
reduced_errors.S_env_C_cneg(1,2) = obj.errout;

%% Save results

save("./results/errors_SPP.mat", "reduced_errors")

%% Scenario 3 - C and B neglectable
% % load model file
% load("modelSPP_C_B_neglectable_full.mat")
% model.multiple.multiple = 0;
% config = repmat("dyn", [1 model.I.nstates]);
% model.I = config2I(model.I, config, model.L);
% 
% % lumping merge both pathways 3 (B), 4 (C)
% options.lumpmat = zeros(4, 5);
% options.lumpmat(1, model.I.A) = 1;
% options.lumpmat(2, model.I.S) = 1;
% options.lumpmat(3, model.I.B) = 1;
% options.lumpmat(3, model.I.C) = 1;
% options.lumpmat(4, model.I.D) = 1;
% reduced_errors.lumping_B_C(1,3) = calculate_lumping_error(model, options);
% 
% % lumping no crosstalk 4 (C), 5 (D)
% options.lumpmat = zeros(4, 5);
% options.lumpmat(1, model.I.A) = 1;
% options.lumpmat(2, model.I.S) = 1;
% options.lumpmat(3, model.I.B) = 1;
% options.lumpmat(4, model.I.C) = 1;
% options.lumpmat(4, model.I.D) = 1;
% reduced_errors.lumping_C_D(1,3) = calculate_lumping_error(model, options);
% 
% % lumping with crosstalk 2 (S), 3 (B)
% options.lumpmat = zeros(4, 5);
% options.lumpmat(1, model.I.A) = 1;
% options.lumpmat(2, model.I.S) = 1;
% options.lumpmat(2, model.I.B) = 1;
% options.lumpmat(3, model.I.C) = 1;
% options.lumpmat(4, model.I.D) = 1;
% reduced_errors.lumping_S_B(1,3) = calculate_lumping_error(model, options);
% 
% % lumping try S C 2 (S), 4 (C)
% options.lumpmat = zeros(4, 5);
% options.lumpmat(1, model.I.A) = 1;
% options.lumpmat(2, model.I.S) = 1;
% options.lumpmat(3, model.I.B) = 1;
% options.lumpmat(2, model.I.C) = 1;
% options.lumpmat(4, model.I.D) = 1;
% reduced_errors.lumping_S_C(1,3) = calculate_lumping_error(model, options);
% 
% % A env
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.A) = "env";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
% reduced_errors.A_env(1,3) = obj.errout;
% 
% % S env
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.S) = "env";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
% reduced_errors.S_env(1,3) = obj.errout;
% 
% % S pss
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.S) = "pss";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
% reduced_errors.S_pss(1,3) = obj.errout;
% 
% % C cneg
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.C) = "cneg";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
% reduced_errors.C_cneg(1,3) = obj.errout;
% 
% % C pss
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.C) = "pss";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
% reduced_errors.C_pss(1,3) = obj.errout;
% 
% % B cneg
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.B) = "cneg";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
% reduced_errors.B_cneg(1,3) = obj.errout;
% 
% % B pss
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.B) = "pss";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
% reduced_errors.B_pss(1,3) = obj.errout;
% 
% % B ssenv
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.B) = "ssenv";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
% reduced_errors.B_ssenv(1,3) = obj.errout;
% 
% % B cneg C pss
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.B) = "cneg";
% config(model.I.C) = "pss";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
% reduced_errors.B_cneg_C_pss(1,3) = obj.errout;
% 
% % B pss C pss
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.B) = "pss";
% config(model.I.C) = "pss";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
% reduced_errors.B_pss_C_pss(1,3) = obj.errout;
% 
% % B cneg C cneg
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.B) = "cneg";
% config(model.I.C) = "cneg";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
% reduced_errors.B_cneg_C_cneg(1,3) = obj.errout;
% 
% % S env C cneg
% config = repmat("dyn", [1 model.I.nstates]);
% config(model.I.S) = "env";
% config(model.I.C) = "cneg";
% obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "rel2NE");
% reduced_errors.S_env_C_cneg(1,3) = obj.errout;

%% reduced indices - no crosstalk
load("modelSPP_no_crosstalk_full.mat")
model.multiple.multiple = 0;

config = repmat("dyn", [model.I.nstates, 1]);
config(model.I.C) = "cneg";
config(model.I.S) = "env";
model.I = config2I(model.I, config, []);

[ir, contr, obs, t_ind] = compute_ir_indices_matlabfun(model);

save("./results/indices_SPP_no_crosstalk_red.mat", "ir", "contr", "obs", "t_ind")

%% reduced indices - with crosstalk
load("modelSPP_with_crosstalk_full.mat")
model.multiple.multiple = 0;

config = repmat("dyn", [model.I.nstates, 1]);
config(model.I.C) = "cneg";
model.I = config2I(model.I, config, []);

[ir, contr, obs, t_ind] = compute_ir_indices_matlabfun(model);

save("./results/indices_SPP_with_crosstalk_red_correct.mat", "ir", "contr", "obs", "t_ind")

config = repmat("dyn", [model.I.nstates, 1]);
config(model.I.B) = "cneg";
% config(model.I.C) = "cneg"
model.I = config2I(model.I, config, []);

[ir, contr, obs, t_ind] = compute_ir_indices_matlabfun(model);

save("./results/indices_SPP_with_crosstalk_red_wrong.mat", "ir", "contr", "obs", "t_ind")