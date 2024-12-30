%% Setup
clear; clc;
addpath(genpath("../../../CPT-Tutorial-ModelReduction"))
reduced_errors = struct;

%% Scenario 1 - no crosstalk
% load model file
load("modelMMEK_Cpss_Eenv_full.mat")
model.multiple.multiple = 0;

% lumping Cpss Eenv 2 (S), 3 (E)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(2, model.I.S) = 1;
lumpmat(2, model.I.E) = 1;
lumpmat(3, model.I.C) = 1;
lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_S_E(1,1) = calculate_lumping_error(model, lumpmat);

% lumping Cpss EpCenv 1 (A), 2 (S)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(1, model.I.S) = 1;
lumpmat(2, model.I.E) = 1;
lumpmat(3, model.I.C) = 1;
lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_A_S(1,1) = calculate_lumping_error(model, lumpmat);

% lumping indices 3 (E), 4 (C)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(2, model.I.S) = 1;
lumpmat(3, model.I.E) = 1;
lumpmat(3, model.I.C) = 1;
lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_E_C(1,1) = calculate_lumping_error(model, lumpmat);

% A env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.A) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.A_env(1,1) = obj.errout;

% S env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.S_env(1,1) = obj.errout;

% S pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.S_pss(1,1) = obj.errout;

% C env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.C_env(1,1) = obj.errout;

% C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.C_pss(1,1) = obj.errout;

% E env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.E) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.E_env(1,1) = obj.errout;

% E pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.E) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.E_pss(1,1) = obj.errout;

%% Scenario 2 - with crosstalk
% load model file
load("modelMMEK_Cpss_EpCenv_full.mat")
model.multiple.multiple = 0;

% lumping Cpss Eenv 2 (S), 3 (E)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(2, model.I.S) = 1;
lumpmat(2, model.I.E) = 1;
lumpmat(3, model.I.C) = 1;
lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_S_E(1,2) = calculate_lumping_error(model, lumpmat);

% lumping Cpss EpCenv 1 (A), 2 (S)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(1, model.I.S) = 1;
lumpmat(2, model.I.E) = 1;
lumpmat(3, model.I.C) = 1;
lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_A_S(1,2) = calculate_lumping_error(model, lumpmat);

% lumping indices 3 (E), 4 (C)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(2, model.I.S) = 1;
lumpmat(3, model.I.E) = 1;
lumpmat(3, model.I.C) = 1;
lumpmat(4, model.I.P) = 1;
reduced_errors.lumping_E_C(1,2) = calculate_lumping_error(model, lumpmat);

% A env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.A) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.A_env(1,2) = obj.errout;

% S env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.S_env(1,2) = obj.errout;

% S pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.S) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.S_pss(1,2) = obj.errout;

% C env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.C_env(1,2) = obj.errout;

% C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.C_pss(1,2) = obj.errout;

% E env
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.E) = "env";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.E_env(1,2) = obj.errout;

% E pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.E) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.E_pss(1,2) = obj.errout;

%% Save results
save("./results/errors_MMEK.mat", "reduced_errors")
