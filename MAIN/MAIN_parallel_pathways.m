%% Setup
clear; clc;
addpath(genpath("../../../CPT-Tutorial-ModelReduction"))
reduced_errors = struct;

%% Scenario 1 - no crosstalk
% load model file
load("modelSPP_no_crosstalk_full.mat")
model.multiple.multiple = 0;
config = repmat("dyn", [1 model.I.nstates]);
model.I = config2I(model.I, config, model.L);

% lumping merge both pathways 3 (B), 4 (C)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(2, model.I.S) = 1;
lumpmat(3, model.I.B) = 1;
lumpmat(3, model.I.C) = 1;
lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_B_C(1,1) = calculate_lumping_error(model, lumpmat);

% lumping no crosstalk 4 (C), 5 (D)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(2, model.I.S) = 1;
lumpmat(3, model.I.B) = 1;
lumpmat(4, model.I.C) = 1;
lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_C_D(1,1) = calculate_lumping_error(model, lumpmat);

% lumping with crosstalk 2 (S), 3 (B)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(2, model.I.S) = 1;
lumpmat(2, model.I.B) = 1;
lumpmat(3, model.I.C) = 1;
lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_S_B(1,1) = calculate_lumping_error(model, lumpmat);

% lumping try S C 2 (S), 4 (C)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(2, model.I.S) = 1;
lumpmat(3, model.I.B) = 1;
lumpmat(2, model.I.C) = 1;
lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_S_C(1,1) = calculate_lumping_error(model, lumpmat);

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

% C cneg
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "cneg";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.C_cneg(1,1) = obj.errout;

% C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.C_pss(1,1) = obj.errout;

% B cneg
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "cneg";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.B_cneg(1,1) = obj.errout;

% B pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.B_pss(1,1) = obj.errout;

%% Scenario 2 - with crosstalk
% load model file
load("modelSPP_with_crosstalk_full.mat")
model.multiple.multiple = 0;
config = repmat("dyn", [1 model.I.nstates]);
model.I = config2I(model.I, config, model.L);

% lumping merge both pathways 3 (B), 4 (C)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(2, model.I.S) = 1;
lumpmat(3, model.I.B) = 1;
lumpmat(3, model.I.C) = 1;
lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_B_C(1,2) = calculate_lumping_error(model, lumpmat);

% lumping no crosstalk 4 (C), 5 (D)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(2, model.I.S) = 1;
lumpmat(3, model.I.B) = 1;
lumpmat(4, model.I.C) = 1;
lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_C_D(1,2) = calculate_lumping_error(model, lumpmat);

% lumping with crosstalk 2 (S), 3 (B)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(2, model.I.S) = 1;
lumpmat(2, model.I.B) = 1;
lumpmat(3, model.I.C) = 1;
lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_S_B(1,2) = calculate_lumping_error(model, lumpmat);

% lumping try S C 2 (S), 4 (C)
lumpmat = zeros(4, 5);
lumpmat(1, model.I.A) = 1;
lumpmat(2, model.I.S) = 1;
lumpmat(3, model.I.B) = 1;
lumpmat(2, model.I.C) = 1;
lumpmat(4, model.I.D) = 1;
reduced_errors.lumping_S_C(1,2) = calculate_lumping_error(model, lumpmat);

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

% C cneg
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "cneg";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.C_cneg(1,2) = obj.errout;

% C pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.C) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.C_pss(1,2) = obj.errout;

% B cneg
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "cneg";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.B_cneg(1,2) = obj.errout;

% B pss
config = repmat("dyn", [1 model.I.nstates]);
config(model.I.B) = "pss";
obj = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");
reduced_errors.B_pss(1,2) = obj.errout;

%% Save results
save("./results/errors_SPP.mat", "reduced_errors")
