%% Setup
clear; clc;
addpath(genpath("../../CPT-Tutorial-ModelReduction"))
reduced_errors = struct;

%% load good reduced model
load('modelEGFR_exh_t120_MRSE_0.1_0.5_linear_dyncnegpnegenvirenv_geompss.mat')
model = redmodel;
config = model.exhaustive_mor.configs(model.exhaustive_mor.objvals(:, 2) == model.redobj.errout, :);

%% code pss states into index-reduced model

% obtain pss states
I_red = config2I(model.I, config, []);

% create symbolic variables
syms t;
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));

% make odefun symbolic
tic
odefun_symbolic = model.odefun(X_sym,par_sym);
toc

nm_X_pss = model.I.nmstate(config == "pss");
X_sym_pss = X_sym(config == "pss");
odefun_symbolic_pss = odefun_symbolic(config == "pss");

% condition ODEs
odefun_symbolic_solved = odefun_symbolic;
odefun_symbolic_solved(config == "pneg") = 0;
odefun_symbolic_solved(config == "cneg") = 0;
odefun_symbolic_solved(config == "env") = 0;
odefun_symbolic_solved(config == "irenv_geom") = 0;
odefun_symbolic_solved(config == "irenv_arith") = 0;
odefun_symbolic_solved = simplify(odefun_symbolic_solved);

% % solve pss states
% tic
% G = solve(odefun_symbolic_pss, X_sym_pss, 'IgnoreAnalyticConstraints', true);
% toc
% 
% % input solved states to ODEs
% for k = 1:length(X_sym_pss)
%     solutions = G.(nm_X_pss{k});
%     odefun_symbolic_solved = subs(odefun_symbolic, X_sym_pss(k), solutions(1));
% end
% 
% % simplify solved ODEs
% odefun_symbolic_solved = simplify(odefun_symbolic_solved);

% solve pss states iteratively
id_states = 1:model.I.nstates;
id_pss = id_states(config == "pss");
for k = 1:length(X_sym_pss)
    disp(k)
    tic
    G = solve(odefun_symbolic_pss(k), X_sym_pss(k), 'IgnoreAnalyticConstraints', false);
    toc
    
    % input solved states to ODEs
    tic
    solutions = G(end);
    odefun_symbolic_solved = simplify(subs(odefun_symbolic_solved, X_sym_pss(k), solutions(1)));
    toc
    
    % remove solved state ODEs
    odefun_symbolic_solved(id_pss(k)) = 0;
end

% differentiate jac with solved ODEs
tic
jacfun_symbolic_solved = sym(zeros(model.I.nstates, model.I.nstates));
for i = 1:model.I.nstates
    for j = 1:model.I.nstates
        jacfun_symbolic_solved(i, j) = diff(odefun_symbolic(i), X_sym(j));
    end
end
toc

% convert to matlabfun and insert to model
model.odefun = matlabFunction(odefun_symbolic_solved,'Vars',{X_sym,par_sym});
model.ode = model.odefun;
model.jacfun = matlabFunction(jacfun_symbolic_solved,'Vars',{X_sym,par_sym,t});
model.jac = model.jacfun;

% update I
model.I.pss_solved = model.I.pss;
model.I.pss = [];

%% check pss solved greedy reduced EGFR model

config(config == "pss") = "pneg";
model.I = config2I(model.I, config, []);
[err_index_solved, ~, tred_solved, Xred_solved] = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, [], model.param, multiple, model.odefun, model.jacfun, config, "MRSE");
