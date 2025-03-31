%% load good reduced model

load('modelEGFR_exh_t120_0.1_0.5_linear_dyncnegpnegirenv_geomenvpss.mat')
model = redmodel;
model.redobj.redconfig = model.exhaustive_mor.configs(model.exhaustive_mor.objvals(:, 2) == model.redobj.errout, :);

%% code pss states into index-reduced model

% obtain pss states
I_red = config2I(model.I, redmodel.exhasutice_mr, []);

% create symbolic variables
syms t;
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));

% make odefun symbolic
tic
odefun_symbolic = model.odefun(X_sym,par_sym);
disp(toc)

nm_X_pss = model.I.nmstate(config == "pss");
X_sym_pss = X_sym(config == "pss");
odefun_symbolic_pss = odefun_symbolic(config == "pss");

% solve pss states
G = solve(odefun_symbolic_pss, X_sym_pss);

% input solved states to ODEs
for k = 1:length(X_sym_pss)
    solutions = G.(nm_X_pss{k});
    odefun_symbolic_solved = subs(odefun_symbolic, X_sym_pss(k), solutions(1));
end

% remove solved state ODEs
odefun_symbolic_solved(config == "pss") = 0;

% simplify solved ODEs
odefun_symbolic_solved = simplify(odefun_symbolic_solved);

% differentiate jac with solved ODEs
tic
jacfun_symbolic_solved = sym(zeros(model.I.nstates, model.I.nstates));
for i = 1:model.I.nstates
    for j = 1:model.I.nstates
        jacfun_symbolic_solved(i, j) = diff(odefun_symbolic(i), X_sym(j));
    end
end
disp(toc)

% convert to matlabfun and insert to model
model.odefun = matlabFunction(odefun_symbolic_solved,'Vars',{X_sym,par_sym});
model.ode = model.odefun;
model.jacfun = matlabFunction(jacfun_symbolic_solved,'Vars',{X_sym,par_sym,t});
model.jac = model.jacfun;

% update I
model.I.pss_solved = model.I.pss;
model.I.pss = [];