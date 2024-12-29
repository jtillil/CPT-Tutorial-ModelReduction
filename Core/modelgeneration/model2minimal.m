function model = model2minimal(bigmodel)

%% basic fields

model = [];
model.name = bigmodel.name;
if isfield(bigmodel, 'scenario')
    model.scenario = bigmodel.scenario;
end
model.I = bigmodel.I;
model.X0 = bigmodel.X0;
model.par = bigmodel.par;
if isfield(bigmodel, 'multiple')
    model.multiple = bigmodel.multiple;
    model.multiple.multiple = true;
    model.u_ref = bigmodel.u_ref;
else
    model.multiple.multiple = false;
end
if isfield(bigmodel, 'L')
    model.L = bigmodel.L;
end
model.t_ref = bigmodel.t_ref;
model.X_ref = bigmodel.X_ref;
% model.yfun = @(model, xt) xt TODO ;
% model.ode = bigmodel.ode;
% model.jac = bigmodel.jac;
% model.odefun = bigmodel.odefun;
% model.jacfun = bigmodel.jacfun;
% model.simODE = bigmodel.simODE;
% model.abserrnorm = bigmodel.abserrnorm;
% model.relerrnorm = bigmodel.relerrnorm;
if isfield(bigmodel, 'analysis')
    model.analysis = bigmodel.analysis;
end
if isfield(bigmodel, 'ir')
    model.ir = bigmodel.ir;
end
if isfield(bigmodel, 'contr')
    model.contr = bigmodel.contr;
end
if isfield(bigmodel, 'obs')
    model.obs = bigmodel.obs;
end
if isfield(bigmodel, 'env')
    model.env = bigmodel.env;
end
if isfield(bigmodel, 'pneg')
    model.pneg = bigmodel.pneg;
end
if isfield(bigmodel, 'cneg')
    model.cneg = bigmodel.cneg;
end
if isfield(bigmodel, 'pss')
    model.pss = bigmodel.pss;
end
if isfield(bigmodel, 'irenv_arith')
    model.I.Irenv_arith = bigmodel.irenv_arith;
else
    model.irenv_arith = [];
end
if isfield(bigmodel, 'irenv_geom')
    model.I.Irenv_geom = bigmodel.irenv_geom;
else
    model.irenv_geom = [];
end
if ~isfield(bigmodel, 'irenv_arith') && ~isfield(bigmodel, 'irenv_geom')
    model = model2irenv(bigmodel);
end
if isfield(bigmodel, 'adjmat')
    model.adjmat = bigmodel.adjmat;
    model.speciesmask = bigmodel.speciesmask;
    model.substratemask = bigmodel.substratemask;
end

%% special old notation fields
if isfield(bigmodel, 'input')
    model.input = bigmodel.input;
end

%% generate matlab functions 
% create symbolic variables
syms t
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));

% generate ode matlabfun
tic
odefun_symbolic = bigmodel.odefun(t,X_sym,par_sym,model);
disp(toc)
% disp(odefun_symbolic)
tic
model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
% model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym,t},'File',['../modelfiles/' model.name '_ode.m'],'Optimize',true);
disp(toc)

% check jac presence
jac_present = true;
if ~isfield(bigmodel, 'jac')
    jac_present = false;
elseif isempty(bigmodel.jac)
    jac_present = false;
end

% generate jac matlabfun
tic
if jac_present
    jacfun_symbolic = bigmodel.jacfun(t,X_sym,par_sym,model);
else
    jacfun_symbolic = sym(zeros(model.I.nstates, model.I.nstates));
    for i = 1:model.I.nstates
        for j = 1:model.I.nstates
            jacfun_symbolic(i, j) = diff(odefun_symbolic(i), X_sym(j));
        end
    end
end
disp(toc)
tic
model.jacfun = matlabFunction(jacfun_symbolic,'Vars',{X_sym,par_sym,t});
% model.jacfun = matlabFunction(jacfun_symbolic,'Vars',{X_sym,par_sym,t},'File',['../modelfiles/' model.name '_jac.m'],'Optimize',true);
disp(toc)

%% param struct
if isfield(bigmodel, 'param')
    model.param = bigmodel.param;
else
    model.param.states2Ipar = {};
    model.param.states2nmpar = {};
    model.param.nreact = model.I.npar;
    
    I = model.I;
    paridx = 1:length(model.par);
    
    tic
    for k = 1:I.nstates
        states2param = [];
        current_symvars = symvar(odefun_symbolic(k));
        for current_symvar = current_symvars
            par_having_current_symvar = has(par_sym, current_symvar);
            if any(par_having_current_symvar)
                states2param = [states2param, paridx(par_having_current_symvar)];
            end
        end
        model.param.states2Ipar{k} = states2param;
        model.param.states2nmpar{k} = [I.nmpar(states2param)];
    end
    disp(toc)
end

%% identify constant states TODO ???
% create boolean vector for constant states
% state_constant = zeros([length(X_sym), 1]);

%% identify unimportant states
% create boolean vector for unnecessary states
state_unimportant = zeros([length(X_sym), 1]);

% step 1: find states that do not contribute to other ODEs at all
tic
for i = 1:length(X_sym)
    if ~any(has(odefun_symbolic(setdiff(1:length(X_sym), i)), X_sym(i)))
        state_unimportant(i) = 1;
    end
end
disp(toc)

% step 2: find states that contribute to other ODEs only with all
% corresponding parameters = 0
tic
for i = 1:length(X_sym)
    % disp(X_sym(i))
    not_present_or_0 = zeros([length(X_sym), 1]);
    not_present_or_0(i) = 1;
    for j = setdiff(1:length(X_sym), i)
        other_ODE = odefun_symbolic(j);
        if ~has(other_ODE, X_sym(i))
            not_present_or_0(j) = 1;
        else
            no_influence_on_all_terms = 1;
            for term = children(other_ODE)
                if has(term{1}, X_sym(i))
                    variables = symvar(term{1});
                    for var = variables
                        if any(has(par_sym, var))
                            par_val = model.par(model.I.(string(var)));
                            if par_val ~= 0
                                no_influence_on_all_terms = 0;
                            end
                        end
                    end
                end
            end
            if no_influence_on_all_terms
                not_present_or_0(j) = 1;
            end
        end
    end
    if all(not_present_or_0)
        state_unimportant(i) = 1;
    end
end
disp(toc)
disp(sum(state_unimportant))

% make sure, output is not unimportant
state_unimportant(model.I.output) = 0;

% append boolean vector to model
model.state_unimportant = state_unimportant;

% add constant vectors to I
model = model2consts(model);

%% re-calculate reference solution

%% calculate gramian(s)
% model = model2gramian(model);

end