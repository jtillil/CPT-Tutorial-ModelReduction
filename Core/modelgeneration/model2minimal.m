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
    if bigmodel.multiple.multiple
        model.multiple = bigmodel.multiple;
        model.multiple.multiple = true;
        model.u_ref = bigmodel.u_ref;
    else
        model.multiple.multiple = false;
    end
else
    model.multiple.multiple = false;
end
if isfield(bigmodel, 'L')
    model.L = bigmodel.L;
else
    model.L = [];
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
if ~isfield(bigmodel, 'irenv_arith') || ~isfield(bigmodel, 'irenv_geom')
    if isfield(bigmodel, 'obs')
        model = model2irenv(model);
    end
end
if isfield(bigmodel, 'adjmat')
    model.adjmat = bigmodel.adjmat;
    model.speciesmask = bigmodel.speciesmask;
    model.substratemask = bigmodel.substratemask;
end

%% special fields for backwards compatibility
if isfield(bigmodel, 'input')
    model.input = bigmodel.input;
end
if isfield(bigmodel, 'simODE')
    model.simODE = bigmodel.simODE;
end
if isfield(bigmodel, 'relerrnorm')
    model.relerrnorm = bigmodel.relerrnorm;
end

%% generate matlab functions
if nargin(bigmodel.odefun) > 2
    [model.odefun, model.ode, model.jacfun, model.jac] = ode2matlabfun(bigmodel);
else
    model.odefun = bigmodel.odefun;
    model.ode = bigmodel.ode;
    model.jacfun = bigmodel.jacfun;
    model.jac = bigmodel.jac;
end

model.X_sym = cell2sym(model.I.nmstate(:));
model.par_sym = cell2sym(model.I.nmpar(:));
model.odefun_sym = model.odefun(model.X_sym, model.par_sym);
model.jacfun_sym = model.jacfun(model.X_sym, model.par_sym);

%% param struct
% if isfield(bigmodel, 'param')
%     model.param = bigmodel.param;
% else
    odefun_sym = model.odefun_sym;

    model.param.states2Ipar = {};
    model.param.states2nmpar = {};
    model.param.nreact = model.I.npar;
    
    I = model.I;
    paridx = 1:length(model.par);
    
    tic
    for k = 1:I.nstates
        states2param = [];
        % current_symvars = symvar(odefun_sym(k));
        current_symvars = identifyMultiplicativeFactors(odefun_sym(k));
        for current_symvar = current_symvars
            par_having_current_symvar = has(model.par_sym, current_symvar);
            if any(par_having_current_symvar)
                states2param = [states2param, paridx(par_having_current_symvar)];
            end
        end
        model.param.states2Ipar{k} = states2param;
        model.param.states2nmpar{k} = [I.nmpar(states2param)];
    end
    disp(toc)
% end

%% identify constant states (based on ODE = 0)

% create boolean vector for constant states
state_constant = zeros([length(model.X_sym), 1]);

% identify constant ODEs
tic
for i = 1:length(model.odefun_sym)
    if model.odefun_sym(i) == 0
        state_constant(i) = 1;
    end
end
disp(toc)

% append boolean vector to model
model.state_constant = state_constant;

%% identify unimportant states

% create boolean vector for unnecessary states
state_unimportant = zeros([length(model.X_sym), 1]);

% step 1: find states that do not contribute to other ODEs at all
tic
for i = 1:length(model.X_sym)
    if ~any(has(model.odefun_sym(setdiff(1:length(model.X_sym), i)), model.X_sym(i)))
        state_unimportant(i) = 1;
    end
end
disp(toc)

% step 2: find states that contribute to other ODEs only with all
% corresponding parameters = 0
tic
for i = 1:length(model.X_sym)
    % disp(X_sym(i))
    not_present_or_0 = zeros([length(model.X_sym), 1]);
    not_present_or_0(i) = 1;
    for j = setdiff(1:length(model.X_sym), i)
        other_ODE = model.odefun_sym(j);
        if ~has(other_ODE, model.X_sym(i))
            not_present_or_0(j) = 1;
        else
            no_influence_on_all_terms = 1;
            for term = children(other_ODE)
                if has(term{1}, model.X_sym(i))
                    variables = symvar(term{1});
                    for var = variables
                        if any(has(model.par_sym, var))
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
if isfield(bigmodel, 'obs')
    model = model2consts(model);
end

%% re-calculate reference solution

config = repmat("dyn", [1 model.I.nstates]);
[model.t_ref, model.X_ref] = simModel_simple(model, config);

%% calculate gramian(s)
% model = model2gramian(model);

end