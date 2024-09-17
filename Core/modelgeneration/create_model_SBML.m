%% load SBML model
modelSBML = TranslateSBML('../lib/libSBML-5/BIOMD0000000613_url.xml');

%% paths

addpath(genpath('modelspecification'))
addpath(genpath('modelfiles'))
addpath(genpath('Indices'))

%% model initialization

model.name = modelSBML.name;
model.scenario = '-';
model.SBML = true;

%% species

model.I.nstates = length({modelSBML.species.id});
model.I.nmstate = {modelSBML.species.id};
for i = 1:length(model.I.nmstate)
    model.I.(model.I.nmstate{i}) = i;
end

model.X0 = cell2mat({modelSBML.species.initialConcentration}');
model.X0(isnan(model.X0)) = 0;

%% parameters

model.I.nmpar = {modelSBML.parameter.id};
for i = 1:length(model.I.nmpar)
    model.I.(model.I.nmpar{i}) = i;
end

model.par = cell2mat({modelSBML.parameter.value}');
model.par(isnan(model.par)) = 0;

%% config

model = config2model(model, repmat("dyn", [1, model.I.nstates]));
model.I.replaceODE = [];

%% initialAssignment --> assign initials and parameters initially

model.initassign.initials_and_params = {modelSBML.initialAssignment.symbol}';
model.initassign.assigns = {modelSBML.initialAssignment.math}';

for i = 1:length(model.I.nmpar)
    assignin('base',model.I.nmpar{i},model.par(i))
end
for i = 1:length(model.I.nmstate)
    assignin('base',model.I.nmstate{i},model.X0(i))
end

tic
for i = 1:length(model.initassign.initials_and_params)
    current_assign = model.initassign.assigns{i};
    if any(contains(model.I.nmstate, model.initassign.initials_and_params{i}))
        model.X0(model.I.(model.initassign.initials_and_params{i})) = eval(current_assign);
        assignin('base', model.initassign.initials_and_params{i}, model.X0(model.I.(model.initassign.initials_and_params{i})))
    else
        model.par(model.I.(model.initassign.initials_and_params{i})) = eval(current_assign);
        assignin('base', model.initassign.initials_and_params{i}, model.par(model.I.(model.initassign.initials_and_params{i})))
    end
end
toc

%% rule --> specify not yet assigned parameters

model.rule.params = {modelSBML.rule.variable}';
model.rule.rules = {modelSBML.rule.formula}';

% apply specific changes
model = changes_to_model(model);

tic
for i = 1:length(model.rule.params)
    current_rule = model.rule.rules{i};
    % try
    model.par(model.I.(model.rule.params{i})) = eval(current_rule);
    assignin('base', model.rule.params{i}, model.par(model.I.(model.rule.params{i})))
    % catch
    %     try
    %         model.par(model.I.(model.rule.params{i})) = subs(str2sym(current_rule));
    %         assignin('base', model.rule.params{i}, model.par(model.I.(model.rule.params{i})))
    %     catch
    %         disp('Error evaluating rules.')
    %         stop()
    %         % model.par(model.I.(model.rule.params{i})) = create_and_execute_sym(model, current_rule);
    %         % assignin('base', model.rule.params{i}, model.par(model.I.(model.rule.params{i})))
    %     end
    % end
end
toc

%% ode function

reactant = {modelSBML.reaction.reactant}';
product = {modelSBML.reaction.product}';
% modifier = {modelSBML.reaction.modifier}';
law = {modelSBML.reaction.kineticLaw}';

X_sym = cell2sym(model.I.nmstate(:));
par_sym = cell2sym(model.I.nmpar(:));

% for i = 1:length(model.I.nmpar)
%     assignin('base',model.I.nmpar{i},par_sym(i))
% end
% for i = 1:length(model.I.nmstate)
%     assignin('base',model.I.nmstate{i},X_sym(i))
% end

odesym = sym(zeros(model.I.nstates, 1));

tic
for i = 1:length(law)
    % count and save reactants
    current_reactants = [];
    for reactantentry = reactant{i}
        current_reactants = [current_reactants, model.I.(reactantentry.species)];
    end
    % count and save products
    current_products = [];
    for productentry = product{i}
        current_products = [current_products, model.I.(productentry.species)];
    end
    % count and save modifiers
    % modifiers probably irrelevant
    % create law
    current_law = str2sym(law{i}.math);
    % append to odesym
    if ~isempty(current_reactants)
        odesym(current_reactants) = odesym(current_reactants) - current_law;
    end
    if ~isempty(current_products)
        odesym(current_products) = odesym(current_products) + current_law;
    end
end
toc

tic
model.odefun = matlabFunction(odesym,'Vars',{X_sym,par_sym});
toc

%% jac function

tic
jacsym = odesym2jacsym(odesym, X_sym, model);
toc

tic
model.jacfun = matlabFunction(jacsym,'Vars',{X_sym,par_sym});
toc

%% param struct for parameter involvement

model.param.states2Ipar = {};
model.param.states2nmpar = {};
model.param.nreact = length(law);

I = model.I;
paridx = 1:length(model.par);

tic
for k = 1:I.nstates
    states2param = [];
    current_symvars = symvar(odesym(k));
    for current_symvar = current_symvars
        par_having_current_symvar = has(par_sym, current_symvar);
        if any(par_having_current_symvar)
            states2param = [states2param, paridx(par_having_current_symvar)];
        end
    end
    model.param.states2Ipar{k} = states2param;
    model.param.states2nmpar{k} = [I.nmpar(states2param)];
end
toc

%% generate reference solution

t = [0 365];
X0 = model.X0;
par = model.par;
I = model.I;
param = model.param;
odefun = model.odefun;
jacfun = model.jacfun;

tic
[model.t_ref, model.X_ref, log] = simModel(t, X0, par, I, param, odefun, jacfun);
toc

%% generate other model components

% model = model2irenv(model);

%% save model

save("modelfiles/modelCHBR_minimal.mat", 'model')

%% helper functions

% function out = create_and_execute_sym(model, current_rule)
%     disp('calling create_and_execute_sym')
% 
%     X_sym = cell2sym(model.I.nmstate(:));
%     par_sym = cell2sym(model.I.nmpar(:));
% 
%     for i = 1:length(model.I.nmpar)
%         feval(@()assignin('caller',model.I.nmpar{i},par_sym(i)))
%     end
%     for i = 1:length(model.I.nmstate)
%         feval(@()assignin('caller',model.I.nmstate{i},X_sym(i)))
%     end
% 
%     % rulefun_symbolic = evalin('caller', current_rule);
%     if contains(current_rule, "piecewise")
%         % rulefun_symbolic = evalin('caller', current_rule);
%         % ATTENTION: only works if there is just one piecewise in the rule
%         out = eval(current_rule);
%     else
%         rulefun_symbolic = evalin(symengine, current_rule);
%         rulefun_matlabfun = matlabFunction(rulefun_symbolic,'Vars',{X_sym,par_sym});
%         out = rulefun_matlabfun(model.X0, model.par);
%     end
% end

function model = changes_to_model(model)
    if model.name == "Peterson2010 - Integrated calcium homeostasis and bone remodelling"
        model.rule.rules{5} = 'exp((log(J14OCgam)/log(J14OCmax*power(OC0,J14OCgam))/T13-power(OC0,J14OCgam))/J14OCgam)';
        model.rule.rules{33} = 'exp((log(kinOCgam)/log(power(M0,kinOCgam))*EmaxMeffOC/(1-E0Meff)-power(M0,kinOCgam))/kinOCgam)';
    end
end

function out = piecewise(val1, cond1, val2)
    if subs(cond1)
        out = val1;
    else
        out = val2;
    end
end

function jacsym = odesym2jacsym(odesym, X_sym, model)

I = model.I;
jacsym = sym(zeros(I.nstates, I.nstates));

for i = 1:I.nstates
    for j = 1:I.nstates
        jacsym(i, j) = diff(odesym(i), X_sym(j));
    end
end

end
