addpath(genpath('../Core'))
addpath(genpath('./results'))

%% standard parameters
eout = 0.2;
eint = 1;
pint = 100;
timeout = 120;
% crit = 'linear';
crit = 'linear';
% errtype = 'MRSE_conditioned';
errtype = 'MRSE';

%% run information
firstrun = false;
pnegrun = true;
conlawrun = false;
% classifs = ["dyn" "cneg" "pneg" "env" "irenv_geom" "pss"];
log_required = false;

%% model selection
% modelname = 'modelBCSV';
modelname = 'modelEGFR';
% modelname = 'modelGlyc';

% load(['./results/' modelname '_exh_t120_MRSE_0.1_0.5_linear_dyncnegpnegenvirenv_geompss.mat'])
% clear log

load(['../Core/modelfiles/' modelname '_minimal.mat'])
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), model.L);

%% variability
variability = 1;
LHS_EOG = 1;
backwards = 0;
variability_input = 0;
var_obj_prctile = 90;
Npop = 100;
if variability
    rng(1)
    X0 = model.X0;
    % Npop has to be even for variability_input
    ref_indv = zeros(1, length(X0));
    ref_indv(1, :) = X0';
    if LHS_EOG
        sigma = 0.16;
        
        virtual_pop_base = lhsnorm(zeros(1,length(X0)), sigma*eye(length(X0)), Npop);
        [U, ~, ~] = svd(model.obs_gramian);
        virtual_pop = (U*virtual_pop_base')';
        virtual_pop = exp(virtual_pop);

        for i = 1:Npop
            virtual_pop(i, :) = virtual_pop(i, :) .* ref_indv;
        end
    elseif variability_input
        virtual_pop = repmat(X0', [Npop, 1]);
        % virtual_pop(:, model.I.input) = 0:(2*X0(model.I.input)/(Npop-1)):(2*X0(model.I.input));
        virtual_pop(:, model.I.input) = (X0(model.I.input)/(0.5*Npop)):(2*X0(model.I.input)/(Npop)):(2*X0(model.I.input));
    else
        virtual_pop = zeros(Npop, length(X0));
        virtual_pop(1, :) = X0';
        for i = 1:length(X0)
            if X0(i) ~= 0
                m = X0(i);
                v = (0.4*X0(i))^2;
                mu = log(m^2 / sqrt(v + m^2));
                sigma = sqrt(log(v / (m^2) + 1));
        
                virtual_pop(1:end, i) = lognrnd(mu, sigma, [1 Npop]);
            end
        end
        % virtual_pop(2:end, model.I.input) = X0(model.I.input);
    end
    % finalize virtual population
    virtual_pop = [ref_indv; virtual_pop];

    X_ref_var = cell([1 Npop + 1]);
    for npop = 2:(Npop + 1)
        [~, X_ref_var{npop}, ~, ~] = simModel(model.t_ref, virtual_pop(npop, :), model.par, model.I, model.param, model.multiple, model.odefun, model.jacfun);
    end
    X_ref_var{1} = model.X_ref;
else
    virtual_pop = 0;
    X_ref_var = 0;
end

if ~backwards
    redmodel = 0;
end

%% mode
mode = 'from_start';
% mode = 'intermediate';
% mode = 'finished';

%%%% classifs
% "dyn" "pneg" "cneg" "irenv_geom" "irenv_arith" "average" "mode" "constant" "constregr" "ssenv" "env" "pss" 

%% run 1
classifs = ["dyn" "cneg" "pneg" "env" "irenv_geom"];

morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, var_obj_prctile, LHS_EOG, variability_input, virtual_pop, X_ref_var, backwards, redmodel, log_required)

%% run 1
% classifs = ["dyn" "cneg" "pneg" "env"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)

%% run 1
% classifs = ["dyn" "cneg" "pneg" "irenv_geom"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)

%% run 1
% classifs = ["dyn" "cneg" "pneg" "irenv_arith"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)

%% run 1
% classifs = ["dyn" "cneg" "pneg" "average"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)

%% run 1
% classifs = ["dyn" "cneg" "pneg" "mode"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)

%% run 1
% classifs = ["dyn" "cneg" "pneg" "ssenv"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)

% %% run 1
% classifs = ["dyn" "cneg" "pneg" "average" "env"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)
% %% run 1
% classifs = ["dyn" "cneg" "pneg" "average" "pss"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)
% %% run 1
% classifs = ["dyn" "cneg" "pneg" "average" "env" "pss"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)

%% run 2
% classifs = ["dyn" "env"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)

%% run 2
% classifs = ["dyn" "irenv_geom"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)

% %% run 2
% eint = 1;
% classifs = ["dyn" "cneg" "pneg" "env" "irenv_geom" "pss"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)
% 
%% run 2
% eint = 1;
% classifs = ["dyn" "cneg" "pneg" "env" "irenv_geom" "irenv_arith"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)
% 
% %% run 3
% eint = 2;
% classifs = ["dyn" "cneg" "pneg" "env" "irenv_arith" "pss"];
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)
% 
% %% run 4
% eint = Inf;
% classifs = ["dyn" "cneg" "pneg" "env" "irenv_geom" "irenv_arith" "pss"];
% errtype = 'MALE';
% 
% morexh(modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, variability_input, virtual_pop, X_ref_var, log_required)

