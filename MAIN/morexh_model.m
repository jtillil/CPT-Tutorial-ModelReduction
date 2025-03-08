%% setup

function morexh_model(model, modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, var_obj_prctile, LHS_EOG, variability_input, virtual_pop, X_ref_var, backwards, redmodel, config, log_required)

% clear; clc; close all;

addpath(genpath('../Core'))
addpath(genpath('./results'))

if isempty(gcp('nocreate'))
    parpool();
end

% mode = 'from_start';
% mode = 'intermediate';
% mode = 'finished';

mor_options.err_out = eout;
mor_options.err_int = eint;
mor_options.prct_int = pint;
mor_options.timeout = timeout; % in seconds, how long to allow a specific redmodel calculation to run. can also be inf for no bound
mor_options.criterion = crit; % one of 'out', 'linear', 'linear_time', 'max', 'remaining', 'quadratic'
mor_options.errtype = errtype; % one of 'MRSE', 'MALE'
mor_options.config = config;

mor_options.variability = variability;
mor_options.var_obj_prctile = var_obj_prctile;
mor_options.backwards = backwards;
mor_options.virtual_pop = virtual_pop;
mor_options.X_ref_var = X_ref_var;

% mor_options.classifs_to_consider = ["dyn" "cneg" "pneg" "env"];
% mor_options.classifs_to_consider = ["dyn" "cneg" "pneg" "irenv_geom" "irenv_arith" "env" "pss"];
% mor_options.classifs_to_consider = ["dyn" "cneg" "pneg" "env" "pss"];
% mor_options.classifs_to_consider = ["dyn" "cneg" "pneg" "env" "irenv_geom" "pss"];
mor_options.classifs_to_consider = classifs;
mor_options.conlawrun = conlawrun; % before firstrun, initially set conlaws
mor_options.firstrun = firstrun; % without qss
mor_options.pnegrun = pnegrun;
mor_options.log = log_required; % if to output logs

mor_options.saveroot = [modelname '_exh_t' num2str(mor_options.timeout) '_' errtype];
if mor_options.variability
    mor_options.saveroot = [mor_options.saveroot '_variability'];
    if variability_input
        mor_options.saveroot = [mor_options.saveroot 'input'];
    end
    if LHS_EOG
        mor_options.saveroot = [mor_options.saveroot 'lhseog'];
    end
    mor_options.saveroot = [mor_options.saveroot char(string(size(virtual_pop, 1) - 1))];
end
if mor_options.backwards
    mor_options.saveroot = [mor_options.saveroot '_backwards'];
end
if mor_options.conlawrun
    mor_options.saveroot = [mor_options.saveroot '_conrun'];
end
if mor_options.firstrun
    mor_options.saveroot = [mor_options.saveroot '_firstrun'];
end
if mor_options.pnegrun
    mor_options.saveroot = [mor_options.saveroot '_pnegrun_greedy'];
end
mor_options.saveroot = [mor_options.saveroot '_' num2str(mor_options.err_out) '_' num2str(mor_options.err_int) '_' mor_options.criterion '_' char(strjoin(mor_options.classifs_to_consider, ''))];

% disp(mor_options)

% If wanted, first run without qss states for efficiency but only until
% half of the error bounds
mor_options_firstrun = mor_options;
mor_options_firstrun.err_out = mor_options.err_out / 2;
mor_options_firstrun.err_int = mor_options.err_int / 2;
mor_options_firstrun.classifs_to_consider = setdiff(mor_options.classifs_to_consider, "pss", 'stable');

%% model reduction
switch mode
    case 'from_start'
        % model = load([modelname '_minimal.mat']).model;
        if mor_options.backwards
            [redmodel, log] = mor_exh_repeated_backwards(model, redmodel, mor_options);
        else
            [redmodel, log] = mor_exh_repeated(model, mor_options_firstrun, mor_options);
        end
    case 'intermediate'
        model = load([modelname '_minimal.mat']).model;
        load(['results/' mor_options.saveroot '_intermediate.mat'])
        [redmodel, log] = mor_exh_repeated_from_iteration(model, exhaustive_mor, log, mor_options_firstrun, mor_options);
    case 'finished'
        model = load([modelname '_minimal.mat']).model;
        finishedmodel = load(['results/' mor_options.saveroot '.mat']).redmodel;
        % exhaustive_mor = load(['results/modelEGFR_red_exh_' num2str(mor_options.err_out) '_' num2str(mor_options.err_int) '_' mor_options.criterion '_gen' num2str(iteration) '.mat']).exhaustive_mor;
        % log = load(['results/modelEGFR_red_exh_' num2str(mor_options.err_out) '_' num2str(mor_options.err_int) '_' mor_options.criterion '_gen' num2str(iteration) '.mat']).log;
        [redmodel, log] = mor_exh_repeated_from_iteration(model, finishedmodel.exhaustive_mor, finishedmodel.log, mor_options_firstrun, mor_options);
end


%% save

save(['results/' mor_options.saveroot '.mat'], 'redmodel', 'log', '-v7.3');

end