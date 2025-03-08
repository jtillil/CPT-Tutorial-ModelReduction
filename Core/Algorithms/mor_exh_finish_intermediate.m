%%% Version: January 24th, 2020
%%%
%%% call by: redmodel  =  model_order_reduction(model,seqofstates,relerrTOL)
%%%
%%% This function returns the classification of the state variables as
%%% environmental, negligible, quasi steady state, mass conserved or
%%% dynamical state variable
%%%
%%% Input:  model               structure specifying the model
%%%         seqofstates         sequence, in which state variables are
%%%                             tested for model order reduction
%%%         relerrTOL           user defined relative error threshold
%%%
%%% Output: redmodel            structure specifying the reduced order model
%%%
%%% Citation:
%%% 
%%% Knoechel, Kloft and Huisinga, "Sensitivity based input-response index to 
%%% analyse and reduce large-scale signalling networks"
%%% PLOS Comp. Biology, 2020 (under review)
%%% 
%%% Authors: Johannes Tillil
%%%

function [redmodel, log] = mor_exh_finish_intermediate(model, exhaustive_mor, mor_options, from_repeated, log)

%% indexing
I = model.I;

%% save final intenal output, log, elapsed time and number of steps

exhaustive_mor.elapsedtime = sum(exhaustive_mor.time);
exhaustive_mor.nsteps = size(exhaustive_mor.objvals, 1);

%% report over reduced model

fprintf('\n\nReturning to last viable model.');

% get reduced model iteration and config
validindices = find(exhaustive_mor.objvals(:, 2) < mor_options.err_out & exhaustive_mor.objvals(:, 3) < mor_options.err_int);
% lastndyn = exhaustive_mor.objvals(validindices(end), 1);
% lastndynindices = find(exhaustive_mor.objvals(:, 1) == lastndyn);
% [~, best_lastndynidx] = min(exhaustive_mor.criterion(lastndynindices));
% best_idx = lastndynindices(best_lastndynidx);
best_idx = validindices(end);
exhaustive_mor.redconfig = exhaustive_mor.configs(best_idx, :);
% exhaustive_mor.redconfig = exhaustive_mor.configs(validindices(end), :);
exhaustive_mor.finaliteration = best_idx;

model.exhaustive_mor = exhaustive_mor;
if mor_options.log
    model.log = log;
end

fprintf(['\nLast viable model found at iteration ' char(num2str(best_idx)) ' with objective values: ' char(num2str(exhaustive_mor.objvals(best_idx, :)))]);

%% perform PNEG testing

if mor_options.pnegrun
    fprintf('\nStart PNEG testing run.')
    % if mor_options.log
    %     [model, pnegconfig, log] = mor_exh_pneg(model, exhaustive_mor.redconfig);
    %     finalconfig = pnegconfig;
    % else
        [model.pnegobj, model.pneg_run, pnegconfig] = mor_exh_pneg_greedy(model, exhaustive_mor.redconfig, mor_options, false);
        finalconfig = pnegconfig;
    % end
else
    finalconfig = exhaustive_mor.redconfig;
end

%% final reporting

% calculate reduced model solution
[model.redobj, model.redlog, model.t_red, model.X_red] = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, finalconfig, mor_options.errtype);
model.Y_red = model.X_red(:,I.output);
model.redobj.redconfig = finalconfig;

% calculate variability reduced model stats
if mor_options.variability
    model.variability.Npop = size(mor_options.virtual_pop);
    model.variability.virtual_pop = mor_options.virtual_pop;
    model.variability.redobj = cell([1 model.variability.Npop]);
    model.variability.errout = zeros([1 model.variability.Npop]);
    model.variability.errdyn100 = zeros([1 model.variability.Npop]);
    for npop = 1:model.variability.Npop
        [model.variability.redobj{npop}, ~, ~, ~] = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, finalconfig, mor_options.errtype);
        model.variability.errout(npop) = model.variability.redobj{npop}.errout;
        model.variability.errdyn100(npop) = model.variability.redobj{npop}.errdyn100;
    end
    model.variability.q90errout = prctile(model.variability.errout, 90);
    model.variability.q90errdyn100 = prctile(model.variability.errdyn100, 90);
end

% finish model
model.mor_options = mor_options;
redmodel = model;
redmodel.I = config2I(model.I, finalconfig, model.L);
rI = redmodel.I;

fprintf('\n\nReduced model consists of');
fprintf('\n %d dynamical state variable(s): ',length(rI.dyn)); 
fprintf('%s, ',I.nmstate{rI.dyn})
fprintf('\n %d environmental state variable(s): ',length(rI.env)); 
fprintf('%s, ',I.nmstate{rI.env})
fprintf('\n %d steady-state environmental state variable(s): ',length(rI.ssenv));
fprintf('%s, ',I.nmstate{rI.ssenv})
fprintf('\n %d obs-mode environmental state variable(s): ',length(rI.mode));
fprintf('%s, ',I.nmstate{rI.mode})
fprintf('\n %d average C-t environmental state variable(s): ',length(rI.average));
fprintf('%s, ',I.nmstate{rI.average})
fprintf('\n %d arithmetic obs-weighted C-t environmental state variable(s): ',length(rI.irenv_arith));
fprintf('%s, ',I.nmstate{rI.irenv_arith})
fprintf('\n %d geometric obs-weighted C-t environmental state variable(s): ',length(rI.irenv_geom));
fprintf('%s, ',I.nmstate{rI.irenv_geom})
fprintf('\n %d partially neglected state variable(s): ',length(rI.pneg)); 
fprintf('%s, ',I.nmstate{rI.pneg})
fprintf('\n %d completely neglected state variable(s): ',length(rI.cneg)); 
fprintf('%s, ',I.nmstate{rI.cneg})
fprintf('\n %d state variable(s) approximated by their quasi-steady state: ',length(rI.pss)); 
fprintf('%s, ',I.nmstate{rI.pss})
fprintf('\n %d state variable(s) eliminated by conservation law(s): ',length(rI.con)); 
fprintf('%s, ',I.nmstate{rI.con})
fprintf('\n')

% save if called outside script
if ~from_repeated
    save(['results/' mor_options.saveroot '.mat'], 'redmodel', 'log', '-v7.3');
end

end
