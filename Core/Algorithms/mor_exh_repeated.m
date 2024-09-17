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

function [redmodel, log] = mor_exh_repeated(model, mor_options_firstrun, mor_options)

%% setup

% start counting run time
mor_startTime = tic;

% indexing
I = model.I;

% initiate first config
exhaustive_mor.configs = repmat("dyn", 1, I.nstates);
exhaustive_mor.configs(logical(model.state_unimportant)) = "pneg";
lastconfig = exhaustive_mor.configs(1, :);

% initiate exhaustive_mor
exhaustive_mor.objvals = [I.nstates 0 0 0];
exhaustive_mor.criterion = 0;
exhaustive_mor.statefromtoreduced = {0 "dyn" "dyn" 0};
exhaustive_mor.time = 0;

% initiate log
if mor_options.log
    log{1}.outflags = 'initial iteration';
else
    log = 0;
end

iteration = 1;

%% start main loop

fprintf('Exhaustive MOR - start procedure');
fprintf(['\n' mor_options.saveroot]);

if mor_options.conlawrun && any(contains(mor_options.classifs_to_consider, 'pss'))
    fprintf('\nStart initial conlaw run.')
    if mor_options.log
        [exhaustive_mor, log] = mor_exh_conlawrun(model, exhaustive_mor, iteration, lastconfig, mor_options, log);
    else
        exhaustive_mor = mor_exh_conlawrun(model, exhaustive_mor, iteration, lastconfig, mor_options);
    end
    
    iteration = size(exhaustive_mor.objvals, 1);
    validindices = find(exhaustive_mor.objvals(:, 2) < mor_options_firstrun.err_out & exhaustive_mor.objvals(:, 3) < mor_options_firstrun.err_int);
    lastndyn = exhaustive_mor.objvals(validindices(end), 1);
    lastndynindices = find(exhaustive_mor.objvals(:, 1) == lastndyn);
    [~, best_lastndynidx] = min(exhaustive_mor.criterion(lastndynindices));
    best_idx = lastndynindices(best_lastndynidx);
    lastconfig = exhaustive_mor.configs(best_idx, :);
    fprintf(['\nLast viable model found at iteration ' char(num2str(best_idx)) ' with objective values: ' char(num2str(exhaustive_mor.objvals(best_idx, :)))]);
    fprintf('\nConservation laws initially searched, now continue with main loop.')
elseif mor_options.conlawrun && ~any(contains(mor_options.classifs_to_consider, 'pss'))
    fprintf('\nInitial conlaw run was requested, but there are no reductions to "pss" considered.')
end

if mor_options.firstrun
    fprintf('\nStart pre-run.')
    if mor_options.log
        [exhaustive_mor, log] = mor_exh_main_loop(model, exhaustive_mor, iteration, lastconfig, mor_options_firstrun, log);
    else
        exhaustive_mor = mor_exh_main_loop(model, exhaustive_mor, iteration, lastconfig, mor_options_firstrun);
    end
    
    iteration = size(exhaustive_mor.objvals, 1);
    validindices = find(exhaustive_mor.objvals(:, 2) < mor_options_firstrun.err_out & exhaustive_mor.objvals(:, 3) < mor_options_firstrun.err_int);
    lastndyn = exhaustive_mor.objvals(validindices(end), 1);
    lastndynindices = find(exhaustive_mor.objvals(:, 1) == lastndyn);
    [~, best_lastndynidx] = min(exhaustive_mor.criterion(lastndynindices));
    best_idx = lastndynindices(best_lastndynidx);
    lastconfig = exhaustive_mor.configs(best_idx, :);
    fprintf('\nReductions without "pss" explored, now include "pss" reductions.')
end

fprintf('\nStart main run.')
if mor_options.log
    [exhaustive_mor, log] = mor_exh_main_loop(model, exhaustive_mor, iteration, lastconfig, mor_options, log);
else
    exhaustive_mor = mor_exh_main_loop(model, exhaustive_mor, iteration, lastconfig, mor_options);
end

%% save elapsed time and number of steps

elapsedtime = toc(mor_startTime); fprintf('\n[elapsed time = %.1fs]',elapsedtime);

exhaustive_mor.elapsedtime = elapsedtime;
exhaustive_mor.nsteps = size(exhaustive_mor.objvals, 1);

%% finish reporting over the model
[redmodel, log] = mor_exh_finish_intermediate(model, exhaustive_mor, mor_options, true, log);

end
