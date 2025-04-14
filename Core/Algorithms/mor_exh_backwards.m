function [backwards_mor, log] = mor_exh_backwards(model, redmodel, backwards_mor, iteration, mor_options, log)

% show options
% fprintf('\n')
% disp(mor_options)

% exhaustive_mor

% arguments
relerrTOL_out = mor_options.err_out;
relerrTOL_int = mor_options.err_int;
timeout = mor_options.timeout;
criterion = mor_options.criterion;
classifs_to_consider = mor_options.classifs_to_consider;
saveroot = mor_options.saveroot;
errtype = mor_options.errtype;
variability = mor_options.variability;
virtual_pop = mor_options.virtual_pop;
X_ref_var = mor_options.X_ref_var;

% LOG
log_required = (nargin > 5);

% indexing
I = model.I;

% X0 values
% X0empty = (model.X0 == 0);

% unnecessary states
% state_unimportant = model.state_unimportant;

%% main loop
while true
    % handle next iteration
    iteration = iteration + 1; fprintf('\n     Iteration: %i', iteration);
    % [test_configs, statefromtoreduced] = generate_configs(lastconfig, I, classifs_to_consider, X0empty, state_unimportant, criterion);
    test_configs = redmodel.exhaustive_mor.configs(backwards_mor.morexh_iteration(end), :);
    statefromtoreduced = {1, "dyn", "dyn", 1};

    % LOG: test_configs AND statefromtoreduced
    if log_required
        log{end+1}.test_configs = test_configs;
        log{end}.statefromtoreduced = statefromtoreduced;
    end

    % calculate all objective values
    objfun_startTime = tic;
    if log_required
        [test_objvals, log] = objfun_vectorized(model, test_configs, statefromtoreduced, timeout, errtype, variability, virtual_pop, X_ref_var, log);
    else
        test_objvals = objfun_vectorized(model, test_configs, statefromtoreduced, timeout, errtype, variability, virtual_pop, X_ref_var);
    end
    objfun_time = toc(objfun_startTime);

    % calculate criteria
    % switch criterion
    %     case 'out'
    %         objval_criteria = test_objvals(:,2);
    %     case 'linear'
    %         objval_criteria = test_objvals(:,2)/relerrTOL_out + test_objvals(:,3)/relerrTOL_int;
    %     case 'linear_time'
    %         objval_criteria = test_objvals(:,2)/relerrTOL_out + test_objvals(:,3)/relerrTOL_int + 0.1 * test_objvals(:,4) / test_objvals(1,4);
    %     % case 'remaining'
    %     %     objval_criteria = (relerrTOL_int - exhaustive_mor.objvals(iteration-1, 3))*test_objvals(:,2) + (relerrTOL_out - exhaustive_mor.objvals(iteration-1, 2))*test_objvals(:,3);
    %     case 'quadratic'
    %         objval_criteria = (test_objvals(:,2)/relerrTOL_out).^2 + (test_objvals(:,3)/relerrTOL_int).^2;
    %     case 'max'
    %         objval_criteria = max(test_objvals(:,2)/relerrTOL_out, test_objvals(:,3)/relerrTOL_int);
    % end

    % if criterion includes time, skip first objval
    % if criterion == "linear_time"
    %     objval_criteria(1) = inf;
    % end

    % LOG: raw criterion
    % if log_required
    %     log{end}.criterion_raw = objval_criteria;
    % end

    % read min
    % [best_val, best_indx] = min(objval_criteria);

    % avoid best_indx = 1 if criterion includes time
    % if criterion == "linear_time" && best_indx == 1
    %     best_indx = 2;
    %     best_val = inf;
    % end

    % handle repeated reductions
    % while statefromtoreduced{best_indx, 4} == 1 && best_val >= exhaustive_mor.criterion(iteration-1)
        % test_objvals(best_indx,:) = [I.nstates+1 1e6 1e6 1e6];
        % switch criterion
        %     case 'out'
        %         objval_criteria(best_indx) = test_objvals(best_indx,2);
        %     case 'linear'
        %         objval_criteria(best_indx) = test_objvals(best_indx,2)/relerrTOL_out + test_objvals(best_indx,3)/relerrTOL_int;
        %     case 'linear_time'
        %         objval_criteria(best_indx) = test_objvals(best_indx,2)/relerrTOL_out + test_objvals(best_indx,3)/relerrTOL_int + 0.1 * test_objvals(best_indx,4) / test_objvals(1,4);
        %     case 'remaining'
        %         objval_criteria(best_indx) = (relerrTOL_int - exhaustive_mor.objvals(iteration-1, 3))*test_objvals(best_indx,2) + (relerrTOL_out - exhaustive_mor.objvals(iteration-1, 2))*test_objvals(best_indx,3);
        %     case 'quadratic'
        %         objval_criteria(best_indx) = (test_objvals(best_indx,2)/relerrTOL_out).^2 + (test_objvals(best_indx,3)/relerrTOL_int).^2;
        %     case 'max'
        %         objval_criteria(best_indx) = max(test_objvals(best_indx,2)/relerrTOL_out, test_objvals(best_indx,3)/relerrTOL_int);
        % end
        % [best_val, best_indx] = min(objval_criteria);
    % end

    % LOG: criterion without reduced states
    % if log_required
    %     log{end}.criterion = objval_criteria;
    % end

    % % if best == env && env == pneg == cneg, take cneg
    % if statefromtoreduced{best_indx, 3} == "env" && (objval_criteria(best_indx) == objval_criteria(best_indx+1) && objval_criteria(best_indx) == objval_criteria(best_indx+2))
    %     best_indx = best_indx + 2;
    % end
    % % if best == env && env == pneg < cneg, take pneg
    % if statefromtoreduced{best_indx, 3} == "env" && objval_criteria(best_indx) == objval_criteria(best_indx+1)
    %     best_indx = best_indx + 1;
    % end
    % % if best == pneg && pneg == cneg, take cneg
    % if statefromtoreduced{best_indx, 3} == "pneg" && objval_criteria(best_indx) == objval_criteria(best_indx+1)
    %     best_indx = best_indx + 1;
    % end
    % % TODO: also works for reduced states?

    % show current objvals
    best_indx = 1;
    fprintf(['\n          ' char(num2str(test_objvals(best_indx, 1))) ' ' char(num2str(test_objvals(best_indx, 2))) ' ' char(num2str(test_objvals(best_indx, 3))) ' ' char(num2str(test_objvals(best_indx, 4))) '; ' char(I.nmstate{statefromtoreduced{best_indx, 1}}) ' (' char(num2str(statefromtoreduced{best_indx, 1})) ') from ' char(statefromtoreduced{best_indx, 2}) ' to ' char(statefromtoreduced{best_indx, 3})])
    
    % save to output
    backwards_mor.configs(iteration, :) = test_configs(best_indx, :);
    backwards_mor.objvals(iteration, :) = test_objvals(best_indx, :);
    % backwards_mor.criterion(iteration) = best_val;
    backwards_mor.statefromtoreduced(iteration, 1:4) = statefromtoreduced(best_indx, 1:4);
    backwards_mor.time(iteration, :) = objfun_time;

    % decrease morexh_iteration
    backwards_mor.morexh_iteration(end+1) = backwards_mor.morexh_iteration(end) - 1;

    % stop if no viable config found
    % if test_objvals(best_indx, 2) > relerrTOL_out || test_objvals(best_indx, 3) > relerrTOL_int
    %     if remaining_additional_dyn_reductions == 0
    %         fprintf('\n     finished at iteration %i ', iteration);
    %         fprintf(['after ' num2str(max_additional_dyn_reductions) ' additional dynamic reductions.'])
    %         break;
    %     elseif statefromtoreduced{best_indx,2} == "dyn"
    %         remaining_additional_dyn_reductions = remaining_additional_dyn_reductions - 1;
    %         fprintf([' remaining dynamic reductions: ' num2str(remaining_additional_dyn_reductions)])
    %     end
    % else
    %     % else set remaining iterations to max
    %     remaining_additional_dyn_reductions = max_additional_dyn_reductions;
    % end

    % stop if test_objvals suddenly below error bounds
    if test_objvals(1,2) <= relerrTOL_out && test_objvals(1,3) <= relerrTOL_int
        fprintf('\nFinished after relative errors back below error bounds.');
        break;
    end

    % save at suitable iterations
    % if log_required
    %     save(['results/' saveroot '_intermediate.mat'], 'backwards_mor', 'mor_options', 'log');
    % else
        save(['results/' saveroot '_bintermediate.mat'], 'backwards_mor', 'mor_options');
    % end

    % set next lastconfig
    % lastconfig = test_configs(best_indx, :);
end

end
