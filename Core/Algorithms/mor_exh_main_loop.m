function [exhaustive_mor, log] = mor_exh_main_loop(model, exhaustive_mor, iteration, lastconfig, mor_options, log)

% show options
% fprintf('\n')
% disp(mor_options)

% arguments
relerrTOL_out = mor_options.err_out;
relerrTOL_int = mor_options.err_int;
timeout = mor_options.timeout;
criterion = mor_options.criterion;
classifs_to_consider = mor_options.classifs_to_consider;
saveroot = mor_options.saveroot;
errtype = mor_options.errtype;
variability = mor_options.variability;
var_obj_prctile = mor_options.var_obj_prctile;
virtual_pop_X0 = mor_options.virtual_pop_X0;
virtual_pop_par = mor_options.virtual_pop_par;
X_ref_var = mor_options.X_ref_var;

% LOG
log_required = (nargin > 5);

% indexing
I = model.I;

% X0 values
X0empty = (model.X0 == 0);

% unnecessary states
state_unimportant = model.state_unimportant;
state_constant = zeros([model.I.nstates, 1]);

%% main loop
while true
    % handle next iteration
    iteration = iteration + 1; fprintf('\n     Iteration: %i', iteration);
    [test_configs, statefromtoreduced] = generate_configs(lastconfig, I, classifs_to_consider, X0empty, state_unimportant, state_constant, criterion);

    % LOG: test_configs AND statefromtoreduced
    if log_required
        log{end+1}.test_configs = test_configs;
        log{end}.statefromtoreduced = statefromtoreduced;
    end

    % calculate all objective values
    objfun_startTime = tic;
    if log_required
        [test_objvals, log] = objfun_vectorized(model, test_configs, statefromtoreduced, timeout, errtype, variability, var_obj_prctile, virtual_pop_X0, virtual_pop_par, X_ref_var, log);
    else
        test_objvals = objfun_vectorized(model, test_configs, statefromtoreduced, timeout, errtype, variability, var_obj_prctile, virtual_pop_X0, virtual_pop_par, X_ref_var);
    end
    objfun_time = toc(objfun_startTime);

    % calculate criteria
    switch criterion
        case 'out'
            objval_criteria = test_objvals(:,2);
        case 'linear'
            objval_criteria = test_objvals(:,2)/relerrTOL_out + test_objvals(:,3)/relerrTOL_int;
        case 'linear_time'
            objval_criteria = test_objvals(:,2)/relerrTOL_out + test_objvals(:,3)/relerrTOL_int + 0.1 * test_objvals(:,4) / test_objvals(1,4);
        case 'remaining'
            objval_criteria = (relerrTOL_int - exhaustive_mor.objvals(iteration-1, 3))*test_objvals(:,2) + (relerrTOL_out - exhaustive_mor.objvals(iteration-1, 2))*test_objvals(:,3);
        case 'quadratic'
            objval_criteria = (test_objvals(:,2)/relerrTOL_out).^2 + (test_objvals(:,3)/relerrTOL_int).^2;
        case 'max'
            objval_criteria = max(test_objvals(:,2)/relerrTOL_out, test_objvals(:,3)/relerrTOL_int);
    end

    % if criterion includes time, skip first objval
    if criterion == "linear_time"
        objval_criteria(1) = inf;
    end

    % LOG: raw criterion
    if log_required
        log{end}.criterion_raw = objval_criteria;
    end

    % read min
    [best_val, best_indx] = min(objval_criteria);

    % avoid best_indx = 1 if criterion includes time
    if criterion == "linear_time" && best_indx == 1
        best_indx = 2;
        best_val = inf;
    end

    % handle repeated reductions
    while statefromtoreduced{best_indx, 4} == 1 && best_val >= exhaustive_mor.criterion(iteration-1)
        if mor_options.variability
            test_objvals(best_indx,:) = [I.nstates+1 1e6 1e6 1e6 1e6 1e6];
        else
            test_objvals(best_indx,:) = [I.nstates+1 1e6 1e6 1e6];
        end
        switch criterion
            case 'out'
                objval_criteria(best_indx) = test_objvals(best_indx,2);
            case 'linear'
                objval_criteria(best_indx) = test_objvals(best_indx,2)/relerrTOL_out + test_objvals(best_indx,3)/relerrTOL_int;
            case 'linear_time'
                objval_criteria(best_indx) = test_objvals(best_indx,2)/relerrTOL_out + test_objvals(best_indx,3)/relerrTOL_int + 0.1 * test_objvals(best_indx,4) / test_objvals(1,4);
            case 'remaining'
                objval_criteria(best_indx) = (relerrTOL_int - exhaustive_mor.objvals(iteration-1, 3))*test_objvals(best_indx,2) + (relerrTOL_out - exhaustive_mor.objvals(iteration-1, 2))*test_objvals(best_indx,3);
            case 'quadratic'
                objval_criteria(best_indx) = (test_objvals(best_indx,2)/relerrTOL_out).^2 + (test_objvals(best_indx,3)/relerrTOL_int).^2;
            case 'max'
                objval_criteria(best_indx) = max(test_objvals(best_indx,2)/relerrTOL_out, test_objvals(best_indx,3)/relerrTOL_int);
        end
        [best_val, best_indx] = min(objval_criteria);
    end

    % LOG: criterion without reduced states
    if log_required
        log{end}.criterion = objval_criteria;
    end

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
    fprintf(['\n          ' char(num2str(test_objvals(best_indx, 1))) ' ' char(num2str(test_objvals(best_indx, 2))) ' ' char(num2str(test_objvals(best_indx, 3))) ' ' char(num2str(test_objvals(best_indx, 4))) '; ' char(I.nmstate{statefromtoreduced{best_indx, 1}}) ' (' char(num2str(statefromtoreduced{best_indx, 1})) ') from ' char(statefromtoreduced{best_indx, 2}) ' to ' char(statefromtoreduced{best_indx, 3})])
    
    % save to output
    exhaustive_mor.configs(iteration, :) = test_configs(best_indx, :);
    exhaustive_mor.objvals(iteration, :) = test_objvals(best_indx, :);
    exhaustive_mor.criterion(iteration) = best_val;
    exhaustive_mor.statefromtoreduced(iteration, 1:4) = statefromtoreduced(best_indx, 1:4);
    exhaustive_mor.time(iteration, :) = objfun_time;

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

    % stop at zero dynamic variables or when all out errors are at least 1e6
    if test_objvals(best_indx, 1) == 0
        fprintf('\nFinished after all dynamic variables were reduced.');
        break;
    elseif test_objvals(best_indx, 2) >= 1e6
        fprintf('\nFinished because no possible reduction left.');
        break;
    end

    % save at suitable iterations
    if log_required
        save(['results/' saveroot '_intermediate.mat'], 'exhaustive_mor', 'mor_options', 'log');
    else
        save(['results/' saveroot '_intermediate.mat'], 'exhaustive_mor', 'mor_options');
    end

    % set next lastconfig
    lastconfig = test_configs(best_indx, :);
end

end

%%%% Helper functions %%%%

function [configs, statefromtoreduced] = generate_configs(lastconfig, I, classifs_to_consider, X0empty, state_unimportant, state_constant, criterion)

if length(classifs_to_consider) < 2
    error('Not enough possible state classifications provided - more than one needed.')
end

% initialize configs, statefromtoreduced (4 cols)
nconfigs = 0;
if criterion == "linear_time"
    nconfigs = nconfigs + 1;
end
for stateID = 1:I.nstates  %([1:51 53 55:78 80:112])
    if state_unimportant(stateID)
        nconfigs = nconfigs + 1;
    elseif state_constant(stateID)
        nconfigs = nconfigs + 0;
    elseif lastconfig(stateID) == "dyn"
        if X0empty(stateID)
            nconfigs = nconfigs + length(classifs_to_consider) - 2;
        else
            nconfigs = nconfigs + length(classifs_to_consider) - 1;
        end
    else
        if X0empty(stateID)
            nconfigs = nconfigs + length(classifs_to_consider) - 3;
        else
            nconfigs = nconfigs + length(classifs_to_consider) - 2;
        end
    end
end
configs = strings(nconfigs, I.nstates);
statefromtoreduced = cell(nconfigs, 4);

% fill configs and statefromtoreduced
if criterion == "linear_time"
    configs(1, :) = "dyn";
    statefromtoreduced{1, 1} = 0;
    statefromtoreduced{1, 2} = "dyn";
    statefromtoreduced{1, 3} = "dyn";
    statefromtoreduced{1, 4} = 0;
    curr_config = 2;
else
    curr_config = 1;
end
for stateID = 1:I.nstates
    if state_unimportant(stateID)
        if lastconfig(stateID) == "pneg"
            configs(curr_config,:) = lastconfig;
            statefromtoreduced{curr_config, 1} = stateID;
            statefromtoreduced{curr_config, 2} = lastconfig(stateID);
            configs(curr_config, stateID) = "cneg";
            statefromtoreduced{curr_config, 3} = "cneg";
            statefromtoreduced{curr_config, 4} = 1;
            curr_config = curr_config + 1;
        else
            configs(curr_config,:) = lastconfig;
            statefromtoreduced{curr_config, 1} = stateID;
            statefromtoreduced{curr_config, 2} = lastconfig(stateID);
            configs(curr_config, stateID) = "pneg";
            statefromtoreduced{curr_config, 3} = "pneg";
            statefromtoreduced{curr_config, 4} = 1;
            curr_config = curr_config + 1;
        end
    elseif state_constant(stateID)
        % ...
    elseif lastconfig(stateID) == "dyn"
        for newstateconfig = setdiff(classifs_to_consider, "dyn", 'stable')
            if newstateconfig == "env" && X0empty(stateID)
                continue
            end
            configs(curr_config,:) = lastconfig;
            statefromtoreduced{curr_config, 1} = stateID;
            statefromtoreduced{curr_config, 2} = lastconfig(stateID);
            configs(curr_config, stateID) = newstateconfig;
            statefromtoreduced{curr_config, 3} = newstateconfig;
            statefromtoreduced{curr_config, 4} = 0;
            curr_config = curr_config + 1;
        end
    else
        for newstateconfig = setdiff(classifs_to_consider, ["dyn" lastconfig(stateID)], 'stable')
            if newstateconfig == "env" && X0empty(stateID)
                continue
            end
            configs(curr_config,:) = lastconfig;
            statefromtoreduced{curr_config, 1} = stateID;
            statefromtoreduced{curr_config, 2} = lastconfig(stateID);
            configs(curr_config, stateID) = newstateconfig;
            statefromtoreduced{curr_config, 3} = newstateconfig;
            statefromtoreduced{curr_config, 4} = 1;
            curr_config = curr_config + 1;
        end
    end
end

end