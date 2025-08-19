%%% Version: October 5th, 2023
%%%
%%% call by: vector_objvals  =  redmodel_objfun_vectorized(model,mat_vector_config)
%%%
%%% This function computes errors from any reduced model. Useful for higher
%%% level algorithms. Vectorized function.
%%%
%%% Input:  model               struct specifying the model
%%%         mat_vector_config   matrix of row vectors specifying reduced model configs
%%%
%%% Output: err                 struct containing output and internal errs
%%%
%%% Citation:
%%% 
%%% ---
%%% 
%%% Authors: Johannes Tillil
%%%

function [vector_objvals, log] = objfun_vectorized(model, mat_vector_config, stateparfromtoreduced, timeout, errtype, variability, var_obj_prctile, virtual_pop_X0, virtual_pop_par, X_ref_var, log)

% LOG
log_required = (nargin > 10);

% number of configs to parallelize
Nconfigs = size(mat_vector_config, 1);

% prepare objects to return and to append to log
if variability
    vector_objvals = zeros(Nconfigs, 6);
else
    vector_objvals = zeros(Nconfigs, 4);
end
if log_required
    err = zeros(Nconfigs, model.I.nstates);
    outflags = strings([Nconfigs, 1]);
end

% prepare arguments for parallelized objfun calls
t_ref = model.t_ref;
X_ref = model.X_ref;
X0 = model.X0;
par = model.par;
param = model.param;
multiple = model.multiple;
I = model.I;
L = model.L;
odefun = model.odefun;
if isfield(model, "jacfun")
    jacfun = model.jacfun;
else
    jacfun = [];
end

% parfor nconfig = 1:Nconfigs
%     [obj, objlog, ~, ~, errentry] = objfun(t_ref, X0, par, model, mat_vector_config(nconfig, :), timeout);
%     vector_objvals(nconfig, :) = [obj.ndyn obj.errout obj.errdyn100];
%     outflags(nconfig) = objlog;
%     err(nconfig, :) = errentry;
% end

% start timing
start_of_parfeval = tic;

% create parfeval futures
if variability
    Npop = size(virtual_pop_X0, 1);
    for nconfig = 1:Nconfigs
        for npop = 1:Npop
            parfevalID = (nconfig - 1) * Npop + npop;
            f(parfevalID) = parfeval(@objfun, 5, t_ref, X_ref_var{npop}, virtual_pop_X0(npop, :)', virtual_pop_par(npop, :)', I, L, param, multiple, odefun, jacfun, mat_vector_config(nconfig, :), errtype);
            % f(parfevalID) = parfeval(@objfun, 5, t_ref, X_ref_var{npop}, virtual_pop_X0(npop, :)', par, I, L, param, multiple, odefun, jacfun, mat_vector_config(nconfig, :), errtype);
        end
    end
    cancelled = zeros(1, Nconfigs * Npop);
else
    for nconfig = 1:Nconfigs
        f(nconfig) = parfeval(@objfun, 5, t_ref, X_ref, X0, par, I, L, param, multiple, odefun, jacfun, mat_vector_config(nconfig, :), errtype);
    end
    cancelled = zeros(1, Nconfigs);
end

% recursively check futures and stop those that ran too long
while 1
    % pause for 0.1 seconds between checks
    pause(0.1)

    % get future states
    states = {f.State};

    % cancel futures that ran too long
    cancelled_in_current_iteration = 0;
    % currenttime = datevec(datetime('now'));
    for futureID = 1:length(states)
        if states{futureID} == "running"
            % if etime(currenttime, datevec(f(futureID).StartDateTime)) > timeout
            if seconds(f(futureID).RunningDuration) > timeout
                cancel(f(futureID))
                cancelled(futureID) = 1;
                cancelled_in_current_iteration = 1;
                % fprintf(['\n          Cancelled future for: ' char(model.I.nmstate{statefromtoreduced{futureID, 1}}) ' (' char(num2str(statefromtoreduced{futureID, 1})) ') from ' char(statefromtoreduced{futureID, 2}) ' to ' char(statefromtoreduced{futureID, 3})]);
            end
        end
    end

    % if all finished, stop
    if all(strcmp(states, 'finished'))
        if cancelled_in_current_iteration
            fprintf('\n    ');
        end
        fprintf(['; time: ' char(string(toc(start_of_parfeval)))])
        break;
    end
end

% collect future outputs
if variability
    for nconfig = 1:Nconfigs
        ndyn = -1;
        % obtain virtual pop results
        errout = zeros([1 Npop-1]);
        errdyn100 = zeros([1 Npop-1]);
        simtime = zeros([1 Npop-1]);
        for npop = 2:Npop
            parfevalID = (nconfig - 1) * Npop + npop;
            if ~cancelled(parfevalID)
                [obj, ~, ~, ~, ~] = f(parfevalID).fetchOutputs();
                errout(npop) = obj.errout;
                if ~isfield(obj, 'errdyn100')
                    obj.errdyn100 = 9e6;
                elseif (isempty(obj.errdyn100) || isnan(obj.errdyn100))
                    obj.errdyn100 = 9e6;
                end
                errdyn100(npop) = obj.errdyn100;
                simtime(npop) = obj.simtime;
                ndyn = obj.ndyn;
            else
                errout(npop) = 2e6;
                errdyn100(npop) = 2e6;
                simtime(npop) = 2e6;
            end
        end
        % obtain reference parametrization results
        errout_ref = 0;
        errdyn100_ref = 0;
        for npop = 1
            parfevalID = (nconfig - 1) * Npop + npop;
            if ~cancelled(parfevalID)
                [obj, ~, ~, ~, ~] = f(parfevalID).fetchOutputs();
                errout_ref = obj.errout;
                if ~isfield(obj, 'errdyn100')
                    obj.errdyn100 = 1e6;
                elseif (isempty(obj.errdyn100) || isnan(obj.errdyn100))
                    obj.errdyn100 = 1e6;
                end
                errdyn100_ref = obj.errdyn100;
            else
                errout_ref = 2e6;
                errdyn100_ref = 2e6;
            end
        end
        if ndyn == -1
            ndyn = model.I.nstates + 1;
        end
    
        % fill bad outputs
        % if ~isfield(obj, 'ndyn')
        %     obj.ndyn = model.I.nstates + 1;
        % end
        % if ~isfield(obj, 'errout')
        %     obj.errout = 1e6;
        % end
    
        % save future outputs
        % vector_objvals(nconfig, :) = [obj.ndyn obj.errout obj.errdyn100 seconds(f(nconfig).RunningDuration)];
        if any(isnan(errout)) || any(isnan(errdyn100))
            vector_objvals(nconfig, :) = [model.I.nstates+1, 3e6, 3e6, 3e6, 3e6, 3e6];
        elseif isnan(errout_ref) || isnan(errdyn100_ref)
            vector_objvals(nconfig, :) = [model.I.nstates+1, 4e6, 4e6, 4e6, 4e6, 4e6];
        else
            vector_objvals(nconfig, :) = [ndyn prctile(errout, var_obj_prctile) prctile(errdyn100, var_obj_prctile) mean(simtime) errout_ref errdyn100_ref];
        end
        if log_required
            outflags = '';
            err = 0;
        end
    end
else
    for nconfig = 1:Nconfigs
        if ~cancelled(nconfig)
            [obj, objlog, ~, ~, errentry] = f(nconfig).fetchOutputs();
        else
            obj.ndyn = model.I.nstates + 1;
            obj.errout = 1e6;
            obj.errdyn100 = 1e6;
            obj.simtime = 1e6;
            if log_required
                objlog = 'future cancelled due to run time limit';
                errentry = zeros(1, model.I.nstates) + 1e6;
            end
        end
    
        % fill bad outputs
        % if ~isfield(obj, 'ndyn')
        %     obj.ndyn = model.I.nstates + 1;
        % end
        % if ~isfield(obj, 'errout')
        %     obj.errout = 1e6;
        % end
        if ~isfield(obj, 'errdyn100')
            obj.errdyn100 = 1e6;
        elseif (isempty(obj.errdyn100) || isnan(obj.errdyn100))
            obj.errdyn100 = 1e6;
        end
    
        % save future outputs
        % vector_objvals(nconfig, :) = [obj.ndyn obj.errout obj.errdyn100 seconds(f(nconfig).RunningDuration)];
        vector_objvals(nconfig, :) = [obj.ndyn obj.errout obj.errdyn100 obj.simtime];
        if log_required
            outflags(nconfig) = objlog;
            err(nconfig, :) = errentry;
        end
    end
end

% remove futures
clear f;

% append to log
if log_required
    log{end}.rel_err = err;
    log{end}.outflags = outflags;
else
    log = 0;
end

end
