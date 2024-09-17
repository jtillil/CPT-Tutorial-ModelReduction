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

function [vector_objvals, log] = objfun_vectorized(model, mat_vector_config, statefromtoreduced, timeout, errtype, variability, virtual_pop, X_ref_var, log)

% LOG
log_required = (nargin > 8);

% number of configs to parallelize
Nconfigs = size(mat_vector_config, 1);

% prepare objects to return and to append to log
vector_objvals = zeros(Nconfigs, 4);
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
jacfun = model.jacfun;

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
    Npop = size(virtual_pop, 1);
    for nconfig = 1:Nconfigs
        for npop = 1:Npop
            parfevalID = (nconfig - 1) * Npop + npop;
            f(parfevalID) = parfeval(@objfun, 5, t_ref, X_ref_var{npop}, virtual_pop(npop, :)', par, I, L, param, multiple, odefun, jacfun, mat_vector_config(nconfig, :), errtype);
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
    % pause for 5 seconds between checks
    pause(5)

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
            fprintf('\n          ');
        end
        fprintf(['; ' char(string(toc(start_of_parfeval)))])
        break;
    end
end

% collect future outputs
if variability
    for nconfig = 1:Nconfigs
        ndyn = -1;
        errout = zeros([1 Npop]);
        errdyn100 = zeros([1 Npop]);
        simtime = zeros([1 Npop]);
        for npop = 1:Npop
            parfevalID = (nconfig - 1) * Npop + npop;
            if ~cancelled(parfevalID)
                [obj, objlog, ~, ~, errentry] = f(parfevalID).fetchOutputs();
                errout(npop) = obj.errout;
                errdyn100(npop) = obj.errdyn100;
                simtime(npop) = obj.simtime;
                ndyn = obj.ndyn;
            else
                errout(npop) = 1e6;
                errdyn100(npop) = 1e6;
                simtime(npop) = 1e6;
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
        if ~isfield(obj, 'errdyn100')
            obj.errdyn100 = 1e6;
        elseif (isempty(obj.errdyn100) || isnan(obj.errdyn100))
            obj.errdyn100 = 1e6;
        end
    
        % save future outputs
        % vector_objvals(nconfig, :) = [obj.ndyn obj.errout obj.errdyn100 seconds(f(nconfig).RunningDuration)];
        vector_objvals(nconfig, :) = [ndyn prctile(errout, 90) prctile(errdyn100, 90) mean(simtime)];
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
