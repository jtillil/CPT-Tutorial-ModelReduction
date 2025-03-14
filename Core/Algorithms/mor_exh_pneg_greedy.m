function [pnegobj, pneg_run, pnegconfig, log] = mor_exh_pneg_greedy(model, redconfig, mor_options, log_required)

% show options
% fprintf('\n')
% disp(mor_options)

% arguments
relerrTOL_out = mor_options.err_out;
relerrTOL_int = mor_options.err_int;
timeout = mor_options.timeout;
criterion = mor_options.criterion;
% classifs_to_consider = mor_options.classifs_to_consider;
% saveroot = mor_options.saveroot;
errtype = mor_options.errtype;

% LOG
% if log_required
    log = 0;
% else
%     log = 0;
% end

% indexing
I = model.I;

% unnecessary states
state_unimportant = model.state_unimportant;

% simulation arguments
t = model.t_ref;
X_ref = model.X_ref;
X0 = model.X0;
par = model.par;
I = model.I;
L = model.L;
param = model.param;
multiple = model.multiple;
odefun = model.odefun;
jacfun = model.jacfun;

% calculate redobj
% [redobj, ~, ~, ~, ~] = objfun(t, X_ref, X0, par, I, L, param, multiple, odefun, jacfun, redconfig, errtype);
f = parfeval(@objfun, 5, t, X_ref, X0, par, I, L, param, multiple, odefun, jacfun, redconfig, errtype);

% monitor redobj futures
while 1
    % pause for 1 second between checks
    pause(1)

    % get future state
    state = f.State;

    % if finished, stop
    if state == "finished"
        break;
    end
end

% evaluate future output
[redobj, ~, ~, ~, ~] = f.fetchOutputs();
clear f;
pnegobj = redobj;
fprintf('\n  Full model objective values:')
fprintf(['\n    ' char(num2str(redobj.errout))])
fprintf(['\n    ' char(num2str(redobj.errdyn100))])

% initiate first pneg_run save
pneg_run.istatefromto = {0 0 "dyn" "dyn"};
pneg_run.configs = redconfig;
pneg_run.objvals = [redobj.ndyn redobj.npss redobj.errout redobj.errdyn100];

% start pneg testing loop
pnegconfig = redconfig;
i = 0;
while 1
    i = i + 1;
    % nr_changed_states = 0;
    fprintf(['\nPNEG iteration nr. ' char(num2str(i))])
    for stateID = 1:length(pnegconfig)
        if pnegconfig(stateID) ~= "pneg"% && ~state_unimportant(stateID)
            tempconfig = pnegconfig;
            tempconfig(stateID) = "pneg";

            % initiate future to cancel long trials
            start_of_parfeval = tic;
            % tempobj = objfun(t, X_ref, X0, par, I, L, param, multiple, odefun, jacfun, tempconfig, errtype);
            f(stateID) = parfeval(@objfun, 5, t, X_ref, X0, par, I, L, param, multiple, odefun, jacfun, tempconfig, errtype);
        else
            f(stateID) = parfeval(@() struct('errout',1e6,'errdyn100',1e6), 1);
        end
    end
    
    % recursively check futures and stop those that ran too long
    cancelled = zeros(1, length(pnegconfig));
    while 1
        % pause for 5 seconds between checks
        pause(1)
    
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

    % evaluate future output
    errout = zeros(1, length(states));
    errdyn100 = zeros(1, length(states));
    for futureID = 1:length(states)
        if ~cancelled(futureID)
            % [tempobj, ~, ~, ~, ~] = f.fetchOutputs();
            tempobj = f(futureID).fetchOutputs();
            errout(futureID) = tempobj.errout;
            errdyn100(futureID) = tempobj.errdyn100;
        else
            errout(futureID) = 1e6;
            errdyn100(futureID) = 1e6;
        end
    end
    switch criterion
        case 'out'
            objval_criteria = errout;
        case 'linear'
            objval_criteria = errout/relerrTOL_out + errdyn100/relerrTOL_int;
        case 'linear_time'
            error("criterion 'linear_time' not implemented")
        case 'remaining'
            error("criterion 'remaining' not implemented")
            % objval_criteria = (relerrTOL_int - exhaustive_mor.objvals(iteration-1, 3))*test_objvals(:,2) + (relerrTOL_out - exhaustive_mor.objvals(iteration-1, 2))*test_objvals(:,3);
        case 'quadratic'
            objval_criteria = (errout/relerrTOL_out).^2 + (errdyn100/relerrTOL_int).^2;
        case 'max'
            objval_criteria = max(errout/relerrTOL_out, errdyn100/relerrTOL_int);
    end

    % find minimum
    [~, idx_min_objval] = min(objval_criteria);

    if (errout(idx_min_objval) <= relerrTOL_out) && (errdyn100(idx_min_objval) <= relerrTOL_int)
        fprintf(['\n  State ' char(num2str(idx_min_objval)) ' (' char(I.nmstate{idx_min_objval}) ') changed from ' char(pnegconfig(idx_min_objval)) ' to pneg.'])
        fprintf(['\n    ' char(num2str(errout(idx_min_objval)))])
        fprintf(['\n    ' char(num2str(errdyn100(idx_min_objval)))])

        % iteration = size(pneg_run.objvals);
        % iteration = iteration(1) + 1;
        pneg_run.istatefromto(i, 1:4) = {i idx_min_objval pnegconfig(idx_min_objval) "pneg"};

        pnegconfig(idx_min_objval) = "pneg";
        pnegobj = objfun(t, X_ref, X0, par, I, L, param, multiple, odefun, jacfun, pnegconfig, errtype);
        % nr_changed_states = nr_changed_states + 1;

        pneg_run.configs(i, :) = pnegconfig;
        pneg_run.objvals(i, :) = [pnegobj.ndyn pnegobj.npss pnegobj.errout pnegobj.errdyn100];
    else
        break
    end

    % remove futures
    clear f;

    % if nr_changed_states == 0
    %     break
    % end
end

% rI = config2I(I, pnegconfig, L);
% 
% fprintf('\n\nPNEG run reduced model consists of');
% fprintf('\n %d dynamical state variable(s): ',length(rI.dyn)); 
% fprintf('%s, ',I.nmstate{rI.dyn})
% fprintf('\n %d environmental state variable(s): ',length(rI.env)); 
% fprintf('%s, ',I.nmstate{rI.env})
% fprintf('\n %d arithmetic ir-weighted environmental state variable(s): ',length(rI.irenv_arith));
% fprintf('%s, ',I.nmstate{rI.irenv_arith})
% fprintf('\n %d geometric ir-weighted environmental state variable(s): ',length(rI.irenv_geom));
% fprintf('%s, ',I.nmstate{rI.irenv_geom})
% fprintf('\n %d partially neglected state variable(s): ',length(rI.pneg)); 
% fprintf('%s, ',I.nmstate{rI.pneg})
% fprintf('\n %d completely neglected state variable(s): ',length(rI.cneg)); 
% fprintf('%s, ',I.nmstate{rI.cneg})
% fprintf('\n %d state variable(s) approximated by their quasi-steady state: ',length(rI.pss)); 
% fprintf('%s, ',I.nmstate{rI.pss})
% fprintf('\n %d state variable(s) eliminated by conservation law(s): ',length(rI.con)); 
% fprintf('%s, ',I.nmstate{rI.con})
% fprintf('\n')

end
