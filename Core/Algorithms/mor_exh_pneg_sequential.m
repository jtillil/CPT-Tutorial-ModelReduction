function [pnegobj, pneg_run, pnegconfig, log] = mor_exh_pneg_sequential(model, redconfig, mor_options, log_required)

% show options
% fprintf('\n')
% disp(mor_options)

% arguments
% relerrTOL_out = mor_options.err_out;
% relerrTOL_int = mor_options.err_int;
timeout = mor_options.timeout;
% criterion = mor_options.criterion;
% classifs_to_consider = mor_options.classifs_to_consider;
% saveroot = mor_options.saveroot;
errtype = mor_options.errtype;

% LOG
if log_required
    log = 0;
end

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

% monitor redobj future
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
    nr_changed_states = 0;
    fprintf(['\nPNEG testing loop nr. ' char(num2str(i))])
    for stateID = 1:length(pnegconfig)
        if pnegconfig(stateID) ~= "pneg" && ~state_unimportant(stateID)
            tempconfig = pnegconfig;
            tempconfig(stateID) = "pneg";

            % initiate future to cancel long trials
            start_of_parfeval = tic;
            % tempobj = objfun(t, X_ref, X0, par, I, L, param, multiple, odefun, jacfun, tempconfig, errtype);
            f = parfeval(@objfun, 5, t, X_ref, X0, par, I, L, param, multiple, odefun, jacfun, tempconfig, errtype);

            % monitor and cancel future
            while 1
                % pause for 5 seconds between checks
                pause(1)
            
                % get future state
                state = f.State;
            
                % cancel futures that ran too long
                cancelled = 0;
                if state == "running"
                    % if etime(currenttime, datevec(f.StartDateTime)) > timeout
                    if seconds(f.RunningDuration) > timeout
                        cancel(f)
                        cancelled = 1;
                        fprintf(['\n          Cancelled future for state ' char(num2str(stateID)) ' (' I.nmstate{stateID} ')']);
                    end
                end
            
                % if all finished, stop
                if state == "finished" || cancelled
                    % if cancelled
                    %     fprintf('\n          ');
                    % end
                    fprintf(['; ' char(string(toc(start_of_parfeval)))])
                    break;
                end
            end

            % evaluate future output
            if ~cancelled
                [tempobj, ~, ~, ~, ~] = f.fetchOutputs();
                if (tempobj.errout <= redobj.errout) && (tempobj.errdyn100 <= redobj.errdyn100)
                    fprintf(['\n  State ' char(num2str(stateID)) ' (' char(I.nmstate{stateID}) ') changed from ' char(pnegconfig(stateID)) ' to pneg.'])
                    fprintf(['\n    ' char(num2str(tempobj.errout))])
                    fprintf(['\n    ' char(num2str(tempobj.errdyn100))])
    
                    iteration = size(pneg_run.objvals);
                    iteration = iteration(1) + 1;
                    pneg_run.istatefromto(iteration, 1:4) = {i stateID pnegconfig(stateID) "pneg"};
    
                    pnegconfig(stateID) = "pneg";
                    pnegobj = tempobj;
                    nr_changed_states = nr_changed_states + 1;
    
                    pneg_run.configs(iteration, :) = pnegconfig;
                    pneg_run.objvals(iteration, :) = [pnegobj.ndyn pnegobj.npss pnegobj.errout pnegobj.errdyn100];
                end
            end

            % remove future
            clear f;
        end
    end
    if nr_changed_states == 0
        break
    end
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
