%%% Version: October 28th, 2023
%%%
%%% call by: [tout, Xout, log] = simModel(t, X0, par, model)
%%%
%%% This function calculates one solution (and returns a log of the
%%% calculation) of a model, specified by 
%%%
%%% Input:  model config
%%%
%%% Output: [tout, Xout, log]   output
%%%
%%% Citation:
%%% 
%%% ---
%%% 
%%% Authors: Johannes Tillil
%%%

% function [tout, Xout, log] = simModel(t, X0, par, model, timeout)
function [tout, Xout, log, simtime] = simModel(t, X0, par, I, param, multiple, odefun, jacfun)

% init log
log = 'Log: ';

%% (1) pre-processing for index analysis

% set negligible state variables to zero (to realise that they are not
% part of the model
X0([I.pneg I.cneg]) = 0;

% set irenv state variables to ir-weighted initial state
% if (isfield(I, 'irenv'))
%     X0(I.irenv) = I.Irenv(I.irenv);
% end
if (isfield(I, 'irenv_arith'))
    if (~isempty(I.irenv_arith))
        X0(I.irenv_arith) = I.Irenv_arith(I.irenv_arith);
    end
end
if (isfield(I, 'irenv_geom'))
    if (~isempty(I.irenv_geom))
        X0(I.irenv_geom) = I.Irenv_geom(I.irenv_geom);
    end
end
if (isfield(I, 'average'))
    if (~isempty(I.average))
        X0(I.average) = I.Average(I.average);
    end
end
if (isfield(I, 'mode'))
    if (~isempty(I.mode))
        X0(I.mode) = I.Mode(I.mode);
    end
end
if (isfield(I, 'ssenv'))
    if (~isempty(I.ssenv))
        X0(I.ssenv) = I.Ssenv(I.ssenv);
    end
end
if (isfield(I, 'constant'))
    if (~isempty(I.constant))
        X0(I.constant) = I.Constant(I.constant);
    end
end
if (isfield(I, 'constregr'))
    if (~isempty(I.constregr))
        X0(I.constregr) = I.Constregrfcn(I.constregr, X0, par);
    end
end

% set all reaction rate constants related to the species in I.cneg (both,
% forward and backward reactions) to zero 
if ~isempty(I.cneg)
    I_par = param.states2Ipar(I.cneg);
    par([I_par{:}]) = 0;
end

% add jacobian, if provided
if ~isempty(jacfun)
    options.Jacobian = @(t,X) jacfunModel(X,par,I,jacfun);
end

% check consistency of states that are part of I.replaceODE
%
% if a state (say the k-th one) is part of I.replaceODE, then its
% initial value X0(k) and its ODE dX(k) are replaced by the sum of
% inital values and ODEs of the states that are specified in
% I.replaceODEby
for p = 1:length(I.replaceODE)

    k = I.replaceODE(p);        % index of ODE to be replaced 
    states = I.replaceODEby{p}; % index of ODEs whose sum is used for replacement  

    % consistency check
    if ~ismember(k,states)
        fprintf('\n\n --> The ODE of state %s does not belong to the ODEs that are used to replace it :-( --- please fix \n\n',I.nmstate{k});
        error('--> Ending here')
    end

    % check, whether any of the other replaceODE states is part of one of the
    % other used replacedbyODEs (this is not allowed)
    % report the first one
    if any( ismember( setdiff(I.replaceODE,k), setdiff(states,k) ) )
        commonstates  = intersect(setdiff(I.replaceODE,k),states);
        fprintf('\n\n --> The state %s in replaceODE does belong to multiple sets of states in replaceODEby :-( --- please fix \n\n',I.nmstate{commonstates(1)});
        error('--> Ending here')      
    end

    % set initial value of replaceODE state
    X0(k) = sum(X0(states));

end   

%% (2) specify options

% specify options
if isempty([I.pss])
    % default, unless changed below
    options.NonNegative = 1:I.nstates;
    
else
    % DAE solver not able to deal with NonNegative option
    options.NonNegative = [];
    
    % so far unresolved problem, if Jacobian is provided for DAE case with
    % simultaneous conservation laws; not providing Jacobian seems to be 
    % an interim solution
    if ~isempty(I.replaceODE) && ~isempty(jacfun)
        options.Jacobian = [];
        %fprintf('\n\n NOTE: for simulatenous pss and con states, analytically provided Jacobian is not used! \n\n ');
    end
    
    % define mass matrix for DAE solver
    M = eye(I.nstates);
    M(I.pss,I.pss) = 0;
    options.Mass = M;
    options.MassSingular = 'yes';
end

%% (3) simulate model ODEs

simtime_start = tic;

if ~multiple.multiple
    % without multiple dosing
    [tout, Xout, log] = odesolver_loop(t, X0, par, I, odefun, options, log);

else
    % handling for multiple dosing
    input_events = multiple.input_events;
    input_doses = multiple.input_doses;
    u_ref = multiple.u_ref;
    tout = [];
    Xout = [];
    % te_out = [];
    tspan_total = t(:);
    X0_ref = X0(:) + u_ref;
    
    % [tout, Xout, log] = odesolver_loop(tspan_local, X0, par, I, odefun, options, log);
    if tspan_total(1) < input_events(1)
        if length(tspan_total)==2
            tspan_local = [tspan_total(1) input_events(1)];
        else
            tspan_local = tspan_total(tspan_total<input_events(1));
        end
        [tout, Xout, log] = odesolver_loop(tspan_local, X0, par, I, odefun, options, log);
        X0_ref = Xout(end,:)' + u_ref;
    end
    for i = 2:length(input_events)
        if input_events(i) <= tspan_total(end) && tspan_total(1) < input_events(i)
            %%% in case of specification of specific timepoints to be solved need
            %%% to account for this
            if length(tspan_total) <= 2
                tspan_local = [input_events(i-1) input_events(i)];
            else
                tspan_local = unique([input_events(i-1);tspan_total(tspan_total>=input_events(i-1) & tspan_total<=input_events(i));input_events(i)]);
            end
            [tout_local, Xout_local, log] = odesolver_loop(tspan_local, X0_ref, par, I, odefun, options, log);
            %%% to avoid double timepoints in the time vector in the multiple
            %%% dosing case, the dosing timepoint is kept and the other discarded
            %%% such that in the state vector only the dosing state is included
            if i < length(input_events)
                tout = [tout; tout_local(1:end-1)];
                Xout = [Xout; Xout_local(1:end-1,:)];
            else
                if tspan_total(end) > input_events(end)
                    tout = [tout; tout_local(1:end-1)];
                    Xout = [Xout; Xout_local(1:end-1,:)];
                else
                    tout = [tout; tout_local];
                    Xout = [Xout; Xout_local];
                end
            end
            if any(input_doses)
                u_ref(I.input) = input_doses(i);
            end
            X0_ref = Xout_local(end,:)' + u_ref;
        elseif tspan_total(1) >= input_events(i)
            continue;
        else
            break;
        end
    end

    if tspan_total(end)>input_events(end)
        if length(tspan_total)==2
            tspan_local = [input_events(end) tspan_total(end)];
        else
            tspan_local = tspan_total(tspan_total >= input_events(end));
        end
        [tout_local, Xout_local, log] = odesolver_loop(tspan_local, X0_ref, par, I, odefun, options, log);
        tout = [tout; tout_local];
        Xout = [Xout; Xout_local];
    end

    Xout = Xout(ismember(tout, tspan_total),:);
    tout = tout(ismember(tout, tspan_total));
end

simtime = toc(simtime_start);

%% (3) check output

if ~isnan(tout)
    if length(tout) < length(t)
        log = [log 'length(tout) ' num2str(length(tout)) ' != ' num2str(length(t)) ' length(t_ref), returning NaN; '];
        tout = NaN; Xout = NaN;
    end
end

%% (4) post-processing for index analysis

% a-posteriori determine value of states eliminated via conservation laws
if ~isnan(tout) 
    for p = 1:length(I.replaceODE)

        % index of state variable whose ODE is to be replaced
        k = I.replaceODE(p); 

        % indices of remaining states that are part of the replaceODEby
        % states
        remstates = setdiff(I.replaceODEby{p},k);

        % back calculated solution of states whose ODE was replaced
        Xout(:,k) = max(0, Xout(:,k) - sum(Xout(:,remstates),2) );

    end
end

end

%% internal helper functions

% odesolver_loop
function [tout, Xout, log] = odesolver_loop(t, X0, par, I, odefun, options, log)
if isempty([I.pss])
    [tout, Xout] = ode15s(@(t,X) odefunModel(X,par,I,odefun), t, X0, options);
    log = [log 'non-pss solving successfull; '];
else
    warning off;
    for X0_init = [1e2 1e4 1 1e-10] % try different initial cond for pss state, in case DAE solver has problems
        
        X0(I.pss) = X0_init;
        options.InitialSlope = odefunModel(X0,par,I,odefun);     
        try
            lastwarn(''); % clear last warning message
            [tout, Xout] = ode15s(@(t,X) odefunModel(X,par,I,odefun), t, X0, options);
            if isempty(lastwarn)
                log = [log 'pss solving successfull with X0_init = ' num2str(X0_init) '; '];
                break; % everything is fine --> exit for X0_init = ... loop
            else
                error(lastwarn)
            end
        catch e
            tout = NaN; Xout = NaN;  %%% ode could not be solved
            log = [log 'pss solving unsuccessfull: ' e.message ', returning NaN; '];
        end
    end
    warning on;
end
end

% added to ode options if timeout parameter supplied
% function status = timeoutFun(t,X,flag,startTime,timeout)
% 
% % read elapsed time
% elapsedTime = toc(startTime);
% 
% if elapsedTime > timeout
%     status = 1; % stop integration
% else
%     status = 0; % continue integration
% end
% 
% end


