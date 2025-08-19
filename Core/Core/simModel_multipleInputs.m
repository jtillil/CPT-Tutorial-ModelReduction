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
function [tout, Xout, log] = simModel_multipleInputs(t, X0, par, I, param, multiple, odefun, jacfun)

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

% set all reaction rate constants related to the species in I.cneg (both,
% forward and backward reactions) to zero 
if ~isempty(I.cneg)
    I_par = param.states2Ipar(I.cneg);
    par([I_par{:}]) = 0;
end

% jacobian, if provided
if isempty(jacfun)
    options.Jacobian = @(t,X) jacfunModel(X,par,I,jacfun);
end

% check consistency of of states that are part of I.replaceODE
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

%% (2) solve ODEs;  (i) no pss states, (ii) with pss states

% nonegative for no pss states
if isempty([I.pss])
    options.NonNegative = 1:I.nstates;
end

% simulate model ODE
[t,X,te] = ode15s(@(t,X) model.odefun(t,X,par,model), t, X0, options);
if ~isempty(model.event) && t(end)~=tspan(end)
    t=[t;tspan(end)];
    X(end+1,:)=0;
end
if isfield(model,'multiple') && model.multiple
    input_events=model.input_events;
    input_doses=model.input_doses;
    u_ref=model.u_ref;
    t_out=[];
    x_out=[];
    te_out=[];
    tspan_local=t(:);
    X0_ref=X0(:);
   
    if tspan_local(1)<input_events(1)
        if length(tspan_local)==2
            tspan = [tspan_local(1) input_events(1)];
        else
            tspan = tspan_local(tspan_local<input_events(1));
        end
        [t_out,x_out,te_out] = ode15s(@(t,X) model.odefun(t,X,par,model), tspan, X0, options);
        X0_ref=x_out(end,:)' + u_ref;
    end
    for i=2:length(input_events)
        if input_events(i)<=tspan_local(end) && tspan_local(1) < input_events(i)
            %%% in case of specification of specific timepoints to be solved need
            %%% to account for this
            if length(tspan_local)<=2
                tspan = [input_events(i-1) input_events(i)];
            else
                tspan = unique([input_events(i-1);tspan_local(tspan_local>=input_events(i-1) & tspan_local<=input_events(i));input_events(i)]);
            end
            [t,x,te]  = ode15s(@(t,X) model.odefun(t,X,par,model), tspan, X0_ref, options);
            %%% to avoid double timepoints in the time vector in the multiple
            %%% dosing case, the dosing timepoint is kept and the other discarded
            %%% such that in the state vector only the dosing state is included
            if i<length(input_events)
                t_out = [t_out;t(1:end-1)];
                x_out = [x_out;x(1:end-1,:)];
                te_out = [te_out,te];
            else
                if tspan_local(end)>input_events(end)
                    t_out = [t_out;t(1:end-1)];
                    x_out = [x_out;x(1:end-1,:)];
                    te_out = [te_out,te];
                else
                    t_out = [t_out;t];
                    x_out = [x_out;x];
                    te_out = [te_out,te];
                end
            end
            if any(input_doses)
                u_ref(I.input)=input_doses(i);
            end
            X0_ref=x(end,:)' + u_ref;
        elseif tspan_local(1)>=input_events(i)
            continue;
        else
            break;
        end
    end

    if tspan_local(end)>input_events(end)
        if length(tspan_local)==2
            tspan = [input_events(end) tspan_local(end)];
        else
            tspan = tspan_local(tspan_local>=input_events(end));
        end
        [t,x,te]  = ode15s(@(t,X) model.odefun(t,X,par,model), tspan, X0_ref, options);
        t_out = [t_out;t];
        x_out = [x_out;x];
        te_out = [te_out,te];
    end

    X=x_out(ismember(t_out,tspan_local),:);
    t=t_out(ismember(t_out,tspan_local));
    te=te_out; 

end

if isempty([I.pss])
    
    % default, unless changed below
    options.NonNegative = 1:I.nstates;
    [tout, Xout] = ode15s(@(t,X) odefunModel(X,par,I,odefun), t, X0, options);
    
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
    options.Mass = M; options.MassSingular = 'yes';
    
    warning off;
    for X0_init = [1e2 1e4 1 1e-10] % try different initial cond for pss state, in case DAE solver has problems
        
        X0(I.pss) = X0_init;
        % options.InitialSlope = odefunModel(t(1),X0,model.par,model); 
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


