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
function [tout, Xout, log] = simModel_model(t, X0, par, model)

% read indexing
I = model.I;

% if pss present, check if model contains solved pss 
if ~isempty(I.pss)
    
end

% init log
log = 'Log: ';

%% (1) pre-processing for index analysis

% set negligible state variables to zero (to realise that they are not
% part of the model
X0([I.pneg I.cneg]) = 0;

% set irenv state variables to ir-weighted initial state
% if (isfield(model, 'irenv'))
%     X0(I.irenv) = model.irenv(I.irenv);
% end
if (isfield(model, 'irenv_arith'))
    X0(I.irenv_arith) = model.irenv_arith(I.irenv_arith);
end
if (isfield(model, 'irenv_geom'))
    X0(I.irenv_geom) = model.irenv_geom(I.irenv_geom);
end
if (isfield(I, 'Irenv_arith'))
    X0(I.irenv_arith) = I.Irenv_arith(I.irenv_arith);
end
if (isfield(I, 'Irenv_geom'))
    X0(I.irenv_geom) = I.Irenv_geom(I.irenv_geom);
end

% set all reaction rate constants related to the species in I.cneg (both,
% forward and backward reactions) to zero 
if ~isempty(I.cneg)
    I_par = model.param.states2Ipar(I.cneg);
    par([I_par{:}]) = 0;
end

% jacobian, if provided
if ~isempty(model.jacfun)
    options.Jacobian = @(t,X) jacfunModel(t,X,par,model);
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

% check if timeout supplied
% if nargin > 4
%     % start timeout counting
%     startTime = tic;
% 
%     % add timeout function to options
%     options.OutputFcn = @(t,X,flag) timeoutFun(t,X,flag,startTime,timeout);
% end

if isempty([I.pss])
    
    % default, unless changed below
    options.NonNegative = 1:I.nstates;
    [tout, Xout] = ode15s(@(t,X) odefunModel(t,X,par,model), t, X0, options);
    
else
    
    % DAE solver not able to deal with NonNegative option
    options.NonNegative = [];

    % so far unresolved problem, if Jacobian is provided for DAE case with
    % simultaneous conservation laws; not providing Jacobian seems to be 
    % an interim solution
    if ~isempty(I.replaceODE) && ~isempty(model.jacfun)
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
        options.InitialSlope = odefunModel(t(1),X0,par,model);     
        try
            lastwarn(''); % clear last warning message
            [tout, Xout] = ode23t(@(t,X) odefunModel(t,X,par,model), t, X0, options);
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
    if length(tout) ~= length(t)
        log = [log 'length(tout) ' num2str(length(tout)) ' != ' num2str(length(model.t_ref)) ' length(t_ref), returning NaN; '];
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
        X(:,k) = max(0, X(:,k) - sum(X(:,remstates),2) );

    end
end

end

%% internal helper functions

% accesses model.ode
function dX = odefunModel(t,X,par,model)

%%% assign model indexing
I  = model.I;

%%% (1) pre-processing for index analysis
%%% set negliglibe states variables to zero and account for intervention
%%%
X([I.pneg I.cneg]) = 0;

% back calculated value of states whose ODE was replaced
for p = 1:length(I.replaceODE)

    % if a state (say the k-th one) is part of I.replaceODE, then its
    % initial value X0(k) and its ODE dX(k) are replaced by the sum of
    % inital values and ODEs of the states that are specified in
    % I.replaceODEby

    % Here, this step is 'undone' by determining the value of the k-th 
    % state from the values of the states in replaceODEby
    % Enforce that difference is non-negative (ODE solver accurracy
    % might otherwise result in negative values)

    k = I.replaceODE(p);    % index of state variable

    remstates = setdiff(I.replaceODEby{p},k);
    X(k) = max(0, X(k) - sum(X(remstates)) );

end

%%% (2) call model ode 
dX = model.odefun(X,par);

%%% (3) post-processing for index analysis
%%% set derivative of environmental & negligible states variables to zero
%%%
dX([I.env I.pneg I.cneg I.irenv_arith, I.irenv_geom]) = 0; 

%%% determine the ODEs for the states in replaceODE
for p = 1:length(I.replaceODE)
    k  = I.replaceODE(p);       % index of state to be replaced
    states = I.replaceODEby{p}; % indices of states used for replacement
    dX(k) = sum(dX(states));
end

end

% accesses model.jac
function DF = jacfunModel(t,X,par,model)

%%% assign model indexing
I  = model.I;

%%% (1) pre-processing for index analysis
%%% set negliglibe states variables to zero and account for intervention
%%%
X([I.pneg I.cneg]) = 0;

% back calculated value of states whose ODE was replaced
for p = 1:length(I.replaceODE)

    % if a state (say the k-th one) is part of I.replaceODE, then its
    % initial value X0(k) and its ODE dX(k) are replaced by the sum of
    % inital values and ODEs of the states that are specified in
    % I.replaceODEby

    % Here, this step is 'undone' by determining the value of the k-th 
    % state from the values of the states in replaceODEby
    % Enforce that difference is non-negative (ODE solver accurracy
    % might otherwise result in negative values)

    k = I.replaceODE(p);    % index of state variable

    remstates = setdiff(I.replaceODEby{p},k);
    X(k) = max(0, X(k) - sum(X(remstates)) );

end

%%% (2) call model jacobian 
DF = model.jacfun(X,par);

%%% (3) post-processing for index analysis
%%% check if there are NaN or Inf entries
if any(isinf(DF),'all') || any(isnan(DF),'all')
    error('\n There is an issue with the Jacobian (contains ''inf'' or ''NaN'' entries).')
end

%%% set jacobian of environmental & negligible states variables to zero
%%%
DF([I.env I.pneg I.cneg I.irenv_arith, I.irenv_geom],:) = 0;

%%% determine the Jacobian for the states in replaceODE
for p = 1:length(I.replaceODE)
    k  = I.replaceODE(p);       % index of state to be replaced
    states = I.replaceODEby{p}; % indices of states used for replacement
    DF(k,:) = sum(DF(states,:));
end

end

% added to ode options of timeout parameter supplied
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


