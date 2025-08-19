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
function [tout, Xout, log, simtime] = simModel(t, X0, par, I, param, multiple, odefun, jacfun, simoptions, odeoptions)

% init log
log = 'Log: ';

% set ode options
if exist("odeoptions", "var")
    options = odeoptions;
end

% check if reduced model provided
if length(I.dyn) < I.nstates % if not all states dynamic
    model_reduced = 1;
else
    model_reduced = 0;
end

% handle lumping [pre and post], ir-index solving
if exist("simoptions", "var")
    if isfield(simoptions, 'prelumpmat')
        prelumping = 1;
    else
        prelumping = 0;
    end
    if isfield(simoptions, 'ir_indices')
        if simoptions.ir_indices
            ir_indices = 1;
        else
            ir_indices = 0;
        end
    else
        ir_indices = 0;
    end
else
    prelumping = 0;
    ir_indices = 0;
end

% init prelumping [and postlumping] matrices
if prelumping
    prelumpmat = simoptions.prelumpmat;
    if ~isfield(simoptions, 'preinvlumpmat')
        preinvlumpmat = pinv(prelumpmat);
    else
        preinvlumpmat = simoptions.preinvlumpmat;
    end
else
    prelumpmat = [];
    preinvlumpmat = [];
end

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

if prelumping
    X0 = prelumpmat * X0;
end

%% (2) specify options

% add jacobian, if provided
if ~isempty(jacfun) && ~prelumping && ~ir_indices
    options.Jacobian = @(t,X) jacfunModel(X,par,I,jacfun);
elseif ir_indices
    options.Jpattern = sparse(extodejacpatfun(I, par, jacfun));
end

% specify options
if isempty([I.pss])
    % default, unless changed below
    if ~prelumping
        options.NonNegative = 1:I.nstates;
    else
        options.NonNegative = 1:size(prelumpmat, 1);
    end
    
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
    if ~ir_indices
        M = eye(I.nstates);
        M(I.pss,I.pss) = 0;
        options.Mass = M;
        options.MassSingular = 'yes';
    else
        M = eye(I.nstates);
        M(I.pss,I.pss) = 0;
        options.Mass = blkdiag(M, kron(M, eye(I.nstates)));
        options.MassSingular = 'yes';

        options.Jacobian = [];
        options.Jpattern = extodejacpatfun(I, par, jacfun);
    end
        % error('pss states not yet implemented for prelumping')
    % end
end

%% (3) simulate model ODEs

simtime_start = tic;

if ~multiple.multiple
    % without multiple dosing
    [tout, Xout, log] = odesolver_loop(t, X0, par, I, odefun, jacfun, options, ir_indices, prelumping, prelumpmat, preinvlumpmat, log);

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
    
    % [tout, Xout, log] = odesolver_loop(tspan_local, X0, par, I, odefun, options, ir_indices, prelumping, prelumpmat, preinvlumpmat, log);
    if tspan_total(1) < input_events(1)
        if length(tspan_total)==2
            tspan_local = [tspan_total(1) input_events(1)];
        else
            tspan_local = tspan_total(tspan_total<input_events(1));
        end
        [tout, Xout, log] = odesolver_loop(tspan_local, X0, par, I, odefun, jacfun, options, ir_indices, prelumping, prelumpmat, preinvlumpmat, log);
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
            [tout_local, Xout_local, log] = odesolver_loop(tspan_local, X0_ref, par, I, odefun, jacfun, options, ir_indices, prelumping, prelumpmat, preinvlumpmat, log);
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
        [tout_local, Xout_local, log] = odesolver_loop(tspan_local, X0_ref, par, I, odefun, jacfun, options, ir_indices, prelumping, prelumpmat, preinvlumpmat, log);
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
        if prelumping
            k = find(prelumpmat(:, k) ~= 0);
        end

        % indices of remaining states that are part of the replaceODEby
        % states
        remstates = setdiff(I.replaceODEby{p},k);
        if prelumping
            for i = 1:length(remstates)
                remstates(i) = find(prelumpmat(:, remstates(i)) ~= 0);
            end
        end

        % back calculated solution of states whose ODE was replaced
        Xout(:,k) = max(0, Xout(:,k) - sum(Xout(:,remstates),2) );

    end
end

end

%% internal helper functions

% odesolver_loop
function [tout, Xout, log] = odesolver_loop(t, X0, par, I, odefun, jacfun, options, ir_indices, prelumping, prelumpmat, preinvlumpmat, log)
if isempty([I.pss])
    if ir_indices
        if ~prelumping
            % odefun == extodefun(which calls odefunModel)
            [tout, Xout] = ode15s(@(t,X) extodefun(X,par,I,odefun,jacfun), t, X0, options);
        else
            % odefun == extodefun(which calls odefunModel)
            [tout, Xout] = ode15s(@(t,X) prelumpmat * extodefun(preinvlumpmat * X,par,I,odefun,jacfun), t, X0, options);
        end
    else
        if ~prelumping
            [tout, Xout] = ode15s(@(t,X) odefunModel(X,par,I,odefun), t, X0, options);
        else
            [tout, Xout] = ode15s(@(t,X) prelumpmat * odefunModel(preinvlumpmat * X,par,I,odefun), t, X0, options);
        end
    end
    log = [log 'non-pss solving successfull; '];
else
    warning off;
    for X0_init = [1 1e-4 1e4 1e-10 1e10] % try different initial conditions for pss state, in case DAE solver has problems
        
        X0(I.pss) = X0_init;
        % options.InitialSlope = odefunModel(X0,par,I,odefun);     
        try
            lastwarn(''); % clear last warning message
            % [tout, Xout] = ode15s(@(t,X) odefunModel(X,par,I,odefun), t, X0, odeoptions);
            if ir_indices
                if ~prelumping
                    % odefun == extodefun(which calls odefunModel)
                    [tout, Xout] = ode15s(@(t,X) extodefun(X,par,I,odefun,jacfun), t, X0, options);
                else
                    % odefun == extodefun(which calls odefunModel)
                    [tout, Xout] = ode15s(@(t,X) prelumpmat * extodefun(preinvlumpmat * X,par,I,odefun,jacfun), t, X0, options);
                end
            else
                if ~prelumping
                    [tout, Xout] = ode15s(@(t,X) odefunModel(X,par,I,odefun), t, X0, options);
                else
                    [tout, Xout] = ode15s(@(t,X) prelumpmat * odefunModel(preinvlumpmat * X,par,I,odefun), t, X0, options);
                end
            end
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

%%% -----------------------------------------------------------------------
%%% extended ODE system to solve the sensitivity equations
%%%
function dextX = extodefun(extX,par,I,odefun,jacfun)

%%% assign model indexing
% I  = model.I;

% decompose extended state vector into states and Wronski matrix
X = extX(1:I.nstates);
W = extX(I.nstates+1:end);

% dX = model.odefun(X,par);
dX = odefunModel(X,par,I,odefun);

%%% calculate wronski matrix
W_matrix = reshape(W,I.nstates,I.nstates);

%%% ODE of Wronski matrix W'(t,s)=df/dx*W(t,s);
% if isempty(jacfun) || prelumping
if isempty(jacfun)
    % determinue finite difference approximation
    dF = numjacfun(X,par,I,odefun);
else
    % use analytical jacobian
    dF = jacfun(X,par);
end

% dF([I.pss I.env I.pneg I.cneg I.irenv_arith, I.irenv_geom I.average I.mode I.constant I.ssenv I.constregr], :) = 0;
dF([I.env I.pneg I.cneg I.irenv_arith, I.irenv_geom I.average I.mode I.constant I.ssenv I.constregr], :) = 0;

dW_matrix = dF * W_matrix;

dextX = [dX;dW_matrix(:)];
end

%%% -----------------------------------------------------------------------
%%% provide numerical approximation of the jacobian, if no analytical
%%% function is provided
%%%
function jac = numjacfun(X,par,I,odefun)

% set relevant parameters for numerical calculation of jacobian
eps = 1e-9;
h0  = sqrt(eps);
n   = length(X);
jac = zeros(n);
E   = eye(n);
% calculate either two sided or one sided numerical jacobian
for k=1:n
    h = max(1,abs(X(k)))*h0;
    if (X(k)-h)>0
        jac(:,k) = (odefunModel(X+h*E(k,:)',par,I,odefun)-odefunModel(X-h*E(k,:)',par,I,odefun))/(2*h);
    else
        jac(:,k) = (odefunModel(X+h*E(k,:)',par,I,odefun)-odefunModel(X,par,I,odefun))/(h);
    end
end

end

%%% -----------------------------------------------------------------------
%%% provide numerical approximation of the jacobian, if no analytical
%%% function is provided
%%%
% function jac = extode_numjacfun(X,par,I,odefun)
% 
% % set relevant parameters for numerical calculation of jacobian
% eps = 1e-9;
% h0  = sqrt(eps);
% n   = length(X);
% jac = zeros(n);
% E   = eye(n);
% % calculate either two sided or one sided numerical jacobian
% for k=1:n
%     h = max(1,abs(X(k)))*h0;
%     if (X(k)-h)>0
%         jac(:,k) = (odefunModel(X+h*E(k,:)',par,I,odefun)-odefunModel(X-h*E(k,:)',par,I,odefun))/(2*h);
%     else
%         jac(:,k) = (odefunModel(X+h*E(k,:)',par,I,odefun)-odefunModel(X,par,I,odefun))/(h);
%     end
% end
% 
% end

%%% -----------------------------------------------------------------------
%%% defines the jacobian sparsity pattern of the extended ODE system 
%%% (only used if no jacobian is provided)
%%%
% function extodejacpat = extodejacpatfun(model)
function extodejacpat = extodejacpatfun(I, par, jacfun)

% number of state variables
nstates = I.nstates;

% initialise output
extodejacpat = zeros(nstates+nstates^2);

% initialize the jacobian 
X_test = ones(nstates,1); % ok for this jacobian
DF = jacfun(X_test,par);

% check if there are NaN or Inf entries
if any(isinf(DF),'all') || any(isnan(DF),'all')
    error('--> There is an issue with the Jacobian: it contains ''inf'' or ''NaN'' entries :-( Please fix.')
end

% determine pattern of non-zero entries and assign it
odejacpat = DF;
odejacpat(odejacpat~=0) = 1;   

extodejacpat(1:nstates,1:nstates) = odejacpat;

% determine pattern of jacobian of Wronski part
for i=1:nstates
    extodejacpat(i*nstates+1:(i+1)*nstates, i*nstates+1:(i+1)*nstates) = odejacpat;
    extodejacpat(i*nstates+1:i*nstates+nstates, 1:nstates) = ones(nstates,nstates);
end

extodejacpat = sparse(extodejacpat);
end


