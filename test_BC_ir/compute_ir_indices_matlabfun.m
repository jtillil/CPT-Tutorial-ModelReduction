%%% Version: 19 Jan 2023
%%%
%%% call by: [ir, contr, obs] = compute_ir_indices(model)
%%%
%%% This function computes the sensitivity based input-response indices for
%%% the given model
%%%
%%% Input:  model               structure specifying the model
%%%         saveresults         1 = yes, 0 = no
%%%
%%% Output: ir                  input-response index 
%%%         contr               controllability index
%%%         obs                 observability index
%%%
%%%
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function [ir, contr, obs] = compute_ir_indices_matlabfun(model,saveresults)

tic;
fprintf('\n Calculate ir, contr & obs indices \n');

% indexing
I = model.I;

% check, whether classification of states is supported
if ~isempty(setxor(I.dyn,1:I.nstates))
    fprintf('\n\n --> ir indices are only computed for models with all states characterised as dynamic; PLEASE FIX! \n\n'); beep; 
    return;
end

% ODE solver options and right hand side of ODE
% options.NonNegative = 1:I.nstates;
% options.AbsTol = [];
% options.RelTol = 1e-2;
% if jacobian specified, also give pattern of jacobian of extODE
% if ~isempty(model.jacfun)
    options.JPattern = extodejacpatfun(model);
% end

% time vectors
t_ref  = model.t_ref; X_ref = model.X_ref;
ntstar = length(t_ref); % number (n) of tstar values

% Initialise variables necessary for analysis
ir.index    = zeros(ntstar,I.nstates);
contr.index = zeros(ntstar,I.nstates);
obs.index   = zeros(ntstar,I.nstates);

% give information about progress of computation
inform = true; 
if inform, fprintf(' Progress report: solve extODE ... '); end

% Set initial value for extended ODE system, containing the ODEs of the
% state variables with indices [1:I.nstates] and the ODEs of the Jacobian 
% (Wronski matrix) with indices [I.nstates+1:(I.nstates+I.nstates^2)].
% In MATLAB, one may alternatively use the indexing [I.nstates+1:end].
% Note: the ODE solver requires a single input state vector, therefor matrices
% like the Jacobian (Wronski) need to be reshaped into vector format. Later
% the vectorised output of the extended ODE system related to the Wronski 
% matrix is again reshaped into matrix form

X0 = X_ref(1,:);        % initial condition of ODE system
W0 = eye(I.nstates);    % initial condition of Wronski matrix equation
extX0 = [X0'; W0(:)];   % initial condition of extended ODE system

% solving the extended ODE system to obtain J_u(0,t*), including reshaping
% the output to obtain Jacobian in matrix form (see above)
[~,extX_ref]  = ode45(@(t,X) extodefun(t,X,model.par,model), t_ref, extX0, options);
% [~,extX_ref]  = ode15s(@(t,X) extodefun(t,X,model.par,model), t_ref, extX0, options);
% [~,extX_ref]  = ode23s(@(t,X) extodefun(t,X,model.par,model), t_ref, extX0, options);

% decompose extended state vector into states and Wronski matrix; initial
% 'e' indicates that X_ref and eX_ref can be expected to differ due to the
% different ODEs used to determine it (original vs. extended ODE)
eX_ref = extX_ref(:,1:I.nstates);
eW_ref = extX_ref(:,I.nstates+1:end);

% extract max
maxeX_ref = max(eX_ref, [], 1);
releX_ref = eX_ref ./ repmat(maxeX_ref, [size(eX_ref, 1), 1]);

%%% define Jacobian J_u(t*,t0)
Jac_u = permute(reshape(eW_ref,length(t_ref),I.nstates,I.nstates),[2,3,1]);

%%% states to determine the indices for
relevantstates = 1:I.nstates;

%%% only for testing purposes
if model.quicktest
    relevantstates = I.output; ntstar = 3; beep;
    fprintf('\n --> quicktest running <-- \n')
end

if inform, fprintf(' and for each t*: '); end
for ts = 1:(ntstar-1)
    
    if inform, fprintf('%d,',ntstar-ts); end
    
    % calculate J_y(t*,t), including reshaping the output to obtain Jacobian 
    % in matrix form (see above) 
    tstarspan  = t_ref(ts:end);
    extX0_tstar = [eX_ref(ts,:)';W0(:)];
    % [t_tstar,extX_tstar]  = ode15s(@(t,X) extodefun(t,X,model.par,model), tstarspan, extX0_tstar, options);
    % [t_tstar,extX_tstar]  = ode23s(@(t,X) extodefun(t,X,model.par,model), tstarspan, extX0_tstar, options);
    [t_tstar,extX_tstar]  = ode45(@(t,X) extodefun(t,X,model.par,model), tstarspan, extX0_tstar, options);

    W_tstar = extX_tstar(:,I.nstates+1:end);
    
    % define Jacobian J_y(T,t*)
    Jac_y  = permute(reshape(W_tstar,length(t_tstar),I.nstates,I.nstates),[2,3,1]);

    % calculate index at timepoint t*, including evaluation of an integral
    % for the observability index via the trapezoidal rule
    if length(t_tstar) == length(tstarspan)
        for k = relevantstates
            obs.index(ts,k)   = sqrt( 1/t_ref(end) * trapz(t_tstar, Jac_y(I.output,k,:).^2) );
            contr.index(ts,k) = abs( Jac_u(k,I.input,ts) );
            ir.index(ts,k)    = contr.index(ts,k) * obs.index(ts,k);
        end
    else
        for k = relevantstates
            obs.index(ts,k)   = NaN;
            contr.index(ts,k) = NaN;
            ir.index(ts,k)    = NaN;
        end
    end
    clear Jac_y extX_tstar
    
end
elapsedtime = toc; fprintf('\n [elapsed time = %.1f]\n\n',elapsedtime);

%%% compute normalized ir index
ir.nindex = diag(1./sum(ir.index,2,'omitnan')) * ir.index;

if saveresults
    fprintf('  \n results saved in %s \n',[model.savenameroot 'ir/contr/obs_index.mat']);
    indx = ir;    save([model.savenameroot '_ir_index.mat'],'indx','model')
    indx = contr; save([model.savenameroot '_contr_index.mat'],'indx','model')
    indx = obs;   save([model.savenameroot '_obs_index.mat'],'indx','model')
end

end


%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% LOCAL SUB-ROUTINES
%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%%% -----------------------------------------------------------------------
%%% extended ODE system to solve the sensitivity equations
%%%
function dextX = extodefun(t,extX,par,model)

%%% assign model indexing
I  = model.I;

% decompose extended state vector into states and Wronski matrix
X = extX(1:I.nstates);
W = extX(I.nstates+1:end);

dX = model.odefun(X,par);

%%% calculate wronski matrix
W_matrix = reshape(W,I.nstates,I.nstates);

%%% ODE of Wronski matrix W'(t,s)=df/dx*W(t,s);
if isempty(model.jacfun)
    % determinue finite difference approximation
    dW_matrix = numjacfun(t,X,model) * W_matrix;
else
    % use analytical jacobian
    dW_matrix = model.jacfun(X,par) * W_matrix;
end

dextX = [dX;dW_matrix(:)];
end

%%% -----------------------------------------------------------------------
%%% provide numerical approximation of the jacobian, if no analytical
%%% function is provided
%%%
function jac = numjacfun(t,X,model)

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
        jac(:,k) = (model.odefun(X+h*E(k,:)',model.par)-model.odefun(X-h*E(k,:)',model.par))/(2*h);
    else
        jac(:,k) = (model.odefun(X+h*E(k,:)',model.par)-model.odefun(X,model.par))/(h);
    end
end

end

%%% -----------------------------------------------------------------------
%%% defines the jacobian sparsity pattern of the extended ODE system 
%%% (only used if no jacobian is provided)
%%%
function extodejacpat = extodejacpatfun(model)

% number of state variables
nstates = model.I.nstates;

% initialise output
extodejacpat = zeros(nstates+nstates^2);

% initialize the jacobian 
X_test = ones(model.I.nstates,1); % ok for this jacobian
DF = model.jacfun(X_test,model.par);

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


