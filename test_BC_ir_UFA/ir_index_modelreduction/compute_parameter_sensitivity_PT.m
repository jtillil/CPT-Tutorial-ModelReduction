%%% Version: March 30th, 2020
%%%
%%% call by: [ir,contr,obs]  =  compute_ir_indices(model)
%%%
%%% This function computes the sensitivity based input-response indices for
%%% the given model
%%%
%%% Input: model                structure specifying the model
%%%
%%% Output: ir                  input-response index 
%%%         contr               controllability index
%%%         obs                 observability index
%%%
%%% Citation:
%%% 
%%% Knoechel, Kloft and Huisinga, "Sensitivity based input-response index to 
%%% analyse and reduce large-scale signalling networks"
%%% PLOS Comp. Biology, 2020 (under review)
%%%
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function [obs]  =  compute_parameter_sensitivity_PT(model)

tic;
fprintf('\n calculate parameter sensitivities \n');

% indexing
I = model.I;
if isfield(model,'options')
    options=model.options;
end

% check, whether classification of states is supported
if ~isempty([I.env, I.neg, I.qss, I.con])
    fprintf(['\n\n --> computation of ir indices not yet supporting ',...
        'anything else then dynamical states. \n    ',...
        'Therefore, all states are assumed to be dynamic for ir computation! \n\n']);
    
    % initialise classification of states for the detailed model
    I.dyn = 1:I.nstates; I.env = []; I.neg = []; I.qss = []; I.con = [];
    
end

% ODE solver options and right hand side of ODE
options.NonNegative = 1:I.nstates;
% if jacobian specified, also give pattern of jacobian of extODE
if ~isempty(model.jacfun)
    options.JPattern = extodejacpatfun(model);
end

% time vectors
t_ref  = model.t_ref;

% Initialise variables necessary for analysis
obs         = zeros(1,I.nstates);

% Set initial value for extended ODE system, containing the ODEs of the
% state variables with indices [1:I.nstates] and the ODEs of the Jacobian 
% (Wronski matrix) with indices [I.nstates+1:(I.nstates+I.nstates^2)].
% In MATLAB, one may alternatively use the indexing [I.nstates+1:end].
% Note: the ODE solver requires a single input state vector, therefor matrices
% like the Jacobian (Wronski) need to be reshaped into vector format. Later
% the vectorised output of the extended ODE system related to the Wronski 
% matrix is again reshaped into matrix form

x0 = model.X_ref(1,:); % initial condition of ODE system
W0 = eye(I.nstates);    % initial condition of Wronski matrix equation

if ~isempty(model.jacfun)
    model.options.Jacobian = @(t,X) model.jacfun(t,X,model.par,model);
end
model.options.NonNegative = 1:I.nstates;

% solving the ODE system
[~,X_ref]  = ode15s(@(t,X) model.odefun(t,X,model.par,model), t_ref, x0, model.options);

nstates=I.nstates;
par=model.par;

% calculate J_y(t*,t), including reshaping the output to obtain Jacobian 
% in matrix form (see above) 
tstarspan  = t_ref;
tstarspan = [tstarspan(tstarspan<0.003195);0.003195];
extX0_tstar = [X_ref(1,:)';W0(:)];
%options.Event=@(t,concentrations) event_fct(t,concentrations,1500/3600,model.I.AUC,1,1);
    [t_tstar,extX_tstar]  = ode15s(@(t,X) extodefun(t,X,par,model), tstarspan, extX0_tstar, options);
    W_tstar = extX_tstar(:,nstates+1:end);        

% define Jacobian J_y(T,t*)
Jac_y  = permute(reshape(W_tstar,length(t_tstar),nstates,nstates),[2,3,1]);

% calculate index at timepoint t*, including evaluation of an integral
% for the observability index via the trapezoidal rule
for k = 1:nstates
    obs(k)   = abs( Jac_y(I.output,k,end) );
end

elapsedtime = toc; fprintf('\n[elapsed time = %.1f]\n\n',elapsedtime);
end


%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% LOCAL SUB-ROUTINES
%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%%% -----------------------------------------------------------------------
%%% extended ODE system to solve the sensitivity equations

function dextX = extodefun(t,extX,par,model)

%%% assign model indexing
I  = model.I;

% decompose extended state vector into states and Wronski matrix
X = extX(1:I.nstates);
W = extX(I.nstates+1:end);

dX = model.odefun(t,X,par,model);

%%% calculate wronski matrix
W_matrix = reshape(W,I.nstates,I.nstates);

%%% ODE of Wronski matrix W'(t,s)=df/dx*W(t,s);
if isempty(model.jacfun)
    % determinue finite difference approximation
    dW_matrix = numjacfun(t,X,model) * W_matrix;
else
    % use analytical jacobian
    dW_matrix = model.jacfun(t,X,par,model) * W_matrix;
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
        jac(:,k) = (model.odefun(t,X+h*E(k,:)',model.par,model)-model.odefun(t,X-h*E(k,:)',model.par,model))/(2*h);
    else
        jac(:,k) = (model.odefun(t,X+h*E(k,:)',model.par,model)-model.odefun(t,X,model.par,model))/(h);
    end
end

end

%%% -----------------------------------------------------------------------
%%% defines the jacobian sparsity pattern of the extended ODE system 
%%%
function extodejacpat = extodejacpatfun(model)

% number of state variables
nstates = model.I.nstates;

% initialise output
extodejacpat = sparse(nstates+nstates^2);

% initialize the jacobian 
X_test = ones(model.I.nstates,1); % ok for this jacobian
DF = model.jacfun(0,X_test,model.par,model);

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

function [position, isterminal, direction] = event_fct(~,concentrations,limit,index,direction,terminal)
% event function: stop simulation when index'th concentration >= limit 
% (<= if direction=-1)
position = concentrations(index)-limit;
isterminal = terminal; %1 if terminal not defined

end