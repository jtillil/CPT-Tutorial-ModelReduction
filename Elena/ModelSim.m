function model = ModelSim(model)
%==========================================================================
%-- sources:
%  
% Andrea Y. Weiße, Diego A. Oyarzún, Vincent Danos, Peter S. Swain:
% 
% "Cellular trade-offs, gene expression, and growth"
% Proceedings of the National Academy of Sciences Mar 2015
% DOI: 10.1073/pnas.1416533112
%==========================================================================

%%%===========================
%%% Step 1: Set up the model                 
%%%===========================

if ~isfield(model,'I')
    model.I = Index();
end
I = model.I;

% Import parameters
if ~isfield(model,'par')
    model.par = Par(I);
end
par = model.par;

%Import initial conditions
if ~isfield(model,'X0')
model.X0 = InitialCon(I);
end 
X0 = model.X0;

%Import timespan
if ~isfield(model,'tspan')
model.tspan = [0 1e6];
end
tspan = model.tspan;

%Right hand side of ODEs
model.odefun = @(t,X,par,model) ODEs(t,X,par,model);

% Function handle to simulate the ODE
model.simODE = @(t,X,par,model) simulateODE(t,X,par,model);


%%%============================
%%% Step 2: compute solution
%%%============================

[t,X,sol] = model.simODE(model.tspan,X0,par,model);
model.sol = sol;
model.t = t;
model.X = X;
model.lam = calculatelam(model)';
% model.Rfr = calculateRfr(model)';

end
 



%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% LOCAL SUB-ROUTINES
%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%--------------------------------------------------------------------------

function [t,X,sol] = simulateODE(tspan,X0,par,model)
I       = model.I; 

% Solver options 
options.NonNegative = 1:I.nstates; 
options.RelTol      = 1e-6; 
options.AbsTol      = 1e-9;

% Simulate ODE
sol = ode15s(@(t,X) model.odefun(t,X,par,model), tspan, X0,options);
t=sol.x';
X=sol.y';
end


function lam = calculatelam(model)

I       = model.I;
par     = model.par;
X       = model.X;

gamma      = par(I.gmax).*X(:,I.a) ./(par(I.Kgamma)+ X(:,I.a));
lam        = (((X(:,I.cq) + X(:,I.cr) + X(:,I.ct) + X(:,I.cm)).*gamma)./(par(I.Mref)));

end

function Rfr = calculateRfr(model)

I       = model.I;
par     = model.par;
X       = model.X;

Rfr     =  par(I.nr).*(X(:,I.r)+X(:,I.cr)+X(:,I.ct)+X(:,I.cm)...
    +X(:,I.cq)+X(:,I.zmr)+X(:,I.zmt)+X(:,I.zmm)+X(:,I.zmq))  ./...
    (par(I.nr)*(X(:,I.r)+X(:,I.cr)+X(:,I.ct)+X(:,I.cm)+X(:,I.cq)...
    +X(:,I.zmr)+X(:,I.zmt)+X(:,I.zmm)+X(:,I.zmq))+ par(I.nx).*...
    (X(:,I.q)+X(:,I.et)+X(:,I.em)));

end


