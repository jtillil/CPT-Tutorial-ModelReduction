function model = set_up_the_model(model)

%%% =======================================================================
%%% 1st part - files to set up the model
%%% =======================================================================

%%% To define a model, there are only six files needed  
%%% (in folder 'modelspecification/'):
%%%
%%%     a) MODELNAME_indexing.m
%%%     b) MODELNAME_initialvalues.m
%%%     c) MODELNAME_parameters.m
%%%     d) MODELNAME_ode.m
%%%     e) MODELNAME_species2params.m
%%%     f) MODELNAME_model_set_up_details.m
%%%
%%% To speed up simulation, a seventh file is helpful 
%%% (in the folder 'modelspecification/'):
%%%
%%%     f) MODELNAME_odejac.m
%%% 

% indexing of states and parameters
filename = [model.name '_indexing'];
check_if_filename_exists(filename);
I = feval(filename);
model.I = I;

% initial values
filename = [model.name '_initialvalues'];
check_if_filename_exists(filename);
X0 = feval(filename,model);
model.X0 = X0;
model.X0prior2input = X0;

% parameter values
filename = [model.name '_parameters'];
check_if_filename_exists(filename);
par = feval(filename,model);
model.par = par;

% determine parameters related to each species
% if exist([model.name '_species2params.mat'], "file") == 2
    param = feval([model.name '_species2params'],model);
    model.param = param;
% end

% define input
I.input  = model.setup.input;
u_ref = zeros(size(X0)); 
u_ref(I.input) = model.setup.u_ref;
model.unit.input = model.setup.unit.input;

% define output 
I.output = model.setup.output;
model.unit.output = model.setup.unit.output;
model.unit.graphic.transfoutput = model.setup.unit.graphic.transfoutput;
model.unit.graphic.output = model.setup.unit.graphic.output;
model.unit.graphic.ylim = model.setup.unit.graphic.ylim;

% simulation time span and unit transformation
model.tspan = model.setup.tspan;  
model.unit.time = model.setup.unit.time;
model.unit.graphic.transftime = model.setup.unit.graphic.transftime;
model.unit.graphic.time = model.setup.unit.graphic.time;
model.unit.graphic.xlim = model.setup.unit.graphic.xlim;

% initial value of reference model (incl. pre-stimulus state and input)
model.X0 = model.X0prior2input + u_ref;

% right hand side of model ODEs
filename = [model.name '_ode'];
check_if_filename_exists(filename);
model.ode = @(t,X,par,model) feval(filename,t,X,par,model);

% ODE function including pre- and post-processes for indices
model.odefun = @(t,X,par,model) odefun4indices(t,X,par,model);

% optional: jacobian of ODE and jacobian pattern of extended ODE 
% false = no jacobian provided 
% true = jacobian provided 
model.jacfun = []; 
if model.setup.jacfunprovided
    filename = [model.name '_odejac'];
    check_if_filename_exists(filename);
    model.jac = @(t,X,par,model) feval(filename,t,X,par,model);

    % Jacobian function including pre- and post-processes for indices
    model.jacfun = @(t,X,par,model) jacfun4indices(t,X,par,model);
end

% conservation law for this small example was set manually
if model.setup.conlaw.n > 0

    % number (n) of conlaws
    L.nconlaw = model.setup.conlaw.n;

    % all fields of model.conlaw are names (nm) of conlaws, except for 'n'
    L.nmconlaw = fields(rmfield(model.setup.conlaw,'n')); 
    for k = 1:L.nconlaw
        nmcl = L.nmconlaw{k};
        L.(nmcl) = [];
        L.(nmcl).states = model.setup.conlaw.(nmcl).states;
    end

model.L = L;
end

% function handle to simulate the ODE, or in case of pss states, the DAE  
model.simODE = @(t,X,par,model) simulateODE(t,X,par,model);

%to compute difference between reference solution and
% approximate solution
model.relerrnorm = @(t_ref,Y_diff,Y_ref) relerrnorm(t_ref,Y_diff,Y_ref);
model.abserrnorm = @(t_ref,Y_diff) abserrnorm(t_ref,Y_diff);

% define lenged lables and colors for all state variables so that for
% state variables the same color is used in every plot
model.state2col = colormap('lines');

if model.setup.legendlabelsprovided 
    filename = [model.name '_legendlabels'];
    check_if_filename_exists(filename);
    I = feval(filename,I);
else
    I = default_legendlabels(I); % 
end

if model.setup.linestyleprovided 
    filename = [model.name '_state2linestyle'];
    check_if_filename_exists(filename);
    model.state2linestyle = feval(filename,I);
else
    model.state2linestyle = default_state2linestyle(I); % 
end


%%% =======================================================================
%%% 2nd part - specify & initialize default state classification   
%%% =======================================================================

% initialise classification of states for the reference model
% input-response (dyn); environmental (env); unimportant (neg);
% short-lived (pss); conserved (con); intervention (int)
I.nm_indices = {'ir','contr','obs','env','pss','tss','pneg','cneg','conlaw'};
I.nm_non_ir_indices = {'env','pss','tss','pneg','cneg'};

I.stateclasses = {'dyn','env','pss','tss','pneg','cneg','conlaw'};
I.nondynstateclasses = {'env','pss','tss','pneg','cneg','conlaw'};
for sc = I.stateclasses
    I.(sc{:}) = [];
end
I.replaceODE = []; % use original system of ODEs
I.dyn = 1:I.nstates; 
model.I = I;


%%% =======================================================================
%%% 3rd part - compute reference solution
%%% =======================================================================

yesno = {'yes','no'};
fprintf('\n  Compute reference solution (jacobian provided: %s)',yesno{isempty(model.jacfun)+1});

% compute reference solution
tic;
[t_ref,X_ref] = model.simODE(model.tspan,model.X0,par,model);
elapsedtime = toc; fprintf(' [simulation of full model = %.1f sec]\n\n',elapsedtime);

% assign values to model structure
model.t_ref = t_ref; model.X_ref = X_ref; model.Y_ref = X_ref(:,I.output);

model.AUCoutput = sqrt( 1/t_ref(end) * trapz(t_ref, (X_ref(:,I.output)).^2) );

% define observed state variables (to be plottet) and plot
model.obsstates = model.setup.obsstates;

%%% plotting for reference output
model.plotrefsolution          = @(figNr,varargin) plotrefsolution(figNr,model,varargin{:});
model.plotsolutionwithoutinput = @(figNr,varargin) plotsolutionwithoutinput(figNr,model,varargin{:});
model.plotredmodelsolution     = @(figNr,Ired,varargin) plotredmodelsolution(figNr,model,Ired,varargin{:});
model.plotredparmodelsolution  = @(figNr,I_states,varargin) plotredparmodelsolution(figNr,model,I_states,varargin{:});

%%% make bold plots
model.makeplotbold = @(figNr) makeplotbold(figNr);

%%% display refrerence output
model.plotrefsolution(1,model.obsstates);

end


%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% LOCAL SUB-ROUTINES
%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -------------------------------------------------------------------------
function dX = odefun4indices(t,X,par,model)

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
dX = model.ode(t,X,par,model);

%%% (3) post-processing for index analysis
%%% set derivative of environmental & negligible states variables to zero
%%%
dX([I.env I.pneg I.cneg]) = 0; 

%%% determine the ODEs for the states in replaceODE
for p = 1:length(I.replaceODE)
    k  = I.replaceODE(p);       % index of state to be replaced
    states = I.replaceODEby{p}; % indices of states used for replacement
    dX(k) = sum(dX(states));
end

end

% -------------------------------------------------------------------------
function DF = jacfun4indices(t,X,par,model)

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
DF = model.jac(t,X,par,model);

%%% (3) post-processing for index analysis
%%% check if there are NaN or Inf entries
if any(isinf(DF),'all') || any(isnan(DF),'all')
    error('\n There is an issue with the Jacobian (contains ''inf'' or ''NaN'' entries).')
end

%%% set jacobian of environmental & negligible states variables to zero
%%%
DF([I.env I.pneg I.cneg],:) = 0; 

%%% determine the Jacobian for the states in replaceODE
for p = 1:length(I.replaceODE)
    k  = I.replaceODE(p);       % index of state to be replaced
    states = I.replaceODEby{p}; % indices of states used for replacement
    DF(k,:) = sum(DF(states,:));
end

end

% -------------------------------------------------------------------------
function [t,X] = simulateODE(tspan,X0,par,model)

% indexing
I = model.I; 

% sharper tolerances (sometimes needed to resolve numerical problems
% (too small step size?) 
% options.RelTol = 1e-6; options.AbsTol = 1e-9;

%%% (1) pre-processing for index analysis

% set negligible state variables to zero (to realise that they are not
% part of the model
X0([I.pneg I.cneg]) = 0;

% set all reaction rate constants related to the species in I.cneg (both,
% forward and backward reactions) to zero 
if ~isempty(I.cneg)
    I_par = model.param.states2Ipar(I.cneg);
    par([I_par{:}]) = 0;
end

% jacobian, if provided
if ~isempty(model.jacfun)
    if model.ode_is_matlabfun
        options.Jacobian = @(t,X) model.jacfun(X,par);
    else
        options.Jacobian = @(t,X) model.jacfun(t,X,par,model);
    end
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
    if any( ismember( setdiff(I.replaceODE,k), setdiff(states,k) ) )
        commonstates  = intersect(setdiff(I.replaceODE,k),states);
        fprintf('\n\n --> The state %s in replaceODE does belong to multiple sets of states in replaceODEby :-( --- please fix \n\n',I.nmstate{commonstates(1)});
        error('--> Ending here')      
    end

    % set initial value of replaceODE state
    X0(k)  = sum(X0(states));

end   

%%% (2) solve ODEs;  (i) no pss states, (ii) with pss states

if isempty([I.pss])
    
    % default, unless changed below
    options.NonNegative = 1:I.nstates;
    if model.ode_is_matlabfun
        [t,X] = ode15s(@(t,X) model.odefun(X,par), tspan, X0, options);
    else
        [t,X] = ode15s(@(t,X) model.odefun(t,X,par,model), tspan, X0, options);
    end
    
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
        if model.ode_is_matlabfun
            options.InitialSlope = model.odefun(X0,model.par);
        else
            options.InitialSlope = model.odefun(tspan(1),X0,model.par,model);
        end
        try    
            lastwarn(''); % clear last warning message
            if model.ode_is_matlabfun
                [t,X] = ode15s(@(t,X) model.odefun(X,par), tspan, X0, options);
            else
                [t,X] = ode15s(@(t,X) model.odefun(t,X,par,model), tspan, X0, options);
            end
            if isempty(lastwarn)
                break; % everything is fine --> exit for X0_init = ... loop
            else
                error('some warning occurred')
            end
        catch
            t = NaN; X = NaN;  %%% ode could not be solved
        end      
    end
    warning on;
       
    
end

%%% (3) post-processing for index analysis

% a-posteriori determine value of states eliminated via conservation laws
if ~isnan(t) 
    for p = 1:length(I.replaceODE)

        % index of state variable whose ODE is to be replaced
        k = I.replaceODE(p); 

        % indices of remaining states that are part of the replaceODEby
        % states
        remstates = setdiff(I.replaceODEby{p},k);

        % back calculated solution of state whose ODE was replaced
        X(:,k) = max(0, X(:,k) - sum(X(:,remstates),2) );

    end
end


end

% -------------------------------------------------------------------------
function err = errorfun(t_ref,Y_ref,Y_appr)

%%% determines the error between Y_ref and Y_appr based on minimum of the
%%% absolute and the relative error as is done, e.g., in a similar way 
%%% for the numerical integrator

abs_err = sqrt( trapz(t_ref,(Y_ref-Y_appr).^2) );
rel_err = abs_err / sqrt( trapz(t_ref,Y_ref.^2) );

err = rel_err; %min(abs_err,rel_err);

end

% -------------------------------------------------------------------------
function abs_err = abserrnorm(t_ref,Y_diff)

%%% determines the time-normalized L2-norm of Y_diff

T_end = t_ref(end);
abs_err = sqrt( 1/T_end * trapz(t_ref,Y_diff.^2) );

end

% -------------------------------------------------------------------------
function rel_err = relerrnorm(t_ref,Y_diff,Y_ref)

%%% determines the L2-norm of Y_diff normalized by the L2-norm of Y_ref
%%% Note: the normalisation constant 1/T_end of the abs_err cancels out

rel_err = sqrt( trapz(t_ref,Y_diff.^2) / trapz(t_ref,Y_ref.^2) );

end

% -------------------------------------------------------------------------
function [] = plotrefsolution(figNr,model,varargin)

I = model.I;

%%% states for graphical output
state2col = model.state2col;
state2linestyle = model.state2linestyle;
states2plot = 1:I.nstates;
if nargin == 3
    states2plot = varargin{1};
end

%%% potential unit conversion
unit = model.unit;
t_ref = unit.graphic.transftime(model.t_ref);
X_ref = unit.graphic.transfoutput(model.X_ref);

figure(figNr); clf
for k = 1:length(states2plot)
    s = states2plot(k);
    plot(t_ref,X_ref(:,states2plot(k)),'LineWidth',2,'LineStyle',state2linestyle{s},'Color',state2col(s,:));
    hold on;
end
xlabel(sprintf('t [%s]',unit.graphic.time)); ylabel(sprintf('concentration [%s]',unit.graphic.output)); 
title('reference solution of original model')
set(gca,'yscale','log'); legend(I.nmstatelegend(states2plot),'Location','eastoutside');
% zoom in/out on the Y axis by changing the command below
xlim(model.unit.graphic.xlim); ylim(model.unit.graphic.ylim);
makeplotbold(figNr); drawnow;

end

% -------------------------------------------------------------------------
function [] = plotsolutionwithoutinput(figNr,model,varargin)

fprintf('\n Compute solution without input');
I = model.I;

% compute reference solution
tic;
[t_ref,X_ref] = model.simODE(model.tspan,model.X0prior2input,model.par,model);
elapsedtime = toc; fprintf('\n [simulation of full model = %.1f sec]\n\n',elapsedtime);

%%% states for graphical output
state2col = model.state2col;
states2plot = 1:I.nstates;
if nargin == 3
    states2plot = varargin{1};
end

%%% potential unit conversion
unit = model.unit;
t_ref = unit.graphic.transftime(t_ref);
X_ref = unit.graphic.transfoutput(X_ref);

% graphical output
figure(figNr); clf
for k = 1:length(states2plot)
    plot(t_ref,X_ref(:,states2plot(k)),'LineWidth',2,'Color',state2col(states2plot(k),:));
    hold on;
end
xlabel(sprintf('t [%s]',unit.graphic.time)); ylabel(sprintf('concentration [%s]',unit.graphic.output)); 
title('solution of original model without input')
set(gca,'yscale','log'); legend(I.nmstatelegend(states2plot),'Location','eastoutside');
% zoom in/out on the Y axis by changing the command below
xlim(model.unit.graphic.xlim); ylim(model.unit.graphic.ylim);
makeplotbold(figNr); drawnow;

%%% brief analysis
%%%
fprintf(' \n check for steady state (rel change  = X_ref(%1.2f)/X_ref(0)) \n',model.tspan(end));

relchange = X_ref(end,:)./X_ref(1,:);
for k = 1:model.I.nstates
    
    if relchange(k)<0.99
       fprintf('\n rel change of %s below 0.99 (not in steady state!)',I.nmstate{k});
    elseif relchange(k)>1.01
       fprintf('\n rel change of %s above 1.01 (not in steady state!)',I.nmstate{k});
    end
    
end
fprintf('\n\n')

end

% -------------------------------------------------------------------------
function [] = plotredmodelsolution(figNr,model,Ired,varargin)

%%% set up perturbed (reduced) model
I = model.I; sumofsizes = 0; I_nondynstates = []; 

for s = I.nondynstateclasses
    class = s{:};
    if isfield(Ired,class)
        % sum of sizes needed to check below for pairwise disjointness
        sumofsizes = sumofsizes + length(Ired.(class));
        I_nondynstates = [I_nondynstates Ired.(class)];
    else
        Ired.(class) = [];
    end   
end
Ired.dyn = setdiff(1:I.nstates,I_nondynstates);

%%% consistency check for Ired: should be pairwise disjoint; use the fact: 
%%% finite sets are pairwise disjoint if size of their union equals 
%%% the sum of their sizes
if ~( length(unique(I_nondynstates)) == sumofsizes )
    fprintf('\n\n --> Set specified in Ired are not pairwise disjoint; PLEASE FIX! \n\n'); beep; 
    return;
end

%%% enrich a Ired by all fields of I that are not yet field of Ired
for fnm = fieldnames(I)'
    if ~isfield(Ired,fnm{:})
        Ired.(fnm{:}) = I.(fnm{:});
    end
end
redmodel = model; redmodel.I = Ired;

fprintf('\n Compute solution of reduced model');

%!% options.RelTol = 1e-6; options.AbsTol = 1e-9;

% compute preturbed model solution
tic;
[~,X_appr] = redmodel.simODE(model.t_ref,redmodel.X0,redmodel.par,redmodel);
elapsedtime = toc; fprintf('\n [simulation of full model = %.1f sec]\n\n',elapsedtime);

%%% states for graphical output
states2plot = 1:I.nstates;
if nargin >= 4
    states2plot = varargin{1};
end
state2col = model.state2col;

%%% potential unit conversion
unit = model.unit;
t_ref = unit.graphic.transftime(model.t_ref);
X_ref = unit.graphic.transfoutput(model.X_ref);
if ~isnan(X_appr)
    X_appr = unit.graphic.transfoutput(X_appr);
    %model.errfun = @(t_ref,Y_ref,Y_appr) sqrt( trapz(t_ref,(X_ref(:,I.output)-X_per(:,I.output)).^2) ) / sqrt( trapz(t_ref,X_ref(:,I.output).^2) );

    rel_err_per = model.relerrnorm(t_ref,X_appr(:,I.output)-X_ref(:,I.output),X_ref(:,I.output));
    fprintf('\n   relative output error of perturbed model = %1.3f\n',rel_err_per);
else
    fprintf('\n\n --> ODE model cannot be integrated! Probably due to infeasible partial steady state(s). NO OUTPUT! \n\n'); beep;
    figure(figNr); clf; plot(1,NaN,'r*')
    legend(sprintf(' ODE model cannot be integrated ! \n Probably due to infeasible partial steady state(s).\n NO OUTPUT!'),'Fontsize',16)
    return;
end

%%% notify, if output is very small (potentiall zero), due to imcomplete
%%% signal transduction
if all(X_appr(:,I.output)<1)
    fprintf('\n\n Note: perturbed output below 1 %s (max = %1.3e)! \n\n',unit.graphic.output,max(X_appr(:,I.output))); beep
end

%%% plot additional quantity: sum of X_ref and X_per for all states
%%% specified as varargin{2}
sumofstates = []; nmsumofstateslegend = '';
if nargin >= 5
    sumofstates = varargin{2};
    if nargin >= 6
        nmsumofstateslegend = varargin{3};
    end
end
sumX_ref = sum(X_ref(:,sumofstates),2);
sumX_per = sum(X_appr(:,sumofstates),2);

%%% graphical output
figure(figNr); clf
for k = 1:length(states2plot)
    plot(t_ref,X_ref(:,states2plot(k)),'Linestyle','-','LineWidth',2,'Color',state2col(states2plot(k),:));
    hold on;
end
for k = 1:length(states2plot)
    plot(t_ref,X_appr(:,states2plot(k)),'Linestyle','--','LineWidth',2,'Color',state2col(states2plot(k),:));
end
plot(t_ref,sumX_ref,'k-','LineWidth',2);
plot(t_ref,sumX_per,'k--','LineWidth',2);
xlabel(sprintf('t [%s]',unit.graphic.time)); ylabel(sprintf('concentration [%s]',unit.graphic.output)); 
title('reduced (--) & original (-) model output')
set(gca,'yscale','log'); legend([I.nmstatelegend(states2plot) nmsumofstateslegend],'Location','eastoutside');
% zoom in/out on the Y axis by changing the command below
xlim(model.unit.graphic.xlim); ylim(model.unit.graphic.ylim);
makeplotbold(figNr); drawnow;

end

% -------------------------------------------------------------------------
function [] = plotredparmodelsolution(figNr,model,I_par,varargin)

%%% set up perturbed (reduced) model
I = model.I; redmodel = model; 

% set parameters of reduced model to zero
redmodel.par(I_par) = 0;

fprintf('\n Compute solution of reduced parameter model');

% compute preturbed model solution
tic;
[~,X_appr] = redmodel.simODE(model.t_ref,redmodel.X0,redmodel.par,redmodel);
elapsedtime = toc; fprintf('\n [simulation of full model = %.1f sec]\n\n',elapsedtime);

%%% states for graphical output
state2col = model.state2col;
states2plot = 1:I.nstates;
if nargin >= 4
    states2plot = varargin{1};
end

%%% potential unit conversion
unit = model.unit;
t_ref = unit.graphic.transftime(model.t_ref);
X_ref = unit.graphic.transfoutput(model.X_ref);
if ~isnan(X_appr)
    X_appr = unit.graphic.transfoutput(X_appr);
else
    fprintf('\n\n --> Reduced ODE model cannot be integrated. NO OUTPUT! \n\n'); beep;
    return;
end

%%% notify, if output is very small (potentiall zero, due to imcomplete
%%% signal transduction
if all(X_appr(:,I.output)<1)
    fprintf('\n\n Note: perturbed output below 1 %s (max = %1.3e)! \n\n',unit.graphic.output,max(X_appr(:,I.output))); beep
end

%%% plot additional quantity: sum of X_ref and X_per for all states
%%% specified as varargin{2}
sumofstates = [];
if nargin >= 5
    sumofstates = varargin{2};
end
sumX_ref = sum(X_ref(:,sumofstates),2);
sumX_per = sum(X_appr(:,sumofstates),2);

%%% graphical output
figure(figNr); clf
for k = 1:length(states2plot)
    plot(t_ref,X_ref(:,states2plot(k)),'Linestyle','-','LineWidth',2,'Color',state2col(states2plot(k),:));
    hold on;
end
for k = 1:length(states2plot)
    plot(t_ref,X_appr(:,states2plot(k)),'Linestyle','--','LineWidth',2,'Color',state2col(states2plot(k),:));
end
plot(t_ref,sumX_ref,'k-',t_ref,sumX_per,'k--','LineWidth',2);
xlabel(sprintf('t [%s]',unit.graphic.time)); ylabel(sprintf('concentration [%s]',unit.graphic.output)); 
title('reduced (--) & original (-) model output')
set(gca,'yscale','log'); legend(I.nmstatelegend(states2plot),'Location','eastoutside');
% zoom in/out on the Y axis by changing the command below
xlim(model.unit.graphic.xlim); ylim(model.unit.graphic.ylim);
makeplotbold(figNr); drawnow;

end

% -------------------------------------------------------------------------
function [] = makeplotbold(figNr)

figure(figNr);
% increase line width (lw) and font size (fs)
lw = 2; fs = 20; fsL = 10;

set(gca,'FontSize',fs);
set(gca,'LineWidth',lw);
set(get(gca,'title'),'Fontsize',fs);
set(get(gca,'xlabel'),'Fontsize',fs);
set(get(gca,'ylabel'),'Fontsize',fs);
set(legend,'Fontsize',fsL);

end

% -------------------------------------------------------------------------
function [] = check_if_filename_exists(filename)

% if ~isfile(['modelspecification/' filename '.m'])
% if ~isfile(['modelspecification/' filename '.m'])
% if ~isfile([filename '.m'])
%     fprintf('\n\n --> No file %s.m found in folder ''modelspecification''! \n\n',filename);
%     error('--> Ending here')
% end

end
