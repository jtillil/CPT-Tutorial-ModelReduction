%%% Version: January 24th, 2020
%%%
%%% call by: redmodel  =  model_order_reduction(model,seqofstates,relerrTOL)
%%%
%%% This function returns the classification of the state variables as
%%% environmental, negligible, quasi steady state, mass conserved or
%%% dynamical state variable
%%%
%%% Input:  model               structure specifying the model
%%%         seqofstates         sequence, in which state variables are
%%%                             tested for model order reduction
%%%         relerrTOL           user defined relative error threshold
%%%
%%% Output: redmodel            structure specifying the reduced order model
%%%
%%% Citation:
%%% 
%%% Knoechel, Kloft and Huisinga, "Sensitivity based input-response index to 
%%% analyse and reduce large-scale signalling networks"
%%% PLOS Comp. Biology, 2020 (under review)
%%% 
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function redmodel = model_order_reduction(model,seqofstates,relerrTOL)

tic
% indexing
I = model.I;

% reference time span and output
t_ref = model.t_ref;
Y_ref = model.Y_ref;

%%% some ODE solver option and right hand side of ODE
options.NonNegative = 1:model.I.nstates;
odefun = @(t,X,model)  model.odefun(t,X,model.par,model);
% jacobian, if provided
jacfun = model.jacfun;
if ~isempty(jacfun)
    jacfun = @(t,X,model)  model.jacfun(t,X,model.par,model);
end

%%% error function to compute difference between full and reduced model
%%% (here: normalized L2-error between full and reduced output)
errfun = @(Y_ref,Y_red) sqrt( trapz(t_ref,(Y_ref-Y_red).^2) ) / sqrt( trapz(t_ref,Y_ref.^2) );


%%% =======================================================================
%%% 1st step: exploit negligible and environmental state variables
%%% =======================================================================

fprintf('\n 1st step - eliminate state variables as negligible or environmental');
fprintf('\n          - tested state variables: ');

% initialise reduced order model structure
redmodel = model;
redmodel.X0_red = model.X0;

for k = seqofstates

    fprintf('%s, ',I.nmstate{k});
    rI = redmodel.I;
    
    %%% first tentatively consider the kth state to belong to the environment
    %%%
    tmodel_env = redmodel;  
    tmodel_env.I.dyn  = setdiff(rI.dyn,k); % remove k from dyn states
    tmodel_env.I.env  = union(rI.env,k);   % add k to env states
    
    % compute solution of tentatively reduced model and corresponding output
    if ~isempty(jacfun)
        options.Jacobian = @(t,X) jacfun(t,X,tmodel_env); 
    end
    [~,tX_red_env] = ode15s(@(t,X) odefun(t,X,tmodel_env),t_ref, tmodel_env.X0_red, options);
    tY_red_env = tX_red_env(:,I.output);
    
    %%% secondly, tentatively consider the kth state to be negligible
    %%%
    tmodel_neg = redmodel;  
    tmodel_neg.I.dyn = setdiff(rI.dyn,k); % remove k from dyn states
    tmodel_neg.I.pneg = union(rI.pneg,k);   % add k to pneg states
    tmodel_neg.X0_red(k) = 0; % set negligible state to zero
    
    % compute solution of tentatively reduced model and corresponding output
    if ~isempty(jacfun)
        options.Jacobian = @(t,X) jacfun(t,X,tmodel_neg);
    end
    [~,tX_red_neg] = ode15s(@(t,X) odefun(t,X,tmodel_neg), t_ref, tmodel_neg.X0_red, options);
    tY_red_neg = tX_red_neg(:,I.output);
    
    %%% determine reduced model dependend of error criterion 
    %%%
    if errfun(Y_ref,tY_red_neg) < relerrTOL
        if errfun(Y_ref,tY_red_env) < errfun(Y_ref,tY_red_neg)
            redmodel = tmodel_env;
            redmodel.X_red = tX_red_env;
        else
            redmodel = tmodel_neg;
            redmodel.X_red = tX_red_neg;
        end
    elseif errfun(Y_ref,tY_red_env) < relerrTOL
        redmodel = tmodel_env;
        redmodel.X_red = tX_red_env;
    end
    
end

% remove identified env and neg states from sequence of states
rI = redmodel.I;
seqofstates = setdiff(seqofstates,[rI.env,rI.pneg],'stable');


%%% =======================================================================
%%% 2n step: exploit quasi steady state approximation
%%% =======================================================================

fprintf('\n 2nd step - eliminate state variables via quasi-steady state approximation');
fprintf('\n          - tested state variables: ');

% mass matrix for the DAE solver
redmodel.M = eye(I.nstates);
options.NonNegative = [];

% If some of the tentatively reduced models results in numerical problem,
% e.g., intergration tolerance could not be meet etc, then such a reduced 
% model is not accepted. In such a case, the approximation error could
% anyway not be computed (due to the numerical issue). 
lastwarn('allfine'); warning off;
for k = seqofstates

    fprintf('%s, ',I.nmstate{k});
    rI = redmodel.I;

    % tentatively consider kth state in pss
    tmodel_pss = redmodel;
    tmodel_pss.I.dyn  = setdiff(rI.dyn,k); % remove k from dyn states
    tmodel_pss.I.pss  = union(rI.pss,k);   % add k to pss states
    tmodel_pss.M(k,k) = 0; % and set corresponding entry in mass matrix to zero
    
    % check whether initial value of tentative con state is zero. This
    % causes problems with the DAE solver. Set to small positive values 
    % to resolve problem. Note: initial value of tentative con state is
    % 'not relevant', since the DAE solver will determine its value
    if tmodel_pss.X0_red(k) == 0 
       tmodel_pss.X0_red(k) = 1e-3;
    end
    
    %%% set DAE solver options
    %%%
    options.Mass = tmodel_pss.M;
    options.MassSingular = 'yes';
    options.InitialSlope = odefun(0,tmodel_pss.X0_red,tmodel_pss);
    
    if ~isempty(jacfun)
        options.Jacobian = @(t,X) jacfun(t,X,tmodel_pss);
    end
    
    % compute solution of tentatively reduced model and corresponding output
    try
        [~,tX_red_pss] = ode15s(@(t,X) odefun(t,X,tmodel_pss), t_ref, tmodel_pss.X0_red, options);
        tY_red_pss = tX_red_pss(:,I.output);
    catch
        lastwarn('error'); % numerical error occoured 
    end
        
    if strcmp(lastwarn,'allfine') % numerics worked out fine, ...     
        if errfun(Y_ref,tY_red_pss) < relerrTOL 
            redmodel = tmodel_pss;
            redmodel.X_red = tX_red_pss;
        end
    else
        % numerics had problems (this could be an error or a warning)
        % ignore tentatively reduced model, reset lastwarn and continue
        lastwarn('allfine');
    end
    
end
warning on; % turn all warnings on again

% remove identified pss states from sequence of states
rI = redmodel.I;
seqofstates = setdiff(seqofstates,rI.pss,'stable');


%%% =======================================================================
%%% 3rd step: exploit conservation laws
%%% =======================================================================

fprintf('\n 3rd step - eliminate state variables via conservation laws');
fprintf('\n          - tested state variables: ');

options.Mass = redmodel.M;
options.MassSingular = 'yes';
options.InitialSlope = odefun(0,redmodel.X0_red,redmodel);

% only test remaining dyn states that are part of some conservation law
seqofstates = intersect(seqofstates,I.conlaw,'stable');

% if there exists already a conserved state, then do not try to reduced
% further (since there is only one conservation law)
if ~isempty(I.con)
    seqofstates = [];
end

for k = seqofstates
  
    fprintf('%s, ',I.nmstate{k});
    rI = redmodel.I;

    % tentatively consider kth state eliminated via a conservation law
    tmodel_con = redmodel;
    tmodel_con.I.dyn  = setdiff(rI.dyn,k); % remove k from dyn states
    tmodel_con.I.con  = union(rI.con,k);   % add k to con states

    % solve DAE system
    options.Jacobian = [];
    [~,tX_red_con] = ode15s(@(t,X) odefun(t,X,tmodel_con), t_ref, tmodel_con.X0_red, options);
    % ... and account for conservation laws
    if ~isempty(rI.con)
        tX_red_con(:,rI.con) = max(redmodel.X0_conlaw - sum(tX_red_con(:,setdiff(rI.conlaw,rI.con)),2),0);
    end
    tY_red_con = tX_red_con(:,I.output);

    
    if errfun(Y_ref,tY_red_con) < relerrTOL
        redmodel = tmodel_con;
        redmodel.X_red = tX_red_con;
        
       break;
    end
    
end


%%% =======================================================================
%%% report about final reduced model
%%% =======================================================================

rI = redmodel.I;

fprintf('\n\n Reduced model consists of');
fprintf('\n %d dynamical state variable(s): ',length(rI.dyn)); 
fprintf('%s, ',I.nmstate{rI.dyn})
fprintf('\n %d environmental state variable(s): ',length(rI.env)); 
fprintf('%s, ',I.nmstate{rI.env})
fprintf('\n %d neglected state variable(s): ',length(rI.pneg)); 
fprintf('%s, ',I.nmstate{rI.pneg})
fprintf('\n %d state variable(s) approximated by their quasi-steady state: ',length(rI.pss)); 
fprintf('%s, ',I.nmstate{rI.pss})
fprintf('\n %d state variable(s) eliminated by conservation law(s): ',length(rI.con)); 
fprintf('%s, ',I.nmstate{rI.con})
fprintf('\n')

elapsedtime = toc; fprintf(' [elapsed time = %.1f]\n\n',elapsedtime);

%%% assign output and approximation error of reduced model
redmodel.t_red  = t_ref;
redmodel.Y_red  = redmodel.X_red(:,I.output);
redmodel.relerr = errfun(Y_ref,redmodel.Y_red);


end
