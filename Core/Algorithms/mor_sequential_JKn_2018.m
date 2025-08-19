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

function redmodel = mor_sequential_JKn_2018(model,seqofstates,relerrTOL,variability,var_obj_prctile)

if nargin < 4
    variability.variability = 0;
end

tic
% indexing
I = model.I;
errtype = "rel2NE";

% reference time span and output
t_ref = model.t_ref;
X_ref = model.X_ref;

%%% some ODE solver option and right hand side of ODE
% options.NonNegative = 1:model.I.nstates;
% odefun = @(t,X,model)  model.odefun(t,X,model.par,model);
% jacobian, if provided
% jacfun = model.jacfun;
% if ~isempty(jacfun)
    % jacfun = @(t,X,model)  model.jacfun(t,X,model.par,model);
% end

%%% error function to compute difference between full and reduced model
%%% (here: normalized L2-error between full and reduced output)
% errfun = @(Y_ref,Y_red) sqrt( trapz(t_ref,(Y_ref-Y_red).^2) ) / sqrt( trapz(t_ref,Y_ref.^2) );


%%% =======================================================================
%%% 1st step: exploit negligible and environmental state variables
%%% =======================================================================

fprintf('\n 1st step - eliminate state variables as negligible or environmental');
fprintf('\n          - tested state variables: ');

% initialise reduced order model structure
redmodel = model;
if ~isfield(redmodel, "redconfig")
    config = repmat("dyn", [1 model.I.nstates]);
else
    config = redmodel.redconfig;
end

for k = seqofstates(ismember(seqofstates, redmodel.I.dyn))

    fprintf('%s',I.nmstate{k});
    
    %%% first tentatively consider the kth state to belong to the environment
    %%%
    config_env = config;
    config_env(k) = "env";
    
    % compute solution of tentatively reduced model and corresponding output
    if ~variability.variability
        [obj_env, ~, ~, Xred_env, ~] = objfun_simple(model, config_env, errtype);
        err_red_env = obj_env.errout;
    else
        obj_env = objfun_vectorized_simple(model, config_env, errtype, variability, var_obj_prctile);
        err_red_env = max(obj_env(2), obj_env(5));
    end

    %%% secondly, tentatively consider the kth state to be negligible
    %%%
    config_pneg = config;
    config_pneg(k) = "pneg";
    
    % compute solution of tentatively reduced model and corresponding output
    if ~variability.variability
        [obj_pneg, ~, ~, Xred_pneg, ~] = objfun_simple(model, config_pneg, errtype);
        err_red_pneg = obj_pneg.errout;
    else
        obj_pneg = objfun_vectorized_simple(model, config_pneg, errtype, variability, var_obj_prctile);
        err_red_pneg = max(obj_pneg(2), obj_pneg(5));
    end
    
    %%% determine reduced model dependend of error criterion 
    %%%
    if err_red_pneg < relerrTOL
        if err_red_env < err_red_pneg
            redmodel.I.dyn = setdiff(redmodel.I.dyn, k);
            redmodel.I.env = union(redmodel.I.env, k);
            if ~variability.variability
                redmodel.X_red = Xred_env;
            end
            config = config_env;
            redmodel.redobj = obj_env;
            fprintf(' (env),')
        else
            redmodel.I.dyn = setdiff(redmodel.I.dyn, k);
            redmodel.I.pneg = union(redmodel.I.pneg, k);
            if ~variability.variability
                redmodel.X_red = Xred_pneg;
            end
            config = config_pneg;
            redmodel.redobj = obj_pneg;
            fprintf(' (pneg),')
        end
    elseif err_red_env < relerrTOL
        redmodel.I.dyn = setdiff(redmodel.I.dyn, k);
        redmodel.I.env = union(redmodel.I.env, k);
        if ~variability.variability
            redmodel.X_red = Xred_env;
        end
        config = config_env;
        redmodel.redobj = obj_env;
        fprintf(' (env),')
    else
        fprintf(' (dyn),')
    end
    
end

rI = redmodel.I;

fprintf('\n\n Reduced model after step 1 consists of');
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

fprintf('\n 2nd step - eliminate remaining dynamic state variables as environmental');
fprintf('\n          - tested state variables: ');

for k = seqofstates(ismember(seqofstates, redmodel.I.dyn))

    fprintf('%s',I.nmstate{k});
    
    %%% tentatively consider every dynamic state again to be environmental
    %%%
    config_env = config;
    config_env(k) = "env";
    
    % compute solution of tentatively reduced model and corresponding output
    if ~variability.variability
        [obj_env, ~, ~, Xred_env, ~] = objfun_simple(model, config_env, errtype);
        err_red_env = obj_env.errout;
    else
        obj_env = objfun_vectorized_simple(model, config_env, errtype, variability, var_obj_prctile);
        err_red_env = max(obj_env(2), obj_env(5));
    end
    
    %%% determine reduced model dependend of error criterion 
    %%%
    if err_red_env < relerrTOL
        redmodel.I.dyn = setdiff(redmodel.I.dyn, k);
        redmodel.I.env = union(redmodel.I.env, k);
        redmodel.I.pss = setdiff(redmodel.I.pss, k);
        redmodel.I.pneg = setdiff(redmodel.I.pneg, k);
        if ~variability.variability
            redmodel.X_red = Xred_env;
        end
        config = config_env;
        redmodel.redobj = obj_env;
        fprintf(' (env),')
    else
        fprintf(' (dyn),')
    end
    
end

rI = redmodel.I;

fprintf('\n\n Reduced model after step 2 consists of');
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

fprintf('\n 3rd step - eliminate remaining dynamic or environmental state variables as negligible');
fprintf('\n          - tested state variables: ');

for k = setdiff(seqofstates, redmodel.I.pneg, 'stable')

    fprintf('%s',I.nmstate{k});
    
    %%% tentatively consider every state again to be negligible
    %%%
    config_pneg = config;
    config_pneg(k) = "pneg";
    
    % compute solution of tentatively reduced model and corresponding output
    if ~variability.variability
        [obj_pneg, ~, ~, Xred_pneg, ~] = objfun_simple(model, config_pneg, errtype);
        err_red_pneg = obj_pneg.errout;
    else
        obj_pneg = objfun_vectorized_simple(model, config_pneg, errtype, variability, var_obj_prctile);
        err_red_pneg = max(obj_pneg(2), obj_pneg(5));
    end
    
    %%% determine reduced model dependend of error criterion 
    %%%
    if err_red_pneg < relerrTOL
        redmodel.I.dyn = setdiff(redmodel.I.dyn, k);
        redmodel.I.env = setdiff(redmodel.I.env, k);
        redmodel.I.pss = setdiff(redmodel.I.pss, k);
        redmodel.I.pneg = union(redmodel.I.pneg, k);
        if ~variability.variability
            redmodel.X_red = Xred_pneg;
        end
        config = config_pneg;
        redmodel.redobj = obj_pneg;
        fprintf(' (pneg),')
    else
        fprintf(' (dyn),')
    end
    
end

% remove identified env and neg states from sequence of states
seqofstates = setdiff(seqofstates, [redmodel.I.env, redmodel.I.pneg], 'stable');


%%% =======================================================================
%%% 2n step: exploit quasi steady state approximation
%%% =======================================================================
do_pss = 0;
if do_pss

fprintf('\n 2nd step - eliminate state variables via quasi-steady state approximation');
fprintf('\n          - tested state variables: ');

% If some of the tentatively reduced models results in numerical problem,
% e.g., intergration tolerance could not be meet etc, then such a reduced 
% model is not accepted. In such a case, the approximation error could
% anyway not be computed (due to the numerical issue). 
lastwarn('allfine'); warning off;
for k = seqofstates

    fprintf('%s, ',I.nmstate{k});
    rI = redmodel.I;

    % tentatively consider kth state in pss
    config_pss = config;
    config_pss(k) = "pss";
    
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
            % redmodel.X_red = tX_red_pss;
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

end

%%% =======================================================================
%%% 3rd step: exploit conservation laws
%%% =======================================================================
do_conlaw = 0;
if do_conlaw

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
        % redmodel.X_red = tX_red_con;
        
       break;
    end
    
end

end

%%% =======================================================================
%%% report about final reduced model
%%% =======================================================================

rI = redmodel.I;

fprintf('\n\n Final reduced model consists of');
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
% redmodel.t_red  = t_ref;
% redmodel.Y_red  = redmodel.X_red(:,I.output);
redmodel.redconfig = config;

redmodel.redobj = objfun_simple(redmodel, redmodel.redconfig, "rel2NE");

end
