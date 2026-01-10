function redmodel = mor_backwards_UFa_2023(model,seqofstates,relerrTOL,variability,var_obj_prctile)

tic

% indexing
I = model.I;
errtype = "rel2NE";
seqofstates = flip(seqofstates);

% reference time span and output
t_ref = model.t_ref;
X_ref = model.X_ref;

%%% some ODE solver options and right hand side of ODE
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

fprintf('\n 1st step - make states dynamical again from negligible or environmental');
fprintf('\n          - tested state variables: ');

% initialise reduced order model structure
redmodel = model;
config = redmodel.redconfig;
% redmodel.X0_red = model.X0;

for k = setdiff(seqofstates, redmodel.I.dyn, 'stable')

    fprintf('%s, ',I.nmstate{k});
    
    %%% tentatively consider the kth state to be dynamic
    %%%
    config_dyn = config;
    config_dyn(k) = "dyn";
    
    % compute solution of tentatively reduced model and corresponding output
    obj_dyn = objfun_vectorized_simple(model, config_dyn, errtype, variability, var_obj_prctile);
    err_red_dyn = max(obj_dyn(2), obj_dyn(5));
    fprintf(['; ' char(num2str(obj_dyn(2))) '; ' char(num2str(obj_dyn(5)))]);
    
    %%% determine reduced model dependend of error criterion 
    %%%
    if err_red_dyn > relerrTOL
        redmodel.I.dyn = union(redmodel.I.dyn, k);
        redmodel.I.env = setdiff(redmodel.I.env, k);
        redmodel.I.pneg = setdiff(redmodel.I.pneg, k);
        redmodel.I.pss = setdiff(redmodel.I.pss, k);
        % redmodel.X_red = Xred_dyn;
        config = config_dyn;
        redmodel.redobj = obj_dyn;
    else
        redmodel.I.dyn = union(redmodel.I.dyn, k);
        redmodel.I.env = setdiff(redmodel.I.env, k);
        redmodel.I.pneg = setdiff(redmodel.I.pneg, k);
        redmodel.I.pss = setdiff(redmodel.I.pss, k);
        % redmodel.X_red = Xred_dyn;
        config = config_dyn;
        redmodel.redobj = obj_dyn;
        break
    end
    
end

%%% =======================================================================
%%% report about final reduced model
%%% =======================================================================

rI = redmodel.I;

fprintf('\n\n Backtracked reduced model consists of');
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

end
