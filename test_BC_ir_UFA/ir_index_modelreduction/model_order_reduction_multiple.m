%%% Version: February 09th, 2020
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

function redmodel = model_order_reduction_multiple(model,seqofstates,relerrTOL)

tic
% indexing
if ~isfield(model.I,'h_red')
    model.I.h_red=model.I.h;
end
I = model.I;

% reference time span, initial value, parameter values and output
t_ref  = model.t_ref;
X0_ref = model.X0_ref;
par    = model.par;
Y_ref  = model.Y_ref(:);


% error function
errfun = model.errfun;

% initialise reduced order model structure
redmodel = model;


%%% perform model order reduction multiple times!
for r = 1:model.nrepeats
    
    %%% =======================================================================
    %%% 1st step: exploit negligible and environmental state variables
    %%% =======================================================================
    
    fprintf('\n 1st step - eliminate state variables as negligible or environmental');
    fprintf('\n          > tested state variables: ');
    
    
    warning off;
    for k = seqofstates
        
        fprintf('%s, ',I.nmstate{k});
        rI = redmodel.I;
        
        %%% first tentatively consider the kth state to belong to the environment
        %%%
        tmodel_env = redmodel;
        tmodel_env.I.dyn  = setdiff(rI.dyn,k,'stable'); % remove k from dyn states
        tmodel_env.I.env  = [rI.env,k];   % add k to env states
        
        try
            
            % compute solution of tentatively reduced model and corresponding output
            [~,tX_red_env,te_red_env] = model.simulateODE(t_ref,X0_ref,par, tmodel_env);
            %tY_red_env = tX_red_env(:,I.output);
            tY_red_env = I.h_red(t_ref,tX_red_env,te_red_env);
            
            if redmodel.scenario=="in_vivo_warfarin"
                tY_red_env=tY_red_env*redmodel.Y_ref(1);
            end
            
            relerr_env = errfun(t_ref,Y_ref,tY_red_env);
        catch
            % numerics had problems (this could be an error or a warning)
            % ignore tentatively reduced model and continue
            fprintf('\b\b (error env), ');
            relerr_env = Inf;
            
        end
        
        %%% secondly, tentatively consider the kth state to be negligible
        %%%
        tmodel_neg = redmodel;
        tmodel_neg.I.dyn = setdiff(rI.dyn,k,'stable'); % remove k from dyn states
        tmodel_neg.I.neg = [rI.neg,k];   % add k to neg states
        
        
        % compute solution of tentatively reduced model and corresponding output
        try
            [~,tX_red_neg,te_red_neg] = model.simulateODE(t_ref,X0_ref,par, tmodel_neg);
            %tY_red_neg = tX_red_neg(:,I.output);
            tY_red_neg = I.h_red(t_ref,tX_red_neg,te_red_neg);
            
            if redmodel.scenario=="in_vivo_warfarin"
                tY_red_neg=tY_red_neg*redmodel.Y_ref(1);
            end
            
            relerr_neg = errfun(t_ref,Y_ref,tY_red_neg);
        catch
            % numerics had problems (this could be an error or a warning)
            % ignore tentatively reduced model and continue
            fprintf('\b\b (error neg), ');
            relerr_neg = Inf;
        end
        
        %%% determine reduced model dependend of error criterion
        %%%
        if isfinite(relerr_neg) && (relerr_neg < relerrTOL)
            if isfinite(relerr_env) && (relerr_env < relerr_neg)
                redmodel = tmodel_env;
                fprintf('\b\b (env), ');
            else
                redmodel = tmodel_neg;
                fprintf('\b\b (neg), ');
            end
        elseif isfinite(relerr_env) && (relerr_env < relerrTOL)
            redmodel = tmodel_env;
            fprintf('\b\b (env), ');
        end
        
        
    end
    warning on;
    
    % remove identified env and neg states from sequence of states
    rI = redmodel.I;
    seqofstates = setdiff(seqofstates,[rI.env,rI.neg],'stable');
    classificationstatus(redmodel.I);
    
    
    %%% =======================================================================
    %%% 2n step: exploit quasi steady state approximation
    %%% =======================================================================
    
    fprintf('\n 2nd step - eliminate state variables via quasi-steady state approximation');
    fprintf('\n          > tested state variables: ');
    
    % If some of the tentatively reduced models results in numerical problem,
    % e.g., intergration tolerance could not be meet etc, then such a reduced
    % model is not accepted. In such a case, the approximation error could
    % anyway not be computed (due to the numerical issue).
    warning off;
    
%     %% qss symbolic
%     X_sym_var=cell2sym(I.nmstate(:));
%     X_sym=X_sym_var;
%     %X_sym(redmodel.I.env)=vpa(redmodel.X0(redmodel.I.env));
%     X_sym(redmodel.I.neg)=0;
%     digits(4)
%     redmodel.symode_red=vpa(simplify(redmodel.odefun(0,X_sym,redmodel.par,redmodel)));
    %%
        
    for k = seqofstates
        
        fprintf('%s, ',I.nmstate{k});
        rI = redmodel.I;
        
        % tentatively consider kth state in qss
        tmodel_qss = redmodel;
        tmodel_qss.I.dyn  = setdiff(rI.dyn,k,'stable'); % remove k from dyn states
        tmodel_qss.I.qss  = [rI.qss,k];   % add k to qss states
        %% qss symbolic
        % compute solution of tentatively reduced model and corresponding output
        try
            [~,tX_red_qss,te_red_qss] = model.simulateODE(t_ref,X0_ref,par, tmodel_qss);
            tY_red_qss = I.h_red(t_ref,tX_red_qss,te_red_qss);
            
            if redmodel.scenario=="in_vivo_warfarin"
                tY_red_qss=tY_red_qss*redmodel.Y_ref(1);
            end
            
            % numerics seems to be fine, so compute the approximation error
            if errfun(t_ref,Y_ref,tY_red_qss) < relerrTOL
                redmodel = tmodel_qss;
                fprintf('\b\b (qss), ');
            end
            
        catch
            % numerics had problems (this could be an error or a warning)
            % ignore tentatively reduced model and continue
            fprintf('\b\b (error qss), ');
            
        end
    end
    warning on; % turn all warnings on again
    
    
    % remove identified qss states from sequence of states
    rI = redmodel.I;
    seqofstates = setdiff(seqofstates,rI.qss,'stable');
    classificationstatus(redmodel.I);
 

%%% =======================================================================
%%% Interim step: try to neglect environmental state variables
%%% =======================================================================

interimseqofstates = [rI.env rI.qss];
fprintf('\n Interim step - try to neglect env and qss state variables');
fprintf('\n              - tested state variables: ');

warning off;
for k = interimseqofstates
    
    fprintf('%s, ',I.nmstate{k});
    rI = redmodel.I;
    
    %%% tentatively consider the kth state to be negligible
    %%%
    tmodel_neg = redmodel;
    if ismember(k,rI.env)
        tmodel_neg.I.env = setdiff(rI.env,k,'stable'); % remove k from env states
    else
        tmodel_neg.I.qss = setdiff(rI.qss,k,'stable'); % remove k from qss states
    end
    tmodel_neg.I.neg = [rI.neg,k];   % add k to neg states
            %% qss symbolic
    % compute solution of tentatively reduced model and corresponding output
    try
        [~,tX_red_neg,te_red_neg] = model.simulateODE(t_ref,X0_ref,par, tmodel_neg);
        %tY_red_neg = tX_red_neg(:,I.output);
        tY_red_neg = I.h_red(t_ref,tX_red_neg,te_red_neg);
        
        if redmodel.scenario=="in_vivo_warfarin"
            tY_red_neg=tY_red_neg*redmodel.Y_ref(1);
        end
            
        
        % numerics seems to be fine, so compute the approximation error
        if errfun(t_ref,Y_ref,tY_red_neg) < relerrTOL
            redmodel = tmodel_neg;
            fprintf('\b\b (neg), ');
        end
        
    catch
        % numerics had problems (this could be an error or a warning)
        % ignore tentatively reduced model and continue
        fprintf('\b\b (error), ');
        
    end
end
warning on;
classificationstatus(redmodel.I);
   
end % of performing steps 1-2 and interim step multiple times


%%% =======================================================================
%%% 3rd step: exploit conservation laws
%%% =======================================================================

fprintf('\n 3rd step - eliminate state variables via conservation laws');

% only perform this step, if number (n) of conservation laws > 0
rL = redmodel.L;
if rL.nconlaws > 0
    
    %%% check, whether conlaws still approximately hold for the reduced
    %%% model
    [~,X_red] = model.simulateODE(t_ref,X0_ref,par, redmodel);
    validconlaws = [];
    redvalconlaw = NaN(rL.nconlaws,1);
    for c = 1:rL.nconlaws
        
        % determine value of conlaw in reduced model
        redvalconlaw(c) = sum(X_red(1,rL.statesofconlaw{c}),2);
        
        % report about relative change in value of the conlaw
        reldiffconlaw = (redvalconlaw(c)-rL.valconlaw(c))/rL.valconlaw(c);
        fprintf('\n Conlaw(%d) rel. change by %+.2f %%',c,reldiffconlaw*100);
        conlaw_env = sum(any(ismember(rI.env,rL.statesofconlaw{c})));
        conlaw_neg = sum(any(ismember(rI.neg,rL.statesofconlaw{c})));
        conlaw_qss = sum(any(ismember(rI.qss,rL.statesofconlaw{c})));
        conlaw_con = sum(any(ismember(rI.con,rL.statesofconlaw{c})));
        fprintf('[%d env, %d neg, %d qss, %d con states]',conlaw_env,conlaw_neg,conlaw_qss,conlaw_con);
        
        
        % still a valid conlaw for the reduced model? Use same threshold as
        % for classifying the state variables
        if max( abs( rL.valconlaw(c) - sum(X_red(:,rL.statesofconlaw{c}),2) ) )/rL.valconlaw(c) < relerrTOL
            validconlaws = [validconlaws,c];
        end
        
        
    end
    
    fprintf('\n Valid conlaws for reduced model: ');
    statesofallvalidconlaws = unique([rL.statesofconlaw{validconlaws}]);
    
    redmodel.L.nvalidconlaws = length(validconlaws);
    redmodel.L.validconlaws = validconlaws;
    redmodel.L.statesofallvalidconlaws = statesofallvalidconlaws;
    
    
    % only test remaining dyn states that are part of some conservation law
    seqofstates = intersect(seqofstates,statesofallvalidconlaws,'stable');
    
    % report about still valid conlaws of the reduced model
    fprintf(' [%d]',validconlaws);
    fprintf('\n          - tested state variables: ');
    
    warning off;
    for k = seqofstates
        
        fprintf('%s ',I.nmstate{k});
        rI = redmodel.I; rL = redmodel.L;
        
        % to avoid conflicts between different conlaws, allow only conlaws
        % to be used that do not contain any state already previously
        % eliminated by a conlaw, i.e., not in rL.con
        possibleconlaws = validconlaws;
        for c = validconlaws
            usedconlaws = rL.con2conlaw;
            statesofusedconlaws = [rL.statesofconlaw{usedconlaws}];
            
            % check, whether cth conlaw has some states in common with
            % the states of already used conlaws
            if any(ismember(rL.statesofconlaw{c},statesofusedconlaws))
                possibleconlaws = setdiff(possibleconlaws,c,'stable');
            end
        end
        
        if ~ismember(k,[rL.statesofconlaw{possibleconlaws}])
            fprintf('(skipped), ');
            continue; % skip the remaining part and continue with next state
        end
        
        for c = possibleconlaws
            
            % check whether kth state belongs to cth conlaw
            if ismember(k, rL.statesofconlaw{c})
                
                fprintf('[%d]',c);
                % tentatively consider kth state eliminated via cth conservation law
                tmodel_con = redmodel;
                tmodel_con.I.dyn  = setdiff(rI.dyn,k,'stable'); % remove k from dyn states
                tmodel_con.I.con  = [rI.con,k];   % add k to con states
                tmodel_con.L.con2conlaw = [rL.con2conlaw,c];
                
                % determine remaining (rem) states of the cth conlaw
                remstates = setdiff(rL.statesofconlaw{c},k,'stable');
                tmodel_con.L.remstatesofconlaw = [rL.remstatesofconlaw,remstates];
                
                % compute solution of tentatively reduced model and corresponding output
                try
                    [~,tX_red_con] = model.simulateODE(t_ref,X0_ref,par, tmodel_con);
                    tY_red_con = I.h_red(t_ref,tX_red_con);
                    
                    % numerics seems to be fine, so compute the approximation error
                    if errfun(t_ref,Y_ref,tY_red_con) < relerrTOL
                        redmodel = tmodel_con;
                        fprintf(' (con)');
                        break; % exit 'for c = unused conlaws' loop and continue with next k in seqofstates
                    end
                    
                catch
                    % numerics had problems (this could be an error or a warning)
                    % ignore tentatively reduced model and continue
                    fprintf('\b error]');
                    
                end
            end
        end
        fprintf(', ')
    end
    warning on; % turn all warnings on again
end % rL.nconlaws > 0
classificationstatus(redmodel.I);


%%% =======================================================================
%%% report about final reduced model
%%% =======================================================================

rI = redmodel.I;

fprintf('\n\n Reduced model consists of');
fprintf('\n %d dynamical state variable(s): ',length(rI.dyn));
fprintf('%s, ',I.nmstate{rI.dyn})
fprintf('\n %d environmental state variable(s): ',length(rI.env));
fprintf('%s, ',I.nmstate{rI.env})
fprintf('\n %d neglected state variable(s): ',length(rI.neg));
fprintf('%s, ',I.nmstate{rI.neg})
fprintf('\n %d state variable(s) approximated by their quasi-steady state: ',length(rI.qss));
fprintf('%s, ',I.nmstate{rI.qss})
fprintf('\n %d state variable(s) eliminated by conservation law(s): ',length(rI.con));
fprintf('%s, ',I.nmstate{rI.con})
fprintf('\n')

elapsedtime = toc; fprintf(' [model order reduction = %.1f sec]\n\n',elapsedtime);

%%% simulate final reduced model and assign output & approximation error
tic;
[t_red,X_red,te] = model.simulateODE(t_ref,X0_ref,par, redmodel);
elapsedtime = toc; fprintf(' [simulation of reduced model = %.1f sec]\n\n',elapsedtime);

Y_red = rI.h_red(t_ref,X_red,te);
redmodel.t_red  = t_red;
redmodel.X_red  = X_red;
redmodel.Y_red  = I.h_red(t_ref,X_red,te);
redmodel.errfun = errfun;
redmodel.relerr = errfun(t_ref,Y_ref,Y_red);

end

%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% LOCAL SUB-ROUTINES
%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -------------------------------------------------------------------------
function [] = classificationstatus(I)

fprintf('\n          > current classification: dyn=%d, env=%d, neg=%d, qss=%d, con=%d ',...
        length(I.dyn),length(I.env),length(I.neg),length(I.qss),length(I.con));

end
