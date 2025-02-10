%%% Version: July 09th, 2021
%%%
%%% call by: redmodel  =  model_order_reduction(model,seqofstates,relerrTOL)
%%%
%%% This function returns the classification of the state variables as
%%% environmental, negligible, quasi steady state, mass conserved or
%%% dynamical state variable. Reduction must respect error threshold for
%%% all parameter sets in the population.
%%%
%%% Input:  model               structure specifying the model, including
%%%                             different parameters and initial values
%%%                             giving the population
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
%%% Authors: Undine Falkenhagen, Jane Knoechel and Wilhelm Huisinga
%%%

function redmodel = model_order_reduction_robust(model,seqofstates,relerrTOL,varargin)
p = inputParser();
addParameter(p, 'share', 1, @isnumeric);
parse(p,varargin{:})
share = p.Results.share;
    
tic
% indexing
if ~isfield(model.I,'h_red')
    model.I.h_red=model.I.h;
end
I = model.I;

% reference time span, initial value, parameter values and output
t_ref  = model.t_ref;
if isfield(model,'X0_refs')
    X0_refs = model.X0_refs;    pars = model.pars;    Y_refs = model.Y_refs;
else
    X0_refs = model.X0_ref;     pars = model.par;     Y_refs = model.Y_ref;
end
h_red = I.h_red;
Y_refs_1 = Y_refs(:,1);

% reduced dose
reduced_dose = model.reduced_dose;

% scenario
scenario = model.scenario;

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
        tmodel_env_reduced_dose = tmodel_env;
        tmodel_env_reduced_dose.u_ref = tmodel_env_reduced_dose.u_ref/4;
        tmodel_env_reduced_dose.input_doses = tmodel_env_reduced_dose.input_doses/4;
        tmodel_env_reduced_dose.X0_refs = tmodel_env_reduced_dose.X0s + tmodel_env_reduced_dose.u_ref';        
        
        nind = size(X0_refs,1);
        relerr_env = zeros(nind,1);
        parfor ind=1:nind
            if reduced_dose(ind)
                tmodel_env_i=tmodel_env_reduced_dose;
            else
                tmodel_env_i=tmodel_env;
            end
            X0_refs(ind,:)=tmodel_env_i.X0_refs(ind,:);
            try
                % compute solution of tentatively reduced model and corresponding output
                [~,tX_red_env,te_red_env] = tmodel_env_i.simulateODE(t_ref,X0_refs(ind,:),pars(ind,:), tmodel_env_i);
                tY_red_env = h_red(t_ref,tX_red_env,te_red_env); %#ok<PFBNS>
                        
                if scenario=="in_vivo_warfarin"
                    tY_red_env=tY_red_env*Y_refs_1(ind);
                end
            
                relerr_env(ind) = errfun(t_ref,Y_refs(ind,:)',tY_red_env); %#ok<PFBNS>
            catch
                % numerics had problems (this could be an error or a warning)
                % ignore tentatively reduced model and continue
                fprintf('\b\b (error env), ');
                relerr_env(ind) = Inf;
    
            end
        end
        
        %%% secondly, tentatively consider the kth state to be negligible
        %%%
        tmodel_neg = redmodel;
        tmodel_neg.I.dyn = setdiff(rI.dyn,k,'stable'); % remove k from dyn states
        tmodel_neg.I.neg = [rI.neg,k];   % add k to neg states
        tmodel_neg_reduced_dose = tmodel_neg;
        tmodel_neg_reduced_dose.u_ref = tmodel_neg_reduced_dose.u_ref/4;
        tmodel_neg_reduced_dose.input_doses = tmodel_neg_reduced_dose.input_doses/4;
        tmodel_neg_reduced_dose.X0_refs = tmodel_neg_reduced_dose.X0s + tmodel_neg_reduced_dose.u_ref'; 
        relerr_neg = zeros(size(X0_refs,1),1);

        parfor ind=1:nind
            if reduced_dose(ind)
                tmodel_neg_i=tmodel_neg_reduced_dose;
            else
                tmodel_neg_i=tmodel_neg;
            end
            X0_refs(ind,:)=tmodel_neg_i.X0_refs(ind,:);
            % compute solution of tentatively reduced model and corresponding output
            try
                [~,tX_red_neg,te_red_neg] = tmodel_neg_i.simulateODE(t_ref,X0_refs(ind,:),pars(ind,:), tmodel_neg_i);
                tY_red_neg = h_red(t_ref,tX_red_neg,te_red_neg); %#ok<PFBNS>

                if scenario=="in_vivo_warfarin"
                    tY_red_neg=tY_red_neg*Y_refs_1(ind);
                end

                relerr_neg(ind) = errfun(t_ref,Y_refs(ind,:)',tY_red_neg); %#ok<PFBNS>
            catch
                % numerics had problems (this could be an error or a warning)
                % ignore tentatively reduced model and continue
                fprintf('\b\b (error neg), ');
                relerr_neg(ind) = Inf;
            end
        end

        %%% determine reduced model dependent on error criterion
        %%%
        if sum(isfinite(relerr_neg))>=share && sum(relerr_neg < relerrTOL)/size(X0_refs,1)>=share
            if all(isfinite(relerr_env)) && any(relerr_env < relerr_neg)
                redmodel = tmodel_env;
                fprintf('\b\b (env), ');
            else
                redmodel = tmodel_neg;
                fprintf('\b\b (neg), ');
            end
        elseif sum(isfinite(relerr_env))>=share && sum(relerr_env < relerrTOL)/size(X0_refs,1)>=share
            redmodel = tmodel_env;
            fprintf('\b\b (env), ');
        end
        
    end
    warning on;
    
    % remove identified env and neg states from sequence of states
    rI = redmodel.I;
    seqofstates = setdiff(seqofstates,[rI.env,rI.neg],'stable');
    classificationstatus(redmodel.I);
    
    
%     %%% =======================================================================
%     %%% 2n step: exploit quasi steady state approximation
%     %%% =======================================================================
%     
%     fprintf('\n 2nd step - eliminate state variables via quasi-steady state approximation');
%     fprintf('\n          > tested state variables: ');
%     
%     % If some of the tentatively reduced models results in numerical problem,
%     % e.g., intergration tolerance could not be meet etc, then such a reduced
%     % model is not accepted. In such a case, the approximation error could
%     % anyway not be computed (due to the numerical issue).
%     warning off;
%     
%         
%     for k = seqofstates
%         
%         fprintf('%s, ',I.nmstate{k});
%         rI = redmodel.I;
%         
%         % tentatively consider kth state in qss
%         tmodel_qss = redmodel;
%         tmodel_qss.I.dyn  = setdiff(rI.dyn,k,'stable'); % remove k from dyn states
%         tmodel_qss.I.qss  = [rI.qss,k];   % add k to qss states
%         tmodel_qss_reduced_dose = tmodel_qss;
%         tmodel_qss_reduced_dose.u_ref = tmodel_qss_reduced_dose.u_ref/4;
%         tmodel_qss_reduced_dose.input_doses = tmodel_qss_reduced_dose.input_doses/4;
%         tmodel_qss_reduced_dose.X0_refs = tmodel_qss_reduced_dose.X0s + tmodel_qss_reduced_dose.u_ref'; 
%         relerr_qss = zeros(size(X0_refs,1),1);
%         % compute solution of tentatively reduced model and corresponding output
%         parfor ind=1:nind
%             if reduced_dose(ind)
%                 tmodel_qss_i=tmodel_qss_reduced_dose;
%             else
%                 tmodel_qss_i=tmodel_qss;
%             end
%             X0_refs(ind,:)=tmodel_qss_i.X0_refs(ind,:);
%             try
%                 [~,tX_red_qss,te_red_qss] = tmodel_qss_i.simulateODE(t_ref,X0_refs(ind,:),pars(ind,:), tmodel_qss_i);
%                 tY_red_qss = h_red(t_ref,tX_red_qss,te_red_qss); %#ok<PFBNS>
% 
%                 if scenario=="in_vivo_warfarin"
%                     tY_red_qss=tY_red_qss*Y_refs_1(ind);
%                 end
% 
%                 % numerics seems to be fine, so compute the approximation error
%                 relerr_qss(ind) = errfun(t_ref,Y_refs(ind,:)',tY_red_qss); %#ok<PFBNS>
% 
%             catch
%                 % numerics had problems (this could be an error or a warning)
%                 % ignore tentatively reduced model and continue
%                 fprintf('\b\b (error qss), ');
%                 relerr_qss(ind) = Inf;
% 
%             end
%         end
%         if sum(relerr_qss < relerrTOL)/length(relerr_qss)>=share
%             redmodel = tmodel_qss;
%             fprintf('\b\b (qss), ');
%         end
%     end
%     warning on; % turn all warnings on again
%     
%     
%     % remove identified qss states from sequence of states
%     rI = redmodel.I;
%     seqofstates = setdiff(seqofstates,rI.qss,'stable');
%     classificationstatus(redmodel.I);
%  

%%% =======================================================================
%%% Interim step: try to neglect environmental state variables
%%% =======================================================================

interimseqofstates = [rI.env rI.qss];
fprintf('\n Interim step - try to neglect env state variables'); %and qss 
fprintf('\n              - tested state variables: ');

warning off;
for k = interimseqofstates
    
    fprintf('%s ',I.nmstate{k});
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
    relerr_neg = zeros(size(X0_refs,1),1);
          
    % compute solution of tentatively reduced model and corresponding output
    for ind=1:nind
        if reduced_dose(ind)
            tmodel_neg.u_ref = tmodel_neg.u_ref/4;
            tmodel_neg.input_doses = tmodel_neg.input_doses/4;
            tmodel_neg.X0_refs = tmodel_neg.X0s + tmodel_neg.u_ref';
            X0_refs(ind,I.Awarf)=1;
        end
        try
            [~,tX_red_neg,te_red_neg] = model.simulateODE(t_ref,X0_refs(ind,:),pars(ind,:), tmodel_neg);
            tY_red_neg = I.h_red(t_ref,tX_red_neg,te_red_neg);

            if redmodel.scenario=="in_vivo_warfarin"
                tY_red_neg=tY_red_neg*redmodel.Y_refs(ind,1);
            end


            % numerics seems to be fine, so compute the approximation error
            relerr_neg(ind) = errfun(t_ref,Y_refs(ind,1:length(tY_red_neg))',tY_red_neg);

        catch
            % numerics had problems (this could be an error or a warning)
            % ignore tentatively reduced model and continue
            fprintf('(error), ');
            relerr_neg(ind) = Inf;

        end
        if reduced_dose(ind)
            tmodel_neg.u_ref = tmodel_neg.u_ref*4;
            tmodel_neg.input_doses = tmodel_neg.input_doses*4;
            tmodel_neg.X0_refs = tmodel_neg.X0s + tmodel_neg.u_ref';
        end
        if sum(relerr_neg>relerrTOL)>nind*(1-share)
            break
        end
    end
    if sum(relerr_neg < relerrTOL)/length(relerr_neg)>=share
        redmodel = tmodel_neg;
        fprintf('(neg), ');
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
    [~,X_red] = model.simulateODE(t_ref,X0_refs,pars, redmodel);
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
                    [~,tX_red_con] = model.simulateODE(t_ref,X0_refs,pars, tmodel_con);
                    tY_red_con = I.h_red(t_ref,tX_red_con);
                    
                    % numerics seems to be fine, so compute the approximation error
                    if errfun(t_ref,Y_refs,tY_red_con) < relerrTOL
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
% fprintf('\n %d state variable(s) approximated by their quasi-steady state: ',length(rI.qss));
% fprintf('%s, ',I.nmstate{rI.qss})
fprintf('\n %d state variable(s) eliminated by conservation law(s): ',length(rI.con));
fprintf('%s, ',I.nmstate{rI.con})
fprintf('\n')

elapsedtime = toc; fprintf(' [model order reduction = %.1f sec]\n\n',elapsedtime);

%%% simulate final reduced model and assign output & approximation error
tic;
Y_reds = zeros(nind,size(Y_refs,2));
X_reds = zeros(nind,1,size(X0_refs,2));
for ind=1:nind
    if reduced_dose(ind)
        redmodel.u_ref = redmodel.u_ref/4;
        redmodel.input_doses = redmodel.input_doses/4;
        redmodel.X0_refs = redmodel.X0s + redmodel.u_ref';
        X0_refs(ind,I.Awarf)=1;
    end
    [~,X_red,te] = model.simulateODE(t_ref,X0_refs(ind,:),pars(ind,:), redmodel);
    Y_reds(ind,:) = rI.h_red(t_ref,X_red,te);
    X_reds(ind,1:size(X_red,1),:) = X_red;
    if redmodel.scenario=="in_vivo_warfarin"
        Y_reds(ind,:)=Y_reds(ind,:)*Y_refs(ind,1);    
    end
    relerr(ind) = errfun(t_ref,Y_refs(ind,1:length(Y_reds(ind,:))),Y_reds(ind,:));
    if reduced_dose(ind)
        redmodel.u_ref = redmodel.u_ref*4;
        redmodel.input_doses = redmodel.input_doses*4;
        redmodel.X0_refs = redmodel.X0s + redmodel.u_ref';
    end
end
elapsedtime = toc; fprintf(' [simulation of reduced model = %.1f sec]\n\n',elapsedtime);
redmodel.Y_reds = Y_reds;
redmodel.X_reds  = X_reds;
redmodel.errfun = errfun;
redmodel.relerr = max(relerr);
end

%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% LOCAL SUB-ROUTINES
%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -------------------------------------------------------------------------
function [] = classificationstatus(I)

fprintf('\n          > current classification: dyn=%d, env=%d, neg=%d, qss=%d, con=%d ',...
        length(I.dyn),length(I.env),length(I.neg),length(I.qss),length(I.con));

end
