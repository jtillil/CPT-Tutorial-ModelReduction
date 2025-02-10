%%% Version: September 17th, 2021
%%%
%%% call by: redmodel  =  model_order_reduction_parameters(model,seqofstates,relerrTOL)
%%%
%%% This function returns the classification of the parameters as
%%% negligible, or normal parameters
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

function redmodel = model_order_reduction_parameters_robust(model,seqofpars,relerrTOL,varargin)
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
X0_refs = model.X0_refs;     
pars = model.pars;     
Y_refs = model.Y_refs;
Y_refs_1 = Y_refs(:,1);
reduced_dose = model.reduced_dose;
scenario = model.scenario;
h_red  = I.h_red;


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
    for k = seqofpars
        
        fprintf('%s, ',I.nmpar{k});
        rI = redmodel.I;
        
        %%% first tentatively consider the kth parameter to be negligible
        %%%
        tmodel_neg = redmodel;
        tmodel_neg.I.p = setdiff(rI.p,k,'stable');
        tmodel_neg.I.pneg = [rI.pneg,k];   % add k to neg states
        tmodel_neg_reduced_dose = tmodel_neg;
        tmodel_neg_reduced_dose.u_ref = tmodel_neg_reduced_dose.u_ref/4;
        tmodel_neg_reduced_dose.input_doses = tmodel_neg_reduced_dose.input_doses/4;
        tmodel_neg_reduced_dose.X0_refs = tmodel_neg_reduced_dose.X0s + tmodel_neg_reduced_dose.u_ref'; 
        
        % compute solution of tentatively reduced model and corresponding output
        nind = size(X0_refs,1);
        relerr_neg = zeros(nind,1);
        %fprintf(1,'(%4d/%4d)',0,nind);
        parfor ind=1:nind
            %fprintf(1,'\b\b\b\b\b\b\b\b\b\b%4d/%4d)',ind,nind);
            if reduced_dose(ind)
                tmodel_neg_i=tmodel_neg_reduced_dose;
            else
                tmodel_neg_i=tmodel_neg;
            end
            X0_refs(ind,:)=tmodel_neg_i.X0_refs(ind,:);
            try
                [~,tX_red_neg,te_red_neg] = tmodel_neg_i.simulateODE(t_ref,X0_refs(ind,:),pars(ind,:), tmodel_neg_i);
                tY_red_neg = h_red(t_ref,tX_red_neg,te_red_neg); %#ok<PFBNS>
                if scenario=="in_vivo_warfarin"
                    tY_red_neg=tY_red_neg*Y_refs_1(ind);
                end
                relerr_neg(ind) = errfun(t_ref,Y_refs(ind,:)',tY_red_neg); %#ok<PFBNS>
            catch
                % numerics had problems -> ignore reduced model and continue
                fprintf('\b\b (error neg), ');
                relerr_neg(ind) = Inf;
            end
%             if sum(relerr_neg>relerrTOL)>nind*(1-share)
%                 break
%             end
        end
        %fprintf('\b\b\b\b\b\b\b\b\b\b\b')
        
        
        if sum(isfinite(relerr_neg))>=share && sum(relerr_neg < relerrTOL)/length(relerr_neg)>=share
            redmodel = tmodel_neg;
            fprintf('\b\b (neg), ');
        else
            tmodel_inf = redmodel;
            tmodel_inf.I.p = setdiff(rI.p,k,'stable');
            tmodel_inf.I.pinf = [rI.pinf,k];
            tmodel_inf_reduced_dose = tmodel_inf;
            tmodel_inf_reduced_dose.u_ref = tmodel_inf_reduced_dose.u_ref/4;
            tmodel_inf_reduced_dose.input_doses = tmodel_inf_reduced_dose.input_doses/4;
            tmodel_inf_reduced_dose.X0_refs = tmodel_inf_reduced_dose.X0s + tmodel_inf_reduced_dose.u_ref'; 
            %fprintf(1,'(%4d/%4d)',0,nind);
            relerr_inf = zeros(nind,1);
            parfor ind=1:nind
                %fprintf(1,'\b\b\b\b\b\b\b\b\b\b%4d/%4d)',ind,nind);         
                if reduced_dose(ind)
                    tmodel_inf_i=tmodel_inf_reduced_dose;
                else
                    tmodel_inf_i=tmodel_inf;
                end
                X0_refs(ind,:)=tmodel_inf_i.X0_refs(ind,:);
                try
                    [~,tX_red_inf,te_red_inf] = tmodel_inf_i.simulateODE(t_ref,X0_refs(ind,:),pars(ind,:), tmodel_inf_i);
                    tY_red_inf = h_red(t_ref,tX_red_inf,te_red_inf); %#ok<PFBNS>
                    if scenario=="in_vivo_warfarin"
                        tY_red_inf=tY_red_inf*Y_refs_1(ind);
                    end
                    relerr_inf(ind) = errfun(t_ref,Y_refs(ind,:)',tY_red_inf); %#ok<PFBNS>
                catch
                    fprintf('\b\b (error inf), ');
                    relerr_inf(ind) = Inf;
                end
%                 if sum(relerr_inf>relerrTOL)>nind*(1-share)
%                     break
%                 end
            end
            %fprintf('\b\b\b\b\b\b\b\b\b\b\b')
            if sum(isfinite(relerr_inf))>=share && sum(relerr_inf < relerrTOL)/length(relerr_inf)>=share
                redmodel = tmodel_inf;
                fprintf('\b\b (inf), ');
            end
        end
        %%% determine reduced model dependend of error criterion
    end
    warning on;
    
    % remove identified env and neg states from sequence of states
    rI = redmodel.I;
    seqofpars = setdiff(seqofpars,[rI.pneg,rI.pinf],'stable'); %,rI.psma
    classificationstatus(redmodel.I);
end % of performing steps 1-2 and interim step multiple times


%%% =======================================================================
%%% report about final reduced model
%%% =======================================================================

rI = redmodel.I;

fprintf('\n\n Reduced model consists of');
fprintf('\n %d remaining parameter(s): ',length(rI.p));
fprintf('%s, ',I.nmpar{rI.p})
fprintf('\n %d neglected parameter(s): ',length(rI.pneg));
fprintf('%s, ',I.nmpar{rI.pneg})
fprintf('\n %d inf-neglected parameter(s): ',length(rI.pinf));
fprintf('%s, ',I.nmpar{rI.pinf})
%fprintf('\n %d fixed-to-0.00001 parameter(s): ',length(rI.psma));
%fprintf('%s, ',I.nmpar{rI.psma})
fprintf('\n')

elapsedtime = toc; fprintf(' [parameter reduction = %.1f sec]\n\n',elapsedtime);

%%% simulate final reduced model and assign output & approximation error
tic;
fprintf(1,'(%4d/%4d)',0,nind);
relerr = zeros(nind,1);
for ind=1:nind
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b%4d/%4d)',ind,nind);
    [~,X_red,te] = model.simulateODE(t_ref,X0_refs(ind,:),pars(ind,:), redmodel);
    Y_reds(ind,:) = rI.h_red(t_ref,X_red,te);
    X_reds(ind,1:size(X_red,1),:) = X_red;
    if redmodel.scenario=="in_vivo_warfarin"
        Y_reds(ind,:)=Y_reds(ind,:)*Y_refs(ind,1);    
    end
    relerr(ind) = errfun(t_ref,Y_refs(ind,:),Y_reds(ind,:));
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
