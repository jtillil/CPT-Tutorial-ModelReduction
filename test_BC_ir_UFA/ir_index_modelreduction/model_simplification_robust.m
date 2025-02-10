%%% Version: September 20th, 2021
%%%
%%% call by: redmodel  =  model_simplification(model,seqofstates,relerrTOL)
%%%
%%% This function returns which states can have simplified reactions 
%%%
%%% Input:  model               structure specifying the model
%%%         seqofstates         sequence, in which state variables are
%%%                             tested for model order reduction
%%%         relerrTOL           user defined relative error threshold
%%%
%%% Output: redmodel            structure specifying the reduced order model
%%%
%%% 

function redmodel = model_simplification_robust(model,seqofstates,relerrTOL,varargin)
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
pars    = model.pars;
Y_refs  = model.Y_refs;
Y_refs_1= Y_refs(:,1);
h_red   = I.h_red;

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
    
    fprintf('\n 1st step - simplify ODEs of state variables');
    fprintf('\n          > tested state variables: ');
    
    par_sym=cell2sym(model.I.nmpar(:));
    par_sym_full=par_sym;
    par_sym(I.pneg)=0;
    par_sym(I.pinf)=Inf;
    X_sym=cell2sym(model.I.nmstate(:));
    X_sym_full=X_sym;
    X_sym(I.neg)=0;
    
    warning off;
    for k = seqofstates
        
        fprintf('%s, ',I.nmstate{k});
        rI = redmodel.I;
        
        %%% first tentatively consider the kth state to be simplifiable
        %%%
        tmodel_sim = redmodel;
        tmodel_sim.I.dyn = setdiff(rI.dyn,k,'stable');
        tmodel_sim.I.sim = [rI.sim,k];   % add k to neg states
        
        % compute solution of tentatively reduced model and corresponding output
        ode_fun=redmodel.odefun(redmodel.t_ref,X_sym,par_sym,redmodel);
        red_fun=taylor(ode_fun(k),X_sym_full,'Order', 3);
        ode_fun(k)=red_fun;
        ode_fun=matlabFunction(ode_fun,'Vars',{X_sym_full,par_sym_full});
        tmodel_sim.odefun=@(t,X,par,model) ode_fun(X,par);
        
        tmodel_sim_reduced_dose = tmodel_sim;
        tmodel_sim_reduced_dose.u_ref = tmodel_sim_reduced_dose.u_ref/4;
        tmodel_sim_reduced_dose.input_doses = tmodel_sim_reduced_dose.input_doses/4;
        tmodel_sim_reduced_dose.X0_refs = tmodel_sim_reduced_dose.X0s + tmodel_sim_reduced_dose.u_ref'; 
        
        nind = size(X0_refs,1);
        %fprintf(1,'(%4d/%4d)',0,nind);
        relerr_sim = zeros(nind,1);
        parfor ind=1:nind
            if reduced_dose(ind)
                tmodel_sim_i=tmodel_sim_reduced_dose;
            else
                tmodel_sim_i=tmodel_sim;
            end
            X0_refs(ind,:)=tmodel_sim_i.X0_refs(ind,:);
            %fprintf(1,'\b\b\b\b\b\b\b\b\b\b%4d/%4d)',ind,nind);
            try
                [~,tX_red_sim,te_red_sim] = tmodel_sim_i.simulateODE(t_ref,X0_refs(ind,:)',pars(ind,:)', tmodel_sim_i);
                tY_red_sim = h_red(t_ref,tX_red_sim,te_red_sim); %#ok<PFBNS>
                if scenario=="in_vivo_warfarin"
                    tY_red_sim=tY_red_sim*Y_refs_1(ind);
                end
                relerr_sim(ind) = errfun(t_ref,Y_refs(ind,:)',tY_red_sim); %#ok<PFBNS>
            catch
                % numerics had problems -> ignore reduced model and continue
                fprintf('\b\b (error lin), ');
                relerr_sim(ind) = Inf;
            end
        end
        fprintf('\b\b\b\b\b\b\b\b\b\b\b')
        
        %%% determine reduced model dependend of error criterion
        %%%
        if sum(isfinite(relerr_sim))>=share && sum(relerr_sim < relerrTOL)/length(relerr_sim)>=share
            redmodel = tmodel_sim;
            fprintf('\b\b (lin), ');
        end
    end
    warning on;
    
    % remove identified env and neg states from sequence of states
    rI = redmodel.I;
    seqofstates = setdiff(seqofstates,[rI.sim],'stable');
    classificationstatus(redmodel.I);
end % of performing steps 1-2 and interim step multiple times


%%% =======================================================================
%%% report about final reduced model
%%% =======================================================================

rI = redmodel.I;

fprintf('\n\n Reduced model consists of');
fprintf('\n %d dynamical nonlinear state variable(s): ',length(rI.dyn));
fprintf('%s, ',I.nmstate{rI.dyn})
fprintf('\n %d environmental state variable(s): ',length(rI.env));
fprintf('%s, ',I.nmstate{rI.env})
fprintf('\n %d neglected state variable(s): ',length(rI.neg));
fprintf('%s, ',I.nmstate{rI.neg})
fprintf('\n %d state variable(s) approximated by their quasi-steady state: ',length(rI.qss));
fprintf('%s, ',I.nmstate{rI.qss})
fprintf('\n %d state variable(s) eliminated by conservation law(s): ',length(rI.con));
fprintf('%s, ',I.nmstate{rI.con})
fprintf('\n %d state variable(s) with linear reaction(s): ',length(rI.sim));
fprintf('%s, ',I.nmstate{rI.sim})
fprintf('\n')


elapsedtime = toc; fprintf(' [model simplification = %.1f sec]\n\n',elapsedtime);

%%% simulate final reduced model and assign output & approximation error
tic;
for ind=1:size(X0_refs,1)
    [~,X_red,te] = model.simulateODE(t_ref,X0_refs(ind,:)',pars(ind,:)', redmodel);
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
