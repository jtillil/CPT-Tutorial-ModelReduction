%% model reduction warfarin
% original results are saved in the folder "original_results/" and can be
% copied from there to the folder "results/"

% choose population size (original n=1000 might take very long)
sample_size = 10;


%% ir-indices and parameter sensitivities are recalculated, only if recalculate = true
% choose if you want to re-calculate ir-indices and parameter sensitivities
recalculate = false;

% calculate parameter sensitivities for in vitro scenario
if recalculate
    addpath('blood_coagulation_modelspecification_parameter_sensitivity/')
    calculate_jacobian
    model_par_sens = simulate_and_reduce_parameter_sensitivity('in_vitro_warfarin',...
        'reduce',true, 'ir',false, 'output_time',(0:30)/3600);
    obs=model_par_sens.obs;
    save('results/parameters_sensitivity_in_vitro_warfarin.mat','obs')

%%% parameter sensitivities for in vivo scenario
% calculated for constant infusion to save runtime
    model_par_sens = simulate_and_reduce_parameter_sensitivity('in_vivo_warfarin',...
        'reduce',true, 'ir',false, 'output_time',0:24:720);
    obs=model_par_sens.obs;
    save('results/parameters_sensitivity_in_vivo_warfarin.mat','obs')
    rmpath('blood_coagulation_modelspecification_parameter_sensitivity/')

% calculate ir-indices for in vitro scenario
    model = simulate_and_reduce('in_vitro_warfarin','reduce',true, 'ir',false,'output_time',0:0.1/3600:15/3600);
    ir=model.ir;
    save('results/ir_in_vitro_warfarin.mat','ir')

%%% calculate ir-indices for in vivo scenario
% calculated for constant infusion to save runtime
    model = simulate_and_reduce('in_vivo_warfarin','reduce',true,'ir',false,...
        'output_time',0:2:720,'output',@(~,x,varargin)...
        (x(:,I.II).*x(:,I.VII).*x(:,I.X)/(x(1,I.II).*x(1,I.VII).*x(1,I.X))).^-0.1944);
    ir=model.ir;
    save('results/ir_in_vivo_warfarin.mat','ir')
end



%% in vitro scenario

%%% simulate warfarin in vivo to generate ensemble of blood samples (steps 1+2)
model = simulate_and_reduce('in_vitro_warfarin');
I=model.I;
model_in_vivo = simulate_and_reduce('in_vivo_warfarin','output_time',...
    [0,24,96,264,720],'multiple',true,'input_events',0:24:720,'sample_size',...
    sample_size,'sigma',0.16,'output', @(~,x,varargin) ...
    (x(:,I.II).*x(:,I.VII).*x(:,I.X)/(x(1,I.II).*x(1,I.VII).*x(1,I.X))).^-0.1975);
pars=repmat(model_in_vivo.pars,4,1);
inits=reshape(model_in_vivo.X_refs(:,2:end,:),[],63);

%%% model order reduction in vitro model ensemble of blood samples (step 3)
% reduce model
red_model=simulate_and_reduce('in_vitro_warfarin','X0s', inits/3,'pars',pars,...
    'reduce',true,'output_time',0:1/3600:90/3600,'share',0.95);

load('results/parameters_sensitivity_in_vitro_warfarin.mat','obs')
[~,ind] = sort( obs(1,64:end) );
ipar = ind;

%%% parameter reduction
red_model.I.p=1:174;
red_model.I.pneg=[];
red_model.I.pinf=[];
red_model.I.sim=[];
parred_model = model_order_reduction_parameters_robust(red_model,ipar,0.1,'share',0.95);

save('results/in_vitro_reduction_parred.mat','red_model','parred_model')

%%% reaction simplification
simred_model=model_simplification_robust(parred_model,parred_model.I.dyn,0.1,'share',0.95);

%%% analytic solution (only for reference, because TF/c30=v36 is relevant)
I=red_model.I;

syms t
X_sym=cell2sym(red_model.I.nmstate(:));
N=63;
for i=1:N
TT = symfun(str2sym([char(parred_model.I.nmstate(i)),'(t)']), t); %declare each element as symbolic handle
X_sym(i)=TT;
end
X_sym(I.neg)=0;
% TF must be known, otherwise division by zero, II, VII and X should remain
% as variables, as they can change dependent on warfarin
X_sym(setdiff(I.env,[I.II,I.VII,I.X,I.Fg]))=red_model.X0_ref(setdiff(I.env,[I.II,I.VII,I.X,I.Fg]));
X_sym(t)=X_sym;

par=cell2sym(red_model.I.nmpar(:));
par(setdiff(1:174,[I.v13,I.v14,I.k13,I.k14,I.v34,I.k34]))=...
    simred_model.par(setdiff(1:174,[I.v13,I.v14,I.k13,I.k14,I.v34,I.k34]));

%par=simred_model.par;
par(simred_model.I.pneg)=0;
%par(simred_model.I.pinf)=Inf;

red_odes=simred_model.odefun(t,X_sym(t),par,simred_model);

X_sym=cell2sym(red_model.I.nmstate(:));
N=63;
for i=1:N
TT = symfun(str2sym([char(parred_model.I.nmstate(i)),'(t)']), t); %declare each element as symbolic handle
X_sym(i)=TT;
end
X_sym(t)=X_sym;

X0=cell2sym(red_model.I.nmstate(:));
for i=1:N
TT = symfun(str2sym([char(simred_model.I.nmstate(i)),'_0']), t); %declare each element as symbolic handle
X0(i)=TT;
end
%X0(red_model.X0_ref==0)=0;
X0(setdiff(1:63,[I.II,I.VII,I.X,I.Fg]))=red_model.X0_ref(setdiff(1:63,[I.II,I.VII,I.X,I.Fg]));
X0(simred_model.I.neg)=0;
cond = X_sym(0) == X0;
solution=dsolve(diff(X_sym,t)==red_odes,cond);
simple_solution=simplify(solution.F,'Steps',50);

save('results/in_vitro_reduction.mat','red_model','parred_model','simred_model','solution','simple_solution')

%% in vivo scenario
% after reduction of in vitro part, we have INR equation as output for in
% vivo part

%%% model order reduction (robust) in vivo with INR equation as output
red_model_in_vivo = simulate_and_reduce('in_vivo_warfarin','reduce',true,...
    'multiple',true,'input_events',0:24:720,'output_time',0:24:720,'share',0.95,...
    'sample_size',sample_size, 'sigma',0.16,'output', @(~,x,varargin) simulate_and_reduce(...
    'in_vitro_warfarin', 'initial_values', x(1,:),'parameters',varargin{2}.par).Y_ref(1)/...
    3.068837*1000*(x(:,I.II).*x(:,I.VII).*x(:,I.X)/(x(1,I.II).*x(1,I.VII).*x(1,I.X))).^-0.1975);

load('results/parameters_sensitivity_in_vivo_warfarin.mat','obs')
[~,ind] = sort( obs(1,64:end) );
ipar = ind;

%%% parameter reduction
red_model_in_vivo.I.p=1:174;
red_model_in_vivo.I.pneg=[];
red_model_in_vivo.I.pinf=[];
red_model_in_vivo.I.sim=[];
parred_model_in_vivo = model_order_reduction_parameters_robust(red_model_in_vivo,ipar,0.1,'share',0.95);

save('results/in_vivo_reduction_parred.mat','red_model_in_vivo','parred_model_in_vivo')

%%% reaction simplification
simred_model_in_vivo=...
    model_simplification_robust(parred_model_in_vivo,parred_model_in_vivo.I.dyn,0.1,'share',0.95);

save('results/in_vivo_reduction.mat','red_model_in_vivo','parred_model_in_vivo','simred_model_in_vivo')

%% simulation results for figures
% Figure 4 loglog-plot
load('results/reference_daily_1000.mat')
for i=1:1000
ref_model.INR_eq(i,1:31)=ref_model.I.h_red(ref_model.t_ref,reshape(ref_model.X_refs(i,:,:),[],63));
end
ref_model.INR_eq=ref_model.INR_eq.*ref_model.Y_refs(:,1);
save('results/loglogplot.mat','ref_model')

% Figure 5a approximation in population

save('results/population_approximation.mat','ref_model','red_model')

% Figure 5b different CYP2C9 genotypes
covariates.cyp=[2,0,0];
model11=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates);
model11_red=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates,'reduced_only',true,'ref_dyn',[I.Cwarf,I.Awarf,I.VKH2,I.II,I.VII,I.X],...
    'ref_qss',[], 'ref_env', [I.VK]);
model11.Y_red=model11_red.Y_red;
covariates.cyp=[1,1,0];
model12=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates);
model12_red=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates,'reduced_only',true,'ref_dyn',[I.Cwarf,I.Awarf,I.VKH2,I.II,I.VII,I.X],...
    'ref_qss',[], 'ref_env', [I.VK]);
model12.Y_red=model12_red.Y_red;
covariates.cyp=[0,2,0];
model22=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates);
model22_red=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates,'reduced_only',true,'ref_dyn',[I.Cwarf,I.Awarf,I.VKH2,I.II,I.VII,I.X],...
    'ref_qss',[], 'ref_env', [I.VK]);
model22.Y_red=model22_red.Y_red;
covariates.cyp=[1,0,1];
model13=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates);
model13_red=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates,'reduced_only',true,'ref_dyn',[I.Cwarf,I.Awarf,I.VKH2,I.II,I.VII,I.X],...
    'ref_qss',[], 'ref_env', [I.VK]);
model13.Y_red=model13_red.Y_red;
covariates.cyp=[0,1,1];
model23=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates);
model23_red=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates,'reduced_only',true,'ref_dyn',[I.Cwarf,I.Awarf,I.VKH2,I.II,I.VII,I.X],...
    'ref_qss',[], 'ref_env', [I.VK]);
model23.Y_red=model23_red.Y_red;
covariates.cyp=[0,0,2];
model33=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates);
model33_red=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates,'reduced_only',true,'ref_dyn',[I.Cwarf,I.Awarf,I.VKH2,I.II,I.VII,I.X],...
    'ref_qss',[], 'ref_env', [I.VK]);
model33.Y_red=model33_red.Y_red;
model33_doseopt=simulate_and_reduce('in_vivo_warfarin','multiple',true,...
    'output_time',0:2:720,'covariates',covariates,'input_doses',ones(1,31));
model33_doseopt_red=simulate_and_reduce('in_vivo_warfarin','multiple',true,...
    'output_time',0:2:720,'covariates',covariates,'input_doses',ones(1,31),...
    'reduced_only',true,'ref_dyn',[I.Cwarf,I.Awarf,I.VKH2,I.II,I.VII,I.X] ,...
    'ref_qss',[], 'ref_env', [I.VK]);
model33_doseopt.Y_red=model33_doseopt_red.Y_red;
save('results/cyp_simulations.mat','model11','model12','model22','model13',...
    'model23','model33','model33_doseopt')

% Figure 5c different VKORC1 genotypes
covariates.cyp=[2,0,0];
covariates.vko=[2,0];
modelv11=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates);
modelv11_red=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates,'reduced_only',true,'ref_dyn',[I.Cwarf,I.Awarf,I.VKH2,I.II,I.VII,I.X],...
    'ref_qss',[], 'ref_env', [I.VK]);
modelv11.Y_red=modelv11_red.Y_red;
covariates.vko=[1,1];
modelv12=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates);
modelv12_red=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates,'reduced_only',true,'ref_dyn',[I.Cwarf,I.Awarf,I.VKH2,I.II,I.VII,I.X],...
    'ref_qss',[], 'ref_env', [I.VK]);
modelv12.Y_red=modelv12_red.Y_red;
covariates.vko=[0,2];
modelv22=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates);
modelv22_red=simulate_and_reduce('in_vivo_warfarin','multiple',true,'output_time',0:2:720,...
    'covariates',covariates,'reduced_only',true,'ref_dyn',[I.Cwarf,I.Awarf,I.VKH2,I.II,I.VII,I.X],...
    'ref_qss',[], 'ref_env', [I.VK]);
modelv22.Y_red=modelv22_red.Y_red;
save('results/vko_simulations.mat','modelv11','modelv12','modelv22')

% Figure Supplement alternative VKORC1 genotype parametrisation
model22=simulate_and_reduce('in_vivo_warfarin','multiple', true, 'output_time',0:2:720,...
    'par_changed',"degVK2",'new_par',0.96/2.05,'init_changed',...
    ["VKH2","PC","PS","II","VII","IX","X"],'new_init',ones(7,1)*0.96/2.05);
model12=simulate_and_reduce('in_vivo_warfarin','multiple', true,'output_time',0:2:720,...
    'par_changed',"degVK2",'new_par',(2.05+0.96)/4.1,'init_changed',...
    ["VKH2","PC","PS","II","VII","IX","X"],'new_init',ones(7,1)*(2.05+0.96)/4.1);
model11=simulate_and_reduce('in_vivo_warfarin','multiple', true,'output_time',0:2:720);
save('results/vko_simulations_alternative.mat','model11','model12','model22')

