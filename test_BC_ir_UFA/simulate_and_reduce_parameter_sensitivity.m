% only needed for calculation of parameter sensitivities
function model = simulate_and_reduce_parameter_sensitivity(scenario,varargin)
% scenario: 'in_vivo_warfarin', 'in_vitro_warfarin',
%   'in_vivo_snakevenom_40h' or 'in_vitro_PTtest_highTF'
% event: either 0 or threshold for AUC (e.g. 1500)
% output_time: tspan for ode-solver
% seed: seed for random parameters, if not 0 (then deterministic)
% sample_size: number of random parameter sets
% reduce: false or true, if reduction should be computed
% nrepepeats: how often the reduction shall be repeated in
%       model_order_reduction_multiple.m
% relerrTOL: relative error tolerance for the reduction
% output: h function to give output
% ir: true or false, if ir indices should be loaded from file
% reduced_only: false or true, if only reduced model shall be simulated,
%       reduced with ref_dyn, ref_qss and ref_env, for output_time
% share: share of population that needs to attain the error for a state to
%       be reduced
% multiple: false or true, if multiple dosing (instead of constant infusion)
% input_events: times of warfarin doses
% input_doses: warfarin doses for the input_events
    p = inputParser();

    addParameter(p, 'save_file', '', @ischar);
    addParameter(p, 'load_file', '', @ischar);
    addParameter(p, 'event', 0, @isnumeric);
    addParameter(p, 'par_changed', '', @isstring);
    addParameter(p, 'new_par', 0, @isnumeric);
    addParameter(p, 'init_changed', '', @isstring);
    addParameter(p, 'new_init', 0, @isnumeric);
    addParameter(p, 'parameters', 0, @isnumeric);
    addParameter(p, 'initial_values', 0, @isnumeric);
    addParameter(p, 'location', 0, @isnumeric);
    addParameter(p, 'output_time', 0, @isnumeric);
    addParameter(p, 'seed', 0, @isnumeric);
    addParameter(p, 'sample_size', 0, @isnumeric);
    addParameter(p, 'factors_par', 0, @isnumeric);
    addParameter(p, 'factors_init', 0, @isnumeric);
    addParameter(p, 'reduce', false, @islogical);
    addParameter(p, 'ref_dyn', [], @isnumeric);
    addParameter(p, 'ref_qss', [], @isnumeric);
    addParameter(p, 'ref_env', [], @isnumeric);
    addParameter(p, 'nrepeats', 1, @isnumeric);
    addParameter(p, 'relerrTOL', 0.1, @isnumeric);
    addParameter(p, 'output', @(t,x,varargin) 0);
    addParameter(p, 'ir', true, @islogical);
    addParameter(p, 'ref_sol', {}, @isstruct);
    addParameter(p, 'reduced_only', false, @islogical);
    addParameter(p, 'share', 1, @isnumeric);
    addParameter(p, 'multiple', false, @islogical);
    addParameter(p, 'input_events', 0:24:720, @isnumeric);
    addParameter(p, 'input_doses', 0, @isnumeric);
    addParameter(p, 'covariates', 0, @isstruct);
    
    parse(p,varargin{:})
 
    save_file = p.Results.save_file;
    load_file = p.Results.load_file;
    event = p.Results.event;
    par_changed = p.Results.par_changed;
    new_par = p.Results.new_par;
    init_changed = p.Results.init_changed;
    new_init = p.Results.new_init;
    parameters = p.Results.parameters;
    initial_values = p.Results.initial_values;
    location = p.Results.location;
    output_time = p.Results.output_time;
    seed = p.Results.seed;
    sample_size = p.Results.sample_size;
    factors_par = p.Results.factors_par;
    factors_init = p.Results.factors_init;
    reduce = p.Results.reduce;
    ref_dyn = p.Results.ref_dyn;
    ref_qss = p.Results.ref_qss;
    ref_env = p.Results.ref_env;
    nrepeats = p.Results.nrepeats;
    relerrTOL = p.Results.relerrTOL;
    output = p.Results.output;
    ir = p.Results.ir;
    ref_sol = p.Results.ref_sol;
    reduced_only = p.Results.reduced_only;
    share = p.Results.share;
    multiple = p.Results.multiple;
    input_events = p.Results.input_events;
    input_doses = p.Results.input_doses;
    covariates = p.Results.covariates;
    
addpath('blood_coagulation_modelspecification_parameter_sensitivity/','ir_index_modelreduction/','blood_coagulation_modelspecification_parameter_sensitivity/optional')
model.scenario = scenario; 
model.multiple=multiple;
model.input_events=input_events;
model.input_doses=input_doses;

%% initialise model
% indexing of states and parameters
I = Wajima2009BloodCoagulation_indexing;
model.I = I;
% initial values and modifications according to the scenarios
X0 = Wajima2009BloodCoagulation_initialvalues(model);
% initial value variabilities
for i=1:length(init_changed)
    X0(I.(init_changed(i)))=new_init(i)*X0(I.(init_changed(i)));
end
% if initial values are provided, use them 
if any(initial_values~=0)
    X0=initial_values;
    if ~any(parameters)
        "Parameters will be adapted to be in steady state with given initial values, if not provided via 'parameters'."
    end
end
model.X0 = X0;

% if covariates defined, use them
if isstruct(covariates)
    model.covariates=covariates;
end

% parameter values
par = Wajima2009BloodCoagulation_parameters(model);
par(I.v43) = 0; 
model.par = par;
% right hand side of ODEs and extended ODEs
model.odefun = @(t,X,par,model) Wajima2009BloodCoagulation_ode(t,X,par,model);
if strcmp(scenario,"in_vivo_warfarin")
    if ismember("jacfun_parameter_sensitivity_in_vivo.mat",arrayfun(@(x) x.name,dir('results'),'uni',false))
        load('results/jacfun_parameter_sensitivity_in_vivo.mat','jac_fun');
        model.jacfun = @(t,X,par,model) jac_fun(X,par); %@(t,X,par,model) Wajima2009BloodCoagulation_odejac(t,X,par,model);
    else
        fprintf("\nJacobian unknown, please run 'calculate_jacobian' first for calculating parameter sensitivity\n\n");
        model.jacfun = [];
    end
elseif strcmp(scenario,"in_vitro_warfarin")
    if ismember("jacfun_parameter_sensitivity_in_vitro.mat",arrayfun(@(x) x.name,dir('results'),'uni',false))
        load('results/jacfun_parameter_sensitivity_in_vitro.mat','jac_fun');
        model.jacfun = @(t,X,par,model) jac_fun(X,par); %@(t,X,par,model) Wajima2009BloodCoagulation_odejac(t,X,par,model);
    else
        fprintf("\nJacobian unknown, please run 'calculate_jacobian' first for calculating parameter sensitivity\n\n");
        model.jacfun = [];
    end
else
    fprintf('\n--> unknown model.scenario :-( Please fix!\n\n');
    error('')
end
% function handle to simulate the ODE, or in case of qss states, the DAE  
model.simulateODE = @(t,X,par,model) simmodelODE(t,X,par,model);
% define colors for all state variables
%model.state2col = Wajima2009BloodCoagulation_state2color(I);
% possible event function
if(event)
    model.event=@(t,concentrations) event_fct(t,concentrations,event/3600,model.I.AUC,1);
else
    model.event=[];
end

%% specify case
u_ref = zeros(size(X0)); 
fprintf('\n Scenario: %s\n',model.scenario);
switch model.scenario
    case 'in_vivo_warfarin'
        if sum(output_time)==0
            tspan = [0 24*30]; % simulation time span (30 days)
        else
            tspan = output_time;
        end
        model.odefun = @(t,X,par,model) Wajima2009BloodCoagulation_ode_const_infusion(t,X,par,model);
        I.input = I.Awarf;
        I.output = I.VII;
        if multiple
            if any(input_doses)
                u_ref(I.input)=input_doses(1);
            else
                u_ref(I.input) = par(I.warf_dose);
            end
            model.odefun=@(t,X,par,model) Wajima2009BloodCoagulation_ode(t,X,par,model);
        end
        %I.h_red = @(t,x) x(:,I.lump).^-
        I.h = @(~,x,varargin) (x(:,I.II).*x(:,I.VII).*x(:,I.X)/(x(1,I.II).*x(1,I.VII).*x(1,I.X))).^-0.1975;
        %I.h = @(~,x,varargin) arrayfun(@(n) (simulate_and_reduce_parameter_sensitivity('in_vitro_warfarin', 'initial_values', x(n,:)/3,'parameters',varargin{2}.par).Y_ref(1)), 1:size(x,1))'/3.068837*1000;%simredmodelODE(x);%simredmodelODE(x);
    case 'in_vivo_ss'
        if sum(output_time)==0
            tspan = [0 24*30]; % simulation time span (30 days)
        else
            tspan = output_time;
        end
        model.odefun = @(t,X,par,model) Wajima2009BloodCoagulation_ode(t,X,par,model);

        I.h_red = @(~,x,varargin) (x(:,I.II).*x(:,I.VII).*x(:,I.X)/(x(1,I.II).*x(1,I.VII).*x(1,I.X))).^-0.1944;
        I.h = @(~,x,varargin) arrayfun(@(n) (simulate_and_reduce('in_vitro_warfarin', 'initial_values', x(n,:)/3,'parameters',varargin{2}.par).Y_ref(1)), 1:size(x,1))'/3.068837*1000;%simredmodelODE(x);%simredmodelODE(x);
    case 'in_vitro_warfarin'
        I.input = I.TF;
        I.output = I.AUC;
        %I.h = @(t,x) repmat([t(find(x(:,I.AUC)>1500/3600,1)) 1],length(t),1);
        I.h = @(t,x,varargin) repmat(max([0; varargin{1}]),length(1),1);
        I.h_red = @(t,x,varargin) repmat(max([0; varargin{1}]),length(1),1);
        model.event=@(t,concentrations) event_fct(t,concentrations,1500/3600,model.I.AUC,1);
        if ~isempty(load_file)
            load(load_file, 'X_invitro_init', 't_invitro_init');
            model.X0=X_invitro_init(location,:)'/3;
        else
            model.X0=model.X0/3; % dilution
        end
        model.X0(I.Tmod)=0;
        u_ref(I.input) = 100; % in [nM]
        if sum(output_time)==0
            tspan = [0 120/3600]; % simulation time span (60sec)
        else
            tspan = output_time;
        end
    case {'in_vivo_snakevenom_1h','in_vivo_snakevenom_40h'}
        % define input/output state variable
        I.input  = I.AVenom;
        I.output = I.Fg;
        I.h = @(t,x,varargin) x(:,I.output);
        I.h_red=I.h;
        dose_snake_venom = 0.0015; % amount in mg
        SF_mg_to_nmol = 1e-3/2e5*1e9; 
        u_ref(I.input) = SF_mg_to_nmol*dose_snake_venom;
        if sum(output_time)==0
            if strcmp(model.scenario, 'in_vivo_snakevenom_1h')
                tspan = [0 1];
            else
                strcmp(model.scenario, 'in_vivo_snakevenom_40h')
                tspan = [0 40];
            end
        else
            tspan = output_time;
        end
    case 'in_vitro_PTtest_highTF'
        % define input/output state variable
        I.input  = I.TF;
        I.output = I.F;
        I.h = @(t,x,varargin) repmat(max([0; varargin{1}]),length(t),1);
        
        model.X0=model.X0/3; % dilution
        model.event=@(t,concentrations) event_fct(t,concentrations,1500/3600,model.I.AUC,1);
        % specify input u_ref for reference solution
        u_ref(I.input) = 100; % in [nM]
        % simulation time span (30sec)
        tspan = [0 30/3600];
    case 'in_vitro_PTtest_highTF_PT'
        % define input/output state variable
        I.input  = I.TF;
        I.output = I.F;
        I.h = @(t,x,varargin) repmat(t(find(x(:,I.AUC)>1500/3600,1)),length(t),1);
        % specify input u_ref for reference solution
        u_ref(I.input) = 100; % in [nM]
        % simulation time span (30sec)
        tspan = [0 30/3600];
    case 'in_vitro_PTtest_lowTF'
        % define input/output state variable
        I.input  = I.TF;
        I.output = I.F;
        I.h = @(t,x,varargin) repmat(max([0; varargin{1}]),length(t),1);
        model.event=@(t,concentrations) event_fct(t,concentrations,1500/3600,model.I.AUC,1);
        % specify input u_ref for reference solution
        u_ref(I.input) = 5e-3; % in [nM]
        % simulation time span (30sec)
        tspan = [0 240/3600];
    otherwise
        fprintf('\n--> unknown model.scenario :-( Please fix!\n\n');
        error('')
end
if sum(output(tspan,X0(:)',10,model)) ~= 0
    I.h=output;
    I.h_red=I.h;
end
% parameter value variabilities
for i=1:length(par_changed)
    par(I.(par_changed(i)))=new_par(i)*par(I.(par_changed(i)));
end
% calculate parameters dependent on other parameters to maintain steady
% state
par(I.degVKH2)  = par(I.degVK2) * (X0(I.VK)./ (X0(I.VKH2)));
par(I.degVKO)   = par(I.degVK2) * (X0(I.VK)./ (X0(I.VKO)));
par(I.aII)      = par(I.degII)  * X0(I.II)  / (X0(I.VKH2));
par(I.aVII)     = par(I.degVII) * X0(I.VII) / (X0(I.VKH2)); 
par(I.aIX)      = par(I.degIX)  * X0(I.IX)  / (X0(I.VKH2));
par(I.aX)       = par(I.degX)   * X0(I.X)   / (X0(I.VKH2));
par(I.aPC)      = par(I.degPC)  * X0(I.PC)  / (X0(I.VKH2)); 
par(I.aPS)      = par(I.degPS)  * X0(I.PS)  / (X0(I.VKH2)); 
par(I.pV)       = par(I.degV)    * X0(I.V); 
par(I.pVIII)    = par(I.degVIII) * X0(I.VIII);
par(I.pFg)      = par(I.degFg)   * X0(I.Fg);  
par(I.pXIII)    = par(I.degXIII) * X0(I.XIII);
par(I.pPg)      = par(I.degPg)   * X0(I.Pg);   
par(I.pTmod)    = par(I.degTmod) * X0(I.Tmod);
par(I.pXI)      = par(I.degXI)   * X0(I.XI);  
par(I.pTFPI)    = par(I.degTFPI) * X0(I.TFPI);
par(I.pVK)      = par(I.degVK)   * X0(I.VK);  
par(I.pXII)     = par(I.degXII)  * X0(I.XII);
par(I.pPk)      = par(I.degPk)   * X0(I.Pk);

% if parameters are provided, use them 
if any(parameters~=0)
    par=parameters;
end

% set synthesis parameters zero in "in vitro" cases
if contains(model.scenario,'in_vitro')
        par(I.aII) = 0;    par(I.aVII)  = 0;
        par(I.aIX) = 0;    par(I.aX)    = 0;
        par(I.aPC) = 0;    par(I.aPS)   = 0;
        par(I.pV)  = 0;    par(I.pVIII) = 0;
        par(I.pFg) = 0;    par(I.pXIII) = 0;
        par(I.pPg) = 0;    par(I.pTmod) = 0;
        par(I.pXI) = 0;    par(I.pTFPI) = 0;
        par(I.pVK) = 0;    par(I.pXII)  = 0;
        par(I.pPk) = 0;
end

% adjust VK_p to ensure steady state
model.X0(I.VK_p)=X0(I.VK)*par(I.VK_V)*par(I.VK_k12)/par(I.VK_k21);

% initial value of reference model
model.X0_ref = model.X0 + u_ref;

model.u_ref=u_ref;
model.par=par;
model.L.nconlaws = 0;
I.dyn = 1:I.nstates; I.env = []; I.neg = []; I.qss = []; I.con = [];
model.I = I;


%% generate population
if sample_size||(any(any(factors_par)) && any(any(factors_init))) 
    if sample_size
        rng(seed)
        factors_par=exp(lhsnorm(zeros(I.npar,1),0.2*eye(I.npar),1000));
        factors_init=exp(lhsnorm(zeros(I.nstates,1),0.2*eye(I.nstates),1000));
        factors_par=factors_par(1:sample_size,:);
        factors_init=factors_init(1:sample_size,:);
    end
    factors_par(:,I.lmax)=1;
    factors_par(:,I.warf_dose)=1;
    model.pars=model.par'.*factors_par;
    model.X0s=model.X0.*factors_init;
    model.pars(:,I.degVKH2)  = model.pars(:,I.degVK2) .* (model.X0s(:,I.VK)./(model.X0s(:,I.VKH2)));
    model.pars(:,I.degVKO)   = model.pars(:,I.degVK2) .* (model.X0s(:,I.VK)./(model.X0s(:,I.VKO)));
    model.pars(:,I.aII)      = model.pars(:,I.degII)  .* model.X0s(:,I.II) ./ (model.X0s(:,I.VKH2));
    model.pars(:,I.aVII)     = model.pars(:,I.degVII) .* model.X0s(:,I.VII)./ (model.X0s(:,I.VKH2)); 
    model.pars(:,I.aIX)      = model.pars(:,I.degIX)  .* model.X0s(:,I.IX) ./ (model.X0s(:,I.VKH2));
    model.pars(:,I.aX)       = model.pars(:,I.degX)   .* model.X0s(:,I.X)  ./ (model.X0s(:,I.VKH2));
    model.pars(:,I.aPC)      = model.pars(:,I.degPC)  .* model.X0s(:,I.PC) ./ (model.X0s(:,I.VKH2)); 
    model.pars(:,I.aPS)      = model.pars(:,I.degPS)  .* model.X0s(:,I.PS) ./ (model.X0s(:,I.VKH2)); 
    model.pars(:,I.pV)       = model.pars(:,I.degV)   .* model.X0s(:,I.V); 
    model.pars(:,I.pVIII)    = model.pars(:,I.degVIII).* model.X0s(:,I.VIII);
    model.pars(:,I.pFg)      = model.pars(:,I.degFg)  .* model.X0s(:,I.Fg);  
    model.pars(:,I.pXIII)    = model.pars(:,I.degXIII).* model.X0s(:,I.XIII);
    model.pars(:,I.pPg)      = model.pars(:,I.degPg)  .* model.X0s(:,I.Pg);   
    model.pars(:,I.pTmod)    = model.pars(:,I.degTmod).* model.X0s(:,I.Tmod);
    model.pars(:,I.pXI)      = model.pars(:,I.degXI)  .* model.X0s(:,I.XI);  
    model.pars(:,I.pTFPI)    = model.pars(:,I.degTFPI).* model.X0s(:,I.TFPI);
    model.pars(:,I.pVK)      = model.pars(:,I.degVK)  .* model.X0s(:,I.VK);  
    model.pars(:,I.pXII)     = model.pars(:,I.degXII) .* model.X0s(:,I.XII);
    model.pars(:,I.pPk)      = model.pars(:,I.degPk)  .* model.X0s(:,I.Pk);
    model.X0_refs = model.X0s + u_ref;
    if isempty(ref_sol) && ~reduced_only
        for ind=1:size(factors_init,1) 
            if sum(output_time)==0
                error('Please define output_time, this is needed, s.t. output time is the same for all individuals');
            end
            [~,X_temp,te_temp] = model.simulateODE(tspan,model.X0_refs(ind,:),model.pars(ind,:),model);
            temp_model=model;
            temp_model.par=model.pars(ind,:);
            Y_refs(ind,:) = model.I.h(tspan,X_temp,te_temp,temp_model);
            X_refs(ind,1:size(X_temp,1),:)=X_temp;
        end
    elseif ~isempty(ref_sol)
        Y_refs=ref_sol.Y_refs;
        X_refs=ref_sol.X_refs;
    else
        model.t_ref = tspan;
        for ind=1:size(factors_init,1)
            temp_model=model;
            temp_model.X0_ref=model.X0_refs(ind,:);
            temp_model.par=model.pars(ind,:);
            temp_model=manualredmodeldesign(temp_model,ref_dyn,ref_qss,ref_env);
            model.Y_reds(ind,:)=temp_model.Y_red;
            model.X_reds(ind,:,:)=temp_model.X_red;
        end
        return
    end
    model.Y_refs= Y_refs;
    model.X_refs= X_refs;
end

%% compute reference solution
fprintf('\n Start model simulation');
% compute reference solution
%load('ref3.mat','X_ref','Y_ref','t_ref')
if isempty(ref_sol) && ~reduced_only
    tic
    [t_ref,X_ref,te] = model.simulateODE(tspan,model.X0_ref,par,model);
    toc
    Y_ref = model.I.h(t_ref,X_ref,te,model);
elseif ~isempty(ref_sol)
    t_ref=ref_sol.t_ref;
    X_ref=ref_sol.X_ref;
    Y_ref=ref_sol.Y_ref;
else % reduced_only
    model.t_ref = tspan;
    model=manualredmodeldesign(model,ref_dyn,ref_qss,ref_env);
    return
end
% assign values to model structure
model.t_ref = t_ref; model.X_ref = X_ref; model.Y_ref= Y_ref;

% further model definition
if strcmp(model.scenario,'in_vivo_warfarin')
    model.errfun = @(t_ref,Y_ref,Y_red) max(abs((Y_ref-Y_red)./Y_ref));
elseif strcmp(model.scenario,'in_vitro_warfarin')
    model.errfun = @(t_ref,Y_ref,Y_red) abs((Y_ref-Y_red)./Y_ref);
else
    model.errfun = @(t_ref,Y_ref,Y_red) sqrt( trapz(t_ref,(Y_ref-Y_red).^2) ) *(1./ sqrt( trapz(t_ref,Y_ref.^2) )');
end
model.relerrTOL = relerrTOL; 
model.nrepeats = nrepeats; 


%% compute model reduction
if reduce
    switch model.scenario
        case 'in_vivo_warfarin'
            ir_file="results/ir_warfarin_in_vivo.mat";
        case 'in_vitro_warfarin'
            ir_file="results/ir_INR_in_vitro.mat";
        otherwise
            fprintf('\n--> unknown model.scenario :-( Please fix!\n\n');
            error('')
    end
    load_ir_indices=ir;
    if ~load_ir_indices
        %load('results/extodejacpat.mat', 'extodejacpat')
        %model.options.JPattern = extodejacpat;
        if strcmp(model.scenario, 'in_vitro_warfarin') || strcmp(model.scenario, 'in_vitro_PTtest_highTF')
            obs  =  compute_parameter_sensitivity_PT(model);
        elseif strcmp(model.scenario, 'in_vivo_warfarin')
            obs  =  compute_parameter_sensitivity_warf(model);
        else
            obs  =  compute_parameter_sensitivity(model);
        end
    else
        load(ir_file,"ir")
    end
    model.obs = obs;
end
end

%% local subroutines


function redmodel = manualredmodeldesign(model,dyn,qss,env)
redmodel = model;
redmodel.I.dyn=dyn;
redmodel.I.qss=qss;
redmodel.I.env=env;
redmodel.I.neg = setdiff(1:model.I.nstates,[dyn,env,qss]);

%%% simulate reduced model
[t_red,X_red,te] = simmodelODE(model.t_ref, redmodel.X0_ref, redmodel.par, redmodel);

%%% assign output and approximation error of full and reduced model
if ~isfield(model.I,'h_red')
    model.I.h_red=model.I.h;
end
redmodel.t_red = t_red;
redmodel.X_red = X_red;
redmodel.Y_red = model.I.h_red(t_red,X_red,te);
if isfield(model,'Y_ref')
    if model.scenario=="in_vivo_warfarin"
        redmodel.Y_red=redmodel.Y_red*model.Y_ref(1);
    end
    %if size(redmodel.Y_red)==size(redmodel.Y_ref)
    redmodel.relerr = model.errfun(t_red,model.Y_ref,redmodel.Y_red);
end
end


function [t,X,te] = simmodelODE(tspan,X0,par,model)

% indexing
I = model.I; L = model.L;

% sharper tolerances (sometimes needed to resolve numerical problems
options.RelTol = 1e-5; options.AbsTol = 1e-7;

% default, unless changed below
options.NonNegative = 1:I.nstates;

% set negligible state variables to zero
if ~isempty(I.neg)
    X0(I.neg) = 0;
end
if isfield(I,'pneg')
    par(I.pneg)=0;
end
if isfield(I,'pinf')
    par(I.pinf)=Inf;
end
if isfield(I,'psma')
    par(I.psma)=0.00001;
end
% define mass matrix for qss states
if ~isempty(I.qss)
    % DAE solver not able to deal with NonNegative option
    options.NonNegative = []; 
    
    % check whether initial value of tentative con state is zero, set to small positive values to resolve problems
    qss0 = I.qss(X0(I.qss)==0);
    X0(qss0) = 1e-3;
    
    % construct mass matrix for DAE solver
    M = eye(I.nstates); 
    M(I.qss,I.qss) = 0;
    options.Mass = M;
    options.MassSingular = 'yes';
    options.InitialSlope = model.odefun(tspan(1),X0,par,model);
end

% jacobian, if provided
if ~isempty(model.jacfun)
    options.Jacobian = @(t,X) model.jacfun(t,X,par,model);
end

% event to stop simulation, if provided
if ~isempty(model.event)
    options.Events = @model.event;
end

% simulate model ODE
[t,X,te] = ode15s(@(t,X) model.odefun(t,X,par,model), tspan, X0, options);
if ~isempty(model.event) && t(end)~=tspan(end)
    t=[t;tspan(end)];
    X(end+1,:)=0;
end
if isfield(model,'multiple') && model.multiple
    input_events=model.input_events;
    input_doses=model.input_doses;
    u_ref=model.u_ref;
    t_out=[];
    x_out=[];
    te_out=[];
    timespan=tspan;
    X0_ref=X0(:)+u_ref;
   
    if timespan(1)<input_events(1)
        if length(timespan)==2
            tspan = [timespan(1) input_events(1)];
        else
            tspan = timespan(timespan<input_events(1));
        end
        [t_out,x_out,te_out] = ode15s(@(t,X) model.odefun(t,X,par,model), tspan, X0, options);
        X0_ref=x_out(end,:)' + u_ref;
    end
    for i=2:length(input_events)
        if input_events(i)<=timespan(end) && timespan(1) < input_events(i)
            %%% in case of specification of specific timepoints to be solved need
            %%% to account for this
            if length(timespan)==2
                tspan = [input_events(i-1) input_events(i)];
            else
                tspan = timespan(timespan>=input_events(i-1) & timespan<=input_events(i));
            end
            [t,x,te]  = ode15s(@(t,X) model.odefun(t,X,par,model), tspan, X0_ref, options);
            %check_whether_ODE_solver_reported_warnings()
            %%% account for perturbed dosing in the case of calculation of indices
            %%% to avoid double timepoints in the time vector in the multiple
            %%% dosing case, the dosing timepoint is kept and the other disgarded
            %%% such that in the state vector only the dosing state is included
            if i<length(input_events)
                t_out = [t_out;t(1:end-1)];
                x_out = [x_out;x(1:end-1,:)];
                te_out = [te_out;te(1:end-1)];
            else
                if timespan(end)>input_events(end)
                    t_out = [t_out;t(1:end-1)];
                    x_out = [x_out;x(1:end-1,:)];
                    te_out = [te_out;te(1:end-1)];
                else
                    t_out = [t_out;t];
                    x_out = [x_out;x];
                    te_out = [te_out;te];
                end
            end
            if any(input_doses)
                u_ref(I.input)=input_doses(i);
            end
            X0_ref=x(end,:)' + u_ref;
        elseif timespan(1)>=input_events(i)
            continue;
        else
            break;
        end
    end

    if timespan(end)>input_events(end)
        if length(timespan)==2
            tspan = [input_events(end) timespan(end)];
        else
            tspan = timespan(timespan>=input_events(end));
        end
        [t,x,te]  = ode15s(@(t,X) model.odefun(t,X,par,model), tspan, X0_ref, options);
        t_out = [t_out;t];
        x_out = [x_out;x];
        te_out = [te_out;te];
    end

    X=x_out;
    t=t_out;
    te=te_out; 
%    model.Y_red = model.I.h_red(model.t_red,model.X_red,te,model);

end

% determine value of states eliminated via conservation laws
for p = 1:length(I.con)
    k = I.con(p);  % index of state variable
    c = L.con2conlaw(p); % number of conservation law
    
    % consistency check
    if ~ismember(k,L.statesofconlaw{c})
        fprintf(['\n\n --> The state ''%s'' does not belong to ',I.nmstate{k},...
            'the conservation law specified by ''%s''',L.conlawspec{c},':-(  Please fix \n\n']);
        error('--> Ending here')
    end
    
    % indices of remaining states that are part of the conlaw
    remstates = L.remstatesofconlaw{p}; 
    X(:,k) = max(0, L.valconlaw(c) - sum(X(:,remstates),2) );
    
end
end

function [position, isterminal, direction] = event_fct(~,concentrations,limit,index,direction)
% event function: stop simulation when index'th concentration >= limit 
% (<= if direction=-1)
position = concentrations(index)-limit;
isterminal = 1;
end
