%%% Version: 19 Jan 2023
%%%
%%% model  = <MODELNAME>_model_set_up_details(model)
%%%
%%% This function specifies all required details that are needed to finally
%%% set of the model
%%% 
%%% Input:  model       model structure (containing the fields model.name 
%%%                     and model.scenario)
%%%                   
%%% Output: model       updated model structure
%%%
%%% Author: Jane Knoechel and Wilhelm Huisinga
%%%

function model = MMEnzymeKinetics_model_set_up_details(model)

%%% define naming for result files 
model.projectfolder = ''; %%% option to provide an absolute path here
model.savenameroot = [model.projectfolder 'results/' model.name];
if ~isempty(model.scenario)
    model.savenameroot = [model.savenameroot '_' model.scenario];
end

%%% for testing purposes, reduced number of computations (i.e., timepoints,
%%% states) when computing the indices
model.quicktest = false;

% load indexing to be able to define input, output, conslaws etc
I = feval([model.name '_indexing']);

fprintf('\n  --> Scenario: %s\n',model.scenario);

% define input
model.setup.u_ref      = 50; % [nM]
% define output
model.setup.unit.graphic.ylim   = [1e-2 5e2]; % in output units

switch model.scenario
    case 'nothing'
        % define input
        model.setup.input  = I.A;
        model.setup.Etot   = 2; % [nM]

        %%% set parameter values
        model.setup.kon    = 50;  % [1/nM/min]
        model.setup.koff   = 100; % [1/min]
        model.setup.kcat   = 1;  % [1/min]
        model.setup.ktrans = 0.1;  % [1/min]
        
    case 'basic_irek'
        % define input
        model.setup.input  = I.S;
        model.setup.Etot   = 2; % [nM]

        %%% set parameter values
        model.setup.kon    = 50;  % [1/nM/min]
        model.setup.koff   = 100; % [1/min]
        model.setup.kcat   = 1;  % [1/min]
        model.setup.kdeg   = 0;  % [1/min]
        model.setup.ktrans = 0;  % [1/min]

    case 'basic_irek_enzyme_excess'
        % define input
        model.setup.input  = I.S;
        model.setup.Etot   = 100; % [nM]

        %%% set parameter values
        model.setup.kon    = 50;  % [1/nM/min]
        model.setup.koff   = 100; % [1/min]
        model.setup.kcat   = 1;  % [1/min]
        model.setup.kdeg   = 0;  % [1/min]
        model.setup.ktrans = 0;  % [1/min]

    case 'Cpss_Eenv'
        % define input
        model.setup.input  = I.A;
        model.setup.Etot   = 25; % [nM]

        %%% set parameter values
        model.setup.kon    = 1;  % [1/nM/min]
        model.setup.koff   = 900; % [1/min]
        model.setup.kcat   = 5;  % [1/min]
        model.setup.ktrans = 0.5;  % [1/min]

    case 'Cpss_EpCenv'
        % define input
        model.setup.input  = I.A;
        model.setup.Etot   = 2; % [nM]

        %%% set parameter values
        model.setup.kon    = 0.2;  %5 [1/nM/min]
        model.setup.koff   = 10; % [1/min]
        model.setup.kcat   = 1;  % [1/min]
        model.setup.ktrans = 0.1;  % [1/min]

    case 'Spss_CpS'

        fac = 5;
        % define input
        model.setup.input  = I.S;
        model.setup.Etot   = 5000; % [nM]

        %%% set parameter values
        model.setup.kon    = 100/fac;% [1/nM/min]
        model.setup.koff   = 1/fac; % [1/min]
        model.setup.kcat   = 1/fac;  % [1/min]
        model.setup.ktrans = 30/fac;  % [1/min]
        
        model.setup.unit.graphic.ylim   = [1e-8 2e5]; % in output units

    otherwise
        fprintf('\n\n  --> unknown scenario---PLEASE FIX! \n\n');
        error(':-)');

end
model.setup.unit.input = 'nM';

% define output
model.setup.output      = I.P;
model.setup.unit.output = 'nM';
model.setup.unit.graphic.transfoutput = @(x) x;
model.setup.unit.graphic.output = 'nM';

% simulation time span and unit transformation (here: sec to min)
model.setup.tspan     = [0 30];  % in [min]
model.setup.unit.time = 'min';
model.setup.unit.graphic.transftime = @(x) x;
model.setup.unit.graphic.time = 'min';
model.setup.unit.graphic.xlim = model.setup.tspan + [-1 1]; % in output units


% define observed state variables (to be plottet) and plot
model.setup.obsstates = 1:I.nstates;

% conservation laws for this small example were set manually
model.setup.conlaw.n = 3;
model.setup.conlaw.EC.states   = [I.E I.C]; 
model.setup.conlaw.ASCP.states = [I.A I.S I.C I.P]; 
model.setup.conlaw.SC.states   = [I.S I.C];  %%% only for testing purposes !!!


% optional: jacobian of ODE and jacobian pattern of extended ODE 
% false = no jacobian provided 
% true = jacobian provided 
model.setup.jacfunprovided = true;

% optional: define your specific legend labels (if false, default values are used)
% use the file 'default_legendlabels.m' in the folder 'generalfiles' as a
% template to define your own legend labels, if desired
model.setup.legendlabelsprovided = false;

% optional: define your specific line styles (if false, default values are used)
% use the file 'default_state2linestyle.m' in the folder 'generalfiles' as a
% template to define your own legend labels, if desired
model.setup.linestyleprovided = false;

