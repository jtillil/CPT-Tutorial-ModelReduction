%%% Version: 19 June 2022
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

function model = Wajima2009BloodCoagulation_model_set_up_details(model)

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
switch model.scenario
    case {'in_vivo_snakevenom_1h','in_vivo_snakevenom_40h','in_vivo_snakevenom_50h'}

        %%% define input
        model.setup.input      = I.AVenom;
        amount_snake_venom     = 0.0015; % amount in mg
        SF_mg_to_nmol          = 1e-3/2e5*1e9;
        model.setup.u_ref      = SF_mg_to_nmol*amount_snake_venom;
        model.setup.unit.input = 'nmol'; % needs to be in nmol due to model ODEs

        %%% define output
        model.setup.output      = I.Fg;
        model.setup.unit.output = 'nM'; % needs to be in nM due to model ODEs
        model.setup.unit.graphic.transfoutput = @(x) x;
        model.setup.unit.graphic.output = 'nM';
        model.setup.unit.graphic.ylim   = [1e-4 2e4]; % in output units

        %%% simulation time span and unit transformation
        model.setup.unit.time = 'h'; % needs to be in h due to model ODEs
        if strcmp(model.scenario, 'in_vivo_snakevenom_1h')
            model.setup.tspan = [0 1]; 
            model.setup.unit.graphic.xlim = [-0.02 1.02]; % in output units
        elseif strcmp(model.scenario, 'in_vivo_snakevenom_40h')
            model.setup.tspan = [0 40];
            model.setup.unit.graphic.xlim = [-1 41]; % in output units
        elseif strcmp(model.scenario, 'in_vivo_snakevenom_50h')
            model.setup.tspan = [0 50];
            model.setup.unit.graphic.xlim = [-1 51]; % in output units
        end
        model.setup.unit.graphic.transftime = @(x) x;
        model.setup.unit.graphic.time = 'h';

        %%% define observed state variables (to be plottet) and plot
        model.setup.obsstates = [I.P,I.Fg,I.IIa,I.APC_PS,I.IIa_Tmod,I.APC,I.AVenom,I.CVenom];

    case {'in_vitro_PTtest_highTF','in_vitro_PTtest_lowTF'}

        %%% define input
        model.setup.input      = I.TF;
        model.setup.unit.input = 'nM'; % needs to be in nM due to model ODEs
        if strcmp(model.scenario, 'in_vitro_PTtest_highTF')
            model.setup.u_ref  = 100;  % in nM
        elseif strcmp(model.scenario, 'in_vitro_PTtest_lowTF')
            model.setup.u_ref  = 5e-3; % in nM
        end

        %%% define output
        model.setup.output      = I.F;
        model.setup.unit.output = 'nM'; % needs to be in nM due to model ODEs
        model.setup.unit.graphic.transfoutput = @(x) x;
        model.setup.unit.graphic.output = 'nM';
        model.setup.unit.graphic.ylim   = [1e-5 1e6]; % in output units

        %%% simulation time span and unit transformation (not applicable here)
        model.setup.unit.time = 'h'; % needs to be in h due to model ODEs
        if strcmp(model.scenario, 'in_vitro_PTtest_highTF')
            model.setup.tspan = [0 30/3600]; % needs to be in h due to model ODEs
            model.setup.unit.graphic.xlim = [-1 31]; % in output units
            model.setup.obsstates = [I.TF, I.VII_TF, I.VIIa_TF, I.Xa, I.IIa, I.Fg, I.F, I.VII, I.X, I.II, I.Xa_Va I.P];
        elseif strcmp(model.scenario, 'in_vitro_PTtest_lowTF')
            model.setup.tspan = [0 240/3600]; % needs to be in h due to model ODEs
            model.setup.unit.graphic.xlim = [-1 241]; % in output units
            model.setup.obsstates = [I.TF,I.VIII,I.VIIIa,I.VII_TF,I.VIIa_TF,I.IXa,I.IXa_VIIIa,I.Xa,I.Va,I.IIa,I.Xa_Va,I.Fg,I.F];
        end
        model.setup.unit.graphic.transftime = @(x) 3600*x;
        model.setup.unit.graphic.time = 'sec';

    otherwise
        fprintf('\n\n  --> unknown scenario---PLEASE FIX! \n\n');
        error(':-)');

end

% conservation laws: none considered
model.setup.conlaw.n = 0;

% optional: jacobian of ODE and jacobian pattern of extended ODE 
% false = no jacobian provided 
% true = jacobian provided 
model.setup.jacfunprovided = true;

% optional: define your specific legend labels (if false, default values are used)
% use the file 'default_legendlabels.m' in the folder 'generalfiles' as a
% template to define your own legend labels, if desired
model.setup.legendlabelsprovided = true;

% optional: define your specific line styles (if false, default values are used)
% use the file 'default_state2linestyle.m' in the folder 'generalfiles' as a
% template to define your own legend labels, if desired
model.setup.linestyleprovided = false;




