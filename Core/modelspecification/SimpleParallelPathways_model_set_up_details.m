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

function model = SimpleParallelPathways_model_set_up_details(model)

%%% define naming for result files 
model.projectfolder = ''; %%% option to provide an absolute path here
model.savenameroot = [model.projectfolder 'results/' model.name];
if ~isempty(model.scenario)
    model.savenameroot = [model.savenameroot '_' model.scenario];
end

switch model.scenario
    case 'no_crosstalk'
        fac = 300;
        model.setup.S     = 10*fac; 
        model.setup.k_ab  = 50/fac;
        model.setup.k_b  = 5;
        
        model.setup.k_ac  = 1/fac;
        model.setup.k_c  = 500;

        model.setup.k_bd  = 0.4;
        model.setup.k_cd  = 11;

        model.setup.tspan = [0 0.05];  % in [min]
    case 'with_crosstalk'
        model.setup.S     = 10; 
        model.setup.k_ab  = 50;
        model.setup.k_b  = 5;
        
        model.setup.k_ac  = 1;
        model.setup.k_c  = 500;

        model.setup.k_bd  = 0.4;
        model.setup.k_cd  = 11;

        model.setup.tspan = [0 0.05];  % in [min]
    case 'C_B_neglectable'
        fac = 300;
        % fac = 10;
        model.setup.S     = 10*fac; 
        % model.setup.S     = 10; 
        model.setup.k_ab  = 10/fac;
        % model.setup.k_ab  = 10;
        model.setup.k_b  = 300;
        
        model.setup.k_ac  = 10/fac;
        % model.setup.k_ac  = 10;
        model.setup.k_c  = 300;

        model.setup.k_bd  = 2;
        model.setup.k_cd  = 2;

        model.setup.tspan = [0 0.05];  % in [min]
    otherwise
        fprintf('\n\n  --> Scenarios %s currently not supported!---PLEASE FIX! \n\n',model.scenario); error(' :-)');
end
%%% for testing purposes, reduced number of computations (i.e., timepoints,
%%% states) when computing the indices
model.quicktest = false;

% load indexing to be able to define input, output, conslaws etc
I = feval([model.name '_indexing']);

% define input
model.setup.input      = I.A;
model.setup.u_ref      = 20;
model.setup.unit.input = 'nM';

% define output 
model.setup.output      = I.D;
model.setup.unit.output = 'nM';
model.setup.unit.graphic.transfoutput = @(x) x;
model.setup.unit.graphic.output = 'nM';
model.setup.unit.graphic.ylim   = [1e-3 1e4]; % in output units

% simulation time span and unit transformation (here: sec to min)
model.setup.unit.time = 'min';
model.setup.unit.graphic.transftime = @(x) x;
model.setup.unit.graphic.time = 'min';
Delta_tspan = model.setup.tspan(2)-model.setup.tspan(1);
model.setup.unit.graphic.xlim = model.setup.tspan + [-1 1]*Delta_tspan*0.05; % in output units

% define observed state variables (to be plottet) and plot
model.setup.obsstates = 1:I.nstates;

% conservation law for this small example was set manually
model.setup.conlaw.n = 1;
model.setup.conlaw.SBC.states = [I.S I.B I.C]; 

% optional: jacobian of ODE and jacobian pattern of extended ODE 
% [] = no jacobian provided 
% [1] = jacobian provided 
model.setup.jacfunprovided = false;

% optional: define your specific legend labels (if false, default values are used)
% use the file 'default_legendlabels.m' in the folder 'generalfiles' as a
% template to define your own legend labels, if desired
model.setup.legendlabelsprovided = false;

% optional: define your specific line styles (if false, default values are used)
% use the file 'default_state2linestyle.m' in the folder 'generalfiles' as a
% template to define your own legend labels, if desired
model.setup.linestyleprovided = false;



