%%% Version: 19 Jan 2023
%%%
%%% par = <MODELNAME>_parameters(model)
%%%
%%% Creates a structure with all parameters value
%%%
%%% Input :  model      model structure containing the index structure of
%%%                     the model and its initial values 
%%%                      
%%% Output : par        parameter vector
%%%
%%% 
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function par = SimpleParallelPathways_parameters(model)

%%% assign model indexing
I = model.I;

%%% initializing parameter vector
par = NaN(I.npar,1);

%%% set parameter values
par(I.k_d)  = 25;

% switch model.scenario
%     case {'no_crosstalk','with_crosstalk'}
        par(I.k_ab) = model.setup.k_ab;
        par(I.k_b)  = model.setup.k_b;
        par(I.k_ac) = model.setup.k_ac;
        par(I.k_c)  = model.setup.k_c;
        par(I.k_bd) = model.setup.k_bd;
        par(I.k_cd) = model.setup.k_cd;
% end


end
