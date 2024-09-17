%%% Version: 19 Jan 2023
%%%
%%% par  =  <MODELNAME>_parameters(model)
%%%
%%% This function creates a structure with all parameters
%%%
%%% Input :  model      model structure containing the index structure of
%%%                     the model and its initial values 
%%%                      
%%% Output : par        parameter vector
%%%
%%%
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function par = MMEnzymeKinetics_parameters(model)

%%% assign model indexing
I = model.I;

%%% initializing parameter vector
par = NaN(I.npar,1);

par(I.ktrans) = model.setup.ktrans;  % [1/min]
par(I.kon)    = model.setup.kon;  % [1/nM/min]
par(I.koff)   = model.setup.koff; % [1/min]
par(I.kcat)   = model.setup.kcat;  % [1/min]

par(I.kdeg) = 0; % [1/min]
        

















