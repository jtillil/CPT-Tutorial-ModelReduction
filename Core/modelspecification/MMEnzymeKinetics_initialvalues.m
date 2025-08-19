%%% Version: 19 Jan 2023
%%%
%%% X0 = <MODELNAME>_initialValues(model)
%%%
%%% This function assigns the initial values for all state variables 
%%%
%%% Input :  model      model structure containing the index structure of
%%%                     the model
%%%
%%% Output : X0         initial values 
%%%
%%%
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%
    
function X0 = MMEnzymeKinetics_initialvalues(model)

%%% assign model indexing
I  = model. I;

%%% initialise state values
X0 = NaN(I.nstates,1);

X0(I.A) = 0; % [nM] input
X0(I.S) = 0; 
X0(I.C) = 0;
X0(I.P) = 0;

%%% set state values. NOTE: the input (here, the number of A molecules) is 
%%% defined in the MAIN script
X0(I.E) = model.setup.Etot;

end