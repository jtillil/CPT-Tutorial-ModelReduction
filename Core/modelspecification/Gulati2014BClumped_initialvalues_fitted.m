%%% Version: July 11th, 2022
%%%
%%% X0 = Wajima2009BloodCoagulation_initialValues(model)
%%%
%%% This function assigns the initial concentrations for all species 
%%% of the blood coagulation model 
%%%
%%% Input :  model      model structure containing the index structure of
%%%                     the model
%%%
%%% Output : X0         initial values 
%%%
%%% 
%%% Authors: Undine Falkenhagen, Jane Knoechel and Wilhelm Huisinga
%%%

function X0 = Gulati2014BClumped_initialvalues_fitted(model)

%%% assign model indexing
I  = model.I;

%%% initialise state values
X0 = NaN(I.nstates,1);

X0(I.AVenom)    = 0.0075;
X0(I.CVenom)    = 0;
% X0(I.lumped)    = 7900;
X0(I.lumped)    = 1394.4;
X0(I.IIa)       = 0;
% X0(I.Fg)        = 8900;
X0(I.Fg)        = 8945.5;

end
