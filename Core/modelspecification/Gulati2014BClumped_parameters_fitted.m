%%% Version: July 11th, 2022
%%%
%%% par = Wajima2009BloodCoagulation_parameters(model)
%%%
%%% This function creates a structure with all parameters for
%%% the blood coagulation model
%%%
%%% Input :  model      model structure containing the index structure of
%%%                     the model and its initial values X0
%%%                      
%%% Output : par        parameter vector
%%%
%%% References:
%%% + Wajima et al. "A Comprehensive Model for the Humoral
%%%   Coagulation Network in Humans " (2009), Supplement material Figure 2
%%%   the presented adjusted parameters were used 
%%% + Gulati et al. "Effect of Australian elapid venom on
%%%   blood coagulation: Australian Snakebite Project (ASP-17)", (2013)
%%%   parameters given in Table 1 and 2 
%%% + Tanos et al. "A model for venom-induced consumptive coagulopathy in
%%%   snake bite", (2008) 
%%%   parameters given in Table 1 and 2
%%%
%%% 
%%% Authors: Undine Falkenhagen, Jane Knoechel and Wilhelm Huisinga
%%%

function par = Gulati2014BClumped_parameters_fitted(model)

%%% assign model indexing
I = model.I;

%%% initializing parameter vector
par = NaN(I.npar,1);

par(I.dIIa)             = 58.3;
par(I.vlumpedVenom2IIa) = 1940;
par(I.klumpedVenom2IIa) = 1.38;
par(I.dFg)              = 0.0174;
par(I.pFg)              = 216;
par(I.vFgIIa2F)         = 8.88;
par(I.kFgIIa2F)         = 22.9;
par(I.dVenom)           = 0.744;
par(I.kaVenom)          = 0.854;
par(I.plumped)          = 147;
par(I.dlumped)          = 0.01;

end

