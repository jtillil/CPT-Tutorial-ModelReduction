%%% Version: July 11th, 2022
%%%
%%% dX  = Wajima2009BloodCoagulation_ode(t,X,par,model)
%%%
%%% This function defines the ode system of the  blood coagulation model
%%% 
%%% Input : t           time
%%%         X           state vector
%%%         par         parameter vector
%%%         model       model structure containing the index structure of
%%%                     the model
%%%                   
%%% Output : dX         right hand side of ODEs
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
%%% Author: Undine Falkenhagen, Jane Knoechel and Wilhelm Huisinga
%%%

function dX = Gulati2014BClumped_ode(~,X,par,model)

%%% assign model indexing
I  = model.I;

%%% -----------------------------------------------------------------------
%%% specify system of ODEs 
dX = 0*X;

%%% -----------------------------------------------------------------------
%%% AVenom, CVenom
abs = par(I.kaVenom)*X(I.AVenom);
d = par(I.dVenom)*X(I.CVenom);
dX(I.AVenom) = - abs;
dX(I.CVenom) = abs - d;

%%% -----------------------------------------------------------------------
%%% lumped
p = par(I.plumped);
d = par(I.dlumped)*X(I.lumped);
act = par(I.vlumpedVenom2IIa)*X(I.CVenom)*X(I.lumped) / ...
    (par(I.klumpedVenom2IIa) + X(I.CVenom));
dX(I.lumped) = p - d - act;

%%% -----------------------------------------------------------------------
%%% IIa
d = par(I.dIIa)*X(I.IIa);
dX(I.IIa) = act - d;

%%% -----------------------------------------------------------------------
%%% Fg
p = par(I.pFg);
d = par(I.dFg)*X(I.Fg);
catd = par(I.vFgIIa2F)*X(I.IIa)*X(I.Fg) / ...
    (par(I.kFgIIa2F) + X(I.IIa));
dX(I.Fg) = p - d - catd;

end