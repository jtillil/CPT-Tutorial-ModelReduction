%%% Version: 19 Jan 2023
%%%
%%% dX  = <MODELNAME>_ode(t,X,par,model)
%%%
%%% This function defines the system of model ODEs 
%%% 
%%% Input : t           time
%%%         X           state vector
%%%         par         parameter vector
%%%         model       model structure containing the index structure of
%%%                     the model
%%%                   
%%% Output : dX         right hand side of ODEs
%%%
%%%
%%% Author: Jane Knoechel and Wilhelm Huisinga
%%%

function dX = MMEnzymeKinetics_ode(~,X,par,model)

%%% assign model indexing
I  = model.I;

%%% initialize rhs vector
dX = 0*X;

%%% -----------------------------------------------------------------------
%%% specify system of ODEs 

dX(I.A)  = -par(I.ktrans) * X(I.A);
dX(I.S)   = par(I.ktrans) * X(I.A) -par(I.kon) * X(I.S) * X(I.E) + par(I.koff) * X(I.C);
dX(I.E)   = -par(I.kon) * X(I.S) * X(I.E) + ( par(I.koff) + par(I.kcat) ) * X(I.C);
dX(I.C)   =  par(I.kon) * X(I.S) * X(I.E) - ( par(I.koff) + par(I.kcat) ) * X(I.C);
dX(I.P)   = par(I.kcat) * X(I.C) - par(I.kdeg) * X(I.P); 

%%% -----------------------------------------------------------------------

end

