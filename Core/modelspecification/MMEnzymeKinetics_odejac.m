%%% Version: 19 Jan 2023
%%%
%%% jac  = <MODELNAME>_odejac(t,X,par,model)
%%%
%%% This function defines the jacobian of the ode system 
%%%
%%% Input   t           time
%%%         X           state vector
%%%         par         parameter vector
%%%         model       model structure containing the index structure of
%%%                     the model
%%%                   
%%%
%%% Output  jac         Jacobian of the right hand side with respect to  
%%%                     the state variables
%%%
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function DF = MMEnzymeKinetics_odejac(~,X,par,model)

%%% assign model indexing
I  = model.I;

%%% initialize the jacobian
DF  = 0*X;

%%% -----------------------------------------------------------------------
%%% define jacobian

DF(I.A,I.A)   = - par(I.ktrans);

DF(I.S,I.A)   =   par(I.ktrans);
DF(I.S,I.S)   = - par(I.kon) * X(I.E);
DF(I.S,I.E)   = - par(I.kon) * X(I.S);
DF(I.S,I.C)   =   par(I.koff);

DF(I.E,I.S)   = - par(I.kon) * X(I.E);
DF(I.E,I.E)   = - par(I.kon) * X(I.S);
DF(I.E,I.C)   =   par(I.koff) + par(I.kcat);

DF(I.C,I.S)   =   par(I.kon) * X(I.E);
DF(I.C,I.E)   =   par(I.kon) * X(I.S);
DF(I.C,I.C)   = -(par(I.koff) + par(I.kcat));

DF(I.P,I.C)   =   par(I.kcat);
DF(I.P,I.P)   = - par(I.kdeg);

%%% -----------------------------------------------------------------------


end
