%%% Version: 19 Jan 2023
%%%
%%% dX  = <MODELNAME>_ode(t,X,par,model)
%%%
%%% This function defines the system of ODEs 
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

function dX = SimpleParallelPathways_ode(~,X,par,model)

%%% assign model indexing
I  = model.I;

%%% initialize rhs vector
dX = 0*X;

%%% -----------------------------------------------------------------------
%%% specify system of ODEs 

dX(I.A)     = 0;
dX(I.S)     = -( par(I.k_ab)+par(I.k_ac) )*X(I.A)*X(I.S) + par(I.k_b)*X(I.B) + par(I.k_c)*X(I.C);
dX(I.B)     = par(I.k_ab)*X(I.A)*X(I.S) - par(I.k_b)*X(I.B);
dX(I.C)     = par(I.k_ac)*X(I.A)*X(I.S) - par(I.k_c)*X(I.C); 
dX(I.D)     = par(I.k_bd)*X(I.B) + par(I.k_cd)*X(I.C) - par(I.k_d)*X(I.D); 

%%% -----------------------------------------------------------------------

end

