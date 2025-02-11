%%% Version: 12 Jul 2023
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

function DF = Wajima2009BloodCoagulation_odejac(~,X,par,model)

%%% assign model indexing
I  = model.I;

%%% initialize the jacobian
% DF  = spalloc(I.nstates,I.nstates,3*I.nstates);
sz = size(X);
if sz(1) == 1
    DF = repmat(0*X, [length(X) 1]);
elseif sz(2) == 1
    DF = repmat(0*X, [1 length(X)]);
else
    disp(sz)
    error('Wrong size of X for jac fun.')
end

%%% -----------------------------------------------------------------------
%%% define jacobian

%%% -----------------------------------------------------------------------
%%% states: XII, XIIa
%%% eqs. 1 and 2: Hageman- Factor (XII) and activation
%%%
DF(I.XII,I.XII)   = -( par(I.v41)*X(I.CA) ) / ( par(I.k41)+X(I.CA) ) - ( par(I.v42)*X(I.K))/(par(I.k42)+X(I.K)) - par(I.degXII);
DF(I.XII,I.CA)    = -(par(I.v41)*par(I.k41))/(par(I.k41)+X(I.CA))^2*X(I.XII);
DF(I.XII,I.K)     = -(par(I.v42)*par(I.k42))/(par(I.k42)+X(I.K))^2*X(I.XII);

DF(I.XIIa,I.XII)  = (par(I.v41)*X(I.CA))/(par(I.k41)+X(I.CA))+(par(I.v42)*X(I.K))/(par(I.k42)+X(I.K));
DF(I.XIIa,I.CA)   = (par(I.v41)*par(I.k41))/(par(I.k41)+X(I.CA))^2*X(I.XII);
DF(I.XIIa,I.K)    = (par(I.v42)*par(I.k42))/(par(I.k42)+X(I.K))^2*X(I.XII);
DF(I.XIIa,I.XIIa) = -par(I.degXIIa);

%%% -----------------------------------------------------------------------
%%% states: VIII, VIIIa
%%% eqs. 3 and 4 : Antihaemophiles Globulin A (VIII) and activation
%%%
DF(I.VIII,I.VIII)    = -(par(I.v1)*X(I.IIa))/(par(I.k1)+X(I.IIa))-par(I.degVIII);
DF(I.VIII,I.IIa)     = -(par(I.v1)*par(I.k1))/(par(I.k1)+X(I.IIa))^2*X(I.VIII);

DF(I.VIIIa,I.VIII)   = (par(I.v1)*X(I.IIa))/(par(I.k1)+X(I.IIa));
DF(I.VIIIa,I.IIa)    = (par(I.v1)*par(I.k1))/(par(I.k1)+X(I.IIa))^2*X(I.VIII);
DF(I.VIIIa,I.VIIIa)  = -(par(I.v2)*X(I.APC_PS))/(par(I.k2)+X(I.APC_PS))-(X(I.IXa))/par(I.c26)-par(I.degVIIIa);
DF(I.VIIIa,I.IXa)    = -(X(I.VIIIa))/par(I.c26);
DF(I.VIIIa,I.APC_PS) = -(par(I.v2)*par(I.k2))/(par(I.k2)+X(I.APC_PS))^2*X(I.VIIIa);

%%% -----------------------------------------------------------------------
%%% states: IX, IXa
%%% eqs. 5 and 6 : Antihaemophiles Globulin B (IX) and activation
%%%
DF(I.IX,I.IX)       = -(par(I.v3)*X(I.XIa))/(par(I.k3)+X(I.XIa))-(par(I.v35)*X(I.VIIa_TF))/(par(I.k35)+X(I.VIIa_TF))-par(I.degIX);
DF(I.IX,I.XIa)      = -(par(I.v3)*par(I.k3))/(par(I.k3)+X(I.XIa))^2*X(I.IX);
DF(I.IX,I.VIIa_TF)  = -(par(I.v35)*par(I.k35))/(par(I.k35)+X(I.VIIa_TF))^2*X(I.IX);
DF(I.IX,I.VKH2)     = par(I.aIX);

DF(I.IXa,I.IX)      = (par(I.v3)*X(I.XIa))/(par(I.k3)+X(I.XIa))+(par(I.v35)*X(I.VIIa_TF))/(par(I.k35)+X(I.VIIa_TF));
DF(I.IXa,I.IXa)     = -(X(I.AT_III_Heparin))/par(I.c46)-(X(I.VIIIa))/par(I.c26)-par(I.degIXa);
DF(I.IXa,I.VIIIa)   = -(X(I.IXa))/par(I.c26);
DF(I.IXa,I.AT_III_Heparin) = -(X(I.IXa))/par(I.c46);
DF(I.IXa,I.XIa)     = (par(I.v3)*par(I.k3))/(par(I.k3)+X(I.XIa))^2*X(I.IX);
DF(I.IXa,I.VIIa_TF) = (par(I.v35)*par(I.k35))/(par(I.k35)+X(I.VIIa_TF))^2*X(I.IX);

%%% -----------------------------------------------------------------------
%%% states: XI, XIa
%%% eqs. 7 and 8 : Rosenthal-factor (XI) and activation
%%%
DF(I.XI,I.XI)    = -(par(I.v4)*X(I.XIIa))/(par(I.k4)+X(I.XIIa))-(par(I.v5)*X(I.IIa))/(par(I.k5)+X(I.IIa))-par(I.degXI);
DF(I.XI,I.XIIa)  = -(par(I.v4)*par(I.k4))/(par(I.k4)+X(I.XIIa))^2*X(I.XI);
DF(I.XI,I.IIa)   = -(par(I.v5)*par(I.k5))/(par(I.k5)+X(I.IIa))^2*X(I.XI);

DF(I.XIa,I.XI)   = (par(I.v4)*X(I.XIIa))/(par(I.k4)+X(I.XIIa))+(par(I.v5)*X(I.IIa))/(par(I.k5)+X(I.IIa));
DF(I.XIa,I.XIIa) = (par(I.v4)*par(I.k4))/(par(I.k4)+X(I.XIIa))^2*X(I.XI);
DF(I.XIa,I.IIa)  = (par(I.v5)*par(I.k5))/(par(I.k5)+X(I.IIa))^2*X(I.XI);
DF(I.XIa,I.XIa)  = -par(I.degXIa);

%%% -----------------------------------------------------------------------
%%% states: VII, VIIa
%%% eqs. 9 and 10: Proconvertin (VII) and activation
%%%
DF(I.VII,I.VII)      = -(X(I.TF))/(par(I.c30))-(par(I.v6)*X(I.IIa))/(par(I.k6)+X(I.IIa))...
                        -(par(I.v38)*X(I.Xa))/(par(I.k38)+X(I.Xa))...
                        -(par(I.v39)*X(I.VIIa_TF))/(par(I.k39)+X(I.VIIa_TF))...
                        -(par(I.v40)*X(I.IXa))/(par(I.k40)+X(I.IXa))-par(I.degVII);
DF(I.VII,I.TF)       = -(X(I.VII))/(par(I.c30));
DF(I.VII,I.IIa)      = -(par(I.v6)*par(I.k6))/(par(I.k6)+X(I.IIa))^2*X(I.VII);
DF(I.VII,I.Xa)       = -(par(I.v38)*par(I.k38))/(par(I.k38)+X(I.Xa))^2*X(I.VII);
DF(I.VII,I.VIIa_TF)  = -(par(I.v39)*par(I.k39))/(par(I.k39)+X(I.VIIa_TF))^2*X(I.VII);
DF(I.VII,I.IXa)      = -(par(I.v40)*par(I.k40))/(par(I.k40)+X(I.IXa))^2*X(I.VII);
DF(I.VII,I.VKH2)     = par(I.aVII);
DF(I.VII,I.delayTaipan2) = - ( par(I.vtaipan)*(par(I.ktaipan)+X(I.delayTaipan2)) - par(I.vtaipan)*X(I.delayTaipan2)  ) / ... 
                            (par(I.ktaipan)+X(I.delayTaipan2))^2 * X(I.VII); %%% added 14 Sep 2021

DF(I.VIIa,I.VII)     = (par(I.v6)*X(I.IIa))/(par(I.k6)+X(I.IIa))...
                        +(par(I.v38)*X(I.Xa))/(par(I.k38)+X(I.Xa))...
                        +(par(I.v39)*X(I.VIIa_TF))/(par(I.k39)+X(I.VIIa_TF))...
                        +(par(I.v40)*X(I.IXa))/(par(I.k40)+X(I.IXa));
DF(I.VIIa,I.TF)      = -(X(I.VIIa))/(par(I.c29));
DF(I.VIIa,I.IIa)     = (par(I.v6)*par(I.k6))/(par(I.k6)+X(I.IIa))^2*X(I.VII);
DF(I.VIIa,I.Xa)      = (par(I.v38)*par(I.k38))/(par(I.k38)+X(I.Xa))^2*X(I.VII);
DF(I.VIIa,I.VIIa_TF) = (par(I.v39)*par(I.k39))/(par(I.k39)+X(I.VIIa_TF))^2*X(I.VII);
DF(I.VIIa,I.IXa)     = (par(I.v40)*par(I.k40))/(par(I.k40)+X(I.IXa))^2*X(I.VII);
DF(I.VIIa,I.VIIa)    = -(X(I.TF))/(par(I.c29))-par(I.degVIIa);

%%% -----------------------------------------------------------------------
%%% eqs. 11 and 12: Stuart-Power-Factor (X) and activation
%%%
DF(I.X,I.X)          = -(par(I.v7)*X(I.IXa))/(par(I.k7)+X(I.IXa))-(par(I.v8)*X(I.IXa_VIIIa))/(par(I.k8)+X(I.IXa_VIIIa))...
                        -(par(I.v9)*X(I.VIIa))/(par(I.k9)+X(I.VIIa))-(par(I.v34)*X(I.VIIa_TF))/(par(I.k34)+X(I.VIIa_TF))-par(I.degX);
DF(I.X,I.IXa)        = -(par(I.v7)*par(I.k7))/(par(I.k7)+X(I.IXa))^2*X(I.X);
DF(I.X,I.IXa_VIIIa)  = -(par(I.v8)*par(I.k8))/(par(I.k8)+X(I.IXa_VIIIa))^2*X(I.X);
DF(I.X,I.VIIa)       = -(par(I.v9)*par(I.k9))/(par(I.k9)+X(I.VIIa))^2*X(I.X);
DF(I.X,I.VIIa_TF)    = -(par(I.v34)*par(I.k34))/(par(I.k34)+X(I.VIIa_TF))^2*X(I.X);
DF(I.X,I.VKH2)       = par(I.aX);

DF(I.Xa,I.X)         = (par(I.v7)*X(I.IXa))/(par(I.k7)+X(I.IXa))+(par(I.v8)*X(I.IXa_VIIIa))/(par(I.k8)+X(I.IXa_VIIIa))...
                        +(par(I.v9)*X(I.VIIa))/(par(I.k9)+X(I.VIIa))+(par(I.v34)*X(I.VIIa_TF))/(par(I.k34)+X(I.VIIa_TF));
DF(I.Xa,I.IXa)       = (par(I.v7)*par(I.k7))/(par(I.k7)+X(I.IXa))^2*X(I.X);
DF(I.Xa,I.IXa_VIIIa) = (par(I.v8)*par(I.k8))/(par(I.k8)+X(I.IXa_VIIIa))^2*X(I.X);
DF(I.Xa,I.VIIa)      = (par(I.v9)*par(I.k9))/(par(I.k9)+X(I.VIIa))^2*X(I.X);
DF(I.Xa,I.VIIa_TF)   = (par(I.v34)*par(I.k34))/(par(I.k34)+X(I.VIIa_TF))^2*X(I.X);
DF(I.Xa,I.Xa)        = -(X(I.Va))/(par(I.c27))-(X(I.TFPI))/(par(I.c32))-((X(I.AT_III_Heparin)))/par(I.c45)-par(I.degXa);
DF(I.Xa,I.Va)        = -(X(I.Xa))/(par(I.c27));
DF(I.Xa,I.TFPI)      = -(X(I.Xa))/(par(I.c32));
DF(I.Xa,I.AT_III_Heparin) = -((X(I.Xa)))/par(I.c45);

%%% -----------------------------------------------------------------------
%%% eqs. 13 and 14 : Proaccelerin (V) and activation
%%%
DF(I.V,I.V)       = -(par(I.v10)*X(I.IIa))/(par(I.k10)+X(I.IIa))-par(I.degV);
DF(I.V,I.IIa)     = -(par(I.v10)*par(I.k10))/(par(I.k10)+X(I.IIa))^2*X(I.V);

DF(I.Va,I.V)      = (par(I.v10)*X(I.IIa))/(par(I.k10)+X(I.IIa));
DF(I.Va,I.IIa)    = (par(I.v10)*par(I.k10))/(par(I.k10)+X(I.IIa))^2*X(I.V);
DF(I.Va,I.Va)     = -(par(I.v11)*X(I.APC_PS))/(par(I.k11)+X(I.APC_PS))-(X(I.Xa))/(par(I.c27))-par(I.degVa);
DF(I.Va,I.Xa)     = -(X(I.Va))/(par(I.c27));
DF(I.Va,I.APC_PS) = -(par(I.v11)*par(I.k11))/(par(I.k11)+X(I.APC_PS))^2*X(I.Va);

%%% -----------------------------------------------------------------------
%%% eq. 15 : Complex Xa_Va
%%%
DF(I.Xa_Va,I.Xa_Va)  = -(par(I.v25)*X(I.APC_PS))/(par(I.k25)+X(I.APC_PS))-par(I.degXaVa);
DF(I.Xa_Va,I.Xa)     = (X(I.Va))/(par(I.c27));
DF(I.Xa_Va,I.Va)     = (X(I.Xa))/(par(I.c27));
DF(I.Xa_Va,I.APC_PS) = -(par(I.v25)*par(I.k25))/(par(I.k25)+X(I.APC_PS))^2*X(I.Xa_Va);

%%% -----------------------------------------------------------------------
%%% eqs. 16 and 17 : Prothrombin (II) and activation
%%%

DF(I.II,I.II)           = -(par(I.v12)* (X(I.Xa_Va)+X(I.CVenom)+X(I.TaipanVenom) ))/(par(I.k12)+(X(I.Xa_Va)+X(I.CVenom)+X(I.TaipanVenom)))-(par(I.v13)*(X(I.Xa)+X(I.CVenom_Tiger)))/(par(I.k13)+(X(I.Xa)+X(I.CVenom_Tiger)))-par(I.degII);
DF(I.II,I.Xa_Va)        = -(par(I.v12)*par(I.k12))/(par(I.k12)+X(I.Xa_Va)+X(I.CVenom)+X(I.TaipanVenom))^2*X(I.II);
DF(I.II,I.CVenom)       = -(par(I.v12)*par(I.k12))/(par(I.k12)+X(I.Xa_Va)+X(I.CVenom)+X(I.TaipanVenom))^2*X(I.II);
DF(I.II,I.TaipanVenom)  = -(par(I.v12)*par(I.k12))/(par(I.k12)+X(I.Xa_Va)+X(I.CVenom)+X(I.TaipanVenom))^2*X(I.II);
DF(I.II,I.Xa)           = -(par(I.v13)*par(I.k13))/(par(I.k13)+X(I.Xa)+X(I.CVenom_Tiger))^2*X(I.II);
DF(I.II,I.CVenom_Tiger) = -(par(I.v13)*par(I.k13))/(par(I.k13)+X(I.Xa)+X(I.CVenom_Tiger))^2*X(I.II);
DF(I.II,I.VKH2)         = par(I.aII);

DF(I.IIa,I.II)              = (par(I.v12)*(X(I.Xa_Va)+X(I.CVenom)+X(I.TaipanVenom)))/(par(I.k12)+X(I.Xa_Va)+X(I.CVenom)+X(I.TaipanVenom))+(par(I.v13)*(X(I.Xa)+X(I.CVenom_Tiger)))/(par(I.k13)+(X(I.Xa)+X(I.CVenom_Tiger)));
DF(I.IIa,I.Xa_Va)           = (par(I.v12)*par(I.k12))/(par(I.k12)+X(I.Xa_Va)+X(I.CVenom)+X(I.TaipanVenom))^2*X(I.II);
DF(I.IIa,I.CVenom)          = (par(I.v12)*par(I.k12))/(par(I.k12)+X(I.Xa_Va)+X(I.CVenom)+X(I.TaipanVenom))^2*X(I.II);
DF(I.IIa,I.TaipanVenom)     = (par(I.v12)*par(I.k12))/(par(I.k12)+X(I.Xa_Va)+X(I.CVenom)+X(I.TaipanVenom))^2*X(I.II);
DF(I.IIa,I.Xa)              = (par(I.v13)*par(I.k13))/(par(I.k13)+X(I.Xa)+X(I.CVenom_Tiger))^2*X(I.II);
DF(I.IIa,I.CVenom_Tiger)    = (par(I.v13)*par(I.k13))/(par(I.k13)+X(I.Xa)+X(I.CVenom_Tiger))^2*X(I.II);
DF(I.IIa,I.IIa)             = -(X(I.Tmod))/(par(I.c28))-((X(I.AT_III_Heparin)))/par(I.c44)-par(I.degIIa);
DF(I.IIa,I.Tmod)            = -(X(I.IIa))/(par(I.c28));
DF(I.IIa,I.AT_III_Heparin)  = -(X(I.IIa))/par(I.c44);

%%% -----------------------------------------------------------------------
%%% eq. 18: thrombin-antithrombin (TAT)
%%%

DF(I.TAT,I.IIa) = par(I.degIIa);
DF(I.TAT,I.TAT) = -par(I.degTAT);

%%% -----------------------------------------------------------------------
%%% eqs. 19, 20, 21, 22 and 23 : Fibrinogen (I), activation and complex/linked
%%%
DF(I.Fg,I.Fg)     = -(par(I.v14)*X(I.IIa))/(par(I.k14)+X(I.IIa))-(par(I.v15)*X(I.P))/(par(I.k15)+X(I.P))-par(I.degFg);
DF(I.Fg,I.IIa)    = -(par(I.v14)*par(I.k14))/(par(I.k14)+X(I.IIa))^2*X(I.Fg);
DF(I.Fg,I.P)      = -(par(I.v15)*par(I.k15))/(par(I.k15)+X(I.P))^2*X(I.Fg);

DF(I.F,I.Fg)      = (par(I.v14)*X(I.IIa))/(par(I.k14)+X(I.IIa));
DF(I.F,I.F)       = -(par(I.v16)*X(I.XIIIa))/(par(I.k16)+X(I.XIIIa))-(par(I.v17)*X(I.P))/(par(I.k17)+X(I.P))-par(I.degF);
DF(I.F,I.IIa)     = (par(I.v14)*par(I.k14))/(par(I.k14)+X(I.IIa))^2*X(I.Fg);
DF(I.F,I.XIIIa)   = -(par(I.v16)*par(I.k16))/(par(I.k16)+X(I.XIIIa))^2*X(I.F);
DF(I.F,I.P)       = -(par(I.v17)*par(I.k17))/(par(I.k17)+X(I.P))^2*X(I.F);

DF(I.XF,I.F)      = (par(I.v16)*X(I.XIIIa))/(par(I.k16)+X(I.XIIIa));
DF(I.XF,I.XF)     = -(par(I.v18)*X(I.P))/(par(I.k18)+X(I.P))-(par(I.v19)*X(I.APC_PS))/(par(I.k19)+X(I.APC_PS))-par(I.degXF);
DF(I.XF,I.P)      = -(par(I.v18)*par(I.k18))/(par(I.k18)+X(I.P))^2*X(I.XF);
DF(I.XF,I.APC_PS) = -(par(I.v19)*par(I.k19))/(par(I.k19)+X(I.APC_PS))^2*X(I.XF);
DF(I.XF,I.XIIIa)  = (par(I.v16)*par(I.k16))/(par(I.k16)+X(I.XIIIa))^2*X(I.F);

DF(I.FDP,I.Fg)    = (par(I.v15)*X(I.P))/(par(I.k15)+X(I.P))+par(I.degFg);
DF(I.FDP,I.F)     = (par(I.v17)*X(I.P))/(par(I.k17)+X(I.P))+par(I.degF);
DF(I.FDP,I.FDP)   = -par(I.degFDP);
DF(I.FDP,I.P)     = (par(I.v15)*par(I.k15))/(par(I.k15)+X(I.P))^2*X(I.Fg)+(par(I.v17)*par(I.k17))/(par(I.k17)+X(I.P))^2*X(I.F);

DF(I.D,I.D)       = -par(I.degD);
DF(I.D,I.XF)      = (par(I.v18)*X(I.P))/(par(I.k18)+X(I.P))+(par(I.v19)*X(I.APC_PS))/(par(I.k19)+X(I.APC_PS))+par(I.degXF);
DF(I.D,I.P)       = (par(I.v18)*par(I.k18))/(par(I.k18)+X(I.P))^2*X(I.XF);
DF(I.D,I.APC_PS)  = (par(I.v19)*par(I.k19))/(par(I.k19)+X(I.APC_PS))^2*X(I.XF);

%%% -----------------------------------------------------------------------
%%% eqs. 24 and 25 : Fibrinstabilisierender Factor (XIII) and activation
%%%
DF(I.XIII,I.XIII)   = -(par(I.v20)*X(I.IIa))/(par(I.k20)+X(I.IIa))-par(I.degXIII);
DF(I.XIII,I.IIa)    = -(par(I.v20)*par(I.k20))/(par(I.k20)+X(I.IIa))^2*X(I.XIII);

DF(I.XIIIa,I.IIa)   = (par(I.v20)*par(I.k20))/(par(I.k20)+X(I.IIa))^2*X(I.XIII);
DF(I.XIIIa,I.XIII)  = (par(I.v20)*X(I.IIa))/(par(I.k20)+X(I.IIa));
DF(I.XIIIa,I.XIIIa) = -par(I.degXIIIa);

%%% -----------------------------------------------------------------------
%%% eqs. 26 and 27: Plasminogen (Pg) and Plasmin (P)
%%%
DF(I.Pg,I.Pg)     = -(par(I.v21)*X(I.IIa))/(par(I.k21)+X(I.IIa))-(par(I.v22)*X(I.F))/(par(I.k22)+X(I.F))-(par(I.v23)*X(I.APC_PS))/(par(I.k23)+X(I.APC_PS))-par(I.degPg);
DF(I.Pg,I.IIa)    = -(par(I.v21)*par(I.k21))/(par(I.k21)+X(I.IIa))^2*X(I.Pg);
DF(I.Pg,I.F)      = -(par(I.v22)*par(I.k22))/(par(I.k22)+X(I.F))^2*X(I.Pg);
DF(I.Pg,I.APC_PS) = -(par(I.v23)*par(I.k23))/(par(I.k23)+X(I.APC_PS))^2*X(I.Pg);

DF(I.P,I.Pg)      = (par(I.v21)*X(I.IIa))/(par(I.k21)+X(I.IIa))+(par(I.v22)*X(I.F))/(par(I.k22)+X(I.F))+(par(I.v23)*X(I.APC_PS))/(par(I.k23)+X(I.APC_PS));
DF(I.P,I.IIa)     = (par(I.v21)*par(I.k21))/(par(I.k21)+X(I.IIa))^2*X(I.Pg);
DF(I.P,I.F)       = (par(I.v22)*par(I.k22))/(par(I.k22)+X(I.F))^2*X(I.Pg);
DF(I.P,I.APC_PS ) = (par(I.v23)*par(I.k23))/(par(I.k23)+X(I.APC_PS))^2*X(I.Pg);
DF(I.P,I.P)       = -par(I.degP);

%%% -----------------------------------------------------------------------
%%% eqs. 28 and 29 : Protein C (PC) and activation (APC)
%%%
DF(I.PC,I.PC)        = -(par(I.v24)*X(I.IIa_Tmod))/(par(I.k24)+X(I.IIa_Tmod))-par(I.degPC);
DF(I.PC,I.IIa_Tmod)  = -(par(I.v24)*par(I.k24))/(par(I.k24)+X(I.IIa_Tmod))^2*X(I.PC);
DF(I.PC,I.VKH2)      = par(I.aPC);

DF(I.APC,I.PC)       = (par(I.v24)*X(I.IIa_Tmod))/(par(I.k24)+X(I.IIa_Tmod));
DF(I.APC,I.IIa_Tmod) = (par(I.v24)*par(I.k24))/(par(I.k24)+X(I.IIa_Tmod))^2*X(I.PC);
DF(I.APC,I.APC)      = -(X(I.PS))/(par(I.c37))-par(I.degAPC);
DF(I.APC,I.PS)       = -(X(I.APC))/(par(I.c37));

%%% -----------------------------------------------------------------------
%%% eq. 30 : Thrombmodulin (Tmod)
%%%
DF(I.Tmod,I.Tmod) = -(X(I.IIa))/(par(I.c28))-par(I.degTmod);
DF(I.Tmod,I.IIa)  = -(X(I.Tmod))/(par(I.c28));

%%% -----------------------------------------------------------------------
%%% eq. 31 : IIa_Tmod complex
%%%
DF(I.IIa_Tmod,I.Tmod)     = (X(I.IIa))/(par(I.c28));
DF(I.IIa_Tmod,I.IIa)      = (X(I.Tmod))/(par(I.c28));
DF(I.IIa_Tmod,I.IIa_Tmod) = -par(I.degIIaTmod);

%%% -----------------------------------------------------------------------
%%% eq. 32 : IXa_VIIIa complex
%%%
DF(I.IXa_VIIIa,I.IXa)       = (X(I.VIIIa))/par(I.c26);
DF(I.IXa_VIIIa,I.VIIIa)     = (X(I.IXa))/par(I.c26);
DF(I.IXa_VIIIa,I.IXa_VIIIa) = -par(I.degIXaVIIIa);

%%% -----------------------------------------------------------------------
%%% eqs. 33, 34 and 35 : tissue factor also called III (TF) and 
%%% its complexes (VII_TF and VIIa_TF)
%%%
DF(I.TF,I.TF)           = -(X(I.VIIa))/(par(I.c29))-(X(I.VII))/(par(I.c30))-par(I.degTF);
DF(I.TF,I.VII)          = -(X(I.TF))/(par(I.c30));
DF(I.TF,I.VIIa)         = -(X(I.TF))/(par(I.c29));

DF(I.VII_TF,I.VII_TF)   = -(par(I.v33)*X(I.Xa))/(par(I.k33)+X(I.Xa))-(par(I.v36)*X(I.TF))/(par(I.k36)+X(I.TF))-par(I.degVIITF);
DF(I.VII_TF,I.VII)      = (X(I.TF))/(par(I.c30));
DF(I.VII_TF,I.TF)       = (X(I.VII))/(par(I.c30))-(par(I.v36)*par(I.k36))/(par(I.k36)+X(I.TF))^2*X(I.VII_TF);
DF(I.VII_TF,I.Xa)       = -(par(I.v33)*par(I.k33))/(par(I.k33)+X(I.Xa))^2*X(I.VII_TF);

DF(I.VIIa_TF,I.Xa)      = (par(I.v33)*par(I.k33))/(par(I.k33)+X(I.Xa))^2*X(I.VII_TF);
DF(I.VIIa_TF,I.TF)      = (X(I.VIIa))/(par(I.c29))+(par(I.v36)*par(I.k36))/(par(I.k36)+X(I.TF))^2*X(I.VII_TF);
DF(I.VIIa_TF,I.VII_TF)  = (par(I.v33)*X(I.Xa))/(par(I.k33)+X(I.Xa))+(par(I.v36)*X(I.TF))/(par(I.k36)+X(I.TF));
DF(I.VIIa_TF,I.VIIa_TF) = -par(I.degVIIaTF)-(X(I.Xa_TFPI))/(par(I.c31));
DF(I.VIIa_TF,I.Xa_TFPI) = -(X(I.VIIa_TF))/(par(I.c31));
DF(I.VIIa_TF,I.VIIa)    = (X(I.TF))/(par(I.c29));

%%% -----------------------------------------------------------------------
%%% eqs. 36, 37 and 38 :tissue factor pathway inhibitor (TFPI) and its
%%% complexes (Xa_TFPI and VIIa_TF_Xa_TFPI)
%%%
DF(I.TFPI,I.TFPI)       = -(X(I.Xa))/(par(I.c32))-par(I.degTFPI);
DF(I.TFPI,I.Xa)         = -(X(I.TFPI))/(par(I.c32));

DF(I.Xa_TFPI,I.Xa_TFPI) = -(X(I.VIIa_TF))/(par(I.c31))-par(I.degXaTFPI);
DF(I.Xa_TFPI,I.Xa)      = (X(I.TFPI))/(par(I.c32));
DF(I.Xa_TFPI,I.TFPI)    = (X(I.Xa))/(par(I.c32));
DF(I.Xa_TFPI,I.VIIa_TF) = -(X(I.Xa_TFPI))/(par(I.c31));

DF(I.VIIa_TF_Xa_TFPI,I.VIIa_TF)         = (X(I.Xa_TFPI))/(par(I.c31));
DF(I.VIIa_TF_Xa_TFPI,I.Xa_TFPI)         = (X(I.VIIa_TF))/(par(I.c31));
DF(I.VIIa_TF_Xa_TFPI,I.VIIa_TF_Xa_TFPI) = -par(I.degVIIaTFXaTFPI);

%%% -----------------------------------------------------------------------
%%% eqs. 39 and 40 Protein S (PS) and complex (APC_PS)
%%%
DF(I.PS,I.PS)   = -(X(I.APC))/(par(I.c37))-par(I.degPS);
DF(I.PS,I.APC)  = -(X(I.PS))/(par(I.c37));
DF(I.PS,I.VKH2) = par(I.aPS);

DF(I.APC_PS,I.APC)    = (X(I.PS))/(par(I.c37));
DF(I.APC_PS,I.PS)     = (X(I.APC))/(par(I.c37));
DF(I.APC_PS,I.APC_PS) = -par(I.degAPCPS);

%%% -----------------------------------------------------------------------
%%% eqs. 41 and 42 :Prekallikrein (PK) and Kallikrein (K)
%%%
DF(I.Pk,I.Pk)   = -(par(I.v43)*X(I.XIIa))/(par(I.k43)+X(I.XIIa))-par(I.degPk);
DF(I.Pk,I.XIIa) = -(par(I.v43)*par(I.k43))/(par(I.k43)+X(I.XIIa))^2*X(I.Pk);

DF(I.K,I.XIIa)  = (par(I.v43)*par(I.k43))/(par(I.k43)+X(I.XIIa))^2*X(I.Pk);
DF(I.K,I.Pk)    = (par(I.v43)*X(I.XIIa))/(par(I.k43)+X(I.XIIa));
DF(I.K,I.K)     = -par(I.degK);

%%% -----------------------------------------------------------------------
%%% eqs. 43,44,45 and 46 : Vitamin K (VK)
%%%

DF(I.VK,I.VK)     = -par(I.degVK)-par(I.VK_k12)-par(I.degVK2)*(1-(par(I.lmax)*X(I.Cwarf)))/(par(I.IC50)+X(I.Cwarf));
DF(I.VK,I.VKO)    = par(I.degVKO)*(1-(par(I.lmax)*X(I.Cwarf))/(par(I.IC50)+X(I.Cwarf)));
DF(I.VK,I.VK_p)   = (par(I.VK_k21)/(par(I.VK_V)));

DF(I.VKH2,I.VK)   = par(I.degVK2)*(1-(par(I.lmax)*X(I.Cwarf))/(par(I.IC50)+X(I.Cwarf)));
DF(I.VKH2,I.VKH2) = -par(I.degVKH2);

DF(I.VKO,I.VKO)   = -par(I.degVKO)*(1-(par(I.lmax)*X(I.Cwarf))/(par(I.IC50)+X(I.Cwarf)));
DF(I.VKO,I.VKH2)  = par(I.degVKH2);

DF(I.VK_p,I.VK)   = par(I.VK_k12)*par(I.VK_V);
DF(I.VK_p,I.VK_p) = -par(I.VK_k21);

%%% -----------------------------------------------------------------------
%%% eqs. 47 and 48 for warfarin PK 

DF(I.Awarf,I.Awarf) = -par(I.ka_Warf);

DF(I.Cwarf,I.Awarf) = par(I.ka_Warf)/par(I.Vd_Warf);
DF(I.Cwarf,I.Cwarf) = -par(I.ke_Warf);

%%% -----------------------------------------------------------------------
%%% eq. 49 :activator for the contact system (CA)
%%%
DF(I.CA,I.CA) = -par(I.degCA);

%%% -----------------------------------------------------------------------
%%% eqs. 50, 51 and 52 for Heparin
%%%
DF(I.AEnox,I.AEnox) = -par(I.ka_Hep);

DF(I.AT_III_Heparin,I.AEnox) = par(I.ka_Hep)/par(I.Vc_Hep);
DF(I.AT_III_Heparin,I.AT_III_Heparin) = -par(I.k12_Hep)-par(I.ke_Hep)-X(I.IIa)/par(I.c44)-X(I.Xa)/par(I.c45)-X(I.IXa)/par(I.c46);
DF(I.AT_III_Heparin,I.ENO_p) = par(I.k21_Hep)*par(I.Vc_Hep);
DF(I.AT_III_Heparin,I.IIa)   = X(I.AT_III_Heparin)/par(I.c44);
DF(I.AT_III_Heparin,I.Xa)    = X(I.AT_III_Heparin)/par(I.c45);
DF(I.AT_III_Heparin,I.IXa)   = X(I.AT_III_Heparin)/par(I.c46);

DF(I.ENO_p,I.AT_III_Heparin) = par(I.k12_Hep)*par(I.Vc_Hep);
DF(I.ENO_p,I.ENO_p)          = - par(I.k21_Hep);

%%% -----------------------------------------------------------------------
%%% eqs. 53 AUC of Fibrin
%%%
DF(I.AUC,I.F) = 1;
%%% -----------------------------------------------------------------------
%%% eqs.54 and 55 for brown snake venom with absorption compartment (A) and
%%% concentration compartment (C)
%%%

DF(I.AVenom,I.AVenom) = -par(I.ka_Brown);

DF(I.CVenom,I.AVenom) = par(I.ka_Brown);
DF(I.CVenom,I.CVenom) = -par(I.d_Brown);

%%% -----------------------------------------------------------------------
%%% eqs. 56,57 and 58  : Taipan snake venom & delay compartment for action of Taipan
%%% on VII
%%%

DF(I.TaipanVenom,I.TaipanVenom) =-par(I.d_Taipan);

DF(I.delayTaipan1,I.TaipanVenom) =par(I.d_Taipan);
DF(I.delayTaipan1,I.delayTaipan1) =- par(I.ktrans_Taipan);

DF(I.delayTaipan2,I.delayTaipan1) = par(I.ktrans_Taipan);
DF(I.delayTaipan2,I.delayTaipan2) = -par(I.ktrans_Taipan);

%%% -----------------------------------------------------------------------
%%% eqs. 59  : ATIII
%%%
DF(I.ATIII,I.ATIII)=0;

%%% -----------------------------------------------------------------------
%%% eqs. 60 and 61 : tiger snake venom with absorption compartment (A) and
%%% concentration compartment (C)
%%%
DF(I.AVenom_Tiger,I.AVenom_Tiger) = -par(I.ka_Tiger);

DF(I.CVenom_Tiger,I.AVenom_Tiger) = par(I.ka_Tiger);
DF(I.CVenom_Tiger,I.CVenom_Tiger) = -par(I.d_Tiger);

%%% -----------------------------------------------------------------------
%%% eqs. 62 : Herparin (UFH) PK
%%%
DF(I.AT_III_UFH,I.AT_III_UFH) = -par(I.ke_Hep)-X(I.IIa)/par(I.c44)-X(I.Xa)/par(I.c45)-X(I.IXa)/par(I.c46);
DF(I.AT_III_UFH,I.IIa)   = -X(I.AT_III_UFH)/par(I.c44);
DF(I.AT_III_UFH,I.Xa)    = -X(I.AT_III_UFH)/par(I.c45);
DF(I.AT_III_UFH,I.IXa)   = -X(I.AT_III_UFH)/par(I.c46);

%%% -----------------------------------------------------------------------

%%% check if there are NaN or Inf entries
if any(isinf(DF),'all') || any(isnan(DF),'all')
    error('\n There is an issue with the Jacobian (contains ''inf'' or ''NaN'' entries).')
end

end
