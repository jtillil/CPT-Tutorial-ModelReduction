%%% Version: 19 Jun 2022
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
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function par = Wajima2009BloodCoagulation_parameters(model)

%%% assign model indexing
I = model.I;

%%% initializing parameter vector
par = NaN(I.npar,1);

%%% -----------------------------------------------------------------------
%%% parameters for Hillfunctions and complex formation
%%% Vmax is denoted by v in [1/h] and Km is denoted by k in [nM]
%%% c denotes the rate constant of complex formation and is given in [(nM*h)]

%%% r1, r2 in VIII
par(I.v1) = 50000;           par(I.k1) = 0.000001; 
par(I.v2) = 50;              par(I.k2) = 1;

%%% r3 in IX
par(I.v3) = 7;               par(I.k3) = 10;

%%% r4,5 in XI
par(I.v4) = 7;               par(I.k4) = 1;
par(I.v5) = 10;              par(I.k5) = 10;

%%% r6 in VII
par(I.v6) = 0.1;             par(I.k6) = 10;

%%% r7,8,9 in X
par(I.v7) = 0.02;            par(I.k7) = 10;
par(I.v8) = 2;               par(I.k8) = 0.1;
par(I.v9) = 0.000000001;     par(I.k9) = 10;

%%% r10,11 in V
par(I.v10) = 50000;          par(I.k10) = 10;
%%% adaptations by Gulati et al for the in-vivo scenario
if contains(model.scenario,'in_vivo')
    par(I.v10)  =  25000;        par(I.k10)  =  1800;
end

par(I.v11) = 50;             par(I.k11) = 1;

%%% r12,13 in II
par(I.v12) = 100;            par(I.k12) = 10;
par(I.v13) = 9;              par(I.k13) = 500;

%%% r14,15,16,17,18,19 in I 
par(I.v14) = 20000;          par(I.k14) = 0.5;
par(I.v15) = 500;            par(I.k15) = 500;
%%% adaptations by Gulati et al for the in-vivo scenario
if contains(model.scenario,'in_vivo')
    par(I.v14) = 21000;          par(I.k14) = 30000;
end

par(I.v16) = 7;              par(I.k16) = 10;
par(I.v17) = 7;              par(I.k17) = 10;
par(I.v18) = 7;              par(I.k18) = 100;
par(I.v19) = 1;              par(I.k19) = 1;

%%% r20 XIII
par(I.v20) = 7;              par(I.k20) = 1;

%%% r21,22,23 in Pg
par(I.v21) = 7;              par(I.k21) = 5000;
par(I.v22) = 5;              par(I.k22) = 10000;
par(I.v23) = 2;              par(I.k23) = 1;

%%% r24 in PC
par(I.v24) = 7;              par(I.k24) = 1;

%%% r25 in Xa_Va
par(I.v25) = 2;              par(I.k25) = 1;

%%% r26 in VIII and IXa_VIIIa
par(I.c26) = 0.01;

%%% r27 in X
par(I.c27) = 0.5;

%%% r28 in II and IIa_Tmod
par(I.c28) = 0.5;

%%% r29,30 in VII
par(I.c29) = 0.5;
par(I.c30) = 0.1;

%%% r31 in TF and VIIa_TF_Xa_TFPI and Xa_TFPI
par(I.c31) = 0.5;

%%% r32 in X and Xa_TFPI and TFPI
par(I.c32) = 0.5;

%%% r33 in TF
par(I.v33) = 70;             par(I.k33) = 1;

%%% r34 in X
par(I.v34) = 900;            par(I.k34) = 200;

%%% r35 in IX
par(I.v35) = 70;             par(I.k35) = 1;

%%% r36 in TF
par(I.v36) = 1000;           par(I.k36) = 1;

%%% r37 in PC and PS and APC_PS
par(I.c37) = 0.5;

%%% r38,39,40 in VII
par(I.v38) = 1;              par(I.k38) = 10;
par(I.v39) = 1;              par(I.k39) = 10;
par(I.v40) = 0.2;            par(I.k40) = 10;

%%% r.41,r.42 in Hageman Factor XII
par(I.v41) = 7;              par(I.k41) = 1;
par(I.v42) = 70;             par(I.k42) = 1;

%%% r43 in Pk
par(I.v43) = 7;              par(I.k43) = 1;


%%% -----------------------------------------------------------------------
%%% parameters for drug/venom/VK dynamics
%%%

%%% -----------------------------------------------------------------------
%%% vitamin k parameters
%%%
par(I.VK_k12)  = 0.0587;
par(I.VK_k21)  = 0.0122;
par(I.VK_V)    = 24;

%%% -----------------------------------------------------------------------
%%% brown snake venom parameters
%%%
par(I.ka_Brown)   =  5.0; % absorption rate constant
par(I.d_Brown)    =  3.5; % degradation rate constant

%%% -----------------------------------------------------------------------
%%% tiger snake venom parameters
%%%
par(I.ka_Tiger)   =  5.0; % absorption rate constant
par(I.d_Tiger)    =  3.5; % degradation rate constant

%%% -----------------------------------------------------------------------
%%% taipan snake venom parameters
%%% 
par(I.d_Taipan)      =  3.5; 
par(I.ktrans_Taipan) = 0.99;

%%% -----------------------------------------------------------------------
%%% taipan snake venom action on blood coagulation parameters
%%% 
par(I.vtaipan)      = 70;  par(I.ktaipan)      =  10; 

%%% -----------------------------------------------------------------------
%%% warfarin PK parameters
%%% 
par(I.ka_Warf) = 1.0;
par(I.Vd_Warf) = 10;
par(I.Cl_Warf) = 0.2;
par(I.ke_Warf) = par(I.Cl_Warf)/par(I.Vd_Warf);

%%% -----------------------------------------------------------------------
%%% warfarin PD parameters
%%% 
par(I.lmax) = 1;
par(I.IC50) = 0.34;

%%% -----------------------------------------------------------------------
%%% heparin PK parameters (unfractioned heparin UFH or low-molecular weight heparin LMWH)
%%%
par(I.ke_Hep)  = 0.693;
par(I.ka_Hep)  = 0.255;
par(I.Vc_Hep)  = 4.567;
par(I.Vp_Hep)  = 29.6;
par(I.Cl_Hep)  = 1.058;
par(I.Q_Hep)   = 0.62;
par(I.ke_Hep)  = par(I.Cl_Hep)/par(I.Vc_Hep);
par(I.k12_Hep) = par(I.Q_Hep)/par(I.Vc_Hep);
par(I.k21_Hep) = par(I.Q_Hep)/par(I.Vp_Hep);

par(I.R1) = 1/7.1; %UFH
%par(I.R1) = 3.9; %enoxaparin

%%% r44 in II
par(I.c44) = 0.85 * par(I.R1);

%%% r45 in X
par(I.c45) = 0.85;
par(I.R2) = 1; %UFH

%%%par(I.R2) = 10; %LMWH
% r46 in IX
par(I.c46) = par(I.c45) * par(I.R2);

par(I.inf_rate_UFH) = 0;%32000/24;

%%% -----------------------------------------------------------------------
%%% degradation rate parameters
%%% All degradation rate constants are in the unit [1/h].
%%%
par(I.degFg) = 0.032;             par(I.degF) = 0.050;          par(I.degXF) = 0.050;             par(I.degII) = 0.010;
par(I.degIIa) = 67.4;             par(I.degTF) = 0.05;          par(I.degV) = 0.043;              par(I.degVa) = 20.0;
par(I.degVII) = 0.12;             par(I.degVIIa) = 20.0;        par(I.degVIII) = 0.058;           par(I.degVIIIa) = 20.0;
par(I.degIX) = 0.029;             par(I.degIXa) = 20.0;         par(I.degX) = 0.018;              par(I.degXa) = 20.0;
par(I.degXI) = 0.10;              par(I.degXIa) = 20.0;         par(I.degXII) = 0.012;            par(I.degXIIa) = 20.0;
par(I.degXIII) = 0.0036;          par(I.degXIIIa) = 0.69;       par(I.degPk) = 0.05;              par(I.degK) = 20.0;
par(I.degPg) = 0.05;              par(I.degP) = 20.0;           par(I.degPC) = 0.050;             par(I.degAPC) = 20.4;
par(I.degPS) = 0.0165;            par(I.degFDP) = 3.5;          par(I.degD) = 0.1;                par(I.degTFPI) = 20.0;
par(I.degVIIaTF) = 20.0;          par(I.degVIITF) = 0.7;        par(I.degAPCPS) = 20.0;           par(I.degXaVa) = 20.0;
par(I.degIXaVIIIa) = 20.0;        par(I.degTmod) = 0.050;       par(I.degIIaTmod) = 20.0;         par(I.degXaTFPI) = 20.0;
par(I.degVIIaTFXaTFPI) = 20.0;    par(I.degTAT) = 0.2;          par(I.degCA) = 0.05;              par(I.degVK) = 0.2052;
par(I.degVK2) = 0.0228;           

%%% degVKH2 and degVKO depend on the initial concentration
%%%
X0 = model.X0prior2input;

par(I.degVKH2) = par(I.degVK2) * (X0(I.VK)./(X0(I.VKH2)));
par(I.degVKO)  = par(I.degVK2) * (X0(I.VK)./(X0(I.VKO)));

%%% to avoid numerical issues if either VKH2 or VK are considered as
%%% neglected species degradation rate constants are assumed to be zero
if isfield(I,'neg') && (ismember('VKH2',I.neg) || ismember('VK',I.neg))
    par(I.degVKH2) = 0;
    par(I.degVKO)  = 0;
end

%%% -----------------------------------------------------------------------
%%% production rate parameters
%%%

par(I.aII)      = par(I.degII)  * X0(I.II)  / (X0(I.VKH2));
par(I.aVII)     = par(I.degVII) * X0(I.VII) / (X0(I.VKH2)); 
par(I.aIX)      = par(I.degIX)  * X0(I.IX)  / (X0(I.VKH2));
par(I.aX)       = par(I.degX)   * X0(I.X)   / (X0(I.VKH2));
par(I.aPC)      = par(I.degPC)  * X0(I.PC)  / (X0(I.VKH2)); 
par(I.aPS)      = par(I.degPS)  * X0(I.PS)  / (X0(I.VKH2)); 

par(I.pV)       = par(I.degV)    * X0(I.V); 
par(I.pVIII)    = par(I.degVIII) * X0(I.VIII);
par(I.pFg)      = par(I.degFg)   * X0(I.Fg);  
par(I.pXIII)    = par(I.degXIII) * X0(I.XIII);
par(I.pPg)      = par(I.degPg)   * X0(I.Pg);   
par(I.pTmod)    = par(I.degTmod) * X0(I.Tmod);
par(I.pXI)      = par(I.degXI)   * X0(I.XI);  
par(I.pTFPI)    = par(I.degTFPI) * X0(I.TFPI);
par(I.pVK)      = par(I.degVK)   * X0(I.VK);  
par(I.pXII)     = par(I.degXII)  * X0(I.XII);
par(I.pPk)      = par(I.degPk)   * X0(I.Pk);


%%% in the case of the in vitro setting (PT or aPTT) all production rates are assumed to
%%% be zero (see Wajima et al)
if contains(model.scenario,'in_vitro')
        par(I.aII) = 0;    par(I.aVII)  = 0;
        par(I.aIX) = 0;    par(I.aX)    = 0;
        par(I.aPC) = 0;    par(I.aPS)   = 0;
        par(I.pV)  = 0;    par(I.pVIII) = 0;
        par(I.pFg) = 0;    par(I.pXIII) = 0;
        par(I.pPg) = 0;    par(I.pTmod) = 0;
        par(I.pXI) = 0;    par(I.pTFPI) = 0;
        par(I.pVK) = 0;    par(I.pXII)  = 0;
        par(I.pPk) = 0;
end


end

