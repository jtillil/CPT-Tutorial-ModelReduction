%%% Version: January 24th, 2020
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
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function X0 = Wajima2009BloodCoagulation_initialvalues(model)

%%% assign model indexing
I  = model. I;

%%% initialise state values
X0 = NaN(I.nstates,1);

X0(I.XII) = 375;      X0(I.XIIa) = 0;               X0(I.VIII) = 0.70;        X0(I.VIIIa) = 0;
X0(I.IX) = 89.6;      X0(I.IXa) = 0;                X0(I.XI) = 30.6;          X0(I.XIa) = 0;
X0(I.VII) = 10.0;     X0(I.VIIa) = 0;               X0(I.X) = 174.3;          X0(I.Xa) = 0;
X0(I.V) = 26.7;       X0(I.Va) = 0;                 X0(I.Xa_Va) = 0;          X0(I.II) = 1394.4;
X0(I.IIa) = 0;        X0(I.TAT) = 0;                X0(I.Fg) = 8945.5;        X0(I.F) = 0;
X0(I.XF) = 0;         X0(I.FDP) = 0;                X0(I.D) = 0;              X0(I.XIII) = 70.3;
X0(I.XIIIa) = 0;      X0(I.Pg) = 2154.3;            X0(I.P) = 0;              X0(I.PC) = 60.0;
X0(I.APC) = 0;        X0(I.Tmod) = 50.0;            X0(I.IIa_Tmod) = 0;       X0(I.IXa_VIIIa) = 0;
X0(I.TF) = 0;         X0(I.VII_TF) = 0;             X0(I.VIIa_TF) = 0;        X0(I.TFPI) = 2.5;
X0(I.Xa_TFPI) = 0;    X0(I.VIIa_TF_Xa_TFPI) = 0;    X0(I.PS) = 300;           X0(I.APC_PS) = 0;
X0(I.Pk) = 450;       X0(I.K) = 0;                  X0(I.VK) = 1.0;           X0(I.VKH2) = 0.1;
X0(I.VKO) = 0.1;      X0(I.VK_p) = 115.47541;       X0(I.Awarf) = 0;          X0(I.Cwarf) = 0;
X0(I.CA) = 0;         X0(I.AEnox) = 0;              X0(I.AT_III_Heparin) = 0; X0(I.ENO_p) = 0; 
X0(I.AUC) = 0;        X0(I.AVenom) = 0;             X0(I.CVenom) = 0;         X0(I.TaipanVenom) = 0;
X0(I.ATIII) = 0;      X0(I.delayTaipan1) = 0;       X0(I.delayTaipan2) = 0;   X0(I.AVenom_Tiger) = 0; 
X0(I.CVenom_Tiger) = 0; X0(I.AT_III_UFH) = 0;       X0(I.lump) = 1;%X0(I.PT) = 0; 


X0(I.v1) = 50000;           X0(I.k1) = 0.000001; 
X0(I.v2) = 50;              X0(I.k2) = 1;

%%% r3 in IX
X0(I.v3) = 7;               X0(I.k3) = 10;

%%% r4,5 in XI
X0(I.v4) = 7;               X0(I.k4) = 1;
X0(I.v5) = 10;              X0(I.k5) = 10;

%%% r6 in VII
X0(I.v6) = 0.1;             X0(I.k6) = 10;

%%% r7,8,9 in X
X0(I.v7) = 0.02;            X0(I.k7) = 10;
X0(I.v8) = 2;               X0(I.k8) = 0.1;
X0(I.v9) = 0.000000001;     X0(I.k9) = 10;

%%% r10,11 in V
X0(I.v10) = 50000;          X0(I.k10) = 10;
%%% adaptations by Gulati et al for the in-vivo scenario
if contains(model.scenario,'in_vivo')
    X0(I.v10)  =  25000;        X0(I.k10)  =  1800;
end
if contains(model.scenario,'warfarin')||contains(model.scenario,'_ss')
    X0(I.v10) = 50000;          X0(I.k10) = 10;
end

X0(I.v11) = 50;             X0(I.k11) = 1;

%%% r12,13 in II
X0(I.v12) = 100;            X0(I.k12) = 10;
X0(I.v13) = 9;              X0(I.k13) = 500;

%%% r14,15,16,17,18,19 in I 
X0(I.v14) = 20000;          X0(I.k14) = 0.5;
X0(I.v15) = 500;            X0(I.k15) = 500;
%%% adaptations by Gulati et al for the in-vivo scenario
if contains(model.scenario,'in_vivo')
    X0(I.v14) = 21000;          X0(I.k14) = 30000;
end
if contains(model.scenario,'warfarin')||contains(model.scenario,'_ss')
    X0(I.v14) = 20000;          X0(I.k14) = 0.5;
end

X0(I.v16) = 7;              X0(I.k16) = 10;
X0(I.v17) = 7;              X0(I.k17) = 10;
X0(I.v18) = 7;              X0(I.k18) = 100;
X0(I.v19) = 1;              X0(I.k19) = 1;

%%% r20 XIII
X0(I.v20) = 7;              X0(I.k20) = 1;

%%% r21,22,23 in Pg
X0(I.v21) = 7;              X0(I.k21) = 5000;
X0(I.v22) = 5;              X0(I.k22) = 10000;
X0(I.v23) = 2;              X0(I.k23) = 1;

%%% r24 in PC
X0(I.v24) = 7;              X0(I.k24) = 1;

%%% r25 in Xa_Va
X0(I.v25) = 2;              X0(I.k25) = 1;

%%% r26 in VIII and IXa_VIIIa
X0(I.c26) = 0.01;

%%% r27 in X
X0(I.c27) = 0.5;

%%% r28 in II and IIa_Tmod
X0(I.c28) = 0.5;

%%% r29,30 in VII
X0(I.c29) = 0.5;
X0(I.c30) = 0.1;

%%% r31 in TF and VIIa_TF_Xa_TFPI and Xa_TFPI
X0(I.c31) = 0.5;

%%% r32 in X and Xa_TFPI and TFPI
X0(I.c32) = 0.5;

%%% r33 in TF
X0(I.v33) = 70;             X0(I.k33) = 1;

%%% r34 in X
X0(I.v34) = 900;            X0(I.k34) = 200;

%%% r35 in IX
X0(I.v35) = 70;             X0(I.k35) = 1;

%%% r36 in TF
X0(I.v36) = 1000;           X0(I.k36) = 1;

%%% r37 in PC and PS and APC_PS
X0(I.c37) = 0.5;

%%% r38,39,40 in VII
X0(I.v38) = 1;              X0(I.k38) = 10;
X0(I.v39) = 1;              X0(I.k39) = 10;
X0(I.v40) = 0.2;            X0(I.k40) = 10;

%%% r.41,r.42 in Hageman Factor XII
X0(I.v41) = 7;              X0(I.k41) = 1;
X0(I.v42) = 70;             X0(I.k42) = 1;

%%% r43 in Pk
X0(I.v43) = 7;              X0(I.k43) = 1;


%%% -----------------------------------------------------------------------
%%% X0ameters for drug/venom/VK dynamics
%%%

%%% -----------------------------------------------------------------------
%%% vitamin k parameters
%%%
X0(I.VK_k12)  = 0.0587;
X0(I.VK_k21)  = 0.0122;
X0(I.VK_V)    = 24;

%%% -----------------------------------------------------------------------
%%% brown snake venom parameters
%%%
X0(I.ka_Brown)   =  5.0; % absorption rate constant
X0(I.d_Brown)    =  3.5; % degradation rate constant

%%% -----------------------------------------------------------------------
%%% tiger snake venom parameters
%%%
X0(I.ka_Tiger)   =  5.0; % absorption rate constant
X0(I.d_Tiger)    =  3.5; % degradation rate constant

%%% -----------------------------------------------------------------------
%%% taipan snake venom parameters
%%% 
X0(I.d_Taipan)      =  3.5; 
X0(I.ktrans_Taipan) = 0.99;

%%% -----------------------------------------------------------------------
%%% taipan snake venom action on blood coagulation parameters
%%% 
X0(I.vtaipan)      = 70;  X0(I.ktaipan)      =  10; 

%%% -----------------------------------------------------------------------
%%% warfarin PK parameters
%%% 
X0(I.warf_dose) = 4.0; %4mg
pk_pars = "Wajima_augmented";
switch pk_pars
    case "Wajima" % original parameters from Wajima et al.
        X0(I.ka_Warf) = 1.0;
        X0(I.Vd_Warf) = 10;
        X0(I.Cl_Warf) = 0.2;
    case "Hamberg" % Hamberg et al. 2010
        X0(I.ka_Warf) = 2.0;
        X0(I.Vd_Warf) = 14.3;
        Cl1 = 0.174;
        Cl2 = 0.0879;
        Cl3 = 0.0422;
        if isfield(model,'covariates')&& isfield(model.covariates,'cyp')
            X0(I.Cl_Warf) = model.covariates.cyp*[Cl1; Cl2; Cl3];
        else
            %fprintf('\n--> unknown CYP-genotype, expectation used \n\n')
            X0(I.Cl_Warf) = 2*[0.815,0.112,0.073]*[Cl1; Cl2; Cl3];
        end
    case "Wajima_augmented" % Wajima parameters as reference, relative changes for genotypes as in Hamberg
        X0(I.ka_Warf) = 1.0;
        X0(I.Vd_Warf) = 10;
        Cl1 = 0.1124;
        Cl2 = 0.0568;
        Cl3 = 0.0273;
        if isfield(model,'covariates')&& isfield(model.covariates,'cyp')
            X0(I.Cl_Warf) = model.covariates.cyp*[Cl1; Cl2; Cl3];
        else
            %fprintf('\n--> unknown CYP-genotype, expectation used \n\n')
            X0(I.Cl_Warf) = 2*[0.815,0.112,0.073]*[Cl1; Cl2; Cl3];
        end
    otherwise
        fprintf('\n--> unknown pk_pars :-( Please fix!\n\n');
        error('')
end
X0(I.ke_Warf) = X0(I.Cl_Warf)/X0(I.Vd_Warf);

%%% -----------------------------------------------------------------------
%%% warfarin PD parameters
%%% 
X0(I.lmax) = 1;
switch pk_pars
    case "Wajima" % original parameters from Wajima et al.
        X0(I.IC50) = 0.34;
    case "Hamberg" % Hamberg et al. 2010
        IC1 = 2.05*0.068;
        IC2 = 0.96*0.068;
        if isfield(model,'covariates')&& isfield(model.covariates,'cyp')
            X0(I.IC50) = model.covariates.vko*[IC1,  IC2]';
        else
            %fprintf('\n--> unknown CYP-genotype, expectation used \n\n')
            X0(I.IC50) = 2*[0.608,0.392]*[IC1,  IC2]'; % allele frequencies from WARG study (Hamberg 2010), parameters chosen s.t. expected value 0.34
        end
    case "Wajima_augmented" % Wajima parameters as reference, relative changes for genotypes as in Hamberg
        IC1 = 0.2148;
        IC2 = 0.1006;
        if isfield(model,'covariates')&& isfield(model.covariates,'cyp')
            X0(I.IC50) = model.covariates.vko*[IC1,  IC2]';
        else
            %fprintf('\n--> unknown CYP-genotype, expectation used \n\n')
            X0(I.IC50) = 2*[0.608,0.392]*[IC1,  IC2]'; % allele frequencies from WARG study (Hamberg 2010), parameters chosen s.t. expected value 0.34
        end
    otherwise
        fprintf('\n--> unknown pk_pars :-( Please fix!\n\n');
        error('')
end


%%% -----------------------------------------------------------------------
%%% heparin PK parameters (unfractioned heparin UFH or low-molecular weight heparin LMWH)
%%%
X0(I.ke_Hep)  = 0.693;
X0(I.ka_Hep)  = 0.255;
X0(I.Vc_Hep)  = 4.567;
X0(I.Vp_Hep)  = 29.6;
X0(I.Cl_Hep)  = 1.058;
X0(I.Q_Hep)   = 0.62;
X0(I.ke_Hep)  = X0(I.Cl_Hep)/X0(I.Vc_Hep);
X0(I.k12_Hep) = X0(I.Q_Hep)/X0(I.Vc_Hep);
X0(I.k21_Hep) = X0(I.Q_Hep)/X0(I.Vp_Hep);

X0(I.R1) = 1/7.1; %UFH
%X0(I.R1) = 3.9; %enoxaparin

%%% r44 in II
X0(I.c44) = 0.85 * X0(I.R1);

%%% r45 in X
X0(I.c45) = 0.85;
X0(I.R2) = 1; %UFH

%%%X0(I.R2) = 10; %LMWH
% r46 in IX
X0(I.c46) = X0(I.c45) * X0(I.R2);

X0(I.inf_rate_UFH) = 0;%32000/24;

%%% -----------------------------------------------------------------------
%%% degradation rate parameters
%%% All degradation rate constants are in the unit [1/h].
%%%
X0(I.degFg) = 0.032;             X0(I.degF) = 0.050;          X0(I.degXF) = 0.050;             X0(I.degII) = 0.010;
X0(I.degIIa) = 67.4;             X0(I.degTF) = 0.05;          X0(I.degV) = 0.043;              X0(I.degVa) = 20.0;
X0(I.degVII) = 0.12;             X0(I.degVIIa) = 20.0;        X0(I.degVIII) = 0.058;           X0(I.degVIIIa) = 20.0;
X0(I.degIX) = 0.029;             X0(I.degIXa) = 20.0;         X0(I.degX) = 0.018;              X0(I.degXa) = 20.0;
X0(I.degXI) = 0.10;              X0(I.degXIa) = 20.0;         X0(I.degXII) = 0.012;            X0(I.degXIIa) = 20.0;
X0(I.degXIII) = 0.0036;          X0(I.degXIIIa) = 0.69;       X0(I.degPk) = 0.05;              X0(I.degK) = 20.0;
X0(I.degPg) = 0.05;              X0(I.degP) = 20.0;           X0(I.degPC) = 0.050;             X0(I.degAPC) = 20.4;
X0(I.degPS) = 0.0165;            X0(I.degFDP) = 3.5;          X0(I.degD) = 0.1;                X0(I.degTFPI) = 20.0;
X0(I.degVIIaTF) = 20.0;          X0(I.degVIITF) = 0.7;        X0(I.degAPCPS) = 20.0;           X0(I.degXaVa) = 20.0;
X0(I.degIXaVIIIa) = 20.0;        X0(I.degTmod) = 0.050;       X0(I.degIIaTmod) = 20.0;         X0(I.degXaTFPI) = 20.0;
X0(I.degVIIaTFXaTFPI) = 20.0;    X0(I.degTAT) = 0.2;          X0(I.degCA) = 0.05;              X0(I.degVK) = 0.2052;
X0(I.degVK2) = 0.0228;           

%%% degVKH2 and degVKO depend on the initial concentration
%%% VK_V adapted for steady state of VK_p


X0(I.degVKH2) = X0(I.degVK2) * (X0(I.VK)./(X0(I.VKH2)));
X0(I.degVKO)  = X0(I.degVK2) * (X0(I.VK)./(X0(I.VKO)));
%X0(I.VK_V)    = model.X0(I.VK_p)/X0(I.VK)*X0(I.VK_k21)/X0(I.VK_k12);


%%% to avoid numerical issues if either VKH2 or VK are considered as
%%% neglected species degradation rate constants are assumed to be zero
if isfield(I,'neg') && (ismember(I.VKH2,I.neg) || ismember(I.VK,I.neg))
    X0(I.degVKH2) = 0;
    X0(I.degVKO)  = 0;
end

%%% -----------------------------------------------------------------------
%%% production rate parameters
%%%

X0(I.aII)      = X0(I.degII)  * X0(I.II)  / (X0(I.VKH2));
X0(I.aVII)     = X0(I.degVII) * X0(I.VII) / (X0(I.VKH2)); 
X0(I.aIX)      = X0(I.degIX)  * X0(I.IX)  / (X0(I.VKH2));
X0(I.aX)       = X0(I.degX)   * X0(I.X)   / (X0(I.VKH2));
X0(I.aPC)      = X0(I.degPC)  * X0(I.PC)  / (X0(I.VKH2)); 
X0(I.aPS)      = X0(I.degPS)  * X0(I.PS)  / (X0(I.VKH2)); 

X0(I.pV)       = X0(I.degV)    * X0(I.V); 
X0(I.pVIII)    = X0(I.degVIII) * X0(I.VIII);
X0(I.pFg)      = X0(I.degFg)   * X0(I.Fg);  
X0(I.pXIII)    = X0(I.degXIII) * X0(I.XIII);
X0(I.pPg)      = X0(I.degPg)   * X0(I.Pg);   
X0(I.pTmod)    = X0(I.degTmod) * X0(I.Tmod);
X0(I.pXI)      = X0(I.degXI)   * X0(I.XI);  
X0(I.pTFPI)    = X0(I.degTFPI) * X0(I.TFPI);
X0(I.pVK)      = X0(I.degVK)   * X0(I.VK);  
X0(I.pXII)     = X0(I.degXII)  * X0(I.XII);
X0(I.pPk)      = X0(I.degPk)   * X0(I.Pk);

%%% in the case of the in vitro setting (PT or aPTT) all production rates are assumed to
%%% be zero (see Wajima et al)
if contains(model.scenario,'in_vitro')
        X0(I.aII) = 0;    X0(I.aVII)  = 0;
        X0(I.aIX) = 0;    X0(I.aX)    = 0;
        X0(I.aPC) = 0;    X0(I.aPS)   = 0;
        X0(I.pV)  = 0;    X0(I.pVIII) = 0;
        X0(I.pFg) = 0;    X0(I.pXIII) = 0;
        X0(I.pPg) = 0;    X0(I.pTmod) = 0;
        X0(I.pXI) = 0;    X0(I.pTFPI) = 0;
        X0(I.pVK) = 0;    X0(I.pXII)  = 0;
        X0(I.pPk) = 0;
end

end
