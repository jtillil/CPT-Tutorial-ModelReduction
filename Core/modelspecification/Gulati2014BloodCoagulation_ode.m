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

function dX = Gulati2014BloodCoagulation_ode(t,X,par,model)

%%% assign model indexing
I  = model.I;

%%% -----------------------------------------------------------------------
%%% account for state variable (in I.con) that are eliminated 
%%% via a conservation law 
%%%

% determine value of states eliminated via conservation laws
% for p = 1:length(I.con)
% 
%     L = model.L;
%     % determine value of the state in I.con by subtracting the sum of the
%     % current values of the remaining (rem) states of the conlaw. Note:
%     % States in I.con are not in the same order as the conservation laws
%     % have been specifed (but in some permuted (p) order).
% 
%     k = I.con(p);  % index of state variable
%     c = L.con2conlaw(p); % number of conservation law
% 
%     % determine value of state using the corresponding conlaw. If resulting
%     % value is ''too negative'' (it should minimally be zero) or 'too large'
%     % (it should maximally be 'valconlaw'), then skip solving
%     % the ODE and report '(error)' for the kth state. Otherwise, ensure
%     % that value is non-negative
%     X(k) = L.valconlaw(c) - sum( X(L.remstatesofconlaw{p}) );
%     conlowTOL = min( L.valconlaw(c)*0.01, 1e-3 ); % rel 1% or abs 1e-3
%     if ( -conlowTOL < X(k) ) && ( X(k) < L.valconlaw(c) + conlowTOL )
%         X(k) = max(0,X(k));
%     else
%         return;
%     end
% 
% end

%%% -----------------------------------------------------------------------
%%% specify system of ODEs 
dX = 0*X;

%%% NOTE: The action of brown snake venom was assumed to be identical to the human
%%% prothrombinase complex (Xa:Va) based on their structural similarity
%%% states: XII, XIIa
%%% eqs. 1 and 2: Hageman- Factor (XII) and activation
%%%
r41 = ( par(I.v41)*X(I.CA) ) / ( par(I.k41)+X(I.CA) ) * X(I.XII);
r42 = ( par(I.v42)*X(I.K) ) / ( par(I.k42)+X(I.K) ) * X(I.XII);

dX(I.XII)  = par(I.pXII)-r41-r42-par(I.degXII) * X(I.XII);
dX(I.XIIa) = r41+r42-par(I.degXIIa) * X(I.XIIa);

%%% -----------------------------------------------------------------------
%%% states: VIII, VIIIa
%%% eqs. 3 and 4 : Antihaemophiles Globulin A (VIII) and activation
%%%
r1  = ( par(I.v1)*X(I.IIa) ) / ( par(I.k1)+X(I.IIa) ) * X(I.VIII);
r2  = ( par(I.v2)*X(I.APC_PS) ) / ( par(I.k2)+X(I.APC_PS) ) * X(I.VIIIa);
r26 = ( X(I.IXa)*X(I.VIIIa) ) / par(I.c26);

dX(I.VIII)  = par(I.pVIII)-r1-par(I.degVIII)*X(I.VIII);
dX(I.VIIIa) = r1-r2-r26-par(I.degVIIIa)*X(I.VIIIa);

%%% -----------------------------------------------------------------------
%%% states: IX, IXa
%%% eqs. 5 and 6 : Antihaemophiles Globulin B (IX) and activation
%%%
r3  = (par(I.v3)*X(I.XIa))/(par(I.k3)+X(I.XIa)) * X(I.IX);
r35 = (par(I.v35)*X(I.VIIa_TF))/(par(I.k35)+X(I.VIIa_TF)) * X(I.IX);
r46 = (X(I.IXa)*(X(I.AT_III_Heparin)))/par(I.c46);
pIX = par(I.aIX) * X(I.VKH2);

dX(I.IX) = pIX-r35-r3-par(I.degIX)*X(I.IX);
dX(I.IXa) = r35+r3-r26-r46-par(I.degIXa)*X(I.IXa);

%%% -----------------------------------------------------------------------
%%% states: XI, XIa
%%% eqs. 7 and 8 : Rosenthal-factor (XI) and activation
%%%
r4 = (par(I.v4)*X(I.XIIa))/(par(I.k4)+X(I.XIIa)) * X(I.XI);
r5 = (par(I.v5)*X(I.IIa))/(par(I.k5)+X(I.IIa)) * X(I.XI);

dX(I.XI)  = par(I.pXI)-r5-r4-par(I.degXI)*X(I.XI);
dX(I.XIa) = r5+r4-par(I.degXIa)*X(I.XIa);

%%% -----------------------------------------------------------------------
%%% states: VII, VIIa
%%% eqs. 9 and 10: Proconvertin (VII) and activation
%%%

% consensus -->  r6   = ( par(I.v6)*X(Ix.IIa) ) / ( par(I.k6)+X(Ix.IIa) ) * X(Ix.VII);

r6   = (par(I.v6)*X(I.IIa))/(par(I.k6)+X(I.IIa)) * X(I.VII);
r29  = (X(I.TF)*X(I.VIIa))/(par(I.c29));
r30  = (X(I.TF)*X(I.VII))/(par(I.c30));
r38  = (par(I.v38)*X(I.Xa))/(par(I.k38)+X(I.Xa)) * X(I.VII);
r39  = (par(I.v39)*X(I.VIIa_TF))/(par(I.k39)+X(I.VIIa_TF)) * X(I.VII);
r40  = (par(I.v40)*X(I.IXa))/(par(I.k40)+X(I.IXa)) * X(I.VII);
rtaipan  = (par(I.vtaipan)*X(I.delayTaipan2))/(par(I.ktaipan)+X(I.delayTaipan2)) * X(I.VII);
pVII = par(I.aVII)*X(I.VKH2);

dX(I.VII)  = pVII-r30-r6-par(I.degVII)*X(I.VII)-r38-r39-r40-rtaipan;
dX(I.VIIa) = -r29+r6-par(I.degVIIa)*X(I.VIIa)+r38+r39+r40+rtaipan;

%%% -----------------------------------------------------------------------
%%% eqs. 11 and 12: Stuart-Power-Factor (X) and activation
%%%
r7  = (par(I.v7)*X(I.IXa))/(par(I.k7)+X(I.IXa)) * X(I.X);
r8  = (par(I.v8)*X(I.IXa_VIIIa))/(par(I.k8)+X(I.IXa_VIIIa)) * X(I.X);
r9  = (par(I.v9)*X(I.VIIa))/(par(I.k9)+X(I.VIIa)) * X(I.X);
r27 = (X(I.Xa)*X(I.Va))/(par(I.c27));
r32 = (X(I.Xa)*X(I.TFPI))/(par(I.c32));
r34 = (par(I.v34)*X(I.VIIa_TF))/(par(I.k34)+X(I.VIIa_TF)) * X(I.X);
r45 = (X(I.Xa)*(X(I.AT_III_Heparin)))/par(I.c45);
pX  = par(I.aX)*X(I.VKH2);

dX(I.X)  = pX-r9-r34-r7-r8-par(I.degX)*X(I.X);
dX(I.Xa) = r9+r34+r7+r8-r27-r32-r45-par(I.degXa)*X(I.Xa);

%%% -----------------------------------------------------------------------
%%% eqs. 13 and 14 : Proaccelerin (V) and activation
%%%
r10 = (par(I.v10)*X(I.IIa))/(par(I.k10)+X(I.IIa)) * X(I.V);
r11 = (par(I.v11)*X(I.APC_PS))/(par(I.k11)+X(I.APC_PS)) * X(I.Va);
r27 = (X(I.Xa)*X(I.Va))/(par(I.c27));

dX(I.V)  = par(I.pV)-r10-par(I.degV)*X(I.V);
dX(I.Va) = r10-r11-r27-par(I.degVa)*X(I.Va);

%%% -----------------------------------------------------------------------
%%% eq. 15 : Complex Xa_Va
%%%
r25 = (par(I.v25)*X(I.APC_PS))/(par(I.k25)+X(I.APC_PS)) * X(I.Xa_Va);
r27 = (X(I.Xa)*X(I.Va))/(par(I.c27));

dX(I.Xa_Va) = r27-r25-par(I.degXaVa)*X(I.Xa_Va);

%%% -----------------------------------------------------------------------
%%% eqs. 16 and 17 : Prothrombin (II) and activation
%%%
r12 =  (par(I.v12)*(X(I.Xa_Va)+X(I.CVenom)+X(I.TaipanVenom))) / ...
          (par(I.k12)+(X(I.Xa_Va)+X(I.CVenom)+X(I.TaipanVenom))) * X(I.II);
r13 = (par(I.v13)*(X(I.Xa)+X(I.CVenom_Tiger)))/...
          (par(I.k13)+(X(I.Xa)+X(I.CVenom_Tiger))) * X(I.II);
r28 = (X(I.IIa)*X(I.Tmod))/(par(I.c28));
r44 = (X(I.IIa)*(X(I.AT_III_Heparin)))/par(I.c44);
pII = par(I.aII)*X(I.VKH2);

dX(I.II)  = pII-r12-r13-par(I.degII)*X(I.II);
dX(I.IIa) = r12+r13-r28-r44-par(I.degIIa)*X(I.IIa);

%%% -----------------------------------------------------------------------
%%% eq. 18: thrombin-antithrombin (TAT)
%%%

dX(I.TAT) = par(I.degIIa)*X(I.IIa)-par(I.degTAT)*X(I.TAT);

%%% -----------------------------------------------------------------------
%%% eqs. 19, 20, 21, 22 and 23 : Fibrinogen (I), activation and complex/linked
%%%
r14 = (par(I.v14)*X(I.IIa))/(par(I.k14)+X(I.IIa)) * X(I.Fg);
r15 = (par(I.v15)*X(I.P))/(par(I.k15)+X(I.P)) * X(I.Fg);
r16 = (par(I.v16)*X(I.XIIIa))/(par(I.k16)+X(I.XIIIa)) * X(I.F);
r17 = (par(I.v17)*X(I.P))/(par(I.k17)+X(I.P)) * X(I.F);
r18 = (par(I.v18)*X(I.P))/(par(I.k18)+X(I.P)) * X(I.XF);
r19 = (par(I.v19)*X(I.APC_PS))/(par(I.k19)+X(I.APC_PS)) * X(I.XF);

dX(I.Fg)  = par(I.pFg)-r14-r15-par(I.degFg)*X(I.Fg);
dX(I.F)   = r14-r16-r17-par(I.degF)*X(I.F);
dX(I.XF)  = r16-r18-r19-par(I.degXF)*X(I.XF);
dX(I.FDP) = r15+par(I.degFg)*X(I.Fg)+r17+par(I.degF)*X(I.F)-par(I.degFDP)*X(I.FDP);
dX(I.D)   = r18+r19+par(I.degXF)*X(I.XF)-par(I.degD)*X(I.D);


%%% -----------------------------------------------------------------------
%%% eqs. 24 and 25 : Fibrinstabilisierender Factor (XIII) and activation
%%%
r20 = (par(I.v20)*X(I.IIa))/(par(I.k20)+X(I.IIa)) * X(I.XIII);

dX(I.XIII)  = par(I.pXIII)-r20-par(I.degXIII)*X(I.XIII);
dX(I.XIIIa) = r20-par(I.degXIIIa)*X(I.XIIIa);

%%% -----------------------------------------------------------------------
%%% eqs. 26 and 27: Plasminogen (Pg) and Plasmin (P)
%%%
r21 = (par(I.v21)*X(I.IIa))/(par(I.k21)+X(I.IIa)) * X(I.Pg);
r22 = (par(I.v22)*X(I.F))/(par(I.k22)+X(I.F)) * X(I.Pg);
r23 = (par(I.v23)*X(I.APC_PS))/(par(I.k23)+X(I.APC_PS)) * X(I.Pg);

dX(I.Pg) = par(I.pPg)-r21-r23-r22-par(I.degPg)*X(I.Pg);
dX(I.P)  = r21+r23+r22-par(I.degP)*X(I.P);

%%% -----------------------------------------------------------------------
%%% eqs. 28 and 29 : Protein C (PC) and activation (APC)
%%%
r24 = (par(I.v24)*X(I.IIa_Tmod))/(par(I.k24)+X(I.IIa_Tmod)) * X(I.PC);
r37 = (X(I.APC)*X(I.PS))/(par(I.c37));
pPC = par(I.aPC)*X(I.VKH2);

dX(I.PC)  = pPC-r24-par(I.degPC)*X(I.PC);
dX(I.APC) = r24-r37-par(I.degAPC)*X(I.APC);

%%% -----------------------------------------------------------------------
%%% eq. 30 : Thrombmodulin (Tmod)
%%%
r28 = (X(I.IIa)*X(I.Tmod))/(par(I.c28));

dX(I.Tmod) = par(I.pTmod)-r28-par(I.degTmod)*X(I.Tmod);

%%% -----------------------------------------------------------------------
%%% eq. 31 : IIa_Tmod complex
%%%

dX(I.IIa_Tmod) = r28-par(I.degIIaTmod)*X(I.IIa_Tmod);

%%% -----------------------------------------------------------------------
%%% eq. 32 : IXa_VIIIa complex
%%%
r26 = (X(I.IXa)*X(I.VIIIa))/par(I.c26);

dX(I.IXa_VIIIa) = r26-par(I.degIXaVIIIa)*X(I.IXa_VIIIa);

%%% -----------------------------------------------------------------------
%%% eqs. 33, 34 and 35 : tissue factor also called III (TF) and 
%%% its complexes (VII_TF and VIIa_TF)
%%%
r33 = (par(I.v33)*X(I.Xa))/(par(I.k33)+X(I.Xa)) * X(I.VII_TF);
r36 = (par(I.v36)*X(I.TF))/(par(I.k36)+X(I.TF)) * X(I.VII_TF);
r31 = (X(I.VIIa_TF)*X(I.Xa_TFPI))/(par(I.c31));
r29 = (X(I.TF)*X(I.VIIa))/(par(I.c29));
r30 = (X(I.TF)*X(I.VII))/(par(I.c30));

dX(I.TF) = -r29-r30-par(I.degTF)*X(I.TF);

dX(I.VII_TF)  = r30-r36-r33-par(I.degVIITF)*X(I.VII_TF);
dX(I.VIIa_TF) = r29+r36+r33-r31-par(I.degVIIaTF)*X(I.VIIa_TF);

%%% -----------------------------------------------------------------------
%%% eqs. 36, 37 and 38 :tissue factor pathway inhibitor (TFPI) and its
%%% complexes (Xa_TFPI and VIIa_TF_Xa_TFPI)
%%%
r31 = (X(I.VIIa_TF)*X(I.Xa_TFPI))/(par(I.c31));
r32 = (X(I.Xa)*X(I.TFPI))/(par(I.c32));

dX(I.TFPI) = par(I.pTFPI)-r32-par(I.degTFPI)*X(I.TFPI);
dX(I.Xa_TFPI) = r32-r31-par(I.degXaTFPI)*X(I.Xa_TFPI);

r31 = (X(I.VIIa_TF)*X(I.Xa_TFPI))/(par(I.c31));

dX(I.VIIa_TF_Xa_TFPI) = r31-par(I.degVIIaTFXaTFPI)*X(I.VIIa_TF_Xa_TFPI);

%%% -----------------------------------------------------------------------
%%% eqs. 39 and 40 Protein S (PS) and complex (APC_PS)
%%%
pPS = par(I.aPS)*X(I.VKH2);
r37 = (X(I.APC)*X(I.PS))/(par(I.c37));

dX(I.PS)     = pPS-r37-par(I.degPS)*X(I.PS);
dX(I.APC_PS) = r37-par(I.degAPCPS)*X(I.APC_PS);

%%% -----------------------------------------------------------------------
%%% eqs. 41 and 42 :Prekallikrein (PK) and Kallikrein (K)
%%%
r43 = (par(I.v43)*X(I.XIIa))/(par(I.k43)+X(I.XIIa))*X(I.Pk);

dX(I.Pk) = par(I.pPk)-r43-par(I.degPk)*X(I.Pk);
dX(I.K)  = r43-par(I.degK)*X(I.K);

%%% -----------------------------------------------------------------------
%%% eqs. 43,44,45 and 46 : Vitamin K (VK)
%%%
r47 = 1-(par(I.lmax)*X(I.Cwarf))/(par(I.IC50)+X(I.Cwarf));
r48 = 1-(par(I.lmax)*X(I.Cwarf))/(par(I.IC50)+X(I.Cwarf));

dX(I.VK)   = par(I.pVK)-par(I.degVK2)*X(I.VK)*r47-par(I.degVK)*X(I.VK) + ...
             par(I.degVKO)*X(I.VKO)*r48-par(I.VK_k12)*X(I.VK)+(par(I.VK_k21)*X(I.VK_p)/(par(I.VK_V)));
dX(I.VKH2) = par(I.degVK2)*X(I.VK)*r47-par(I.degVKH2)*X(I.VKH2);
dX(I.VKO)  = -par(I.degVKO)*X(I.VKO)*r48+par(I.degVKH2)*X(I.VKH2);
dX(I.VK_p) = par(I.VK_k12)*X(I.VK)*par(I.VK_V)-par(I.VK_k21)*X(I.VK_p);

%%% -----------------------------------------------------------------------
%%% eqs. 47,48 : Warfarin PK (oral)
%%%

dX(I.Awarf)=-par(I.ka_Warf)*X(I.Awarf);
dX(I.Cwarf)=par(I.ka_Warf)*(X(I.Awarf)/par(I.Vd_Warf))-par(I.ke_Warf)*X(I.Cwarf);

%%% -----------------------------------------------------------------------
%%% eq. 49 :activator for the contact system (CA)
%%%
dX(I.CA) = -par(I.degCA)*X(I.CA);

%%% -----------------------------------------------------------------------
%%% eqs. 50,51 and 52 : Heparin (LMWH) PK 
%%%

dX(I.AEnox)             = -par(I.ka_Hep)*X(I.AEnox);
dX(I.AT_III_Heparin)    = par(I.ka_Hep)*(X(I.AEnox)/par(I.Vc_Hep))-par(I.k12_Hep)*X(I.AT_III_Heparin)...
                            +par(I.k21_Hep)*X(I.ENO_p)*par(I.Vc_Hep)-r44-r45-r46...
                            -par(I.ke_Hep)*X(I.AT_III_Heparin);
dX(I.ENO_p)             = par(I.k12_Hep)*X(I.AT_III_Heparin)*par(I.Vc_Hep)-par(I.k21_Hep)*X(I.ENO_p); 

%ENO_p: amount of enoxaparin in peripheral compartment
   

%%% -----------------------------------------------------------------------
%%% eqs. 53 : AUC of fibrin concentration 
%%%

dX(I.AUC) = X(I.F);

%%% -----------------------------------------------------------------------
%%% eqs. 54 and 55 : brown snake venom with absorption compartment (A) and
%%% concentration compartment (C)
%%%

dX(I.AVenom)  =  -par(I.ka_Brown)*X(I.AVenom);
dX(I.CVenom)  =  par(I.ka_Brown)*X(I.AVenom) -par(I.d_Brown)*X(I.CVenom);

%%% -----------------------------------------------------------------------
%%% eqs. 56,57 and 58  : Taipan snake venom & delay compartment for action of Taipan
%%% on VII
%%%

dX(I.TaipanVenom)  =  -par(I.d_Taipan)*X(I.TaipanVenom);
dX(I.delayTaipan1)  =  par(I.d_Taipan)*X(I.TaipanVenom)- par(I.ktrans_Taipan)*X(I.delayTaipan1);
dX(I.delayTaipan2)  =  par(I.ktrans_Taipan)*X(I.delayTaipan1)- par(I.ktrans_Taipan)*X(I.delayTaipan2);

%%% -----------------------------------------------------------------------
%%% eqs. 59  : ATIII
%%%

dX(I.ATIII)  =  0;
%%% -----------------------------------------------------------------------
%%% eqs. 60 and 61 : tiger snake venom with absorption compartment (A) and
%%% concentration compartment (C)
%%%

dX(I.AVenom_Tiger)  =  -par(I.ka_Tiger)*X(I.AVenom_Tiger);
dX(I.CVenom_Tiger)  =  par(I.ka_Tiger)*X(I.AVenom_Tiger) -par(I.d_Tiger)*X(I.CVenom_Tiger);

%%% -----------------------------------------------------------------------
%%% eqs. 62 : Herparin (UFH) PK
%%%
dX(I.AT_III_UFH) = par(I.ke_Hep)*X(I.AT_III_UFH)-r44-r45-r46+par(I.inf_rate_UFH);

%%% -----------------------------------------------------------------------
% dX(I.lump) = 0;


%%% set derivative of environmental, negligible and conservation law 
%%% state variables to zero
%%%
%%%------------------------------------------------

dX(I.env) = 0; dX(I.neg) = 0; 
dX=dX(:);

end