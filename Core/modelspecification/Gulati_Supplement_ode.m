
% ORDINARY DIFFERENTIAL EQUATIONS FOR THE ORIGINAL 62-STATE COAGULATION SYSTEMS MODEL

% Initial conditions are provided in Supplementary Table 1. Parameter values can be obtained from the following publications:
%
% Wajima, T., Isbister, G.K. & Duffull, S.B. A comprehensive model for the humoral coagulation network in humans. Clin Pharmacol Ther 86, 290-298 (2009)
% Gulati, A., Faed, J.M., Isbister, G.K. & Duffull, S.B. Development and evaluation of a prototype of a novel clotting time test to monitor enoxaparin. Pharm Res 29, 225-235 (2012)
% Gulati, A., Isbister, G.K. & Duffull, S.B. Effect of Australian elapid venoms on blood coagulation: Australian Snakebite Project (ASP-17). Toxicon 61, 94-104 (2013)

function dxdt=Gulati_Supplement_ode(t,x,pX,pV,pII,pVIII,pIX,pFg,pXIII,pPg,pPC,pTmod,pXI,pVII,pTFPI,pPS,pVK,pXII,pPK,dX,dXa,dV,dVa,dVaXa,dII,dIIa,dVIII,dVIIIa,dIXa,dIXaVIIIa,dIX,dXIa,dFg,dFDP,dF,dXF,dD,dXIII,dXIIIa,dPg,dP,dPC,dAPC,dTmod,dIIaTmod,dTaipan,dXI,dXIIa,dVII,dVIIa,dTAT,dBrown,dTF,dVIITF,dVIIaTF,dTFPI,dXaTFPI,dVIIaTFXaTFPI,dPS,dAPCPS,dHeparin,dVK_1,dVK_2,dVK_3,dVKH2,dVKO,dXII,dPK,dK,dca,vIXa,kIXa,vII2IIaXa,kII2IIaXa,vV2VaIIa,kV2VaIIa,vVIII2VIIIaIIa,kVIII2VIIIaIIa,vIXaVIIIa,kIXaVIIIa,vXIa,kXIa,vVaXa,kVaXa,vFg2F,kFg2F,vFg2FDP,kFg2FDP,vF2FDP,kF2FDP,vF2XF,kF2XF,vXF2D,kXF2D,vXIII2XIIIa,kXIII2XIIIa,vPg2PF,kPg2PF,vPg2PIIa,kPg2PIIa,vPC2APC,kPC2APC,vVIIIaLoss,kVIIIaLoss,vVaLoss,kVaLoss,vXaVaLoss,kXaVaLoss,vXF2DAPC,kXF2DAPC,vPg2PAPC,kPg2PAPC,vXI2XIaIIa,kXI2XIaIIa, vXI2XIaXIIa,kXI2XIaXIIa,vX2XaVIIa,kX2XaVIIa,vVII2VIIaIIa,kVII2VIIaIIa,vVII2VIIaTaipan,kVII2VIIaTaipan,vVIITF_Xa,kVIITF_Xa,vX_VIIaTF,kX_VIIaTF,vIX_VIIaTF,kIX_VIIaTF,vVIITF_TF,kVIITF_TF,vVII_Xa,kVII_Xa,vVII_VIIaTF,kVII_VIIaTF,vVII_IXa,kVII_IXa,vXII_ca,kXII_ca,vXII_K,kXII_K,vPK_XIIa,kPK_XIIa,cVaXa,cIXaVIIIa,cIIaTmod,cVIIaTF,cVIITF,cVIIaTFXaTFPI,cXaTFPI,cAPCPS,cXaHeparin,cIIaHeparin,cIXaHeparin,cIXaUFH,cXaUFH,cIIaUFH,Imax,IC50,ka,v,ke,lagt,ka_e,k10_e,k12_e,k21_e,v2_e,k12_vk,k21_vk,vc_vk,ke_UFH,R_UFH,T_UFH,time_of_antivenom,kr,FLAG1,FLAG2,KDXa,KDE,ka_brown,ka_tiger,dTiger)

effective_XaVa = x(5) + x(62);
effective_Xa = x(2) + x(61);

Enox_conc = 1/100/4500*10^9;
Bmax = x(34);
f_KD=1-((Enox_conc/v2_e)/(KDE+(Enox_conc/v2_e)));
binding=Bmax*x(2)/(f_KD*KDXa+x(2));

% for anitivenom treatment
% if t>time_of_antivenom
%     antivenom=100;
% else
%     antivenom=1;
% end

% for heparin (UFH) infusion
% if t>T_UFH
%     INF_UFH=0;
% else
%     INF_UFH=1;
% end

% FOR SYMBOLIC ODE CALCULATION
antivenom = 1;
INF_UFH = 0;

    %(1) X
dxdt=[pX*x(52)-dX*x(1)-vIXa*x(10)*x(1)/(kIXa+x(10))-vIXaVIIIa*x(11)*x(1)/(kIXaVIIIa+x(11))-vX2XaVIIa*x(32)*x(1)/(kX2XaVIIa+x(32))-vX_VIIaTF*x(39)*x(1)/(kX_VIIaTF+x(39))
    %(2) Xa
    vIXa*x(10)*x(1)/(kIXa+x(10))+vIXaVIIIa*x(11)*x(1)/(kIXaVIIIa+x(11))-dXa*x(2)-effective_Xa*x(4)/cVaXa+vX2XaVIIa*x(32)*x(1)/(kX2XaVIIa+x(32))+vX_VIIaTF*x(39)*x(1)/(kX_VIIaTF+x(39))-x(2)*x(40)/cXaTFPI-x(2)*x(49)*(1-FLAG1)/cXaHeparin-x(2)*x(59)/cXaUFH-binding*FLAG2
    %(3) V
    pV-dV*x(3)-vV2VaIIa*x(3)*x(7)/(kV2VaIIa+x(7))
    %(4) Va
    vV2VaIIa*x(3)*x(7)/(kV2VaIIa+x(7))-dVa*x(4)-effective_Xa*x(4)/cVaXa-vVaLoss*x(46)*x(4)/(kVaLoss+x(46))
    %(5) VaXa
    effective_Xa*x(4)/cVaXa-dVaXa*x(5)-vXaVaLoss*x(46)*x(5)/(kXaVaLoss+x(46))
    %(6) II
    pII*x(52)-dII*x(6)-vVaXa*effective_XaVa*x(6)/(kVaXa+effective_XaVa)-vII2IIaXa*x(6)*effective_Xa/(kII2IIaXa+effective_Xa)-vVaXa*x(21)*x(6)/(kVaXa+x(21))
    %(7) IIa
    vII2IIaXa*x(6)*effective_Xa/(kII2IIaXa+effective_Xa)+vVaXa*effective_XaVa*x(6)/(kVaXa+effective_XaVa)-dIIa*x(7)-x(7)*x(26)/cIIaTmod+vVaXa*x(21)*x(6)/(kVaXa+x(21))-x(7)*x(49)*(1-FLAG1)/cIIaHeparin-x(7)*x(59)/cIIaUFH-x(7)*Enox_conc*FLAG1/cIIaHeparin
    %(8) VIII
    pVIII-dVIII*x(8)-vVIII2VIIIaIIa*x(7)*x(8)/(kVIII2VIIIaIIa+x(7))
    %(9) VIIIa
    vVIII2VIIIaIIa*x(7)*x(8)/(kVIII2VIIIaIIa+x(7))-dVIIIa*x(9)-x(9)*x(10)/cIXaVIIIa-vVIIIaLoss*x(46)*x(9)/(kVIIIaLoss+x(46))
    %(10) IXa
    vXIa*x(12)*x(13)/(kXIa+x(13))-dIXa*x(10)-x(9)*x(10)/cIXaVIIIa+vIX_VIIaTF*x(39)*x(12)/(kIX_VIIaTF+x(39))-x(10)*x(49)*(1-FLAG1)/cIXaHeparin-x(10)*x(59)/cIXaUFH-x(10)*Enox_conc*FLAG1/cIXaHeparin
    %(11) IXaVIIIa
    x(9)*x(10)/cIXaVIIIa-dIXaVIIIa*x(11)
    %(12) IX
    pIX*x(52)-dIX*x(12)-vXIa*x(12)*x(13)/(kXIa+x(13))-vIX_VIIaTF*x(39)*x(12)/(kIX_VIIaTF+x(39))
    %(13) XIa
    vXI2XIaIIa*x(7)*x(29)/(kXI2XIaIIa+x(7))+vXI2XIaXIIa*x(30)*x(29)/(kXI2XIaXIIa+x(30))-dXIa*x(13)
    %(14) Fibrinogen
    pFg-dFg*x(14)-vFg2F*x(7)*x(14)/(kFg2F+x(7))-vFg2FDP*x(23)*x(14)/(kFg2FDP+x(23))
    %(15) FDPs
    dFg*x(14)+vF2FDP*x(23)*x(16)/(kF2FDP+x(23))+vFg2FDP*x(23)*x(14)/(kFg2FDP+x(23))+dF*x(16)-dFDP*x(15)
    %(16) Fibrin
    vFg2F*x(7)*x(14)/(kFg2F+x(7))-dF*x(16)-vF2FDP*x(23)*x(16)/(kF2FDP+x(23))-vF2XF*x(20)*x(16)/(kF2XF+x(20))
    %(17) X-linked Fibrin
    vF2XF*x(20)*x(16)/(kF2XF+x(20))-dXF*x(17)-vXF2D*x(23)*x(17)/(kXF2D+x(23))-vXF2DAPC*x(46)*x(17)/(kXF2DAPC+x(46))
    %(18) D-dimers
    dXF*x(17)+vXF2D*x(23)*x(17)/(kXF2D+x(23))+vXF2DAPC*x(46)*x(17)/(kXF2DAPC+x(46))-dD*x(18)
    %(19) XIII
    pXIII-dXIII*x(19)-vXIII2XIIIa*x(7)*x(19)/(kXIII2XIIIa+x(7))
    %(20) XIIIa
    vXIII2XIIIa*x(7)*x(19)/(kXIII2XIIIa+x(7))-dXIIIa*x(20)
    %(21) Taipan venom
    -dTaipan*x(21)*antivenom
    %(22) Plasminogen
    pPg-dPg*x(22)-vPg2PF*x(16)*x(22)/(kPg2PF+x(16))-vPg2PIIa*x(7)*x(22)/(kPg2PIIa+x(7))- vPg2PAPC*x(46)*x(22)/(kPg2PAPC +x(46))
    %(23) Plasmin
    vPg2PF*x(16)*x(22)/(kPg2PF+x(16))+vPg2PIIa*x(7)*x(22)/(kPg2PIIa+x(7))+vPg2PAPC*x(46)*x(22)/(kPg2PAPC +x(46))-dP*x(23)
    %(24) Protein C
    pPC*x(52)-dPC*x(24)-vPC2APC*x(27)*x(24)/(kPC2APC+x(27))
    %(25) Activated Protein C
    vPC2APC*x(27)*x(24)/(kPC2APC+x(27))-x(25)*x(45)/cAPCPS-dAPC*x(25)
    %(26) Thrombmodulin
    pTmod-dTmod*x(26)-x(7)*x(26)/cIIaTmod
    %(27) IIa:Tmod
    x(7)*x(26)/cIIaTmod-dIIaTmod*x(27)
    %(28) Brown snake venom absorption site
    -ka_brown*x(28)
    %(29) XI
    pXI-dXI*x(29)-vXI2XIaIIa*x(7)*x(29)/(kXI2XIaIIa+x(7))-vXI2XIaXIIa*x(30)*x(29)/(kXI2XIaXIIa+x(30))
    %(30) XIIa
    -dXIIa*x(30)+vXII_ca*x(58)*x(55)/(kXII_ca+x(58))+vXII_K*x(57)*x(55)/(kXII_K+x(57))
    %(31) VII
    pVII*x(52)-dVII*x(31)-vVII2VIIaIIa*x(7)*x(31)/(kVII2VIIaIIa+x(7))-vVII2VIIaTaipan*x(36)*x(31)/(kVII2VIIaTaipan+x(36))-x(31)*x(37)/cVIITF-vVII_Xa*x(2)*x(31)/(kVII_Xa+x(2))-vVII_VIIaTF*x(39)*x(31)/(kVII_VIIaTF+x(39))-vVII_IXa*x(10)*x(31)/(kVII_IXa+x(10))
    %(32) VIIa
    vVII2VIIaIIa*x(7)*x(31)/(kVII2VIIaIIa+x(7))+vVII2VIIaTaipan*x(36)*x(31)/(kVII2VIIaTaipan+x(36))-dVIIa*x(32)-x(37)*x(32)/cVIIaTF+vVII_Xa*x(2)*x(31)/(kVII_Xa+x(2))+vVII_VIIaTF*x(39)*x(31)/(kVII_VIIaTF+x(39))+vVII_IXa*x(10)*x(31)/(kVII_IXa+x(10))
    %(33) TAT
    dIIa*x(7)-dTAT*x(33)
    %(34) ATIII
    0
    %(35) delay for taipan on VII
    dTaipan*x(21)*antivenom-kr*x(35)
    %(36) another delay compartment for taipan on VII
    kr*x(35)-kr*x(36)
    %(37) TF
    -x(37)*x(32)/cVIIaTF-x(31)*x(37)/cVIITF-dTF*x(37)
    %(38) VII-TF
    x(37)*x(31)/cVIITF-vVIITF_TF*x(37)*x(38)/(kVIITF_TF+x(37))-vVIITF_Xa*x(2)*x(38)/(kVIITF_Xa+x(2))-dVIITF*x(38)
    %(39) VIIa-TF
    x(37)*x(32)/cVIIaTF+vVIITF_TF*x(37)*x(38)/(kVIITF_TF+x(37))+vVIITF_Xa*x(2)*x(38)/(kVIITF_Xa+x(2))-x(39)*x(41)/cVIIaTFXaTFPI-dVIIaTF*x(39)
    %(40) TFPI
    pTFPI-x(2)*x(40)/cXaTFPI-dTFPI*x(40)
    %(41) Xa-TFPI
    x(2)*x(40)/cXaTFPI-x(39)*x(41)/cVIIaTFXaTFPI-dXaTFPI*x(41)
    %(42) VIIa-TF-Xa-TFPI
    x(39)*x(41)/cVIIaTFXaTFPI-dVIIaTFXaTFPI*x(42)
    %(43) warfarin absorption site
    -ka*x(43)
    %(44) warfarin plasma compartment
    ka*x(43)/v-ke*x(44)
    %(45) Protein S
    pPS*x(52)-x(25)*x(45)/cAPCPS-dPS*x(45)
    %(46) APC:PS
    x(25)*x(45)/cAPCPS-dAPCPS*x(46)
    %(47) integral fibrin
    x(16)*60*60
    %(48) LMWH absorption site
    -ka_e*x(48)
    %(49) LMWH plasma compartment
    -x(2)*x(49)/cXaHeparin+ka_e*x(48)-k12_e*x(49)-k10_e*x(49)+k21_e*x(50)-dHeparin*x(49)-x(7)*x(49)/cIIaHeparin-x(10)*x(49)/cIXaHeparin
    %(50) LMWH peripheral compartment
    +k12_e*x(49)-k21_e*x(50)
    %(51) vitamin K
    -dVK_3*x(51)*(1-Imax*x(44)/(IC50+x(44)))-dVK_1*x(51)-dVK_2*x(51)+pVK+dVKO*x(53)*(1-Imax*x(44)/(IC50+x(44)))-k12_vk*x(51)+k21_vk*x(54)/vc_vk
    %(52) VKH2
    dVK_3*x(51)*(1-Imax*x(44)/(IC50+x(44)))+dVK_2*x(51)-dVKH2*x(52)
    %(53) VKO
    dVKH2*x(52)-dVKO*x(53)*(1-Imax*x(44)/(IC50+x(44)))
    %(54) VK_peripheral
    k12_vk*x(51)*vc_vk-k21_vk*x(54)
    %(55) XII
    pXII-dXII*x(55)-vXII_ca*x(58)*x(55)/(kXII_ca+x(58))-vXII_K*x(57)*x(55)/(kXII_K+x(57))
    %(56) Prekallikrein
    -vPK_XIIa*x(30)*x(56)/(kPK_XIIa+x(30))+pPK-dPK*x(56)
    %(57) Kallikrein
    vPK_XIIa*x(30)*x(56)/(kPK_XIIa+x(30))-dK*x(57)
    %(58) CA
    -dca*x(58)
    %(59) UFH
    -dHeparin*x(59)-ke_UFH*x(59)+R_UFH*INF_UFH-x(10)*x(59)/cIXaUFH-x(2)*x(59)/cXaUFH-x(7)*x(59)/cIIaUFH
    %(60) Tiger snake venom absorption site
    -ka_tiger*x(60)
    %(61) Tiger snake venom plasma compartment
    ka_tiger*x(60)-dTiger*x(61)*antivenom
    %(62) Brown snake venom plasma compartment
    ka_brown*x(28)-dBrown*x(62)*antivenom];