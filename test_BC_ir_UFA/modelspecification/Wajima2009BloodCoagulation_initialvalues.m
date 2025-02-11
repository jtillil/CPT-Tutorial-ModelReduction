%%% Version: 19 Jun 2022
%%%
%%% X0 = <MODELNAME>_initialValues(model)
%%%
%%% This function assigns the initial values for all state variables 
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
X0(I.CVenom_Tiger) = 0; X0(I.AT_III_UFH) = 0;             

 if contains(model.scenario,'in_vitro')
    %%% in the invitro setting Tmod was assumed to be zero Wajima (2009)
    X0(I.Tmod) = 0;
    
    %%% and the delution of the blood was considered by computing 1/3 of the
    %%% orignal initial value
    X0 = 1/3 * X0;
end

end
