%% Setup
clear; clc;
addpath(genpath("../../CPT-Tutorial-ModelReduction"))
reduced_errors = struct;

%% Calculate ir matlabfun

% load model
% load("modelBC_SV40_from_JKn_2024.mat")
% config = repmat("dyn", [1, model.I.nstates]);
% model.I = config2I(model.I, config, []);
% 
% [ir, contr, obs] = compute_ir_indices_matlabfun(model);
% 
% plot(model.t_ref, ir.nindex(:, model.I.output))

%% Reduce model index analysis

% load model
load("modelBC_SV40_from_JKn_2024.mat")

% initialize config and threshold
config = repmat("dyn", [1, model.I.nstates]);
t = 0.05;

% assign config based on index analysis results
for i = 1:model.I.nstates
neg = 0;
ir = max(model.ir.nindex(:, i));
env = max(model.env.nindex(:, i));
envrel = max(model.env.relstateerr(:, i));
pss = max(model.pss.nindex(:, i));
pssrel = max(model.pss.relstateerr(:, i));
cneg = max(model.cneg.nindex(:, i));
pneg = max(model.pneg.nindex(:, i));
if ir <  t
if env < t || pss < t
    if env <= pss && envrel < t
        config(i) = "env";
    elseif env <= pss && pss < t && pssrel < t
        config(i) = "pss";
    elseif pss <= env && pssrel < t
        config(i) = "pss";
    elseif pss <= env && env < t && envrel < t
        config(i) = "env";
    else
        neg = 1;
    end
else
    neg = 1;
end
if neg == 1
    if cneg < t || pneg < t
        if pneg <= cneg
            config(i) = "pneg";
        elseif cneg <= pneg
            config(i) = "cneg";
        end
    end
end
end
end

% calculate error
multiple.multiple = 0;
[err_index, ~, tred, Xred] = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, [], model.param, multiple, model.odefun, model.jacfun, config, "MRSE");

% show results
disp(err_index.errout)
disp(sum(config == "dyn"))
disp(sum(config == "env"))
disp(sum(config == "pss"))
disp(sum(config == "pneg"))
disp(sum(config == "cneg"))

% calculate reduced model ir indices
model.I = config2I(model.I, config, []);
[irred.ir, irred.contr, irred.obs, irred.t_ir] = compute_ir_indices_matlabfun(model);

% save indices
save("results_irred_BC_40h_JKn_t005.mat", "irred")

% nir-indices
size = 12;
lw = 1;
lwt = 0.5;

figure
plot(irred.t_ir, irred.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold), 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
% xlim([-0.14 4.15])
xlim([-0.01 1])
ylim([-0.01 1])
legend([model.analysis.ir.nmstates_above_nindex_threshold; 'threshold'], 'Location','east')
xlabel("t [h]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_JKn_ir_index_red_0_1.pdf")


%% Blood Coagulation 40h: lumping Aarons
load("modelBC_SV40_from_JKn_2024.mat")

[lump_matrices,inv_lump_matrices,errors,out_states] = lumping_Aarons(model);

save("lumpingres_BC_SV40_from_JKn.mat", "out_states", "errors", "inv_lump_matrices", "lump_matrices")

%% analyze lumping
idx = 1:length(errors);
idx = idx(errors < 0.05);
idx = idx(end);

lumpmat = lump_matrices{idx};

for k = 1:size(lumpmat, 1)
    if sum(lumpmat(k, :)) == 1
        disp(k)
        disp(model.I.nmstate{logical(lumpmat(k, :))})
    end
end

disp(model.I.nmstate{logical(lumpmat(out_states(idx), :))})
% disp(model.I.nmstate{logical(lumpmat(out_states(idx), :))})

disp(errors(idx))

%% save reduced solution plot
size = 12;
lw = 1;
lwt = 0.5;

figure
grid on
semilogy(tred, Xred(:, model.analysis.ir.I_sorted_max_nindex_above_threshold), 'LineWidth', lw) %DisplayName', plotnames(i))
xlim([-2 42])
ylim([1e-7 5e4])
legend(model.analysis.ir.nmstates_above_nindex_threshold, 'Location','southeast')
xlabel("t [h]")
ylabel("concentration [g/L]")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_ref_sol_indices.pdf")

%% Full model relevant states

X_full = model.X_ref(:, [I.Fg, I.II, I.IIa, I.AVenom, I.CVenom, I.Xa, I.Xa_Va, I.Tmod, I.AT_III_Heparin, I.TaipanVenom, I.CVenom_Tiger]);

%% Lumping as in Gulati 2014
% load model file
% load("modelBC_Gulati2014_in_vivo_full.mat")
% load("modelBCSV_minimal.mat")
load("modelBC_temp.mat")
model.X0 = Gulati2014BloodCoagulation_initialvalues(model);
model.multiple.multiple = 0;
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);

idxFg = model.I.Fg;
idxII = model.I.II;
idxIIa = model.I.IIa;
idxAvenom = model.I.AVenom;
idxPvenom = model.I.CVenom;

lumpmat_Gulati = zeros(5, model.I.nstates);
lumpmat_Gulati(5, :) = 1;
lumpmat_Gulati(1, idxFg) = 1;
lumpmat_Gulati(5, idxFg) = 0;
lumpmat_Gulati(2, idxIIa) = 1;
lumpmat_Gulati(5, idxIIa) = 0;
lumpmat_Gulati(3, idxAvenom) = 1;
lumpmat_Gulati(5, idxAvenom) = 0;
lumpmat_Gulati(4, idxPvenom) = 1;
lumpmat_Gulati(5, idxPvenom) = 0;

opt.lumpmat = lumpmat_Gulati;

% lumpmat_Gulati = zeros(6, model.I.nstates);
% lumpmat_Gulati(5, :) = 1;
% lumpmat_Gulati(1, idxFg) = 1;
% lumpmat_Gulati(5, idxFg) = 0;
% lumpmat_Gulati(6, idxII) = 1;
% lumpmat_Gulati(5, idxII) = 0;
% lumpmat_Gulati(2, idxIIa) = 1;
% lumpmat_Gulati(5, idxIIa) = 0;
% lumpmat_Gulati(3, idxAvenom) = 1;
% lumpmat_Gulati(5, idxAvenom) = 0;
% lumpmat_Gulati(4, idxPvenom) = 1;
% lumpmat_Gulati(5, idxPvenom) = 0;

[reduced_errors.lumping_Gulati, X_red] = calculate_lumping_error(model, opt);

AUC = trapz(model.t_ref, model.X_ref, 1);
opt.invlumpmat = opt.lumpmat';
for col = 1:size(opt.invlumpmat, 2)
    opt.invlumpmat(:, col) = (AUC.*opt.invlumpmat(:, col)') ./ (AUC*opt.invlumpmat(:, col));
end
[reduced_errors.scaled_lumping_Gulati, X_red] = calculate_lumping_error(model, opt);

% Save results
save("./results/errors_BC_Gulati2014.mat", "reduced_errors")

%% Diagnostics
I = model.I;
par = model.par;

disp(par(I.degIIa))

disp(par(I.v14))
disp(par(I.k14))

disp(par(I.v15))
disp(par(I.k15))

% disp("v15")
% disp(par(I.v15))

disp("v12")
disp(par(I.v12))
disp(par(I.k12))

disp("v13")
disp(par(I.v13))
disp(par(I.k13))

disp("v14")
disp(par(I.v14))
disp(par(I.k14))

% disp(par(I.pFg)/par(I.degFg))

disp("Brown")
disp(par(I.ka_Brown))
disp(par(I.d_Brown))

disp("Fg")
disp(par(I.pFg))
disp(par(I.degFg))

disp(model.X0(I.Fg))

%% Plots
load("modelBC_temp.mat")

t = tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
nexttile(t)
plot(model.t_ref, model.X_ref(:, [model.I.AVenom model.I.CVenom model.I.II model.I.IIa model.I.Fg]), 'LineWidth', 2)
ylim([1e-4 1e4])
grid on
set(gca, 'YScale', 'log')
legend('AVenom', 'CVenom', 'II', 'IIa', 'Fibrinogen', 'Location', 'southeast')
title('Full model')

nexttile(t)
% plot(model.t_ref, X_red(:, [3 4 5 2 1]).*repmat([1 1 1/58 1 1], [length(model.t_ref), 1]), 'LineWidth', 2)
% plot(t_ref, X_red(:, [3 4 5 2 1]).*repmat([1 1 1/58 1 1], [length(t_ref), 1]), 'LineWidth', 2)
plot(t_red, X_red, 'LineWidth', 2)
ylim([1e-4 1e4])
grid on
set(gca, 'YScale', 'log')
legend('AVenom', 'CVenom', 'lumped', 'IIa', 'Fibrinogen', 'Location', 'southeast')
% title('Lumped model from our implementation')
title('Lumped model as reported in Gulati2014')

nexttile(t)
plot(t_red, (X_red(:, 5) - X_ref(:, model.I.Fg))./X_ref(:, model.I.Fg), 'LineWidth', 2)
% plot(t_ref, (X_red(:, 1) - X_ref(:, model.I.Fg))./X_ref(:, model.I.Fg), 'LineWidth', 2)
% ylim([1e-4 1e4])
grid on
% set(gca, 'YScale', 'log')
legend('relative error')
title('Relative error')

exportgraphics(t, "./figures/BClumped_reported_Gulati_ODEs.pdf")
% exportgraphics(t, "./figures/BClumped_our_implementation.pdf")

%% Get ODE form
% Uni Potsdam implementation
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));

ODE_red = simplify(lumpmat_Gulati*model.odefun(X_sym, par_sym));

%% Get ODE form
% Gulati implementation
x = sym('x', [62 1]);
syms t
syms pX pV pII pVIII pIX pFg pXIII pPg pPC pTmod pXI pVII pTFPI pPS pVK pXII pPK dX dXa dV dVa dVaXa dII dIIa dVIII dVIIIa dIXa dIXaVIIIa dIX dXIa dFg dFDP dF dXF dD dXIII dXIIIa dPg dP dPC dAPC dTmod dIIaTmod dTaipan dXI dXIIa dVII dVIIa dTAT dBrown dTF dVIITF dVIIaTF dTFPI dXaTFPI dVIIaTFXaTFPI dPS dAPCPS dHeparin dVK_1 dVK_2 dVK_3 dVKH2 dVKO dXII dPK dK dca vIXa kIXa vII2IIaXa kII2IIaXa vV2VaIIa kV2VaIIa vVIII2VIIIaIIa kVIII2VIIIaIIa vIXaVIIIa kIXaVIIIa vXIa kXIa vVaXa kVaXa vFg2F kFg2F vFg2FDP kFg2FDP vF2FDP kF2FDP vF2XF kF2XF vXF2D kXF2D vXIII2XIIIa kXIII2XIIIa vPg2PF kPg2PF vPg2PIIa kPg2PIIa vPC2APC kPC2APC vVIIIaLoss kVIIIaLoss vVaLoss kVaLoss vXaVaLoss kXaVaLoss vXF2DAPC kXF2DAPC vPg2PAPC kPg2PAPC vXI2XIaIIa kXI2XIaIIa vXI2XIaXIIa kXI2XIaXIIa vX2XaVIIa kX2XaVIIa vVII2VIIaIIa kVII2VIIaIIa vVII2VIIaTaipan kVII2VIIaTaipan vVIITF_Xa kVIITF_Xa vX_VIIaTF kX_VIIaTF vIX_VIIaTF kIX_VIIaTF vVIITF_TF kVIITF_TF vVII_Xa kVII_Xa vVII_VIIaTF kVII_VIIaTF vVII_IXa kVII_IXa vXII_ca kXII_ca vXII_K kXII_K vPK_XIIa kPK_XIIa cVaXa cIXaVIIIa cIIaTmod cVIIaTF cVIITF cVIIaTFXaTFPI cXaTFPI cAPCPS cXaHeparin cIIaHeparin cIXaHeparin cIXaUFH cXaUFH cIIaUFH Imax IC50 ka v ke lagt ka_e k10_e k12_e k21_e v2_e k12_vk k21_vk vc_vk ke_UFH R_UFH T_UFH time_of_antivenom kr FLAG1 FLAG2 KDXa KDE ka_brown ka_tiger dTiger

lumpmat_Gulati_supp = zeros(5, 62);
lumpmat_Gulati_supp(5, :) = 1;
lumpmat_Gulati_supp(1, 14) = 1;              % Fg
lumpmat_Gulati_supp(5, 14) = 0;
lumpmat_Gulati_supp(2, 7) = 1;               % IIa
lumpmat_Gulati_supp(5, 7) = 0;
lumpmat_Gulati_supp(3, 28) = 1;              % AVenom
lumpmat_Gulati_supp(5, 28) = 0;
lumpmat_Gulati_supp(4, 62) = 1;              % CVenom
lumpmat_Gulati_supp(5, 62) = 0;

ODE_red_Gulati = simplify(lumpmat_Gulati_supp * Gulati_Supplement_ode(t,x,pX,pV,pII,pVIII,pIX,pFg,pXIII,pPg,pPC,pTmod,pXI,pVII,pTFPI,pPS,pVK,pXII,pPK,dX,dXa,dV,dVa,dVaXa,dII,dIIa,dVIII,dVIIIa,dIXa,dIXaVIIIa,dIX,dXIa,dFg,dFDP,dF,dXF,dD,dXIII,dXIIIa,dPg,dP,dPC,dAPC,dTmod,dIIaTmod,dTaipan,dXI,dXIIa,dVII,dVIIa,dTAT,dBrown,dTF,dVIITF,dVIIaTF,dTFPI,dXaTFPI,dVIIaTFXaTFPI,dPS,dAPCPS,dHeparin,dVK_1,dVK_2,dVK_3,dVKH2,dVKO,dXII,dPK,dK,dca,vIXa,kIXa,vII2IIaXa,kII2IIaXa,vV2VaIIa,kV2VaIIa,vVIII2VIIIaIIa,kVIII2VIIIaIIa,vIXaVIIIa,kIXaVIIIa,vXIa,kXIa,vVaXa,kVaXa,vFg2F,kFg2F,vFg2FDP,kFg2FDP,vF2FDP,kF2FDP,vF2XF,kF2XF,vXF2D,kXF2D,vXIII2XIIIa,kXIII2XIIIa,vPg2PF,kPg2PF,vPg2PIIa,kPg2PIIa,vPC2APC,kPC2APC,vVIIIaLoss,kVIIIaLoss,vVaLoss,kVaLoss,vXaVaLoss,kXaVaLoss,vXF2DAPC,kXF2DAPC,vPg2PAPC,kPg2PAPC,vXI2XIaIIa,kXI2XIaIIa, vXI2XIaXIIa,kXI2XIaXIIa,vX2XaVIIa,kX2XaVIIa,vVII2VIIaIIa,kVII2VIIaIIa,vVII2VIIaTaipan,kVII2VIIaTaipan,vVIITF_Xa,kVIITF_Xa,vX_VIIaTF,kX_VIIaTF,vIX_VIIaTF,kIX_VIIaTF,vVIITF_TF,kVIITF_TF,vVII_Xa,kVII_Xa,vVII_VIIaTF,kVII_VIIaTF,vVII_IXa,kVII_IXa,vXII_ca,kXII_ca,vXII_K,kXII_K,vPK_XIIa,kPK_XIIa,cVaXa,cIXaVIIIa,cIIaTmod,cVIIaTF,cVIITF,cVIIaTFXaTFPI,cXaTFPI,cAPCPS,cXaHeparin,cIIaHeparin,cIXaHeparin,cIXaUFH,cXaUFH,cIIaUFH,Imax,IC50,ka,v,ke,lagt,ka_e,k10_e,k12_e,k21_e,v2_e,k12_vk,k21_vk,vc_vk,ke_UFH,R_UFH,T_UFH,time_of_antivenom,kr,FLAG1,FLAG2,KDXa,KDE,ka_brown,ka_tiger,dTiger));

%% modify Wajima2009 equations to match Gulati2014 paper ODE diagram

load("modelBC_temp.mat")
% turn off P deactivation of Fg
% model.par(model.I.v15) = 0;
% model.par(model.I.v13) = 0;
% model.par(model.I.v12) = 1940;
% model.par(model.I.k12) = 1.38;
% load modified odefun with autoactivating compounds to IIa turned off
% model.odefun = @(t,X,par,model) Wajima2009BloodCoagulation_ode_Gulati_compatibility(t,X,par,model);
% convert odefun to (t,X) matlabfun
% syms t
% par_sym = cell2sym(model.I.nmpar(:));
% X_sym = cell2sym(model.I.nmstate(:));
% odefun_symbolic = model.odefun(t,X_sym,par_sym,model);
% model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
% make lumping matrix
I = model.I;
lumpmat_Gulati = zeros(5, model.I.nstates);
lumpmat_Gulati(5, :) = 1;
lumpmat_Gulati(1, I.Fg) = 1;
lumpmat_Gulati(5, I.Fg) = 0;
lumpmat_Gulati(2, I.IIa) = 1;
lumpmat_Gulati(5, I.IIa) = 0;
lumpmat_Gulati(3, I.AVenom) = 1;
lumpmat_Gulati(5, I.AVenom) = 0;
lumpmat_Gulati(4, I.CVenom) = 1;
lumpmat_Gulati(5, I.CVenom) = 0;
% check lumping error
[reduced_errors.lumping_Gulati, X_red] = calculate_lumping_error(model, lumpmat_Gulati);

load("modelBC_temp.mat")
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
[t_ref, X_ref] = simModel(model.t_ref, model.X0, model.par, model.I, [], model.multiple, model.odefun, []);
relativeErrorL2(t_ref, X_ref(:, model.I.Fg), X_red(:, 5))

%% explicitly written lumped Gulati2014 model - fitted parameters

% make explicitly lumped model
model = struct;
model.I = Gulati2014BClumped_indexing();
model.X0 = Gulati2014BClumped_initialvalues_fitted(model);
model.par = Gulati2014BClumped_parameters_fitted(model);
model.odefun = @(t,X,par,model) Gulati2014BClumped_ode(t,X,par,model);
syms t
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));
odefun_symbolic = model.odefun(t,X_sym,par_sym,model);
model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
model.multiple.multiple = 0;
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
% simulate
[t_red, X_red] = simModel([0 40], model.X0, model.par, model.I, [], model.multiple, model.odefun, []);

%% explicitly written lumped Gulati2014 model - original parameters (except for lumped state)

% make explicitly lumped model
model = struct;
model.I = Gulati2014BClumped_indexing();
model.X0 = Gulati2014BClumped_initialvalues_original(model);
model.par = Gulati2014BClumped_parameters_original(model);
model.odefun = @(t,X,par,model) Gulati2014BClumped_ode(t,X,par,model);
syms t
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));
odefun_symbolic = model.odefun(t,X_sym,par_sym,model);
model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
model.multiple.multiple = 0;
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
% simulate
[t_red, X_red] = simModel([0 40], model.X0, model.par, model.I, [], model.multiple, model.odefun, []);

load("modelBC_temp.mat")
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
[t_ref, X_ref] = simModel(t_red, model.X0, model.par, model.I, [], model.multiple, model.odefun, []);
relativeErrorL2(t_ref, X_ref(:, model.I.Fg), X_red(:, 5))

%% ?????????????

% make explicitly lumped model
model = struct;
model.I = Gulati2014BClumped_indexing();
model.X0 = Gulati2014BClumped_initialvalues_original(model);
model.par = Gulati2014BClumped_parameters_original(model);
model.odefun = @(t,X,par,model) Gulati2014BClumped_ode(t,X,par,model);
syms t
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));
odefun_symbolic = model.odefun(t,X_sym,par_sym,model);
model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
model.multiple.multiple = 0;
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
% simulate
[t_red, X_red] = simModel([0 40], model.X0, model.par, model.I, [], model.multiple, model.odefun, []);

load("modelBC_temp.mat")
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
[t_ref, X_ref] = simModel(t_red, model.X0, model.par, model.I, [], model.multiple, model.odefun, []);
relativeErrorL2(t_ref, X_ref(:, model.I.Fg), X_red(:, 5))

%% Error fun
function Error=relativeErrorL2(t,X,Y)
    Error=sqrt(trapz(t,(X-Y).^2)/trapz(t,X.^2));
end
