%% Setup
clear; clc;
addpath(genpath("../../CPT-Tutorial-ModelReduction"))
reduced_errors = struct;

%% Lumping as in Gulati 2014
% load model file
% load("modelBC_Gulati2014_in_vivo_full.mat")
% load("modelBCSV_minimal.mat")
load("modelBC_temp.mat")
model.X0 = Gulati2014BloodCoagulation_initialvalues(model);
model.multiple.multiple = 0;

idxFg = model.I.Fg;
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

[reduced_errors.lumping_Gulati, X_red] = calculate_lumping_error(model, lumpmat_Gulati);

%% Diagnostics

disp(model.par(model.I.degIIa))

%% Plots

t = tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile(t)
plot(model.t_ref, model.X_ref(:, [model.I.AVenom model.I.CVenom model.I.II model.I.IIa model.I.Fg]), 'LineWidth', 2)
ylim([1e-4 1e4])
grid on
set(gca, 'YScale', 'log')
legend('AVenom', 'CVenom', 'II', 'IIa', 'Fibrinogen')

nexttile(t)
plot(model.t_ref, X_red(:, [3 4 5 2 1]).*repmat([1 1 1/58 1 1], [length(model.t_ref), 1]), 'LineWidth', 2)
ylim([1e-4 1e4])
grid on
set(gca, 'YScale', 'log')
legend('AVenom', 'CVenom', 'II', 'IIa', 'Fibrinogen')

%% Get ODE form
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));

ODE_red = simplify(lumpmat_Gulati*model.odefun(X_sym, par_sym))


%% Save results
save("./results/errors_BC_Gulati2014.mat", "reduced_errors")
