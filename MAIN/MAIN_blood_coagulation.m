%% Setup
clear; clc;
addpath(genpath("../../../CPT-Tutorial-ModelReduction"))
reduced_errors = struct;

%% Lumping as in Gulati 2014
% load model file
load("modelBC_Gulati2014_in_vivo_full.mat")
model.multiple.multiple = 0;

idxFg = model.I.Fg;
idxIIa = model.I.IIa;
idxAvenom = model.I.AVenom;
idxPvenom = model.I.CVenom;

lumpmat_Gulati = zeros(5, model.I.nstates-1);
lumpmat_Gulati(5, :) = 1;
lumpmat_Gulati(1, idxFg) = 1;
lumpmat_Gulati(5, idxFg) = 0;
lumpmat_Gulati(2, idxIIa) = 1;
lumpmat_Gulati(5, idxIIa) = 0;
lumpmat_Gulati(3, idxAvenom) = 1;
lumpmat_Gulati(5, idxAvenom) = 0;
lumpmat_Gulati(4, idxPvenom) = 1;
lumpmat_Gulati(5, idxPvenom) = 0;

reduced_errors.lumping_Gulati = calculate_lumping_error(model, lumpmat_Gulati);

%% Save results
save("./results/errors_BC_Gulati2014.mat", "reduced_errors")
