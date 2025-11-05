%% Setup
addpath(genpath("../../CPT-Tutorial-ModelReduction"))

% run("./setup_models.m")

%% Analyze two small example models

run("./MAIN_parallel_pathways.m")
run("./MAIN_enzyme_kinetics.m")

%% Analyze large QSP model of blood coagulation

run("./MAIN_blood_coagulation.m")

%% Generate figures

run("./figures_example_models.m")
run("./figures_blood_coagulation.m")