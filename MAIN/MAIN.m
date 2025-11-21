%% Setup
addpath(genpath("../."))

%% Analyze small parallel pathways example model

run("./MAIN_parallel_pathways.m")

%% Analyze large QSP model of blood coagulation

run("./MAIN_blood_coagulation.m")

%% Generate figures

run("./figures_example_models.m")
run("./figures_blood_coagulation.m")