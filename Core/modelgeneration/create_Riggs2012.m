addpath(genpath('../../Core'))

% run original script
run("Riggs2012_calciumBoneModel_MAIN.m")

% convert old nomenclature to new
model.X_ref = model.x_ref;
model.I.nstates = model.I.nrOfStates;
model.I = rmfield(model.I, "nrOfStates");
model.I.nmstate = model.I.stateName;
model.I = rmfield(model.I, "stateName");
model.I.npar = model.I.nrOfPar;
model.I = rmfield(model.I, "nrOfPar");
model.I.nmpar = model.I.parName;
model.I = rmfield(model.I, "parName");

% run minimal model script
model.odefun = @(t,X,par,model) Riggs2012_calciumBoneModel_ode(t,X,par,model);
model_minimal = model2minimal(model);