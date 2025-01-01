%% Setup
addpath(genpath("../../CPT-Tutorial-ModelReduction"))

config.ir       = true;
config.non_ir   = true;
config.analyze  = true;

%% parallel pathways
name = 'SimpleParallelPathways';
namesimple = 'SPP';
scenario = 'no_crosstalk';

setupModel(name, namesimple, scenario, config)

scenario = 'with_crosstalk';

setupModel(name, namesimple, scenario, config)

%% enzyme kinetics
name = 'MMEnzymeKinetics';
namesimple = 'MMEK';
scenario = 'Cpss_Eenv';

setupModel(name, namesimple, scenario, config)

scenario = 'Cpss_EpCenv';

setupModel(name, namesimple, scenario, config)

%% blood coagulation
model = struct;
model.name = 'Gulati2014BloodCoagulation';
model.namesimple = 'BC_Gulati2014';
model.scenario = 'in_vivo';

model.I = Gulati2014BloodCoagulation_indexing();
model.I = config2I(model.I, repmat("dyn", [1, model.I.nstates]), 0);
model.X0 = Gulati2014BloodCoagulation_initialvalues(model);
model.par = Gulati2014BloodCoagulation_parameters(model);
model.ode = @(t, X) Gulati2014BloodCoagulation_ode(t,X,par,model);

% calculate indices
if config.ir
    [ir, contr, obs] = compute_ir_indices(model, false);
    model.ir = ir;
    model.contr = contr;
    model.obs = obs;
end
if config.non_ir
    model.env = compute_non_ir_indices(model, "env", false);
    model.pneg = compute_non_ir_indices(model, "pneg", false);
    model.cneg = compute_non_ir_indices(model, "cneg", false);
    model.pss = compute_non_ir_indices(model, "pss", false);
end
if config.analyze
    % model = compute_and_analyse_indices(model, 'compute');
    model.threshold = 0.1;
    model.analysis = analyse_all_indices(model);
end

% add relevant components
model = model2minimal(model);

% save model
save(['../Core/modelfiles/model' namesimple '_' scenario '_full.mat'], 'model');
