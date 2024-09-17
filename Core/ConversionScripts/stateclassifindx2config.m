addpath(genpath('modelfiles'))

% load model
model = load('modelEGFR_full.mat').model;
I = model.I;
analysis = model.analysis;

% collect indices in single matrix
indices = zeros(I.nstates, 4);
indx = ["env" "pss" "pneg" "cneg"];
for i = 1:4
    indices(:, i) = analysis.(indx(i)).max_nindex;
end

% convert indices to config
config = strings(1, I.nstates);
for state = 1:I.nstates
    [val, best_classif] = min(indices(state, :));
    if best_classif == 1
        if ismember(state, analysis.env.I_states_below_nindex_and_state_threshold)
            config(state) = "env";
        else
            config(state) = "dyn";
        end
    elseif best_classif == 2
        if ismember(state, analysis.pss.I_states_below_nindex_and_state_threshold)
            config(state) = "pss";
        else
            config(state) = "dyn";
        end
    elseif best_classif == 3
        config(state) = "pneg";
    elseif best_classif == 4
        config(state) = "cneg";
    else
        config(state) = "dyn";
    end
end

% simulate reduced model and collect objvals
model = config2model(model, config);
obj = objfun(model.t_ref, model.X0, model.par, model, config, 100);
obj
