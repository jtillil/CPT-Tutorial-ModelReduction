function model = model2consts(model)

%% check if everything necessary is present in the model
if ~isfield(model, 'obs')
    error('No model.obs index present.')
end

%% calculate averages of C-t profile

% extract C-t profile
t_ref = model.t_ref;
X_ref = model.X_ref;

% calculate arithmetic average
model.I.Average = trapz(t_ref, X_ref, 1) ./ (t_ref(end) - t_ref(1));

% calculate geometric average
% TODO

%% calculate C at mode of obs index

% extract obs index
Modes = [];
obs = model.obs.index;

% calculate mode for every state
for stateID = 1:model.I.nstates
    obs_state = obs(:, stateID);
    [~, mode_index] = max(obs_state);
    Modes = [Modes X_ref(mode_index(1), stateID)];
end
model.I.Mode = Modes;

%% calculate steady-state constant concentration

% extract last concentration of simulation
model.I.Ssenv = X_ref(end,:);

end
