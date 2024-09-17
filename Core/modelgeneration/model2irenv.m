function model = model2irenv(model)

%% check if everything necessary is present in the model
if ~isfield(model, 'obs')
    error('No model.obs index present.')
end

%% calculate irenv initial states from obs index

% calculate timepoint contribution to trapz
% this ignores the last timepoint which is anyways irrelevant
trapzcontrib = diff(model.t_ref);

% weigh obs index by timepoint contribution
wobs = model.obs.index(1:(end-1), :) .* trapzcontrib;

% normalize obs index to obtain weights that sum to 1
nobs = wobs ./ sum(wobs);
nobs_geom = wobs(2:end, :) ./ sum(wobs(2:end, :));
nobs_geom_shorter = wobs(3:end, :) ./ sum(wobs(3:end, :));

% set NaN columns to normalized trapz contribution (that indicates that the obs index is 0 everywhere)
ntrapzcontrib = trapzcontrib ./ sum(trapzcontrib);
ntrapzcontrib_geom = trapzcontrib(2:end, :) ./ sum(trapzcontrib(2:end, :));
ntrapzcontrib_geom_shorter = trapzcontrib(3:end, :) ./ sum(trapzcontrib(3:end, :));
for col = 1:size(nobs, 2)
    if any(isnan(nobs(:, col)))
        nobs(:, col) = ntrapzcontrib;
        nobs_geom(:, col) = ntrapzcontrib_geom;
        nobs_geom_shorter(:, col) = ntrapzcontrib_geom_shorter;
    end
end

% obtain arithmetic mean
model.I.Irenv_arith = sum(model.X_ref(1:(end-1), :) .* nobs);

% obtain geometric mean
irenv_geom_full = prod(model.X_ref(1:(end-1), :) .^ nobs);
irenv_geom_short = prod(model.X_ref(2:(end-1), :) .^ nobs_geom);
irenv_geom_shorter = prod(model.X_ref(3:(end-1), :) .^ nobs_geom_shorter);
irenv_geom_short(irenv_geom_short == 0) = irenv_geom_shorter(irenv_geom_short == 0);
irenv_geom_full(irenv_geom_full == 0) = irenv_geom_short(irenv_geom_full == 0);
model.I.Irenv_geom = irenv_geom_full;

end
