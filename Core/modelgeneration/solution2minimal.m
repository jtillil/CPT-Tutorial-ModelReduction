function model = solution2minimal(redmodel)

% read standard components from model
model.name = redmodel.name;
model.scenario = redmodel.scenario;
model.I = redmodel.I;
model.t_ref = redmodel.t_ref;
model.X_ref = redmodel.X_ref;
for idx = {'ir','contr','obs','env','pss','pneg','cneg'}
    model.(idx{1}) = redmodel.(idx{1});
end

% convert mor for R
model.nsteps = redmodel.exhaustive_mor.nsteps;
model.classifs = {'dyn','env','cneg','pneg','pss','irenv','irenv_arith','irenv_geom'};
model.configs = zeros(size(redmodel.exhaustive_mor.configs));
for i = 1:size(model.configs, 1)
for j = 1:size(model.configs, 2)
    model.configs(i,j) = find(strcmp(model.classifs, redmodel.exhaustive_mor.configs(i, j)));
end
end

% read graph components
model.adjmat = redmodel.adjmat;
model.speciesmask = redmodel.speciesmask;
model.substratemask = redmodel.substratemask;

save("minimalsolution.mat", 'model')

end