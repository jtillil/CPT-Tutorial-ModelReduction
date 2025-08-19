% addpath(genpath('modelfiles'))

name = 'BC';

scenario_present = false;
scenario = 'no_crosstalk';

tic
if scenario_present
    bigmodel = load(['model' name '_' scenario '_full.mat']).model;
else
    bigmodel = load(['model' name '_full.mat']).model;
end
disp(toc)

bigmodel.name = 'Wajima2009BloodCoagulation';

model = model2minimal(bigmodel);

if scenario_present
    save(['modelfiles/model' name '_' scenario '_minimal.mat'], "model");
else
    save(['modelfiles/model' name '_minimal.mat'], "model");
end


