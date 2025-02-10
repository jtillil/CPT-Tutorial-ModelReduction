model.name = 'Wajima2009BloodCoagulation';
% namesimple = 'BC';
% model.scenario = 'in_vivo_snakevenom_1h';
model.scenario = 'in_vivo_snakevenom_40h';
% model.scenario = 'in_vitro_PTtest_lowTF';
% model.scenario = 'in_vitro_PTtest_highTF';

model = Wajima2009BloodCoagulation_model_set_up_details(model);
model = set_up_the_model(model);

model = ode2matlabfun(model);

[model.ir, model.contr, model.obs] = compute_ir_indices_matlabfun(model,saveresults);

% model = compute_and_analyse_indices_matlabfun(model,'compute');

% saveresults = false;
% model.env  = compute_non_ir_indices(model,'env',saveresults);
% model.pneg = compute_non_ir_indices(model,'pneg',saveresults);
% model.cneg = compute_non_ir_indices(model,'cneg',saveresults);
% model.pss  = compute_non_ir_indices(model,'pss',saveresults);