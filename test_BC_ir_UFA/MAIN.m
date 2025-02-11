% scenario = 'in_vivo_warfarin';
% 
% model = simulate_and_reduce(scenario);

model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
[model.ir, model.contr, model.obs]  =  compute_ir_indices(model);