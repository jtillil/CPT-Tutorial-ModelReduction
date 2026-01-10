function [tout, Xout, log, simtime] = simModel_simple(model, config)

if nargin == 1
    config = repmat("dyn", [1, model.I.nstates]);
end

[tout, Xout, log, simtime] = simModel(model.t_ref, model.X0, model.par, config2I(model.I, config, model.L), model.param, model.multiple, model.odefun, model.jacfun);

end
