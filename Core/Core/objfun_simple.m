function [obj, objlog, tred, Xred, err] = objfun_simple(model, config, errtype)

[obj, objlog, tred, Xred, err] = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, errtype);

end
