%%% Version: October 28th, 2023
%%%
%%% call by: [tout, Xout, log] = simModel(t, X0, par, model)
%%%
%%% This function calculates one solution (and returns a log of the
%%% calculation) of a model, specified by 
%%%
%%% Input:  model config
%%%
%%% Output: [tout, Xout, log]   output
%%%
%%% Citation:
%%% 
%%% ---
%%% 
%%% Authors: Johannes Tillil
%%%

% function [tout, Xout, log] = simModel(t, X0, par, model, timeout)
function [obj, objlog, tred, Xred, err] = objfun_simple(model, config, errtype)

[obj, objlog, tred, Xred, err] = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, model.L, model.param, model.multiple, model.odefun, model.jacfun, config, errtype);

end
