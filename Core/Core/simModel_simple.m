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
function [tout, Xout, log, simtime] = simModel_simple(model, config)

if nargin == 1
    config = repmat("dyn", [1, model.I.nstates]);
end

[tout, Xout, log, simtime] = simModel(model.t_ref, model.X0, model.par, config2I(model.I, config, model.L), model.param, model.multiple, model.odefun, model.jacfun);

end
