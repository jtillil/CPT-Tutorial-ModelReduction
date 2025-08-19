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

function [tred, Xred, obj, objlog] = objfun(t, X0, par, model, config, timeout)

% write config to model indexing
model = config2model(model, config);
I = model.I;

% calculate re-configured model solution
[tred, Xred, redlog] = simModel(t, X0, par, model, timeout);

% check re-configured model solution and write to log
% TODO

% define error function
errfun = @(t_ref,X_ref,X_red) sqrt( trapz(t_ref,(X_ref-X_red).^2,1) ) ./ sqrt( trapz(t_ref,X_ref.^2,1) );

% calculate re-configured model error compared to full model
err = errfun(model.t_ref, model.X_ref, Xred);

% DEBUG
% disp(err);

% save obj values
obj.ndyn = length(I.dyn);
obj.npss = length(I.pss);
obj.err = err;
obj.errout = err(I.output);

obj.errdyn90 = prctile(err(I.dyn), 90);
obj.errdyn95 = prctile(err(I.dyn), 95);
obj.errdyn100 = max(err(I.dyn));

obj.errdynenvpss90 = prctile(err([I.dyn I.env I.pss]), 90);
obj.errdynenvpss95 = prctile(err([I.dyn I.env I.pss]), 95);
obj.errdynenvpss100 = max(err([I.dyn I.env I.pss]));

% finish log
objlog = redlog;

end
