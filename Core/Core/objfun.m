function [obj, objlog, tred, Xred, err] = objfun(t, X_ref, X0, par, I, L, param, multiple, odefun, jacfun, config, errtype)

% write config to indexing
I = config2I(I, config(1:I.nstates), L);

% write simplification of parameters
if length(config) > I.nstates
    par(config((I.nstates+1):end) == "0") = 0;
    par(config((I.nstates+1):end) == "1") = 1;
    par(config((I.nstates+1):end) == "Inf") = Inf;
end

% check integrity of conservation laws
for p = 1:length(I.replaceODE)
    k = I.replaceODE(p);        % index of ODE to be replaced 
    states = I.replaceODEby{p}; % index of ODEs whose sum is used for replacement  

    % check, whether any of the other replaceODE states is part of one of the
    % other used replacedbyODEs (this is not allowed)
    % report the first one
    if any( ismember( setdiff(I.replaceODE,k), setdiff(states,k) ) )
        commonstates  = intersect(setdiff(I.replaceODE,k),states);
        fprintf('\n --> The state %s in replaceODE does belong to multiple sets of states in replaceODEby, continue with next conservation law.',I.nmstate{commonstates(1)});
        
        % save obj values
        obj.ndyn = length(I.dyn);
        obj.npss = length(I.pss);
        obj.err = 1e6;
        obj.errout = 1e6;
        
        obj.errdyn90 = 1e6;
        obj.errdyn95 = 1e6;
        obj.errdyn100 = 1e6;
        
        obj.errdynenvpss90 = 1e6;
        obj.errdynenvpss95 = 1e6;
        obj.errdynenvpss100 = 1e6;
        
        % finish log
        objlog = 'current conservation law overlapped with others';

        % red output
        tred = NaN;
        Xred = NaN;

        % err
        err = 1e6;

        return
    end
end

% calculate re-configured model solution
% [tred, Xred, redlog] = simModel(t, X0, par, model, timeout);
[tred, Xred, redlog, simtime] = simModel(t, X0, par, I, param, multiple, odefun, jacfun);

% check re-configured model solution and write to log
% TODO

% define error function
if errtype == "rel2NE"
    errfun = @(t_ref,X_ref,X_red) sqrt( trapz(t_ref, (X_ref - X_red).^2, 1) ) ./ sqrt( trapz(t_ref,X_ref.^2,1) );
elseif errtype == "rel2NE_conditioned"
    errfun = @(t_ref,X_ref,X_red) sqrt( trapz(t_ref, (X_ref - X_red).^2 + 1e-100, 1) ) ./ sqrt( trapz(t_ref, (X_ref).^2 + 1e-100,1) );
elseif errtype == "relAE"
    errfun = @(t_ref,X_ref,X_red) trapz(t_ref, abs(X_red - X_ref), 1) ./ trapz(t_ref, abs(X_ref), 1);
elseif errtype == "relRE"
    errfun = @(t_ref,X_ref,X_red) sqrt( trapz(t_ref, (X_ref - X_red), 1) ) ./ sqrt( trapz(t_ref,X_ref,1) );
elseif errtype == "relALE"
    errfun = @(t_ref,X_ref,X_red) trapz(t_ref, abs(log10( (X_red + 1e-100) ./ (X_ref + 1e-100) )), 1) / t_ref(end);
else
    disp(errtype)
    error("Unknown error type!")
end

% calculate re-configured model error compared to full model
err = errfun(t, X_ref, Xred);

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

obj.simtime = simtime;

% finish log
objlog = redlog;

end
