function [error, X_current] = calculate_lumping_error(model, options)
%%% take a model and a lumping matrix to calculate the lumping error
options.prelumpmat = options.lumpmat;
if isfield(options, 'invlumpmat')
    options.preinvlumpmat = options.invlumpmat;
end

% par = model.par;
% model.lumping.lumpmat = lump_matrix;
% model.lumping.invlumpmat = pinv(lump_matrix);
out_state = find(options.lumpmat(:, model.I.output) ~= 0);

% X0_new = lumpmat * model.X0;

% options = odeset;
% % options.Jacobian = @(t,X) jac_lumping(X, par, model);
% options.AbsTol = 1;
% options.RelTol = 1e-3;
% % options.InitialStep = 1e-2;
% options.NonNegative = 1:size(lumpmat, 1);

% [~,X_current] = ode15s(@(t,X) ode_lumping(X,par,model),model.t_ref,X0_new,options);
[~,X_current] = simModel(model.t_ref, model.X0, model.par, model.I, model.param, model.multiple, model.odefun, model.jacfun, options);

try
    error = relativeErrorL2(model.t_ref, model.X_ref(:, model.I.output), X_current(:, out_state));
catch
    error = Inf;
end

function Error=relativeErrorL2(t,X,Y)
Error=sqrt(trapz(t,(X-Y).^2)/trapz(t,X.^2));
end

% function dX = ode_lumping(X, par, model)
%     dX = model.lumping.lumpmat * model.odefun(model.lumping.invlumpmat * X,par);
% end
% 
% function dX2 = jac_lumping(X, par, model)
%     dX2 = model.lumping.lumpmat * model.jacfun(model.lumping.invlumpmat * X,par);
% end

end