function error = calculate_lumping_error(model, lump_matrix)
%%% take a model and a lumping matrix to calculate the lumping error

par = model.par;
model.lumping.lumpmat = lump_matrix;
model.lumping.invlumpmat = pinv(lump_matrix);
out_state = find(lump_matrix(:, model.I.output) ~= 0);

X0_new = lump_matrix * model.X0;

options = odeset;
% options.Jacobian = @(t,X) jac_lumping(X, par, model);
options.AbsTol = 1;
options.RelTol = 1e-3;
% options.InitialStep = 1e-2;
options.NonNegative = 1:size(lump_matrix, 1);

[~,X_current] = ode15s(@(t,X) ode_lumping(X,par,model),model.t_ref,X0_new,options);
try
    error = relativeErrorL2(model.t_ref, model.X_ref(:, model.I.output), X_current(:, out_state));
catch
    error = Inf;
end

function Error=relativeErrorL2(t,X,Y)
Error=sqrt(trapz(t,(X-Y).^2)/trapz(t,X.^2));
end

function dX = ode_lumping(X, par, model)
    dX = model.lumping.lumpmat * model.odefun(model.lumping.invlumpmat * X,par);
end

function dX2 = jac_lumping(X, par, model)
    dX2 = model.lumping.lumpmat * model.jacfun(model.lumping.invlumpmat * X,par);
end

end