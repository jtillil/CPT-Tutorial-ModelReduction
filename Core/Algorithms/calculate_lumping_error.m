function [error, X_current] = calculate_lumping_error(model, options)
%%% take a model and a lumping matrix to calculate the lumping error

options.prelumpmat = options.lumpmat;
if isfield(options, 'invlumpmat')
    options.preinvlumpmat = options.invlumpmat;
end

out_state = find(options.lumpmat(:, model.I.output) ~= 0);

[~,X_current] = simModel(model.t_ref, model.X0, model.par, model.I, model.param, model.multiple, model.odefun, model.jacfun, options);

try
    error = relativeErrorL2(model.t_ref, model.X_ref(:, model.I.output), X_current(:, out_state));
catch
    error = Inf;
end

function Error=relativeErrorL2(t,X,Y)
Error=sqrt(trapz(t,(X-Y).^2)/trapz(t,X.^2));
end

end