%% setup

model.multiple.multiple = false;

par = model.par;
num_par = 0;

% conf_idx = find(redmodel.exhaustive_mor.objvals(:, 1) == redmodel.redobj.ndyn, 1, 'last');
% conf_idx = find(redmodel.exhaustive_mor.objvals(:, 1) == 25, 1, 'last');
% config_red = redmodel.exhaustive_mor.configs(conf_idx, :);
% I_red = config2I(model.I, config_red, model.L);

lumping = 1;

if lumping
    simoptions.prelumpmat = lumpmat;
    config = repmat("dyn", [model.I.nstates, 1]);
    model.I = config2I(model.I, config, []);
else
    simoptions = struct;
end

[tred, X_red, log] = simModel(model.t_ref, model.X0, model.par, model.I, model.param, model.multiple, model.odefun, model.jacfun, simoptions);

ndyn = sum(config == "dyn");
npneg = sum(config == "pneg");
ncneg = sum(config == "cneg");
nenv = sum(config == "env");
ngeom = sum(config == "irenv_geom");
narith = sum(config == "irenv_arith");
npss = sum(config == "pss");

ntstart = 1;

%% calc

% TODO show which parameters are affected and how much
% TODO also count env and irenv for parameters and remove parameters that
%       are double due to env*par

for i = 1:length(par)
    disp(i)
    parval = par(i);
    if parval ~= 0
        % 1/1000x
        par_change = model.par;
        par_change(i) = (1e-3)*par_change(i);
        [tout, X_out, log] = simModel(model.t_ref, model.X0, par_change, model.I, model.param, model.multiple, model.odefun, model.jacfun, simoptions);
        rel_err = abs(X_out(ntstart:end, :) - X_red(ntstart:end, :)) ./ X_red(ntstart:end, :);
        max_rel_err_div1e3 = max(rel_err(~isinf(rel_err)), [], 'all');
        % half
        par_change = model.par;
        par_change(i) = 0.5*par_change(i);
        [tout, X_out, log] = simModel(model.t_ref, model.X0, par_change, model.I, model.param, model.multiple, model.odefun, model.jacfun, simoptions);
        rel_err = abs(X_out(ntstart:end, :) - X_red(ntstart:end, :)) ./ X_red(ntstart:end, :);
        max_rel_err_half = max(rel_err(~isinf(rel_err)), [], 'all');
        % double
        par_change = model.par;
        par_change(i) = 2*par_change(i);
        [tout, X_out, log] = simModel(model.t_ref, model.X0, par_change, model.I, model.param, model.multiple, model.odefun, model.jacfun, simoptions);
        rel_err = abs(X_out(ntstart:end, :) - X_red(ntstart:end, :)) ./ X_red(ntstart:end, :);
        max_rel_err_double = max(rel_err(~isinf(rel_err)), [], 'all');
        % 1000x
        par_change = model.par;
        par_change(i) = (1e3)*par_change(i);
        [tout, X_out, log] = simModel(model.t_ref, model.X0, par_change, model.I, model.param, model.multiple, model.odefun, model.jacfun, simoptions);
        rel_err = abs(X_out(ntstart:end, :) - X_red(ntstart:end, :)) ./ X_red(ntstart:end, :);
        max_rel_err_mul1e3 = max(rel_err(~isinf(rel_err)), [], 'all');
        if max([max_rel_err_half, max_rel_err_double, max_rel_err_div1e3, max_rel_err_mul1e3]) >= 1e-16
            num_par = num_par + 1;
        end
        disp(max_rel_err_half)
        disp(max_rel_err_double)
    end
end
