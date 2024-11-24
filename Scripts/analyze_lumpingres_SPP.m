%% own lumping

res = lumpingres_BCSnake;

errors = res.errors;
lumpmats = res.lump_matrices;
inv_lumpmats = res.inv_lump_matrices;
out_states = res.out_states;

idx = find(errors < 0.1, 1, 'last');

error = errors(idx);
lumpmat = lumpmats{idx};
inv_lumpmat = inv_lumpmats{idx};
out_state = out_states(idx);

for i = 1:size(lumpmat, 1)
    disp(sum(lumpmat(i, :)))
end

lumpedstates = cell([size(lumpmat, 1), 1]);
for i = 1:size(lumpmat, 2)
    idx_lumpmat = find(lumpmat(:, i) == 1, 1);
    lumpedstates{idx_lumpmat} = [lumpedstates{idx_lumpmat} i];
end

options = odeset;
% options.Jacobian = @(t,X) jac_lumping(X, par, model);
options.AbsTol = 1;
options.RelTol = 1e-3;
% options.InitialStep = 1e-2;
options.NonNegative = 1:size(lumpmat, 1);
% [~,lumped_out] = ode15s(@(t,X) ode_lumping(X,model.par,model,lumpmat,inv_lumpmat), model.t_ref, lumpmat*model.X0, options);

% plot(model.t_ref, lumped_out(:, out_state))

%% Duffull lumping result

% check ode of BC model
Xsym = sym(model.I.nmstate)';
psym = sym(model.I.nmpar)';
disp(model.odefun(Xsym, psym))

idxFg = model.I.Fg;
idxIIa = model.I.IIa;
idxAvenom = model.I.AVenom;
idxPvenom = model.I.CVenom;

lumpmat_Duffull = zeros(5, model.I.nstates-1);
lumpmat_Duffull(5, :) = 1;
lumpmat_Duffull(1, idxFg) = 1;
lumpmat_Duffull(5, idxFg) = 0;
lumpmat_Duffull(2, idxIIa) = 1;
lumpmat_Duffull(5, idxIIa) = 0;
lumpmat_Duffull(3, idxAvenom) = 1;
lumpmat_Duffull(5, idxAvenom) = 0;
lumpmat_Duffull(4, idxPvenom) = 1;
lumpmat_Duffull(5, idxPvenom) = 0;

invlumpmat_Duffull = pinv(lumpmat_Duffull);

options.NonNegative = 1:5;
X0_Duffull = model.X0(1:(end-1));
init = lumpmat_Duffull*X0_Duffull;
% init = init(1:(end-1));
[~,lumped_out_Duffull] = ode15s(@(t,X) ode_lumping(X,model.par,model,lumpmat_Duffull,invlumpmat_Duffull), model.t_ref, init, options);

hold on;
grid on;
plot(model.t_ref, model.X_ref(:, idxFg))
plot(model.t_ref, lumped_out_Duffull(:, 1))
hold off;

%% helper functions

function dX = ode_lumping(X, par, model, lumpmat, invlumpmat)
    % X = X(1:(end-1));
    odeout = model.odefun(invlumpmat * X, par);
    odeout = odeout(1:(end-1));
    dX = lumpmat * odeout;
end