function odejacstr = odesym2odejacstr(handle_odesym, model)

I = model.I;
X_sym = transpose(sym(I.nmstate));
par_sym = transpose(sym(I.nmpar));

fprintf("Symbolic ode\n")
% tic
dX_sym = handle_odesym(X_sym, par_sym, model);
% toc

fprintf("Symbolic jacobian\n")
% tic
odejacsym = odesym2odejacsym(dX_sym, X_sym, model);
% toc

fprintf("String jacobian\n")
% tic
odejacstr = odejacsym2odejacstr(odejacsym, model);
% toc

end