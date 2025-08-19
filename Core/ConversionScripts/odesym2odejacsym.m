function odejacsym = odesym2odejacsym(dX_sym, X_sym, model)

I = model.I;
odejacsym = sym(zeros(I.nstates, I.nstates));

for i = 1:I.nstates
    for j = 1:I.nstates
        odejacsym(i, j) = diff(dX_sym(i), X_sym(j));
    end
end

end