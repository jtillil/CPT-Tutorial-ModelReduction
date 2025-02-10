%% produce Jacobian for parameter sensitivity symbolically
% in vitro
model = simulate_and_reduce_parameter_sensitivity('in_vitro_warfarin');
I=model.I;

syms t
par_sym=cell2sym(model.I.nmpar(:));
X_sym=cell2sym(model.I.nmstate(:));
symbolic_odes=model.odefun(t,X_sym,par_sym,model);

for i=1:I.nstates
    parfor j=1:I.nstates
        jac_sym(i,j)=diff(symbolic_odes(i),X_sym(j));
    end
end
jac_fun=matlabFunction(jac_sym,'Vars',{X_sym,par_sym});
save('results/jacfun_parameter_sensitivity_in_vitro.mat','jac_fun')


% in vivo
model = simulate_and_reduce_parameter_sensitivity('in_vivo_warfarin');
I=model.I;

syms t
par_sym=cell2sym(model.I.nmpar(:));
X_sym=cell2sym(model.I.nmstate(:));
% simplify a little for numerical reasons
X_sym(all(model.X_ref==0))=0;
symbolic_odes=model.odefun(t,X_sym,par_sym,model);

for i=1:I.nstates
    parfor j=1:I.nstates
        jac_sym(i,j)=diff(symbolic_odes(i),X_sym(j));
    end
end
X_sym=cell2sym(model.I.nmstate(:));
jac_fun=matlabFunction(jac_sym,'Vars',{X_sym,par_sym});
save('results/jacfun_parameter_sensitivity_in_vivo.mat','jac_fun')