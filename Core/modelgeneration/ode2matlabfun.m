function [odefun, ode, jacfun, jac] = ode2matlabfun(bigmodel)

%% read bigmodel
model = bigmodel;

%% generate matlab functions 
% create symbolic variables
syms t
par_sym = cell2sym(bigmodel.I.nmpar(:));
X_sym = cell2sym(bigmodel.I.nmstate(:));

% generate ode matlabfun
tic
odefun_symbolic = bigmodel.odefun(t,X_sym,par_sym,bigmodel);
disp(toc)
% disp(odefun_symbolic)
tic
odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
ode = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
% model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym,t},'File',['../modelfiles/' model.name '_ode.m'],'Optimize',true);
disp(toc)

% check jac presence
% jac_present = true;
% if ~isfield(bigmodel, 'jac')
%     jac_present = false;
% elseif isempty(bigmodel.jac)
%     jac_present = false;
% end

% always calculate jacfun
jac_present = false;

% generate jac matlabfun
tic
if jac_present
    jacfun_symbolic = bigmodel.jacfun(t,X_sym,par_sym,bigmodel);
else
    jacfun_symbolic = sym(zeros(model.I.nstates, model.I.nstates));
    for i = 1:model.I.nstates
        for j = 1:model.I.nstates
            jacfun_symbolic(i, j) = diff(odefun_symbolic(i), X_sym(j));
        end
    end
end
disp(toc)
tic
jacfun = matlabFunction(jacfun_symbolic,'Vars',{X_sym,par_sym,t});
jac = matlabFunction(jacfun_symbolic,'Vars',{X_sym,par_sym,t});
% model.jacfun = matlabFunction(jacfun_symbolic,'Vars',{X_sym,par_sym,t},'File',['../modelfiles/' model.name '_jac.m'],'Optimize',true);
disp(toc)

end