addpath(genpath('../../Core'))

Dano = false;

if ~Dano
    % run original script
    run("Hynne2001Glycolysis_MAIN.m")

    % convert old nomenclature to new
    model.X_ref = model.X;
    model.t_ref = model.t;
    model.I.nstates = model.I.nrOfStates;
    model.I = rmfield(model.I, "nrOfStates");
    model.I.nmstate = model.I.stateName;
    model.I = rmfield(model.I, "stateName");
    model.I.npar = model.I.nrOfPar;
    model.I = rmfield(model.I, "nrOfPar");
    model.I.nmpar = model.I.parName;
    model.I = rmfield(model.I, "parName");
    
    % add needed fields
    model.I.input = model.I.('GlcX0');
    model.I.output = model.I.('Pyr');
    model.L.nconlaw = 2;
    model.L.nmconlaw = {"ATP"; "NADH"};
    model.L.ATP.states = [3 5 22];
    model.L.NADH.states = [10 12];
    
    % run minimal model script
    model.odefun = @(t,X,par,model) Hynne2001Glycolysis_ode(t,X,par,model);
    model_minimal = model2minimal(model);
    
    % add config to I
    model_minimal.I = config2I(model_minimal.I, repmat("dyn", [1 model_minimal.I.nstates]), 0);
    
    % add irenv
    [ir, contr, obs] = compute_ir_indices(model_minimal,false);
    model_minimal.ir = ir;
    model_minimal.contr = contr;
    model_minimal.obs = obs;
    model_irenv = model2irenv(model_minimal);
else
    load("../modelfiles/modelGlyc_full.mat")
    model.odefun = @(t,X,par,model) Hynne2001Glycolysis_ode_20E6D(t,X,par,model);%% generate matlab functions 
    % create symbolic variables
    syms t
    par_sym = cell2sym(model.I.nmpar(:));
    X_sym = cell2sym(model.I.nmstate(:));
    
    % generate ode matlabfun
    tic
    odefun_symbolic = model.odefun(t,X_sym,par_sym,model);
    disp(toc)
    % disp(odefun_symbolic)
    tic
    model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
    % model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym,t},'File',['../modelfiles/' model.name '_ode.m'],'Optimize',true);
    disp(toc)
    
    % generate jac matlabfun
    tic
    jacfun_symbolic = sym(zeros(model.I.nstates, model.I.nstates));
    for i = 1:model.I.nstates
        for j = 1:model.I.nstates
            jacfun_symbolic(i, j) = diff(odefun_symbolic(i), X_sym(j));
        end
    end
    disp(toc)
    tic
    model.jacfun = matlabFunction(jacfun_symbolic,'Vars',{X_sym,par_sym,t});
    % model.jacfun = matlabFunction(jacfun_symbolic,'Vars',{X_sym,par_sym,t},'File',['../modelfiles/' model.name '_jac.m'],'Optimize',true);
    disp(toc)
    
    % param struct
    model.param.states2Ipar = {};
    model.param.states2nmpar = {};
    model.param.nreact = model.I.npar;
    
    I = model.I;
    paridx = 1:length(model.par);
    
    tic
    for k = 1:I.nstates
        states2param = [];
        current_symvars = symvar(odefun_symbolic(k));
        for current_symvar = current_symvars
            par_having_current_symvar = has(par_sym, current_symvar);
            if any(par_having_current_symvar)
                states2param = [states2param, paridx(par_having_current_symvar)];
            end
        end
        model.param.states2Ipar{k} = states2param;
        model.param.states2nmpar{k} = [I.nmpar(states2param)];
    end
    disp(toc)
end

% save model
if Dano
    save("../modelfiles/modelGlyc_full_20E6D.mat", "model");
    model = rmfield(model, "ir");
    model = rmfield(model, "contr");
    model = rmfield(model, "obs");
    save("../modelfiles/modelGlyc_minimal_20E6D.mat", "model");
else
    model = model_irenv;
    save("../modelfiles/modelGlyc_full.mat", "model");
    model = rmfield(model, "ir");
    model = rmfield(model, "contr");
    model = rmfield(model, "obs");
    save("../modelfiles/modelGlyc_minimal.mat", "model");
end

