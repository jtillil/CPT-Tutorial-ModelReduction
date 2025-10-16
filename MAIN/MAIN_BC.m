%% Setup
clear; clc;
addpath(genpath("../."))
reduced_errors = struct;

%% Calculate ir matlabfun

% load model
% load("modelBC_SV40_from_JKn_2024.mat")
% config = repmat("dyn", [1, model.I.nstates]);
% model.I = config2I(model.I, config, []);
% 
% [ir, contr, obs] = compute_ir_indices_matlabfun(model);
% 
% plot(model.t_ref, ir.nindex(:, model.I.output))

%% Reduce model with index analysis (Jane 2024 scheme including state-approximation errors)

% load model
load("modelBC_SV40_from_JKn_2024.mat")

% initialize config and threshold
config = repmat("dyn", [1, model.I.nstates]);
t = 0.1;

errfun = @(t_ref,X_ref,X_red) sqrt( trapz(t_ref,(X_ref - X_red).^2,1) ) ./ sqrt( trapz(t_ref,X_ref.^2,1) );

% assign config based on index analysis results
for i = 1:model.I.nstates
neg = 0;
ir = max(model.ir.nindex(:, i));
env = max(model.env.nindex(:, i));

envrel = max(model.env.relstateerr(:, i));
% calc_config = repmat("dyn", [1, model.I.nstates]);
% calc_config(i) = "env";
% calc_I = config2I(model.I, calc_config, []);
% [~, X_red, ~, ~] = simModel(model.t_ref, model.X0, model.par, calc_I, model.param, model.multiple, model.odefun, model.jacfun);
% envrel = errfun(model.t_ref, model.X_ref, X_red);
% envrel = envrel(i);

pss = max(model.pss.nindex(:, i));

pssrel = max(model.pss.relstateerr(:, i));
% calc_config = repmat("dyn", [1, model.I.nstates]);
% calc_config(i) = "pss";
% calc_I = config2I(model.I, calc_config, []);
% [~, X_red, ~, ~] = simModel(model.t_ref, model.X0, model.par, calc_I, model.param, model.multiple, model.odefun, model.jacfun);
% pssrel = errfun(model.t_ref, model.X_ref, X_red);
% pssrel = pssrel(i);

cneg = max(model.cneg.nindex(:, i));
pneg = max(model.pneg.nindex(:, i));
if ir <  t
if env < t || pss < t
    if env <= pss && envrel < t
        config(i) = "env";
    elseif env <= pss && pss < t && pssrel < t
        config(i) = "pss";
    elseif pss <= env && pssrel < t
        config(i) = "pss";
    elseif pss <= env && env < t && envrel < t
        config(i) = "env";
    else
        neg = 1;
    end
else
    neg = 1;
end
if neg == 1
    if cneg < t || pneg < t
        if pneg <= cneg
            config(i) = "pneg";
        elseif cneg <= pneg
            config(i) = "cneg";
        end
    else
        disp("unclassified")
        disp(model.I.nmstate(i))
    end
end
end
end

model.I.nmstate(config == "dyn")

%%
% correct reduction
% config(model.I.PS) = "pneg";
% config(model.I.PC) = "pneg";
% config(model.I.APC) = "pneg";
% config(model.I.Tmod) = "pneg";

% calculate error
multiple.multiple = 0;
[err_index, ~, tred, Xred] = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, [], model.param, multiple, model.odefun, model.jacfun, config, "MRSE");

% show results
disp(err_index.errout)
disp(sum(config == "dyn"))
disp(sum(config == "env"))
disp(sum(config == "pss"))
disp(sum(config == "pneg"))
disp(sum(config == "cneg"))

model.I.nmstate(config == "pss")
max(model.pneg.nindex(:, model.I.TFPI))
max(model.pneg.nindex(:, model.I.VK))

%% Reduce model with index analysis (Tutorial scheme where only lowest index important + rel state err)

% load model
load("modelBC_SV40_from_JKn_2024.mat")

% initialize config and threshold
config = repmat("dyn", [1, model.I.nstates]);
t = 0.109;

errfun = @(t_ref,X_ref,X_red) sqrt( trapz(t_ref,(X_ref - X_red).^2,1) ) ./ sqrt( trapz(t_ref,X_ref.^2,1) );

% assign config based on index analysis results
idx = 1:model.I.nstates;
% idx = idx(max(model.ir.index) > 0);
for i = idx
ir = max(model.ir.nindex(:, i));
env = max(model.env.nindex(:, i));
pss = max(model.pss.nindex(:, i));
cneg = max(model.cneg.nindex(:, i));
pneg = max(model.pneg.nindex(:, i));

envrel = max(model.env.relstateerr(:, i));
pssrel = max(model.pss.relstateerr(:, i));

scindices = [pneg cneg env pss];
% if envrel >= t
%     scindices(1) = 1e6;
% elseif pssrel >= t
%     scindices(2) = 1e6;
% end
[minscindices, idxminscindices] = min(scindices);

% if ir <  t && minscindices < t
if minscindices < t
    if idxminscindices == 1
        config(i) = "pneg";
    elseif idxminscindices == 2
        config(i) = "cneg";
    elseif idxminscindices == 3
        % config(i) = "env";
        if envrel >= t
            config(i) = "env_bad_approx";
        else
            config(i) = "env";
        end
    elseif idxminscindices == 4
        % config(i) = "pss";
        if pssrel >= t
            config(i) = "pss_bad_approx";
        else
            config(i) = "pss";
        end
    end
end

% simplification iteration 1
if config(i) == "env_bad_approx"
    scindices = [pneg cneg pss];
    [minscindices, idxminscindices] = min(scindices);
    if minscindices < t
        if idxminscindices == 1
            config(i) = "pneg";
        elseif idxminscindices == 2
            config(i) = "cneg";
        elseif idxminscindices == 3
            if pssrel >= t
                config(i) = "pss_bad_approx";
            else
                config(i) = "pss";
            end
        end
    else
        config(i) = "unclassified";
    end
end

if config(i) == "pss_bad_approx"
    scindices = [pneg cneg env];
    [minscindices, idxminscindices] = min(scindices);
    if minscindices < t
        if idxminscindices == 1
            config(i) = "pneg";
        elseif idxminscindices == 2
            config(i) = "cneg";
        elseif idxminscindices == 3
            if envrel >= t
                config(i) = "env_bad_approx";
            else
                config(i) = "env";
            end
        end
    else
        config(i) = "unclassified";
    end
end

% simplification iteration 2
if config(i) == "env_bad_approx" || config(i) == "pss_bad_approx"
    scindices = [pneg cneg];
    [minscindices, idxminscindices] = min(scindices);
    if minscindices < t
        if idxminscindices == 1
            config(i) = "pneg";
        elseif idxminscindices == 2
            config(i) = "cneg";
        end
    else
        config(i) = "unclassified";
    end
end

end

% correct reduction
% config(model.I.PS) = "pneg";
% config(model.I.PC) = "pneg";
% config(model.I.APC) = "pneg";
% config(model.I.Tmod) = "pneg";

% calculate error
% multiple.multiple = 0;
% [err_index, ~, tred, Xred] = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, [], model.param, multiple, model.odefun, model.jacfun, config, "MRSE");

% show results
% disp(err_index.errout)
disp(sum(config == "dyn"))
model.I.nmstate(config == "dyn")
disp(sum(config == "env"))
model.I.nmstate(config == "env")
disp(sum(config == "env_bad_approx"))
disp(sum(config == "pneg"))
disp(sum(config == "cneg"))
model.I.nmstate(config == "cneg")
disp(sum(config == "pss"))
model.I.nmstate(config == "pss")
disp(sum(config == "pss_bad_approx"))
disp(sum(config == "unclassified"))
model.I.nmstate(config == "unclassified")

% pneg indices of env and pss classified states
max(model.pneg.nindex(:, config == "env"))
max(model.pneg.nindex(:, config == "pss"))

% number of states where approaches are appropriate
sum(max(model.pneg.nindex) < 0.1)
sum(max(model.cneg.nindex) < 0.1)
sum(max(model.env.nindex) < 0.1)
sum(max(model.pss.nindex) < 0.1)
sum(max(model.env.nindex) < 0.1 & max(model.env.relstateerr) < 0.1)
model.I.nmstate(max(model.env.nindex) < 0.1)
model.I.nmstate(max(model.env.nindex) < 0.1 & max(model.env.relstateerr) < 0.1)
sum(max(model.pss.nindex) < 0.1 & max(model.pss.relstateerr) < 0.1)
model.I.nmstate(max(model.pss.nindex) < 0.1 & max(model.pss.relstateerr) < 0.1)

i = model.I.TFPI;
i = model.I.X;
i = model.I.D;
a = max(model.ir.nindex(:, i))
b = max(model.env.nindex(:, i))
c = max(model.pss.nindex(:, i))
d = max(model.cneg.nindex(:, i))
e = max(model.pneg.nindex(:, i))

ir = max(model.ir.index);
sum(ir > 0)
sum(ir > 1e-10)

ir = max(model.ir.nindex, [], 1);
env = max(model.env.nindex, [], 1);
pss = max(model.pss.nindex, [], 1);
cneg = max(model.cneg.nindex, [], 1);
pneg = max(model.pneg.nindex, [], 1);

envrel = max(model.env.relstateerr, [], 1);
pssrel = max(model.pss.relstateerr, [], 1);

pcneg = config == "pneg" | config == "cneg";
diffcpneg = cneg(pcneg) - pneg(pcneg);
sum(isnan(diffcpneg) | diffcpneg > 1e-4)

%% plot index-reduced model

size = 12;
lw = 1;
lwt = 0.5;

model.I.nmstatelegend = cellfun(@(x) strrep(x, '_', ':'), model.I.nmstatelegend, 'UniformOutput', false);

% reference solution
figure
hold on
for i = 1:7
    semilogy(model.t_ref, Xred(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    semilogy(model.t_ref, Xred(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
xlim([-2 42])
ylim([1e-7 5e4])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'southeast')
xlabel("t [h]")
ylabel("concentration [nmol/L]")
yscale('log')
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_index_reduced_005_24_ref_sol.pdf")

% reference solution 0 - 0.15 h
figure
hold on
for i = 1:7
    semilogy(model.t_ref, Xred(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    semilogy(model.t_ref, Xred(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
xlim([-0.005 0.155])
ylim([1e-7 5e4])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'southeast')
xlabel("t [h]")
ylabel("concentration [nmol/L]")
yscale('log')
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_index_reduced_005_24_ref_sol_015h.pdf")

%% code pss states into index-reduced model

% obtain pss states
I_red = config2I(model.I, config, []);

% create symbolic variables
syms t;
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));

% make odefun symbolic
tic
odefun_symbolic = model.odefun(X_sym,par_sym);
disp(toc)

nm_X_pss = model.I.nmstate(config == "pss");
X_sym_pss = X_sym(config == "pss");
odefun_symbolic_pss = odefun_symbolic(config == "pss");

% solve pss states
G = solve(odefun_symbolic_pss, X_sym_pss);

% input solved states to ODEs
for k = 1:length(X_sym_pss)
    solutions = G.(nm_X_pss{k});
    odefun_symbolic_solved = subs(odefun_symbolic, X_sym_pss(k), solutions(1));
end

% remove solved state ODEs
odefun_symbolic_solved(config == "pss") = 0;

% simplify solved ODEs
odefun_symbolic_solved = simplify(odefun_symbolic_solved);

% differentiate jac with solved ODEs
tic
jacfun_symbolic_solved = sym(zeros(model.I.nstates, model.I.nstates));
for i = 1:model.I.nstates
    for j = 1:model.I.nstates
        jacfun_symbolic_solved(i, j) = diff(odefun_symbolic(i), X_sym(j));
    end
end
disp(toc)

% convert to matlabfun and insert to model
model.odefun = matlabFunction(odefun_symbolic_solved,'Vars',{X_sym,par_sym});
model.ode = model.odefun;
model.jacfun = matlabFunction(jacfun_symbolic_solved,'Vars',{X_sym,par_sym,t});
model.jac = model.jacfun;

% update I
model.I.pss_solved = model.I.pss;
model.I.pss = [];

%% check pss solved index reduced model

config(config == "pss") = "pneg";
model.I = config2I(model.I, config, []);
[err_index_solved, ~, tred_solved, Xred_solved] = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, [], model.param, multiple, model.odefun, model.jacfun, config, "MRSE");

%% ir-indices of pss solved index reduced model

% model.t_ref = model.t_ref([1:80:end]);
% model.X_ref = model.X_ref([1:80:end], :);

config = repmat("dyn", [model.I.nstates, 1]);
model.I = config2I(model.I, config, []);

[ir, contr, obs] = compute_ir_indices_matlabfun(model);

save("indices_irred_0.05_BCSV40_from_JKn_pss_solved.mat", "ir", "contr", "obs")

%% pneg run of index-reduced model
mor_options.err_out         = err_index.errout + 0.01;
% mor_options.err_out         = 0.05;
mor_options.err_int         = 10;
mor_options.timeout         = 120;
mor_options.criterion       = 'linear';
mor_options.errtype         = 'MRSE';

model.L = [];

[model.pnegobj, model.pneg_run, model.pnegconfig, ~] = mor_exh_pneg_greedy(model, config, mor_options, false);

disp(sum(model.pnegconfig == "dyn"))
disp(sum(model.pnegconfig == "env"))
disp(sum(model.pnegconfig == "pss"))
disp(sum(model.pnegconfig == "pneg"))
disp(sum(model.pnegconfig == "cneg"))

%% further optimize index-reduced model

model.L = [];

modelname = 'BCSV40_from_JKn';

mode = 'from_start';
eout = 0.05;
eint = 10;
pint = 100;
timeout = 120;
crit = 'linear';
errtype = 'MRSE';
firstrun = 0;
pnegrun = 1;
conlawrun = 0;
% classifs = ["dyn" "cneg" "pneg" "env" "pss"];
classifs = ["dyn" "cneg" "pneg" "env"];

variability = 0;
var_obj_prctile = 90;
LHS_EOG = 0;
backwards = 0;
variability_input = 0;
virtual_pop = [];
X_ref_var = [];
redmodel = [];

log_required = 0;

morexh_model(model, modelname, mode, eout, eint, pint, timeout, crit, errtype, firstrun, pnegrun, conlawrun, classifs, variability, var_obj_prctile, LHS_EOG, variability_input, virtual_pop, X_ref_var, backwards, redmodel, config, log_required)

%% plot optimized model
size = 12;
lw = 1;
lwt = 0.5;

load("modelBC_SV40_from_JKn_2024.mat")

% reference solution
figure
grid on
semilogy(model.t_ref, model.X_ref(:, model.analysis.ir.I_sorted_max_nindex_above_threshold), 'LineWidth', lw) %DisplayName', plotnames(i))
xlim([-2 42])
ylim([1e-7 5e4])
legend(model.analysis.ir.nmstates_above_nindex_threshold, 'Location', 'southeast')
xlabel("t [h]")
ylabel("concentration [g/L]")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_ref_sol.pdf")

load("BCSV40_from_JKn_exh_t120_MRSE_pnegrun_0.05_10_linear_dyncnegpnegenvpss.mat")
% I_red = config2I(redmodel.I, redmodel.redobj.redconfig, []);
[~, X_red, ~, ~] = simModel(redmodel.t_ref, redmodel.X0, redmodel.par, redmodel.I, redmodel.param, redmodel.multiple, redmodel.odefun, redmodel.jacfun);

% reference solution
figure
grid on
semilogy(redmodel.t_ref, X_red(:, redmodel.analysis.ir.I_sorted_max_nindex_above_threshold), 'LineWidth', lw) %DisplayName', plotnames(i))
xlim([-2 42])
ylim([1e-7 5e4])
legend(model.analysis.ir.nmstates_above_nindex_threshold, 'Location', 'southeast')
xlabel("t [h]")
ylabel("concentration [g/L]")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_redmodel_sol.pdf")

%% indices for ir-reduced and optimized 8-state model

% load model
load("BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv.mat")
redmodel.I = config2I(redmodel.I, redmodel.pneg_run.configs(end, :), []);

[ir, contr, obs] = compute_ir_indices_matlabfun(redmodel);

save("indices_irred_0.05_BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv.mat", "ir", "contr", "obs")

%% plot
size = 12;
lw = 1;
lwt = 0.5;

% reference solution
figure
semilogy(redmodel.t_ref, redmodel.X_red(:, redmodel.analysis.ir.I_sorted_max_nindex_above_threshold), 'LineWidth', lw) %DisplayName', plotnames(i))
xlim([-2 42])
ylim([1e-7 5e4])
legend(redmodel.analysis.ir.nmstates_above_nindex_threshold, 'Location', 'southeast')
xlabel("t [h]")
ylabel("concentration [g/L]")
% grid on

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_indices_irred_0.05_BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv_ref_sol.pdf")

% nir indices
figure
plot(redmodel.t_ref, ir.nindex(:, redmodel.analysis.ir.I_sorted_max_nindex_above_threshold), 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
% xlim([-0.14 4.15])
xlim([-0.3 20])
ylim([-0.01 1])
legend([redmodel.analysis.ir.nmstates_above_nindex_threshold; 'threshold'], 'Location','east')
xlabel("t [h]")
ylabel("nir-index")
% grid on

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_indices_irred_0.05_BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv_ir.pdf")

%% analyze env-pruned relstateerr

load("modelBC_SV40_from_JKn_2024.mat")

max(model.env.nindex(:, model.analysis.ir.I_sorted_max_nindex))
max(model.pss.nindex(:, model.analysis.ir.I_sorted_max_nindex))
max(model.pneg.nindex(:, model.analysis.ir.I_sorted_max_nindex))
max(model.cneg.nindex(:, model.analysis.ir.I_sorted_max_nindex))

max(model.env.relstateerr(:, model.analysis.ir.I_sorted_max_nindex))
max(model.pss.relstateerr(:, model.analysis.ir.I_sorted_max_nindex))

env_relstateerr_pruned = model.env.relstateerr;
env_relstateerr_pruned(model.env.nindex < 0.01) = 0;
env_relstateerr_pruned = env_relstateerr_pruned(:, model.analysis.ir.I_sorted_max_nindex);
max(env_relstateerr_pruned);

%% indices for ir-reduced model

% calculate reduced model ir indices
model.I = config2I(model.I, config, []);
[irred.ir, irred.contr, irred.obs, irred.t_ir] = compute_ir_indices_matlabfun(model);

% save indices
save("results_irred_BC_40h_JKn_t005.mat", "irred")

% nir-indices
size = 12;
lw = 1;
lwt = 0.5;

figure
plot(irred.t_ir, irred.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold), 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
% xlim([-0.14 4.15])
xlim([-0.01 1])
ylim([-0.01 1])
legend([model.analysis.ir.nmstates_above_nindex_threshold; 'threshold'], 'Location','east')
xlabel("t [h]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_JKn_ir_index_red_0_1.pdf")


%% Blood Coagulation 40h: lumping Aarons
load("modelBC_SV40_from_JKn_2024.mat")

[lump_matrices,inv_lump_matrices,errors,out_states] = lumping_Aarons(model);

save("lumpingres_BC_SV40_from_JKn.mat", "out_states", "errors", "inv_lump_matrices", "lump_matrices")

%% analyze lumping
idx = 1:length(errors);
idx = idx(errors < 0.05);
idx = idx(end);

lumpmat = lump_matrices{idx};

for k = 1:size(lumpmat, 1)
    if sum(lumpmat(k, :)) == 1
        disp(k)
        disp(model.I.nmstate{logical(lumpmat(k, :))})
    end
end

disp(model.I.nmstate{logical(lumpmat(out_states(idx), :))})
% disp(model.I.nmstate{logical(lumpmat(out_states(idx), :))})

disp(errors(idx))

%% save reduced solution plot
size = 12;
lw = 1;
lwt = 0.5;

figure
grid on
semilogy(tred, Xred(:, model.analysis.ir.I_sorted_max_nindex_above_threshold), 'LineWidth', lw) %DisplayName', plotnames(i))
xlim([-2 42])
ylim([1e-7 5e4])
legend(model.analysis.ir.nmstates_above_nindex_threshold, 'Location','southeast')
xlabel("t [h]")
ylabel("concentration [g/L]")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_ref_sol_indices.pdf")

%% Full model relevant states

X_full = model.X_ref(:, [I.Fg, I.II, I.IIa, I.AVenom, I.CVenom, I.Xa, I.Xa_Va, I.Tmod, I.AT_III_Heparin, I.TaipanVenom, I.CVenom_Tiger]);

%% Lumping as in Gulati 2014
% load model file
% load("modelBC_Gulati2014_in_vivo_full.mat")
% load("modelBCSV_minimal.mat")
load("modelBC_temp.mat")
model.X0 = Gulati2014BloodCoagulation_initialvalues(model);
model.multiple.multiple = 0;
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);

idxFg = model.I.Fg;
idxII = model.I.II;
idxIIa = model.I.IIa;
idxAvenom = model.I.AVenom;
idxPvenom = model.I.CVenom;

lumpmat_Gulati = zeros(5, model.I.nstates);
lumpmat_Gulati(5, :) = 1;
lumpmat_Gulati(1, idxFg) = 1;
lumpmat_Gulati(5, idxFg) = 0;
lumpmat_Gulati(2, idxIIa) = 1;
lumpmat_Gulati(5, idxIIa) = 0;
lumpmat_Gulati(3, idxAvenom) = 1;
lumpmat_Gulati(5, idxAvenom) = 0;
lumpmat_Gulati(4, idxPvenom) = 1;
lumpmat_Gulati(5, idxPvenom) = 0;

opt.lumpmat = lumpmat_Gulati;

% lumpmat_Gulati = zeros(6, model.I.nstates);
% lumpmat_Gulati(5, :) = 1;
% lumpmat_Gulati(1, idxFg) = 1;
% lumpmat_Gulati(5, idxFg) = 0;
% lumpmat_Gulati(6, idxII) = 1;
% lumpmat_Gulati(5, idxII) = 0;
% lumpmat_Gulati(2, idxIIa) = 1;
% lumpmat_Gulati(5, idxIIa) = 0;
% lumpmat_Gulati(3, idxAvenom) = 1;
% lumpmat_Gulati(5, idxAvenom) = 0;
% lumpmat_Gulati(4, idxPvenom) = 1;
% lumpmat_Gulati(5, idxPvenom) = 0;

[reduced_errors.lumping_Gulati, X_red] = calculate_lumping_error(model, opt);

AUC = trapz(model.t_ref, model.X_ref, 1);
opt.invlumpmat = opt.lumpmat';
for col = 1:size(opt.invlumpmat, 2)
    opt.invlumpmat(:, col) = (AUC.*opt.invlumpmat(:, col)') ./ (AUC*opt.invlumpmat(:, col));
end
[reduced_errors.scaled_lumping_Gulati, X_red] = calculate_lumping_error(model, opt);

% Save results
save("./results/errors_BC_Gulati2014.mat", "reduced_errors")

%% Diagnostics
I = model.I;
par = model.par;

disp(par(I.degIIa))

disp(par(I.v14))
disp(par(I.k14))

disp(par(I.v15))
disp(par(I.k15))

% disp("v15")
% disp(par(I.v15))

disp("v12")
disp(par(I.v12))
disp(par(I.k12))

disp("v13")
disp(par(I.v13))
disp(par(I.k13))

disp("v14")
disp(par(I.v14))
disp(par(I.k14))

% disp(par(I.pFg)/par(I.degFg))

disp("Brown")
disp(par(I.ka_Brown))
disp(par(I.d_Brown))

disp("Fg")
disp(par(I.pFg))
disp(par(I.degFg))

disp(model.X0(I.Fg))

%% Plots
load("modelBC_temp.mat")

t = tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
nexttile(t)
plot(model.t_ref, model.X_ref(:, [model.I.AVenom model.I.CVenom model.I.II model.I.IIa model.I.Fg]), 'LineWidth', 2)
ylim([1e-4 1e4])
grid on
set(gca, 'YScale', 'log')
legend('AVenom', 'CVenom', 'II', 'IIa', 'Fibrinogen', 'Location', 'southeast')
title('Full model')

nexttile(t)
% plot(model.t_ref, X_red(:, [3 4 5 2 1]).*repmat([1 1 1/58 1 1], [length(model.t_ref), 1]), 'LineWidth', 2)
% plot(t_ref, X_red(:, [3 4 5 2 1]).*repmat([1 1 1/58 1 1], [length(t_ref), 1]), 'LineWidth', 2)
plot(t_red, X_red, 'LineWidth', 2)
ylim([1e-4 1e4])
grid on
set(gca, 'YScale', 'log')
legend('AVenom', 'CVenom', 'lumped', 'IIa', 'Fibrinogen', 'Location', 'southeast')
% title('Lumped model from our implementation')
title('Lumped model as reported in Gulati2014')

nexttile(t)
plot(t_red, (X_red(:, 5) - X_ref(:, model.I.Fg))./X_ref(:, model.I.Fg), 'LineWidth', 2)
% plot(t_ref, (X_red(:, 1) - X_ref(:, model.I.Fg))./X_ref(:, model.I.Fg), 'LineWidth', 2)
% ylim([1e-4 1e4])
grid on
% set(gca, 'YScale', 'log')
legend('relative error')
title('Relative error')

exportgraphics(t, "./figures/BClumped_reported_Gulati_ODEs.pdf")
% exportgraphics(t, "./figures/BClumped_our_implementation.pdf")

%% Get ODE form
% Uni Potsdam implementation
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));

ODE_red = simplify(lumpmat_Gulati*model.odefun(X_sym, par_sym));

%% Get ODE form
% Gulati implementation
x = sym('x', [62 1]);
syms t
syms pX pV pII pVIII pIX pFg pXIII pPg pPC pTmod pXI pVII pTFPI pPS pVK pXII pPK dX dXa dV dVa dVaXa dII dIIa dVIII dVIIIa dIXa dIXaVIIIa dIX dXIa dFg dFDP dF dXF dD dXIII dXIIIa dPg dP dPC dAPC dTmod dIIaTmod dTaipan dXI dXIIa dVII dVIIa dTAT dBrown dTF dVIITF dVIIaTF dTFPI dXaTFPI dVIIaTFXaTFPI dPS dAPCPS dHeparin dVK_1 dVK_2 dVK_3 dVKH2 dVKO dXII dPK dK dca vIXa kIXa vII2IIaXa kII2IIaXa vV2VaIIa kV2VaIIa vVIII2VIIIaIIa kVIII2VIIIaIIa vIXaVIIIa kIXaVIIIa vXIa kXIa vVaXa kVaXa vFg2F kFg2F vFg2FDP kFg2FDP vF2FDP kF2FDP vF2XF kF2XF vXF2D kXF2D vXIII2XIIIa kXIII2XIIIa vPg2PF kPg2PF vPg2PIIa kPg2PIIa vPC2APC kPC2APC vVIIIaLoss kVIIIaLoss vVaLoss kVaLoss vXaVaLoss kXaVaLoss vXF2DAPC kXF2DAPC vPg2PAPC kPg2PAPC vXI2XIaIIa kXI2XIaIIa vXI2XIaXIIa kXI2XIaXIIa vX2XaVIIa kX2XaVIIa vVII2VIIaIIa kVII2VIIaIIa vVII2VIIaTaipan kVII2VIIaTaipan vVIITF_Xa kVIITF_Xa vX_VIIaTF kX_VIIaTF vIX_VIIaTF kIX_VIIaTF vVIITF_TF kVIITF_TF vVII_Xa kVII_Xa vVII_VIIaTF kVII_VIIaTF vVII_IXa kVII_IXa vXII_ca kXII_ca vXII_K kXII_K vPK_XIIa kPK_XIIa cVaXa cIXaVIIIa cIIaTmod cVIIaTF cVIITF cVIIaTFXaTFPI cXaTFPI cAPCPS cXaHeparin cIIaHeparin cIXaHeparin cIXaUFH cXaUFH cIIaUFH Imax IC50 ka v ke lagt ka_e k10_e k12_e k21_e v2_e k12_vk k21_vk vc_vk ke_UFH R_UFH T_UFH time_of_antivenom kr FLAG1 FLAG2 KDXa KDE ka_brown ka_tiger dTiger

lumpmat_Gulati_supp = zeros(5, 62);
lumpmat_Gulati_supp(5, :) = 1;
lumpmat_Gulati_supp(1, 14) = 1;              % Fg
lumpmat_Gulati_supp(5, 14) = 0;
lumpmat_Gulati_supp(2, 7) = 1;               % IIa
lumpmat_Gulati_supp(5, 7) = 0;
lumpmat_Gulati_supp(3, 28) = 1;              % AVenom
lumpmat_Gulati_supp(5, 28) = 0;
lumpmat_Gulati_supp(4, 62) = 1;              % CVenom
lumpmat_Gulati_supp(5, 62) = 0;

ODE_red_Gulati = simplify(lumpmat_Gulati_supp * Gulati_Supplement_ode(t,x,pX,pV,pII,pVIII,pIX,pFg,pXIII,pPg,pPC,pTmod,pXI,pVII,pTFPI,pPS,pVK,pXII,pPK,dX,dXa,dV,dVa,dVaXa,dII,dIIa,dVIII,dVIIIa,dIXa,dIXaVIIIa,dIX,dXIa,dFg,dFDP,dF,dXF,dD,dXIII,dXIIIa,dPg,dP,dPC,dAPC,dTmod,dIIaTmod,dTaipan,dXI,dXIIa,dVII,dVIIa,dTAT,dBrown,dTF,dVIITF,dVIIaTF,dTFPI,dXaTFPI,dVIIaTFXaTFPI,dPS,dAPCPS,dHeparin,dVK_1,dVK_2,dVK_3,dVKH2,dVKO,dXII,dPK,dK,dca,vIXa,kIXa,vII2IIaXa,kII2IIaXa,vV2VaIIa,kV2VaIIa,vVIII2VIIIaIIa,kVIII2VIIIaIIa,vIXaVIIIa,kIXaVIIIa,vXIa,kXIa,vVaXa,kVaXa,vFg2F,kFg2F,vFg2FDP,kFg2FDP,vF2FDP,kF2FDP,vF2XF,kF2XF,vXF2D,kXF2D,vXIII2XIIIa,kXIII2XIIIa,vPg2PF,kPg2PF,vPg2PIIa,kPg2PIIa,vPC2APC,kPC2APC,vVIIIaLoss,kVIIIaLoss,vVaLoss,kVaLoss,vXaVaLoss,kXaVaLoss,vXF2DAPC,kXF2DAPC,vPg2PAPC,kPg2PAPC,vXI2XIaIIa,kXI2XIaIIa, vXI2XIaXIIa,kXI2XIaXIIa,vX2XaVIIa,kX2XaVIIa,vVII2VIIaIIa,kVII2VIIaIIa,vVII2VIIaTaipan,kVII2VIIaTaipan,vVIITF_Xa,kVIITF_Xa,vX_VIIaTF,kX_VIIaTF,vIX_VIIaTF,kIX_VIIaTF,vVIITF_TF,kVIITF_TF,vVII_Xa,kVII_Xa,vVII_VIIaTF,kVII_VIIaTF,vVII_IXa,kVII_IXa,vXII_ca,kXII_ca,vXII_K,kXII_K,vPK_XIIa,kPK_XIIa,cVaXa,cIXaVIIIa,cIIaTmod,cVIIaTF,cVIITF,cVIIaTFXaTFPI,cXaTFPI,cAPCPS,cXaHeparin,cIIaHeparin,cIXaHeparin,cIXaUFH,cXaUFH,cIIaUFH,Imax,IC50,ka,v,ke,lagt,ka_e,k10_e,k12_e,k21_e,v2_e,k12_vk,k21_vk,vc_vk,ke_UFH,R_UFH,T_UFH,time_of_antivenom,kr,FLAG1,FLAG2,KDXa,KDE,ka_brown,ka_tiger,dTiger));

%% modify Wajima2009 equations to match Gulati2014 paper ODE diagram

load("modelBC_temp.mat")
% turn off P deactivation of Fg
% model.par(model.I.v15) = 0;
% model.par(model.I.v13) = 0;
% model.par(model.I.v12) = 1940;
% model.par(model.I.k12) = 1.38;
% load modified odefun with autoactivating compounds to IIa turned off
% model.odefun = @(t,X,par,model) Wajima2009BloodCoagulation_ode_Gulati_compatibility(t,X,par,model);
% convert odefun to (t,X) matlabfun
% syms t
% par_sym = cell2sym(model.I.nmpar(:));
% X_sym = cell2sym(model.I.nmstate(:));
% odefun_symbolic = model.odefun(t,X_sym,par_sym,model);
% model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
% make lumping matrix
I = model.I;
lumpmat_Gulati = zeros(5, model.I.nstates);
lumpmat_Gulati(5, :) = 1;
lumpmat_Gulati(1, I.Fg) = 1;
lumpmat_Gulati(5, I.Fg) = 0;
lumpmat_Gulati(2, I.IIa) = 1;
lumpmat_Gulati(5, I.IIa) = 0;
lumpmat_Gulati(3, I.AVenom) = 1;
lumpmat_Gulati(5, I.AVenom) = 0;
lumpmat_Gulati(4, I.CVenom) = 1;
lumpmat_Gulati(5, I.CVenom) = 0;
% check lumping error
[reduced_errors.lumping_Gulati, X_red] = calculate_lumping_error(model, lumpmat_Gulati);

load("modelBC_temp.mat")
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
[t_ref, X_ref] = simModel(model.t_ref, model.X0, model.par, model.I, [], model.multiple, model.odefun, []);
relativeErrorL2(t_ref, X_ref(:, model.I.Fg), X_red(:, 5))

%% explicitly written lumped Gulati2014 model - fitted parameters

% make explicitly lumped model
model = struct;
model.I = Gulati2014BClumped_indexing();
model.X0 = Gulati2014BClumped_initialvalues_fitted(model);
model.par = Gulati2014BClumped_parameters_fitted(model);
model.odefun = @(t,X,par,model) Gulati2014BClumped_ode(t,X,par,model);
syms t
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));
odefun_symbolic = model.odefun(t,X_sym,par_sym,model);
model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
model.multiple.multiple = 0;
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
% simulate
[t_red, X_red] = simModel([0 40], model.X0, model.par, model.I, [], model.multiple, model.odefun, []);

%% explicitly written lumped Gulati2014 model - original parameters (except for lumped state)

% make explicitly lumped model
model = struct;
model.I = Gulati2014BClumped_indexing();
model.X0 = Gulati2014BClumped_initialvalues_original(model);
model.par = Gulati2014BClumped_parameters_original(model);
model.odefun = @(t,X,par,model) Gulati2014BClumped_ode(t,X,par,model);
syms t
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));
odefun_symbolic = model.odefun(t,X_sym,par_sym,model);
model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
model.multiple.multiple = 0;
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
% simulate
[t_red, X_red] = simModel([0 40], model.X0, model.par, model.I, [], model.multiple, model.odefun, []);

load("modelBC_temp.mat")
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
[t_ref, X_ref] = simModel(t_red, model.X0, model.par, model.I, [], model.multiple, model.odefun, []);
relativeErrorL2(t_ref, X_ref(:, model.I.Fg), X_red(:, 5))

%% ?????????????

% make explicitly lumped model
model = struct;
model.I = Gulati2014BClumped_indexing();
model.X0 = Gulati2014BClumped_initialvalues_original(model);
model.par = Gulati2014BClumped_parameters_original(model);
model.odefun = @(t,X,par,model) Gulati2014BClumped_ode(t,X,par,model);
syms t
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));
odefun_symbolic = model.odefun(t,X_sym,par_sym,model);
model.odefun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
model.multiple.multiple = 0;
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
% simulate
[t_red, X_red] = simModel([0 40], model.X0, model.par, model.I, [], model.multiple, model.odefun, []);

load("modelBC_temp.mat")
model.I = config2I(model.I, repmat("dyn", [1 model.I.nstates]), []);
[t_ref, X_ref] = simModel(t_red, model.X0, model.par, model.I, [], model.multiple, model.odefun, []);
relativeErrorL2(t_ref, X_ref(:, model.I.Fg), X_red(:, 5))

%% Error fun
function Error=relativeErrorL2(t,X,Y)
    Error=sqrt(trapz(t,(X-Y).^2)/trapz(t,X.^2));
end

%% investigate TaipanVenom pss

load("modelBC_SV40_from_JKn_2024.mat")
% model.L = [];

ir = max(model.ir.index);

config = repmat("dyn", [1 model.I.nstates]);
config(ir == 0) = "pneg";
% config(model.I.TaipanVenom) = "pss";
% model.I.TaipanVenom

[~, Xout] = simModel_simple(model, config);

model.relerrnorm(model.t_ref, Xout(:, model.I.output), model.X_ref(:, model.I.output))

objfun_simple(model, config, "rel2NE")

f = figure;
hold on
plot(model.t_ref, model.X_ref(:, model.I.output))
plot(model.t_ref, Xout(:, model.I.output))

%% reduce via ir-indices of 1h scenario for 40h model (using pneg, env) - derive 8-state model

model_BC_SV1 = load("modelBC_SV1_from_JKn_2024_irmatlabfun.mat").model;
model = load("modelBC_SV40_from_JKn_2024.mat").model;

redmodel = mor_sequential_JKn_2018(model, flip(model_BC_SV1.analysis.ir.I_sorted_max_nindex), 0.1);

save("modelBC_SV40_from_JKn_2024_reduced_8_state.mat", "redmodel")

%% ir-indices of 8-state model

% load model
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")

[ir, contr, obs, t_ind] = compute_ir_indices_matlabfun(redmodel);

save("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat", "ir", "contr", "obs", "t_ind")

%% reduce via ir-indices of 40h scenario for 40h model (using pneg, env) - derive 13-state model

model = load("modelBC_SV40_from_JKn_2024.mat").model;

redmodel = mor_sequential_JKn_2018(model, flip(model.analysis.ir.I_sorted_max_nindex), 0.1);

save("modelBC_SV40_from_JKn_2024_reduced_13_state.mat", "redmodel")

%% ir-indices of 13-state model

% load model
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")

[ir, contr, obs, t_ind] = compute_ir_indices_matlabfun(redmodel);

save("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat", "ir", "contr", "obs", "t_ind")

%% make small Venom-II-Fg model

redmodel = load("modelBC_SV40_from_JKn_2024.mat").model;

redmodel.redconfig = repmat("pneg", [1 redmodel.I.nstates]);
redmodel.redconfig(redmodel.I.AVenom) = "dyn";
redmodel.redconfig(redmodel.I.CVenom) = "dyn";
redmodel.redconfig(redmodel.I.Fg) = "dyn";
redmodel.redconfig(redmodel.I.II) = "env";
redmodel.redconfig(redmodel.I.IIa) = "dyn";
redmodel.redconfig(redmodel.I.Pg) = "env";
redmodel.redconfig(redmodel.I.P) = "dyn";
[redmodel.redobj, ~, redmodel.t_red, redmodel.X_red, ~] = objfun_simple(redmodel, redmodel.redconfig, "rel2NE");

save("modelBC_SV40_from_JKn_2024_reduced_small.mat", "redmodel")

hold on
plot(redmodel.t_ref, redmodel.X_ref(:, redmodel.I.Fg))
plot(redmodel.t_ref, redmodel.X_red(:, redmodel.I.Fg))
hold off

%% reduce model for n=1000 virtual population

model_BC_SV1 = load("modelBC_SV1_from_JKn_2024_irmatlabfun.mat").model;
model = load("modelBC_SV40_from_JKn_2024.mat").model;
load("modelBC_SV40_population_CV40.mat")
% variability = load("modelBC_population_from_UFa_2023.mat").variability;

redmodel_seq1 = mor_sequential_JKn_2018(model, flip(model_BC_SV1.analysis.ir.I_sorted_max_nindex), 0.1);

redmodel_back = mor_backwards_UFa_2023(redmodel_seq1, flip(model.analysis.ir.I_sorted_max_nindex), 0.1, variability, 95);

redmodel_seq_var = mor_sequential_JKn_2018(redmodel_back, flip(model.analysis.ir.I_sorted_max_nindex), 0.1, variability, 95);

% UFa 40% CV/90% of pop: 25/4
% JTi 40% CV/95% of pop: 
% JTi 40% CV/90% of pop: 24/3
% JTi 20% CV/90% of pop: 24/3

save("modelBC_SV40_redvar_CV40_popprct95.mat", "redmodel_seq_var")

%% check 25-state reduced model for UFa 2023 virtual population

% 8-state model
% load("modelBC_population_from_UFa_2023")
load("modelBC_SV40_population_CV40.mat")
load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")

n_8_state = 0;
tic
for i = 2:1001
    model.X0 = variability.X0_pop(i, :)';
    model.par = variability.par_pop(i, :)';
    model.X_ref = variability.X_ref_pop{i};
    [redobj, ~, tred, Xred, err] = objfun_simple(model, redmodel.redconfig, "rel2NE");
    if redobj.errout <= 0.1
        n_8_state = n_8_state+1;
    end
end
disp(n_8_state)
toc

% 25-state model
load("modelBC_SV40_from_JKn_2024.mat")
load("modelBC_SV40_redvar_CV40_popprct95.mat")

% I = model.I;
% redconfig = repmat("pneg", [1 model.I.nstates]);
% % redconfig([I.AVenom, I.CVenom, I.Fg, I.F, I.II, I.IIa, I.V, I.Va, I.VIII, I.VIIIa, I.IX, I.IXa, I.IXa_VIIIa, I.Xa, I.Xa_Va, ...
% % I.XI, I.XIa, I.XIII, I.XIIIa, I.Tmod, I.IIa_Tmod, I.APC, I.APC_PS, I.Pg, I.P]) = "dyn";
% redconfig([I.AVenom, I.CVenom, I.Fg, I.II, I.IIa, I.V, I.Va, I.VIII, I.VIIIa, I.IX, I.IXa, I.IXa_VIIIa, I.X, I.Xa, I.Xa_Va, ...
% I.XI, I.XIa, I.Tmod, I.IIa_Tmod, I.APC, I.APC_PS, I.Pg, I.P]) = "dyn";
% % redconfig([I.X, I.PS, I.PC, I.TFPI, I.VKH2]) = "env";
% redconfig([I.PS, I.TFPI, I.VKH2]) = "env";
redconfig = redmodel_seq_var.redconfig;

n_25_state = 0;
errors = [];
tic
for i = 2:1001
    model.X0 = variability.X0_pop(i, :)';
    model.par = variability.par_pop(i, :)';
    model.X_ref = variability.X_ref_pop{i};
    [redobj, ~, tred, Xred, err] = objfun_simple(model, redconfig, "rel2NE");
    errors = [errors, redobj.errout];
    if redobj.errout <= 0.1
        n_25_state = n_25_state+1;
    end
end
disp(n_25_state)
toc

%% make virtual population

% load model
load("modelBC_SV40_from_JKn_2024.mat")
par = model.par;
X0 = model.X0;

% generate virtual pop
seed = 1234;
rng(seed, "twister");
% Npop = 1000;
Npop = 10;

CV = 0.4;

virtual_pop_par = zeros(Npop+1, length(par));
% lognpdf_par = zeros(Npop, length(par));
for i = 1:length(par)
    if par(i) ~= 0
        m = par(i);
        v = (CV*par(i))^2;
        mu = log(m^2 / sqrt(v + m^2));
        sigma = sqrt(log(v / (m^2) + 1));
        % disp(m)
        % disp(sqrt(v))
        disp(sigma)

        virtual_pop_par(2:end, i) = lognrnd(mu, sigma, [Npop 1]);
        % lognpdf_par(:, i) = lognpdf(virtual_pop_par(:, i), mu, sigma);
    else
        % lognpdf_par(:, i) = 1;
    end
end

% steady-state initials and reference solution for virtual pop
X0_no_input = X0;
X0_no_input(model.I.input) = 0;
virtual_pop_par_initials_ss = zeros(Npop+1, length(X0));
X_ref_var_ss = cell([1 Npop+1]);
for npop = 2:(Npop + 1)
    [~, curr_X_ref, ~, ~] = simModel([0 2000], X0_no_input, virtual_pop_par(npop, :)', model.I, model.param, model.multiple, model.odefun, model.jacfun);
    virtual_pop_par_initials_ss(npop, :) = curr_X_ref(end, :);
    virtual_pop_par_initials_ss(npop, model.I.input) = X0(model.I.input);
    virtual_pop_par_initials_ss(npop, (virtual_pop_par_initials_ss(npop, :) < 1e-10)) = 0;
    [~, X_ref_var_ss{npop}, ~, ~] = simModel(model.t_ref, virtual_pop_par_initials_ss(npop, :)', virtual_pop_par(npop, :)', model.I, model.param, model.multiple, model.odefun, model.jacfun);
end

X_ref_var_ss{1} = model.X_ref;
virtual_pop_par_initials_ss(1, :) = X0';
virtual_pop_par(1, :) = par';

variability.X0_pop = virtual_pop_par_initials_ss;
variability.par_pop = virtual_pop_par;
variability.X_ref_pop = X_ref_var_ss;
variability.variability = 1;

% save(['../Core/modelfiles/modelBC_SV40_population_CV' char(num2str(100*CV)) '.mat'], "variability")