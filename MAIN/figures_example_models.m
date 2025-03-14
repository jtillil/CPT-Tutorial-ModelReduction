%% Setup
clear; clc;
addpath(genpath("../../CPT-Tutorial-ModelReduction"))
size = 8;
lw = 1;
lwt = 0.5;

%% Parallel Pathways: Scenario 1 - no crosstalk
load("modelSPP_no_crosstalk_full.mat")

% reference solution
figure
semilogy(model.t_ref, model.X_ref, 'LineWidth', lw) %DisplayName', plotnames(i))
xlim([-0.002 0.052])
ylim([1e-3 1e4])
legend('A', 'S', 'B', 'C', 'D', 'Location','northeast')
xlabel("t [min]")
ylabel("concentration [nM]")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_no_crosstalk_ref_sol.pdf")

% nir-indices
figure
plot(model.t_ref, model.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-0.002 0.052])
ylim([-0.01 1])
legend('A', 'S', 'B', 'C', 'D', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_no_crosstalk_ir.pdf")

% non ir-indices B
figure
hold on
semilogy(model.t_ref, model.env.nindex(:, model.I.B), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pss.nindex(:, model.I.B), 'LineStyle', '--', 'LineWidth', lw)
semilogy(model.t_ref, model.cneg.nindex(:, model.I.B), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pneg.nindex(:, model.I.B), 'LineStyle', '--', 'LineWidth', lw)
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-0.002 0.052])
ylim([5e-4 1e1])
set(gca, 'YScale', 'log')
box on
legend('env', 'pss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_no_crosstalk_non_ir_B.pdf")

% non ir-indices C
figure
hold on
semilogy(model.t_ref, model.env.nindex(:, model.I.C), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pss.nindex(:, model.I.C), 'LineStyle', '--', 'LineWidth', lw)
semilogy(model.t_ref, model.cneg.nindex(:, model.I.C), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pneg.nindex(:, model.I.C), 'LineStyle', '--', 'LineWidth', lw)
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-0.002 0.052])
ylim([5e-4 1e1])
set(gca, 'YScale', 'log')
box on
legend('env', 'pss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
% hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_no_crosstalk_non_ir_C.pdf")

% non ir-indices S
figure
hold on
semilogy(model.t_ref, model.env.nindex(:, model.I.S), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pss.nindex(:, model.I.S), 'LineStyle', '--', 'LineWidth', lw)
semilogy(model.t_ref, model.cneg.nindex(:, model.I.S), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pneg.nindex(:, model.I.S), 'LineStyle', '--', 'LineWidth', lw)
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-0.002 0.052])
ylim([5e-4 1e1])
set(gca, 'YScale', 'log')
box on
legend('env', 'pss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
% hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_no_crosstalk_non_ir_S.pdf")

%% Parallel Pathways: Scenario 2 - with crosstalk
load("modelSPP_with_crosstalk_full.mat")

% reference solution
figure
semilogy(model.t_ref, model.X_ref, 'LineWidth', lw) %DisplayName', plotnames(i))
xlim([-0.002 0.052])
ylim([1e-3 1e4])
legend('A', 'S', 'B', 'C', 'D', 'Location','northeast')
xlabel("t [min]")
ylabel("concentration [nM]")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_with_crosstalk_ref_sol.pdf")

% nir-indices
figure
plot(model.t_ref, model.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-0.002 0.052])
ylim([-0.01 1])
legend('A', 'S', 'B', 'C', 'D', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_with_crosstalk_ir.pdf")

% non ir-indices B
figure
hold on
semilogy(model.t_ref, model.env.nindex(:, model.I.B), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pss.nindex(:, model.I.B), 'LineStyle', '--', 'LineWidth', lw)
semilogy(model.t_ref, model.cneg.nindex(:, model.I.B), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pneg.nindex(:, model.I.B), 'LineStyle', '--', 'LineWidth', lw)
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-0.002 0.052])
ylim([5e-4 1e1])
set(gca, 'YScale', 'log')
box on
legend('env', 'pss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_with_crosstalk_non_ir_B.pdf")

% non ir-indices C
figure
hold on
semilogy(model.t_ref, model.env.nindex(:, model.I.C), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pss.nindex(:, model.I.C), 'LineStyle', '--', 'LineWidth', lw)
semilogy(model.t_ref, model.cneg.nindex(:, model.I.C), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pneg.nindex(:, model.I.C), 'LineStyle', '--', 'LineWidth', lw)
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-0.002 0.052])
ylim([5e-4 1e1])
set(gca, 'YScale', 'log')
box on
legend('env', 'pss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_with_crosstalk_non_ir_C.pdf")


%% Enzyme Kinetics: Scenario 1 - Cpss Eenv
load("modelMMEK_Cpss_Eenv_full.mat")

% reference solution
figure
semilogy(model.t_ref, model.X_ref, 'LineWidth', lw) %DisplayName', plotnames(i))
xlim([-1 31])
ylim([1e-2 5e3])
legend('A', 'S', 'E', 'C', 'P', 'Location','northeast')
xlabel("t [min]")
ylabel("concentration [nM]")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_Eenv_ref_sol.pdf")

% nir-indices
figure
plot(model.t_ref, model.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-1 31])
ylim([-0.01 1])
legend('A', 'S', 'E', 'C', 'P', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_Eenv_ir.pdf")

% non ir-indices E
figure
hold on
semilogy(model.t_ref, model.env.nindex(:, model.I.E), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pss.nindex(:, model.I.E), 'LineStyle', '--', 'LineWidth', lw)
semilogy(model.t_ref, model.cneg.nindex(:, model.I.E), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pneg.nindex(:, model.I.E), 'LineStyle', '--', 'LineWidth', lw)
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-1 31])
ylim([5e-4 1e1])
set(gca, 'YScale', 'log')
box on
legend('env', 'pss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_Eenv_non_ir_E.pdf")

% non ir-indices C
figure
hold on
semilogy(model.t_ref, model.env.nindex(:, model.I.C), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pss.nindex(:, model.I.C), 'LineStyle', '--', 'LineWidth', lw)
semilogy(model.t_ref, model.cneg.nindex(:, model.I.C), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pneg.nindex(:, model.I.C), 'LineStyle', '--', 'LineWidth', lw)
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-1 31])
ylim([5e-4 1e1])
set(gca, 'YScale', 'log')
box on
legend('env', 'pss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_Eenv_non_ir_C.pdf")

%% Enzyme Kinetics: Scenario 2 - Cpss EpCenv
load("modelMMEK_Cpss_EpCenv_full.mat")

% reference solution
figure
semilogy(model.t_ref, model.X_ref, 'LineWidth', lw) %DisplayName', plotnames(i))
xlim([-1 31])
ylim([1e-2 5e3])
legend('A', 'S', 'E', 'C', 'P', 'Location','northeast')
xlabel("t [min]")
ylabel("concentration [nM]")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_EpCenv_ref_sol.pdf")

% nir-indices
figure
plot(model.t_ref, model.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-1 31])
ylim([-0.01 1])
legend('A', 'S', 'E', 'C', 'P', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_EpCenv_ir.pdf")

% non ir-indices E
figure
hold on
semilogy(model.t_ref, model.env.nindex(:, model.I.E), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pss.nindex(:, model.I.E), 'LineStyle', '--', 'LineWidth', lw)
semilogy(model.t_ref, model.cneg.nindex(:, model.I.E), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pneg.nindex(:, model.I.E), 'LineStyle', '--', 'LineWidth', lw)
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-1 31])
ylim([5e-4 1e1])
set(gca, 'YScale', 'log')
box on
legend('env', 'pss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_EpCenv_non_ir_E.pdf")

% non ir-indices C
figure
hold on
semilogy(model.t_ref, model.env.nindex(:, model.I.C), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pss.nindex(:, model.I.C), 'LineStyle', '--', 'LineWidth', lw)
semilogy(model.t_ref, model.cneg.nindex(:, model.I.C), 'LineStyle', '-', 'LineWidth', lw)
semilogy(model.t_ref, model.pneg.nindex(:, model.I.C), 'LineStyle', '--', 'LineWidth', lw)
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-1 31])
ylim([5e-4 1e1])
set(gca, 'YScale', 'log')
box on
legend('env', 'pss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_EpCenv_non_ir_C.pdf")

