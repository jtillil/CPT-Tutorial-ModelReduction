%% Setup
clear; clc;
addpath(genpath("../../../CPT-Tutorial-ModelReduction"))
size = 12;
lw = 1;
lwt = 0.5;

%% Blood Coagulation: 40h Wajima 2009
load("modelBC_SV40_from_JKn_2024.mat")
model.I.nmstatelegend = cellfun(@(x) strrep(x, '_', ':'), model.I.nmstatelegend, 'UniformOutput', false);

% reference solution
figure
hold on
for i = 1:7
    semilogy(model.t_ref, model.X_ref(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    semilogy(model.t_ref, model.X_ref(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
xlim([-2 42])
ylim([1e-7 5e4])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'southeast')
xlabel("t [h]")
ylabel("concentration [nmol/L]")
yscale('log')
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_ref_sol.pdf")

% reference solution 0 - 0.15 h
figure
hold on
for i = 1:7
    semilogy(model.t_ref, model.X_ref(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    semilogy(model.t_ref, model.X_ref(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
xlim([-0.005 0.155])
ylim([1e-7 5e4])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'southeast')
xlabel("t [h]")
ylabel("concentration [nmol/L]")
yscale('log')
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_ref_sol_015h.pdf")

% nir-indices
figure
hold on
for i = 1:7
    plot(model.t_ref, model.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    plot(model.t_ref, model.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
yline(0.1, 'k--', 'LineWidth', lwt)
% xlim([-0.14 4.15])
xlim([-0.3 10.3])
ylim([-0.01 1])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'east')
xlabel("t [h]")
ylabel("nir-index")
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_ir.pdf")

% nir-indices 0 - 0.15 h
figure
hold on
for i = 1:7
    plot(model.t_ref, model.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    plot(model.t_ref, model.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
yline(0.1, 'k--', 'LineWidth', lwt)
% xlim([-0.14 4.15])
xlim([-0.005 0.155])
ylim([-0.01 1])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'northeast')
xlabel("t [h]")
ylabel("nir-index")
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_ir_015h.pdf")

%% non-ir

load("modelBC_SV40_from_JKn_2024_t0_01.mat")

% non ir-indices
for i = 1:length(model.analysis.ir.nmstates_above_nindex_threshold)
    state = model.analysis.ir.nmstates_above_nindex_threshold{i};
    idx = model.I.(state);
    figure
    hold on
    semilogy(model.t_ref, model.env.nindex(:, idx), 'LineStyle', '-', 'LineWidth', lw)
    semilogy(model.t_ref, model.pss.nindex(:, idx), 'LineStyle', '--', 'LineWidth', lw)
    semilogy(model.t_ref, model.cneg.nindex(:, idx), 'LineStyle', '-', 'LineWidth', lw)
    semilogy(model.t_ref, model.pneg.nindex(:, idx), 'LineStyle', '--', 'LineWidth', lw)
    yline(0.1, 'k--', 'LineWidth', lwt)
    xlim([-0.14 4.15])
    ylim([5e-4 1e1])
    set(gca, 'YScale', 'log')
    % legend('env', 'pss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
    legend('env', 'pss', 'cneg', 'pneg', 'Location','northeast')
    xlabel("t [h]")
    ylabel("normalised index")
    box on
    hold off
    
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

    exportgraphics(gcf, ['./figures/BC_SV40_ir_' state '.pdf'])
end

%% index-reduced and optimized 8-state model

load("BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv.mat")
model = redmodel;
model.I.nmstatelegend = cellfun(@(x) strrep(x, '_', ':'), model.I.nmstatelegend, 'UniformOutput', false);

% reference solution
figure
hold on
for i = 1:7
    semilogy(model.t_ref, model.X_red(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    semilogy(model.t_ref, model.X_red(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
xlim([-2 42])
ylim([1e-7 5e4])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'southeast')
xlabel("t [h]")
ylabel("concentration [nmol/L]")
yscale('log')
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv_ref_sol.pdf")

% reference solution 0 - 0.15 h
figure
hold on
for i = 1:7
    semilogy(model.t_ref, model.X_red(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    semilogy(model.t_ref, model.X_red(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
xlim([-0.005 0.155])
ylim([1e-7 5e4])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'southeast')
xlabel("t [h]")
ylabel("concentration [nmol/L]")
yscale('log')
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv_ref_sol_015h.pdf")

load("indices_irred_0.05_BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv.mat")

% nir-indices
figure
hold on
for i = 1:7
    plot(model.t_ref, ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    plot(model.t_ref, ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
yline(0.1, 'k--', 'LineWidth', lwt)
% xlim([-0.14 4.15])
xlim([-0.3 10.3])
ylim([-0.01 1])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'east')
xlabel("t [h]")
ylabel("nir-index")
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv_ir.pdf")

% nir-indices 0 - 0.15 h
figure
hold on
for i = 1:7
    plot(model.t_ref, ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    plot(model.t_ref, ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
yline(0.1, 'k--', 'LineWidth', lwt)
% xlim([-0.14 4.15])
xlim([-0.005 0.155])
ylim([-0.01 1])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'northeast')
xlabel("t [h]")
ylabel("nir-index")
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv_ir_015h.pdf")

%% index-reduced model - TODO UPDATE IR TO READ FROM MODEL ONCE SAVED

model.I.nmstatelegend = cellfun(@(x) strrep(x, '_', ':'), model.I.nmstatelegend, 'UniformOutput', false);

% nir-indices
figure
hold on
for i = 1:7
    plot(model.t_ref, ir.nindex(1:(end-2), model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    plot(model.t_ref, ir.nindex(1:(end-2), model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
yline(0.1, 'k--', 'LineWidth', lwt)
% xlim([-0.14 4.15])
xlim([-0.3 10.3])
ylim([-0.01 1])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'east')
xlabel("t [h]")
ylabel("nir-index")
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/irred_0.05_BCSV40_from_JKn_pss_solved_ir.pdf")

% nir-indices 0 - 0.15 h
figure
hold on
for i = 1:7
    plot(model.t_ref, ir.nindex(1:(end-2), model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    plot(model.t_ref, ir.nindex(1:(end-2), model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
yline(0.1, 'k--', 'LineWidth', lwt)
% xlim([-0.14 4.15])
xlim([-0.005 0.155])
ylim([-0.01 1])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'northeast')
xlabel("t [h]")
ylabel("nir-index")
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/irred_0.05_BCSV40_from_JKn_pss_solved_ir_015h.pdf")