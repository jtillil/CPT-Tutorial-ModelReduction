%% Setup
clear; clc;
addpath(genpath("../../../CPT-Tutorial-ModelReduction"))
size = 12;
lw = 1;
lwt = 0.5;

%% Blood Coagulation: 40h Wajima 2009
load("modelBC_SV40_full.mat")

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

% nir-indices
figure
plot(model.t_ref, model.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold), 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
% xlim([-0.14 4.15])
xlim([-0.3 20])
ylim([-0.01 1])
legend([model.analysis.ir.nmstates_above_nindex_threshold; 'threshold'], 'Location','east')
xlabel("t [h]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV40_ir.pdf")

load("modelBC_SV40_t0_01.mat")

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
    legend('env', 'pss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
    xlabel("t [h]")
    ylabel("normalised index")
    hold off
    
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

    exportgraphics(gcf, ['./figures/BC_SV40_ir_' state '.pdf'])
end

%% new ir indices

load("ir_functioning_jac_from_UFa.mat")

% what to plot
threshold = 0.1;
maxnir = max(model.nir, [], 1);
[~, idx_sorted] = sort(maxnir, "descend");
idx_sorted = idx_sorted(maxnir(idx_sorted) > threshold);

% reference solution
figure
grid on
semilogy(model.t_ref, model.X_ref(:, idx_sorted), 'LineWidth', lw) %DisplayName', plotnames(i))
xlim([-2 42])
ylim([1e-7 5e4])
legend(model.I.nmstatelegend(idx_sorted), 'Location', 'southeast')
xlabel("t [h]")
ylabel("concentration [g/L]")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/functioning_jac_from_UFa_ref_sol.pdf")

% ir-indices
figure
plot(model.t_ref, model.nir(:, idx_sorted), 'LineStyle', '-', 'LineWidth', lw)
legend(model.I.nmstatelegend(idx_sorted))
xlim([-0.3 20])
xlabel("t [h]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, ['./figures/functioning_jac_from_UFa_ir_indices.pdf'])





