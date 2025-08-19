%% Setup
clear; clc;
addpath(genpath("../../jti-code"))
sz = 12;
lw = 1;
lwt = 0.5;
interpreter = 'latex';
% Set default interpreter to LaTeX
set(groot, 'defaultTextInterpreter', interpreter);
set(groot, 'defaultAxesTickLabelInterpreter', interpreter);
set(groot, 'defaultLegendInterpreter', interpreter);
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

%% Blood Coagulation: 40h Wajima 2009
% load("modelBC_SV40_from_JKn_2024.mat")
load("res.mat")
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
xlim([-0.05 1])
ylim([-0.01 1])
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'east')
xlabel("t [h]")
ylabel("nir-index")
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV1_ir.pdf")

% nir-indices 0015
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
legend(model.I.nmstatelegend(model.analysis.ir.I_sorted_max_nindex_above_threshold), 'Location', 'east')
xlabel("t [h]")
ylabel("nir-index")
box on
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size, size]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BC_SV1_ir_0015.pdf")

% nir-indices 0 - 0.15 h
figure
hold on
for i = 1:7
    plot(model.t_ref, model.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), 'LineWidth', lw) %DisplayName', plotnames(i))
end
for i = 8:length(model.analysis.ir.I_sorted_max_nindex_above_threshold)
    plot(model.t_ref, model.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
end
plot(model.t_ref, model.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex_above_threshold(i)), '--', 'LineWidth', lw) %DisplayName', plotnames(i))
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

%% Blood Coagulation: Combined Figure (Adjusted for consistent styling)

% ADJUSTED: Define constants for clarity and easy modification
lw = 1.5; % Line width for plots
lwt = 1;  % Line width for threshold lines

% Create figure
figure('Units', 'centimeters', 'Position', [0, 0, 20, 5.5]);

% Create tiled layout
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% ADJUSTED: Pre-define title and label strings for cleaner code
% Define title strings for each subplot
title_strings = {
    "A) Concentration-time profiles", "B) Normalized ir-indices", "C) Normalized ir-indices in first 0.15h";
    "D) C-t profiles (reduced)", "E) Norm. ir-indices (reduced)", "F) Norm. ir-indices in first 0.15h (reduced)"
};

% Define y-labels (only need one for each type)
ylabel_strings = {
    "concentration [nM]",
    "nir-index"
};

% --- Loop to create all plots ---
% for row = 1:2
for row = 1
    for col = 1:3
        % Load data based on the row
        if row == 1
            load("modelBC_SV40_from_JKn_2024.mat")
            nmstatelegend = cellfun(@(x) strrep(x, '_', ':'), model.I.nmstatelegend, 'UniformOutput', false);
        else
            load("BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv.mat")
            load("indices_irred_0.05_BCSV40_from_JKn_exh_t120_MRSE_pnegrun_greedy_0.05_10_linear_dyncnegpnegenv.mat")
            % We need nmstatelegend for the second row too if we want a legend
            % Assuming the states are the same, we can load the full model just for this
            full_model = load("modelBC_SV40_from_JKn_2024.mat");
            nmstatelegend = cellfun(@(x) strrep(x, '_', ':'), full_model.model.I.nmstatelegend, 'UniformOutput', false);
        end

        % Create next tile
        ax = nexttile;
        hold(ax, 'on');
        
        % idx to plot
        idx_sorted_to_plot = setdiff(model.analysis.ir.I_sorted_max_nindex_above_threshold, [model.I.VIII, model.I.VIIIa], 'stable');
        idx_sorted_to_plot = [idx_sorted_to_plot model.I.IIa];

        % Plot data based on column
        if col == 1 % Concentration profiles
            if row == 1
                for i = 1:7
                    semilogy(ax, model.t_ref, model.X_ref(:, idx_sorted_to_plot(i)), 'LineWidth', lw);
                end
                for i = 8:length(idx_sorted_to_plot)
                    semilogy(ax, model.t_ref, model.X_ref(:, idx_sorted_to_plot(i)), '--', 'LineWidth', lw);
                end
            elseif row == 2
                for i = 1:7
                    semilogy(ax, redmodel.t_ref, redmodel.X_red(:, idx_sorted_to_plot(i)), 'LineWidth', lw);
                end
                for i = 8:length(idx_sorted_to_plot)
                    semilogy(ax, redmodel.t_ref, redmodel.X_red(:, idx_sorted_to_plot(i)), '--', 'LineWidth', lw);
                end
            end
            % xlim(ax, [-2 42]);
            xlim(ax, [-0.3 10.3]);
            ylim(ax, [1e-6 1e5]);
            yticks([1e-5, 1e0, 1e5]);
            set(ax, 'YScale', 'log');
            
        elseif col == 2 % ir-indices (full time)
            if row == 1
                for i = 1:7
                    plot(ax, model.t_ref, model.ir.nindex(:, idx_sorted_to_plot(i)), 'LineWidth', lw);
                end
                for i = 8:length(idx_sorted_to_plot)
                    plot(ax, model.t_ref, model.ir.nindex(:, idx_sorted_to_plot(i)), '--', 'LineWidth', lw);
                end
            elseif row == 2
                for i = 1:7
                    plot(ax, model.t_ref, ir.nindex(:, idx_sorted_to_plot(i)), 'LineWidth', lw);
                end
                for i = 8:length(idx_sorted_to_plot)
                    plot(ax, model.t_ref, ir.nindex(:, idx_sorted_to_plot(i)), '--', 'LineWidth', lw);
                end
            end
            yline(ax, 0.1, 'k--', 'LineWidth', lwt);
            xlim(ax, [-0.3 10.3]);
            ylim(ax, [-0.01 1]);
            
        elseif col == 3 % ir-indices (zoomed in)
            if row == 1
                for i = 1:7
                    plot(ax, model.t_ref, model.ir.nindex(:, idx_sorted_to_plot(i)), 'LineWidth', lw);
                end
                for i = 8:length(idx_sorted_to_plot)
                    plot(ax, model.t_ref, model.ir.nindex(:, idx_sorted_to_plot(i)), '--', 'LineWidth', lw);
                end
            elseif row == 2
                for i = 1:7
                    plot(ax, model.t_ref, ir.nindex(:, idx_sorted_to_plot(i)), 'LineWidth', lw);
                end
                for i = 8:length(idx_sorted_to_plot)
                    plot(ax, model.t_ref, ir.nindex(:, idx_sorted_to_plot(i)), '--', 'LineWidth', lw);
                end
            end
            yline(ax, 0.1, 'k--', 'LineWidth', lwt);
            xlim(ax, [-0.005 0.155]);
            ylim(ax, [-0.01 1]);
        end
        
        box(ax, 'on');
        hold(ax, 'off');
        
        % ADJUSTED: Use text() for titles for precise, left-aligned placement
        % This replaces the old title() function. The coordinates are normalized,
        % with (0,1) being the top-left corner of the axes.
        % We adjust the horizontal position for a consistent look across columns.
        if col == 1
            text(ax, -0.28, 1.13, title_strings{row, col}, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');
        else
            text(ax, -0.20, 1.13, title_strings{row, col}, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');
        end

        % Y-labels only for left column
        if col == 1
            ylabel(ylabel_strings{1});
        elseif col == 2
            ylabel(ylabel_strings{2});
        end
        
        % X-labels only for the bottom row
        if row == 1
            xlabel("t [h]");
        else
            % Remove x-tick labels for non-bottom rows for a cleaner look
            set(ax, 'XTickLabel', []);
        end

        grid on
        box on
        
        % Legends only for the right-most column
        if col == 3
             % We only add a legend to the first row to avoid redundancy
             if row == 1
                lg = legend(ax, nmstatelegend(idx_sorted_to_plot), 'Location', 'eastoutside');
                lg.FontSize = 8;
             end
             % For the second row, a legend would be redundant if plotting the same species
        end
    end
end

% ADJUSTED: Removed manual position adjustment.
% The 'tiledlayout' with 'compact' spacing and 'eastoutside' legends
% generally handles alignment well, making this manual step unnecessary and
% sometimes problematic.

% Export
% ADJUSTED: Added 'ContentType','vector' for a high-quality, scalable PDF
exportgraphics(gcf, "./figures/BC_complete.pdf", 'ContentType', 'vector');

%% combined layout analysis

% Load data first, as it's needed for both plot sections
load("modelBC_SV40_from_JKn_2024.mat")

% --- Prepare Data for All Plots ---
Xdat = [...
    max(model.env.nindex(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.env.relstateerr(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.pss.nindex(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.pss.relstateerr(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.cneg.nindex(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.pneg.nindex(:, model.analysis.ir.I_sorted_max_nindex))'];
Xdat(isnan(Xdat)) = 10;
names = model.I.nmstate(model.analysis.ir.I_sorted_max_nindex);
names = replace(names, '_', ':');
names{46} = 'VKp';
names{54} = 'ENOp';
names{60} = 'AVenom\_Tiger';
names{61} = 'CVenom\_Tiger';
dat = max(model.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex))';
dat([30 31]) = 0;

AUC_t_idx = 298;
dat_AUC = trapz(model.t_ref(1:AUC_t_idx), model.ir.nindex(1:AUC_t_idx, model.analysis.ir.I_sorted_max_nindex), 1);
sum_AUC = sum(dat_AUC);
dat_AUC_norm = dat_AUC./sum_AUC;
dat_AUC_norm([30 31]) = 0;

% --- 1. Setup Figure and Panels ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 15]);

mainLayout = tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- 2. Create the Top Scatter Plot inside the Top Panel ---
ax1 = nexttile(mainLayout, [1,1]);

myblue = [0.0660, 0.4430, 0.7450];
myorange = [0.8660, 0.3290, 0.0000];
scatter(ax1, 1:length(dat), dat, 36, 'filled', 'MarkerFaceColor', 'black');
hold(ax1, 'on');
scatter(ax1, 1:length(dat), dat_AUC_norm, 'ko', 'MarkerFaceColor', 'white');
hold(ax1, 'off');
legend(["maximum", "fraction of AUC"], 'Location', 'northeast')
box(ax1, 'on');
xlim(ax1, [0 length(dat)+1]);
ylim(ax1, [1e-4 2e0]);
ylabel(ax1, 'nir-index');
title_text_scatter = "A) Ordering of states by decreasing dynamic relevance";
ax1.YTick = [1e-4 1e-3 1e-2 1e-1 1e0];
ax1.XTick = 1:length(names);
ax1.XTickLabel = names(1:end);
ax1.XTickLabelRotation = 60;
ax1.XAxis.FontSize = 8;
ax1.TickLength = [0.005 0.005];
set(ax1, 'YScale', 'log');
grid(ax1, 'on');
text(ax1, -0.08, 1.19, title_text_scatter, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');


% --- 3. Create the Bottom Heatmap Plots inside the Bottom Panel ---
bottomLayout = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
bottomLayout.Layout.Tile = 2;
bottomLayout.Layout.TileSpan = [2, 1];

% --- Define colors and column names for the heatmaps ---
colNames = {'env', 'err (env)', 'pss', 'err (pss)', 'cneg', 'pneg'};
mypurple = [0.4940, 0.1840, 0.5560];
mygold   = [0.9290, 0.6940, 0.1250];
white    = [1.0000, 1.0000, 1.0000];
cmap = [white; myblue; myorange; mygold; mypurple];
row_splits = [1, 22, 43, length(names)+1];

% Loop to create the three heatmaps
for i = 1:3
    row_indices = row_splits(i) : row_splits(i+1)-1;
    Xdat_subset = Xdat(row_indices, :);
    names_subset = names(row_indices);
    [num_rows, num_cols] = size(Xdat_subset);

    % --- This part remains the same: Determine the color index for each cell ---
    color_index_matrix = ones(num_rows, num_cols);
    color_index_matrix(Xdat_subset(:, 1) < 0.1, 1) = 2;
    color_index_matrix(Xdat_subset(:, 2) < 0.1, 2) = 2;
    color_index_matrix(Xdat_subset(:, 3) < 0.1, 3) = 3;
    color_index_matrix(Xdat_subset(:, 4) < 0.1, 4) = 3;
    color_index_matrix(Xdat_subset(:, 5) < 0.1, 5) = 4;
    color_index_matrix(Xdat_subset(:, 6) < 0.1, 6) = 5;

    ax_heat = nexttile(bottomLayout);
    
    % --- MODIFIED SECTION: Draw Colored Checkmarks on White Grid ---
    % We will create an empty plot with a white background and draw our own content.
    % Set axes limits manually, as the 'image' function is no longer doing it.
    xlim(ax_heat, [0.5, num_cols + 0.5]);
    ylim(ax_heat, [0.5, num_rows + 0.5]);
    % Invert the Y-axis to match the behavior of 'image' (row 1 at the top)
    set(ax_heat, 'YDir', 'reverse');

    hold(ax_heat, 'on');

    % --- NEW LOGIC: Loop through non-white cells and plot colored checkmarks ---
    % Find the coordinates of all cells that should not be white (index ~= 1)
    [rows, cols] = find(color_index_matrix ~= 1);
    
    % Loop through each of these coordinates
    for k = 1:length(rows)
        % Get the specific row and column for this checkmark
        r = rows(k);
        c = cols(k);
        
        % Get the color index (e.g., 2 for blue, 3 for orange, etc.)
        color_idx = color_index_matrix(r, c);
        
        % Look up the actual RGB color from our custom colormap
        checkmark_color = cmap(color_idx, :);
        
        % Place the checkmark with the specific color
        text(ax_heat, c, r, '✓', ...
            'Color', checkmark_color, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 12, ...
            'FontWeight', 'bold');
    end
    
    hold(ax_heat, 'off');
    
    % --- The rest of the formatting code remains the same ---
    ax_heat.XTick = 1:num_cols;
    ax_heat.XTickLabel = colNames;
    ax_heat.XTickLabelRotation = 45;
    
    ax_heat.YTick = 1:num_rows;
    ax_heat.YTickLabel = names_subset;
    
    ax_heat.TickLength = [0 0];
    x_grid_lines = 0.5 : 1 : num_cols + 0.5;
    y_grid_lines = 0.5 : 1 : num_rows + 0.5;
    ax_heat.XAxis.MinorTickValues = x_grid_lines;
    ax_heat.YAxis.MinorTickValues = y_grid_lines;
    grid(ax_heat, 'off'); % Turn off major grid
    grid(ax_heat); % Turn on major grid again to apply settings
    ax_heat.GridAlpha = 0; % Make major grid invisible
    ax_heat.XMinorGrid = 'on';
    ax_heat.YMinorGrid = 'on';
    ax_heat.MinorGridColor = [0 0 0]; % Black grid lines
    ax_heat.MinorGridAlpha = 1; % Fully opaque
    ax_heat.MinorGridLineStyle = '-'; 
    ax_heat.Layer = 'top'; % Ensure grid is drawn over any background color

    if i == 1
        text(ax_heat, -0.33, 1.07, "B) Analysis of reduction approaches", ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');
    end
end

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_analysis_combined.pdf", 'ContentType', 'vector');

%% combined layout analysis with colored background

% Load data first, as it's needed for both plot sections
load("modelBC_SV40_from_JKn_2024.mat")

% --- Prepare Data for All Plots ---
Xdat = [...
    max(model.env.nindex(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.env.relstateerr(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.pss.nindex(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.pss.relstateerr(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.cneg.nindex(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.pneg.nindex(:, model.analysis.ir.I_sorted_max_nindex))'];
Xdat(isnan(Xdat)) = 10;
names = model.I.nmstate(model.analysis.ir.I_sorted_max_nindex);
names = replace(names, '_', ':');
names{46} = 'VKp';
names{54} = 'ENOp';
names{60} = 'AVenom\_Tiger';
names{61} = 'CVenom\_Tiger';
dat = max(model.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex))';
dat([30 31]) = 0;

AUC_t_idx = 298;
dat_AUC = trapz(model.t_ref(1:AUC_t_idx), model.ir.nindex(1:AUC_t_idx, model.analysis.ir.I_sorted_max_nindex), 1);
sum_AUC = sum(dat_AUC);
dat_AUC_norm = dat_AUC./sum_AUC;
dat_AUC_norm([30 31]) = 0;

% --- 1. Setup Figure and Panels ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 15]);

mainLayout = tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- 2. Create the Top Scatter Plot inside the Top Panel ---
ax1 = nexttile(mainLayout, [1,1]);

myblue = [0.0660, 0.4430, 0.7450];
myorange = [0.8660, 0.3290, 0.0000];
scatter(ax1, 1:length(dat), dat, 36, 'filled', 'MarkerFaceColor', 'black');
hold(ax1, 'on');
% scatter(ax1, 1:length(dat), dat_AUC_norm, 'ko', 'MarkerFaceColor', 'white');
hold(ax1, 'off');
% legend(["maximum", "fraction of AUC"], 'Location', 'northeast')
box(ax1, 'on');
xlim(ax1, [0 length(dat)+1]);
ylim(ax1, [1e-4 2e0]);
ylabel(ax1, 'nir-index');
title_text_scatter = "A) Ordering of states by decreasing dynamic relevance";
ax1.YTick = [1e-4 1e-3 1e-2 1e-1 1e0];
ax1.XTick = 1:length(names);
ax1.XTickLabel = names(1:end);
ax1.XTickLabelRotation = 60;
ax1.XAxis.FontSize = 8;
ax1.TickLength = [0.005 0.005];
set(ax1, 'YScale', 'log');
grid(ax1, 'on');
text(ax1, -0.08, 1.19, title_text_scatter, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');


% --- 3. Create the Bottom Heatmap Plots inside the Bottom Panel ---
bottomLayout = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
bottomLayout.Layout.Tile = 2;
bottomLayout.Layout.TileSpan = [2, 1];

% --- Define colors and column names for the heatmaps ---
colNames = {'env', 'err (env)', 'pss', 'err (pss)', 'cneg', 'pneg'};
mypurple = [0.4940, 0.1840, 0.5560];
mygold   = [0.9290, 0.6940, 0.1250];
white    = [1.0000, 1.0000, 1.0000];
cmap = [white; myblue; myorange; mygold; mypurple];
row_splits = [1, 22, 43, length(names)+1];

% --- Define rows and colors for background highlighting ---
gray_rows_global = [1, 2, 3, 7, 9, 10, 13, 23];
blue_rows_global = [5, 15, 19, 21, 26];
% Calculate a 20% gray (i.e., 80% white)
light_gray = [0.8, 0.8, 0.8]; 
% Calculate a 20% tint of 'myblue' mixed with 80% white
light_blue = 0.2 * myblue + 0.8 * white;

% Loop to create the three heatmaps
for i = 1:3
    row_indices = row_splits(i) : row_splits(i+1)-1;
    Xdat_subset = Xdat(row_indices, :);
    names_subset = names(row_indices);
    [num_rows, num_cols] = size(Xdat_subset);

    color_index_matrix = ones(num_rows, num_cols);
    color_index_matrix(Xdat_subset(:, 1) < 0.1, 1) = 2;
    color_index_matrix(Xdat_subset(:, 2) < 0.1, 2) = 2;
    color_index_matrix(Xdat_subset(:, 3) < 0.1, 3) = 3;
    color_index_matrix(Xdat_subset(:, 4) < 0.1, 4) = 3;
    color_index_matrix(Xdat_subset(:, 5) < 0.1, 5) = 4;
    color_index_matrix(Xdat_subset(:, 6) < 0.1, 6) = 5;

    ax_heat = nexttile(bottomLayout);
    
    xlim(ax_heat, [0.5, num_cols + 0.5]);
    ylim(ax_heat, [0.5, num_rows + 0.5]);
    set(ax_heat, 'YDir', 'reverse');

    hold(ax_heat, 'on');

    % --- ADDED CODE: Draw background highlights for specific rows --- %%
    current_gray_rows = intersect(gray_rows_global, row_indices);
    current_blue_rows = intersect(blue_rows_global, row_indices);

    % Draw gray background patches
    for k = 1:length(current_gray_rows)
        global_row = current_gray_rows(k);
        local_row = global_row - row_indices(1) + 1;
        patch(ax_heat, [0.5, num_cols + 0.5, num_cols + 0.5, 0.5], ...
                       [local_row - 0.5, local_row - 0.5, local_row + 0.5, local_row + 0.5], ...
                       light_gray, 'EdgeColor', 'none');
    end
    
    % Draw blue background patches
    for k = 1:length(current_blue_rows)
        global_row = current_blue_rows(k);
        local_row = global_row - row_indices(1) + 1;
        patch(ax_heat, [0.5, num_cols + 0.5, num_cols + 0.5, 0.5], ...
                       [local_row - 0.5, local_row - 0.5, local_row + 0.5, local_row + 0.5], ...
                       light_blue, 'EdgeColor', 'none');
    end
    % --- END OF ADDED CODE --- %%

    [rows, cols] = find(color_index_matrix ~= 1);
    
    for k = 1:length(rows)
        r = rows(k);
        c = cols(k);
        color_idx = color_index_matrix(r, c);
        checkmark_color = cmap(color_idx, :);
        text(ax_heat, c, r, '✓', ...
            'Color', checkmark_color, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 12, ...
            'FontWeight', 'bold');
    end
    
    hold(ax_heat, 'off');
    
    ax_heat.XTick = 1:num_cols;
    ax_heat.XTickLabel = colNames;
    ax_heat.XTickLabelRotation = 45;
    
    ax_heat.YTick = 1:num_rows;
    ax_heat.YTickLabel = names_subset;
    
    ax_heat.TickLength = [0 0];
    x_grid_lines = 0.5 : 1 : num_cols + 0.5;
    y_grid_lines = 0.5 : 1 : num_rows + 0.5;
    ax_heat.XAxis.MinorTickValues = x_grid_lines;
    ax_heat.YAxis.MinorTickValues = y_grid_lines;
    grid(ax_heat, 'off');
    grid(ax_heat);
    ax_heat.GridAlpha = 0;
    ax_heat.XMinorGrid = 'on';
    ax_heat.YMinorGrid = 'on';
    ax_heat.MinorGridColor = [0 0 0];
    ax_heat.MinorGridAlpha = 1;
    ax_heat.MinorGridLineStyle = '-'; 
    ax_heat.Layer = 'top';

    if i == 1
        text(ax_heat, -0.33, 1.07, "B) Analysis of reduction approaches", ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');
    end
end

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_analysis_combined.pdf", 'ContentType', 'vector');

%% combined layout analysis with colored background, separate columns for models

% Load data first, as it's needed for both plot sections
load("modelBC_SV40_from_JKn_2024.mat")

% --- Prepare Data for All Plots ---
Xdat = [...
    max(model.env.nindex(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.env.relstateerr(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.pss.nindex(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.pss.relstateerr(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.cneg.nindex(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.pneg.nindex(:, model.analysis.ir.I_sorted_max_nindex))'];
Xdat(isnan(Xdat)) = 10;
names = model.I.nmstate(model.analysis.ir.I_sorted_max_nindex);
names = replace(names, '_', ':');
names{46} = 'VKp';
names{54} = 'ENOp';
names{60} = 'AVenom\_Tiger';
names{61} = 'CVenom\_Tiger';
dat = max(model.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex))';
dat([30 31]) = 0;

AUC_t_idx = 298;
dat_AUC = trapz(model.t_ref(1:AUC_t_idx), model.ir.nindex(1:AUC_t_idx, model.analysis.ir.I_sorted_max_nindex), 1);
sum_AUC = sum(dat_AUC);
dat_AUC_norm = dat_AUC./sum_AUC;
dat_AUC_norm([30 31]) = 0;

% --- 1. Setup Figure and Panels ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 15]);

mainLayout = tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- 2. Create the Top Scatter Plot inside the Top Panel ---
ax1 = nexttile(mainLayout, [1,1]);

myblue = [0.0660, 0.4430, 0.7450];
myorange = [0.8660, 0.3290, 0.0000];
scatter(ax1, 1:length(dat), dat, 36, 'filled', 'MarkerFaceColor', 'black');
% hold(ax1, 'on');
% scatter(ax1, 1:length(dat), dat_AUC_norm, 'ko', 'MarkerFaceColor', 'white');
% hold(ax1, 'off');
% legend(["maximum", "fraction of AUC"], 'Location', 'northeast')
box(ax1, 'on');
xlim(ax1, [0 length(dat)+1]);
ylim(ax1, [1e-4 2e0]);
ylabel(ax1, 'nir-index');
title_text_scatter = "A) Ordering of states by decreasing dynamic relevance";
ax1.YTick = [1e-4 1e-3 1e-2 1e-1 1e0];
ax1.XTick = 1:length(names);
ax1.XTickLabel = names(1:end);
ax1.XTickLabelRotation = 60;
ax1.XAxis.FontSize = 8;
ax1.TickLength = [0.005 0.005];
set(ax1, 'YScale', 'log');
grid(ax1, 'on');
text(ax1, -0.08, 1.19, title_text_scatter, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');


% --- 3. Create the Bottom Heatmap Plots inside the Bottom Panel ---
bottomLayout = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
bottomLayout.Layout.Tile = 2;
bottomLayout.Layout.TileSpan = [2, 1];

% --- Define colors and column names for the heatmaps ---
colNames = {'env', 'err (env)', 'pss', 'err (pss)', 'cneg', 'pneg', '8-state', '13-state'};
mypurple = [0.4940, 0.1840, 0.5560];
mygold   = [0.9290, 0.6940, 0.1250];
white    = [1.0000, 1.0000, 1.0000];
cmap = [white; myblue; myorange; mygold; mypurple];
row_splits = [1, 22, 43, length(names)+1];

% --- Define rows and colors for background highlighting ---
% For '8-state' column (column 7)
gray_rows_8 = [1, 2, 3, 7, 9, 10, 13, 23]; 
blue_rows_8 = [5, 15, 19, 21, 26];       

% For '13-state' column (column 8)
gray_rows_13 = [1, 2, 3, 5, 6, 8, 11, 12, 13, 14, 17, 18, 20];
blue_rows_13 = [4, 16, 22, 24];

% Calculate tinted colors
light_gray = [0.8, 0.8, 0.8]; 
light_gold = 0.2 * mygold + 0.8 * white;
light_blue = 0.2 * myblue + 0.8 * white;

% Loop to create the three heatmaps
for i = 1:3
    row_indices = row_splits(i) : row_splits(i+1)-1;
    Xdat_subset = Xdat(row_indices, :);
    names_subset = names(row_indices);
    [num_rows, num_data_cols] = size(Xdat_subset);
    num_total_cols = num_data_cols + 2; % 6 data columns + 2 model columns

    color_index_matrix = ones(num_rows, num_data_cols);
    color_index_matrix(Xdat_subset(:, 1) < 0.1, 1) = 2;
    color_index_matrix(Xdat_subset(:, 2) < 0.1, 2) = 2;
    color_index_matrix(Xdat_subset(:, 3) < 0.1, 3) = 3;
    color_index_matrix(Xdat_subset(:, 4) < 0.1, 4) = 3;
    color_index_matrix(Xdat_subset(:, 5) < 0.1, 5) = 4;
    color_index_matrix(Xdat_subset(:, 6) < 0.1, 6) = 5;

    ax_heat = nexttile(bottomLayout);
    
    xlim(ax_heat, [0.5, num_total_cols + 0.5]);
    ylim(ax_heat, [0.5, num_rows + 0.5]);
    set(ax_heat, 'YDir', 'reverse');

    hold(ax_heat, 'on');

    % --- MODIFIED CODE: Separate logic for column 7 and column 8 --- %%
    all_local_rows = 1:num_rows;

    % --- Coloring for Column 7 ('8-state') ---
    local_blue_rows_8 = find(ismember(row_indices, blue_rows_8));
    local_gray_rows_8 = find(ismember(row_indices, gray_rows_8));
    local_gold_rows_8 = setdiff(all_local_rows, [local_blue_rows_8(:)' local_gray_rows_8(:)']);
    
    for k = 1:length(local_gold_rows_8); patch(ax_heat, [6.5, 7.5, 7.5, 6.5], [local_gold_rows_8(k) - 0.5, local_gold_rows_8(k) - 0.5, local_gold_rows_8(k) + 0.5, local_gold_rows_8(k) + 0.5], light_gold, 'EdgeColor', 'none'); end
    for k = 1:length(local_blue_rows_8); patch(ax_heat, [6.5, 7.5, 7.5, 6.5], [local_blue_rows_8(k) - 0.5, local_blue_rows_8(k) - 0.5, local_blue_rows_8(k) + 0.5, local_blue_rows_8(k) + 0.5], light_blue, 'EdgeColor', 'none'); end
    for k = 1:length(local_gray_rows_8); patch(ax_heat, [6.5, 7.5, 7.5, 6.5], [local_gray_rows_8(k) - 0.5, local_gray_rows_8(k) - 0.5, local_gray_rows_8(k) + 0.5, local_gray_rows_8(k) + 0.5], light_gray, 'EdgeColor', 'none'); end

    % --- Coloring for Column 8 ('13-state') ---
    local_blue_rows_13 = find(ismember(row_indices, blue_rows_13));
    local_gray_rows_13 = find(ismember(row_indices, gray_rows_13));
    local_gold_rows_13 = setdiff(all_local_rows, [local_blue_rows_13(:)' local_gray_rows_13(:)']);

    for k = 1:length(local_gold_rows_13); patch(ax_heat, [7.5, 8.5, 8.5, 7.5], [local_gold_rows_13(k) - 0.5, local_gold_rows_13(k) - 0.5, local_gold_rows_13(k) + 0.5, local_gold_rows_13(k) + 0.5], light_gold, 'EdgeColor', 'none'); end
    for k = 1:length(local_blue_rows_13); patch(ax_heat, [7.5, 8.5, 8.5, 7.5], [local_blue_rows_13(k) - 0.5, local_blue_rows_13(k) - 0.5, local_blue_rows_13(k) + 0.5, local_blue_rows_13(k) + 0.5], light_blue, 'EdgeColor', 'none'); end
    for k = 1:length(local_gray_rows_13); patch(ax_heat, [7.5, 8.5, 8.5, 7.5], [local_gray_rows_13(k) - 0.5, local_gray_rows_13(k) - 0.5, local_gray_rows_13(k) + 0.5, local_gray_rows_13(k) + 0.5], light_gray, 'EdgeColor', 'none'); end
    % --- END OF MODIFIED CODE --- %%

    [rows, cols] = find(color_index_matrix ~= 1);
    
    for k = 1:length(rows)
        r = rows(k);
        c = cols(k);
        color_idx = color_index_matrix(r, c);
        checkmark_color = cmap(color_idx, :);
        text(ax_heat, c, r, '✓', ...
            'Color', checkmark_color, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 12, ...
            'FontWeight', 'bold');
    end
    
    hold(ax_heat, 'off');
    
    ax_heat.XTick = 1:num_total_cols;
    ax_heat.XTickLabel = colNames;
    ax_heat.XTickLabelRotation = 45;
    
    ax_heat.YTick = 1:num_rows;
    ax_heat.YTickLabel = names_subset;
    
    ax_heat.TickLength = [0 0];
    x_grid_lines = 0.5 : 1 : num_total_cols + 0.5;
    y_grid_lines = 0.5 : 1 : num_rows + 0.5;
    ax_heat.XAxis.MinorTickValues = x_grid_lines;
    ax_heat.YAxis.MinorTickValues = y_grid_lines;
    grid(ax_heat, 'off');
    grid(ax_heat);
    ax_heat.GridAlpha = 0;
    ax_heat.XMinorGrid = 'on';
    ax_heat.YMinorGrid = 'on';
    ax_heat.MinorGridColor = [0 0 0];
    ax_heat.MinorGridAlpha = 1;
    ax_heat.MinorGridLineStyle = '-'; 
    ax_heat.Layer = 'top';

    if i == 1
        text(ax_heat, -0.33, 1.07, "B) Analysis of reduction approaches", ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');
    end
end

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_analysis_combined.pdf", 'ContentType', 'vector');

%% combined layout analysis with colored background, separate columns for models

% Load data first, as it's needed for both plot sections
load("modelBC_SV40_from_JKn_2024.mat")

% --- Prepare Data for All Plots ---
Xdat = [...
    max(model.env.nindex(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.env.relstateerr(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.pss.nindex(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.pss.relstateerr(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.cneg.nindex(:, model.analysis.ir.I_sorted_max_nindex))',...
    max(model.pneg.nindex(:, model.analysis.ir.I_sorted_max_nindex))'];
Xdat(isnan(Xdat)) = 10;
names = model.I.nmstate(model.analysis.ir.I_sorted_max_nindex);
names = replace(names, '_', ':');
names{46} = 'VKp';
names{54} = 'ENOp';
names{60} = 'AVenom\_Tiger';
names{61} = 'CVenom\_Tiger';
dat = max(model.ir.nindex(:, model.analysis.ir.I_sorted_max_nindex))';
dat([30 31]) = 0;

AUC_t_idx = 298;
dat_AUC = trapz(model.t_ref(1:AUC_t_idx), model.ir.nindex(1:AUC_t_idx, model.analysis.ir.I_sorted_max_nindex), 1);
sum_AUC = sum(dat_AUC);
dat_AUC_norm = dat_AUC./sum_AUC;
dat_AUC_norm([30 31]) = 0;

% --- 1. Setup Figure and Panels ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 20]);

mainLayout = tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- 2. Create the Top Scatter Plot with a Broken Y-Axis ---
% This is more complex than a single call to nexttile. We create a
% placeholder tile to get its position, then delete it and create two
% new axes in its place.

% Get position for the top plot area
ax_placeholder = nexttile(mainLayout, [1,1]);
pos = get(ax_placeholder, 'Position');
delete(ax_placeholder);

% Define heights for the two axes and the gap between them
bottom_ax_height = pos(4) * 0.10; % Bottom 15% for the linear scale
top_ax_height = pos(4) * 0.58;    % Top 78% for the log scale
gap_height = pos(4) * 0.04;       % A 7% gap

% Create the bottom axes (linear scale)
ax_bottom_pos = [pos(1), pos(2), pos(3), bottom_ax_height];
ax_bottom = axes('Position', ax_bottom_pos);

scatter(ax_bottom, 1:length(dat), dat, 36, 'filled', 'MarkerFaceColor', 'black');
box(ax_bottom, 'on');
ylim(ax_bottom, [-0.1e-10 1e-10]); % Set linear Y-limits
ax_bottom.YTick = [0 1e-10];
ax_bottom.YTickLabel = {'0', '1e-10'}; % Custom labels for this small section
grid(ax_bottom, 'on');
ax_bottom.XAxis.FontSize = 8;
ax_bottom.XTick = 1:length(names);
ax_bottom.XTickLabel = names(1:end);
ax_bottom.XTickLabelRotation = 60;
ax_bottom.TickLength = [0.005 0.005];

% Create the top axes (log scale)
ax_top_pos = [pos(1), pos(2) + bottom_ax_height + gap_height, pos(3), top_ax_height];
ax_top = axes('Position', ax_top_pos);

scatter(ax_top, 1:length(dat), dat, 36, 'filled', 'MarkerFaceColor', 'black');
box(ax_top, 'on');
ylim(ax_top, [1e-4 2e0]); % Set log Y-limits
set(ax_top, 'YScale', 'log');
set(ax_top, 'XTick', []); % Remove x-ticks from top plot
ylabel(ax_top, 'nir-index');
title_text_scatter = "A) Ordering of states by decreasing dynamic relevance";
ax_top.YTick = [1e-4 1e-3 1e-2 1e-1 1e0];
grid(ax_top, 'on');
ax_top.TickLength = [0.005 0.005];

% Link the X-axes of the two plots so they zoom and pan together
linkaxes([ax_top, ax_bottom], 'x');
xlim(ax_top, [0 length(dat)+1]); % Set x-limit on one to affect both

% Add the title to the top axes
text(ax_top, -0.08, 1.19, title_text_scatter, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');

% Draw the break marks in the gap (requires normalized figure units)
drawnow; % Ensure positions are finalized before drawing annotations
y_break_pos = ax_top.Position(2) - 0.015; % Y position for break marks
break_width = 0.007; % How wide the 'tick' is
% Left break mark
annotation('line', [ax_top.Position(1) ax_top.Position(1)], [y_break_pos-0.01 y_break_pos+0.01], 'LineWidth', 0.75);
annotation('line', [ax_top.Position(1)-break_width ax_top.Position(1)+break_width], [y_break_pos-0.01-0.003 y_break_pos-0.01+0.003], 'LineWidth', 0.75);
annotation('line', [ax_top.Position(1)-break_width ax_top.Position(1)+break_width], [y_break_pos+0.01-0.003 y_break_pos+0.01+0.003], 'LineWidth', 0.75);
% Right break mark
x_right_break = ax_top.Position(1) + ax_top.Position(3);
annotation('line', [x_right_break x_right_break], [y_break_pos-0.01 y_break_pos+0.01], 'LineWidth', 0.75);
annotation('line', [x_right_break-break_width x_right_break+break_width], [y_break_pos-0.01-0.003 y_break_pos-0.01+0.003], 'LineWidth', 0.75);
annotation('line', [x_right_break-break_width x_right_break+break_width], [y_break_pos+0.01-0.003 y_break_pos+0.01+0.003], 'LineWidth', 0.75);


% --- 3. Create the Bottom Heatmap Plots inside the Bottom Panel ---
bottomLayout = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
bottomLayout.Layout.Tile = 2;
bottomLayout.Layout.TileSpan = [2, 1];

% --- Define colors and column names for the heatmaps ---
colNames = {'env', 'err (env)', 'pss', 'err (pss)', 'cneg', 'pneg', '8-state', '13-state'};
mypurple = [0.4940, 0.1840, 0.5560];
mygold   = [0.9290, 0.6940, 0.1250];
white    = [1.0000, 1.0000, 1.0000];
cmap = [white; myblue; myorange; mygold; mypurple];
row_splits = [1, 22, 43, length(names)+1];

% --- Define rows and colors for background highlighting ---
% For '8-state' column (column 7)
gray_rows_8 = [1, 2, 3, 7, 9, 10, 13, 23]; 
blue_rows_8 = [5, 15, 19, 21, 26];       

% For '13-state' column (column 8)
gray_rows_13 = [1, 2, 3, 5, 6, 8, 11, 12, 13, 14, 17, 18, 20];
blue_rows_13 = [4, 16, 22, 24];

% Calculate tinted colors
light_gray = [0.8, 0.8, 0.8]; 
light_gold = 0.2 * mygold + 0.8 * white;
light_blue = 0.2 * myblue + 0.8 * white;

% Loop to create the three heatmaps
for i = 1:3
    row_indices = row_splits(i) : row_splits(i+1)-1;
    Xdat_subset = Xdat(row_indices, :);
    names_subset = names(row_indices);
    [num_rows, num_data_cols] = size(Xdat_subset);
    num_total_cols = num_data_cols + 2; % 6 data columns + 2 model columns

    color_index_matrix = ones(num_rows, num_data_cols);
    color_index_matrix(Xdat_subset(:, 1) < 0.1, 1) = 2;
    color_index_matrix(Xdat_subset(:, 2) < 0.1, 2) = 2;
    color_index_matrix(Xdat_subset(:, 3) < 0.1, 3) = 3;
    color_index_matrix(Xdat_subset(:, 4) < 0.1, 4) = 3;
    color_index_matrix(Xdat_subset(:, 5) < 0.1, 5) = 4;
    color_index_matrix(Xdat_subset(:, 6) < 0.1, 6) = 5;

    ax_heat = nexttile(bottomLayout);
    
    xlim(ax_heat, [0.5, num_total_cols + 0.5]);
    ylim(ax_heat, [0.5, num_rows + 0.5]);
    set(ax_heat, 'YDir', 'reverse');

    hold(ax_heat, 'on');

    % --- MODIFIED CODE: Separate logic for column 7 and column 8 --- %%
    all_local_rows = 1:num_rows;

    % --- Coloring for Column 7 ('8-state') ---
    local_blue_rows_8 = find(ismember(row_indices, blue_rows_8));
    local_gray_rows_8 = find(ismember(row_indices, gray_rows_8));
    local_gold_rows_8 = setdiff(all_local_rows, [local_blue_rows_8(:)' local_gray_rows_8(:)']);
    
    for k = 1:length(local_gold_rows_8); patch(ax_heat, [6.5, 7.5, 7.5, 6.5], [local_gold_rows_8(k) - 0.5, local_gold_rows_8(k) - 0.5, local_gold_rows_8(k) + 0.5, local_gold_rows_8(k) + 0.5], light_gold, 'EdgeColor', 'none'); end
    for k = 1:length(local_blue_rows_8); patch(ax_heat, [6.5, 7.5, 7.5, 6.5], [local_blue_rows_8(k) - 0.5, local_blue_rows_8(k) - 0.5, local_blue_rows_8(k) + 0.5, local_blue_rows_8(k) + 0.5], light_blue, 'EdgeColor', 'none'); end
    for k = 1:length(local_gray_rows_8); patch(ax_heat, [6.5, 7.5, 7.5, 6.5], [local_gray_rows_8(k) - 0.5, local_gray_rows_8(k) - 0.5, local_gray_rows_8(k) + 0.5, local_gray_rows_8(k) + 0.5], light_gray, 'EdgeColor', 'none'); end

    % --- Coloring for Column 8 ('13-state') ---
    local_blue_rows_13 = find(ismember(row_indices, blue_rows_13));
    local_gray_rows_13 = find(ismember(row_indices, gray_rows_13));
    local_gold_rows_13 = setdiff(all_local_rows, [local_blue_rows_13(:)' local_gray_rows_13(:)']);

    for k = 1:length(local_gold_rows_13); patch(ax_heat, [7.5, 8.5, 8.5, 7.5], [local_gold_rows_13(k) - 0.5, local_gold_rows_13(k) - 0.5, local_gold_rows_13(k) + 0.5, local_gold_rows_13(k) + 0.5], light_gold, 'EdgeColor', 'none'); end
    for k = 1:length(local_blue_rows_13); patch(ax_heat, [7.5, 8.5, 8.5, 7.5], [local_blue_rows_13(k) - 0.5, local_blue_rows_13(k) - 0.5, local_blue_rows_13(k) + 0.5, local_blue_rows_13(k) + 0.5], light_blue, 'EdgeColor', 'none'); end
    for k = 1:length(local_gray_rows_13); patch(ax_heat, [7.5, 8.5, 8.5, 7.5], [local_gray_rows_13(k) - 0.5, local_gray_rows_13(k) - 0.5, local_gray_rows_13(k) + 0.5, local_gray_rows_13(k) + 0.5], light_gray, 'EdgeColor', 'none'); end
    % --- END OF MODIFIED CODE --- %%

    [rows, cols] = find(color_index_matrix ~= 1);
    
    for k = 1:length(rows)
        r = rows(k);
        c = cols(k);
        color_idx = color_index_matrix(r, c);
        checkmark_color = cmap(color_idx, :);
        text(ax_heat, c, r, '✓', ...
            'Color', checkmark_color, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 12, ...
            'FontWeight', 'bold');
    end
    
    hold(ax_heat, 'off');
    
    ax_heat.XTick = 1:num_total_cols;
    ax_heat.XTickLabel = colNames;
    ax_heat.XTickLabelRotation = 45;
    
    ax_heat.YTick = 1:num_rows;
    ax_heat.YTickLabel = names_subset;
    
    ax_heat.TickLength = [0 0];
    x_grid_lines = 0.5 : 1 : num_total_cols + 0.5;
    y_grid_lines = 0.5 : 1 : num_rows + 0.5;
    ax_heat.XAxis.MinorTickValues = x_grid_lines;
    ax_heat.YAxis.MinorTickValues = y_grid_lines;
    grid(ax_heat, 'off');
    grid(ax_heat);
    ax_heat.GridAlpha = 0;
    ax_heat.XMinorGrid = 'on';
    ax_heat.YMinorGrid = 'on';
    ax_heat.MinorGridColor = [0 0 0];
    ax_heat.MinorGridAlpha = 1;
    ax_heat.MinorGridLineStyle = '-'; 
    ax_heat.Layer = 'top';

    if i == 1
        text(ax_heat, -0.33, 1.07, "B) Analysis of reduction approaches", ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');
    end
end

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_analysis_combined.pdf", 'ContentType', 'vector');

%% sum of ir-indices

% Load data first, as it's needed for both plot sections
load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")

% --- 1. Setup Figure and Panels ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 5]);
hold on;
plot(model.t_ref, sum(model.ir.index, 2), 'LineWidth', 1.5, 'LineStyle', '--');
load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
plot(model.t_ref, sum(ir.index, 2), 'LineWidth', 1.5);
load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
plot(model.t_ref, sum(ir.index, 2), 'LineWidth', 1.5);
% xlim([-0.2, 10])
xlim([-1, 41])
ylim([-0.1e4, 2.5e4])
xlabel("t [h]")
ylabel("sum of ir-indices")
% title_text_maxnir = "Sum of ir-indices";
legend(["original", "8-state", "13-state"])
grid("on")
box on;
% text(-0.17, 1.27, title_text_maxnir, ...
%                 'Units', 'normalized', ...
%                 'HorizontalAlignment', 'left', ...
%                 'VerticalAlignment', 'top', ...
%                 'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');


% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_sum_ir_indices.pdf", 'ContentType', 'vector');

% Load data first, as it's needed for both plot sections
load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")

% --- 1. Setup Figure and Panels ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 5]);
hold on;
plot(model.t_ref, sum(model.ir.index, 2), 'LineWidth', 1.5, 'LineStyle', '--');
load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
plot(model.t_ref, sum(ir.index, 2), 'LineWidth', 1.5);
load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
plot(model.t_ref, sum(ir.index, 2), 'LineWidth', 1.5);
% xlim([-0.2, 10])
% xlim([-1, 41])
xlim([-0.005 0.155]);
ylim([-0.1e4, 2.5e4])
xlabel("t [h]")
ylabel("sum of ir-indices")
% title_text_maxnir = "Sum of ir-indices";
legend(["original", "8-state", "13-state"])
grid("on")
box on;
% text(-0.17, 1.27, title_text_maxnir, ...
%                 'Units', 'normalized', ...
%                 'HorizontalAlignment', 'left', ...
%                 'VerticalAlignment', 'top', ...
%                 'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');


% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_sum_ir_indices_015.pdf", 'ContentType', 'vector');

%% Fibrinogen timecourse

% Load data first, as it's needed for both plot sections
load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
default_colors = get(groot, 'DefaultAxesColorOrder');

% --- 1. Setup Figure and Panels ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 12, 7]);
hold on;
idx = model.I.Fg;
plot(model.t_ref, model.X_ref(:, idx), 'LineWidth', 1.5, 'Color', default_colors(2, :));
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
plot(model.t_ref, redmodel.X_red(:, idx), 'LineWidth', 1.5, 'LineStyle', '--', 'Color', default_colors(2, :));
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
plot(model.t_ref, redmodel.X_red(:, idx), 'LineWidth', 1.5, 'LineStyle', '-.', 'Color', default_colors(2, :));
% xlim([-0.2, 10])
xlim([-1, 41])
ylim([-300, 10000])
% ylim([-0.1e4, 2.5e4])
xlabel("t [h]")
ylabel("Fibrinogen concentration [nM]")
% title_text_maxnir = "Sum of ir-indices";
legend(["original", "8-state", "13-state"], 'Location', 'eastoutside')
grid("on")
% yscale('log')
box on;
% text(-0.17, 1.27, title_text_maxnir, ...
%                 'Units', 'normalized', ...
%                 'HorizontalAlignment', 'left', ...
%                 'VerticalAlignment', 'top', ...
%                 'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');


% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_Fg_timecourses.pdf", 'ContentType', 'vector');

%% Timecourses and ir-indices --- original model


% SINGLE PLOT --- timecourses
default_colors = get(groot, 'DefaultAxesColorOrder');
new_colors = [
    0.90, 0.10, 0.15; 0.60, 0.90, 0.10; 0.00, 0.60, 0.50;
    0.95, 0.30, 0.70; 0.60, 0.40, 0.20; 0.20, 0.20, 0.60;
];
colors = [default_colors; new_colors];
colors = colors(1:11, :);

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
names = model.I.nmstate;
names = replace(names, '_', ':');

% --- 1. Setup Figure and Tiled Layout ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 12, 7]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- LEFT TILE: Original .index data, NO Legend ---
ax1 = gca;
colororder(ax1, colors);
hold(ax1, 'on');
idx = model.analysis.ir.I_sorted_max_nindex([1:10, 13]);
plot(ax1, model.t_ref, model.X_ref(:, idx), 'LineWidth', 1.5);
% plot(ax1, model.t_ref, model.ir.index(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax1, [-1, 41])
ylim(ax1, [1e-6, 1e5])
xlabel(ax1, "t [h]")
ylabel(ax1, "concentration [nM]")
legend(ax1, names{idx}, 'Location', 'eastoutside')
yscale(ax1, 'log')
grid(ax1, "on")
box(ax1, 'on');
% Title removed

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_timecourses_original_model.pdf", 'ContentType', 'vector');


% SINGLE INSET PLOT
default_colors = get(groot, 'DefaultAxesColorOrder');
new_colors = [
    0.90, 0.10, 0.15; 0.60, 0.90, 0.10; 0.00, 0.60, 0.50;
    0.95, 0.30, 0.70; 0.60, 0.40, 0.20; 0.20, 0.20, 0.60;
];
colors = [default_colors; new_colors];
colors = colors(1:11, :);

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
names = model.I.nmstate;
names = replace(names, '_', ':');

% --- 1. Setup Figure and Tiled Layout ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 12, 7]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- LEFT TILE: Original .index data, NO Legend ---
ax1 = gca;
colororder(ax1, colors);
hold(ax1, 'on');
% plot(ax1, model.t_ref, ir.index(:, idx), 'LineWidth', 1.5);
idx = model.analysis.ir.I_sorted_max_nindex([1:10, 13]);
plot(ax1, model.t_ref, model.ir.index(:, idx), 'LineWidth', 1.5);
xlim(ax1, [-1, 41])
ylim(ax1, [-300, 10000])
xlabel(ax1, "t [h]")
ylabel(ax1, "ir-index")
legend(ax1, names{idx}, 'Location', 'eastoutside')
grid(ax1, "on")
box(ax1, 'on');
% Title removed

inset_axes = axes('Position', [0.42, 0.5, 0.25, 0.4]);
colororder(inset_axes, colors);
hold(inset_axes, 'on');
idx = model.analysis.ir.I_sorted_max_nindex([1:10, 13]);
% plot(inset_axes, model.t_ref, ir.index(:, idx), 'LineWidth', 1.5);
plot(inset_axes, model.t_ref, model.ir.index(:, idx), 'LineWidth', 1.5);
xlim(inset_axes, [-0.005 0.155])
ylim(inset_axes, [-300, 10000])
% xlabel(inset_axes, "t [h]")
% ylabel(inset_axes, "ir-index")
grid(inset_axes, "on")
box(inset_axes, 'on');

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_ir_indices_original_model_inset.pdf", 'ContentType', 'vector');

%% Individual ir-indices --- 8-state model


% SINGLE INSET PLOT
default_colors = get(groot, 'DefaultAxesColorOrder');
new_colors = [
    0.90, 0.10, 0.15; 0.60, 0.90, 0.10; 0.00, 0.60, 0.50;
    0.95, 0.30, 0.70; 0.60, 0.40, 0.20; 0.20, 0.20, 0.60;
];
colors = [default_colors; new_colors];
colors = colors(1:8, :);

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
names = model.I.nmstate;
names = replace(names, '_', ':');

% --- 1. Setup Figure and Tiled Layout ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 12, 7]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- LEFT TILE: Original .index data, NO Legend ---
ax1 = gca;
colororder(ax1, colors);
hold(ax1, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
plot(ax1, model.t_ref, ir.index(:, idx), 'LineWidth', 1.5);
plot(ax1, model.t_ref, model.ir.index(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax1, [-1, 41])
ylim(ax1, [-300, 10000])
xlabel(ax1, "t [h]")
ylabel(ax1, "ir-index")
legend(ax1, names{idx}, 'Location', 'eastoutside')
grid(ax1, "on")
box(ax1, 'on');
% Title removed

inset_axes = axes('Position', [0.42, 0.5, 0.25, 0.4]);
colororder(inset_axes, colors);
hold(inset_axes, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
plot(inset_axes, model.t_ref, ir.index(:, idx), 'LineWidth', 1.5);
plot(inset_axes, model.t_ref, model.ir.index(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(inset_axes, [-0.005 0.155])
ylim(inset_axes, [-300, 10000])
% xlabel(inset_axes, "t [h]")
% ylabel(inset_axes, "ir-index")
grid(inset_axes, "on")
box(inset_axes, 'on');

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_ir_indices_8_state_reduced_model_inset.pdf", 'ContentType', 'vector');


% --- Setup code (kept from original) ---
default_colors = get(groot, 'DefaultAxesColorOrder');
new_colors = [
    0.90, 0.10, 0.15; 0.60, 0.90, 0.10; 0.00, 0.60, 0.50;
    0.95, 0.30, 0.70; 0.60, 0.40, 0.20; 0.20, 0.20, 0.60;
];
colors = [default_colors; new_colors];
colors = colors(1:8, :);

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
names = model.I.nmstate;
names = replace(names, '_', ':');

% --- 1. Setup Figure and Tiled Layout ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 8]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- LEFT TILE: Original .index data, NO Legend ---
ax1 = nexttile;
colororder(ax1, colors);
hold(ax1, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
plot(ax1, model.t_ref, ir.index(:, idx), 'LineWidth', 1.5);
plot(ax1, model.t_ref, model.ir.index(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax1, [-1, 41])
ylim(ax1, [-300, 10000])
xlabel(ax1, "t [h]")
ylabel(ax1, "ir-index")
grid(ax1, "on")
box(ax1, 'on');
% Title removed

% --- RIGHT TILE: Normalized .nindex data, Legend Outside ---
ax2 = nexttile;
colororder(ax2, colors);
hold(ax2, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
plot(ax2, model.t_ref, ir.nindex(:, idx), 'LineWidth', 1.5);
plot(ax2, model.t_ref, model.ir.nindex(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax2, [-1, 41])
ylim(ax2, [-0.035, 1.035])
xlabel(ax2, "t [h]")
% --- CHANGE: Y-label updated ---
ylabel(ax2, "nir-index")
legend(ax2, names{idx}, 'Location', 'eastoutside')
grid(ax2, "on")
box(ax2, 'on');
% Title removed

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_ir_indices_8_state_reduced_model.pdf", 'ContentType', 'vector');


% Script 2 Final Version: No Titles, Updated Y-Label (Zoomed)

% --- Setup code (kept from original) ---
default_colors = get(groot, 'DefaultAxesColorOrder');
new_colors = [
    0.90, 0.10, 0.15; 0.60, 0.90, 0.10; 0.00, 0.60, 0.50;
    0.95, 0.30, 0.70; 0.60, 0.40, 0.20; 0.20, 0.20, 0.60;
];
colors = [default_colors; new_colors];
colors = colors(1:8, :);

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
names = model.I.nmstate;
names = replace(names, '_', ':');

% --- 1. Setup Figure and Tiled Layout ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 8]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- LEFT TILE: Original .index data, NO Legend ---
ax1 = nexttile;
colororder(ax1, colors);
hold(ax1, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
plot(ax1, model.t_ref, ir.index(:, idx), 'LineWidth', 1.5);
plot(ax1, model.t_ref, model.ir.index(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax1, [-0.005 0.155]);
ylim(ax1, [-300, 10000])
xlabel(ax1, "t [h]")
ylabel(ax1, "ir-index")
grid(ax1, "on")
box(ax1, 'on');
% Title removed

% --- RIGHT TILE: Normalized .nindex data, Legend Outside ---
ax2 = nexttile;
colororder(ax2, colors);
hold(ax2, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
plot(ax2, model.t_ref, ir.nindex(:, idx), 'LineWidth', 1.5);
plot(ax2, model.t_ref, model.ir.nindex(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax2, [-0.005 0.155]);
ylim(ax2, [-0.035, 1.035])
xlabel(ax2, "t [h]")
% --- CHANGE: Y-label updated ---
ylabel(ax2, "nir-index")
legend(ax2, names{idx}, 'Location', 'eastoutside')
grid(ax2, "on")
box(ax2, 'on');
% Title removed

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_ir_indices_8_state_reduced_model_015.pdf", 'ContentType', 'vector');

% 
% % Load data first, as it's needed for both plot sections
% load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
% 
% % --- 1. Setup Figure and Panels ---
% f = figure('Units', 'centimeters', 'Position', [0, 0, 12, 7]);
% ax = gca;
% colors = colors(1:8, :);
% colororder(ax, colors);
% hold on;
% load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
% idx = redmodel.redconfig == "dyn";% | redmodel.redconfig == "env";
% numbers = 1:redmodel.I.nstates;
% idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
% load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
% plot(model.t_ref, ir.index(:, idx), 'LineWidth', 1.5);
% load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
% plot(model.t_ref, model.ir.index(:, idx), 'LineWidth', 1, 'LineStyle', '--');
% % xlim([-0.2, 10])
% % xlim([-1, 41])
% xlim([-0.005 0.155]);
% ylim([-300, 10000])
% xlabel("t [h]")
% ylabel("ir-index")
% % title_text_maxnir = "Sum of ir-indices";
% legend(names{idx}, 'Location', 'eastoutside')
% grid("on")
% box on;
% % text(-0.17, 1.27, title_text_maxnir, ...
% %                 'Units', 'normalized', ...
% %                 'HorizontalAlignment', 'left', ...
% %                 'VerticalAlignment', 'top', ...
% %                 'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');
% 
% 
% % --- 4. Export the Figure ---
% drawnow;
% exportgraphics(f, "./figures/BC_ir_indices_8_state_reduced_model_015_quadratic.pdf", 'ContentType', 'vector');

%% 2 types of normalized indices --- 8-state model

% --- Setup code (kept from original) ---
default_colors = get(groot, 'DefaultAxesColorOrder');
new_colors = [
    0.90, 0.10, 0.15; 0.60, 0.90, 0.10; 0.00, 0.60, 0.50;
    0.95, 0.30, 0.70; 0.60, 0.40, 0.20; 0.20, 0.20, 0.60;
];
colors = [default_colors; new_colors];
colors = colors(1:8, :);

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
names = model.I.nmstate;
names = replace(names, '_', ':');

% --- 1. Setup Figure and Tiled Layout ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 8]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- LEFT TILE: Original .index data, NO Legend ---
ax1 = nexttile;
colororder(ax1, colors);
hold(ax1, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
plot(ax1, model.t_ref, ir.nindex(:, idx), 'LineWidth', 1.5);
plot(ax1, model.t_ref, model.ir.nindex(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax1, [-1, 41])
ylim(ax1, [-0.035, 1.035])
xlabel(ax1, "t [h]")
ylabel(ax1, "nir-index")
grid(ax1, "on")
box(ax1, 'on');
title(ax1, "nir-indices (normalized by all states)")
% Title removed

% --- RIGHT TILE: Normalized .nindex data, Legend Outside ---
ax2 = nexttile;
colororder(ax2, colors);
hold(ax2, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
plot(ax2, model.t_ref, ir.index(:, idx) ./ sum(ir.index(:, idx), 2), 'LineWidth', 1.5);
plot(ax2, model.t_ref, model.ir.index(:, idx) ./ sum(model.ir.index(:, idx), 2), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax2, [-1, 41])
ylim(ax2, [-0.035, 1.035])
xlabel(ax2, "t [h]")
% --- CHANGE: Y-label updated ---
% ylabel(ax2, "nir-index")
legend(ax2, names{idx}, 'Location', 'eastoutside')
grid(ax2, "on")
box(ax2, 'on');
title(ax2, "nir-indices (normalized by remaining dynamic states)")
% Title removed

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_nir_indices_8_state_reduced_model.pdf", 'ContentType', 'vector');


% Script 2 Final Version: No Titles, Updated Y-Label (Zoomed)

% --- Setup code (kept from original) ---
default_colors = get(groot, 'DefaultAxesColorOrder');
new_colors = [
    0.90, 0.10, 0.15; 0.60, 0.90, 0.10; 0.00, 0.60, 0.50;
    0.95, 0.30, 0.70; 0.60, 0.40, 0.20; 0.20, 0.20, 0.60;
];
colors = [default_colors; new_colors];
colors = colors(1:8, :);

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
names = model.I.nmstate;
names = replace(names, '_', ':');

% --- 1. Setup Figure and Tiled Layout ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 8]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- LEFT TILE: Original .index data, NO Legend ---
ax1 = nexttile;
colororder(ax1, colors);
hold(ax1, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
plot(ax1, model.t_ref, ir.nindex(:, idx), 'LineWidth', 1.5);
plot(ax1, model.t_ref, model.ir.nindex(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax1, [-0.005 0.155]);
ylim(ax1, [-0.035, 1.035])
xlabel(ax1, "t [h]")
ylabel(ax1, "ir-index")
grid(ax1, "on")
box(ax1, 'on');
title(ax1, "nir-indices (normalized by all states)")
% Title removed

% --- RIGHT TILE: Normalized .nindex data, Legend Outside ---
ax2 = nexttile;
colororder(ax2, colors);
hold(ax2, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_8_state_full.mat")
plot(ax2, model.t_ref, ir.index(:, idx) ./ sum(ir.index(:, idx), 2), 'LineWidth', 1.5);
plot(ax2, model.t_ref, model.ir.index(:, idx) ./ sum(model.ir.index(:, idx), 2), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax2, [-0.005 0.155]);
ylim(ax2, [-0.035, 1.035])
xlabel(ax2, "t [h]")
% --- CHANGE: Y-label updated ---
ylabel(ax2, "nir-index")
legend(ax2, names{idx}, 'Location', 'eastoutside')
grid(ax2, "on")
box(ax2, 'on');
title(ax2, "nir-indices (normalized by remaining dynamic states)")
% Title removed

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_nir_indices_8_state_reduced_model_015.pdf", 'ContentType', 'vector');


%% Individual ir-indices --- 13-state model

% SINGLE INSET PLOT
default_colors = get(groot, 'DefaultAxesColorOrder');
new_colors = [
    0.90, 0.10, 0.15; 0.60, 0.90, 0.10; 0.00, 0.60, 0.50;
    0.95, 0.30, 0.70; 0.60, 0.40, 0.20; 0.20, 0.20, 0.60;
];
colors = [default_colors; new_colors];
% colors = colors(1:8, :);

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
names = model.I.nmstate;
names = replace(names, '_', ':');

% --- 1. Setup Figure and Tiled Layout ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 12, 7]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- LEFT TILE: Original .index data, NO Legend ---
ax1 = gca;
colororder(ax1, colors);
hold(ax1, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
plot(ax1, model.t_ref, ir.index(:, idx), 'LineWidth', 1.5);
plot(ax1, model.t_ref, model.ir.index(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax1, [-1, 41])
ylim(ax1, [-300, 10000])
xlabel(ax1, "t [h]")
ylabel(ax1, "ir-index")
legend(ax1, names{idx}, 'Location', 'eastoutside')
grid(ax1, "on")
box(ax1, 'on');
% Title removed

inset_axes = axes('Position', [0.42, 0.5, 0.25, 0.4]);
colororder(inset_axes, colors);
hold(inset_axes, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
plot(inset_axes, model.t_ref, ir.index(:, idx), 'LineWidth', 1.5);
plot(inset_axes, model.t_ref, model.ir.index(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(inset_axes, [-0.005 0.155])
ylim(inset_axes, [-300, 10000])
% xlabel(inset_axes, "t [h]")
% ylabel(inset_axes, "ir-index")
grid(inset_axes, "on")
box(inset_axes, 'on');

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_ir_indices_13_state_reduced_model_inset.pdf", 'ContentType', 'vector');


% --- Setup code (kept from original) ---
default_colors = get(groot, 'DefaultAxesColorOrder');
new_colors = [
    0.90, 0.10, 0.15; 0.60, 0.90, 0.10; 0.00, 0.60, 0.50;
    0.95, 0.30, 0.70; 0.60, 0.40, 0.20; 0.20, 0.20, 0.60;
];
colors = [default_colors; new_colors];
% colors = colors(1:8, :);

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
names = model.I.nmstate;
names = replace(names, '_', ':');

% --- 1. Setup Figure and Tiled Layout ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 8]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- LEFT TILE: Original .index data, NO Legend ---
ax1 = nexttile;
colororder(ax1, colors);
hold(ax1, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
plot(ax1, model.t_ref, ir.index(:, idx), 'LineWidth', 1.5);
plot(ax1, model.t_ref, model.ir.index(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax1, [-1, 41])
ylim(ax1, [-300, 10000])
xlabel(ax1, "t [h]")
ylabel(ax1, "ir-index")
grid(ax1, "on")
box(ax1, 'on');
% Title removed

% --- RIGHT TILE: Normalized .nindex data, Legend Outside ---
ax2 = nexttile;
colororder(ax2, colors);
hold(ax2, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
plot(ax2, model.t_ref, ir.nindex(:, idx), 'LineWidth', 1.5);
plot(ax2, model.t_ref, model.ir.nindex(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax2, [-1, 41])
ylim(ax2, [-0.035, 1.035])
xlabel(ax2, "t [h]")
% --- CHANGE: Y-label updated ---
ylabel(ax2, "nir-index")
legend(ax2, names{idx}, 'Location', 'eastoutside')
grid(ax2, "on")
box(ax2, 'on');
% Title removed

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_ir_indices_13_state_reduced_model.pdf", 'ContentType', 'vector');


% Script 2 Final Version: No Titles, Updated Y-Label (Zoomed)

% --- Setup code (kept from original) ---
default_colors = get(groot, 'DefaultAxesColorOrder');
new_colors = [
    0.90, 0.10, 0.15; 0.60, 0.90, 0.10; 0.00, 0.60, 0.50;
    0.95, 0.30, 0.70; 0.60, 0.40, 0.20; 0.20, 0.20, 0.60;
];
colors = [default_colors; new_colors];
% colors = colors(1:8, :);

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
names = model.I.nmstate;
names = replace(names, '_', ':');

% --- 1. Setup Figure and Tiled Layout ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 8]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- LEFT TILE: Original .index data, NO Legend ---
ax1 = nexttile;
colororder(ax1, colors);
hold(ax1, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
plot(ax1, model.t_ref, ir.index(:, idx), 'LineWidth', 1.5);
plot(ax1, model.t_ref, model.ir.index(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax1, [-0.005 0.155]);
ylim(ax1, [-300, 10000])
xlabel(ax1, "t [h]")
ylabel(ax1, "ir-index")
grid(ax1, "on")
box(ax1, 'on');
% Title removed

% --- RIGHT TILE: Normalized .nindex data, Legend Outside ---
ax2 = nexttile;
colororder(ax2, colors);
hold(ax2, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
plot(ax2, model.t_ref, ir.nindex(:, idx), 'LineWidth', 1.5);
plot(ax2, model.t_ref, model.ir.nindex(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax2, [-0.005 0.155]);
ylim(ax2, [-0.035, 1.035])
xlabel(ax2, "t [h]")
% --- CHANGE: Y-label updated ---
ylabel(ax2, "nir-index")
legend(ax2, names{idx}, 'Location', 'eastoutside')
grid(ax2, "on")
box(ax2, 'on');
% Title removed

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_ir_indices_13_state_reduced_model_015.pdf", 'ContentType', 'vector');

% 
% % Load data first, as it's needed for both plot sections
% load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
% 
% % --- 1. Setup Figure and Panels ---
% f = figure('Units', 'centimeters', 'Position', [0, 0, 12, 7]);
% ax = gca;
% colors = colors(1:13, :);
% colororder(ax, colors);
% hold on;
% load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
% idx = redmodel.redconfig == "dyn";% | redmodel.redconfig == "env";
% numbers = 1:redmodel.I.nstates;
% idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
% load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
% plot(model.t_ref, ir.index(:, idx), 'LineWidth', 1.5);
% load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
% plot(model.t_ref, model.ir.index(:, idx), 'LineWidth', 1, 'LineStyle', '--');
% % xlim([-0.2, 10])
% % xlim([-1, 41])
% xlim([-0.005 0.155]);
% ylim([-300, 10000])
% xlabel("t [h]")
% ylabel("ir-index")
% % title_text_maxnir = "Sum of ir-indices";
% leg = legend(names{idx}, 'Location', 'eastoutside');
% set(leg.EntryContainer.NodeChildren, 'LineStyle', '-');
% grid("on")
% box on;
% hold off;
% % text(-0.17, 1.27, title_text_maxnir, ...
% %                 'Units', 'normalized', ...
% %                 'HorizontalAlignment', 'left', ...
% %                 'VerticalAlignment', 'top', ...
% %                 'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');
% 
% 
% % --- 4. Export the Figure ---
% drawnow;
% exportgraphics(f, "./figures/BC_ir_indices_13_state_reduced_model_015_quadratic.pdf", 'ContentType', 'vector');

%% 2 types of normalized indices --- 8-state model

% --- Setup code (kept from original) ---
default_colors = get(groot, 'DefaultAxesColorOrder');
new_colors = [
    0.90, 0.10, 0.15; 0.60, 0.90, 0.10; 0.00, 0.60, 0.50;
    0.95, 0.30, 0.70; 0.60, 0.40, 0.20; 0.20, 0.20, 0.60;
];
colors = [default_colors; new_colors];
% colors = colors(1:8, :);

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
names = model.I.nmstate;
names = replace(names, '_', ':');

% --- 1. Setup Figure and Tiled Layout ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 8]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- LEFT TILE: Original .index data, NO Legend ---
ax1 = nexttile;
colororder(ax1, colors);
hold(ax1, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
plot(ax1, model.t_ref, ir.nindex(:, idx), 'LineWidth', 1.5);
plot(ax1, model.t_ref, model.ir.nindex(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax1, [-1, 41])
ylim(ax1, [-0.035, 1.035])
xlabel(ax1, "t [h]")
ylabel(ax1, "nir-index")
grid(ax1, "on")
box(ax1, 'on');
title(ax1, "nir-indices (normalized by all states)")
% Title removed

% --- RIGHT TILE: Normalized .nindex data, Legend Outside ---
ax2 = nexttile;
colororder(ax2, colors);
hold(ax2, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
plot(ax2, model.t_ref, ir.index(:, idx) ./ sum(ir.index(:, idx), 2), 'LineWidth', 1.5);
plot(ax2, model.t_ref, model.ir.index(:, idx) ./ sum(model.ir.index(:, idx), 2), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax2, [-1, 41])
ylim(ax2, [-0.035, 1.035])
xlabel(ax2, "t [h]")
% --- CHANGE: Y-label updated ---
% ylabel(ax2, "nir-index")
legend(ax2, names{idx}, 'Location', 'eastoutside')
grid(ax2, "on")
box(ax2, 'on');
title(ax2, "nir-indices (normalized by remaining dynamic states)")
% Title removed

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_nir_indices_13_state_reduced_model.pdf", 'ContentType', 'vector');


% Script 2 Final Version: No Titles, Updated Y-Label (Zoomed)

% --- Setup code (kept from original) ---
default_colors = get(groot, 'DefaultAxesColorOrder');
new_colors = [
    0.90, 0.10, 0.15; 0.60, 0.90, 0.10; 0.00, 0.60, 0.50;
    0.95, 0.30, 0.70; 0.60, 0.40, 0.20; 0.20, 0.20, 0.60;
];
colors = [default_colors; new_colors];
% colors = colors(1:8, :);

load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
names = model.I.nmstate;
names = replace(names, '_', ':');

% --- 1. Setup Figure and Tiled Layout ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 8]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- LEFT TILE: Original .index data, NO Legend ---
ax1 = nexttile;
colororder(ax1, colors);
hold(ax1, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
plot(ax1, model.t_ref, ir.nindex(:, idx), 'LineWidth', 1.5);
plot(ax1, model.t_ref, model.ir.nindex(:, idx), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax1, [-0.005 0.155]);
ylim(ax1, [-0.035, 1.035])
xlabel(ax1, "t [h]")
ylabel(ax1, "ir-index")
grid(ax1, "on")
box(ax1, 'on');
title(ax1, "nir-indices (normalized by all states)")
% Title removed

% --- RIGHT TILE: Normalized .nindex data, Legend Outside ---
ax2 = nexttile;
colororder(ax2, colors);
hold(ax2, 'on');
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
idx = redmodel.redconfig == "dyn";
numbers = 1:redmodel.I.nstates;
idx = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, numbers(idx)));
load("modelBC_SV40_from_JKn_2024_indices_13_state_full.mat")
plot(ax2, model.t_ref, ir.index(:, idx) ./ sum(ir.index(:, idx), 2), 'LineWidth', 1.5);
plot(ax2, model.t_ref, model.ir.index(:, idx) ./ sum(model.ir.index(:, idx), 2), 'LineWidth', 1, 'LineStyle', '--');
xlim(ax2, [-0.005 0.155]);
ylim(ax2, [-0.035, 1.035])
xlabel(ax2, "t [h]")
% --- CHANGE: Y-label updated ---
ylabel(ax2, "nir-index")
legend(ax2, names{idx}, 'Location', 'eastoutside')
grid(ax2, "on")
box(ax2, 'on');
title(ax2, "nir-indices (normalized by remaining dynamic states)")
% Title removed

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_nir_indices_13_state_reduced_model_015.pdf", 'ContentType', 'vector');


%% approximation quality in 8-state reduced model

% Load data first, as it's needed for both plot sections
load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")

% --- Prepare Data for All Plots ---
names = model.I.nmstate;
names = replace(names, '_', ':');
id_names = 1:model.I.nstates;
id_names = id_names(redmodel.redconfig == "env");
for id = id_names
    names{id} = [names{id}, ' (env)'];
end
% names{46} = 'VKp';
% names{54} = 'ENOp';
% names{60} = 'AVenom\_Tiger';
% names{61} = 'CVenom\_Tiger';

idx_in_redmodel = 1:model.I.nstates;
idx_in_redmodel = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, idx_in_redmodel(redmodel.redconfig == "dyn" | redmodel.redconfig == "env")));

% --- 1. Setup Figure and Panels ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 16]);

mainLayout = tiledlayout(f, 4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(idx_in_redmodel)
    idx = idx_in_redmodel(i);
    ax = nexttile(mainLayout);
    hold on;
    plot(redmodel.t_ref, redmodel.X_ref(:, idx), 'LineWidth', 1.5, 'LineStyle', '--');
    plot(redmodel.t_ref, redmodel.X_red(:, idx), 'LineWidth', 1.5);
    % xlim([-0.2, 10])
    xlim([-1, 41])
    ylim([1e-8 1e5]);
    yticks([1e-5, 1e0, 1e5]);
    % ylim([-0.1e4, 2.5e4])
    if i>3*4
        xlabel("t [h]")
    end
    if mod(i, 4) == 1
        ylabel("concentration [nM]")
    end
    % title_text_maxnir = "Sum of ir-indices";
    % legend(["original model", "8-state reduced model"])
    title(names{idx})
    grid("on")
    set(ax, 'YScale', 'log');
    box on

    % Get the current axis limits
    ax_lims_x = xlim;
    ax_lims_y = ylim;
    
    % Format the text string using the value from your error vector
    % '%.2e' formats the number in scientific notation with 2 decimal places
    error_text = sprintf('Error: %.3g%%', 100*redmodel.redobj.err(idx)); 
    % error_text = sprintf(['Error: ', char(num2str(100*redmodel.redobj.err(idx))), '\%']); 
    
    % Add the text to the plot
    % We place it at the bottom-right corner and use alignment properties
    % to ensure it sits neatly inside the plot box.
    text(ax_lims_x(2), ax_lims_y(1), error_text, ...
         'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'bottom', ...
         'FontSize', 8); % Optional: add a white background for readability
end

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_approximation_8_state_reduced_model.pdf", 'ContentType', 'vector');

%% approximation quality in 13-state reduced model

% Load data first, as it's needed for both plot sections
load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")

% --- Prepare Data for All Plots ---
names = model.I.nmstate;
names = replace(names, '_', ':');
id_names = 1:model.I.nstates;
id_names = id_names(redmodel.redconfig == "env");
for id = id_names
    names{id} = [names{id}, ' (env)'];
end
% names{46} = 'VKp';
% names{54} = 'ENOp';
% names{60} = 'AVenom\_Tiger';
% names{61} = 'CVenom\_Tiger';

idx_in_redmodel = 1:model.I.nstates;
idx_in_redmodel = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, idx_in_redmodel(redmodel.redconfig == "dyn" | redmodel.redconfig == "env")));

% --- 1. Setup Figure and Panels ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 35]);

mainLayout = tiledlayout(f, 5, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(idx_in_redmodel)
    idx = idx_in_redmodel(i);
    ax = nexttile(mainLayout);
    hold on;
    plot(redmodel.t_ref, redmodel.X_ref(:, idx), 'LineWidth', 1.5, 'LineStyle', '--');
    plot(redmodel.t_ref, redmodel.X_red(:, idx), 'LineWidth', 1.5);
    % xlim([-0.2, 10])
    xlim([-1, 41])
    ylim([1e-8 1e5]);
    yticks([1e-5, 1e0, 1e5]);
    % ylim([-0.1e4, 2.5e4])
    if i>4*4
        xlabel("t [h]")
    end
    if mod(i, 4) == 1
        ylabel("concentration [nM]")
    end
    % title_text_maxnir = "Sum of ir-indices";
    % legend(["original model", "8-state reduced model"])
    title(names{idx})
    grid("on")
    set(ax, 'YScale', 'log');
    box on

    % Get the current axis limits
    ax_lims_x = xlim;
    ax_lims_y = ylim;
    
    % Format the text string using the value from your error vector
    % '%.2e' formats the number in scientific notation with 2 decimal places
    error_text = sprintf('Error: %.3g%%', 100*redmodel.redobj.err(idx)); 
    % error_text = sprintf(['Error: ', char(num2str(100*redmodel.redobj.err(idx))), '\%']); 
    
    % Add the text to the plot
    % We place it at the bottom-right corner and use alignment properties
    % to ensure it sits neatly inside the plot box.
    text(ax_lims_x(2), ax_lims_y(1), error_text, ...
         'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'bottom', ...
         'FontSize', 8); % Optional: add a white background for readability
end

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_approximation_13_state_reduced_model.pdf", 'ContentType', 'vector');

%% approximation quality in 25-state reduced model

% Load data first, as it's needed for both plot sections
load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
% load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")

I = model.I;
% redconfig = repmat("pneg", [1 model.I.nstates]);
% redconfig([I.AVenom, I.CVenom, I.Fg, I.F, I.II, I.IIa, I.V, I.Va, I.VIII, I.VIIIa, I.IX, I.IXa, I.IXa_VIIIa, I.Xa, I.Xa_Va, ...
% I.XI, I.XIa, I.XIII, I.XIIIa, I.Tmod, I.IIa_Tmod, I.APC, I.APC_PS, I.Pg, I.P]) = "dyn";
% redconfig([I.X, I.PS, I.PC, I.TFPI, I.VKH2]) = "env";

load("modelBC_SV40_redvar_CV40_popprct95.mat")
redconfig = redmodel_seq_var.redconfig;

[redobj, ~, tred, Xred, err] = objfun_simple(model, redconfig, "rel2NE");

% --- Prepare Data for All Plots ---
names = model.I.nmstate;
names = replace(names, '_', ':');
for id_names = [I.PS, I.TFPI, I.VKH2]
    names{id_names} = [names{id_names}, ' (env)'];
end
% names{46} = 'VKp';
% names{54} = 'ENOp';
% names{60} = 'AVenom\_Tiger';
% names{61} = 'CVenom\_Tiger';

idx_in_redmodel = 1:model.I.nstates;
idx_in_redmodel = model.analysis.ir.I_sorted_max_nindex(ismember(model.analysis.ir.I_sorted_max_nindex, idx_in_redmodel(redconfig == "dyn" | redconfig == "env")));

% --- 1. Setup Figure and Panels ---
f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 35]);

mainLayout = tiledlayout(f, 6, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(idx_in_redmodel)
    idx = idx_in_redmodel(i);
    ax = nexttile(mainLayout);
    hold on;
    plot(model.t_ref, model.X_ref(:, idx), 'LineWidth', 1.5, 'LineStyle', '--');
    plot(tred, Xred(:, idx), 'LineWidth', 1.5);
    % xlim([-0.2, 10])
    xlim([-1, 41])
    ylim([1e-8 1e5]);
    yticks([1e-5, 1e0, 1e5]);
    % ylim([-0.1e4, 2.5e4])
    if i>5*5
        xlabel("t [h]")
    end
    if mod(i, 5) == 1
        ylabel("concentration [nM]")
    end
    % title_text_maxnir = "Sum of ir-indices";
    % legend(["original model", "8-state reduced model"])
    title(names{idx})
    grid("on")
    set(ax, 'YScale', 'log');
    box on

    % Get the current axis limits
    ax_lims_x = xlim;
    ax_lims_y = ylim;
    
    % Format the text string using the value from your error vector
    % '%.2e' formats the number in scientific notation with 2 decimal places
    error_text = sprintf('Error: %.3g%%', 100*redobj.err(idx)); 
    % error_text = sprintf(['Error: ', char(num2str(100*redmodel.redobj.err(idx))), '\%']); 
    
    % Add the text to the plot
    % We place it at the bottom-right corner and use alignment properties
    % to ensure it sits neatly inside the plot box.
    text(ax_lims_x(2), ax_lims_y(1), error_text, ...
         'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'bottom', ...
         'FontSize', 8); % Optional: add a white background for readability
end

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/BC_approximation_25_state_reduced_model.pdf", 'ContentType', 'vector');

%% Reduced models perturbation analysis

default_colors = get(groot, 'DefaultAxesColorOrder');

load("modelBC_SV40_from_JKn_2024.mat")

config = repmat("dyn", [1 model.I.nstates]);

% SIMULATION

% t_star = 1e-20;
% t_star = 0.01;
% t_star = 1/60;
% t_star = 5/60;
t_star = 10/60;
% t_star = 20/60;
% t_star = 30/60;
t_star = 60/60;
% t_star = 120/60;
t_star = 180/60;
t_star = 5;

% model.X0(model.I.AVenom) = 0.001*model.X0(model.I.AVenom);
Gulati = 1;
save_gulati = '';
if Gulati
    model.par(model.I.d_Brown) = 0.744;
    model.par(model.I.ka_Brown) = 0.854;
    save_gulati = '_gulati2014';
end

model_sim_init = model;

model_sim_init.t_ref = [0, t_star];
[t_1_full, X_1_full, ~, ~] = simModel_simple(model_sim_init, config);
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
[t_1_8_state, X_1_8_state, ~, ~] = simModel_simple(model_sim_init, redmodel.redconfig);
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
[t_1_13_state, X_1_13_state, ~, ~] = simModel_simple(model_sim_init, redmodel.redconfig);

t_2 = [t_star; model.t_ref(end);];
model_sim_init.t_ref = t_2;
model_sim_init.X0 = X_1_full(end, :)';
[t_2_full_base, X_2_full_base, ~, ~] = simModel_simple(model_sim_init, config);
load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
model_sim_init.X0 = X_1_8_state(end, :)';
[t_2_8_state_base, X_2_8_state_base, ~, ~] = simModel_simple(model_sim_init, redmodel.redconfig);
load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
model_sim_init.X0 = X_1_13_state(end, :)';
[t_2_13_state_base, X_2_13_state_base, ~, ~] = simModel_simple(model_sim_init, redmodel.redconfig);

model_sim = model;
model_sim.t_ref = t_2;
a = 0;
c = 0;
kb = 10;

X0_sim_full = X_1_full(end, :)';
% X0_sim_full(model.I.IIa) = 0;
X0_sim_full(model.I.AVenom) = a*X0_sim_full(model.I.AVenom);
X0_sim_full(model.I.CVenom) = c*X0_sim_full(model.I.CVenom);

X0_sim_8_state = X_1_8_state(end, :)';
% X0_sim_8_state(model.I.IIa) = 0;
X0_sim_8_state(model.I.AVenom) = a*X0_sim_8_state(model.I.AVenom);
X0_sim_8_state(model.I.CVenom) = c*X0_sim_8_state(model.I.CVenom);

X0_sim_13_state = X_1_13_state(end, :)';
% X0_sim_13_state(model.I.IIa) = 0;
X0_sim_13_state(model.I.AVenom) = a*X0_sim_13_state(model.I.AVenom);
X0_sim_13_state(model.I.CVenom) = c*X0_sim_13_state(model.I.CVenom);

model_sim.X0 = X0_sim_full;
[t_full_2, X_full_2, ~, ~] = simModel_simple(model_sim, config);

load("modelBC_SV40_from_JKn_2024_reduced_8_state.mat")
model_sim.X0 = X0_sim_8_state;
[t_8_state_2, X_8_state_2, ~, ~] = simModel_simple(model_sim, redmodel.redconfig);

load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")
model_sim.X0 = X0_sim_13_state;
[t_13_state_2, X_13_state_2, ~, ~] = simModel_simple(model_sim, redmodel.redconfig);

% PLOT

f = figure('Units', 'centimeters', 'Position', [0, 0, 20, 10]);
tlay = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
lw = 1.1;
if Gulati
    xlimits = [0 10];
else
    xlimits = [0 5];
end

ax = nexttile(tlay);
hold on;
plot(t_full_2, X_full_2(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(1, :))
plot(t_8_state_2, X_8_state_2(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(2, :))
plot(t_13_state_2, X_13_state_2(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(3, :))
plot(t_2_full_base, X_2_full_base(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(1, :), 'LineStyle', '--')
plot(t_2_8_state_base, X_2_8_state_base(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(2, :), 'LineStyle', '--')
plot(t_2_13_state_base, X_2_13_state_base(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(3, :), 'LineStyle', '--')
plot(t_1_full, X_1_full(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(1, :))
plot(t_1_8_state, X_1_8_state(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(2, :))
plot(t_1_13_state, X_1_13_state(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(3, :))
xline(t_star, '-k', 'LineWidth', lw, 'LineStyle', ':')
xlim([-1, 41])
ylim([-300, 10000])
ylabel("concentration [nM]")
% legend(["original", "8-state", "13-state"])
title("Fibrinogen")
box on
grid on

ax = nexttile(tlay);
hold on;
plot(t_full_2, X_full_2(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(1, :))
plot(t_8_state_2, X_8_state_2(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(2, :))
plot(t_13_state_2, X_13_state_2(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(3, :))
plot(t_2_full_base, X_2_full_base(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(1, :), 'LineStyle', '--')
plot(t_2_8_state_base, X_2_8_state_base(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(2, :), 'LineStyle', '--')
plot(t_2_13_state_base, X_2_13_state_base(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(3, :), 'LineStyle', '--')
plot(t_1_full, X_1_full(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(1, :))
plot(t_1_8_state, X_1_8_state(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(2, :))
plot(t_1_13_state, X_1_13_state(:, [model.I.Fg]), 'LineWidth', lw, 'Color', default_colors(3, :))
xline(t_star, '-k', 'LineWidth', lw, 'LineStyle', ':')
xlim(xlimits)
legend(["original", "8-state", "13-state", "original (no antivenom)", "8-state (no antivenom)", "13-state (no antivenom)"], 'Location', 'eastoutside')
title("Fibrinogen (log scale)")
yscale('log')
box on
grid on

ax = nexttile(tlay);
hold on;
plot(t_full_2, X_full_2(:, [model.I.AVenom]), 'LineWidth', lw, 'Color', default_colors(1, :))
% plot(t_8_state_2, X_8_state_2(:, [model.I.AVenom]), 'LineWidth', lw, 'Color', default_colors(2, :), 'LineStyle', '--')
% plot(t_13_state_2, X_13_state_2(:, [model.I.AVenom]), 'LineWidth', lw, 'Color', default_colors(3, :), 'LineStyle', ':')
plot(t_2_full_base, X_2_full_base(:, [model.I.AVenom]), 'LineWidth', lw, 'Color', default_colors(1, :), 'LineStyle', '--')
plot(t_1_full, X_1_full(:, [model.I.AVenom]), 'LineWidth', lw, 'Color', default_colors(1, :))
xline(t_star, '-k', 'LineWidth', lw, 'LineStyle', ':')
xlim(xlimits)
ylim([-0.03*8e-3, 8e-3])
ylabel("concentration [nM]")
xlabel("time [h]")
% legend(["original", "8-state", "13-state"])
% legend(["antivenom", "no antivenom"])
title("AVenom")
box on
grid on

ax = nexttile(tlay);
hold on;
plot(t_full_2, X_full_2(:, [model.I.CVenom]), 'LineWidth', lw, 'Color', default_colors(1, :))
% plot(t_8_state_2, X_8_state_2(:, [model.I.CVenom]), 'LineWidth', lw, 'Color', default_colors(2, :), 'LineStyle', '--')
% plot(t_13_state_2, X_13_state_2(:, [model.I.CVenom]), 'LineWidth', lw, 'Color', default_colors(3, :), 'LineStyle', ':')
plot(t_2_full_base, X_2_full_base(:, [model.I.CVenom]), 'LineWidth', lw, 'Color', default_colors(1, :), 'LineStyle', '--')
plot(t_1_full, X_1_full(:, [model.I.CVenom]), 'LineWidth', lw, 'Color', default_colors(1, :))
xline(t_star, '-k', 'LineWidth', lw, 'LineStyle', ':')
xlim(xlimits)
ylim([-0.03*4e-3, 4e-3])
xlabel("time [h]")
% legend(["original", "8-state", "13-state"])
legend(["antivenom", "no antivenom"], 'Location', 'eastoutside')
title("CVenom")
box on
grid on

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, ['./figures/BC_antivenom_sim_', char(num2str(t_star)), save_gulati, '.pdf'], 'ContentType', 'vector');

%% visualize n=1000 UFa Fg trajectories

load("modelBC_SV40_from_JKn_2024.mat")
% load("modelBC_SV40_population_CV20.mat")
load("modelBC_SV40_population_CV40.mat")
% load("modelBC_population_from_UFa_2023")

X = variability.X_ref_pop{1};
historical = [model.t_ref-model.t_ref(2) X(:, model.I.Fg)];
forecast = zeros([length(model.t_ref), 1001]);
forecast(:, 1) = model.t_ref;
for i = 1:1000
    X = variability.X_ref_pop{i+1};
    forecast(:, i+1) = X(:, model.I.Fg);
end

hold on
fanplot(historical, forecast)
yscale log