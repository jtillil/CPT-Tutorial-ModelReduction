%% Setup
clear; clc;
addpath(genpath("../."))
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

%% Figure 5 - Fg timecourses

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

%% Figure 5 - state trajectories; ir indices


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

%% Figure 6

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
title_text_scatter = "a   Ordering of states by decreasing dynamic relevance";
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
light_gray = [0.5 0.5 0.5]; %[0.8, 0.8, 0.8]; 
light_gold = 0.5 * mygold + 0.5 * white;
light_blue = 0.5 * myblue + 0.5 * white;

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
        text(ax_heat, -0.33, 1.07, "b   Analysis of reduction approaches", ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');
    end
end

% --- 4. Export the Figure ---
drawnow;
exportgraphics(f, "./figures/Figure_6.pdf", 'ContentType', 'vector');

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

%% 2 types of normalized indices --- 13-state model

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


%% Supplement: approximation quality in 8-state reduced model

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

%% Supplement: approximation quality in 13-state reduced model

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

%% Supplement: approximation quality in 24-state reduced model

% Load data first, as it's needed for both plot sections
load("../Core/modelfiles/modelBC_SV40_from_JKn_2024.mat")
% load("modelBC_SV40_from_JKn_2024_reduced_13_state.mat")

I = model.I;
% redconfig = repmat("pneg", [1 model.I.nstates]);
% redconfig([I.AVenom, I.CVenom, I.Fg, I.F, I.II, I.IIa, I.V, I.Va, I.VIII, I.VIIIa, I.IX, I.IXa, I.IXa_VIIIa, I.Xa, I.Xa_Va, ...
% I.XI, I.XIa, I.XIII, I.XIIIa, I.Tmod, I.IIa_Tmod, I.APC, I.APC_PS, I.Pg, I.P]) = "dyn";
% redconfig([I.X, I.PS, I.PC, I.TFPI, I.VKH2]) = "env";

load("./results/modelBC_SV40_from_JKn_2024_reduced_24_state.mat")
redconfig = redmodel.redconfig;

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
exportgraphics(f, "./figures/BC_approximation_24_state_reduced_model.pdf", 'ContentType', 'vector');

%% Supplement: antivenom simulation

for Gulati = [0, 1]
if Gulati == 0
for t_star = [10/60 60/60]
    default_colors = get(groot, 'DefaultAxesColorOrder');
    
    load("modelBC_SV40_from_JKn_2024.mat")
    
    config = repmat("dyn", [1 model.I.nstates]);
    
    % SIMULATION
    % t_star = 1e-20;
    % t_star = 0.01;
    % t_star = 1/60;
    % t_star = 5/60;
    % t_star = 10/60;
    % t_star = 20/60;
    % t_star = 30/60;
    % t_star = 60/60;
    % t_star = 120/60;
    % t_star = 180/60;
    % t_star = 5;
    
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
end
else
for t_star = [180/60 5]
    default_colors = get(groot, 'DefaultAxesColorOrder');
    
    load("modelBC_SV40_from_JKn_2024.mat")
    
    config = repmat("dyn", [1 model.I.nstates]);
    
    % SIMULATION
    % t_star = 1e-20;
    % t_star = 0.01;
    % t_star = 1/60;
    % t_star = 5/60;
    % t_star = 10/60;
    % t_star = 20/60;
    % t_star = 30/60;
    % t_star = 60/60;
    % t_star = 120/60;
    % t_star = 180/60;
    % t_star = 5;
    
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
end
end
end
