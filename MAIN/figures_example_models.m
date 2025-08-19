%% Setup
clear; clc;
addpath(genpath("../../jti-code"))
sizex = 10;
sizey = 5;
lw = 1;
lwt = 0.5;
interpreter = 'latex';
% Set default interpreter to LaTeX
set(groot, 'defaultTextInterpreter', interpreter);
set(groot, 'defaultAxesTickLabelInterpreter', interpreter);
set(groot, 'defaultLegendInterpreter', interpreter);
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

%% Parallel Pathways: Scenario 1 - no crosstalk
load("modelSPP_no_crosstalk_full.mat")

% reference solution
figure
p1 = semilogy(model.t_ref, model.X_ref, 'LineWidth', lw); %DisplayName', plotnames(i))
xlim([-0.002 0.052])
ylim([1e-3 1e4])
xlabel("t [min]")
ylabel("concentration [nM]")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]
exportgraphics(gcf, "./figures/SPP_no_crosstalk_ref_sol.pdf")

legend('A', 'S', 'B', 'C', 'D', 'Location','eastoutside')
exportgraphics(gcf, "./figures/SPP_no_crosstalk_ref_sol_legend.pdf")

% set(p1, 'visible', 'off');
% set(gca, 'visible', 'off');
% exportgraphics(gcf, "./figures/SPP_no_crosstalk_ref_sol_legend.pdf");
% exportgraphics(leg, "./figures/SPP_no_crosstalk_ref_sol_legend.pdf")
% saveas(leg, 'my_legend.png');

% nir-indices
figure
plot(model.t_ref, model.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-0.002 0.052])
ylim([-0.01 1])
xlabel("t [min]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]
exportgraphics(gcf, "./figures/SPP_no_crosstalk_ir.pdf")

legend('A', 'S', 'B', 'C', 'D', 'threshold', 'Location','eastoutside')
exportgraphics(gcf, "./figures/SPP_no_crosstalk_ir_legend.pdf")

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
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]
exportgraphics(gcf, "./figures/SPP_no_crosstalk_non_ir_B.pdf")

legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','eastoutside')
exportgraphics(gcf, "./figures/SPP_no_crosstalk_non_ir_B_legend.pdf")

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
xlabel("t [min]")
ylabel("normalised index")
% hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]
exportgraphics(gcf, "./figures/SPP_no_crosstalk_non_ir_C.pdf")

legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','eastoutside')
exportgraphics(gcf, "./figures/SPP_no_crosstalk_non_ir_C_legend.pdf")

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
xlabel("t [min]")
ylabel("normalised index")
% hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]
exportgraphics(gcf, "./figures/SPP_no_crosstalk_non_ir_S.pdf")

legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','eastoutside')
exportgraphics(gcf, "./figures/SPP_no_crosstalk_non_ir_S_legend.pdf")

%% Parallel Pathways: Scenario 2 - with crosstalk
load("modelSPP_with_crosstalk_full.mat")

% reference solution
figure
semilogy(model.t_ref, model.X_ref, 'LineWidth', lw) %DisplayName', plotnames(i))
xlim([-0.002 0.052])
ylim([1e-3 1e4])
xlabel("t [min]")
ylabel("concentration [nM]")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]
exportgraphics(gcf, "./figures/SPP_with_crosstalk_ref_sol.pdf")

legend('A', 'S', 'B', 'C', 'D', 'Location','eastoutside')
exportgraphics(gcf, "./figures/SPP_with_crosstalk_ref_sol_legend.pdf")

% nir-indices
figure
plot(model.t_ref, model.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-0.002 0.052])
ylim([-0.01 1])
xlabel("t [min]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]
exportgraphics(gcf, "./figures/SPP_with_crosstalk_ir.pdf")

legend('A', 'S', 'B', 'C', 'D', 'threshold', 'Location','eastoutside')
exportgraphics(gcf, "./figures/SPP_with_crosstalk_ir_legend.pdf")

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
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]
exportgraphics(gcf, "./figures/SPP_with_crosstalk_non_ir_S.pdf")

legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','eastoutside')
exportgraphics(gcf, "./figures/SPP_with_crosstalk_non_ir_S_legend.pdf")

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
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]
exportgraphics(gcf, "./figures/SPP_with_crosstalk_non_ir_B.pdf")

legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','eastoutside')
exportgraphics(gcf, "./figures/SPP_with_crosstalk_non_ir_B_legend.pdf")

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
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]
exportgraphics(gcf, "./figures/SPP_with_crosstalk_non_ir_C.pdf")

legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','eastoutside')
exportgraphics(gcf, "./figures/SPP_with_crosstalk_non_ir_C_legend.pdf")

%% Parallel Pathways: Combined Figure

% Create figure
figure('Units', 'centimeters', 'Position', [0, 0, 18, 20]);

% --- MODIFIED: Reverted to 'compact' spacing ---
t = tiledlayout(6, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Define constants and labels for clarity ---
lw = 1.5; % Line width for plots
lwt = 1;  % Line width for threshold lines

% Define title strings in a 5x2 cell array for unique letters
title_strings = {
    "B) Concentration-time profiles",   "C) Concentration-time profiles";
    "C) Normalized ir-indices",         "E) Normalized ir-indices";
    "D) State-classification indices for S", "G) State-classification indices for S";
    "E) State-classification indices for B", "I) State-classification indices for B";
    "F) State-classification indices for C", "K) State-classification indices for C";
    "G) Normalized ir-indices (reduced models)",         "M) Normalized ir-indices (reduced models)"
};

% Define y-labels in a cell array
ylabel_strings = {
    "concentration [nM]",
    "nir-index",
    "index / error",
    "index / error",
    "index / error",
    "nir-index"
};

% --- Loop to create all plots ---
for row = 1:6
    for col = 1:3
        if col == 2 && row ~= 6
            % ax = nexttile;
            continue
        end
        % Load the appropriate data file for the column
        if col == 1
            load("modelSPP_no_crosstalk_full.mat")
        else
            load("modelSPP_with_crosstalk_full.mat")
        end

        % Create the next tile and get its handle
        ax = nexttile(3*(row-1)+col);
        hold(ax, 'on'); % Use hold on the specific axes
        co = get(ax, 'ColorOrder');
        
        % Plot data based on the current row
        switch row
            case 1 % Concentration-time profiles
                semilogy(ax, model.t_ref, model.X_ref, 'LineWidth', lw);
                ylim([1e-3 1e4]);
                yticks([1e-2, 1e1, 1e4]);
                set(ax, 'YScale', 'log');

            case 2 % Normalised ir-indices (no threshold)
                plot(ax, model.t_ref, model.ir.nindex, 'LineWidth', lw);
                ylim([-0.04 1.04]);

            case {3, 4, 5} % State-classification indices (S, B, C)
                % Determine which species to plot based on the row
                switch row
                    case 3, idx = model.I.S;
                    case 4, idx = model.I.B;
                    case 5, idx = model.I.C;
                end
                
                % Plot nindex AND relstateerr for env and pss
                semilogy(ax, model.t_ref, model.env.nindex(:, idx), 'LineStyle', '-', 'Color', co(1,:), 'LineWidth', lw);
                semilogy(ax, model.t_ref, model.env.relstateerr(:, idx), 'LineStyle', ':', 'Color', co(1,:), 'LineWidth', lw);
                
                semilogy(ax, model.t_ref, model.pss.nindex(:, idx), 'LineStyle', '--', 'Color', co(2,:), 'LineWidth', lw);
                semilogy(ax, model.t_ref, model.pss.relstateerr(:, idx), 'LineStyle', ':', 'Color', co(2,:), 'LineWidth', lw);
                
                % Plot ONLY nindex for cneg and pneg
                semilogy(ax, model.t_ref, model.cneg.nindex(:, idx), 'LineStyle', '-', 'Color', co(3,:), 'LineWidth', lw);
                semilogy(ax, model.t_ref, model.pneg.nindex(:, idx), 'LineStyle', '--', 'Color', co(4,:), 'LineWidth', lw);

                % Common settings for these rows
                yline(ax, 0.1, 'k--', 'LineWidth', lwt);
                ylim([5e-4 1e1]);
                yticks([1e-3, 1e-1, 1e1]);
                set(ax, 'YScale', 'log');

            case 6 % Normalised ir-indices (no threshold)
                switch col
                    case 1
                        load("indices_SPP_no_crosstalk_red.mat")
                        plot(ax, t_ind, ir.nindex(:, [1]), 'LineWidth', lw, 'Color', co([1],:));
                        plot(ax, t_ind, ir.nindex(:, [3]), 'LineWidth', lw, 'Color', co([3],:));
                        plot(ax, t_ind, ir.nindex(:, [5]), 'LineWidth', lw, 'Color', co([5],:));
                        % plot(ax, t_ind, ir.nindex(:, [1 3 5]), 'LineWidth', lw, 'Color', co([1 3 5],:));
                        ylim([-0.04 1.04]);
                    case 2
                        load("indices_SPP_with_crosstalk_red_correct.mat")
                        plot(ax, t_ind, ir.nindex(:, [1]), 'LineWidth', lw, 'Color', co([1],:));
                        plot(ax, t_ind, ir.nindex(:, [2]), 'LineWidth', lw, 'Color', co([2],:));
                        plot(ax, t_ind, ir.nindex(:, [3]), 'LineWidth', lw, 'Color', co([3],:));
                        plot(ax, t_ind, ir.nindex(:, [5]), 'LineWidth', lw, 'Color', co([5],:));
                        % plot(ax, t_ind, ir.nindex(:, [1 3 5]), 'LineWidth', lw, 'Color', co([1 2 3 5],:));
                        ylim([-0.04 1.04]);
                    case 3
                        load("indices_SPP_with_crosstalk_red_wrong_pss_solved.mat")
                        plot(ax, t_ind, ir.nindex(:, [1]), 'LineWidth', lw, 'Color', co([1],:));
                        plot(ax, t_ind, ir.nindex(:, [2]), 'LineWidth', lw, 'Color', co([2],:));
                        % plot(ax, t_ind, ir.nindex(:, [4]), 'LineWidth', lw, 'Color', co([4],:));
                        plot(ax, t_ind, ir.nindex(:, [5]), 'LineWidth', lw, 'Color', co([5],:));
                        % plot(ax, t_ind, ir.nindex(:, [1 2 4 5]), 'LineWidth', lw, 'Color', co([1 2 4 5],:));
                        ylim([-0.04 1.04]);
                end
        end
        
        % Common settings for all plots
        grid(ax, 'on')
        box(ax, 'on');
        xlim(ax, [-0.002 0.052]);
        
        % --- Decorations (Titles, Labels, Legends) ---

        % --- MODIFIED: Reverted to text() inside the plot for perfect alignment ---
        % text(ax, 0.02, 0.98, title_strings{row, col}, ...
        switch col
            case 1
                text(ax, -0.23, 1.13, title_strings{row, col}, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', ...
                'FontSize', 10, ...
                'Interpreter', 'none');%, ...
                % 'BackgroundColor', 'w'); % Add background for readability
            case 3
                col2 = col-1;
                text(ax, -0.13, 1.13, title_strings{row, col2}, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontWeight', 'bold', ...
                'FontSize', 10, ...
                'Interpreter', 'none');%, ...
                % 'BackgroundColor', 'w'); % Add background for readability
        end

        hold(ax, 'off');

        % Y-labels only for the left column
        if col == 1
            ylabel(ylabel_strings{row});
        else
            % Right column now has y-tick labels
        end
        
        % X-labels only for the bottom row
        if row == 6
            xlabel("t [min]");
        else
            % Remove x-tick labels for all other rows
            set(ax, 'XTickLabel', []);
        end
        
        % Legends only for the right column
        if col == 3
            if row == 1
                lg = legend(ax, 'A', 'S', 'B', 'C', 'D', 'Location','eastoutside');
                % lg.Layout.Tile = 2;
            elseif row == 3 % Rows 3, 4, 5 have the same legend
                legend_labels = {
                    'env (index)', 'env (error)', ...
                    'pss (index)', 'pss (error)', ...
                    'cneg (index)', ...
                    'pneg (index)', ...
                    'threshold'
                };
                lg = legend(ax, legend_labels, 'Location','eastoutside');
                % lg.Layout.Tile = 8;
            elseif row == 6
                % lg = legend(ax, 'A (\#1)', 'S (\#1)', 'B (\#1)', 'D (\#1)', 'A (\#2)', 'S (\#2)', 'D (\#2)', 'Location','eastoutside');
            end
            lg.FontSize = 8;
        end
    end
end

% Export the final figure
exportgraphics(gcf, "./figures/SPP_complete.pdf", 'ContentType', 'vector');


%% SPP complete new

% Create figure
figure('Units', 'centimeters', 'Position', [0, 0, 18, 20]);

% Using 'compact' spacing. Titles added with text() won't auto-adjust layout.
t = tiledlayout(6, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Define constants and labels for clarity ---
lw = 1.5; % Line width for plots
lwt = 1;  % Line width for threshold lines

% Titles for the first column of plots only.
title_strings_col1 = {
    ["Concentration-time"; "profiles"],
    ["Normalized"; "ir-indices"],
    "", % Row 3 has no main title, only the 'S' annotation
    "State-classification indices", % Shared title for rows 3-5
    "", % Row 5 has no main title, only the 'C' annotation
    ["Normalized ir-indices"; "(reduced models)"]
};

% Annotations for state-classification rows
state_labels = {"S", "B", "C"};

% Define y-labels in a cell array
ylabel_strings = {
    "concentration [nM]",
    "nir-index",
    "index / error",
    "index / error",
    "index / error",
    "nir-index"
};

% Initialize a counter for subplot labels (B, C, D...)
plot_counter = 0;

% --- Loop to create all plots ---
for row = 1:6
    for col = 1:3
        % Skip the empty middle-column plots
        if col == 2 && row ~= 6
            continue
        end
        % Load the appropriate data file for the column
        if col == 1
            load("modelSPP_no_crosstalk_full.mat")
        else
            load("modelSPP_with_crosstalk_full.mat")
        end

        % Create the next tile and get its handle
        ax = nexttile(3*(row-1)+col);
        hold(ax, 'on'); % Use hold on the specific axes
        co = get(ax, 'ColorOrder');
        
        % Plot data based on the current row
        switch row
            case 1 % Concentration-time profiles
                semilogy(ax, model.t_ref, model.X_ref, 'LineWidth', lw);
                ylim([1e-3 1e4]);
                yticks([1e-2, 1e1, 1e4]);
                set(ax, 'YScale', 'log');

            case 2 % Normalised ir-indices (no threshold)
                plot(ax, model.t_ref, model.ir.nindex, 'LineWidth', lw);
                ylim([-0.04 1.04]);

            case {3, 4, 5} % State-classification indices (S, B, C)
                % Determine which species to plot based on the row
                switch row
                    case 3, idx = model.I.S;
                    case 4, idx = model.I.B;
                    case 5, idx = model.I.C;
                end
                
                semilogy(ax, model.t_ref, model.env.nindex(:, idx), 'LineStyle', '-', 'Color', co(1,:), 'LineWidth', lw);
                semilogy(ax, model.t_ref, model.env.relstateerr(:, idx), 'LineStyle', ':', 'Color', co(1,:), 'LineWidth', lw);
                semilogy(ax, model.t_ref, model.pss.nindex(:, idx), 'LineStyle', '--', 'Color', co(2,:), 'LineWidth', lw);
                semilogy(ax, model.t_ref, model.pss.relstateerr(:, idx), 'LineStyle', ':', 'Color', co(2,:), 'LineWidth', lw);
                semilogy(ax, model.t_ref, model.cneg.nindex(:, idx), 'LineStyle', '-', 'Color', co(3,:), 'LineWidth', lw);
                semilogy(ax, model.t_ref, model.pneg.nindex(:, idx), 'LineStyle', '--', 'Color', co(4,:), 'LineWidth', lw);

                yline(ax, 0.1, 'k--', 'LineWidth', lwt);
                ylim([5e-4 1e1]);
                yticks([1e-3, 1e-1, 1e1]);
                set(ax, 'YScale', 'log');

            case 6 % Normalised ir-indices (reduced models)
                switch col
                    case 1
                        load("indices_SPP_no_crosstalk_red.mat")
                        plot(ax, t_ind, ir.nindex(:, [1]), 'LineWidth', lw, 'Color', co([1],:));
                        plot(ax, t_ind, ir.nindex(:, [3]), 'LineWidth', lw, 'Color', co([3],:));
                        plot(ax, t_ind, ir.nindex(:, [5]), 'LineWidth', lw, 'Color', co([5],:));
                        % plot(ax, t_ind, ir.nindex(:, [1 3 5]), 'LineWidth', lw, 'Color', co([1 3 5],:));
                        ylim([-0.04 1.04]);
                    case 2
                        load("indices_SPP_with_crosstalk_red_correct.mat")
                        plot(ax, t_ind, ir.nindex(:, [1]), 'LineWidth', lw, 'Color', co([1],:));
                        plot(ax, t_ind, ir.nindex(:, [2]), 'LineWidth', lw, 'Color', co([2],:));
                        plot(ax, t_ind, ir.nindex(:, [3]), 'LineWidth', lw, 'Color', co([3],:));
                        plot(ax, t_ind, ir.nindex(:, [5]), 'LineWidth', lw, 'Color', co([5],:));
                        % plot(ax, t_ind, ir.nindex(:, [1 3 5]), 'LineWidth', lw, 'Color', co([1 2 3 5],:));
                        ylim([-0.04 1.04]);
                    case 3
                        load("indices_SPP_with_crosstalk_red_wrong_pss_solved.mat")
                        plot(ax, t_ind, ir.nindex(:, [1]), 'LineWidth', lw, 'Color', co([1],:));
                        plot(ax, t_ind, ir.nindex(:, [2]), 'LineWidth', lw, 'Color', co([2],:));
                        % plot(ax, t_ind, ir.nindex(:, [4]), 'LineWidth', lw, 'Color', co([4],:));
                        plot(ax, t_ind, ir.nindex(:, [5]), 'LineWidth', lw, 'Color', co([5],:));
                        % plot(ax, t_ind, ir.nindex(:, [1 2 4 5]), 'LineWidth', lw, 'Color', co([1 2 4 5],:));
                        ylim([-0.04 1.04]);
                end
        end
        
        grid(ax, 'on')
        box(ax, 'on');
        xlim(ax, [-0.002 0.052]);
        
        % --- Decorations (Titles, Labels, Legends) ---
        
        % --- START OF MODIFIED SECTION =======================================
        % Add panel labels (B, C, D...) as text ABOVE each plot
        letter_label = sprintf('%c)', 'B' + plot_counter);
        
        % Use text() with normalized coordinates to place it above the axes.
        % x=0 places it at the left edge. y=1.05 places it just above the top.
        % 'VerticalAlignment','bottom' aligns the bottom of the text to this y-coordinate.
        text(ax, 0, 1.05, letter_label, ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', ... 
            'FontWeight', 'bold', ...
            'FontName', 'Helvetica', ...
            'FontSize', 10, 'Interpreter', 'none');
        
        plot_counter = plot_counter + 1; % Increment for the next plot
        
        % Place vertical titles and annotations only for the first column
        if col == 1
            % --- Main Vertical Titles ---
            the_title = title_strings_col1{row};
            if ~isempty(the_title)
                % Ensure font consistency with the panel labels
                text(ax, -0.45, 0.5, the_title, ...
                    'Units', 'normalized', 'Rotation', 90, ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold', 'FontName', 'Helvetica', ...
                    'FontSize', 10, 'Interpreter', 'none');
            end
            
            % --- "S", "B", "C" Annotations ---
            if ismember(row, [3, 4, 5])
                annotation_text = state_labels{row - 2};
                text(ax, -0.25, 0.5, annotation_text, ...
                    'Units', 'normalized',...% 'Rotation', 90, ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold', 'FontName', 'Helvetica', ...
                    'FontSize', 12);
            end
        end
        % --- END OF MODIFIED SECTION =========================================

        hold(ax, 'off');

        % Y-labels only for the left column
        if col == 1
            ylabel(ylabel_strings{row});
        end
        
        % X-labels only for the bottom row
        if row == 6
            xlabel("t [min]");
        else
            set(ax, 'XTickLabel', []);
        end
        
        % Legends only for the right column
        if col == 3
            if row == 1
                lg = legend(ax, 'A', 'S', 'B', 'C', 'D', 'Location','eastoutside');
            elseif row == 3
                legend_labels = {
                    'env (index)', 'env (error)', ...
                    'pss (index)', 'pss (error)', ...
                    'cneg (index)', ...
                    'pneg (index)', ...
                    'threshold'
                };
                lg = legend(ax, legend_labels, 'Location','eastoutside');
            end
            if exist('lg', 'var') && isvalid(lg), lg.FontSize = 8; end
        end
    end
end

% Export the final figure
exportgraphics(gcf, "./figures/SPP_complete.pdf", 'ContentType', 'vector');

%% code pss states into wrong reduced model scenario 2

load("modelSPP_with_crosstalk_full.mat")

% obtain pss states
config = repmat("dyn", [model.I.nstates, 1]);
config(model.I.B) = "cneg";
config(model.I.C) = "pss";
% I_red = config2I(model.I, config, []);

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
if length(X_sym_pss) > 1
    for k = 1:length(X_sym_pss)
        solutions = G.(nm_X_pss{k});
        odefun_symbolic_solved = subs(odefun_symbolic, X_sym_pss(k), solutions(1));
    end
else
    odefun_symbolic_solved = subs(odefun_symbolic, X_sym_pss(k), G(1));
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
        jacfun_symbolic_solved(i, j) = diff(odefun_symbolic_solved(i), X_sym(j));
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

%% check pss solved wrong reduced model scenario 2

% config(config == "pss") = "pneg";
% model.I = config2I(model.I, config, []);
model.multiple.multiple = 0;
config(model.I.C) = "pneg";
[err_index_solved, ~, tred_solved, Xred_solved] = objfun(model.t_ref, model.X_ref, model.X0, model.par, model.I, [], model.param, model.multiple, model.odefun, model.jacfun, config, "MRSE");

%% ir-indices of pss solved wrong reduced model scenario 2

% model.t_ref = model.t_ref([1:80:end]);
% model.X_ref = model.X_ref([1:80:end], :);

% config = repmat("dyn", [model.I.nstates, 1]);
model.I = config2I(model.I, config, []);

[ir, contr, obs, t_ind] = compute_ir_indices_matlabfun(model);

save("indices_SPP_with_crosstalk_red_wrong_pss_solved.mat", "ir", "contr", "obs", "t_ind")

%% Parallel Pathways: Scenario 3 - C B neglectable
load("modelSPP_C_B_neglectable_full.mat")

% reference solution
figure
semilogy(model.t_ref, model.X_ref, 'LineWidth', lw) %DisplayName', plotnames(i))
xlim([-0.002 0.052])
ylim([1e-3 1e4])
legend('A', 'S', 'B', 'C', 'D', 'Location','northeast')
xlabel("t [min]")
ylabel("concentration [nM]")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_C_B_neglectable_ref_sol.pdf")

% nir-indices
figure
plot(model.t_ref, model.ir.nindex, 'LineWidth', lw) %DisplayName', plotnames(i))
yline(0.1, 'k--', 'LineWidth', lwt)
xlim([-0.002 0.052])
ylim([-0.01 1])
legend('A', 'S', 'B', 'C', 'D', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("nir-index")

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_C_B_neglectable_ir.pdf")

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
legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_C_B_neglectable_non_ir_B.pdf")

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
legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
% hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_C_B_neglectable_non_ir_C.pdf")

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
legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
% hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

exportgraphics(gcf, "./figures/SPP_C_B_neglectable_non_ir_S.pdf")

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

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

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

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

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
legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

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
legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

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

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

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

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

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
legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

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
legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("t [min]")
ylabel("normalised index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, sizex, sizey]); % [x, y, width, height]

exportgraphics(gcf, "./figures/MMEK_Cpss_EpCenv_non_ir_C.pdf")

%% figures for BioRender graphics

% ordered ir indices
states = 1:12;
% ind = ((13 - states)/12).^4;
ind = [0 0 0 0.01 0.03 0.04 0.07 0.08 0.42 0.55 1 1];

figure
hold on
b = scatter(states, ind, "filled");%, 'MarkerColor', "#0072BD");
% yline(0.15, 'k--', 'LineWidth', 1)
xlim([0.3 12.7])
ylim([-0.06 1.06])
% set(gca, 'YScale', 'log')
box on
grid on
% legend('const', 'qss', 'cneg', 'pneg', 'threshold', 'Location','northeast')
xlabel("ordered states")
ylabel("max normalized ir-index")
hold off

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 10, 8]); % [x, y, width, height]

exportgraphics(gcf, "./figures/BioRender_ordered_ir_indices.png")

