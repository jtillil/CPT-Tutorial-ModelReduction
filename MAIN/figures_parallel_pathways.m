%% Setup
clear; clc;
addpath(genpath("../."))
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

%% Figure 4

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
                        load("indices_SPP_with_crosstalk_red_wrong.mat")
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
        letter_label = sprintf('%c', 'b' + plot_counter);
        
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
exportgraphics(gcf, "./figures/Figure_4.pdf", 'ContentType', 'vector');
