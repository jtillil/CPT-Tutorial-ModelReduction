simple_calculator_gui()

function simple_calculator_gui()
    % Create main figure window
    fig = uifigure('Name', 'Simple Calculator', 'Position', [100 100 400 300]);
    
    % Create input fields
    num1_field = uieditfield(fig, 'numeric', 'Position', [50 220 100 22], ...
                            'Value', 0, 'ValueDisplayFormat', '%.2f');
    
    num2_field = uieditfield(fig, 'numeric', 'Position', [250 220 100 22], ...
                            'Value', 0, 'ValueDisplayFormat', '%.2f');
    
    % Create labels
    uilabel(fig, 'Position', [50 245 100 22], 'Text', 'First Number:');
    uilabel(fig, 'Position', [250 245 100 22], 'Text', 'Second Number:');
    uilabel(fig, 'Position', [175 220 50 22], 'Text', '+', ...
            'HorizontalAlignment', 'center', 'FontSize', 16);
    
    % Create result display
    result_field = uieditfield(fig, 'text', 'Position', [150 150 100 22], ...
                              'Editable', 'off', 'Value', '0');
    uilabel(fig, 'Position', [150 175 100 22], 'Text', 'Result:');
    
    % Create operation buttons
    add_btn = uibutton(fig, 'push', 'Position', [50 100 60 30], ...
                       'Text', 'Add (+)', 'ButtonPushedFcn', @(btn, event) calculate('add'));
    
    sub_btn = uibutton(fig, 'push', 'Position', [120 100 60 30], ...
                       'Text', 'Subtract (-)', 'ButtonPushedFcn', @(btn, event) calculate('subtract'));
    
    mul_btn = uibutton(fig, 'push', 'Position', [190 100 60 30], ...
                       'Text', 'Multiply (×)', 'ButtonPushedFcn', @(btn, event) calculate('multiply'));
    
    div_btn = uibutton(fig, 'push', 'Position', [260 100 60 30], ...
                       'Text', 'Divide (÷)', 'ButtonPushedFcn', @(btn, event) calculate('divide'));
    
    % Create clear button
    clear_btn = uibutton(fig, 'push', 'Position', [330 100 60 30], ...
                         'Text', 'Clear', 'ButtonPushedFcn', @(btn, event) clear_all());
    
    % Create plot button for demonstration
    plot_btn = uibutton(fig, 'push', 'Position', [150 50 100 30], ...
                        'Text', 'Plot Sin Wave', 'ButtonPushedFcn', @(btn, event) plot_demo());
    
    % Callback function for calculations
    function calculate(operation)
        num1 = num1_field.Value;
        num2 = num2_field.Value;
        
        switch operation
            case 'add'
                result = num1 + num2;
            case 'subtract'
                result = num1 - num2;
            case 'multiply'
                result = num1 * num2;
            case 'divide'
                if num2 == 0
                    result_field.Value = 'Error: Division by zero';
                    return;
                else
                    result = num1 / num2;
                end
            otherwise
                result = 0;
        end
        
        result_field.Value = sprintf('%.4f', result);
    end
    
    % Callback function to clear all fields
    function clear_all()
        num1_field.Value = 0;
        num2_field.Value = 0;
        result_field.Value = '0';
    end
    
    % Callback function for plot demonstration
    function plot_demo()
        % Create new figure for plot
        plot_fig = figure('Name', 'Sin Wave Demo', 'Position', [500 200 500 300]);
        x = 0:0.1:4*pi;
        y = sin(x);
        plot(x, y, 'b-', 'LineWidth', 2);
        grid on;
        title('Sine Wave');
        xlabel('x');
        ylabel('sin(x)');
    end
end

% Alternative approach using GUIDE (older method)
% To create a GUI using GUIDE:
% 1. Type 'guide' in MATLAB command window
% 2. Select "Blank GUI"
% 3. Drag and drop components from the palette
% 4. Right-click components to set properties
% 5. Save the .fig file
% 6. MATLAB will generate a .m file with callback functions

% Example of programmatic GUI with more advanced features
function advanced_gui_example()
    % Create main figure
    fig = uifigure('Name', 'Advanced GUI Example', 'Position', [200 200 600 500]);
    
    % Create tab group
    tabgroup = uitabgroup(fig, 'Position', [10 10 580 480]);
    
    % Tab 1: Data Input
    tab1 = uitab(tabgroup, 'Title', 'Data Input');
    
    % Create table for data input
    data_table = uitable(tab1, 'Position', [20 200 540 200], ...
                         'Data', magic(5), ...
                         'ColumnEditable', true);
    
    % Load data button
    load_btn = uibutton(tab1, 'push', 'Position', [20 150 100 30], ...
                        'Text', 'Load Data', ...
                        'ButtonPushedFcn', @(btn, event) load_data());
    
    % Save data button
    save_btn = uibutton(tab1, 'push', 'Position', [140 150 100 30], ...
                        'Text', 'Save Data', ...
                        'ButtonPushedFcn', @(btn, event) save_data());
    
    % Tab 2: Visualization
    tab2 = uitab(tabgroup, 'Title', 'Visualization');
    
    % Create axes for plotting
    ax = uiaxes(tab2, 'Position', [20 100 540 300]);
    
    % Plot type dropdown
    plot_dropdown = uidropdown(tab2, 'Position', [20 50 150 22], ...
                               'Items', {'Line Plot', 'Bar Chart', 'Histogram', 'Scatter'}, ...
                               'Value', 'Line Plot');
    
    % Plot button
    plot_btn = uibutton(tab2, 'push', 'Position', [190 50 100 30], ...
                        'Text', 'Create Plot', ...
                        'ButtonPushedFcn', @(btn, event) create_plot());
    
    % Callback functions
    function load_data()
        [file, path] = uigetfile('*.mat', 'Select MAT file');
        if file ~= 0
            load(fullfile(path, file));
            % Assume loaded data is in variable 'data'
            if exist('data', 'var')
                data_table.Data = data;
            end
        end
    end
    
    function save_data()
        data = data_table.Data;
        [file, path] = uiputfile('*.mat', 'Save data as');
        if file ~= 0
            save(fullfile(path, file), 'data');
        end
    end
    
    function create_plot()
        data = data_table.Data;
        plot_type = plot_dropdown.Value;
        
        switch plot_type
            case 'Line Plot'
                plot(ax, data);
            case 'Bar Chart'
                bar(ax, mean(data, 1));
            case 'Histogram'
                histogram(ax, data(:));
            case 'Scatter'
                if size(data, 2) >= 2
                    scatter(ax, data(:,1), data(:,2));
                end
        end
        
        title(ax, plot_type);
        grid(ax, 'on');
    end
end

% To run these GUIs, save each function in separate .m files and call:
% simple_calculator_gui()
% advanced_gui_example()