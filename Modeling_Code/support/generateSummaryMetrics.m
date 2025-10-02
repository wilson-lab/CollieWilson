% generateSummaryMetrics.m
%
% This function generates summary plots for various metrics across steering gain values (k values),
% comparing multiple conditions side by side. Each metric is displayed in a separate tile for easier
% comparison, with axes and labels appropriately set for clear visualization.
%
% INPUTS:
%   kValues           - Vector of steering gain values (k) used in the experiment.
%   metrics_prob      - Matrix containing probability near setpoint for each condition (columns).
%   metrics_var       - Matrix containing circular variance values for each condition (columns).
%   metrics_ISE       - Matrix containing Integral of Squared Error (ISE) values for each condition (columns).
%   metrics_IAE       - Matrix containing Integral of Absolute Error (IAE) values for each condition (columns).
%   comparisonLabel   - Cell array of labels for all conditions being compared.
%   comparisonType    - Name of comparison
%   folder            - Struct with fields for saving figures (e.g., `final` and `vectors`).
%
% OUTPUTS:
%   Generates a figure with four summary plots showing metrics for each condition across k values.
%
% CREATED: 10/30/2024 - MC
% UPDATED: 11/16/2024 - MC refactored for multiple conditions
%

function generateSummaryMetrics(kValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, comparisonLabel, comparisonType, folder)

    % Determine number of conditions
    numConditions = size(metrics_prob, 2);
    colorList = lines(numConditions);  % Use MATLAB's default color set

    % Generate summary plots
    figure;
    set(gcf, 'Position', [100 100 600 1000]);  % Adjust figure size for better layout
    tiledlayout(4, 1, 'TileSpacing', 'compact');  % Four rows: one for each metric

    % Plot 1: Probability near setpoint across k values
    nexttile;
    hold on;
    for condIdx = 1:numConditions
        plot(kValues, metrics_prob(:, condIdx), '.', ...
            'Color', colorList(condIdx, :), ...
            'MarkerSize', 12, ...
            'LineStyle', '-', ...
            'LineWidth', 1, ...
            'DisplayName', comparisonLabel{condIdx});
    end
    xlabel('Steering Gain (k)');
    ylabel('Probability Near Setpoint');
    legend('show', 'Location', 'best');
    title('Probability Near Setpoint');
    axis padded;
    ylim([0 1]);
    grid on;

    % Plot 2: Circular variance across k values
    nexttile;
    hold on;
    for condIdx = 1:numConditions
        plot(kValues, metrics_var(:, condIdx), '_', ...
            'Color', colorList(condIdx, :), ...
            'MarkerSize', 12, ...
            'LineStyle', '-', ...
            'LineWidth', 1, ...
            'DisplayName', comparisonLabel{condIdx});
    end
    xlabel('Steering Gain (k)');
    ylabel('1 - Circular Variance');
    title('Circular Variance');
    axis padded;
    ylim([0 1]);
    grid on;

    % Plot 3: Integral of Squared Error (ISE) across k values
    nexttile;
    hold on;
    for condIdx = 1:numConditions
        plot(kValues, metrics_ISE(:, condIdx), '_', ...
            'Color', colorList(condIdx, :), ...
            'MarkerSize', 12, ...
            'LineStyle', '-', ...
            'LineWidth', 1, ...
            'DisplayName', comparisonLabel{condIdx});
    end
    xlabel('Steering Gain (k)');
    ylabel('ISE');
    title('Integral of Squared Error (ISE)');
    axis padded;
    grid on;

    % Plot 4: Integral of Absolute Error (IAE) across k values
    nexttile;
    hold on;
    for condIdx = 1:numConditions
        plot(kValues, metrics_IAE(:, condIdx), '_', ...
            'Color', colorList(condIdx, :), ...
            'MarkerSize', 12, ...
            'LineStyle', '-', ...
            'LineWidth', 1, ...
            'DisplayName', comparisonLabel{condIdx});
    end
    xlabel('Steering Gain (k)');
    ylabel('IAE');
    title('Integral of Absolute Error (IAE)');
    axis padded;
    grid on;

    % Overall title
    sgtitle(['Summary of Model Performance Across Metrics (' strjoin(comparisonLabel, ' vs ') ')']);

    % Save summary metrics plot as PNG and SVG
    saveas(gcf, fullfile(folder.final, [comparisonType '_SummaryMetrics.png']));
    set(gcf, 'renderer', 'Painters');  % Set renderer to 'Painters' for better vector graphic rendering
    saveas(gcf, fullfile(folder.vectors, [comparisonType '_SummaryMetrics.svg']));

end
