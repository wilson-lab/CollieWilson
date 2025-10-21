% GENERATESETTLINGSUMMARY Generate a summary plot of settling times across k values.
%
% This function plots the settling times across k values for multiple conditions
% and saves the figure with a title and filenames based on the provided labels.
%
% INPUTS:
%   k               - Array of k values corresponding to rows in avgSettlingTime.
%   avgSettlingTime - 2D array of settling times, where rows are k values
%                     and columns are the conditions.
%   comparisonLabel - Cell array containing strings for each condition.
%   comparisonType  - Name of comparison
%   folder          - Struct with fields 'final' and 'vectors' specifying save locations.
%
% OUTPUTS:
%   None. The function saves the figure in the specified folders.
%
% CREATED: 11/16/2024 - MC
% UPDATED: 11/16/2024 - MC refactored for multiple conditions

function generateSettlingSummary(k, avgSettlingTime, comparisonLabel, comparisonType, folder)
    % Create a figure
    figure;
    hold on;

    % Number of conditions
    numConditions = size(avgSettlingTime, 2);
    colorList = lines(numConditions);  % Distinct colors per condition

    % Plot the settling times for all conditions with connected lines and '_' markers
    for condIdx = 1:numConditions
        plot(k, avgSettlingTime(:, condIdx), ...
            'Color', colorList(condIdx, :), ...
            'Marker', '.', ...
            'MarkerSize', 12, ...
            'LineStyle', '-', ...
            'LineWidth', 1.5, ...
            'DisplayName', comparisonLabel{condIdx});
    end

    % Add labels, legend, and grid
    xlabel('Starting Pos (deg)');
    ylabel('Settling Time (s)');
    legend('Location', 'best');
    axis padded
    ylim([0 15]);  % Adjust as needed
    grid on;

    % Overall title for the figure
    sgtitle(['Settling Times (' strjoin(comparisonLabel, ' vs ') ')']);

    % Save figure as PNG and SVG in specified folders
    saveas(gcf, fullfile(folder.final, [comparisonType '_SettlingSummary.png']));
    set(gcf, 'renderer', 'Painters');  % Use 'Painters' renderer to maintain vector quality
    saveas(gcf, fullfile(folder.vectors, [comparisonType '_SettlingSummary.svg']));
end
