% generateSummaryEVT.m
%
% This function generates summary plots for average angular velocity versus object position
% across different steering gain (k) values, allowing a comparison across multiple conditions.
% Each k value is displayed in a separate tile, and all conditions are plotted together.
%
% INPUTS:
%   nK                - Number of steering gain values (k) in the experiment.
%   kValues           - Vector of k values used in the experiment.
%   evt               - 3D array of average angular velocity (k x position bins x conditions).
%   posBins           - Vector defining the position bins (degrees).
%   comparisonLabel   - Cell array with labels for all conditions being compared.
%   comparisonType    - Name of comparison
%   folder            - Struct with fields 'final' and 'vectors' specifying save locations.
%
% OUTPUTS:
%   Generates a figure with plots showing the relationship between angular velocity and position for each k value.
%
% CREATED: 10/30/2024 - MC
% UPDATED: 11/16/2024 - MC refactored for multiple conditions
%

function generateSummaryEVT(nK, kValues, evt, posBins, comparisonLabel, comparisonType, folder)

    % Determine number of conditions
    numConditions = size(evt, 3);

    % Generate summary plot for Angular Velocity vs Object Position
    figure;
    set(gcf, 'Position', [100 100 1500 500]);  % Set figure size
    tiledlayout(1, nK, 'TileSpacing', 'compact');  % One row for each k value

    % Loop over k values for plotting
    for kIdx = 1:nK
        nexttile; hold on;
        
        % Plot the relationship for all conditions
        for condIdx = 1:numConditions
            plot(posBins, evt(kIdx, :, condIdx), '-', 'linewidth', 1, 'DisplayName', comparisonLabel{condIdx});
        end

        % Customize the plot
        xlabel('Position (deg)');
        ylabel('Avg. Angular Velocity (deg/s)');
        title(['k=' num2str(kValues(kIdx))]);
        if kIdx == 1
            legend('show', 'Location', 'northwest');
        end
        grid on;
        ylim([-150 150]);
        xlim([-80 80]);
    end

    % Add overall title
    sgtitle(['Binned Angular Velocity vs Object Position Across Conditions (' ...
             strjoin(comparisonLabel, ' vs ') ')']);

    % Save EVT plot in both PNG and SVG formats
    saveas(gcf, fullfile(folder.final, [comparisonType '_Binned_EVT.png']));
    set(gcf, 'renderer', 'Painters');  % Ensure high-resolution SVG output by setting renderer
    saveas(gcf, fullfile(folder.vectors, [comparisonType '_Binned_EVT.svg']));

    % Apply zoomed x-axis limits to each tile
    for tileIdx = 1:nK
        nexttile(tileIdx);  % Select each tile in the current figure
        xlim([-30 30]);     % Set x-axis limits for zoomed view
    end

    % Save zoomed EVT plot in both PNG and SVG formats
    saveas(gcf, fullfile(folder.final, [comparisonType '_Binned_EVT_Zoomed.png']));
    set(gcf, 'renderer', 'Painters');  % Set renderer for vector quality in zoomed version
    saveas(gcf, fullfile(folder.vectors, [comparisonType '_Binned_EVT_Zoomed.svg']));

end
