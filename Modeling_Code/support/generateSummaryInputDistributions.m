% generateSummaryInputDistributions.m
%
% This function generates summary plots comparing the distributions of input history 
% across multiple conditions and steering gain (k) values. Each k value is displayed in 
% a separate tile, and an input-output curve for DNa02 is plotted above each distribution comparison.
%
% INPUTS:
%   nK                - Number of steering gain values (k) in the experiment.
%   kValues           - Vector of k values used in the experiment.
%   inputDist         - 3D array of input history distributions (k x bins x conditions).
%   bins              - Vector defining the bins for input values.
%   comparisonLabel   - Cell array with labels for all conditions being compared.
%   runSettings       - Struct containing the parameters, including DNa02input and shift.
%   folder            - Folder directories for saving figures.
%
% OUTPUTS:
%   Generates a figure with plots showing the DNa02 input-output curves and 
%   distributions for each k value.
%
% CREATED: 11/05/2024 - MC
% UPDATED: 11/16/2024 - MC refactored for multiple conditions
%

function generateSummaryInputDistributions(nK, kValues, inputDist, bins, comparisonLabel, runSettings, folder)

    % Motor input-output for DNa02
    DNa02input = runSettings.DNa02input;  % Input range for downstream neurons
    alpha = runSettings.alpha;           % Specify alpha based on model requirements
    DNa02output = adjELU(DNa02input, alpha, runSettings.shift); % Output range for downstream neurons

    % Generate summary plot for DNa02 Input-Output Curve and Input History Distributions
    figure;
    set(gcf, 'Position', [100 100 1400 600]);  % Set figure size
    tiledlayout(2, nK, 'TileSpacing', 'compact');  % Two rows for each k value
    xrange = [runSettings.minIn, runSettings.maxIn];
    numConditions = size(inputDist, 3);  % Determine number of conditions

    % Plot DNa02 input-output curve in the first row
    for kIdx = 1:nK
        nexttile(kIdx); hold on;
        
        % Plot DNa02 input-output curve
        plot(DNa02input, DNa02output, 'DisplayName', 'DNa02 Input-Output Curve');

        % Customize the plot
        xlabel('Input Value');
        ylabel('Output Value');
        title(['k=' num2str(kValues(kIdx))]);
        xline(0, 'k');
        axis tight;
        xlim(xrange);
        grid on;
    end

    % Plot input history distributions for all conditions in the second row
    for kIdx = 1:nK
        nexttile(kIdx + nK); hold on;
        
        % Plot the distributions for all conditions
        for condIdx = 1:numConditions
            plot(bins, inputDist(kIdx, :, condIdx), '-', 'DisplayName', comparisonLabel{condIdx});
        end

        % Customize the plot
        xlabel('Input Value Bins');
        ylabel('Normalized Frequency');
        if kIdx == 1
            legend('show', 'Location', 'northwest');
        end
        grid on;
        xline(0, 'k');
        xlim(xrange);
        ylim([0 1]);  % Set y-axis to normalized range
    end

    % Overall title for the figure
    sgtitle(['DNa02 Input-Output and Input History Distribution Comparison Across Conditions (' ...
             strjoin(comparisonLabel, ' vs ') ')']);

    % Save the plot in both PNG and SVG formats
    saveas(gcf, fullfile(folder.final, [strjoin(comparisonLabel, 'v') '_Input_Distributions_with_DNa02_Curve.png']));
    set(gcf, 'renderer', 'Painters');  % Ensure high-resolution SVG output by setting renderer
    saveas(gcf, fullfile(folder.vectors, [strjoin(comparisonLabel, 'v') '_Input_Distributions_with_DNa02_Curve.svg']));

end
