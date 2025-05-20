function plotOscSlopes(kir_oscAngRL, wt_oscAngRL, na_oscAngRL, oscSweepPos, folder)
    %% Initialize
    % Combine inputs into a cell array for processing
    genotypes = {'Kir', 'WT', 'Na'};
    dataArrays = {kir_oscAngRL, wt_oscAngRL, na_oscAngRL};
    
    % Check if dimensions match for all arrays
    for i = 1:length(dataArrays)
        if size(oscSweepPos, 1) ~= size(dataArrays{i}, 1)
            error('Position data and response data dimensions must align.');
        end
    end
    
%% Find Index Range for Middle to Closest to +30
midPointIdx = round(length(oscSweepPos) / 2); % Find the midpoint index

% Find the index of the position closest to +30 after the midpoint
[~, maxIdx] = min(abs(oscSweepPos(midPointIdx:end) - 30)); 
maxIdx = midPointIdx + maxIdx - 1; % Adjust to get the actual index in oscSweepPos

if isempty(maxIdx) || maxIdx <= midPointIdx
    error('No valid indices found from the middle to the closest position to +30.');
end

idxRange = midPointIdx:maxIdx; % Set the range from midpoint to closest to +30
    
    %% Fit Slopes
    % Initialize output for slopes
    slopes = cell(1, length(dataArrays));
    medians = zeros(1, length(dataArrays));
    allSlopes = [];
    groupLabels = [];

    % Loop over genotypes
    for g = 1:length(dataArrays)
        % Extract data for this genotype
        data = dataArrays{g};
        nAnimals = size(data, 2);

        % Initialize slopes for this genotype
        slopes{g} = nan(1, nAnimals);

        % Fit a linear slope for each animal
        for a = 1:nAnimals
            x = oscSweepPos(idxRange);
            y = data(idxRange, a);
            if ~any(isnan(y)) % Check for valid data
                p = polyfit(x, y, 1); % Linear fit
                slopes{g}(a) = p(1); % Store the slope
            end
        end

        % Append slopes and labels for ANOVA
        allSlopes = [allSlopes, slopes{g}];
        groupLabels = [groupLabels, repelem(g, nAnimals)];

        % Calculate the median slope for this genotype
        medians(g) = median(slopes{g}, 'omitnan');
    end
    
    %% Perform ANOVA
    [p, tbl, stats] = anova1(allSlopes, groupLabels, 'off'); % Perform one-way ANOVA
    genotypeComparison = multcompare(stats, 'Display', 'off'); % Tukey's post-hoc test
    
    % Convert ANOVA table to cell and save to Excel
    anovaTbl = cell2table(tbl(2:end,:), 'VariableNames', tbl(1,:));
    genotypeComparisonTbl = array2table(genotypeComparison, ...
        'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'pValue'});
    
    % Define the Excel file path and save
    anovaXlsxFile = fullfile(folder.summary, 'ANOVA_oscslope.xlsx');
    writetable(anovaTbl, anovaXlsxFile, 'Sheet', 'ANOVA');
    writetable(genotypeComparisonTbl, anovaXlsxFile, 'Sheet', 'Genotype Comparisons');
    
    %% Plot
    % Generate the plot
    figure;
    hold on;
    jitterAmount = 0.1; % Define the jitter amount

    for g = 1:length(dataArrays)
        % Generate jitter for x-coordinates
        jitter = (rand(1, numel(slopes{g})) - 0.5) * jitterAmount;

        % Plot slopes for each animal as gray dots with jitter
        scatter(g + jitter, slopes{g}, '.', 'MarkerEdgeColor', [0.5, 0.5, 0.5]);

        % Plot median for this genotype as a line
        plot([g - 0.2, g + 0.2], [medians(g), medians(g)], 'k-', 'LineWidth', 2);
    end

    % Add p-value to plot
    text(2.5, max(allSlopes) * 0.9, ['p(genotype) = ' num2str(p, '%.3g')], ...
        'HorizontalAlignment', 'right', 'FontSize', 12);

    % Customize plot
    ylim([0 8])
    xlim([0.5, length(dataArrays) + 0.5]);
    xticks(1:length(dataArrays));
    xticklabels(genotypes); % Use genotype names
    ylabel('Linear Slope');
    xlabel('Genotype');
    title('Slopes for Each Genotype');
    hold off;

end