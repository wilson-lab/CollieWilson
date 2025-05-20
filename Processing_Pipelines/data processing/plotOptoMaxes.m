function plotOptoMaxes(kir_optoRL, wt_optoRL, na_optoRL, folder)
    %% Initialize
    % Combine inputs into a cell array for processing
    genotypes = {'Kir', 'WT', 'Na'};
    dataArrays = {kir_optoRL, wt_optoRL, na_optoRL};
    
    % Initialize output for maxima
    maxResponses = cell(1, length(dataArrays));
    medians = zeros(1, length(dataArrays));
    allMaxes = [];
    groupLabels = [];

    %% Find Maximum Responses
    % Loop over genotypes
    for g = 1:length(dataArrays)
        % Extract data for this genotype
        data = dataArrays{g};
        nAnimals = size(data, 2);

        % Initialize maxima for this genotype
        maxResponses{g} = nan(1, nAnimals);

        % Find the maximum response for each animal
        for a = 1:nAnimals
            y = data(:, a);
            if ~any(isnan(y)) % Check for valid data
                maxResponses{g}(a) = max(y); % Store the maximum response
            end
        end

        % Append maxima and labels for ANOVA
        allMaxes = [allMaxes, maxResponses{g}];
        groupLabels = [groupLabels, repelem(g, nAnimals)];

        % Calculate the median maximum response for this genotype
        medians(g) = median(maxResponses{g}, 'omitnan');
    end

    %% Perform ANOVA
    [p, tbl, stats] = anova1(allMaxes, groupLabels, 'off'); % Perform one-way ANOVA
    genotypeComparison = multcompare(stats, 'Display', 'off'); % Tukey's post-hoc test
    
    % Convert ANOVA table to cell and save to Excel
    anovaTbl = cell2table(tbl(2:end,:), 'VariableNames', tbl(1,:));
    genotypeComparisonTbl = array2table(genotypeComparison, ...
        'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'pValue'});
    
    % Define the Excel file path and save
    anovaXlsxFile = fullfile(folder.summary, 'ANOVA_optoMax.xlsx');
    writetable(anovaTbl, anovaXlsxFile, 'Sheet', 'ANOVA');
    writetable(genotypeComparisonTbl, anovaXlsxFile, 'Sheet', 'Genotype Comparisons');

%% Plot
% Generate the plot
figure;
hold on;
jitterAmount = 0.1; % Define the jitter amount

for g = 1:length(dataArrays)
    % Generate jitter for x-coordinates
    jitter = (rand(1, numel(maxResponses{g})) - 0.5) * jitterAmount;

    % Plot maxima for each animal as gray dots with jitter
    scatter(g + jitter, maxResponses{g}, '.', 'MarkerEdgeColor', [0.5, 0.5, 0.5]);

    % Plot median for this genotype as a dashed line
    plot([g - 0.2, g + 0.2], [medians(g), medians(g)], 'k', 'LineWidth', 2);
end

% Add p-value to plot
text(2.5, max(allMaxes) * 0.9, ['p(genotype) = ' num2str(p, '%.3g')], ...
    'HorizontalAlignment', 'right', 'FontSize', 12);

% Customize plot
xlim([0.5, length(dataArrays) + 0.5]);
xticks(1:length(dataArrays));
xticklabels(genotypes); % Use genotype names
ylabel('Maximum Turn Response');
xlabel('Genotype');
title('Maximum Turn Responses for Each Genotype');
hold off;

end
