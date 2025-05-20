% run_genotype_anova_repeated
% This function performs a mixed-effects model analysis on a dependent variable based on genotype and condition.
% It handles NaNs in the data by removing them before processing. After the model, it runs post-hoc
% Tukey-adjusted pairwise comparisons on both genotype and condition factors. The function saves the model table 
% and comparison results to a single Excel file with separate sheets.
%
% INPUTS:
%   dependentVar_kir   - Matrix of the dependent variable for KIR genotype (animals x conditions)
%   dependentVar_wt    - Matrix of the dependent variable for WT genotype (animals x conditions)
%   dependentVar_na    - Matrix of the dependent variable for NA genotype (animals x conditions)
%   dependent_var_name - Name of the dependent variable (e.g., 'Lag Time')
%   folder             - Structure containing folder paths for saving files
%
% OUTPUTS:
%   p     - p-values from the mixed-effects model for genotype and condition
%   tbl   - Model table
%
% CREATED: 11/01/2024 - MC adapted for mixed-effects model using fitlme with Tukey-adjusted post-hoc comparisons
% UPDATED: 01/07/2025 - MC added posthoc for interaction term
%
function [p, tbl] = run_genotype_anova_repeated(dependentVar_kir, dependentVar_wt, dependentVar_na, dependent_var_name, folder)

% Replace any '-' or ' ' in dependent_var_name with '_'
dependent_var_name = strrep(dependent_var_name, '-', '_');
dependent_var_name = strrep(dependent_var_name, ' ', '_');

% Define condition names for table
nConditions = size(dependentVar_kir, 2);
conditionNames = strcat('Cond', string(1:nConditions));

% Convert data into long format for mixed-effects model
kir_data = array2table(dependentVar_kir, 'VariableNames', conditionNames);
wt_data = array2table(dependentVar_wt, 'VariableNames', conditionNames);
na_data = array2table(dependentVar_na, 'VariableNames', conditionNames);

% Add genotype and animal identifiers
kir_data.Genotype = categorical(repmat({'KIR'}, size(kir_data, 1), 1));
wt_data.Genotype = categorical(repmat({'WT'}, size(wt_data, 1), 1));
na_data.Genotype = categorical(repmat({'NA'}, size(na_data, 1), 1));
kir_data.AnimalID = (1:size(kir_data, 1))';
wt_data.AnimalID = (1:size(wt_data, 1))';
na_data.AnimalID = (1:size(na_data, 1))';

% Stack tables and reshape to long format
all_data = [kir_data; wt_data; na_data];
long_data = stack(all_data, conditionNames, 'NewDataVariableName', 'Response', 'IndexVariableName', 'Condition');

% Remove rows with any NaNs in the Response variable
long_data = rmmissing(long_data);

% Fit a linear mixed-effects model with random intercept for AnimalID
lme = fitlme(long_data, 'Response ~ Genotype * Condition + (1|AnimalID)');

% Display the results
tbl = anova(lme, 'DFMethod', 'satterthwaite');  % ANOVA table with p-values
disp([dependent_var_name ' Mixed-Effects Model Results:']);
disp(tbl);

% Extract p-values for Genotype and Condition
pGenotype = tbl.pValue(2);    % For Genotype
pCondition = tbl.pValue(3);    % For Condition
pInteraction = tbl.pValue(4);  % For Genotype * Condition interaction
p = [pGenotype, pCondition, pInteraction];

% Define Tukey adjustment function
tukey_adjust = @(p_val, n_comparisons) 1 - (1 - p_val).^n_comparisons;

% Define levels of Genotype and Condition for comparisons
genotype_levels = unique(long_data.Genotype);
condition_levels = unique(long_data.Condition);

% Calculate predicted means for each level of Genotype
genotype_comparisons = {};
n_genotype_comparisons = nchoosek(length(genotype_levels), 2); % Number of comparisons
genotype_predictions = zeros(length(genotype_levels), 1);

for i = 1:length(genotype_levels)
    % Create a temporary dataset for prediction with each Genotype level
    temp_data = long_data(1, :);
    temp_data.Genotype(:) = genotype_levels(i);
    temp_data.Condition(:) = condition_levels(1); % Hold Condition constant
    
    % Predict conditional mean for this level of Genotype
    genotype_predictions(i) = predict(lme, temp_data);
end

% Perform Tukey-adjusted pairwise comparisons for Genotype
for i = 1:length(genotype_levels)
    for j = i+1:length(genotype_levels)
        % Calculate difference, standard error, t-statistic, and adjusted p-value
        diff_value = genotype_predictions(i) - genotype_predictions(j);
        se_diff = sqrt(lme.CoefficientCovariance(i, i) + lme.CoefficientCovariance(j, j) - 2 * lme.CoefficientCovariance(i, j));
        t_stat = diff_value / se_diff;
        p_value = 2 * (1 - tcdf(abs(t_stat), lme.DFE));
        p_adj = tukey_adjust(p_value, n_genotype_comparisons);
        
        % Store results
        genotype_comparisons = [genotype_comparisons; {genotype_levels(i), genotype_levels(j), diff_value, se_diff, t_stat, p_adj}];
    end
end
genotype_comparison_tbl = cell2table(genotype_comparisons, 'VariableNames', {'Group1', 'Group2', 'Difference', 'SE', 'TStatistic', 'AdjPValue'});

% Calculate predicted means for each level of Condition
condition_comparisons = {};
n_condition_comparisons = nchoosek(length(condition_levels), 2); % Number of comparisons
condition_predictions = zeros(length(condition_levels), 1);

for i = 1:length(condition_levels)
    % Create a temporary dataset for prediction with each Condition level
    temp_data = long_data(1, :);
    temp_data.Condition(:) = condition_levels(i);
    temp_data.Genotype(:) = genotype_levels(1); % Hold Genotype constant
    
    % Predict conditional mean for this level of Condition
    condition_predictions(i) = predict(lme, temp_data);
end

% Perform Tukey-adjusted pairwise comparisons for Condition
for i = 1:length(condition_levels)
    for j = i+1:length(condition_levels)
        % Calculate difference, standard error, t-statistic, and adjusted p-value
        diff_value = condition_predictions(i) - condition_predictions(j);
        se_diff = sqrt(lme.CoefficientCovariance(i, i) + lme.CoefficientCovariance(j, j) - 2 * lme.CoefficientCovariance(i, j));
        t_stat = diff_value / se_diff;
        p_value = 2 * (1 - tcdf(abs(t_stat), lme.DFE));
        p_adj = tukey_adjust(p_value, n_condition_comparisons);
        
        % Store results
        condition_comparisons = [condition_comparisons; {condition_levels(i), condition_levels(j), diff_value, se_diff, t_stat, p_adj}];
    end
end
condition_comparison_tbl = cell2table(condition_comparisons, 'VariableNames', {'Group1', 'Group2', 'Difference', 'SE', 'TStatistic', 'AdjPValue'});

% Post-hoc tests for interaction term
interaction_comparisons = {};
n_interaction_comparisons = nchoosek(length(genotype_levels) * length(condition_levels), 2);

% Generate predicted means for each genotype-condition combination
interaction_predictions = zeros(length(genotype_levels), length(condition_levels));
for g = 1:length(genotype_levels)
    for c = 1:length(condition_levels)
        temp_data = long_data(1, :);
        temp_data.Genotype(:) = genotype_levels(g);
        temp_data.Condition(:) = condition_levels(c);
        interaction_predictions(g, c) = predict(lme, temp_data);
    end
end

% Perform pairwise comparisons for interaction term
for g1 = 1:length(genotype_levels)
    for c1 = 1:length(condition_levels)
        for g2 = g1:length(genotype_levels)
            for c2 = (g1 == g2) * (c1 + 1) + (g1 ~= g2) * 1:length(condition_levels)
                % Calculate difference, standard error, t-statistic, and adjusted p-value
                diff_value = interaction_predictions(g1, c1) - interaction_predictions(g2, c2);
                se_diff = sqrt(lme.CoefficientCovariance(g1, g1) + lme.CoefficientCovariance(g2, g2) - ...
                               2 * lme.CoefficientCovariance(g1, g2));
                t_stat = diff_value / se_diff;
                p_value = 2 * (1 - tcdf(abs(t_stat), lme.DFE));
                p_adj = 1 - (1 - p_value).^n_interaction_comparisons; % Tukey adjustment
                
                % Store results
                interaction_comparisons = [interaction_comparisons; ...
                    {genotype_levels(g1), condition_levels(c1), ...
                     genotype_levels(g2), condition_levels(c2), ...
                     diff_value, se_diff, t_stat, p_adj}];
            end
        end
    end
end

% Convert interaction comparisons to table
interaction_comparison_tbl = cell2table(interaction_comparisons, ...
    'VariableNames', {'Genotype1', 'Condition1', 'Genotype2', 'Condition2', 'Difference', 'SE', 'TStatistic', 'AdjPValue'});

% Save model and comparison results to a single Excel file with separate sheets
model_tbl_xlsx_file = fullfile(folder.summary, ['anova_' dependent_var_name '.xlsx']);
writetable(dataset2table(tbl), model_tbl_xlsx_file, 'Sheet', 'Model');
writetable(genotype_comparison_tbl, model_tbl_xlsx_file, 'Sheet', 'Genotype Comparisons');
writetable(condition_comparison_tbl, model_tbl_xlsx_file, 'Sheet', 'Condition Comparisons');
writetable(interaction_comparison_tbl, model_tbl_xlsx_file, 'Sheet', 'Interaction Comparisons');

disp(['Model and Pairwise Comparison results saved to: ' model_tbl_xlsx_file]);
end
