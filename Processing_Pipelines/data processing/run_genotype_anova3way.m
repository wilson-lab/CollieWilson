% run_genotype_anova3way
% This function performs a three-way ANOVA on a dependent variable based on genotype, condition, and animal as a factor.
% It includes post-hoc Tukey-adjusted pairwise comparisons for main effects (genotype and condition).
%
% INPUTS:
%   dependentVar_kir   - Matrix of the dependent variable for KIR genotype (animals x conditions)
%   dependentVar_wt    - Matrix of the dependent variable for WT genotype (animals x conditions)
%   dependentVar_na    - Matrix of the dependent variable for NA genotype (animals x conditions)
%   dependent_var_name - Name of the dependent variable (e.g., 'Lag Time')
%   folder             - Structure containing folder paths for saving files
%
% OUTPUTS:
%   p     - p-values from the three-way ANOVA for genotype and condition
%   tbl   - ANOVA table
%
% CREATED: 11/08/2024 - MC adapted for three-way ANOVA with Tukey-adjusted post-hoc comparisons

function [p, tbl] = run_genotype_anova3way(dependentVar_kir, dependentVar_wt, dependentVar_na, dependent_var_name, folder)

% Replace any '-' or ' ' in dependent_var_name with '_'
dependent_var_name = strrep(dependent_var_name, '-', '_');
dependent_var_name = strrep(dependent_var_name, ' ', '_');

% Define condition names for table
nConditions = size(dependentVar_kir, 2);
conditionNames = strcat('Cond', string(1:nConditions));

% Convert data into long format for ANOVA
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

% Perform three-way ANOVA for main effects only
[p, tbl, stats] = anovan(long_data.Response, {long_data.Genotype, long_data.Condition, long_data.AnimalID}, ...
    'model', 'linear', 'varnames', {'Genotype', 'Condition', 'AnimalID'}, 'display', 'off');

% Extract p-values for main effects directly from the p output
pGenotype = p(1);    % Genotype
pCondition = p(2);   % Condition
pAnimalID = p(3);    % AnimalID
p = [pGenotype, pCondition, pAnimalID];

% Perform Tukey HSD post-hoc tests for main effects only
genotype_comparison_tbl = multcompare(stats, 'Dimension', 1, 'CType', 'tukey-kramer', 'Display', 'off');
condition_comparison_tbl = multcompare(stats, 'Dimension', 2, 'CType', 'tukey-kramer', 'Display', 'off');

% Convert comparison results to tables
genotype_comparison_tbl = array2table(genotype_comparison_tbl, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'pValue'});
condition_comparison_tbl = array2table(condition_comparison_tbl, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'pValue'});

% Save ANOVA table and comparisons to Excel
anova_xlsx_file = fullfile(folder.summary, ['anova_' dependent_var_name '.xlsx']);
writetable(cell2table(tbl), anova_xlsx_file, 'Sheet', 'ANOVA');
writetable(genotype_comparison_tbl, anova_xlsx_file, 'Sheet', 'Genotype Comparisons');
writetable(condition_comparison_tbl, anova_xlsx_file, 'Sheet', 'Condition Comparisons');

disp(['ANOVA and Pairwise Comparison results saved to: ' anova_xlsx_file]);
end
