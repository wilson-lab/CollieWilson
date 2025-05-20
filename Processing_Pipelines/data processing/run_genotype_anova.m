% run_genotype_anova
% This function performs a two-way ANOVA on a dependent variable based on genotype and condition.
% It handles NaNs in the data by removing them before processing. After ANOVA, it runs post-hoc
% Tukey-Kramer tests on both genotype and condition factors. The function saves the ANOVA table 
% and Tukey-Kramer results to a single Excel file with separate sheets.
%
% INPUTS:
%   dependentVar_kir   - Matrix of the dependent variable for KIR genotype (animals x conditions)
%   dependentVar_wt    - Matrix of the dependent variable for WT genotype (animals x conditions)
%   dependentVar_na    - Matrix of the dependent variable for NA genotype (animals x conditions)
%   dependent_var_name - Name of the dependent variable (e.g., 'Lag Time')
%   folder             - Structure containing folder paths for saving files
%
% OUTPUTS:
%   p     - p-values from the two-way ANOVA for genotype and condition
%   tbl   - ANOVA table
%
% CREATED: 10/??/2024 - MC
%
function [p, tbl] = run_genotype_anova(dependentVar_kir, dependentVar_wt, dependentVar_na, dependent_var_name, folder)
% Replace any '-' or ' ' in dependent_var_name with '_'
dependent_var_name = strrep(dependent_var_name, '-', '_');
dependent_var_name = strrep(dependent_var_name, ' ', '_');

% Initialize variables to store valid data and labels
all_data = [];
genotype_labels = [];
condition_labels = [];

% Helper function to process each genotype array, handling NaNs
    function process_genotype_data(dep_var, genotype_label)
        [nAnimals, nConditions] = size(dep_var);

        for iCondition = 1:nConditions
            valid_idx = ~isnan(dep_var(:, iCondition));  % Find non-NaN entries for each condition
            current_data = dep_var(valid_idx, iCondition);  % Extract valid data

            % Append valid data
            all_data = [all_data; current_data];

            % Append genotype and condition labels
            genotype_labels = [genotype_labels; repmat(genotype_label, sum(valid_idx), 1)];  % Genotype label
            condition_labels = [condition_labels; repmat(iCondition, sum(valid_idx), 1)];  % Condition label
        end
    end

% Process data for each genotype
process_genotype_data(dependentVar_kir, 1);  % Process KIR data
process_genotype_data(dependentVar_wt, 2);  % Process WT data
process_genotype_data(dependentVar_na, 3);  % Process NA data

% Ensure that the data and labels have matching lengths
if length(all_data) ~= length(genotype_labels) || length(all_data) ~= length(condition_labels)
    error('Mismatch between the length of data and the labels.');
end

% Run two-way ANOVA with genotype and condition as factors
[p, tbl, stats] = anovan(all_data, {genotype_labels, condition_labels}, 'model', 'interaction', ...
    'varnames', {'Genotype', 'Gain'}, 'display', 'off');
% Return p-values
p = p(1:2);  % Return p-values for Genotype and Condition

% Display ANOVA results in the command window
disp([dependent_var_name ' ANOVA Results:']);
disp(tbl);

% Run Tukey-Kramer post-hoc test for Genotype with display off
[c_genotype, ~, ~, ~] = multcompare(stats, 'Dimension', 1, 'ctype', 'tukey-kramer', 'display', 'off');

% Run Tukey-Kramer post-hoc test for Condition with display off
[c_gain, ~, ~, ~] = multcompare(stats, 'Dimension', 2, 'ctype', 'tukey-kramer', 'display', 'off');

% Prepare Tukey-Kramer results for saving
tukey_genotype_tbl = array2table(c_genotype, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'PValue'});
tukey_gain_tbl = array2table(c_gain, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'PValue'});

% Save ANOVA and Tukey-Kramer results to a single Excel file with separate sheets
anova_tbl_xlsx_file = fullfile(folder.summary, ['anova_' dependent_var_name '.xlsx']);
writetable(cell2table(tbl), anova_tbl_xlsx_file, 'Sheet', 'ANOVA');
writetable(tukey_genotype_tbl, anova_tbl_xlsx_file, 'Sheet', 'Tukey-Kramer-Genotype');
writetable(tukey_gain_tbl, anova_tbl_xlsx_file, 'Sheet', 'Tukey-Kramer-Gain');

disp(['ANOVA and Post-Hoc Tukey-Kramer results saved to: ' anova_tbl_xlsx_file]);
end
