% run_genotype_anova3_repeated
% This function performs a three-way mixed-effects model analysis on a dependent variable based on genotype,
% angular velocity bins (continuous), and gain conditions (categorical). It handles NaNs in the data by removing
% them before processing. After the model, it runs post-hoc Tukey-adjusted pairwise comparisons for genotype,
% gain, and bins. The function saves the model table and comparison results to a single Excel file.
%
% INPUTS:
%   dependentVar_kir   - Matrix of the dependent variable for KIR genotype (bins x gains x animals)
%   dependentVar_wt    - Matrix of the dependent variable for WT genotype (bins x gains x animals)
%   dependentVar_na    - Matrix of the dependent variable for NA genotype (bins x gains x animals)
%   bins               - Array of angular velocity bins (1 x N)
%   dependent_var_name - Name of the dependent variable (e.g., 'Lag Time')
%   folder             - Structure containing folder paths for saving files
%
% OUTPUTS:
%   p     - p-values for the mixed-effects model (Genotype, Gain, Bins, and their interactions)
%   tbl   - Model table
%
% CREATED: 11/2024 MC adapted for three-way mixed-effects model using fitlme with Tukey-adjusted post-hoc comparisons
%
function [p, tbl] = run_genotype_anova3_repeated(dependentVar_kir, dependentVar_wt, dependentVar_na, bins, dependent_var_name, folder)

    % Replace any '-' or ' ' in dependent_var_name with '_'
    dependent_var_name = strrep(dependent_var_name, '-', '_');
    dependent_var_name = strrep(dependent_var_name, ' ', '_');

    % Initialize variables to store valid data and labels
    all_data = [];
    genotype_labels = categorical();  % Initialize as an empty categorical array
    gain_labels = [];
    bin_labels = [];
    animal_labels = [];

    % Function to process each genotype's data and append labels
    function process_genotype_data(dep_var, genotype_label)
        [nBins, nGains, nAnimals] = size(dep_var);
        for iAnimal = 1:nAnimals
            for iGain = 1:nGains
                for iBin = 1:nBins
                    current_data = dep_var(iBin, iGain, iAnimal);
                    % Append valid data, ignoring NaNs
                    if ~isnan(current_data)
                        all_data = [all_data; current_data];
                        genotype_labels = [genotype_labels; genotype_label];
                        gain_labels = [gain_labels; iGain];
                        bin_labels = [bin_labels; bins(iBin)];
                        animal_labels = [animal_labels; iAnimal];
                    end
                end
            end
        end
    end

    % Process data for each genotype
    process_genotype_data(dependentVar_kir, categorical({'KIR'}));
    process_genotype_data(dependentVar_wt, categorical({'WT'}));
    process_genotype_data(dependentVar_na, categorical({'NA'}));

    % Create a table for the mixed-effects model
    data_tbl = table(all_data, genotype_labels, categorical(gain_labels), bin_labels, categorical(animal_labels), ...
                     'VariableNames', {'Response', 'Genotype', 'Gain', 'Bins', 'AnimalID'});

    % Fit a linear mixed-effects model with random intercept for AnimalID
    lme = fitlme(data_tbl, 'Response ~ Genotype * Gain * Bins + (1|AnimalID)');

    % Display the results
    tbl = anova(lme, 'DFMethod', 'satterthwaite');  % ANOVA table with p-values
    disp([dependent_var_name ' Three-Way Mixed-Effects Model Results:']);
    disp(tbl);

    % Extract p-values for Genotype, Gain, Bins, and their interactions
    pGenotype = tbl.pValue(2);        % For Genotype
    pGain = tbl.pValue(3);             % For Gain
    pBins = tbl.pValue(4);             % For Bins
    pInteraction = tbl.pValue(5:end);  % For interactions
    p = [pGenotype, pGain, pBins, pInteraction'];

    % Define Tukey adjustment function
    tukey_adjust = @(p_val, n_comparisons) 1 - (1 - p_val).^n_comparisons;

    % Levels of Genotype, Gain, and Bins for post-hoc comparisons
    genotype_levels = unique(data_tbl.Genotype);
    gain_levels = unique(data_tbl.Gain);
    bin_levels = unique(data_tbl.Bins);

    % Tukey-adjusted post-hoc comparisons for Genotype
    n_genotype_comparisons = nchoosek(length(genotype_levels), 2);
    genotype_comparisons = perform_tukey_posthoc(lme, genotype_levels, n_genotype_comparisons, 'Genotype');

    % Tukey-adjusted post-hoc comparisons for Gain
    n_gain_comparisons = nchoosek(length(gain_levels), 2);
    gain_comparisons = perform_tukey_posthoc(lme, gain_levels, n_gain_comparisons, 'Gain');

    % Tukey-adjusted post-hoc comparisons for Bins
    n_bin_comparisons = nchoosek(length(bin_levels), 2);
    bin_comparisons = perform_tukey_posthoc(lme, bin_levels, n_bin_comparisons, 'Bins');

    % Convert comparison results to tables for saving
    tukey_genotype_tbl = array2table(genotype_comparisons, 'VariableNames', {'Group1', 'Group2', 'Difference', 'SE', 'TStatistic', 'AdjPValue'});
    tukey_gain_tbl = array2table(gain_comparisons, 'VariableNames', {'Group1', 'Group2', 'Difference', 'SE', 'TStatistic', 'AdjPValue'});
    tukey_bins_tbl = array2table(bin_comparisons, 'VariableNames', {'Group1', 'Group2', 'Difference', 'SE', 'TStatistic', 'AdjPValue'});

    % Save model and comparison results to an Excel file
    anova_tbl_xlsx_file = fullfile(folder.summary, ['anova_' dependent_var_name '.xlsx']);
    writetable(dataset2table(tbl), anova_tbl_xlsx_file, 'Sheet', 'Model');
    writetable(tukey_genotype_tbl, anova_tbl_xlsx_file, 'Sheet', 'Tukey-Genotype');
    writetable(tukey_gain_tbl, anova_tbl_xlsx_file, 'Sheet', 'Tukey-Gain');
    writetable(tukey_bins_tbl, anova_tbl_xlsx_file, 'Sheet', 'Tukey-Bins');

    disp(['Model and Pairwise Comparison results saved to: ' anova_tbl_xlsx_file]);
end

% Helper function for Tukey-adjusted post-hoc tests
function comparisons = perform_tukey_posthoc(lme, levels, n_comparisons, factor_name)
    comparisons = {};
    predictions = zeros(length(levels), 1);
    for i = 1:length(levels)
        temp_data = lme.Variables(1, :);
        temp_data.(factor_name)(:) = levels(i);
        predictions(i) = predict(lme, temp_data);
    end
    for i = 1:length(levels)
        for j = i+1:length(levels)
            diff_value = predictions(i) - predictions(j);
            se_diff = sqrt(lme.CoefficientCovariance(i, i) + lme.CoefficientCovariance(j, j) - 2 * lme.CoefficientCovariance(i, j));
            t_stat = diff_value / se_diff;
            p_value = 2 * (1 - tcdf(abs(t_stat), lme.DFE));
            p_adj = 1 - (1 - p_value).^n_comparisons;
            comparisons = [comparisons; {levels(i), levels(j), diff_value, se_diff, t_stat, p_adj}];
        end
    end
end
