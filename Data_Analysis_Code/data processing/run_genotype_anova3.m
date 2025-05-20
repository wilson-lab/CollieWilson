% run_genotype_anova3
% This function performs a three-way ANOVA on a dependent variable based on genotype,
% angular velocity bins (continuous), and gain conditions (categorical).
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
%   p     - p-values for the ANOVA (Genotype, Gain, Bins, and their interactions)
%   tbl   - ANOVA table
%
% CREATED: 10/2024 MC
%
function [p, tbl] = run_genotype_anova3(dependentVar_kir, dependentVar_wt, dependentVar_na, bins, dependent_var_name, folder)
    % Replace any '-' or ' ' in dependent_var_name with '_'
    dependent_var_name = strrep(dependent_var_name, '-', '_');
    dependent_var_name = strrep(dependent_var_name, ' ', '_');
    
    % Initialize variables to store valid data and labels
    all_data = [];
    genotype_labels = [];
    gain_labels = [];
    bin_labels = [];

    % Process KIR genotype data
    [nBins_kir, nGains_kir, nAnimals_kir] = size(dependentVar_kir);
    for iAnimal = 1:nAnimals_kir
        for iGain = 1:nGains_kir
            for iBin = 1:nBins_kir
                current_data = dependentVar_kir(iBin, iGain, iAnimal);  % Extract data for current bin, gain, and animal
                
                % Append valid data (ignoring NaNs)
                if ~isnan(current_data)
                    all_data = [all_data; current_data];  
                    genotype_labels = [genotype_labels; 1];  % Genotype label for KIR
                    gain_labels = [gain_labels; iGain];  % Gain label
                    bin_labels = [bin_labels; bins(iBin)];  % Bin label
                end
            end
        end
    end

    % Process WT genotype data
    [nBins_wt, nGains_wt, nAnimals_wt] = size(dependentVar_wt);
    for iAnimal = 1:nAnimals_wt
        for iGain = 1:nGains_wt
            for iBin = 1:nBins_wt
                current_data = dependentVar_wt(iBin, iGain, iAnimal);  % Extract data for current bin, gain, and animal
                
                % Append valid data (ignoring NaNs)
                if ~isnan(current_data)
                    all_data = [all_data; current_data];  
                    genotype_labels = [genotype_labels; 2];  % Genotype label for WT
                    gain_labels = [gain_labels; iGain];  % Gain label
                    bin_labels = [bin_labels; bins(iBin)];  % Bin label
                end
            end
        end
    end

    % Process NA genotype data
    [nBins_na, nGains_na, nAnimals_na] = size(dependentVar_na);
    for iAnimal = 1:nAnimals_na
        for iGain = 1:nGains_na
            for iBin = 1:nBins_na
                current_data = dependentVar_na(iBin, iGain, iAnimal);  % Extract data for current bin, gain, and animal
                
                % Append valid data (ignoring NaNs)
                if ~isnan(current_data)
                    all_data = [all_data; current_data];  
                    genotype_labels = [genotype_labels; 3];  % Genotype label for NA
                    gain_labels = [gain_labels; iGain];  % Gain label
                    bin_labels = [bin_labels; bins(iBin)];  % Bin label
                end
            end
        end
    end

    % Ensure that the data and labels have matching lengths
    if length(all_data) ~= length(genotype_labels) || length(all_data) ~= length(gain_labels) || length(all_data) ~= length(bin_labels)
        error('Mismatch between the length of data and the labels.');
    end

    % Run three-way ANOVA with genotype (categorical), bins (categorical), and gain (categorical) as factors
    [p, tbl, stats] = anovan(all_data, {genotype_labels, gain_labels, bin_labels}, 'model', 'full', ...
                             'varnames', {'Genotype', 'Gain', 'Bins'}, 'display', 'off');
    
    % Display ANOVA results in the command window
    disp([dependent_var_name ' ANOVA Results:']);
    disp(tbl);

    % Return p-values for Genotype, Gain, and Bins
    p = p(1:4);  % p-values for Genotype, Gain, Bins, and their interactions
    
    % Run Tukey-Kramer post-hoc tests
    % Post-hoc for Genotype
    [c_genotype, ~, ~, ~] = multcompare(stats, 'Dimension', 1, 'ctype', 'tukey-kramer', 'display', 'off');

    % Post-hoc for Gain
    [c_gain, ~, ~, ~] = multcompare(stats, 'Dimension', 2, 'ctype', 'tukey-kramer', 'display', 'off');

    % Post-hoc for Bins
    [c_bins, ~, ~, ~] = multcompare(stats, 'Dimension', 3, 'ctype', 'tukey-kramer', 'display', 'off');

    % Prepare Tukey-Kramer results for saving
    tukey_genotype_tbl = array2table(c_genotype, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'PValue'});
    tukey_gain_tbl = array2table(c_gain, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'PValue'});
    tukey_bins_tbl = array2table(c_bins, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'PValue'});

    % Save ANOVA and Tukey-Kramer results to a single Excel file with separate sheets
    anova_tbl_xlsx_file = fullfile(folder.summary, ['anova_' dependent_var_name '.xlsx']);
    writetable(cell2table(tbl), anova_tbl_xlsx_file, 'Sheet', 'ANOVA');
    writetable(tukey_genotype_tbl, anova_tbl_xlsx_file, 'Sheet', 'Tukey-Kramer-Genotype');
    writetable(tukey_gain_tbl, anova_tbl_xlsx_file, 'Sheet', 'Tukey-Kramer-Gain');
    writetable(tukey_bins_tbl, anova_tbl_xlsx_file, 'Sheet', 'Tukey-Kramer-Bins');

    disp(['ANOVA and Post-Hoc Tukey-Kramer results saved to: ' anova_tbl_xlsx_file]);
end
