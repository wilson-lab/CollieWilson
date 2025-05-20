% run_anova_on_jumps
% This function performs a two-way ANOVA on zero-time points above a threshold for each genotype and gain condition,
% and a one-way ANOVA comparing genotypes across all gain conditions combined.
% It saves the ANOVA and post-hoc Tukey-Kramer results to an Excel file.
%
% INPUTS:
%   kirZeroTAboveThresh  - Zero-time points above the threshold for KIR genotype (cell array by gain condition)
%   wtZeroTAboveThresh   - Zero-time points above the threshold for WT genotype (cell array by gain condition)
%   naZeroTAboveThresh   - Zero-time points above the threshold for NA genotype (cell array by gain condition)
%   dependent_var_name   - Name of the dependent variable (for saving results)
%   folder               - Folder path to save the results
%
% OUTPUTS:
%   pSeparate   - p-values for genotype, gain, and interaction from the two-way ANOVA (structure)
%   pCombined   - p-value for the genotype effect across all gain conditions combined (one-way ANOVA)
%
% CREATED: [Date] MC
%
function [pSeparate, pCombined] = run_anova_on_jumps(kirZeroTAboveThresh, wtZeroTAboveThresh, naZeroTAboveThresh, dependent_var_name, folder)
% Perform two-way ANOVA on zero-time points > threshold for each genotype and gain condition,
nGain = size(naZeroTAboveThresh, 2);
% Preallocate arrays to store data for two-way ANOVA
allZeroT = [];
genotypes = [];
gains = [];

% Combine across all gain conditions for one-way ANOVA
combinedZeroT = [];
combinedGenotypes = [];

% Loop over each genotype
for g = 1:3
    % Fetch data based on genotype
    switch g
        case 1
            zeroTAboveThresh = kirZeroTAboveThresh;
            genotypeLabel = 'KIR';
        case 2
            zeroTAboveThresh = wtZeroTAboveThresh;
            genotypeLabel = 'WT';
        case 3
            zeroTAboveThresh = naZeroTAboveThresh;
            genotypeLabel = 'NA';
    end

    % Loop over each gain condition
    for c = 1:nGain
        % Get the zeroT points > threshold for the current gain condition
        selectedZeroT = zeroTAboveThresh{c};

        % Prepare data for two-way ANOVA
        allZeroT = [allZeroT, selectedZeroT];
        genotypes = [genotypes, repmat({genotypeLabel}, 1, length(selectedZeroT))];
        gains = [gains, repmat(c, 1, length(selectedZeroT))];

        % Combine data across gain conditions for one-way ANOVA
        combinedZeroT = [combinedZeroT, selectedZeroT];
        combinedGenotypes = [combinedGenotypes, repmat({genotypeLabel}, 1, length(selectedZeroT))];
    end
end

% Convert to categorical variables for ANOVA
genotypes = categorical(genotypes);
gains = categorical(gains);
combinedGenotypes = categorical(combinedGenotypes);

%% Two-way ANOVA: Genotype and Gain
% Perform two-way ANOVA
[p, tbl, stats] = anovan(allZeroT, {genotypes, gains}, 'model', 'interaction', 'varnames', {'Genotype', 'Gain'}, 'display', 'off');

% Extract p-values for genotype, gain, and interaction
pSeparate.Genotype = p(1);
pSeparate.Gain = p(2);
pSeparate.Interaction = p(3);

% Display ANOVA results for the two-way ANOVA
disp([dependent_var_name ' Two-Way ANOVA Results:']);
disp(tbl);

% Save ANOVA results to Excel
anova_tbl_xlsx_file = fullfile(folder.summary, ['anova_' dependent_var_name '.xlsx']);
writetable(cell2table(tbl), anova_tbl_xlsx_file, 'Sheet', 'Two-Way ANOVA');

%% One-way ANOVA: Genotype only (across combined gain conditions)
% Perform one-way ANOVA on combined data across all gain conditions
[pCombined, tblCombined, statsCombined] = anovan(combinedZeroT, {combinedGenotypes}, 'varnames', {'Genotype'}, 'display', 'off');

% Display ANOVA results for the one-way ANOVA
disp([dependent_var_name ' One-Way ANOVA (Combined Gains) Results:']);
disp(tblCombined);

% Save ANOVA results for one-way ANOVA to Excel
writetable(cell2table(tblCombined), anova_tbl_xlsx_file, 'Sheet', 'One-Way ANOVA');

% Save Tukey-Kramer post-hoc results (two-way)
% Post-hoc for Genotype
[c_genotype, ~, ~, ~] = multcompare(stats, 'Dimension', 1, 'ctype', 'tukey-kramer', 'display', 'off');

% Post-hoc for Gain
[c_gain, ~, ~, ~] = multcompare(stats, 'Dimension', 2, 'ctype', 'tukey-kramer', 'display', 'off');

% Prepare Tukey-Kramer results for saving
tukey_genotype_tbl = array2table(c_genotype, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'PValue'});
tukey_gain_tbl = array2table(c_gain, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'PValue'});

% Save Tukey-Kramer results
writetable(tukey_genotype_tbl, anova_tbl_xlsx_file, 'Sheet', 'Tukey-Kramer-Genotype');
writetable(tukey_gain_tbl, anova_tbl_xlsx_file, 'Sheet', 'Tukey-Kramer-Gain');

% Display save location
disp(['ANOVA and Post-Hoc Tukey-Kramer results saved to: ' anova_tbl_xlsx_file]);

end
