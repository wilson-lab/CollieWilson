% bootstrapDifferenceComparisonCluster
% This function compares slopes between genotypes (KIR vs WT, KIR vs NA) using bootstrapping 
% to assess the significance of differences. It calculates the difference between the bootstrapped 
% means of the two groups based on the specified assumption and calculates p-values based on 
% the position of 0 within the bootstrapped difference distribution.
%
% INPUTS:
%   posBins        - Vector containing the centers of each bin (corresponding to rows in the binned data)
%   kirBinnedData  - 2D array of KIR binned data (bins x conditions)
%   wtBinnedData   - 2D array of WT binned data (bins x conditions)
%   naBinnedData   - 2D array of NA binned data (bins x conditions)
%   n_bootstraps   - Number of bootstrap iterations to perform
%   assumption      - String to specify the hypothesis ('kir<control' or 'kir>control')
%
% OUTPUTS:
%   pValueRaw      - Array of raw p-values for each comparison
%   pValueBonf     - Array of Bonferroni-corrected p-values
%   pValueFDR      - Array of FDR-corrected p-values
%   kirWtDiffs     - Bootstrapped distribution of differences between KIR and WT slopes
%   kirNaDiffs     - Bootstrapped distribution of differences between KIR and NA slopes
%   naWtDiffs      - Bootstrapped distribution of differences between NA and WT slopes
%
% CREATED: [Date] MC
%
function [p_vals, CI_wt_kir, CI_na_kir, median_wt_kir, median_na_kir] = bootstrapDifferenceComparison(kirBinnedData, wtBinnedData, naBinnedData, n_bootstraps, assumption)
    %% Initialize
    % Find the number of bins
    n_bins = size(kirBinnedData, 1);

    % Pre-allocate for p-values, confidence intervals, and medians
    p_vals = nan(n_bins, 2);  % [kir vs wt, kir vs na]
    CI_wt_kir = nan(n_bins, 2);  % Confidence intervals for wt vs kir differences
    CI_na_kir = nan(n_bins, 2);  % Confidence intervals for na vs kir differences
    median_wt_kir = nan(n_bins, 1);  % Median bootstrap difference for wt vs kir
    median_na_kir = nan(n_bins, 1);  % Median bootstrap difference for na vs kir

    %% Loop through each bin
    for bin = 1:n_bins
        % Remove NaNs from each group before bootstrapping
        kir_data_clean = kirBinnedData(bin, :);
        kir_data_clean = kir_data_clean(~isnan(kir_data_clean));  % Remove NaNs

        wt_data_clean = wtBinnedData(bin, :);
        wt_data_clean = wt_data_clean(~isnan(wt_data_clean));  % Remove NaNs

        na_data_clean = naBinnedData(bin, :);
        na_data_clean = na_data_clean(~isnan(na_data_clean));  % Remove NaNs

        % Skip bins where there are fewer than 5 data points in any group
        if numel(kir_data_clean) < 5 || numel(wt_data_clean) < 5 || numel(na_data_clean) < 5
            continue;
        end

        % Bootstrapping the difference for control - kir based on assumption
        if strcmp(assumption, 'kir<control')
            boot_diff_wt_kir = bootstrapDifference(wt_data_clean, kir_data_clean, n_bootstraps);  % wt - kir
            boot_diff_na_kir = bootstrapDifference(na_data_clean, kir_data_clean, n_bootstraps);  % na - kir
        elseif strcmp(assumption, 'kir>control')
            boot_diff_wt_kir = bootstrapDifference(kir_data_clean, wt_data_clean, n_bootstraps);  % kir - wt
            boot_diff_na_kir = bootstrapDifference(kir_data_clean, na_data_clean, n_bootstraps);  % kir - na
        else
            error('Invalid assumption. Use "kir<control" or "kir>control".');
        end

        % Calculate confidence intervals (percentiles) for both comparisons
        CI_wt_kir(bin, :) = prctile(boot_diff_wt_kir, [5 95]);  % 95% CI for wt vs kir
        CI_na_kir(bin, :) = prctile(boot_diff_na_kir, [5 95]);  % 95% CI for na vs kir

        % Calculate the median bootstrap difference
        median_wt_kir(bin) = median(boot_diff_wt_kir);
        median_na_kir(bin) = median(boot_diff_na_kir);

        % Calculate p-values for both comparisons
        % Proportion of bootstraps greater or less than 0
        prop_wt_kir_greater = mean(boot_diff_wt_kir > 0);
        prop_wt_kir_less = mean(boot_diff_wt_kir < 0);
        
        prop_na_kir_greater = mean(boot_diff_na_kir > 0);
        prop_na_kir_less = mean(boot_diff_na_kir < 0);

        % Two-tailed p-value: if 0 is near the center of the distribution, p should be large
        p_val_wt_kir = 2 * min(prop_wt_kir_greater, prop_wt_kir_less);
        p_val_na_kir = 2 * min(prop_na_kir_greater, prop_na_kir_less);

        % Handle the case where 0 is at the exact center (would otherwise return 0 p-value)
        if p_val_wt_kir == 0
            p_val_wt_kir = 1;  % Non-significant if 0 is central
        end

        if p_val_na_kir == 0
            p_val_na_kir = 1;  % Non-significant if 0 is central
        end

        % Store p-values
        p_vals(bin, 1) = p_val_wt_kir;  % wt vs kir
        p_vals(bin, 2) = p_val_na_kir;  % na vs kir
    end
end
