% bootstrapDifferenceComparisonCluster
% This function performs bootstrap comparisons of slopes between genotypes (KIR vs WT, KIR vs NA) 
% using bootstrapping to assess the significance of differences while applying cluster-based 
% test statistics within a specified range (default 0-40 degrees). It calculates raw p-values, 
% Bonferroni-corrected p-values, and FDR-corrected p-values.
%
% INPUTS:
%   posBins      - Vector containing the center of each bin (corresponding to rows in the binned data)
%   kirBinnedData - 2D array of KIR binned data (bins x conditions)
%   wtBinnedData  - 2D array of WT binned data (bins x conditions)
%   naBinnedData  - 2D array of NA binned data (bins x conditions)
%   n_bootstraps  - Number of bootstrap iterations to perform
%   assumption     - String to specify the hypothesis ('kir<control' or 'kir>control')
%
% OUTPUTS:
%   pValueRaw     - Array of raw p-values for each comparison
%   pValueBonf    - Array of Bonferroni-corrected p-values
%   pValueFDR     - Array of FDR-corrected p-values
%   kirWtDiffs    - Bootstrapped distribution of differences between KIR and WT slopes
%   kirNaDiffs    - Bootstrapped distribution of differences between KIR and NA slopes
%   naWtDiffs     - Bootstrapped distribution of differences between NA and WT slopes
%
% CREATED: [Date] MC
%
function [p_vals_corrected, CI_wt_kir, CI_na_kir, median_wt_kir, median_na_kir] = bootstrapDifferenceComparisonCluster(posBins, kirBinnedData, wtBinnedData, naBinnedData, n_bootstraps, assumption)
    % Define the bin range to test and automatically set the cluster size to all bins in the range
    binRange = [0 40];  % Test bins between 0-40 by default
    clusterBins = posBins >= binRange(1) & posBins <= binRange(2);  % Logical vector to select bins in the range

    % Get indices of the bins within the specified range
    binIndices = find(clusterBins);
    if isempty(binIndices)
        error('No bins within the specified binRange.');
    end

    % Pre-allocate for p-values, confidence intervals, and medians (only one cluster)
    p_vals = nan(1, 2);  % [kir vs wt, kir vs na]
    CI_wt_kir = nan(1, 2);  % Confidence intervals for wt vs kir (1 cluster)
    CI_na_kir = nan(1, 2);  % Confidence intervals for na vs kir (1 cluster)
    median_wt_kir = nan(1, 1);  % Median bootstrap difference for wt vs kir (1 cluster)
    median_na_kir = nan(1, 1);  % Median bootstrap difference for na vs kir (1 cluster)

    % Aggregate the data for the entire cluster (all bins within the binRange)
    kir_cluster = kirBinnedData(binIndices, :);
    wt_cluster = wtBinnedData(binIndices, :);
    na_cluster = naBinnedData(binIndices, :);

    % Flatten the data across all bins in the cluster
    kir_data_clean = kir_cluster(:);
    wt_data_clean = wt_cluster(:);
    na_data_clean = na_cluster(:);

    % Remove NaNs from each group before bootstrapping
    kir_data_clean = kir_data_clean(~isnan(kir_data_clean));
    wt_data_clean = wt_data_clean(~isnan(wt_data_clean));
    na_data_clean = na_data_clean(~isnan(na_data_clean));

    % Skip if there are fewer than 5 data points in any group
    if numel(kir_data_clean) < 5 || numel(wt_data_clean) < 5 || numel(na_data_clean) < 5
        error('Not enough data points in one or more groups.');
    end

    % Calculate the observed difference (cluster-level statistic)
    if strcmp(assumption, 'kir<control')
        observed_diff_wt_kir = mean(wt_data_clean) - mean(kir_data_clean);  % wt - kir
        observed_diff_na_kir = mean(na_data_clean) - mean(kir_data_clean);  % na - kir
    elseif strcmp(assumption, 'kir>control')
        observed_diff_wt_kir = mean(kir_data_clean) - mean(wt_data_clean);  % kir - wt
        observed_diff_na_kir = mean(kir_data_clean) - mean(na_data_clean);  % kir - na
    else
        error('Invalid assumption. Use "kir<control" or "kir>control".');
    end

    % Bootstrap the difference for control - kir
    boot_diff_wt_kir = bootstrapDifference(wt_data_clean, kir_data_clean, n_bootstraps);
    boot_diff_na_kir = bootstrapDifference(na_data_clean, kir_data_clean, n_bootstraps);

    % Calculate confidence intervals (percentiles) for both comparisons
    CI_wt_kir(1, :) = prctile(boot_diff_wt_kir, [5 95]);  % 95% CI for wt vs kir
    CI_na_kir(1, :) = prctile(boot_diff_na_kir, [5 95]);  % 95% CI for na vs kir

    % Calculate the median bootstrap difference
    median_wt_kir(1) = median(boot_diff_wt_kir);
    median_na_kir(1) = median(boot_diff_na_kir);

    % Calculate continuous p-values based on where 0 falls in the bootstrap distribution
    % wt vs kir
    p_val_wt_kir = mean(boot_diff_wt_kir >= 0);  % Proportion of bootstrap differences >= 0

    % na vs kir
    p_val_na_kir = mean(boot_diff_na_kir >= 0);  % Proportion of bootstrap differences >= 0

    % Combine p-values for FDR correction
    p_vals(1, 1) = p_val_wt_kir;
    p_vals(1, 2) = p_val_na_kir;

    % Apply FDR correction to control for multiple comparisons
    p_vals_corrected = applyFDR(p_vals, 0.05);  % Apply FDR correction to p-values with a significance level of 0.05
end
