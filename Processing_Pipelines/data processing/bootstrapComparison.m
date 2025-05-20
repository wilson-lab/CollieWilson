% bootstrapDifferenceComparison
% This function performs bootstrapping to compare the significance of differences 
% in means between KIR, WT, and NA groups at each binned object position. 
% It calculates 95% confidence intervals for each bootstrapped mean and checks 
% whether the means of one group fall within the confidence interval of another 
% group to assess significance.
%
% INPUTS:
%   kirBinnedData   - 2D array where rows correspond to binned object positions 
%                     and columns correspond to measured output for each animal in the KIR group.
%   wtBinnedData    - 2D array where rows correspond to binned object positions 
%                     and columns correspond to measured output for each animal in the WT group.
%   naBinnedData    - 2D array where rows correspond to binned object positions 
%                     and columns correspond to measured output for each animal in the NA group.
%   n_bootstraps    - Number of bootstrap iterations.
%   assumption       - String input to specify the assumption ('kir<control' or 'kir>control').
%
% OUTPUTS:
%   p_vals          - Matrix of p-values for comparisons across bins.
%                     Each row corresponds to a binned object position and contains
%                     [p_val_kir_vs_wt, p_val_kir_vs_na, p_val_na_vs_wt].
%   CI_wt_kir      - Confidence intervals for WT vs KIR differences across bins.
%   CI_na_kir      - Confidence intervals for NA vs KIR differences across bins.
%   median_wt_kir  - Median bootstrap difference for WT vs KIR.
%   median_na_kir  - Median bootstrap difference for NA vs KIR.
%
% CREATED: [Date] MC
%
function [p_vals, kir_CI, wt_CI, na_CI] = bootstrapComparison(kirBinnedData, wtBinnedData, naBinnedData, n_bootstraps)
    %% initialize
    % Find the number of bins
    n_bins = size(kirBinnedData, 1);
    % Find the number of animals per group
    nKIR = size(kirBinnedData,2);
    nWT = size(wtBinnedData,2);
    nNA = size(naBinnedData,2);
    binReq = 2; % if < n/binReq bins are available, do not calculate significance

    % Pre-allocate for p-values and confidence intervals
    p_vals = nan(n_bins, 3);  % [kir vs wt, kir vs na, na vs wt]
    kir_CI = nan(n_bins, 2);  % Confidence intervals for kir group [lower_CI upper_CI]
    wt_CI = nan(n_bins, 2);   % Confidence intervals for wt group [lower_CI upper_CI]
    na_CI = nan(n_bins, 2);   % Confidence intervals for na group [lower_CI upper_CI]

    %% Loop through each bin
    for bin = 1:n_bins
        % Remove NaNs from each group before bootstrapping
        kir_data_clean = kirBinnedData(bin, :);
        kir_data_clean = kir_data_clean(~isnan(kir_data_clean));  % Remove NaNs

        wt_data_clean = wtBinnedData(bin, :);
        wt_data_clean = wt_data_clean(~isnan(wt_data_clean));  % Remove NaNs

        na_data_clean = naBinnedData(bin, :);
        na_data_clean = na_data_clean(~isnan(na_data_clean));  % Remove NaNs

        % Bootstrapping
        kir_boot_means = bootstrapMeans(kir_data_clean, n_bootstraps);
        wt_boot_means = bootstrapMeans(wt_data_clean, n_bootstraps);
        na_boot_means = bootstrapMeans(na_data_clean, n_bootstraps);

        % Compute the 95% confidence intervals (percentiles)
        kir_CI(bin, :) = prctile(kir_boot_means, [2.5 97.5]);
        wt_CI(bin, :) = prctile(wt_boot_means, [2.5 97.5]);
        na_CI(bin, :) = prctile(na_boot_means, [2.5 97.5]);

        % Mean for each group
        kir_mean = mean(kir_boot_means);
        wt_mean = mean(wt_boot_means);
        na_mean = mean(na_boot_means);

        % Skip bins where there are fewer than 5 data points in any group
        if numel(kir_data_clean) < nKIR/binReq || numel(wt_data_clean) < nWT/binReq || numel(na_data_clean) < nNA/binReq
            continue;
        end

        % Check if means fall within CIs of other groups
        % kir vs wt
        if kir_mean >= wt_CI(bin, 1) && kir_mean <= wt_CI(bin, 2)
            p_val_kir_vs_wt = 1;  % No significant difference
        else
            p_val_kir_vs_wt = 0;  % Significant difference
        end

        % kir vs na
        if kir_mean >= na_CI(bin, 1) && kir_mean <= na_CI(bin, 2)
            p_val_kir_vs_na = 1;  % No significant difference
        else
            p_val_kir_vs_na = 0;  % Significant difference
        end

        % na vs wt
        if na_mean >= wt_CI(bin, 1) && na_mean <= wt_CI(bin, 2)
            p_val_na_vs_wt = 1;  % No significant difference
        else
            p_val_na_vs_wt = 0;  % Significant difference
        end

        % Store p-values
        p_vals(bin, 1) = p_val_kir_vs_wt;  % kir vs wt
        p_vals(bin, 2) = p_val_kir_vs_na;  % kir vs na
        p_vals(bin, 3) = p_val_na_vs_wt;   % na vs wt
    end

    %% Optional Plot
    optPlot=0;
    if optPlot
        figure;
        for bin = 1:n_bins
            subplot(ceil(n_bins / 3), 3, bin);
            hold on;
            plot([1 2 3], [mean(kir_boot_means), mean(wt_boot_means), mean(na_boot_means)], 'o');
            errorbar(1, mean(kir_boot_means), kir_CI(bin, 1), kir_CI(bin, 2), 'r');
            errorbar(2, mean(wt_boot_means), wt_CI(bin, 1), wt_CI(bin, 2), 'g');
            errorbar(3, mean(na_boot_means), na_CI(bin, 1), na_CI(bin, 2), 'b');
            xticks([1 2 3]);
            xticklabels({'KIR', 'WT', 'NA'});
            title(['Bin ' num2str(bin)]);
            ylabel('Bootstrapped Means');
            hold off;
        end
        sgtitle('Bootstrapped Means with 95% Confidence Intervals');
    end
end