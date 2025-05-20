% applyFDR
%
% Applies the False Discovery Rate (FDR) correction to a set of p-values
% using the Benjamini-Hochberg procedure to control for the rate of
% false discoveries in multiple hypothesis testing.
%
% INPUTS:
% - p_vals : Array of p-values to be corrected (e.g., results from
%             statistical tests).
% - q_level: Desired FDR level (e.g., 0.05), representing the maximum
%             allowable proportion of false discoveries.
%
% OUTPUT:
% - p_vals_corrected : Array of corrected p-values (FDR-adjusted), where
%                      significant p-values are adjusted to account for
%                      multiple comparisons.
%
% PROCESS:
% The function first determines the number of comparisons and sorts the
% p-values while keeping track of their original indices. It then applies
% the Benjamini-Hochberg procedure by calculating the corresponding
% thresholds. The largest p-value that is below the threshold is found,
% and all p-values up to that index are marked as significant in the
% output array.
%
% The result is an array of corrected p-values where significant
% discoveries are preserved based on the specified FDR level.
%
function p_vals_corrected = applyFDR(p_vals, q_level)
    % Number of comparisons
    m = numel(p_vals);
    
    % Sort p-values and store the original indices
    [sorted_p_vals, sort_idx] = sort(p_vals);
    
    % Apply the Benjamini-Hochberg procedure
    bh_thresholds = (1:m) / m * q_level;
    
    % Find the largest p-value that is smaller than the threshold
    significant_idx = find(sorted_p_vals <= bh_thresholds, 1, 'last');
    
    % Mark p-values as significant
    p_vals_corrected = ones(size(p_vals));
    if ~isempty(significant_idx)
        p_vals_corrected(sort_idx(1:significant_idx)) = sorted_p_vals(1:significant_idx);
    end
end