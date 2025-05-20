% bootstrapSlopeComparison
% This function performs a bootstrap analysis to compare slopes between different genotypes 
% (KIR, WT, NA) by generating bootstrap samples and calculating p-values for the comparisons. 
% It provides raw p-values, Bonferroni-corrected p-values, and FDR-corrected p-values, along 
% with bootstrapped distributions of differences between the slopes.
%
% INPUTS:
%   kirFit     - Slope values for KIR genotype (vector)
%   wtFit      - Slope values for WT genotype (vector)
%   naFit      - Slope values for NA genotype (vector)
%   nBootstraps - Number of bootstrap samples to generate
%
% OUTPUTS:
%   pValueRaw  - Array of raw p-values for each comparison (KIR vs. WT, KIR vs. NA, NA vs. WT)
%   pValueBonf - Array of Bonferroni-corrected p-values for each comparison
%   pValueFDR  - Array of FDR-corrected p-values for each comparison
%   kirWtDiffs - Bootstrapped distribution of differences between KIR and WT slopes
%   kirNaDiffs - Bootstrapped distribution of differences between KIR and NA slopes
%   naWtDiffs  - Bootstrapped distribution of differences between NA and WT slopes
%
% CREATED: [Date] MC
%
function [pValueRaw, pValueBonf, pValueFDR, kirWtDiffs, kirNaDiffs, naWtDiffs] = bootstrapSlopeComparison(kirFit, wtFit, naFit, nBootstraps)   
    % Initialize bootstrap difference distributions
    kirWtDiffs = zeros(nBootstraps, 1);
    kirNaDiffs = zeros(nBootstraps, 1);
    naWtDiffs = zeros(nBootstraps, 1);
    
    % Function to compute bootstrap differences between two genotypes
    function bootDiffs = computeBootstrapDifferences(group1, group2, nBootstraps)
        n1 = length(group1);
        n2 = length(group2);
        bootDiffs = zeros(nBootstraps, 1);
        for i = 1:nBootstraps
            % Sample with replacement
            sample1 = group1(randi(n1, n1, 1));
            sample2 = group2(randi(n2, n2, 1));
            % Compute the mean difference
            bootDiffs(i) = mean(sample2) - mean(sample1);
        end
    end
    
    % Compute bootstrapped differences between KIR and WT, KIR and NA, and NA and WT
    kirWtDiffs = computeBootstrapDifferences(kirFit, wtFit, nBootstraps);
    kirNaDiffs = computeBootstrapDifferences(kirFit, naFit, nBootstraps);
    naWtDiffs = computeBootstrapDifferences(naFit, wtFit, nBootstraps);
    
    % Compute the raw p-values
    function pValue = computePValue(bootDiffs)
        % Calculate the proportion of differences that are below or above 0
        pValue = mean(bootDiffs > 0);
        % Adjust p-value for two-tailed test
        pValue = 2 * min(pValue, 1 - pValue);
    end
    
    % Calculate raw p-values
    pValueRaw = [computePValue(kirWtDiffs), computePValue(kirNaDiffs), computePValue(naWtDiffs)];
    
    % Bonferroni correction
    nComparisons = length(pValueRaw);  % Number of comparisons (3 in this case)
    pValueBonf = min(pValueRaw * nComparisons, 1);  % Apply Bonferroni correction
    
    % FDR correction using Benjamini-Hochberg
    [pValuesSorted, sortIdx] = sort(pValueRaw);
    bhThresholds = (1:nComparisons) / nComparisons * 0.05;  % FDR at alpha = 0.05
    pValueFDR = zeros(1, nComparisons);
    
    % Determine significant p-values under FDR
    for i = 1:nComparisons
        if pValuesSorted(i) <= bhThresholds(i)
            pValueFDR(sortIdx(i)) = pValuesSorted(i);
        else
            pValueFDR(sortIdx(i)) = 1;  % Not significant under FDR
        end
    end
    
end
