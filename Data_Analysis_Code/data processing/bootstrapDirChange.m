function [pValueRaw, pValueFDR, kirWtDiffs, kirNaDiffs, naWtDiffs] = bootstrapDirChange(kirData, wtData, naData, nBootstraps)
    % Initialize outputs
    nBins = size(kirData, 1); % Number of bins (rows)
    kirWtDiffs = cell(nBins, 1);
    kirNaDiffs = cell(nBins, 1);
    naWtDiffs = cell(nBins, 1);
    pValueRaw = zeros(nBins, 3); % Store raw p-values for each bin and comparison

    % Function to compute bootstrap differences between two genotypes
    function bootDiffs = computeBootstrapDifferences(group1, group2, nBootstraps)
        n1 = size(group1, 2); % Number of animals in group 1
        n2 = size(group2, 2); % Number of animals in group 2
        bootDiffs = zeros(nBootstraps, 1);
        for i = 1:nBootstraps
            % Sample with replacement
            sample1 = group1(:, randi(n1, n1, 1));
            sample2 = group2(:, randi(n2, n2, 1));
            % Compute the mean difference
            bootDiffs(i) = mean(sample2(:)) - mean(sample1(:));
        end
    end

    % Compute bootstrap differences and raw p-values for each bin
    for bin = 1:nBins
        kirWtDiffs{bin} = computeBootstrapDifferences(kirData(bin, :), wtData(bin, :), nBootstraps);
        kirNaDiffs{bin} = computeBootstrapDifferences(kirData(bin, :), naData(bin, :), nBootstraps);
        naWtDiffs{bin} = computeBootstrapDifferences(naData(bin, :), wtData(bin, :), nBootstraps);

        % Compute raw p-values for each comparison
        pValueRaw(bin, 1) = computePValue(kirWtDiffs{bin});
        pValueRaw(bin, 2) = computePValue(kirNaDiffs{bin});
        pValueRaw(bin, 3) = computePValue(naWtDiffs{bin});
    end

    % FDR correction across bins for each comparison
    pValueFDR = zeros(size(pValueRaw));
    for comparison = 1:3
        [pValuesSorted, sortIdx] = sort(pValueRaw(:, comparison));
        bhThresholds = (1:nBins)' / nBins * 0.05; % FDR at alpha = 0.05
        isSignificant = pValuesSorted <= bhThresholds;
        significantIdx = find(isSignificant, 1, 'last');

        % Assign FDR-corrected p-values
        pValueFDR(:, comparison) = 1; % Default to not significant
        if ~isempty(significantIdx)
            pValueFDR(sortIdx(1:significantIdx), comparison) = pValuesSorted(1:significantIdx);
        end
    end

    % Nested function to compute raw p-value from bootstrapped differences
    function pValue = computePValue(bootDiffs)
        % Calculate the proportion of differences greater or less than 0
        pValue = mean(bootDiffs > 0);
        % Adjust p-value for two-tailed test
        pValue = 2 * min(pValue, 1 - pValue);
    end
end
