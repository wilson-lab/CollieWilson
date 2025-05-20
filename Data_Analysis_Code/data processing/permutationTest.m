% permutationTest
% This function performs a permutation test to compare the means of two groups.
% It calculates the observed difference in means and generates a distribution of differences 
% from permuted data to evaluate the significance of the observed difference.
%
% INPUTS:
%   group1        - Numeric vector of data for the first group
%   group2        - Numeric vector of data for the second group
%   n_permutations - Number of permutations to perform for the test
%
% OUTPUT:
%   p_val         - P-value indicating the significance of the difference between the two groups' means
%
% CREATED: [Date] MC
%
function p_val = permutationTest(group1, group2, n_permutations)
    % Permutation test to compare the means of two groups
    observed_diff = abs(mean(group1) - mean(group2));
    combined_data = [group1; group2];
    n_group1 = numel(group1);
    perm_diffs = nan(n_permutations, 1);

    for i = 1:n_permutations
        perm_labels = randperm(numel(combined_data));
        perm_group1 = combined_data(perm_labels(1:n_group1));
        perm_group2 = combined_data(perm_labels(n_group1+1:end));
        perm_diffs(i) = abs(mean(perm_group1) - mean(perm_group2));
    end

    p_val = mean(perm_diffs >= observed_diff);
end