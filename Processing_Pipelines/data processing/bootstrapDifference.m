% bootstrapDifference
% This helper function performs bootstrap resampling to calculate the difference between 
% two groups. The difference is calculated as the absolute mean of group1 minus the absolute 
% mean of group2. It allows for assessing the significance of differences based on the 
% specified assumption (e.g., control vs. kir).
%
% INPUTS:
%   group1        - Data from the first group (control or kir based on assumption)
%   group2        - Data from the second group (kir or control based on assumption)
%   n_bootstraps  - Number of bootstrap iterations to perform
%
% OUTPUT:
%   boot_diff     - Vector of bootstrapped differences between group1 and group2
%
% CREATED: [Date] MC
%
function boot_diff = bootstrapDifference(group1, group2, n_bootstraps)
    n1 = numel(group1);
    n2 = numel(group2);
    boot_diff = zeros(n_bootstraps, 1);

    % Bootstrap resampling and calculating differences using absolute values
    for i = 1:n_bootstraps
        resample_group1 = randsample(group1, n1, true);  % Resample with replacement
        resample_group2 = randsample(group2, n2, true);  % Resample with replacement
        boot_diff(i) = abs(mean(resample_group1)) - abs(mean(resample_group2));  % Calculate absolute difference
    end
end