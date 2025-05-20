% bootstrap_difference
% This function performs bootstrapping to calculate the difference in medians 
% between two genotypes. It returns the p-value and the bootstrap difference 
% distribution, allowing for statistical comparison of the two groups.
%
% INPUTS:
%   genotype1  - Data for the first genotype (vector)
%   genotype2  - Data for the second genotype (vector)
%   nBootstrap  - Number of bootstrap iterations (scalar)
%
% OUTPUTS:
%   p_value    - P-value based on the percentile of 0 in the bootstrapped 
%                difference distribution (right-tailed test)
%   boot_diff   - Distribution of bootstrapped differences between the medians 
%                 of the two genotypes
%
% CREATED: [Date] MC
%
function [p_value, boot_diff] = bootstrap_difference(genotype1, genotype2, nBootstrap)
    % Get the number of samples in each genotype
    nGenotype1 = length(genotype1);
    nGenotype2 = length(genotype2);
    
    % Initialize array to store bootstrap differences
    boot_diff = zeros(nBootstrap, 1);
    
    % Perform bootstrap sampling
    for i = 1:nBootstrap
        % Resample with replacement from each genotype
        genotype1_boot = genotype1(randi(nGenotype1, nGenotype1, 1));
        genotype2_boot = genotype2(randi(nGenotype2, nGenotype2, 1));
        
        % Calculate the difference in medians for this bootstrap sample
        boot_diff(i) = median(genotype1_boot) - median(genotype2_boot);
    end
    
    % Calculate the p-value based on the percentile of 0 in the difference distribution
    p_value = mean(boot_diff >= 0); % Right-tailed test
    
end
