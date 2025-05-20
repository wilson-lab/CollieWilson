% calculateCI
% This helper function calculates the mean and 95% confidence intervals for slope 
% and intercept estimates from a linear model. It uses the t-distribution to determine 
% the confidence intervals based on the standard error of the mean.
%
% INPUTS:
%   slopes      - Array of slope estimates (numeric vector)
%   intercepts  - Array of intercept estimates (numeric vector)
%
% OUTPUTS:
%   fitCI       - Structure containing the following fields:
%                 - slope: 95% confidence interval for the slope
%                 - intercept: 95% confidence interval for the intercept
%                 - meanSlope: Mean value of the slopes
%                 - meanIntercept: Mean value of the intercepts
%
% CREATED: [Date] MC
%
function fitCI = calculateCI(slopes, intercepts)
meanSlope = mean(slopes);
meanIntercept = mean(intercepts);

% Compute 95% CI using t-distribution
SEM_slope = std(slopes) / sqrt(length(slopes));
SEM_intercept = std(intercepts) / sqrt(length(intercepts));

tVal = tinv(0.975, length(slopes) - 1);  % 95% CI

slopeCI = [meanSlope - tVal * SEM_slope, meanSlope + tVal * SEM_slope];
interceptCI = [meanIntercept - tVal * SEM_intercept, meanIntercept + tVal * SEM_intercept];

fitCI.slope = slopeCI;
fitCI.intercept = interceptCI;
fitCI.meanSlope = meanSlope;
fitCI.meanIntercept = meanIntercept;
end