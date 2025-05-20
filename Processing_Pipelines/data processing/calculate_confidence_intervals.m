% calculate_confidence_intervals
% This function calculates the confidence intervals for the overall mean of the input data.
% It computes the mean and standard deviation while accounting for NaN values, and returns
% the lower and upper bounds of the confidence interval based on the specified confidence level.
%
% INPUTS:
%   data              - A matrix of data (rows: observations, columns: groups)
%   confidence_level  - Desired confidence level (e.g., 0.95 for 95% confidence interval)
%
% OUTPUTS:
%   ci_lower          - Lower bound of the confidence interval for the overall mean
%   ci_upper          - Upper bound of the confidence interval for the overall mean
%
% CREATED: [Date] MC
%
function [ci_lower, ci_upper] = calculate_confidence_intervals(data, confidence_level)
    % Calculate overall mean and standard deviation
    overall_mean = mean(data, 'all', 'omitnan');  % Overall mean, ignoring NaNs
    overall_std = std(data, 0, 'all', 'omitnan'); % Overall standard deviation, ignoring NaNs
    n = sum(~isnan(data), 'all');  % Total number of valid observations

    % Calculate the z-score for the desired confidence level
    z = norminv((1 + confidence_level) / 2, 0, 1);  % Z-score

    % Calculate the standard error
    se = overall_std / sqrt(n);  % Standard Error

    % Calculate margin of error
    margin_of_error = z * se;

    % Calculate confidence interval
    ci_lower = overall_mean - margin_of_error;  % Lower CI
    ci_upper = overall_mean + margin_of_error;  % Upper CI
end
