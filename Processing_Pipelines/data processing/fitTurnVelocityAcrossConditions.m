% fitTurnVelocityAcrossConditions
% This function fits a linear model to the relationship between object position bins and 
% turning velocity across different conditions for each genotype (KIR, WT, NA). It calculates 
% slope values, confidence intervals, and correlation coefficients for each animal.
%
% INPUTS:
%   posBins      - Array containing the bin centers for object positions (degrees)
%   kirEVT      - 3D array of KIR angular velocities (bins x conditions x animals)
%   wtEVT       - 3D array of WT angular velocities (bins x conditions x animals)
%   naEVT       - 3D array of NA angular velocities (bins x conditions x animals)
%   rangeValue   - Range of interest for analysis (e.g., -20 to 20 degrees)
%
% OUTPUTS:
%   kirFit       - Slope values for KIR genotype across conditions
%   wtFit        - Slope values for WT genotype across conditions
%   naFit        - Slope values for NA genotype across conditions
%   kirCI        - Confidence intervals for KIR slopes across animals
%   wtCI         - Confidence intervals for WT slopes across animals
%   naCI         - Confidence intervals for NA slopes across animals
%   kirR         - Correlation coefficients for KIR genotype across conditions
%   wtR          - Correlation coefficients for WT genotype across conditions
%   naR          - Correlation coefficients for NA genotype across conditions
%
% CREATED: [Date] MC
%
function [kirFit, wtFit, naFit, kirCI, wtCI, naCI, kirR, wtR, naR] = fitTurnVelocityAcrossConditions(posBins, kirEVT, wtEVT, naEVT, rangeValue)
% Get the number of gain conditions and animals for each genotype
[~, numConditions, nKIR] = size(kirEVT);
[~, ~, nWT] = size(wtEVT);
[~, ~, nNA] = size(naEVT);

% Initialize output structures
kirFit = zeros(nKIR, 1);
wtFit = zeros(nWT, 1);
naFit = zeros(nNA, 1);

kirR = zeros(nKIR, 1);  % Initialize r values for KIR group
wtR = zeros(nWT, 1);    % Initialize r values for WT group
naR = zeros(nNA, 1);    % Initialize r values for NA group

% Range of interest (e.g., -20 to 20)
rangeMask = (posBins >= -rangeValue) & (posBins <= rangeValue);
xData = posBins(rangeMask);

% Fit linear models for each animal in kir group (fixed intercept = 0)
for animal = 1:nKIR
    % Combine data across all conditions
    yData_kir = reshape(kirEVT(rangeMask, :, animal), [], 1);  % Reshape yData to a column vector
    xData_all = repmat(xData, 1, numConditions);  % Repeat xData for each condition
    xData_all = reshape(xData_all, [], 1);  % Reshape xData to match yData_kir

    % Fit model with fixed intercept at 0
    mdl_kir = fitlm(xData_all, yData_kir, 'Intercept', false);
    kirFit(animal) = mdl_kir.Coefficients.Estimate(1);  % Slope only

    % Store the r value (correlation coefficient)
    kirR(animal) = mdl_kir.Rsquared.Ordinary;  % Use R-squared as r value
end

% Fit linear models for each animal in wt group (fixed intercept = 0)
for animal = 1:nWT
    % Combine data across all conditions
    yData_wt = reshape(wtEVT(rangeMask, :, animal), [], 1);  % Reshape yData to a column vector
    xData_all = repmat(xData, 1, numConditions);  % Repeat xData for each condition
    xData_all = reshape(xData_all, [], 1);  % Reshape xData to match yData_wt

    % Fit model with fixed intercept at 0
    mdl_wt = fitlm(xData_all, yData_wt, 'Intercept', false);
    wtFit(animal) = mdl_wt.Coefficients.Estimate(1);  % Slope only

    % Store the r value (correlation coefficient)
    wtR(animal) = mdl_wt.Rsquared.Ordinary;  % Use R-squared as r value
end

% Fit linear models for each animal in na group (fixed intercept = 0)
for animal = 1:nNA
    % Combine data across all conditions
    yData_na = reshape(naEVT(rangeMask, :, animal), [], 1);  % Reshape yData to a column vector
    xData_all = repmat(xData, 1, numConditions);  % Repeat xData for each condition
    xData_all = reshape(xData_all, [], 1);  % Reshape xData to match yData_na

    % Fit model with fixed intercept at 0
    mdl_na = fitlm(xData_all, yData_na, 'Intercept', false);
    naFit(animal) = mdl_na.Coefficients.Estimate(1);  % Slope only

    % Store the r value (correlation coefficient)
    naR(animal) = mdl_na.Rsquared.Ordinary;  % Use R-squared as r value
end

% Calculate confidence intervals across animals
kirCI = calculateCI(kirFit, []);
wtCI = calculateCI(wtFit, []);
naCI = calculateCI(naFit, []);
end
