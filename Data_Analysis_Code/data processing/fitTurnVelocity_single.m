% fitTurnVelocity_single
% This function fits linear models to the relationship between object position bins
% and turning velocity across different conditions for a single genotype.
% It calculates slope values, confidence intervals, and correlation coefficients
% for each animal across specified conditions.
%
% INPUTS:
%   posBins     - Array containing the bin centers for object positions (degrees)
%   allEVT      - 3D array of angular velocities (bins x conditions x animals)
%   rangeValue  - Range of interest for analysis (e.g., -20 to 20 degrees)
%
% OUTPUTS:
%   allSlopes   - Slope values across conditions (cell array, one cell per condition)
%   allCI       - Confidence intervals for slopes (cell array, one cell per condition)
%   allR        - R-squared values across conditions (cell array, one cell per condition)
%
% CREATED: 09/02/2025 - MC
%
function [allSlopes, allCI, allR] = fitTurnVelocity_single(posBins, allEVT, rangeValue)

% Get the number of gain conditions and animals
[~, numConditions, nAnimals] = size(allEVT);

% Initialize outputs
allSlopes = cell(numConditions, 1);
allCI     = cell(numConditions, 1);
allR      = cell(numConditions, 1);

% Range of interest (e.g., -20 to 20 deg)
rangeMask = (posBins >= -rangeValue) & (posBins <= rangeValue);

% Loop over each condition
for condition = 1:numConditions
    slopesCond = [];
    RsCond     = [];

    % Fit linear models for each animal (fit intercept)
    for animal = 1:nAnimals
        xData = posBins(rangeMask);
        yData = allEVT(rangeMask, condition, animal);
        
        % Check if there are enough finite data points
        if sum(isfinite(yData)) < 2
            slopesCond(end+1) = NaN;
            RsCond(end+1)     = NaN;
        else
            mdl = fitlm(xData, yData);  % slope + intercept
            slopesCond(end+1) = mdl.Coefficients.Estimate(2);   % slope
            RsCond(end+1)     = mdl.Rsquared.Ordinary;          % R²
        end

    end

    % Store results
    allSlopes{condition} = slopesCond;
    allCI{condition}     = calculateCI(slopesCond, []);
    allR{condition}      = RsCond;
end
end
