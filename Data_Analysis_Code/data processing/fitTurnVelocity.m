% fitTurnVelocity
% This function fits linear models to the relationship between object position bins 
% and turning velocity across different conditions for each genotype (KIR, WT, NA). 
% It calculates slope values, confidence intervals, and correlation coefficients for 
% each animal across specified conditions.
%
% INPUTS:
%   posBins    - Array containing the bin centers for object positions (degrees)
%   kirEVT     - 3D array of KIR angular velocities (bins x conditions x animals)
%   wtEVT      - 3D array of WT angular velocities (bins x conditions x animals)
%   naEVT      - 3D array of NA angular velocities (bins x conditions x animals)
%   rangeValue  - Range of interest for analysis (e.g., -20 to 20 degrees)
%
% OUTPUTS:
%   kirSlopes  - Slope values for KIR genotype across conditions
%   wtSlopes   - Slope values for WT genotype across conditions
%   naSlopes   - Slope values for NA genotype across conditions
%   kirCI      - Confidence intervals for KIR slopes across animals
%   wtCI       - Confidence intervals for WT slopes across animals
%   naCI       - Confidence intervals for NA slopes across animals
%   kirR       - Correlation coefficients for KIR genotype across conditions
%   wtR        - Correlation coefficients for WT genotype across conditions
%   naR        - Correlation coefficients for NA genotype across conditions
%
% CREATED: [Date] MC
%
function [kirSlopes, wtSlopes, naSlopes, kirCI, wtCI, naCI, kirR, wtR, naR] = fitTurnVelocity(posBins, kirEVT, wtEVT, naEVT, rangeValue)
% Get the number of gain conditions and animals for each genotype
[~, numConditions, nKIR] = size(kirEVT);
[~, ~, nWT] = size(wtEVT);
[~, ~, nNA] = size(naEVT);

% Initialize output structures
kirSlopes = cell(numConditions, 1);
wtSlopes = cell(numConditions, 1);
naSlopes = cell(numConditions, 1);

kirCI = cell(numConditions, 1);  % Confidence intervals for KIR group
wtCI = cell(numConditions, 1);   % Confidence intervals for WT group
naCI = cell(numConditions, 1);   % Confidence intervals for NA group

kirR = cell(numConditions, 1);  % R-squared values for KIR group
wtR = cell(numConditions, 1);   % R-squared values for WT group
naR = cell(numConditions, 1);   % R-squared values for NA group

% Range of interest (e.g., -20 to 20)
rangeMask = (posBins >= -rangeValue) & (posBins <= rangeValue);

% Loop over each condition
for condition = 1:numConditions
    kirSlopesCond = [];
    wtSlopesCond = [];
    naSlopesCond = [];

    kirRsCond = [];  % R-squared for KIR group
    wtRsCond = [];   % R-squared for WT group
    naRsCond = [];   % R-squared for NA group

    % Fit linear models for each animal in kir group (fit intercept)
    for animal = 1:nKIR
        xData = posBins(rangeMask);
        yData_kir = kirEVT(rangeMask, condition, animal);
        mdl_kir = fitlm(xData, yData_kir);  % Fit both slope and intercept
        kirSlopesCond(end+1) = mdl_kir.Coefficients.Estimate(2);  % Slope
        kirRsCond(end+1) = mdl_kir.Rsquared.Ordinary;  % R-squared (r value)
    end

    % Fit linear models for each animal in wt group (fit intercept)
    for animal = 1:nWT
        xData = posBins(rangeMask);
        yData_wt = wtEVT(rangeMask, condition, animal);
        mdl_wt = fitlm(xData, yData_wt);  % Fit both slope and intercept
        wtSlopesCond(end+1) = mdl_wt.Coefficients.Estimate(2);  % Slope
        wtRsCond(end+1) = mdl_wt.Rsquared.Ordinary;  % R-squared (r value)
    end

    % Fit linear models for each animal in na group (fit intercept)
    for animal = 1:nNA
        xData = posBins(rangeMask);
        yData_na = naEVT(rangeMask, condition, animal);
        mdl_na = fitlm(xData, yData_na);  % Fit both slope and intercept
        naSlopesCond(end+1) = mdl_na.Coefficients.Estimate(2);  % Slope
        naRsCond(end+1) = mdl_na.Rsquared.Ordinary;  % R-squared (r value)
    end

    % Store slopes, confidence intervals, and R-squared values for each condition
    kirSlopes{condition} = kirSlopesCond;
    wtSlopes{condition} = wtSlopesCond;
    naSlopes{condition} = naSlopesCond;

    kirCI{condition} = calculateCI(kirSlopesCond, []);
    wtCI{condition} = calculateCI(wtSlopesCond, []);
    naCI{condition} = calculateCI(naSlopesCond, []);

    kirR{condition} = kirRsCond;
    wtR{condition} = wtRsCond;
    naR{condition} = naRsCond;
end
end
