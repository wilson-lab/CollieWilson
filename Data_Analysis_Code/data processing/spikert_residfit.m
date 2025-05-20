% SPIKERT_RESIDFIT - Perform residual-based partitioning of variance explained for firing rate
%
% This function calculates the unique variance explained in firing rate (spikert) 
% by forward velocity, angular velocity (only for ipsilateral turns), and 
% a Gaussian fit based on visual object position, using a residual-based approach.
% Two fitting orders are tested:
%   - Order A: Fit forward and angular velocity, then fit the Gaussian model on the residuals.
%   - Order B: Fit the Gaussian model first, then fit forward and angular velocity on the residuals.
%
% INPUTS:
%   spikert        - Array of spike rates (firing rate) (time x trials)
%   panelps        - Array of visual object positions (time x trials)
%   forward        - Array of forward velocities (time x trials)
%   angular        - Array of angular velocities (time x trials)
%   ttime          - Time vector (in seconds) for lag shifting the velocities
%   peak_position  - The fixed peak location for the Gaussian model (in degrees)
%   spread         - The fixed spread (standard deviation) for the Gaussian model (in degrees)
%
% OUTPUTS:
%   R2             - Array of R² values for each model fit in the following order:
%                      R2(1) - R² for forward + angular velocity only (Order A)
%                      R2(2) - R² for Gaussian model on residuals (Order A)
%                      R2(3) - R² for Gaussian model only (Order B)
%                      R2(4) - R² for forward + angular velocity on residuals (Order B)
%   mdlfit         - Structure containing final model fits on the full dataset:
%                      mdlfit.forward_angular           - Full model for forward + angular velocity
%                      mdlfit.gaussian_on_residuals     - Gaussian fit on residuals from forward + angular velocity
%                      mdlfit.gaussian_only             - Gaussian model alone
%                      mdlfit.forward_angular_on_residuals - Forward + angular velocity fit on residuals from Gaussian model
%
% CREATED: 11/07/2024 - MC
%
function [R2, mdlfit] = spikert_residfit(spikert, panelps, forward, angular, ttime, peak_position, spread)

    % Fetch processing settings
    settings = processSettings;

    % Adjust forward and angular velocity for lag
    forward_lag = lagShift(forward, ttime, settings.fwdLag);
    angular_lag = lagShift(angular, ttime, settings.angLag);

    % Replace NaNs in visual object position with a placeholder value
    placeholder_value = -999;
    panelps(isnan(panelps)) = placeholder_value;

    % Reshape variables to ensure they are vectors
    firing_rate = reshape(spikert, [], 1);
    forward_velocity = reshape(forward_lag, [], 1);
    angular_velocity = reshape(angular_lag, [], 1);
    visual_object_position = reshape(panelps, [], 1);

    % Define valid indices for data (ipsilateral turns only, no NaNs)
    valid_idx = ~any(isnan([forward_velocity, angular_velocity]), 2) & visual_object_position ~= placeholder_value & angular_velocity > 0;

    % Separate valid data for each variable
    X_forward = forward_velocity(valid_idx);
    X_angular = angular_velocity(valid_idx);
    X_gaussian = visual_object_position(valid_idx);
    y = firing_rate(valid_idx);

    % Calculate baseline as the mean firing rate
    baseline_estimate = mean(spikert(:), 'omitnan');

    % Define Gaussian model with fixed baseline, peak, and spread
    gaussian_model = @(b, X) baseline_estimate + b * exp(-((X - peak_position).^2) / (2 * spread^2)) .* (X >= -15 & X <= 180);

    %% Order A: Fit forward + angular first, then Gaussian on residuals

    % Step 1: Fit forward + angular velocity model
    mdl_forward_angular = fitlm([X_forward, X_angular], y);
    residuals_forward_angular = mdl_forward_angular.Residuals.Raw; % Residuals after forward + angular fit

    % Step 2: Fit Gaussian model to residuals from forward + angular model
    mdl_gaussian_on_residuals = fitnlm(X_gaussian, residuals_forward_angular, gaussian_model, 1);

    % Store R² values for Order A
    R2_orderA_forward_angular = mdl_forward_angular.Rsquared.Ordinary;
    R2_orderA_gaussian = mdl_gaussian_on_residuals.Rsquared.Ordinary;

    %% Order B: Fit Gaussian model first, then forward + angular on residuals

    % Step 1: Fit Gaussian model alone
    mdl_gaussian_only = fitnlm(X_gaussian, y, gaussian_model, 1);
    residuals_gaussian = mdl_gaussian_only.Residuals.Raw; % Residuals after Gaussian model fit

    % Step 2: Fit forward + angular model to residuals from Gaussian model
    mdl_forward_angular_on_residuals = fitlm([X_forward, X_angular], residuals_gaussian);

    % Store R² values for Order B
    R2_orderB_gaussian = mdl_gaussian_only.Rsquared.Ordinary;
    R2_orderB_forward_angular = mdl_forward_angular_on_residuals.Rsquared.Ordinary;

    %% Store R² values in an array and model fits in a structure

    % Final R² output array
    R2 = [R2_orderA_forward_angular, R2_orderA_gaussian, R2_orderB_gaussian, R2_orderB_forward_angular];

    % Final model fits in a structure
    mdlfit = struct();
    mdlfit.forward_angular = mdl_forward_angular;
    mdlfit.gaussian_on_residuals = mdl_gaussian_on_residuals;
    mdlfit.gaussian_only = mdl_gaussian_only;
    mdlfit.forward_angular_on_residuals = mdl_forward_angular_on_residuals;
end
