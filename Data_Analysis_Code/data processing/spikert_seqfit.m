% SPIKERT_SEQFIT - Perform sequential partitioning of variance explained for firing rate
%
% This function calculates the variance explained in firing rate (spikert)
% by forward velocity, angular velocity (only for ipsilateral turns), and
% a Gaussian fit based on visual object position. The function performs
% sequential model fitting using cross-validation to determine the unique
% and incremental contribution of each predictor to the variance in firing rate.
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
%   R2             - Array of R² values, in the following order:
%                      R2(1) - R² for forward velocity only (average across folds)
%                      R2(2) - Incremental R² for angular velocity (ipsilateral turns only)
%                      R2(3) - Incremental R² for Gaussian model (object position)
%                      R2(4) - Total R² for the combined model
%   mdlfit         - Structure containing final model fits on the full dataset:
%                      mdlfit.forward_full           - Final model for forward velocity only
%                      mdlfit.forward_angular_full   - Final model for forward + angular velocity (ipsilateral turns only)
%                      mdlfit.full_final             - Final model for combined predictors (forward, angular, Gaussian)
%
% CREATED: 11/07/2024 - MC
% UPDATED: 11/09/2024 - MC removed start/stops
%
function [R2,mdlfit] = spikert_seqfit(spikert, panelps, forward, angular, ttime, peak_position, spread)
%% Initialize

% Fetch processing settings
settings = processSettings;

% Adjust behavior for lag estimate
forward_lag = lagShift(forward, ttime, settings.fwdLag);
angular_lag = lagShift(angular, ttime, settings.angLag);

% Replace NaNs in visual_object_position with a placeholder value
placeholder_value = -999;
panelps(isnan(panelps)) = placeholder_value;

%% Pre-process data

% Optional: Exclude Start/Stop Transitions
% Set flag to exclude start/stop transitions using transition window settings
ex_startstop = 1;
postStartWin = 0.1; % Time window after start (in seconds)
preStopWin = 0.2;   % Time window before stop (in seconds)

% Number of trials (columns in the cell activity array)
nTrials = size(spikert, 2);

if ex_startstop
    % Convert post-start and pre-stop windows to indices based on time array
    postStartIdx = fetchTimeIdx(ttime, postStartWin);
    preStopIdx = fetchTimeIdx(ttime, preStopWin);

    % Loop over each trial
    for trial = 1:nTrials
        % Calculate run index using Schmitt Trigger
        runIdx = schmittTrigger(forward(:, trial), settings.runThreshE, 0.1);

        % Identify start and stop transitions in runIdx for the current trial
        runTransitions = diff(runIdx);    % Calculate transitions in run state
        startTrans = find(runTransitions == 1); % 0 to 1 (start running)
        stopTrans = find(runTransitions == -1); % 1 to 0 (stop running)

        % Loop over each start transition to set post-start period as NaN
        for st = 1:length(startTrans)
            tStart = startTrans(st); % Start index
            tEnd = min(size(spikert, 1), tStart + postStartIdx); % End index, within bounds
            spikert(tStart:tEnd, trial) = nan; % Set post-start window to NaN in cell activity data
        end

        % Loop over each stop transition to set pre-stop period as NaN
        for sp = 1:length(stopTrans)
            tStop = stopTrans(sp); % Stop index
            tStart = max(1, tStop - preStopIdx); % Start index, within bounds
            spikert(tStart:tStop, trial) = nan; % Set pre-stop window to NaN in cell activity data
        end
    end
end

% Reshape variables
firing_rate = reshape(spikert, [], 1);
forward_velocity = reshape(forward_lag, [], 1);
angular_velocity = reshape(angular_lag, [], 1);
visual_object_position = reshape(panelps, [], 1);

% Define valid indices (non-NaN) for sequential partitioning and positive angular velocity
valid_idx = ~any(isnan([forward_velocity, angular_velocity]), 2) & visual_object_position ~= placeholder_value & angular_velocity > 0;

% Separate valid data for the sequential partitioning with ipsilateral turns only
X_forward = forward_velocity(valid_idx);
X_angular = angular_velocity(valid_idx);
X_gaussian = visual_object_position(valid_idx);
y = firing_rate(valid_idx);

% Calculate a fixed baseline as the mean of spikert
baseline_estimate = mean(spikert(:), 'omitnan');

% Define the Truncated Gaussian model with fixed baseline, peak, and spread, truncated within [-15, 180]
gaussian_model = @(b, X) baseline_estimate + b * exp(-((X - peak_position).^2) / (2 * spread^2)) .* (X >= -15 & X <= 180);

% Fit the Gaussian model to estimate object position's influence
mdl_gaussian = fitnlm(X_gaussian, y, gaussian_model, 1);  % Initial guess for amplitude = 1
fitted_gaussian = mdl_gaussian.Fitted;

%% Cross-validation for Sequential Partitioning of Variance Explained

numFolds = 5;
cv = cvpartition(length(y), 'KFold', numFolds);

R2_folds_forward = zeros(numFolds, 1);
R2_folds_forward_angular = zeros(numFolds, 1);
R2_folds_full = zeros(numFolds, 1);

parfor i = 1:numFolds
    trainIdx = training(cv, i);
    testIdx = test(cv, i);

    % Step 1: Forward velocity only
    mdl_forward = fitlm(X_forward(trainIdx), y(trainIdx));
    R2_folds_forward(i) = mdl_forward.Rsquared.Ordinary;

    % Step 2: Forward + Angular velocity (ipsilateral turns only)
    mdl_forward_angular = fitlm([X_forward(trainIdx), X_angular(trainIdx)], y(trainIdx));
    R2_folds_forward_angular(i) = mdl_forward_angular.Rsquared.Ordinary;

    % Step 3: Forward + Angular velocity + Gaussian model (object position)
    mdl_full = fitlm([X_forward(trainIdx), X_angular(trainIdx), fitted_gaussian(trainIdx)], y(trainIdx));
    R2_folds_full(i) = mdl_full.Rsquared.Ordinary;
end

% Calculate average R-squared values across folds
R2_forward = mean(R2_folds_forward);
R2_forward_angular = mean(R2_folds_forward_angular);
R2_full = mean(R2_folds_full);

% Calculate incremental R² values
incremental_R2_angular = R2_forward_angular - R2_forward;
incremental_R2_gaussian = R2_full - R2_forward_angular;

%% Final Fit on Full Dataset
% Step 1: Forward velocity only
mdl_forward_full = fitlm(X_forward, y);

% Step 2: Forward + Angular velocity (ipsilateral turns only)
mdl_forward_angular_full = fitlm([X_forward, X_angular], y);

% Step 3: Forward + Angular velocity + Gaussian model (object position)
mdl_full_final = fitlm([X_forward, X_angular, fitted_gaussian], y);

% Store R² values directly in an array
R2 = [R2_forward, R2_forward_angular, R2_full];

% Store final model fits in a structure
mdlfit = struct();
mdlfit.forward_full = mdl_forward_full;
mdlfit.forward_angular_full = mdl_forward_angular_full;
mdlfit.full_final = mdl_full_final;
end