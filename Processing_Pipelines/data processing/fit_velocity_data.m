% fit_velocity_data
% This function fits linear models to raw forward, angular, and sideways velocity data across all trials for one animal.
% It applies preprocessing steps, including start/stop exclusion and optional lag adjustments,
% before combining data across trials and fitting a single linear model to each velocity type (forward, angular, and sideways).
%
% Inputs:
%   forward      - array of forward velocities (rows = time points, columns = trials for a single animal)
%   angular      - array of angular velocities (rows = time points, columns = trials for a single animal)
%   sideway      - array of sideways velocities (rows = time points, columns = trials for a single animal)
%   cellactivity - cell activity data, can be firing rate or voltage data (same dimensions as velocities)
%   ttime        - array of trial times (in seconds)
%   lagOpt       - flag (1 = apply lag shifts, 0 = omit lag shifts)
%
% Outputs:
%   slope_fwd - slope of the linear fit for forward velocity across all trials
%   slope_ang - slope of the linear fit for angular velocity (ipsilateral only) across all trials
%   slope_sid - slope of the linear fit for sideways velocity (ipsilateral only) across all trials

function [slope_fwd, slope_ang, slope_sid, r2_fwd, r2_ang, r2_sid] = fit_velocity_data(forward, angular, sideway, cellactivity, ttime, lagOpt)

    % Retrieve settings and initialize output variables
    settings = processSettings();
    nTrials = size(cellactivity, 2);
    
    %% Optional: Exclude Start/Stop Transitions
    ex_startstop = 1;
    postStartWin = 0.1; % Time window after start (in seconds)
    preStopWin = 0.2;   % Time window before stop (in seconds)

    if ex_startstop
        % Convert post-start and pre-stop windows to indices based on time array
        postStartIdx = fetchTimeIdx(ttime, postStartWin);
        preStopIdx = fetchTimeIdx(ttime, preStopWin);

        % Loop over each trial
        for trial = 1:nTrials
            % Calculate run index using Schmitt Trigger
            runIdx = schmittTrigger(forward(:, trial), settings.runThreshE, 0.1);

            % Identify start and stop transitions in runIdx for the current trial
            runTransitions = diff(runIdx);
            startTrans = find(runTransitions == 1); % 0 to 1 (start running)
            stopTrans = find(runTransitions == -1); % 1 to 0 (stop running)

            % Loop over each start transition to set post-start period as NaN
            for st = 1:length(startTrans)
                tStart = startTrans(st);
                tEnd = min(size(cellactivity, 1), tStart + postStartIdx);
                cellactivity(tStart:tEnd, trial) = nan;
            end

            % Loop over each stop transition to set pre-stop period as NaN
            for sp = 1:length(stopTrans)
                tStop = stopTrans(sp);
                tStart = max(1, tStop - preStopIdx);
                cellactivity(tStart:tStop, trial) = nan;
            end

            % Set cellactivity to NaN where runIdx is 0 (not running)
            cellactivity(runIdx == 0, trial) = nan;
        end
    end

    %% Optional: Apply Lag Shifts
    if lagOpt
        idxf = fetchTimeIdx(ttime, 0) - 1;
        idxa = fetchTimeIdx(ttime, settings.angLag) - 1;
        idxs = fetchTimeIdx(ttime, settings.sidLag) - 1;

        % Apply shifts and mark shifted portions as NaN
        forward = circshift(forward, -idxf, 1);
        forward(end-idxf+1:end, :) = nan;
        angular = circshift(angular, -idxa, 1);
        angular(end-idxa+1:end, :) = nan;
        sideway = circshift(sideway, -idxs, 1);
        sideway(end-idxs+1:end, :) = nan;
    end

    %% Combine All Trials and Fit Linear Models
    % Concatenate data across all trials for forward, angular, and sideways velocities
    forward_all = forward(:);
    forward_all(forward_all<2.5) = nan;
    angular_all = angular(:);
    sideway_all = sideway(:);
    cellactivity_all = cellactivity(:);

    % Fit linear model to combined forward velocity data
    validIdx_fwd = ~isnan(forward_all) & ~isnan(cellactivity_all);
    mdl_fwd = fitlm(forward_all(validIdx_fwd), cellactivity_all(validIdx_fwd));
    slope_fwd = mdl_fwd.Coefficients.Estimate(2); % Extract slope for forward fit
    r2_fwd = mdl_fwd.Rsquared.Ordinary;           % Extract R² for forward fit

    % Fit linear model to combined angular velocity data (ipsilateral only: angular >= 0)
    validIdx_ang = ~isnan(angular_all) & ~isnan(cellactivity_all) & (angular_all >= 0);
    mdl_ang = fitlm(angular_all(validIdx_ang), cellactivity_all(validIdx_ang));
    slope_ang = mdl_ang.Coefficients.Estimate(2); % Extract slope for angular fit
    r2_ang = mdl_ang.Rsquared.Ordinary;           % Extract R² for angular fit

    % Fit linear model to combined sideways velocity data (ipsilateral only: sideway >= 0)
    validIdx_sid = ~isnan(sideway_all) & ~isnan(cellactivity_all) & (sideway_all >= 0);
    mdl_sid = fitlm(sideway_all(validIdx_sid), cellactivity_all(validIdx_sid));
    slope_sid = mdl_sid.Coefficients.Estimate(2); % Extract slope for sideways fit
    r2_sid = mdl_sid.Rsquared.Ordinary;           % Extract R² for sideways fit

end
