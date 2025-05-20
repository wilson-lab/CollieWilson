% spikert_xcorr
% This function determines the optimal lag between spike rates and behavior 
% by calculating the cross-correlation between spike rates and forward, angular, 
% and sideways velocities. The analysis allows for NaNs to handle gaps in the data.
%
% INPUT
% spikert - array of spike rates (time x trials)
% forward - array of forward velocities (time x trials)
% angular - array of angular velocities (time x trials)
% sideway - array of sideways velocities (time x trials)
% ttime   - time vector (in seconds) for the xcorr window
%
% OUTPUT
% r_val   - structure containing r values across the xcorr window for each velocity
% lag_t   - times (in milliseconds) corresponding to the cross-correlation lags
%
% CREATED 06/27/23 MC
% UPDATED 09/07/23 MC - Implemented Helen's xcorr analysis to handle NaNs
% UPDATED 10/16/23 MC - Concatenated data matrix instead of running xcorr per trial
% UPDATED 11/09/24 MC - Removing start/stop transitions
%

function [r_val, lag_t] = spikert_xcorr(spikert, forward, angular, sideway, ttime)
%% Initialize variables

% Fetch settings
settings = processSettings();

% Number of trials (columns in the spike rate array)
nTrials = size(spikert, 2);

% Define the total size of the xcorr window (e.g., ±1 second around zero lag)
xc_t = 2; % Cross-correlation window length in seconds
[xc_window] = fetchTimeIdx(ttime, xc_t); % Fetch indices corresponding to the xcorr window

%% Prepare data for cross-correlation

% Optional: Exclude Start/Stop Transitions
% Set flag to exclude start/stop transitions using transition window settings
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
        runTransitions = diff(runIdx);    % Calculate transitions in run state
        startTrans = find(runTransitions == 1); % 0 to 1 (start running)
        stopTrans = find(runTransitions == -1); % 1 to 0 (stop running)

        % Loop over each start transition to set post-start period as NaN
        for st = 1:length(startTrans)
            tStart = startTrans(st); % Start index
            tEnd = min(size(spikert, 1), tStart + postStartIdx); % End index, within bounds
            spikert(tStart:tEnd, trial) = nan; % Set post-start window to NaN in spike data
        end

        % Loop over each stop transition to set pre-stop period as NaN
        for sp = 1:length(stopTrans)
            tStop = stopTrans(sp); % Stop index
            tStart = max(1, tStop - preStopIdx); % Start index, within bounds
            spikert(tStart:tStop, trial) = nan; % Set pre-stop window to NaN in spike data
        end

        % Set spikert to NaN where runIdx is 0 (not running)
        spikert(runIdx == 0, trial) = nan;
    end
end

% % Optional: Fetch timepoints where forward velocity variance was high
% % Define parameters
% wsize = 500; % Window size for variance computation (adjust based on data)
% thresh = 1; % Threshold for significant forward velocity variance
% 
% % Loop through trials
% for trial = 1:nTrials
%     % Compute moving variance of forward velocity
%     varF = movvar(forward(:, trial), wsize,'omitnan');
% 
%     % Set forward velocity to NaN where variance is low
%     forward(varF <= thresh, trial) = NaN;
% end

% Filter for ipsilateral turns by excluding negative angular and sideways velocities
angular(angular < 0) = nan;
sideway(sideway < 0) = nan;

% Add a buffer of NaNs to the velocity and spike rate data to handle window edges
% The buffer ensures that xcorr computations do not include out-of-bounds data
forward_buff = [forward; nan(xc_window * 2, nTrials)];
angular_buff = [angular; nan(xc_window * 2, nTrials)];
sideway_buff = [sideway; nan(xc_window * 2, nTrials)];
spikert_buff = [spikert; nan(xc_window * 2, nTrials)];

% Concatenate the buffered data for each trial into a single vector
% This approach enables the cross-correlation to be computed across all trials at once
forward_cc = reshape(forward_buff, [], 1);
angular_cc = reshape(angular_buff, [], 1);
sideway_cc = reshape(sideway_buff, [], 1);
spikert_cc = reshape(spikert_buff, [], 1);

%% Run cross-correlation for each velocity type

% Perform cross-correlation between spike rate and each velocity type
% xcorrWGaps accounts for missing data (NaNs) in the inputs

disp('Running xcorr...');
[r_val.fwd, lag_idx] = xcorrWGaps(forward_cc, spikert_cc, xc_window); % Forward velocity xcorr
disp('Forward complete.');
[r_val.ang, ~] = xcorrWGaps(angular_cc, spikert_cc, xc_window); % Angular velocity xcorr
disp('Angular complete.');
[r_val.sid, ~] = xcorrWGaps(sideway_cc, spikert_cc, xc_window); % Sideways velocity xcorr
disp('Sideways complete.');

% Convert the lag index into time (in milliseconds) using the time step in ttime
lag_t = (lag_idx .* ttime(2)) * 1000; % Convert to milliseconds

%% Determine the optimal lag (optional, commented)

% Uncomment this section if you need to extract the optimal lag time (e.g., where r is maximal)
% The optimal lag corresponds to the peak of the cross-correlation function for each velocity type

% % Define a minimum peak prominence threshold for peak detection
% pk_prom = 0.08;
% 
% % Find peaks and their indices in the cross-correlation results
% [~, peak_idx_fwd] = findpeaks(r_val.fwd, 'SortStr', 'ascend', 'NPeaks', 1, 'MinPeakProminence', pk_prom);
% [~, peak_idx_ang] = findpeaks(r_val.ang, 'SortStr', 'ascend', 'NPeaks', 1, 'MinPeakProminence', pk_prom);
% [~, peak_idx_sid] = findpeaks(r_val.sid, 'SortStr', 'ascend', 'NPeaks', 1, 'MinPeakProminence', pk_prom);
% 
% % Extract the optimal lag time (in milliseconds) for each velocity type
% lag_opt = struct();
% if peak_idx_fwd
%     lag_opt.fwd = lag_t(peak_idx_fwd); % Optimal lag for forward velocity
% else
%     lag_opt.fwd = nan;
% end
% if peak_idx_ang
%     lag_opt.ang = lag_t(peak_idx_ang); % Optimal lag for angular velocity
% else
%     lag_opt.ang = nan;
% end
% if peak_idx_sid
%     lag_opt.sid = lag_t(peak_idx_sid); % Optimal lag for sideways velocity
% else
%     lag_opt.sid = nan;
% end

end
