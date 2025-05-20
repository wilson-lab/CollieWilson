% BEHAVIOR_V_PANELPS_LAG - This function estimates the optimal lag between an object's position 
% (panelps) and the corresponding behavior (angular) for a given time window. It does so by 
% calculating the cross-correlation (xcorr) between the two signals, identifying the most prominent 
% peak in the cross-correlation to estimate the ideal lag. The lag is expressed in milliseconds (ms).

% Inputs:
%   panelps   - 3D array containing object position data, where each slice along the 3rd dimension
%               represents a different condition, and each column represents a trial.
%   angular   - 3D array containing the corresponding angular velocity data, with dimensions matching
%               panelps.
%   ttime     - Time vector, representing the time points at which the data were sampled.
%   settings  - Structure containing parameters, including 'minXCorrProm', the minimum peak prominence
%               threshold for cross-correlation.

% Outputs:
%   opt_lag   - Estimated optimal lag (in milliseconds) for each condition based on the time of the
%               most prominent peak in the cross-correlation.
%   r_out     - Cross-correlation results for each condition, stored as columns where each column 
%               represents the cross-correlation for a different condition.
%
function [opt_lag,r_out] = behavior_v_panelps_lag(panelps,angular,ttime,settings)
%% initialize
nCond = size(panelps,3);
nTrial = size(panelps,2);

% set size of xcorr window
xc_t = 0.75; %s, total (e.g., 1/2 on each side of 0)
[xc_window] = fetchTimeIdx(ttime,xc_t);

%% prepare data

% add bufer window
buffPanelps = [panelps; nan(xc_window*2,nTrial,nCond)];
buffAngular = [angular; nan(xc_window*2,nTrial,nCond)];

%% estimate ideal lag between object position and behavior
peakLocation = [];

% for each condition
parfor c = 1:nCond
    thisPanelps = reshape(buffPanelps(:,:,c),[],1);
    thisAngular = reshape(buffAngular(:,:,c),[],1);
    [r_val, r_lag] = xcorrWGaps(thisPanelps,thisAngular,xc_window);

    % find most prominent peak
    [~, locs, ~, p] = findpeaks(r_val, r_lag, 'Annotate', 'extents','MinPeakProminence',settings.minXCorrProm);
    if ~isempty(p)
        [~, maxPromIdx] = max(p);
        peakLocation(c) = locs(maxPromIdx);
        r_out(:,c) = r_val;
    else
        peakLocation(c) = nan;
        r_out(:,c) = nan(1,length(r_val));
    end
    
end
opt_lag = (ttime(2)*peakLocation).*1000; %ms

end