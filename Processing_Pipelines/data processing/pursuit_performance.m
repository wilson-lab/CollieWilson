% pursuit_performance
% This function assesses the fly's pursuit performance across repeated presentations of the same stimulus,
% inspired by Sten Ruta (2021). It computes the fly's fidelity (correlation between target and turning)
% and vigor (net amount of turning in the direction of the target) over a sliding time window.
%
% INPUTS:
%   panelps   - Panel positions (degrees)
%   angular   - Angular velocity (degrees/second)
%   spikert   - Spike rate (Hz)
%   ttime     - Time vector (seconds)
%
% OUTPUTS:
%   fidelity  - Correlation between target position and turning (per time window)
%   vigor     - Net amount of turning in the direction of the target (per time window)
%   w_sr      - Spike rate per time window
%   w_time    - Time per window (mean time for each window)
%
% CREATED: 07/15/2024 - MC (adapted from osc_v_output)
%
function [fidelity,vigor,w_sr,w_time] = pursuit_performance(panelps,angular,spikert,ttime)
%% initialize
% remove panel data noise
panelps_r = round((panelps*2))/2; %round to nearest 0.5

% parameters
windowTime = 5; % s, size of the sliding window
stepTime = 0.05; % s, size of steps

[windowSize] = fetchTimeIdx(ttime,windowTime);
[stepSize] = fetchTimeIdx(ttime,stepTime);

%% calculate performance in sliding time-bin

% fetch number of trials
nTrials = size(panelps,2);

% initialize the vector to store correlation coefficients
numWindows = floor((length(angular) - windowSize) / stepSize) + 1;
corrCoeffs = zeros(numWindows,nTrials);
turnStrength = zeros(numWindows,nTrials);

% for each trial
for t = 1:nTrials
    % compute sliding window correlation
    for i = 1:numWindows
        % define the start and end indices of the window
        startIdx = (i - 1) * stepSize + 1;
        endIdx = startIdx + windowSize - 1;

        % extract the window data
        window_pos = panelps_r(startIdx:endIdx,t);
        window_ang = angular(startIdx:endIdx,t);
        window_sr = spikert(startIdx:endIdx,t);
        if t==1
            w_time(i) = mean(ttime(startIdx:endIdx));
        end

        % compute the correlation coefficient
        R = corrcoef(window_pos, window_ang);
        corrCoeffs(i,t) = R(1, 2);

        % compute the turn strength and fetch spikerate
        % if any timepoints where the target was on the right
        if any(window_pos>0)
            window_angR = window_ang;
            window_angR(window_ang<=0) = nan; % fetch only right turns
            strength_right = sum(window_angR(window_pos>0),'omitnan'); %fetch right turns when on right

            window_sr(window_pos<=0) = nan;
            w_sr(i,t) = mean(window_sr,'omitnan');
        else
            strength_right = 0;
        end
        % if any timepoints where the target was on the left
        if any(window_pos<0)
            window_angL = window_ang;
            window_angL(window_ang>=0) = nan; % fetch only left turns
            window_angL = window_angL *-1; %invert
            strength_left = sum(window_angL(window_pos<0),'omitnan'); %fetch left turns when on left
        else
            strength_left = 0;
        end
        turnStrength(i,t) = strength_right + strength_left;
    end
end

%% prepare outputs

% cap trackers and center according to window size
fidelity = corrCoeffs;
vigor = turnStrength;

end

