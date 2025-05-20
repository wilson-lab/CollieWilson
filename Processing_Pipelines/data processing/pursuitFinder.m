% pursuitFinder.m
% Function used to pull only pursuit epochs by first thresholding with a
% Schmitt trigger and then by thresholding based on epoch duration
%
% INPUTS:
% allForward - forward velocity array containing all data
% allAngular - angular velocity array containing all data
% allSideway - sidways velocity array containing all data
% allSpikeRt - spike rate array containing all data
% expttime - time array for 1 trial
% pursuitThresh - minimum forward velocity for "pursuit"
% 
% INPUTS:
% purForward - forward velocity array containing only pursuit, others nan
% purAngular - angular velocity array containing only pursuit, others nan
% purSideway - sidways velocity array containing only pursuit, others nan
% purSpikeRt - spike rate array containing only pursuit, others nan
%
% CREATED: 12/12/2022 MC
% UPDATED: 12/13/2022 MC made pipeline components optional
%

function [purForward,purAngular,purSideway,purSpikeRt] = pursuitFinder(allForward,allAngular,allSideway,allSpikeRt,expttime,pursuitThresh)
%% initialize

% check which data is included
includeBehave = any(allAngular,'all'); %true if behavior not zero
includeSpikeRt = any(allSpikeRt,'all'); %true if spikerate not zero

% copy data
purForward = allForward;
purAngular = allAngular;
purSideway = allSideway;
purSpikeRt = allSpikeRt;

nTrials = size(allForward,2);
trialDur = size(allForward,1);

% set min pursuit time
minPursuitTime = 5; %sec


%% threshold for pursuit epochs based on forward speed

% find pursuit indices using schmitt trigger
pursuit_idx = schmittTrigger(allForward,pursuitThresh,0.1);
% blank out non pursuit behavior
if includeBehave
    purForward(~pursuit_idx) = NaN;
    purAngular(~pursuit_idx) = NaN;
    purSideway(~pursuit_idx) = NaN;
end
if includeSpikeRt
    purSpikeRt(~pursuit_idx) = NaN;
end


%% threshold for pursuit epochs based on run duration

% convert min pursuit time to sample rate
minPursuitN = minPursuitTime * find(expttime==1)-1;

% for each trial...
for t = 1:nTrials
    % find pursuit start/stops for this trial
    pursuit_diff = diff(pursuit_idx(:,t));
    start_idx = find(pursuit_diff==1);
    stop_idx = find(pursuit_diff==-1);

    % if fly ends trial with pursuit (no last stop)
    if length(start_idx)>length(stop_idx)
        stop_idx = [stop_idx ; trialDur];
    % if fly starts trial with pursuit (no first start)
    elseif length(start_idx)<length(stop_idx)
        start_idx = [1 ; start_idx];
    end

    % for each stop/start pair
    pursuit_duration = stop_idx-start_idx; 
    for p = 1:length(pursuit_duration)
        % if this epoch is shorter than minimum
        if pursuit_duration(p)<minPursuitN
            % blank out epochs that are too short
            if includeBehave
                purForward(start_idx(p):stop_idx(p),t) = NaN;
                purAngular(start_idx(p):stop_idx(p),t) = NaN;
                purSideway(start_idx(p):stop_idx(p),t) = NaN;
            end
            if includeSpikeRt
                purSpikeRt(start_idx(p):stop_idx(p),t) = NaN;
            end
        end
    end
end

end