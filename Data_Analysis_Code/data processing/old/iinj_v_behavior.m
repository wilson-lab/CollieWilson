% iinj_v_behavior
% This analysis function retrieves and processes the behavioral responses from repeated current 
% injections. It calculates the mean responses for depolarizing and hyperpolarizing steps, allowing
% for analysis of the relationship between current injection and various output variables (e.g., 
% spike rate, velocity).
%
% INPUTS:
%   iinject     - Current stimulation data (amperes)
%   spikert     - Spike rate data (spikes/sec)
%   forward     - Forward velocity data (mm/s)
%   angular     - Angular velocity data (degrees/second)
%   sideway     - Sideways velocity data (degrees/second)
%   ttime       - Time vector (seconds)
%   runSelect   - Behavioral parameters (0 for all trials, >0 for thresholding)
%
% OUTPUTS:
%   depolMean   - Mean responses for depolarizing steps (structure containing various metrics)
%   hypolMean   - Mean responses for hyperpolarizing steps (structure containing various metrics)
%
% CREATED: 07/10/2023 - MC (adapted from spikerate_v_fullsweep2)
% UPDATED: 07/09/2024 - MC (simplified and combined with threshold version)
%
function [depolMean,hypolMean] = iinj_v_behavior(iinject,spikert,forward,angular,sideway,ttime,runSelect)
%% initialize
nTrials = size(iinject,2);

% set analysis window
preSize = 0.5; %size of window before
pstSize = 1.5; %size of window including and post
[~,preIdx] = min(abs(ttime - preSize)); %find nearest index
[~,pstIdx] = min(abs(ttime - pstSize)); %find nearest index

% for each pulse trial, set the percentage of time behavior must have been
% sufficient for said trial to be included in the average
minInclusion = 0.25;

%% fetch behavioral index depending on run select
if runSelect>0 %running only
    runIdx = schmittTrigger(forward,runSelect,0.1);
else %all
    runIdx = ones(size(forward));
end

%% find iinject pulses

% initialize
dpIinject = [];
dpSpikeRt = [];
dpForward = [];
dpAngular = [];
dpSideway = [];
dpRunIdx = [];
dp = 1;
hpIinject = [];
hpSpikeRt = [];
hpForward = [];
hpAngular = [];
hpSideway = [];
hpRunIdx = [];
hp = 1;

% for each trial
for nt = 1:nTrials
    % find all start/stops based on iinj output changes
    startstopidx = find(ischange(iinject(:,nt)));
    nSteps = length(startstopidx);

    % for each start step
    for ns = 1:2:nSteps
        thisStepVal = round(iinject(startstopidx(ns)+1,nt),-1); %step value
        thisStepIdx = startstopidx(ns)-preIdx:startstopidx(ns)+pstIdx; %step idx
        
        % for + pulse
        if thisStepVal>0
            dpIinject(:,dp) = iinject(thisStepIdx,nt);
            dpSpikeRt(:,dp) = spikert(thisStepIdx,nt);
            dpForward(:,dp) = forward(thisStepIdx,nt);
            dpAngular(:,dp) = angular(thisStepIdx,nt);
            dpSideway(:,dp) = sideway(thisStepIdx,nt);
            dpRunIdx(:,dp) = runIdx(thisStepIdx,nt);
            dp = dp+1;
        % for - pulse
        else
            hpIinject(:,hp) = iinject(thisStepIdx,nt);
            hpSpikeRt(:,hp) = spikert(thisStepIdx,nt);
            hpForward(:,hp) = forward(thisStepIdx,nt);
            hpAngular(:,hp) = angular(thisStepIdx,nt);
            hpSideway(:,hp) = sideway(thisStepIdx,nt);
            hpRunIdx(:,hp) = runIdx(thisStepIdx,nt);
            hp = hp+1;
        end
    end
end
% fetch step duration
stepDur = size(thisStepIdx,2);

%% determine which trials met behavioral requirements (if any) for this run
minIdx = round(stepDur*minInclusion);
goodTrials_dp = sum(dpRunIdx,1)>minIdx;
goodTrials_hp = sum(hpRunIdx,1)>minIdx;

%% calculate average responses and store for output
% for + pulse
depolMean.iinject = mean(dpIinject(:,goodTrials_dp),2,'omitnan');
depolMean.spikert = mean(dpSpikeRt(:,goodTrials_dp),2,'omitnan');
depolMean.forward = mean(dpForward(:,goodTrials_dp),2,'omitnan');
depolMean.angular = mean(dpAngular(:,goodTrials_dp),2,'omitnan');
depolMean.sideway = mean(dpSideway(:,goodTrials_dp),2,'omitnan');

% for - pulse
hypolMean.iinject = mean(hpIinject(:,goodTrials_hp),2,'omitnan');
hypolMean.spikert = mean(hpSpikeRt(:,goodTrials_hp),2,'omitnan');
hypolMean.forward = mean(hpForward(:,goodTrials_hp),2,'omitnan');
hypolMean.angular = mean(hpAngular(:,goodTrials_hp),2,'omitnan');
hypolMean.sideway = mean(hpSideway(:,goodTrials_hp),2,'omitnan');

end