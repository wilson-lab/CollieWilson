% iinj_v_behavior
% This function analyzes the relationship between injected current pulses 
% and behavioral responses. It categorizes trials based on depolarizing, 
% hyperpolarizing, and control pulses, filtering trials based on behavioral 
% criteria, and computing mean response values for each category.
%
% INPUTS:
% - iinject: Injected current data (time x trials)
% - spikert: Spike rate data (time x trials)
% - forward: Forward velocity data (time x trials)
% - angular: Angular velocity data (time x trials)
% - sideway: Sideway velocity data (time x trials)
% - ttime: Time vector
% - runSelect: Threshold for selecting running trials (0 for all trials)
%
% OUTPUTS:
% - depolMean: Mean response values for depolarizing pulses
% - hypolMean: Mean response values for hyperpolarizing pulses
% - controlMean: Mean response values for control trials
% - responseStats: Summary statistics for frequency and max response
%
% CREATED: 03/04/2025 - MC
% UPDATED: 03/10/2025 - MC updated for optional hyperpolarization
%          03/11/2025 - MC added pre-pulse bins
%          03/17/2925 - MC added high inclusion analyses
%
function [depolMean, hypolMean, controlMean, responseStats] = iinj_v_behavior(iinject, spikert, forward, angular, sideway, ttime, runSelect)
%% Initialize variables
nTrials = size(iinject, 2);

% Define analysis window
buffDur = 0.5; % Time buffer before and after pulse
pulseDur = 1; % Time duration of pulse
ctrStart = 3; % Time buffer after pulse to start control block

% Set analysis window relative to pulse onset
preIdx = fetchTimeIdx(ttime, buffDur); % Index for pre-pulse window
pstIdx = fetchTimeIdx(ttime, pulseDur + buffDur); % Index for post-pulse window
ctrIdx = fetchTimeIdx(ttime, ctrStart); % Index for control window after a pulse window

%% Determine run index based on selection criteria
if runSelect > 0 % Running trials only
    runIdx = schmittTrigger(forward, runSelect, 0.1);
else % Include all trials
    runIdx = ones(size(forward));
end

%% Identify iinject pulses and categorize trials
% Initialize storage variables
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
ctrIinject = [];
ctrSpikeRt = [];
ctrForward = [];
ctrAngular = [];
ctrSideway = [];
ctrRunIdx = [];
ctr = 1;

% Process each trial
for nt = 1:nTrials
    % Identify changes in iinject to find pulse start and stop
    startstopidx = find(ischange(iinject(:, nt)));
    nSteps = length(startstopidx);

    % Process each pulse event
    for ns = 1:2:nSteps
        thisStepVal = round(iinject(startstopidx(ns) + 1, nt), -1); % Pulse value
        thisStepIdx = startstopidx(ns) - preIdx : startstopidx(ns) + pstIdx; % Time window
        thisCtrlIdx = thisStepIdx + ctrIdx; % Control time window
        
        % Store data based on pulse polarity
        if thisStepVal > 0 % Depolarizing pulse
            dpIinject(:, dp) = iinject(thisStepIdx, nt);
            dpSpikeRt(:, dp) = spikert(thisStepIdx, nt);
            dpForward(:, dp) = forward(thisStepIdx, nt);
            dpAngular(:, dp) = angular(thisStepIdx, nt);
            dpSideway(:, dp) = sideway(thisStepIdx, nt);
            dpRunIdx(:, dp) = runIdx(thisStepIdx, nt);
            dp = dp + 1;
        else % Hyperpolarizing pulse
            hpIinject(:, hp) = iinject(thisStepIdx, nt);
            hpSpikeRt(:, hp) = spikert(thisStepIdx, nt);
            hpForward(:, hp) = forward(thisStepIdx, nt);
            hpAngular(:, hp) = angular(thisStepIdx, nt);
            hpSideway(:, hp) = sideway(thisStepIdx, nt);
            hpRunIdx(:, hp) = runIdx(thisStepIdx, nt);
            hp = hp + 1;
        end
        % Control pulse
        if max(thisCtrlIdx)<length(iinject)
            ctrIinject(:, ctr) = iinject(thisCtrlIdx, nt);
            ctrSpikeRt(:, ctr) = spikert(thisCtrlIdx, nt);
            ctrForward(:, ctr) = forward(thisCtrlIdx, nt);
            ctrAngular(:, ctr) = angular(thisCtrlIdx, nt);
            ctrSideway(:, ctr) = sideway(thisCtrlIdx, nt);
            ctrRunIdx(:, ctr) = runIdx(thisCtrlIdx, nt);
            ctr = ctr + 1;
        end
    end
end

%% Determine which trials met behavioral requirements
% Set inclusion parameters
minInclusion = 0.3;
minIdx = round(size(dpRunIdx, 1) * minInclusion);
% Find trials for each condition that meet inclusion parameters
goodTrials_dp = sum(dpRunIdx, 1) > minIdx;
goodTrials_hp = sum(hpRunIdx, 1) > minIdx;
goodTrials_ctr = sum(ctrRunIdx, 1) > minIdx;

% Fetch trials for each condition that meet inclusion parameters
dpIinject_inc = dpIinject(:, goodTrials_dp);
dpSpikeRt_inc = dpSpikeRt(:, goodTrials_dp);
dpForward_inc = dpForward(:, goodTrials_dp);
dpAngular_inc = dpAngular(:, goodTrials_dp);
dpSideway_inc = dpSideway(:, goodTrials_dp);

hpIinject_inc = hpIinject(:, goodTrials_hp);
hpSpikeRt_inc = hpSpikeRt(:, goodTrials_hp);
hpForward_inc = hpForward(:, goodTrials_hp);
hpAngular_inc = hpAngular(:, goodTrials_hp);
hpSideway_inc = hpSideway(:, goodTrials_hp);

ctrIinject_inc = ctrIinject(:, goodTrials_ctr);
ctrSpikeRt_inc = ctrSpikeRt(:, goodTrials_ctr);
ctrForward_inc = ctrForward(:, goodTrials_ctr);
ctrAngular_inc = ctrAngular(:, goodTrials_ctr);
ctrSideway_inc = ctrSideway(:, goodTrials_ctr);

%% Partition data according to change in firing rate for hyperpolarization trials

% Set quartile to include
qVal= 0.25;

% Filter depolarizing pulses
% Compute change in firing rate for each trial
fr_change = mean(dpSpikeRt_inc(preIdx:pstIdx,:),1) - mean(dpSpikeRt_inc(1:preIdx,:),1);

% Determine top quartile threshold
qSelect = quantile(fr_change, qVal);

% Logical index of trials in top quartile
top_quartile_trials = fr_change >= qSelect;

% Pull corresponding trials for all variables
dpIinject_inc = dpIinject_inc(:, top_quartile_trials);
dpSpikeRt_inc = dpSpikeRt_inc(:, top_quartile_trials);
dpForward_inc = dpForward_inc(:, top_quartile_trials);
dpAngular_inc = dpAngular_inc(:, top_quartile_trials);
dpSideway_inc = dpSideway_inc(:, top_quartile_trials);

% If available, filter hyperpolarizing pulses
if ~isempty(hpIinject_inc)
    % Compute change in firing rate for each trial
    fr_change = mean(hpSpikeRt_inc(1:preIdx,:),1) - mean(hpSpikeRt_inc(preIdx:pstIdx,:),1);

    % Determine top quartile threshold
    qSelect = quantile(fr_change, qVal);

    % Logical index of trials in top quartile
    top_quartile_trials = fr_change >= qSelect;

    % Pull corresponding trials for all variables
    hpIinject_inc = hpIinject_inc(:, top_quartile_trials);
    hpSpikeRt_inc = hpSpikeRt_inc(:, top_quartile_trials);
    hpForward_inc = hpForward_inc(:, top_quartile_trials);
    hpAngular_inc = hpAngular_inc(:, top_quartile_trials);
    hpSideway_inc = hpSideway_inc(:, top_quartile_trials);
end

%% Determine which trials met HIGH behavioral requirements
% Set HIGH inclusion parameters
minInclusion_high = 0.8;
minIdx_high = round(size(dpRunIdx, 1) * minInclusion_high);

% Find trials for each condition that meet HIGH inclusion parameters
goodTrials_dph = sum(dpRunIdx, 1) > minIdx_high;

% Fetch trials for each condition that meet HIGH inclusion parameters
dpIinject_incH = dpIinject(:, goodTrials_dph);
dpSpikeRt_incH = dpSpikeRt(:, goodTrials_dph);
dpForward_incH = dpForward(:, goodTrials_dph);
dpAngular_incH = dpAngular(:, goodTrials_dph);
dpSideway_incH = dpSideway(:, goodTrials_dph);

%% Partition data according to forward velocity being above or below median
% Determine average forward velocity prior to pulse
mean_prepulse_forward = mean(dpForward_inc(1:preIdx, :), 1, 'omitnan');

% Compute the median of average forward velocities
median_prepulse_forward = median(mean_prepulse_forward, 'omitnan');

% Find trials where forward velocity is below or above the median
below_median_trials = mean_prepulse_forward < median_prepulse_forward;
above_median_trials = mean_prepulse_forward >= median_prepulse_forward;

% Fetch trials for each condition where forward velocity is below median
dpIinject_below = dpIinject_inc(:, below_median_trials);
dpSpikeRt_below = dpSpikeRt_inc(:, below_median_trials);
dpForward_below = dpForward_inc(:, below_median_trials);
dpAngular_below = dpAngular_inc(:, below_median_trials);
dpSideway_below = dpSideway_inc(:, below_median_trials);

% Fetch trials for each condition where forward velocity is above median
dpIinject_above = dpIinject_inc(:, above_median_trials);
dpSpikeRt_above = dpSpikeRt_inc(:, above_median_trials);
dpForward_above = dpForward_inc(:, above_median_trials);
dpAngular_above = dpAngular_inc(:, above_median_trials);
dpSideway_above = dpSideway_inc(:, above_median_trials);

%% Analyze responses for each pulse subset
% Run analyses for each condition
% Assign remaining variables for output
responseStats = struct;

% Analyze depolarizing pulses
[responseStats.dpFreq, dpTurnStats] = analyze_iinj(dpIinject_inc, dpSpikeRt_inc, dpForward_inc, dpAngular_inc, dpSideway_inc, buffDur, pulseDur, ttime);
responseStats.dpSpikert = dpTurnStats.spikert;
responseStats.dpAngular = dpTurnStats.angular;

% Analyze hyperpolarizing pulses (if available)
if isempty(hpIinject_inc)
    responseStats.hpFreq = nan;
    responseStats.hpSpikert = nan;
    responseStats.hpAngular = nan;
else
    [responseStats.hpFreq, hpTurnStats] = analyze_iinj(hpIinject_inc, hpSpikeRt_inc, hpForward_inc, -hpAngular_inc, hpSideway_inc, buffDur, pulseDur, ttime);
    responseStats.hpSpikert = hpTurnStats.spikert;
    responseStats.hpAngular = -hpTurnStats.angular;
end

% Analyze control segments
[responseStats.ctrFreq, ctrTurnStats] = analyze_iinj(ctrIinject_inc, ctrSpikeRt_inc, ctrForward_inc, ctrAngular_inc, ctrSideway_inc, buffDur, pulseDur, ttime);
responseStats.ctrSpikert = ctrTurnStats.spikert;
responseStats.ctrAngular = ctrTurnStats.angular;

%% Additional analyses for binned pulses

% Analyze depolarizing pulses with HIGH inclusion criteria
[responseStats.dphFreq, dphTurnStats] = analyze_iinj(dpIinject_incH, dpSpikeRt_incH, dpForward_incH, dpAngular_incH, dpSideway_incH, buffDur, pulseDur, ttime);
responseStats.dphSpikert = dphTurnStats.spikert;
responseStats.dphAngular = dphTurnStats.angular;

% Analyze depolarizing pulses for forward velocity BELOW median
[responseStats.dpbFreq, dpbTurnStats] = analyze_iinj(dpIinject_below, dpSpikeRt_below, dpForward_below, dpAngular_below, dpSideway_below, buffDur, pulseDur, ttime);
responseStats.dpbSpikert = dpbTurnStats.spikert;
responseStats.dpbAngular = dpbTurnStats.angular;
% Analyze depolarizing pulses for forward velocity ABOVE median
[responseStats.dpaFreq, dpaTurnStats] = analyze_iinj(dpIinject_above, dpSpikeRt_above, dpForward_above, dpAngular_above, dpSideway_above, buffDur, pulseDur, ttime);
responseStats.dpaSpikert = dpaTurnStats.spikert;
responseStats.dpaAngular = dpaTurnStats.angular;

%% Compute mean responses for included trials
% Depolarizing pulse
depolMean.iinject = mean(dpIinject_inc, 2, 'omitnan');
depolMean.spikert = mean(dpSpikeRt_inc, 2, 'omitnan');
depolMean.forward = mean(dpForward_inc, 2, 'omitnan');
depolMean.angular = mean(dpAngular_inc, 2, 'omitnan');
depolMean.sideway = mean(dpSideway_inc, 2, 'omitnan');

% Hyperpolarizing pulse
if isempty(hpIinject_inc)
    hypolMean.iinject = nan(size(dpIinject_inc, 1), 1);
    hypolMean.spikert = nan(size(dpSpikeRt_inc, 1), 1);
    hypolMean.forward = nan(size(dpForward_inc, 1), 1);
    hypolMean.angular = nan(size(dpAngular_inc, 1), 1);
    hypolMean.sideway = nan(size(dpSideway_inc, 1), 1);
else
    hypolMean.iinject = mean(hpIinject_inc, 2, 'omitnan');
    hypolMean.spikert = mean(hpSpikeRt_inc, 2, 'omitnan');
    hypolMean.forward = mean(hpForward_inc, 2, 'omitnan');
    hypolMean.angular = mean(hpAngular_inc, 2, 'omitnan');
    hypolMean.sideway = mean(hpSideway_inc, 2, 'omitnan');
end

% Control pulse
controlMean.iinject = mean(ctrIinject_inc, 2, 'omitnan');
controlMean.spikert = mean(ctrSpikeRt_inc, 2, 'omitnan');
controlMean.forward = mean(ctrForward_inc, 2, 'omitnan');
controlMean.angular = mean(ctrAngular_inc, 2, 'omitnan');
controlMean.sideway = mean(ctrSideway_inc, 2, 'omitnan');

%% Compute mean responses for binned included trials
% Depolarizing pulse w/ HIGH threshold
depolMean.iinjectH = mean(dpIinject_incH, 2, 'omitnan');
depolMean.spikertH = mean(dpSpikeRt_incH, 2, 'omitnan');
depolMean.forwardH = mean(dpForward_incH, 2, 'omitnan');
depolMean.angularH = mean(dpAngular_incH, 2, 'omitnan');
depolMean.sidewayH = mean(dpSideway_incH, 2, 'omitnan');

% Depolarizing pulse w/ BELOW forward median threshold
depolMean.iinjectB = mean(dpIinject_below, 2, 'omitnan');
depolMean.spikertB = mean(dpSpikeRt_below, 2, 'omitnan');
depolMean.forwardB = mean(dpForward_below, 2, 'omitnan');
depolMean.angularB = mean(dpAngular_below, 2, 'omitnan');
depolMean.sidewayB = mean(dpSideway_below, 2, 'omitnan');
% Depolarizing pulse w/ ABOVE forward median threshold
depolMean.iinjectA = mean(dpIinject_above, 2, 'omitnan');
depolMean.spikertA = mean(dpSpikeRt_above, 2, 'omitnan');
depolMean.forwardA = mean(dpForward_above, 2, 'omitnan');
depolMean.angularA = mean(dpAngular_above, 2, 'omitnan');
depolMean.sidewayA = mean(dpSideway_above, 2, 'omitnan');

%% Bin according to pre-pulse behavior
% Fetch trials where pre-pulse fly turned ipsi
ipsiTrials = mean(dpAngular_inc(1:preIdx,:))>10;

depolMean.iinject_ipsi = mean(dpIinject_inc(:,ipsiTrials), 2, 'omitnan');
depolMean.spikert_ipsi = mean(dpSpikeRt_inc(:,ipsiTrials), 2, 'omitnan');
depolMean.forward_ipsi = mean(dpForward_inc(:,ipsiTrials), 2, 'omitnan');
depolMean.angular_ipsi = mean(dpAngular_inc(:,ipsiTrials), 2, 'omitnan');
depolMean.sideway_ipsi = mean(dpSideway_inc(:,ipsiTrials), 2, 'omitnan');

% Fetch trials where pre-pulse fly turned ipsi
contraTrials = mean(dpAngular_inc(1:preIdx,:))<-10;

depolMean.iinject_contra = mean(dpIinject_inc(:,contraTrials), 2, 'omitnan');
depolMean.spikert_contra = mean(dpSpikeRt_inc(:,contraTrials), 2, 'omitnan');
depolMean.forward_contra = mean(dpForward_inc(:,contraTrials), 2, 'omitnan');
depolMean.angular_contra = mean(dpAngular_inc(:,contraTrials), 2, 'omitnan');
depolMean.sideway_contra = mean(dpSideway_inc(:,contraTrials), 2, 'omitnan');

end