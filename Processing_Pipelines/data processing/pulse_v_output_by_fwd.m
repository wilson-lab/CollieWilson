% pulse_v_output_by_fwd
% This function analyzes the relationship between a motion pulse and an output variable
% (e.g., angular velocity or spike rate), grouping pulse-aligned data based on the
% forward velocity of the animal into three behavioral categories: stationary, slow forward, and fast forward.
%
% The function identifies motion pulse events in each trial, extracts the corresponding data,
% categorizes timepoints by forward velocity, and averages responses across trials that meet
% a minimum inclusion threshold for each velocity group.
%
% INPUTS:
%   panelps      - Panel positions (degrees), time x trial
%   forward      - Forward velocity (mm/s), time x trial
%   varOut       - Output variable to analyze (same size as panelps)
%   ttime        - Time vector (seconds)
%   pulseSelect  - Index of the pulse to analyze (e.g., 1 or 2)
%   pulseOptions - Possible motion pulses (e.g., [25 75])
%   nSweep       - Number of sweeps to extract per trial
%
% OUTPUTS:
%   pulseTrials  - Structure containing pulse-aligned trial data for varOut, split by velocity group
%   pulseMeans   - Structure containing averaged varOut traces for each velocity group:
%                   - pulseMeans.stationary
%                   - pulseMeans.slow
%                   - pulseMeans.fast
%
% FORWARD VELOCITY GROUPS:
%   Group 1 - Stationary     : |forward| ≤ 0.5 mm/s
%   Group 2 - Slow Forward   : 0.5 < forward ≤ median(forward > 0)
%   Group 3 - Fast Forward   : forward > median(forward > 0)
%
% CREATED:
%   03/30/2025 - MC
%
function [pulseTrials, pulseMeans] = pulse_v_output_by_fwd(panelps,forward,varOut,ttime,pulseSelect,pulseOptions,nSweep,minFwd)
%% initialize
nTrials = size(panelps,2);
nPulses = length(pulseOptions);

% set analysis window
windowSize = 2; %size of full window
[winIdx] = fetchTimeIdx(ttime,windowSize);

% set number of pulse trials for average to be taken
minTrials = 2;

%% calculate categorical runIdx for behavioral state
runIdx = zeros(size(forward)); % initialize as 0 (uncategorized)
fwdNoise = 0.1; %noise floor

% Define stationary: abs(forward) <= fwdNoise
runIdx(abs(forward) <= fwdNoise) = 1;

% Get median of all positive forward velocities
medFwd = median(forward(forward > fwdNoise));%
%medFwd = minFwd;

% Define slow: 0 < forward <= median
runIdx((forward > fwdNoise) & (forward <= medFwd)) = 2;

% Define fast
%pursuitIdx = schmittTrigger(forward,3,0.1);
runIdx((forward > medFwd)) = 3;

% Minimum proportion of time points that must meet criteria to include trial
minInclusion = 0.4;

%% pull motion pulses from pseudorandomized dataset
% for each trial
c = 1; %counter
for t = 1:nTrials
    % pull this trial data
    thisPanelps = panelps(:,t);
    thisVarOut = varOut(:,t);
    try
        % run ordering function to find motion pulses in order
        [thisOrderR,thisOrderL] = order_motion_pulse(thisPanelps,nPulses,pulseSelect,winIdx);
        % if any indices fall outside of trial duration, fix
        thisOrderR(thisOrderR>size(thisPanelps,1)) = size(thisPanelps,1);
        thisOrderL(thisOrderL>size(thisPanelps,1)) = size(thisPanelps,1);
        % check to make sure analysis ran correctly
        if size(thisOrderR,2)==nSweep
            % pull and store data based on ordered motion pulse indices
            % rows = data, columns = trials, z = sweeps
            for p = 1:nSweep
                % store panel data
                pulsePanelpsR(:,c,p) = thisPanelps(thisOrderR(:,p));
                pulsePanelpsL(:,c,p) = thisPanelps(thisOrderL(:,p));
                % store output variable data
                pulseVarR(:,c,p) = thisVarOut(thisOrderR(:,p));
                pulseVarL(:,c,p) = thisVarOut(thisOrderL(:,p));
                % store index
                runIdxR(:,c,p) = runIdx(thisOrderR(:,p),t);
                runIdxL(:,c,p) = runIdx(thisOrderL(:,p),t);
            end
            c = c+1; %update counter
        end
    catch
        % pull and store data based on ordered motion pulse indices
        % rows = data, columns = trials, z = sweeps
        for p = 1:nSweep
            % blank panel data
            pulsePanelpsR(1:winIdx,c,p) = nan(winIdx,1);
            pulsePanelpsL(1:winIdx,c,p) = nan(winIdx,1);
            % blank output variable data
            pulseVarR(1:winIdx,c,p) = nan(winIdx,1);
            pulseVarL(1:winIdx,c,p) = nan(winIdx,1);
            % blank index
            runIdxR(1:winIdx,c,p) = zeros(winIdx,1);
            runIdxL(1:winIdx,c,p) = zeros(winIdx,1);
        end
        c = c+1; %update counter
    end
end


%% determine which trials met behavioral requirements (if any) for each forward group
pulseMeans = [];

% fetch pulse duration
groupLabels = {'stationary','slow','fast'};
groupCodes = [1, 2, 3];

% Define middle two-thirds indices
pulseDur = size(pulsePanelpsR,1);
midStart = round(pulseDur / 6);
midEnd   = round(5 * pulseDur / 6);
midIdx   = midStart:midEnd;

for g = 1:3
    group = groupLabels{g};
    code = groupCodes(g);

    for p = 1:nSweep
        % Check runIdx only within middle third
        idxR = sum(runIdxR(midIdx,:,p)==code,1) > round(length(midIdx) * minInclusion);
        idxL = sum(runIdxL(midIdx,:,p)==code,1) > round(length(midIdx) * minInclusion);

        % Store trial data
        pulseTrials.(group).varOutR{:,p} = pulseVarR(:,idxR,p);
        pulseTrials.(group).varOutL{:,p} = pulseVarL(:,idxL,p);

        % Compute means if enough trials
        if sum(idxR) >= minTrials
            pulseMeans.(group).varOutR(:,1,p) = mean(pulseVarR(:,idxR,p), 2);
        else
            pulseMeans.(group).varOutR(:,1,p) = nan(pulseDur,1);
        end

        if sum(idxL) >= minTrials
            pulseMeans.(group).varOutL(:,1,p) = mean(pulseVarL(:,idxL,p), 2);
        else
            pulseMeans.(group).varOutL(:,1,p) = nan(pulseDur,1);
        end
    end

    % Calculate RL average
    meanVarOutR = pulseMeans.(group).varOutR;
    meanVarOutL = pulseMeans.(group).varOutL;

    % Motion pulse check
    isItMoving = (max(meanVarOutR(:,1,1),[],'omitnan') - min(meanVarOutR(:,1,1),[],'omitnan')) > 5;

    if isItMoving
        pulseMeans.(group).varOutRL = mean([meanVarOutR -flip(meanVarOutL,3)],2);
    else
        pulseMeans.(group).varOutRL = mean([meanVarOutR -flip(meanVarOutL,3)],2);
        nSweepMid = round(nSweep/2);
        pulseMeans.(group).varOutRL(:,:,nSweepMid) = meanVarOutR(:,:,nSweepMid);
        pulseMeans.(group).varOutRL(:,:,1:nSweepMid-1) = nan;
    end
end

% calculate mean panel position and store
pulseMeans.panelpsR = mean(pulsePanelpsR,2,"omitnan");
pulseMeans.panelpsL = mean(pulsePanelpsL,2,"omitnan");

end