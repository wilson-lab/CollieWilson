% pulse_v_output
% This data processing function analyzes the relationship between each motion pulse 
% and the corresponding output variable (e.g., angular velocity of the fly or spike rate of a neuron). 
% It extracts and averages the data for each pulse based on behavioral parameters.
%
% INPUTS:
%   panelps      - Panel positions (degrees)
%   forward      - Forward velocity (mm/s)
%   varOut       - Output variable to compare against (e.g., angular velocity, spike rate)
%   ttime        - Time vector (seconds)
%   pulseSelect  - Index of the pulse to analyze (e.g., 1 or 2)
%   pulseOptions - Possible motion pulses (e.g., 25 and 75)
%   nSweep       - Number of sweeps to analyze
%   runSelect    - Behavioral parameter selection (0: quiescent, 1: running, -1: all)
%
% OUTPUTS:
%   pulseTrials  - Trial data for each pulse (structure containing output variable data)
%   pulseMeans   - Mean data for each pulse (structure containing mean values)
%
% CREATED:
%   07/05/2024 - MC (adapted from pulse_v_behavior_overlay)
%   08/01/2024 - MC (adjusted buffer window determination)
%
function [pulseTrials, pulseMeans] = pulse_v_output(panelps,forward,varOut,ttime,pulseSelect,pulseOptions,nSweep,runSelect)
%% initialize
nTrials = size(panelps,2);
nPulses = length(pulseOptions);

% set analysis window
windowSize = 2; %size of full window
[winIdx] = fetchTimeIdx(ttime,windowSize);

% set number of pulse trials for average to be taken
minTrials = 3;

%% fetch behavioral index depending on run select
if runSelect==0 %quiescent only
    moveThresh = 0.5; %anything above/below considered moving
    runIdx = ~(forward<-moveThresh | forward>moveThresh);
    % for each pulse trial, set the percentage of time behavior must have been
    % sufficient for said trial to be included in the average
    minInclusion = 0.5;
elseif runSelect>0 %running only
    runIdx = schmittTrigger(forward,runSelect,0.1);
    minInclusion = 0.75;
else %all
    runIdx = ones(size(forward));
    minInclusion = 0.75;
end

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
% fetch pulse duration
pulseDur = size(pulsePanelpsR,1);

%% determine which trials met behavioral requirements (if any) for this run
minIdx = round(pulseDur*minInclusion);
goodTrialsR = sum(runIdxR,1)>minIdx;
goodTrialsL = sum(runIdxL,1)>minIdx;

%% prepare outputs
% calculate mean panel position and store
pulseMeans.panelpsR = mean(pulsePanelpsR,2,"omitnan");
pulseMeans.panelpsL = mean(pulsePanelpsL,2,"omitnan");
% check if its moving
isItMoving = (max(pulseMeans.panelpsR(:,1,1)) - min(pulseMeans.panelpsR(:,1,1)))>5;

% for each sweep position, calculate mean of included trials and store
for p = 1:nSweep
    % for rightward sweeps
    idxR = goodTrialsR(1,:,p);
    % fetch trials
    trialsVarOutR{:,p} = pulseVarR(:,idxR,p);
    % if sufficient, calculate mean
    if sum(idxR)>=minTrials
        meanVarOutR(:,1,p) = mean(pulseVarR(:,idxR,p),2);
    else
        meanVarOutR(:,1,p) = nan(pulseDur,1);
    end
    % for leftward sweeps
    idxL = goodTrialsL(1,:,p);
    % fetch trials
    trialsVarOutL{:,p} = pulseVarL(:,idxR,p);
    % if sufficient, calculate mean
    if sum(idxL)>=minTrials
        meanVarOutL(:,1,p) = mean(pulseVarL(:,idxL,p),2);
    else
        meanVarOutL(:,1,p) = nan(pulseDur,1);
    end
end
% store trials
pulseTrials.varOutR = trialsVarOutR;
pulseTrials.varOutL = trialsVarOutL;

% store means
pulseMeans.varOutR = meanVarOutR;
pulseMeans.varOutL = meanVarOutL;
if isItMoving %motion pulse
    pulseMeans.varOutRL = mean([meanVarOutR -flip(meanVarOutL,3)],2);
else %stationary pulse
    pulseMeans.varOutRL = mean([meanVarOutR -flip(meanVarOutL,3)],2);
    % middle sweep at 0 should be omitted from R+L
    nSweepMid = round(nSweep/2);
    pulseMeans.varOutRL(:,:,nSweepMid) = meanVarOutR(:,:,nSweepMid);
    % sweeps on L in R+L are identical to sweeps on R, so blank
    pulseMeans.varOutRL(:,:,1:nSweepMid-1) = nan;
end

end