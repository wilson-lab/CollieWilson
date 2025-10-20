% pulse_v_output
%
% data processing function, for analyzing the relationship between each
% motion pulse and the corresponding output (e.g., angular velocity of the
% fly or spikerate of the cell)
%
% INPUT
% panelps - panel position
% forward - forward velocity
% varOut  - output variable to compare against
% pulseSelect - select index of pulse to analyze (e.g., 1 or 2)
% pulseOptions - possible motion pulses (e.g., 25 and 75)
% runSelect
%
% OUTPUT
% pulseTrials - trial data for each pulse
% pulseMeans  - mean data for each pulse
%
% CREATED   07/05/2024 - MC adapted from pulse_v_behavior_overlay
%

function [pulseTrials, pulseMeans] = pulse_v_output(panelps,forward,varOut,pulseSelect,pulseOptions,runSelect)
%% initialize
nTrials = size(panelps,2);
nPulses = length(pulseOptions);

% buffer before/after each pulse
bufferOptions = [1650 1000];
bufferWindow = bufferOptions(pulseSelect);

% for each pulse trial, set the percentage of time behavior must have been
% sufficient for said trial to be included in the average
minInclusion = 0.8;
% set number of pulse trials for average to be taken
minTrials = 3;

%% fetch behavioral index depending on run select
if runSelect==0 %quiescent only
    moveThresh = 0.1; %anything above/below considered moving
    runIdx = ~(forward<-moveThresh | forward>moveThresh);
elseif runSelect>0 %running only
    runIdx = schmittTrigger(forward,runSelect,0.1);
else %all
    runIdx = ones(size(forward));
end

%% pull motion pulses from pseudorandomized dataset
% for each trial
for t = 1:nTrials
    % pull this panel data
    thisTrial = panelps(:,t);
    try
        % run ordering function to find motion pulses in order
        [thisOrderR,thisOrderL] = order_motion_pulse(thisTrial,nPulses,pulseSelect,bufferWindow);
        % if any indices fall outside of trial duration, fix
        thisOrderR(thisOrderR>size(panelps,1)) = size(panelps,1);
        thisOrderL(thisOrderL>size(panelps,1)) = size(panelps,1);
        nSweep = size(thisOrderR,2);

        % pull and store data based on ordered motion pulse indices
        % rows = data, columns = trials, z = sweeps
        for p = 1:nSweep
            % store panel data
            pulsePanelpsR(:,t,p) = panelps(thisOrderR(:,p),t);
            pulsePanelpsL(:,t,p) = panelps(thisOrderL(:,p),t);
            % store output variable data
            pulseVarR(:,t,p) = varOut(thisOrderR(:,p),t);
            pulseVarL(:,t,p) = varOut(thisOrderL(:,p),t);
            % store index
            runIdxR(:,t,p) = runIdx(thisOrderR(:,p),t);
            runIdxL(:,t,p) = runIdx(thisOrderL(:,p),t);
        end
    catch
        % pull and store data based on ordered motion pulse indices
        % rows = data, columns = trials, z = sweeps
        for p = 1:nSweep
            % blank panel data
            pulsePanelpsR(:,t,p) = nan;
            pulsePanelpsL(:,t,p) = nan;
            % blank output variable data
            pulseVarR(:,t,p) = nan;
            pulseVarL(:,t,p) = nan;
            % blank index
            runIdxR(:,t,p) = 0;
            runIdxL(:,t,p) = 0;
        end
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
if nPulses>1 %motion pulse
    pulseMeans.varOutRL = mean([meanVarOutR -flip(meanVarOutL,3)],2,'omitnan');
else %stationary pulse
    pulseMeans.varOutRL = mean([meanVarOutR -flip(meanVarOutL,3)],2,'omitnan');
    % middle sweep at 0 should be omitted from R+L
    nSweepMid = round(nSweep/2);
    pulseMeans.varOutRL(:,:,nSweepMid) = meanVarOutR(:,:,nSweepMid);
    % sweeps on L in R+L are identical to sweeps on R, so blank
    pulseMeans.varOutRL(:,:,1:nSweepMid-1) = nan;
end

end