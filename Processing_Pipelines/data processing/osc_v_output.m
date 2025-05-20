% osc_v_output
% This analysis function generates a summary plot of panel position versus 
% a given output variable (e.g., velocity or firing rate). It processes 
% repeated presentations of target sweeps, calculating mean data across 
% sweeps while accounting for behavioral parameters.
%
% INPUTS:
%   panelps   - Downsampled panel positions (degrees)
%   forward   - Downsampled forward velocities (mm/s)
%   varOut    - Output variable to compare against (e.g., velocity or firing rate)
%   runSelect  - Set behavioral parameters (0 for quiescent, >0 for running, -1 for all)
%
% OUTPUTS:
%   sweepPos  - Position of the target sweep (mean position)
%   sweepMean - Mean data across target sweeps (mean of output variable)
%
% CREATED:
%   07/08/2024 - MC (adapted from spikerate_v_fullsweep2)
%
function [sweepPos,sweepMean] = osc_v_output(panelps,forward,varOut,runSelect,optLag,ttime,settings)
%% initialize
nTrials = size(panelps,2);

% for each repeated sweep, set the percentage of time behavior must have been
% sufficient for said repetition to be included in the average
minInclusion = 0.75;

% remove panel data noise
panelps_r = round((panelps*2))/2; %round to nearest 0.5

%% fetch behavioral index depending on run select
if runSelect==0 %quiescent only
    moveThresh = 0.25; %anything above/below considered moving
    runIdx = ~(forward<-moveThresh | forward>moveThresh);
elseif runSelect>0 %running only
    runIdx = schmittTrigger(forward,runSelect,0.1);
else %all
    runIdx = ones(size(forward));
end

%% (optional) shift according to lag estimates
% if lag estimates were provided, shift
if optLag
    % fetch shift indices for each lag
    [idx_vm] = fetchTimeIdx(ttime,settings.visuomotorLag);
    idx_vm = idx_vm-1;
    
    % shift and exclude data at start/stop of trial
    varOut = circshift(varOut,+idx_vm,1); %shift
    varOut(end-idx_vm+1:end,:) = nan; %exclude shifts
end

%% pull output variable for each repeated sweep sweep

% for each trial
c = 1; %counter
for nt = 1:nTrials
    thisPanelps = panelps_r(:,nt);
    thisVarOut = varOut(:,nt);
    % find where the oscillating sweep crosses the midline (0)
    crossIdx = find(thisPanelps==0);
    % check this trial is valid before fetching data
    if ~isempty(crossIdx)
        % check no crosses have been indexed repeatedly by mistake
        crossIdx = crossIdx([1; find(diff(crossIdx)>1)+1]);

        % determine if sweep goes right or left after first crossing
        firstSweepDir = mean(panelps(crossIdx(1):crossIdx(2),1));
        % either way, sweep index for left-right-left sweep
        if firstSweepDir>0 % if right first, start at second crossing
            sweepIdx = crossIdx(2:2:end);
        else % if left first, start at first crossing
            sweepIdx = crossIdx(1:2:end);
        end
        nSweeps = length(sweepIdx)-1; %ignore last (incomplete) sweep

        % initialize
        if c==1
            pxDwellTime = 22/2;
            sweepDur = round(mean(sweepIdx(2:end)-sweepIdx(1:end-1)))+pxDwellTime;
            varBinned = NaN(sweepDur,nTrials*nSweeps);
            runBinned = NaN(sweepDur,nTrials*nSweeps);
            is=0; %sweep index
        end
        % for each sweep, pull and store data
        for ns = 1:nSweeps
            dIdx = ns + is; %data index
            % store panel data
            sweepPos_single(:,ns,c) = thisPanelps(sweepIdx(ns):sweepIdx(ns)+sweepDur-1);
            % store output data
            varBinned(1:sweepDur,dIdx) = thisVarOut(sweepIdx(ns):sweepIdx(ns)+sweepDur-1);
            % store index
            runBinned(1:sweepDur,dIdx) = runIdx(sweepIdx(ns):sweepIdx(ns)+sweepDur-1,c);
        end
        is = is+nSweeps; %update sweep index
        c = c+1; %update counter
    end
end
% store mean sweep position
sweepPos = mean(sweepPos_single,[2 3]);

%% determine which trials met behavioral requirements (if any) for this run
minIdx = round(sweepDur*minInclusion);
goodSweeps = sum(runBinned,1)>minIdx;

% check that enough sweeps can be included for the mean estimate
minGoodSweeps = 3;
if sum(goodSweeps)>=minGoodSweeps
    % calculate and store mean output variable
    sweepMean= mean(varBinned(:,goodSweeps),2);
else
    sweepMean= nan(length(varBinned),1);
end


