% pulse_v_output_snr
%
% Calculates the mean, variance, and signal-to-noise ratio (SNR) of firing rate
% responses during motion pulses across object sweep positions.
%
% For each sweep position, the function:
%   1) Computes the mean firing rate during the middle third of the pulse for
%      each trial that meets behavioral inclusion criteria.
%   2) Calculates the variance across trials of these per-trial mean responses.
%   3) Calculates the mean across trials of these per-trial mean responses.
%   4) Computes SNR as variance divided by mean for each sweep position.
%
% INPUTS:
%   panelps       - Panel position trace (time × trials)
%   forward       - Forward velocity trace (time × trials)
%   varOut        - Output variable (e.g., firing rate; time × trials)
%   ttime         - Time vector corresponding to rows of input data
%   pulseSelect   - Pulse type selector passed to order_motion_pulse
%   pulseOptions  - Available pulse types
%   nSweep        - Number of sweep positions
%   runSelect     - 0 = quiescent only, >0 = running only, <0 = all data
%
% OUTPUTS:
%   posVar        - Variance across trials of mean firing rates per position
%   posMean       - Mean across trials of mean firing rates per position
%   posSNR        - Signal-to-noise ratio (variance / mean) per position
%
% CREATED: 10/05/2025 - MC

function [posVar, posMean, posSNR] = pulse_v_output_snr(panelps,forward,varOut,ttime,pulseSelect,pulseOptions,nSweep,runSelect)

%% initialize
nTrials  = size(panelps,2);
nPulses  = length(pulseOptions);

% set analysis window
windowSize = 2; % size of full window (sec)
winIdx = fetchTimeIdx(ttime, windowSize);

% set number of pulse trials required
minTrials = 3;

%% fetch behavioral index depending on run select
if runSelect==0 % quiescent only
    moveThresh = 0.5;                    % anything above/below is moving
    runIdx = ~(forward<-moveThresh | forward>moveThresh);
    minInclusion = 0.5;                  % fraction of window that must be quiescent
elseif runSelect>0 % running only
    runIdx = schmittTrigger(forward, runSelect, 0.1);
    minInclusion = 0.75;
else % all
    runIdx = ones(size(forward));
    minInclusion = 0.75;
end

%% pull motion pulses from pseudorandomized dataset
c = 1; % counter over successful trials
for t = 1:nTrials
    thisPanelps = panelps(:,t);
    thisVarOut  = varOut(:,t);
    try
        [thisOrderR, thisOrderL] = order_motion_pulse(thisPanelps, nPulses, pulseSelect, winIdx);
        thisOrderR(thisOrderR>size(thisPanelps,1)) = size(thisPanelps,1);
        thisOrderL(thisOrderL>size(thisPanelps,1)) = size(thisPanelps,1);

        if size(thisOrderR,2)==nSweep
            for p = 1:nSweep
                pulsePanelpsR(:,c,p) = thisPanelps(thisOrderR(:,p));
                pulsePanelpsL(:,c,p) = thisPanelps(thisOrderL(:,p));
                pulseVarR(:,c,p)     = thisVarOut(thisOrderR(:,p));  % firing rate
                pulseVarL(:,c,p)     = thisVarOut(thisOrderL(:,p));
                runIdxR(:,c,p)       = runIdx(thisOrderR(:,p), t);
                runIdxL(:,c,p)       = runIdx(thisOrderL(:,p), t);
            end
            c = c+1;
        end
    catch
        nWin = numel(winIdx);
        for p = 1:nSweep
            pulsePanelpsR(1:nWin,c,p) = nan(nWin,1);
            pulsePanelpsL(1:nWin,c,p) = nan(nWin,1);
            pulseVarR(1:nWin,c,p)     = nan(nWin,1);
            pulseVarL(1:nWin,c,p)     = nan(nWin,1);
            runIdxR(1:nWin,c,p)       = false(nWin,1);
            runIdxL(1:nWin,c,p)       = false(nWin,1);
        end
        c = c+1;
    end
end

% pulse duration (# time points in window)
pulseDur = size(pulsePanelpsR,1);
% calculate mean panel position and store
panelpsR = mean(pulsePanelpsR,2,"omitnan");
panelpsR(panelpsR==0) = nan;
% find when the object is actually presented
sweepIdx = find(~isnan(panelpsR(:,1,1)));

%% determine trials meeting behavior criterion
minIdx      = round(pulseDur*minInclusion);
goodTrialsR = sum(runIdxR,1) >= minIdx;    % [1 x trials x sweep]

if ~isempty(sweepIdx)
    %% middle third indices (analysis window)
    t1     = sweepIdx(1);
    t2     = sweepIdx(end);
    midIdx = t1:t2;

    %% outputs across object positions (rightward sweeps)
    posVar  = nan(1, nSweep);   % variance across trials of per-trial means
    posMean = nan(1, nSweep);   % mean     across trials of per-trial means
    posSNR  = nan(1, nSweep);   % SNR = posVar ./ posMean

    for p = 1:nSweep
        idxR = squeeze(goodTrialsR(1,:,p));                  % logical [trials]
        if ~isempty(idxR) && sum(idxR) >= minTrials
            % Data for this position & direction: [time_mid x trials]
            dat_mid = pulseVarR(midIdx, idxR, p);

            % Per-trial mean firing rate in the window (1 x trials)
            trialMeans = mean(dat_mid, 1, 'omitnan');

            % Across-trial summary
            mu  = mean(trialMeans, 'omitnan');
            sg2 = var( trialMeans, 0, 'omitnan');            % sample variance across trials

            posMean(p) = mu;
            posVar(p)  = sg2;
            posSNR(p)  = mu/sg2;
        end
    end
else
    posMean(1:nSweep) = nan;
    posVar(1:nSweep)  = nan;
    posSNR(1:nSweep)  = nan;
end
end
