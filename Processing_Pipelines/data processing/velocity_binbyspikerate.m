% velocity_binbyspikerate
% This function generates a summary plot of average behavioral velocities 
% (forward, angular, and sideways) binned according to spike rate. 
% Spike rate bins are defined in 10 spk/s increments up to 120 spk/s.
% The function optionally excludes transition periods (start/stop movements)
% and applies lag shifts to align behavior with neural activity.
%
% INPUTS:
% forward      - array of forward velocities (mm/s)
% angular      - array of angular velocities (deg/s)
% sideway      - array of sideways velocities (mm/s)
% cellactivity - spike rate data (spikes/s), time x trial
% ttime        - array of trial times (in seconds)
% lagOpt       - flag (1 = apply lag shifts, 0 = omit lag shifts)
%
% OUTPUTS:
% binOut       - structure containing binned behavioral means:
%                  .spikeRateBin - spike rate bin centers (0:10:120)
%                  .fwdMean      - mean forward velocity per spike rate bin
%                  .angMean      - mean angular velocity per spike rate bin
%                  .sidMean      - mean sideways velocity per spike rate bin
%
% CREATED: 04/02/2025 - MC

function binOut = velocity_binbyspikerate(forward, angular, sideway, cellactivity, ttime, lagOpt)

%% Set Analysis Parameters
settings = processSettings();

% Define spike rate binning parameters
srMax = 100;     % Max spike rate (spikes/s)
srStep = 5;     % Spike rate bin size
spkEdges = -srStep/2:srStep:srMax+srStep/2;
spkBins = 0:srStep:srMax;

% Minimum number of samples per bin
minCount = fetchTimeIdx(ttime, 5);

%% Optional: Exclude Start/Stop Transitions
ex_startstop = 1;
postStartWin = 0.1;
preStopWin = 0.2;
nTrials = size(cellactivity, 2);

if ex_startstop
    postStartIdx = fetchTimeIdx(ttime, postStartWin);
    preStopIdx = fetchTimeIdx(ttime, preStopWin);

    for trial = 1:nTrials
        runIdx = schmittTrigger(forward(:, trial), settings.runThreshE, 0.1);
        runTransitions = diff(runIdx);
        startTrans = find(runTransitions == 1);
        stopTrans = find(runTransitions == -1);

        for st = 1:length(startTrans)
            tStart = startTrans(st);
            tEnd = min(size(cellactivity,1), tStart + postStartIdx);
            cellactivity(tStart:tEnd, trial) = nan;
        end

        for sp = 1:length(stopTrans)
            tStop = stopTrans(sp);
            tStart = max(1, tStop - preStopIdx);
            cellactivity(tStart:tStop, trial) = nan;
        end

        cellactivity(runIdx == 0, trial) = nan;
    end
end

%% Optional: Apply Lag Shifts
if lagOpt
    idxf = fetchTimeIdx(ttime, settings.fwdLag) - 1;
    idxa = fetchTimeIdx(ttime, settings.angLag) - 1;
    idxs = fetchTimeIdx(ttime, settings.sidLag) - 1;

    forward = circshift(forward, -idxf, 1); forward(end-idxf+1:end, :) = nan;
    angular = circshift(angular, -idxa, 1); angular(end-idxa+1:end, :) = nan;
    sideway = circshift(sideway, -idxs, 1); sideway(end-idxs+1:end, :) = nan;
end

%% Reshape Datasets for Bin Processing
cellactivity_r = reshape(cellactivity, [], 1);
forward_r = reshape(forward, [], 1);
angular_r = reshape(angular, [], 1);
sideway_r = reshape(sideway, [], 1);

%% Discretize Spike Rate into Bins
spk_disc = discretize(cellactivity_r, spkEdges, spkBins);

%% Compute Average Behavior per Spike Rate Bin
fwdBySR = arrayfun(@(b) mean(forward_r(spk_disc == b), 'omitnan'), spkBins);
angBySR = arrayfun(@(b) mean(angular_r(spk_disc == b), 'omitnan'), spkBins);
sidBySR = arrayfun(@(b) mean(sideway_r(spk_disc == b), 'omitnan'), spkBins);

%% Threshold Bins with Too Few Samples
spkCounts = arrayfun(@(b) sum(spk_disc == b), spkBins);
fwdBySR(spkCounts < minCount) = NaN;
angBySR(spkCounts < minCount) = NaN;
sidBySR(spkCounts < minCount) = NaN;

%% Store Results in Output Structure
binOut.spikeRateBin = spkBins;
binOut.fwdMean = fwdBySR;
binOut.angMean = angBySR;
binOut.sidMean = sidBySR;

end
