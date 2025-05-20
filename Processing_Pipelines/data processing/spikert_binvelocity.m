% spikert_binvelocity
% This function generates a summary plot of cell activity data (firing rate 
% or membrane potential) binned according to directional velocities: 
% forward, angular, and sideways. The cell activity is averaged for each 
% velocity bin, optionally excluding transition periods (start/stop 
% movements) and applying lag shifts.
%
% INPUTS:
% forward      - array of forward velocities
% angular      - array of angular velocities
% sideway      - array of sideways velocities
% cellactivity - cell activity data, can be firing rate or voltage data (vm data)
% ttime        - array of trial times (in seconds)
% lagOpt       - flag (1 = apply lag shifts, 0 = omit lag shifts)
%
% OUTPUTS:
% binOut       - structure containing binned means for each directional velocity:
%                  .fwdBin  - forward velocity bin centers
%                  .angBin  - angular velocity bin centers
%                  .sidBin  - sideways velocity bin centers
%                  .fwdMean - mean cell activity in each forward velocity bin
%                  .angMean - mean cell activity in each angular velocity bin
%                  .sidMean - mean cell activity in each sideways velocity bin
%
% ORIGINAL: 12/12/2022 - MC
% UPDATED:  02/07/2023 - MC added output variables
%           06/07/2023 - MC updated to focus on pursuit distributions
%           08/07/2024 - MC added start/stop exclusion and lag shift
%           11/09/2024 - MC changed start/stop exclusion
%           03/23/2025 - MC added bin inclusion requirement
%
function binOut = spikert_binvelocity(forward, angular, sideway, cellactivity, ttime, lagOpt)

%% Set Analysis Parameters
% Retrieve settings and trial duration
settings = processSettings();

% Define binning parameters
fwdMax = 10;   % Forward max velocity (mm/s)
angMax = 260;  % Angular max velocity (deg/s)
sidMax = 3.5;  % Sideways max velocity (mm/s)
fs = 1;        % Forward bin size
as = 20;       % Angular bin size
ss = 0.25;     % Sideways bin size

% Set threshold for minimum number of samples per bin
minCount = fetchTimeIdx(ttime, 5);

%% Optional: Exclude Start/Stop Transitions
% Set flag to exclude start/stop transitions using transition window settings
ex_startstop = 1;
postStartWin = 0.1; % Time window after start (in seconds)
preStopWin = 0.2;   % Time window before stop (in seconds)

% Number of trials (columns in the cell activity array)
nTrials = size(cellactivity, 2);

if ex_startstop
    % Convert post-start and pre-stop windows to indices based on time array
    postStartIdx = fetchTimeIdx(ttime, postStartWin);
    preStopIdx = fetchTimeIdx(ttime, preStopWin);

    % Loop over each trial
    for trial = 1:nTrials
        % Calculate run index using Schmitt Trigger
        runIdx = schmittTrigger(forward(:, trial), settings.runThreshE, 0.1);

        % Identify start and stop transitions in runIdx for the current trial
        runTransitions = diff(runIdx);    % Calculate transitions in run state
        startTrans = find(runTransitions == 1); % 0 to 1 (start running)
        stopTrans = find(runTransitions == -1); % 1 to 0 (stop running)

        % Loop over each start transition to set post-start period as NaN
        for st = 1:length(startTrans)
            tStart = startTrans(st); % Start index
            tEnd = min(size(cellactivity, 1), tStart + postStartIdx); % End index, within bounds
            cellactivity(tStart:tEnd, trial) = nan; % Set post-start window to NaN in cell activity data
        end

        % Loop over each stop transition to set pre-stop period as NaN
        for sp = 1:length(stopTrans)
            tStop = stopTrans(sp); % Stop index
            tStart = max(1, tStop - preStopIdx); % Start index, within bounds
            cellactivity(tStart:tStop, trial) = nan; % Set pre-stop window to NaN in cell activity data
        end

        % Set cellactivity to NaN where runIdx is 0 (not running)
        cellactivity(runIdx == 0, trial) = nan;
    end
end

%% Optional: Apply Lag Shifts
% Shift directional velocity data if lag option is enabled
if lagOpt
    idxf = fetchTimeIdx(ttime, settings.fwdLag) - 1; % Lag shift for forward
    idxa = fetchTimeIdx(ttime, settings.angLag) - 1; % Lag shift for angular
    idxs = fetchTimeIdx(ttime, settings.sidLag) - 1; % Lag shift for sideways
    
    % Apply shifts and mark shifted portions as NaN
    forward = circshift(forward, -idxf, 1); forward(end-idxf+1:end, :) = nan;
    angular = circshift(angular, -idxa, 1); angular(end-idxa+1:end, :) = nan;
    sideway = circshift(sideway, -idxs, 1); sideway(end-idxs+1:end, :) = nan;
end

%% Reshape Datasets for Bin Processing
% Flatten all data arrays to ensure compatibility with binning functions
cellactivity_r = reshape(cellactivity, [], 1);
forward_r = reshape(forward, [], 1);
angular_r = reshape(angular, [], 1);
sideway_r = reshape(sideway, [], 1);

%% Discretize Velocity Data into Bins
% Create bins for each directional velocity and assign data to bins

% Forward velocity binning
f_edge = -fs/2:fs:fwdMax+fs/2; % Bin edges for forward
f_bins = 0:fs:fwdMax; % Bin centers
forward_disc = discretize(forward_r, f_edge, f_bins);

% Angular velocity binning
a_edge = -angMax-as/2:as:angMax+as/2; % Bin edges for angular
a_bins = -angMax:as:angMax; % Bin centers
angular_disc = discretize(angular_r, a_edge, a_bins);

% Sideways velocity binning
s_edge = -sidMax-ss/2:ss:sidMax+ss/2; % Bin edges for sideways
s_bins = -sidMax:ss:sidMax; % Bin centers
sideway_disc = discretize(sideway_r, s_edge, s_bins);

%% Compute Mean Cell Activity for Each Bin
% Calculate average cell activity for each directional velocity bin

fwdAllMean = arrayfun(@(b) mean(cellactivity_r(forward_disc == b), 'omitnan'), f_bins);
angAllMean = arrayfun(@(b) mean(cellactivity_r(angular_disc == b), 'omitnan'), a_bins);
sidAllMean = arrayfun(@(b) mean(cellactivity_r(sideway_disc == b), 'omitnan'), s_bins);

%% Threshold Bins
% Count number of samples in each bin
fwdCounts = arrayfun(@(b) sum(forward_disc == b), f_bins);
angCounts = arrayfun(@(b) sum(angular_disc == b), a_bins);
sidCounts = arrayfun(@(b) sum(sideway_disc == b), s_bins);

% Set bins with too few samples to NaN
fwdAllMean(fwdCounts < minCount) = NaN;
angAllMean(angCounts < minCount) = NaN;
sidAllMean(sidCounts < minCount) = NaN;

%% Store Results in Output Structure
% Compile bin centers and mean cell activity into output structure
binOut.fwdBin = f_bins;
binOut.angBin = a_bins;
binOut.sidBin = s_bins;
binOut.fwdMean = fwdAllMean;
binOut.angMean = angAllMean;
binOut.sidMean = sidAllMean;

%% NEW: Angular velocity relationship split by forward velocity median
% Compute forward velocity median (ignoring NaNs)
fwd_median = median(forward_r(forward_r>=0), 'omitnan');

% Create logical indices for fast and slow forward velocity
isFast = forward_r > fwd_median;
isSlow = forward_r <= fwd_median;

% Compute angular bin means for fast vs slow forward velocity
angMean_fastFwd = arrayfun(@(b) ...
    mean(cellactivity_r(angular_disc == b & isFast), 'omitnan'), a_bins);
angMean_slowFwd = arrayfun(@(b) ...
    mean(cellactivity_r(angular_disc == b & isSlow), 'omitnan'), a_bins);

% Apply minimum sample count threshold
angCounts_fast = arrayfun(@(b) ...
    sum(angular_disc == b & isFast), a_bins);
angCounts_slow = arrayfun(@(b) ...
    sum(angular_disc == b & isSlow), a_bins);

angMean_fastFwd(angCounts_fast < minCount) = NaN;
angMean_slowFwd(angCounts_slow < minCount) = NaN;

% Store in output
binOut.angMean_fastFwd = angMean_fastFwd;
binOut.angMean_slowFwd = angMean_slowFwd;
end
