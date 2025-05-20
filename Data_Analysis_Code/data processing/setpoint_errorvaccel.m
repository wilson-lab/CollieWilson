% setpoint_errorvaccel
% This function calculates the relationship between setpoint error 
% (object position) and the fly's angular acceleration behavior. 
% It first converts angular velocity to angular acceleration using the 
% velocity2acceleration function, then bins the angular acceleration 
% based on the object position to estimate the average angular 
% acceleration for different setpoint errors.
%
% Inputs:
% - panelps: matrix of object positions (panel position) over time.
% - angular: matrix of angular velocity over time.
% - ttime: vector of time points corresponding to the velocities.
% - settings: structure containing parameters like visuomotorLag and pursuitGain.
% - optLag: flag indicating whether to shift angular velocity based on visuomotor lag.
% - optPlot: flag to optionally generate plots of binned angular acceleration vs setpoint error.
%
% Outputs:
% - accelPosvang: matrix of binned angular acceleration.
% - accelPosvangRL: matrix of binned angular acceleration, combined to account for bias.
% - posBins: vector of position bin labels.

function [accelPosvang, accelPosvangRL, posBins] = setpoint_errorvaccel(panelps,angular,ttime,settings,optLag,optPlot)
%% initialize
% fetch number of conditions
nCond = size(panelps, 3);

% set bin parameters
posMax = 117; %+/- deg
posBin = 9; %deg

% create bins
posEdge = -posMax - posBin/2 : posBin : posMax + posBin/2; % bin edges
posBins = -posMax : posBin : posMax; % bin labels (center)
nPosBins = length(posBins);

%% (optional) shift according to lag estimates
if optLag
    % fetch shift indices for each lag
    [idx_vm] = fetchTimeIdx(ttime, settings.visuomotorLag);
    idx_vm = idx_vm - 1;
    
    % shift and exclude data at start/stop of trial
    angular = circshift(angular, -idx_vm, 1); % shift
    angular(end-idx_vm+1:end, :) = nan; % exclude shifts
end

%% calculate angular acceleration using velocity2acceleration
accel = velocity2acceleration([], angular, [], ttime); % only use angular velocity

%% bin angular acceleration according to setpoint error (panel position)
% initialize
accelBinned = nan(nPosBins, nCond);
minBin = 500;
nBin = [];

% for each condition
for c = 1:nCond
    x = 1;
    % estimate HD bias
    biasHD = mean(panelps(:,:,c), 'all', 'omitnan');

    % fetch and reshape data
    thisPanelps = reshape(panelps(:,:,c) - biasHD, [], 1);
    thisAccel = reshape(accel.angular(:,:,c), [], 1);

    % discretize panel data
    discPanelps = discretize(thisPanelps, posEdge, posBins);

    % calculate mean for each position bin
    for p = 1:nPosBins
        thisBin = posBins(p);
        thisBinIdx = find(discPanelps == thisBin);
        nBin(x, c) = length(thisBinIdx);
        x = x + 1;
        
        % only include the bin if it has at least 3 points
        if length(thisBinIdx) >= minBin
            accelBinned(p, c) = mean(thisAccel(thisBinIdx), 'omitnan');
        else
            accelBinned(p, c) = NaN; % exclude bin with fewer than 3 points
        end
    end
end

% no combine
raw_accelPosvang = accelBinned;

% combine R and L to account for bias
RL_accelPosvang = (accelBinned + flip(-accelBinned, 1)) / 2;

% clean isolated points
minSequenceLength = 5;
rawclean_accelPosvang = removeIsolatedPoints(raw_accelPosvang, minSequenceLength);
RLclean_accelPosvang = removeIsolatedPoints(RL_accelPosvang, minSequenceLength);

% store
accelPosvang = rawclean_accelPosvang;
accelPosvangRL = RLclean_accelPosvang;

%% optional plot
if optPlot
    % initialize
    figure; set(gcf, 'Position', [100 100 1500 600])
    tiledlayout(1, nCond, 'TileSpacing', 'compact')
    
    % for each condition
    for c = 1:nCond
        nexttile
        plot(posBins, accelPosvang(:, c), 'Color', settings.HDColor)
        xline(0); yline(0); ylim([-100 100])
        title([num2str(settings.pursuitGain(c)) 'X'])
    end
end

end
