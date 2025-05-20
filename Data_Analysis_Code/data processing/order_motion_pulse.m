% order_motion_pulse
% This function extracts the sweep indices for motion pulse experiments based on panel positions.
% It identifies the start and end of motion pulses, categorizing them by speed and direction,
% and adjusting for specified analysis windows. It handles stationary and moving pulses separately.
%
% INPUTS:
%   panelps      - Panel positions (degrees)
%   nSpeeds      - Number of speed conditions present
%   thisSpeed    - Select speed condition to pull (e.g., 1 or 2)
%   windowIdx    - Size of the analysis window (in indices)
%
% OUTPUTS:
%   orderSweepsR - Indices of rightward motion sweeps
%   orderSweepsL - Indices of leftward motion sweeps
%
% CREATED: 05/02/2023 - MC
% UPDATED: 08/08/2023 - MC (added stationary pulse ordering)
% UPDATED: 07/30/2024 - MC (added catch for partial sweeps if needed)
%
function  [orderSweepsR,orderSweepsL] = order_motion_pulse(panelps,nSpeeds,thisSpeed,windowIdx)
%% make sure no outliers were missed
panelps(abs(diff(panelps))>10) = nan;

% find when motion pulse based on where sweeps start/stop
findSweeps = ~isnan(panelps); %find sweeps (1) vs blanks (0)
findStartStop = diff(findSweeps); %find both sweep starts (1) and stops (-1)
findStarts = find(findStartStop>0)+1;
findStops = find(findStartStop<0);
% ensure there are no incomplete sweeps
if length(findStarts)>length(findStops)
    findStarts = findStarts(1:length(findStops));
end
% store
idxStartStop=[];
idxStartStop(:,1) = findStarts;
idxStartStop(:,2) = findStops;

% for each motion pulse, bin by sweep speed (using length as proxy)
findSweepLengths = round(idxStartStop(:,2) - idxStartStop(:,1), -1);
checkBadSweeps = isoutlier(findSweepLengths);
if nSpeeds>1
    if any(checkBadSweeps)
        idxStartStop(checkBadSweeps,:) = [];
        findSweepLengths(checkBadSweeps) = [];
    end
end
[binnedSpeed,~] = discretize(findSweepLengths,nSpeeds);
for s = 1:nSpeeds
    lengthBins(s) = median(findSweepLengths(binnedSpeed==s));
end


% set sweep window buffers
win = ceil((windowIdx - lengthBins)./2);
% check if window is even split
if win==floor(win)
    preWin = win;
    pstWin = win;
else
    preWin = floor(win);
    pstWin = ceil(win);
end


%% determine if these pulses are moving or stationary
checkMotion = 2 < abs(panelps(idxStartStop(1,1))-panelps(idxStartStop(1,2))); %is sweep smaller than 1 px
if checkMotion
    % for each motion pulse, bin by sweep direction
    % 1 is leftward, 2 is rightward
    binnedDirection = double(panelps(idxStartStop(:,1))>panelps(idxStartStop(:,2))) + 1;

    % for each motion pulse, bin by start location
    nStart = round(length(binnedDirection)/(nSpeeds*2)); %number of possible start positions
    [binnedStart,~] = discretize(panelps(idxStartStop(:,1)),nStart);

    % put it all together now...
    s = thisSpeed;
    for a = 1:nStart
        % pull start positions (from left to right) going right (1) or
        % left (2) for each sweep speed (from slow to fast)
        try
            thisIdxR(a) = intersect(intersect(find(binnedStart==a),find(binnedDirection==1)),find(binnedSpeed==s)); %rightward
            thisIdxL(a) = intersect(intersect(find(binnedStart==a),find(binnedDirection==2)),find(binnedSpeed==s)); %leftward

            % pull sweep
            orderSweepsR(:,a) = idxStartStop(thisIdxR(a))-preWin(s):idxStartStop(thisIdxR(a))+pstWin(s)+lengthBins(s)-1;
            orderSweepsL(:,a) = idxStartStop(thisIdxL(a))-preWin(s):idxStartStop(thisIdxL(a))+pstWin(s)+lengthBins(s)-1;
        end
    end

    % ensure buffer didnt fall outside of available indices
    orderSweepsR(orderSweepsR<1) = 1;
    orderSweepsL(orderSweepsL<1) = 1;
else
    % for each motion pulse, bin by start location
    nStart = length(binnedSpeed); %number of possible start positions
    [binnedStart,~] = discretize(panelps(idxStartStop(:,1)),nStart);
    % put it all together now...
    s = thisSpeed;
    for a = 1:nStart
        % pull start positions (from left to right)
        try
            thisIdxR(a) = intersect(find(binnedStart==a),find(binnedSpeed==s)); %rightward
            thisIdxL(a) = intersect(find(binnedStart==a),find(binnedSpeed==s)); %leftward

            orderSweepsR(:,a) = idxStartStop(thisIdxR(a))-preWin:idxStartStop(thisIdxR(a))+pstWin+lengthBins(s)-1;
            orderSweepsL(:,a) = idxStartStop(thisIdxL(a))-preWin:idxStartStop(thisIdxL(a))+pstWin+lengthBins(s)-1;
        end
    end

    % ensure buffer didnt fall outside of available indices
    orderSweepsR(orderSweepsR<1) = 1;
    orderSweepsL(orderSweepsL<1) = 1;
end
end