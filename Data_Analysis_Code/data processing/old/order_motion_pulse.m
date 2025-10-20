%order_motion_pulse
%
% used for pulling out the sweep indices for motion pulse experiments
%
% INPUT
% panelps
% nSpeeds - number of speed conditions present
% thisSpeed - select speed condition to pull
% bufferWindow - indices on either side of sweep to pull
%
% OUTPUT
% orderSweepsR/L
%
% CREATED   05/02/2023 - MC
% UPDATED   08/08/2023 - MC added stationary pulse ordering
%
function  [orderSweepsR,orderSweepsL] = order_motion_pulse(panelps,nSpeeds,thisSpeed,bufferWindow)
%% make sure no outliers were missed
panelps(abs(diff(panelps))>10) = nan;

% find when motion pulse based on where sweeps start/stop
findSweeps = ~isnan(panelps); %find sweeps (1) vs blanks (0)
findStartStop = diff(findSweeps); %find both sweep starts (1) and stops (-1)
idxStartStop=[];
idxStartStop(:,1) = find(findStartStop>0)+1;
idxStartStop(:,2) = find(findStartStop<0);

% for each motion pulse, bin by sweep speed (using length as proxy)
findSweepLengths = round(idxStartStop(:,2) - idxStartStop(:,1), -1);
if any(findSweepLengths<350)
    idxBadSweep = find(findSweepLengths<350);
    idxStartStop(idxBadSweep,:) = [];
    findSweepLengths(idxBadSweep) = [];
end
[binnedSpeed,~] = discretize(findSweepLengths,nSpeeds);
for s = 1:nSpeeds
    lengthBins(s) = max(findSweepLengths(binnedSpeed==s));
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
        thisIdxR(a) = intersect(intersect(find(binnedStart==a),find(binnedDirection==1)),find(binnedSpeed==s)); %rightward
        thisIdxL(a) = intersect(intersect(find(binnedStart==a),find(binnedDirection==2)),find(binnedSpeed==s)); %leftward

        % pull sweep
        win = bufferWindow; %add buffer window before/after sweep

        orderSweepsR(:,a) = idxStartStop(thisIdxR(a))-win:idxStartStop(thisIdxR(a))+win+lengthBins(s)-1;
        orderSweepsL(:,a) = idxStartStop(thisIdxL(a))-win:idxStartStop(thisIdxL(a))+win+lengthBins(s)-1;
    end

    % ensure buffer didnt fall outside of available indices
    orderSweepsR(find(orderSweepsR<1)) = 1;
    orderSweepsL(find(orderSweepsL<1)) = 1;
else
    % for each motion pulse, bin by start location
    nStart = length(binnedSpeed); %number of possible start positions
    [binnedStart,~] = discretize(panelps(idxStartStop(:,1)),nStart);
    % put it all together now...
    s = thisSpeed;
    for a = 1:nStart
        % pull start positions (from left to right) going right (1) or
        % left (2) for each sweep speed (from slow to fast)
        thisIdxR(a) = intersect(find(binnedStart==a),find(binnedSpeed==s)); %rightward
        thisIdxL(a) = intersect(find(binnedStart==a),find(binnedSpeed==s)); %leftward

        % pull sweep
        win = bufferWindow; %add buffer window before/after sweep

        orderSweepsR(:,a) = idxStartStop(thisIdxR(a))-win:idxStartStop(thisIdxR(a))+win+lengthBins(s)-1;
        orderSweepsL(:,a) = idxStartStop(thisIdxL(a))-win:idxStartStop(thisIdxL(a))+win+lengthBins(s)-1;
    end

    % ensure buffer didnt fall outside of available indices
    orderSweepsR(find(orderSweepsR<1)) = 1;
    orderSweepsL(find(orderSweepsL<1)) = 1;
end
end