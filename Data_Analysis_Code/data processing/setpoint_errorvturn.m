% setpoint_errorvturn
% This function calculates the relationship between the setpoint error 
% (object position) and the fly's turn behavior (angular velocity). It bins 
% the angular velocity based on the object position (panel position) to estimate 
% the turn magnitude for different setpoint errors. The function also provides an 
% option to shift the angular data according to visuomotor lag estimates.
%
% Inputs:
% - panelps: matrix of object positions (panel position) over time.
% - angular: matrix of angular velocity over time.
% - ttime: vector of time points corresponding to the velocities.
% - settings: structure containing parameters like visuomotorLag and pursuitGain.
% - optLag: flag indicating whether to shift angular velocity based on visuomotor lag.
% - optPlot: flag to optionally generate plots of binned angular velocity vs setpoint error.
%
% Outputs:
% - posvang: matrix of binned angular velocity without accounting for bias.
% - posvangRL: matrix of binned angular velocity, combined to account for bias.
% - posBins: vector of position bin labels.
%
function [posvang,posvangRL,posBins] = setpoint_errorvturn(panelps,angular,ttime,settings,optLag,optPlot)
%% initialize
% fetch number of conditions
nCond = size(panelps,3);

% set bin parameters
posMax = 117; %+/- deg
posBin = 9; %deg

% create bins
posEdge = -posMax-posBin/2:posBin:posMax+posBin/2; %bin edges
posBins = -posMax:posBin:posMax; %bin labels (center)
nPosBins = length(posBins);

%% (optional) shift according to lag estimates
% if lag estimates were provided, shift
if optLag
    % fetch shift indices for each lag
    [idx_vm] = fetchTimeIdx(ttime,settings.visuomotorLag);
    idx_vm = idx_vm-1;
    
    % shift and exclude data at start/stop of trial
    angular = circshift(angular,-idx_vm,1); %shift
    angular(end-idx_vm+1:end,:) = nan; %exclude shifts
end

%% bin turn magnitude according to setpoint error (panel position)

% initialize
posBinnedAngular = nan(nPosBins,nCond);
minBin = 250;
nBin = [];

% for each condition
for c = 1:nCond
    x=1;
    % estimate HD bias
    biasHD = mean(panelps(:,:,c),'all','omitnan');

    % fetch and reshape data
    thisPanelps = reshape(panelps(:,:,c)-biasHD,[],1);
    thisAngular = reshape(angular(:,:,c),[],1);

    % discretize panel data
    binIdx = discretize(thisPanelps, posEdge);

    % calculate mean for each position bin
    for p = 1:nPosBins
        thisBinIdx = (binIdx == p);
        nBin(p,c) = sum(thisBinIdx);
        if nBin(p,c) >= minBin
            posBinnedAngular(p,c) = mean(thisAngular(thisBinIdx), 'omitnan');
        else
            posBinnedAngular(p,c) = NaN;
        end
    end

end

% no combine
raw_posvang = posBinnedAngular;

% combine R and L to account for bias
RL_posvang = (posBinnedAngular + flip(-posBinnedAngular,1))./2;

% clean isolated points
minSequenceLength = 5;
rawclean_posvang = removeIsolatedPoints(raw_posvang, minSequenceLength);
RLclean_posvang = removeIsolatedPoints(RL_posvang, minSequenceLength);

% check to make sure no errors were made
zeroidx = find(posBins==0);
for c = 1:nCond
    if isnan(RLclean_posvang(zeroidx,c))
        RLclean_posvang(:,c) = nan;
        rawclean_posvang(:,c) = nan;
    end
end

% store
posvang = rawclean_posvang;
posvangRL = RLclean_posvang;

%% optional plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1500 600])
    tiledlayout(1,nCond,'TileSpacing','compact')
    % for each condition
    for c = 1:nCond
        nexttile
        plot(posBins,posvang(:,c),'Color',settings.HDColor)
        xline(0); yline(0); ylim([-100 100])
        title([num2str(settings.pursuitGain(c)) 'X'])
    end
end

end
