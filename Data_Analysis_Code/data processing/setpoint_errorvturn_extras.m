% setpoint_errorvturn_extras
% This function determines the relationship between setpoint error (object position) 
% and the fly's turn behavior. It bins the object position and calculates the mean turn 
% velocity (angular) for different behavioral categorizations, such as object velocity 
% (fast/slow) and fly's forward velocity (fast/slow).
%
% INPUTS:
% panelps   - panel positions (degrees) representing object position
% panelvel  - panel velocity (degrees/second) representing object motion
% forward   - forward velocities of the fly
% angular   - angular velocities of the fly
% ttime     - trial time (in seconds)
% settings  - structure containing analysis settings
% optLag    - 1 to apply lag shift based on visuomotor delay, 0 to omit
% optPlot   - 1 to generate plots, 0 to skip plotting
%
% OUTPUTS:
% posvang   - binned turn velocity data for different behavioral conditions (fast/slow)
% posBins   - position bins (degrees)
%
% CREATED: [Date] MC
%
function [posvang,posBins] = setpoint_errorvturn_extras(panelps,panelvel,forward,angular,ttime,settings,optLag,optPlot)
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

% fetch only forward forward velocities
forward(forward<0) = nan;

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

%% add extra categorizations

% break up when object was moving right vs left
rightvel = panelvel>0;
leftvel = panelvel<0;

panelps_right = panelps;
panelps_right(~rightvel) = nan;
panelps_left = panelps;
panelps_left(~leftvel) = nan;

% break up when object was moving fast vs slow
velocityThreshold = 100;

fastobj = abs(panelvel)>velocityThreshold;
slowobj = abs(panelvel)<velocityThreshold;

panelps_objfast = panelps;
panelps_objfast(~fastobj) = nan;
panelps_objslow = panelps;
panelps_objslow(~slowobj) = nan;

% break up when fly was running fast v slow
forwardThreshold = prctile(forward(:),80);

fastforward = forward>forwardThreshold;
slowforward = forward<forwardThreshold;

panelps_fwdfast = panelps;
panelps_fwdfast(~fastforward) = nan;
panelps_fwdslow = panelps;
panelps_fwdslow(~slowforward) = nan;


%% bin turn magnitude according to setpoint error (panel position)

% initialize
posvang = nan(nPosBins,nCond,4);

binRightward = nan(nPosBins,nCond);
binLeftward = nan(nPosBins,nCond);

binObjFast = nan(nPosBins,nCond);
binObjSlow = nan(nPosBins,nCond);

binFwdFast = nan(nPosBins,nCond);
binFwdSlow = nan(nPosBins,nCond);

% for each condition
for c = 1:nCond
    % estimate HD bias
    biasHD = mean(panelps(:,:,c),'all','omitnan');

    % fetch and reshape data
    thisPosRight = reshape(panelps_right(:,:,c)-biasHD,[],1);
    thisPosLeft = reshape(panelps_left(:,:,c)-biasHD,[],1);

    thisPosFast = reshape(panelps_objfast(:,:,c)-biasHD,[],1);
    thisPosSlow = reshape(panelps_fwdslow(:,:,c)-biasHD,[],1);

    thisFwdFast = reshape(panelps_fwdfast(:,:,c)-biasHD,[],1);
    thisFwdSlow = reshape(panelps_objslow(:,:,c)-biasHD,[],1);

    thisAngular = reshape(angular(:,:,c),[],1);

    % discretize panel data
    discPosRight = discretize(thisPosRight,posEdge,posBins);
    discPosLeft = discretize(thisPosLeft,posEdge,posBins);

    discPosFast = discretize(thisPosFast,posEdge,posBins);
    discPosSlow = discretize(thisPosSlow,posEdge,posBins);

    discFwdFast = discretize(thisFwdFast,posEdge,posBins);
    discFwdSlow = discretize(thisFwdSlow,posEdge,posBins);

    % calculate mean for each bin
    for p = 1:nPosBins
        thisBin = posBins(p);

        % Rightward object movement
        rIdx = find(discPosRight==thisBin);
        if length(rIdx) >= 3
            binRightward(p,c) = mean(thisAngular(rIdx),'omitnan');
        else
            binRightward(p,c) = NaN; % exclude bin with fewer than 3 points
        end

        % Leftward object movement
        lIdx = find(discPosLeft==thisBin);
        if length(lIdx) >= 3
            binLeftward(p,c) = mean(thisAngular(lIdx),'omitnan');
        else
            binLeftward(p,c) = NaN; % exclude bin with fewer than 3 points
        end

        % Fast-moving object
        fIdx = find(discPosFast==thisBin);
        if length(fIdx) >= 3
            binObjFast(p,c) = mean(thisAngular(fIdx),'omitnan');
        else
            binObjFast(p,c) = NaN; % exclude bin with fewer than 3 points
        end

        % Slow-moving object
        sIdx = find(discPosSlow==thisBin);
        if length(sIdx) >= 3
            binObjSlow(p,c) = mean(thisAngular(sIdx),'omitnan');
        else
            binObjSlow(p,c) = NaN; % exclude bin with fewer than 3 points
        end

        % Fast-forward velocity
        fIdx = find(discFwdFast==thisBin);
        if length(fIdx) >= 3
            binFwdFast(p,c) = mean(thisAngular(fIdx),'omitnan');
        else
            binFwdFast(p,c) = NaN; % exclude bin with fewer than 3 points
        end

        % Slow-forward velocity
        sIdx = find(discFwdSlow==thisBin);
        if length(sIdx) >= 3
            binFwdSlow(p,c) = mean(thisAngular(sIdx),'omitnan');
        else
            binFwdSlow(p,c) = NaN; % exclude bin with fewer than 3 points
        end
    end

end

% combine R and L to account for bias
RL_fast = (binFwdFast + flip(-binFwdFast,1))./2;
RL_slow = (binFwdSlow + flip(-binFwdSlow,1))./2;
% clean isolated points
minSequenceLength = 5;
RLclean_fast = removeIsolatedPoints(RL_fast, minSequenceLength);
RLclean_slow = removeIsolatedPoints(RL_slow, minSequenceLength);
% check to make sure no errors were made
zeroidx = find(posBins==0);
for c = 1:nCond
    if isnan(RLclean_fast(zeroidx,c))
        RLclean_fast(:,c) = nan;
        RLclean_slow(:,c) = nan;
    end
end

% store
posvang(:,:,1) = RLclean_fast;
posvang(:,:,2) = RLclean_slow;

%posvang(:,:,1) = binRightward;
%posvang(:,:,2) = binLeftward;

% optionally combine R and L to account for bias
%posvang(:,:,1) = (binRightward + flip(-binRightward,1))./2;
%posvang(:,:,2) = (binLeftward + flip(-binLeftward,1))./2;

%posvang(:,:,3) = (binObjFast + flip(-binObjFast,1))./2;
%posvang(:,:,4) = (binObjSlow + flip(-binObjSlow,1))./2;

%% optional plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1500 600])
    tiledlayout(2,nCond,'TileSpacing','compact')
    % for each condition
    for c = 1:nCond
        nexttile
        plot(posBins,posvang(:,c),'Color',settings.HDColor)
        xline(0); yline(0); ylim([-100 100])
        title([num2str(settings.pursuitGain(c)) 'X'])
    end
    for c = 1:nCond
        nexttile
        plot(velBins,velvang(:,c),'Color',settings.HDColor)
        xline(0); yline(0); ylim([-100 100])
        title([num2str(settings.pursuitGain(c)) 'X'])
    end
end

end