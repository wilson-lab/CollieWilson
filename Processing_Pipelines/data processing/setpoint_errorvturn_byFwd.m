% setpoint_errorvturn_byFwd
% This function calculates the relationship between setpoint error (object position) 
% and the fly's angular velocity, separately for low, medium, and high forward velocity. 
% For time points where forward velocity > 0, it splits the data into thirds (tertiles) 
% and bins angular velocity by setpoint error within each group. 
% Optionally, angular and forward velocity data can be shifted to account for visuomotor lag.
%
% Inputs:
% - panelps: matrix of object positions over time.
% - angular: matrix of angular velocity over time.
% - forward: matrix of forward velocity over time.
% - ttime: vector of time points corresponding to the data.
% - settings: structure containing parameters like visuomotorLag and pursuitGain.
% - optLag: flag indicating whether to shift data based on visuomotor lag.
% - optPlot: flag to optionally generate plots of binned angular velocity vs setpoint error.
%
% Outputs:
% - posvang: 3D matrix of binned angular velocity [position bin, forward group, condition].
% - posvangRL: same as posvang, but left and right bins are combined to account for bias.
% - posBins: vector of position bin centers.
%
% Created: 03/26/2025 MC
%
function [posvang, posvangRL, posBins] = setpoint_errorvturn_byFwd(panelps, angular, forward, ttime, settings, optLag, optPlot)
%% initialize
nCond = size(panelps,3);

% set bin parameters
posMax = 117; % +/- deg
posBin = 9; % deg
posEdge = -posMax-posBin/2:posBin:posMax+posBin/2;
posBins = -posMax:posBin:posMax;
nPosBins = length(posBins);

% optional shift based on visuomotor lag
if optLag
    [idx_vm] = fetchTimeIdx(ttime,settings.visuomotorLag);
    idx_vm = idx_vm - 1;

    angular = circshift(angular,-idx_vm,1);
    angular(end-idx_vm+1:end,:,:) = nan;

    forward = circshift(forward,-idx_vm,1);
    forward(end-idx_vm+1:end,:,:) = nan;
end

% output matrices: [posBin, condition, forward group]
posvang = nan(nPosBins, nCond, 3);
posvangRL = nan(nPosBins, nCond, 3);

minBin = 500;

for c = 1:nCond
    % center position to remove bias
    biasHD = mean(panelps(:,:,c), 'all', 'omitnan');
    thisPanelps = reshape(panelps(:,:,c) - biasHD, [], 1);
    thisAngular = reshape(angular(:,:,c), [], 1);
    thisForward = reshape(forward(:,:,c), [], 1);

    % only keep forward velocity > 0
    fwdIdx = thisForward > 0;
    forwardNonzero = thisForward(fwdIdx);

    if isempty(forwardNonzero)
        continue
    end

    % define tertile cutoffs
    cutoff1 = quantile(forwardNonzero, 1/3);
    cutoff2 = quantile(forwardNonzero, 2/3);

    % categorize into low (1), medium (2), high (3)
    fwdGroup = nan(size(thisForward));
    fwdGroup(thisForward > 0 & thisForward <= cutoff1) = 1;
    fwdGroup(thisForward > cutoff1 & thisForward <= cutoff2) = 2;
    fwdGroup(thisForward > cutoff2) = 3;

    % discretize panel position
    discPanelps = discretize(thisPanelps, posEdge, posBins);

    for f = 1:3
        for p = 1:nPosBins
            thisBin = posBins(p);
            idx = (discPanelps == thisBin) & (fwdGroup == f);

            if sum(idx) >= minBin
                posvang(p,c,f) = mean(thisAngular(idx), 'omitnan');
            end
        end
    end
end

% combine R and L
for f = 1:3
    RL = (posvang(:,:,f) + flip(-posvang(:,:,f),1)) ./ 2;
    for c = 1:nCond
        % clean isolated points
        minSequenceLength = 5;
        raw = posvang(:,c,f);
        rl = RL(:,c);

        rawclean = removeIsolatedPoints(raw, minSequenceLength);
        rlclean = removeIsolatedPoints(rl, minSequenceLength);

        % if 0 bin is missing, mark the whole condition
        zeroidx = find(posBins==0);
        if isnan(rlclean(zeroidx))
            rawclean(:) = nan;
            rlclean(:) = nan;
        end

        posvang(:,c,f) = rawclean;
        posvangRL(:,c,f) = rlclean;
    end
end

%% optional plot
if optPlot
    fwdLabel = {'Low Fwd','Mid Fwd','High Fwd'};
    figure; set(gcf,'Position',[100 100 1500 600])
    tiledlayout(3,nCond,'TileSpacing','compact')

    for c = 1:nCond
        for f = 1:3
            nexttile
            plot(posBins, posvang(:,c,f), 'Color', settings.HDColor)
            xline(0); yline(0); ylim([-100 100])
            title([fwdLabel{f} ', ' num2str(settings.pursuitGain(c)) 'X'])
        end
    end
end

end
