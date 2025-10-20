% spikerate_v_heading
%
% analysis function for binning the average spikerate response across
% heading positions (aka x position of the panels)
%
% INPUT
% allPanelps - panel positions
% allSpikert - spike rate
%
% OUTPUT
% headingData - summary
%
% CREATED 11/16/2023 - MC
%

function [headingData] = spikerate_v_heading(allPanelps,allSpikert,optPlot)
%% settings

% set the bin size for panel positions
arenaSize = 360; %deg
nPixels = 192; %px
stepSize = arenaSize/nPixels;
binEdges = (0:stepSize:arenaSize)-180;
binValue = (stepSize/2:stepSize:360 - stepSize/2)-180;

% reshape datasets
rshPanelps = reshape(allPanelps,[],1)-180;
rshSpikert = reshape(allSpikert,[],1);

% plot details
lw = 2;

%% discretize heading
% heading position
discPanelps = discretize(rshPanelps,binEdges,binValue);
% store
headingData.bins = discPanelps;

% calculate average spikerate for each bin
for pbin = 1:length(binValue)
    thisBin = binValue(pbin);
    thisIdx = find(discPanelps==thisBin);
    
    meanSpikert(pbin) = mean(rshSpikert(thisIdx),'omitnan');
    stdSpikert(pbin) = std(rshSpikert(thisIdx),'omitnan');
end

% store
headingData.spikert = meanSpikert;
headingData.std = stdSpikert;

%% plot
if optPlot
    figure; set(gcf,'Position',[100 100 1500 500])
    tiledlayout(1,2,'TileSpacing','compact')

    % plot heading distribution
    nexttile
    histogram(rshPanelps,'BinEdges',-180:10:180,'FaceColor','#77AC30')
    ylabel('Distribution (counts)')
    xlabel('Target Position (deg)')
    xlim([-180 180])
    xline(0)

    % plot spikerate relationship
    nexttile
    % plot SEM band
    s = patch([binValue'; flipud(binValue')],[(meanSpikert-stdSpikert)'; flipud((meanSpikert+stdSpikert)')], 'r', 'FaceAlpha',0.1, 'EdgeColor','none');
    s.FaceColor = '#77AC30';
    % plot mean
    hold on
    plot(binValue,meanSpikert,'Color','#77AC30','LineWidth',lw)
    ylabel('Average Spike Rate')
    xlabel('Target Position (deg)')
    xlim([-180 180])
    xline(0)
end

end