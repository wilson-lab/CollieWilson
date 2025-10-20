% spikerate_v_sweep_psth
% analysis function for generating a summary plot of spike raster binned by
% object position and sweep direction
%
% OUTPUTS:
% sweep_time - sweep time for both L and R
% spike_psth - spike freqency as a function of sweep time, 1 L, 2 R
%
% INPUTS:
% panelpos - panel positions
% spikerst - spike raster
% time - trial time
% optPlot - whether data is plotted
%
% ORIGINAL: 01/06/2023 - MC
%

function [sweep_time,spike_psth] = spikerate_v_sweep_psth(panelpos,spikerst,time,optPlot)
%% initialize
maxPx = max(panelpos(:,2))-2;
minPx = min(panelpos(:,2))+2;

nTrial = size(panelpos,2);

% pxOptions = 1.875:1.875:360;
% panelpos_adj=panelpos;
% for px = 1:length(pxOptions)
%     pxSearch = abs(round(panelpos-pxOptions(px),1));
%     panelpos_adj(pxSearch<=0.4) = pxOptions(px);
% end


%% pull spikerate for each sweep

% initialize
leftwardSpikeRst = [];
leftwardPanelPos = [];
rightwardSpikeRst = [];
rightwardPanelPos = [];

% for each trial

for nt = 1:nTrial
    % find the high peak (start of leftward) and low peak (start of rightward)
    [~,pIdxH] = findpeaks(panelpos(:,nt),'MinPeakHeight',maxPx,'MinPeakDistance',20e3);
    [~,pIdxL] = findpeaks(-panelpos(:,nt),'MinPeakHeight',-minPx,'MinPeakDistance',20e3);
    %plot(panelpos(:,nt),'k'); hold on; xline(pIdxH,'r'); xline(pIdxL,'g'); hold off

    % if sweep start was included, remove
    if pIdxL(1)<pIdxH
        pIdxL(1) = [];
    end

    % if first trial, find average sweep length
    if nt==1
        sweepDur = round(mean([pIdxH(2:end)-pIdxL(1:end) ; pIdxL(1:end)-pIdxH(1:end-1)]))-1;
    end
    
    % for each sweep
    % note: incomplete first and last sweep omitted
    for p=1:length(pIdxH)-1
        % pull this sweep indices for both left and right
        leftIdx = pIdxH(p):pIdxH(p)+sweepDur;
        rightIdx = pIdxL(p):pIdxL(p)+sweepDur;
        % pull data for both left and right
        leftwardSpikeRst = [leftwardSpikeRst spikerst(leftIdx,nt)];
        rightwardSpikeRst = [rightwardSpikeRst spikerst(rightIdx,nt)];
%         leftwardPanelPos = [leftwardPanelPos panelpos(leftIdx,nt)];
%         rightwardPanelPos = [rightwardPanelPos panelpos(rightIdx,nt)];
    end
    sweepTime = time(1:sweepDur+1)*1000; %msec
end


%% generate psth

% initialize
spike_psth = [];
leftwardSpikeTimes = [];
rightwardSpikeTimes = [];

% convert spike raster into spike times
for ns = 1:size(leftwardSpikeRst,2)
    leftwardSpikeTimes = [leftwardSpikeTimes ; sweepTime(find(leftwardSpikeRst(:,ns)==1))];
    rightwardSpikeTimes = [rightwardSpikeTimes ; sweepTime(find(rightwardSpikeRst(:,ns)==1))];
end

% generate histogram
tBins = 50; %size of bins, in msec.
tEdges = 0:tBins:round(sweepTime(end),-1);

% for leftward sweeps
[leftValues,bins] = histcounts(leftwardSpikeTimes,tEdges);
sweep_time = mean([bins(1:end-1) ; bins(2:end)],1);
spike_psth(:,1) = leftValues;
% for rightward sweeps
[rightValues,bins] = histcounts(rightwardSpikeTimes,tEdges);
spike_psth(:,2) = rightValues;


%% plot

if optPlot
    figure; set(gcf,'Position',[100 100 1000 500])
    lw = 2;
    plot(sweep_time,spike_psth(:,1),'Color', '#0032A0','LineWidth',lw)
    hold on
    plot(sweep_time,spike_psth(:,2),'Color', '#ff0080','LineWidth',lw)
    hold off
    ylabel('Frequency')
    xlabel('Time (msec.)')
end

end
