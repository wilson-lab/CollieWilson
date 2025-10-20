%pulse_v_spikert_overlay
%
% data processing function, for plotting the relationship between each
% motion pulse and the corresponding spike response of the cell
%
% INPUT
% panelps
% spikert
% ttime
% thisSpeed
% trialTypes - pulse speeds
%
% OUTPUT
% pulseTrials - trial data for each pulse
% pulseMeans - mean data for each pulse
%
% CREATED
% 05/02/2023 - MC
%
function pulse_v_spikert_overlayacross(panelps,spikert,ttime,trialTypes,nFly)
%% initialize
nTrial = size(panelps,2);
nSpeeds = length(trialTypes);

%% pull motion pulses from pseudorandomized dataset

% initialize
bufferWindow = [165 100];

% for each speed
for sp = 1:nSpeeds
    % initialize
    pulsePanelpsR = [];
    pulsePanelpsL = [];
    pulseSpikertR = [];
    pulseSpikertL = [];
    % for each trial
    for t = 1:nTrial
        % pull this panel data
        thisTrial = panelps(:,t);
        try
            % run ordering function to find motion pulses in order
            [thisOrderR,thisOrderL] = order_motion_pulse(thisTrial,nSpeeds,sp,bufferWindow(sp));
            % if any indices fall outside of trial duration, fix
            thisOrderR(thisOrderR>size(panelps,1)) = size(panelps,1);
            thisOrderL(thisOrderL>size(panelps,1)) = size(panelps,1);

            % pull and store data based on ordered motion pulse indices
            % rows = data, columns = trials, z = sweeps
            nSweeps = size(thisOrderR,2);
            sweepDur(sp) = size(thisOrderR,1);
            for p = 1:nSweeps
                pulsePanelpsR(:,t,p) = panelps(thisOrderR(:,p),t);
                pulsePanelpsL(:,t,p) = panelps(thisOrderL(:,p),t);

                pulseSpikertR(:,t,p) = spikert(thisOrderR(:,p),t);
                pulseSpikertL(:,t,p) = spikert(thisOrderL(:,p),t);

                % pull peak SR for each trial to calculate error
                sweepIdx = find(~isnan(pulsePanelpsR(:,t,p)));
                peakSpikertR(t,p,sp) = max(pulseSpikertR(sweepIdx,t,p));
                peakSpikertL(t,p,sp) = max(pulseSpikertL(sweepIdx,t,p));
            end
        catch
            for p = 1:nSweeps
                pulsePanelpsR(:,t,p) = nan;
                pulsePanelpsL(:,t,p) = nan;

                pulseSpikertR(:,t,p) = nan;
                pulseSpikertL(:,t,p) = nan;

                % pull peak SR for each trial to calculate error
                peakSpikertR(t,p,sp) = nan;
                peakSpikertL(t,p,sp) = nan;
            end
        end
    end
    % calculate means for each
    meanPanelpsR{sp} = reshape(mean(pulsePanelpsR,2,"omitnan"),[],nSweeps);
    meanPanelpsL{sp} = reshape(mean(pulsePanelpsL,2,"omitnan"),[],nSweeps);
    meanSpikertR{sp} = reshape(mean(pulseSpikertR,2,"omitnan"),[],nSweeps);
    meanSpikertL{sp} = reshape(mean(pulseSpikertL,2,"omitnan"),[],nSweeps);

    % pull peak SR for each mean
    peakMeanSpikeRtR(sp,:) = max(meanSpikertR{sp}(sweepIdx,:));
    peakMeanSpikeRtL(sp,:) = max(meanSpikertL{sp}(sweepIdx,:));
end

% calculate peak SEM
semSpikertR = reshape(std(peakSpikertR,[],1)/sqrt(nFly),2,[]);
semSpikertL = reshape(std(peakSpikertL,[],1)/sqrt(nFly),2,[]);

%% reshape for plotting
maxLength = max(sweepDur);
for sp = 1:nSpeeds
    meanPanelpsR_rs(:,:,sp) = imresize(meanPanelpsR{sp},[maxLength nSweeps], 'nearest');
    meanPanelpsL_rs(:,:,sp) = imresize(meanPanelpsL{sp},[maxLength nSweeps], 'nearest');
    meanSpikertR_rs(:,:,sp) = imresize(meanSpikertR{sp},[maxLength nSweeps], 'nearest');
    meanSpikertL_rs(:,:,sp) = imresize(meanSpikertL{sp},[maxLength nSweeps], 'nearest');
end

%% plot

% initialize
figure; set(gcf,'Position',[100 100 1600 800])
subx = nSweeps+1; % separate by left and right
suby = 6; % panels, spikert

% store start positions
startPos = round(reshape(meanPanelpsR_rs(bufferWindow(2)+1,:,2),[],1));
nSweep = size(startPos,1);
nSweep_ipsi = sum(startPos>=0)+1;
nSweep_cntr = nSweep-nSweep_ipsi;

x = 1:size(meanPanelpsR_rs,1);
xrange = [0 max(x)];
panelmax = ceil(max(meanPanelpsR_rs,[],'all'));
panelrange = [-panelmax panelmax];
%spkmin = floor(min(meanSpike,[],'all'));
spkmax = ceil(max(meanSpikertR_rs,[],'all'));
if spkmax>25
    spkrange = [0 spkmax];
elseif spkmax<10
    spkrange = [0 10];
else
    spkrange = [0 25];
end

% set color properties
lw = 1.5;
linecolors = {"#0072BD";"#7E2F8E";"#4DBEEE"};

for p = 1:nSweeps
    % plot rightward panel position
    subplot(suby,subx,p)
    for sp = 1:nSpeeds
        plot(x,meanPanelpsR_rs(:,p,sp),'Color',linecolors{sp},'LineWidth',lw)
        hold on
    end
    if p == 1
        ylabel('Rightward Target (deg)')
    end
    ylim(panelrange)
    yline(0)
    xlim(xrange)
    xticklabels([])

    % plot rightward firing rate
    subplot(suby,subx,[p+subx p+subx*2])
    for sp = 1:nSpeeds
        plot(x,meanSpikertR_rs(:,p,sp),'Color',linecolors{sp},'LineWidth',lw)
        hold on
    end
    if p == 1
        ylabel('Rightward Avg SR (spikes/sec)')
    end
    ylim(spkrange)
    xlim(xrange)
    xticklabels([])


    % plot leftward panel position
    subplot(suby,subx,p+subx*3)
    for sp = 1:nSpeeds
        plot(x,meanPanelpsL_rs(:,p,sp),'Color',linecolors{sp},'LineWidth',lw)
        hold on
    end
    if p == 1
        ylabel('Leftward Target (deg)')
    end
    ylim(panelrange)
    yline(0)
    xlim(xrange)
    xticklabels([])

    % plot rightward firing rate
    subplot(suby,subx,[p+subx*4 p+subx*5])
    for sp = 1:nSpeeds
        plot(x,meanSpikertL_rs(:,p,sp),'Color',linecolors{sp},'LineWidth',lw)
        hold on
    end
    if p == 1
        ylabel('Leftward Avg SR (spikes/sec)')
    end
    ylim(spkrange)
    xlim(xrange)
    xticklabels([])
end

% plot peak FR
subplot(suby,subx,[subx subx*2 subx*3])
for sp = 1:nSpeeds
    plot(startPos,peakMeanSpikeRtR(sp,:),'--o','Color',linecolors{sp})
    %errorbar(1:nSweeps,peakMeanSpikeRtR(sp,:),semSpikertR(sp,:),'o','Color',linecolors{sp})
    hold on
end
ylabel('Peak Rightward SR (spikes/sec)')
xlim(panelrange)
ylim(spkrange)
xticklabels([])

subplot(suby,subx,[subx*4 subx*5 subx*6])
for sp = 1:nSpeeds
    plot(startPos,peakMeanSpikeRtL(sp,:),'--o','Color',linecolors{sp})
    %errorbar(1:nSweeps,peakMeanSpikeRtL(sp,:),semSpikertL(sp,:),'o','Color',linecolors{sp})
    hold on
end
ylabel('Peak Leftward SR (spikes/sec)')
xlim(panelrange)
ylim(spkrange)
xticklabels([])


end