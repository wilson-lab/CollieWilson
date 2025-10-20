%pulse_v_spikert
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
function [pulseTrials, pulseMeans] = pulse_v_spikert(panelps,spikert,ttime,thisSpeed,trialTypes)
%% initialize
nTrial = size(panelps,2);
nSpeeds = length(trialTypes);

%% pull motion pulses from pseudorandomized dataset

% initialize
bufferWindow = 100;

% for each trial
for t = 1:nTrial
    % pull this panel data
    thisTrial = panelps(:,t);
    try
        % run ordering function to find motion pulses in order
        [thisOrderR,thisOrderL] = order_motion_pulse(thisTrial,nSpeeds,thisSpeed,bufferWindow);

        % pull and store data based on ordered motion pulse indices
        % rows = data, columns = trials, z = sweeps
        for p = 1:size(thisOrderR,2)
            pulsePanelpsR(:,t,p) = panelps(thisOrderR(:,p),t);
            pulsePanelpsL(:,t,p) = panelps(thisOrderL(:,p),t);

            pulseSpikeR(:,t,p) = spikert(thisOrderR(:,p),t);
            pulseSpikeL(:,t,p) = spikert(thisOrderL(:,p),t);
        end
    catch
        for p = 1:size(thisOrderR,2)
            pulsePanelpsR(:,t,p) = nan;
            pulsePanelpsL(:,t,p) = nan;

            pulseSpikeR(:,t,p) = nan;
            pulseSpikeL(:,t,p) = nan;
        end
    end
end

% calculate means for each
meanPanelps(:,1,:) = mean(pulsePanelpsR,2,"omitnan");
meanPanelps(:,2,:) = mean(pulsePanelpsL,2,"omitnan");
meanSpike(:,1,:) = mean(pulseSpikeR,2,"omitnan");
meanSpike(:,2,:) = mean(pulseSpikeL,2,"omitnan");


%% motion check
% determine if bar swept or stationary
motionCheck = abs(sum(thisOrderR(:,1)-thisOrderL(:,1)))>0;
if motionCheck
    reSweep = 2;
else
    reSweep = 1;
end


%% plot

% initialize
figure; set(gcf,'Position',[100 100 1800 600])
colorlabels = {'#77AC30','k'};
lw = 1.5;
nSweep = size(pulsePanelpsR,3);
subx = nSweep*reSweep; % all left and right sweeps
suby = 2; % panels, spikert
c = 1; %counter

tpulse = ttime(1:size(pulsePanelpsR,1))*1000;
trange = [0 max(tpulse)];
panelmax = round(max(meanPanelps,[],'all'),-1);
panelrange = [-panelmax panelmax];
spkmax = ceil(max(meanSpike,[],'all'));
if spkmax>25
    spkrange = [0 spkmax];
elseif spkmax<10
    spkrange = [0 10];
else
    spkrange = [0 25];
end

% plot rightward target
for s = 1:nSweep
    subplot(suby,subx,c)
    plot(tpulse,meanPanelps(:,1,s),'Color',colorlabels{1},'LineWidth',lw)
    xlim(trange)
    ylim(panelrange)
    yline(0)
    c = c+1; %update counter
end
if motionCheck
    % plot leftward target
    for s = nSweep:-1:1
        subplot(suby,subx,c)
        plot(tpulse,meanPanelps(:,2,s),'Color',colorlabels{1},'LineWidth',lw)
        xlim(trange)
        ylim(panelrange)
        yline(0)
        c = c+1; %update counter
    end
end
xlabel('Time (msec)')
subplot(suby,subx,1); ylabel('Target Pos (deg)')

% plot rightward spikert
for s = 1:nSweep
    subplot(suby,subx,c)
    plot(tpulse,meanSpike(:,1,s),'Color',colorlabels{2},'LineWidth',lw)
    xlim(trange)
    ylim(spkrange)
    yline(0)
    c = c+1; %update counter
end
if motionCheck
    % plot leftward spikert
    for s = nSweep:-1:1
        subplot(suby,subx,c)
        plot(tpulse,meanSpike(:,2,s),'Color',colorlabels{2},'LineWidth',lw)
        xlim(trange)
        ylim(spkrange)
        yline(0)
        c = c+1; %update counter
    end
end
xlabel('Time (msec)')
subplot(suby,subx,subx+1); ylabel('Spikerate (spk/sec)')


%% store variables for output

% store trials
pulseTrials.spikertR = pulseSpikeR;
pulseTrials.spikertL = pulseSpikeL;

% store means
pulseMeans.panelps = meanPanelps;
pulseMeans.spikert = meanSpike;


end