%pulse_v_behavior
%
% data processing function, for plotting the relationship between each
% motion pulse and the corresponding behavioral response of the fly
%
% INPUT
% panelps
% forward
% angular
% sideway
% ttime
% thisSpeed
% trialTypes - pulse speeds
%
% OUTPUT
% pulseTrials - trial data for each pulse
% pulseMeans - mean data for each pulse
%
% CREATED   05/01/2023 - MC
% UPDATED   08/08/2023 - MC added stationary pulses
%
function [pulseTrials, pulseMeans] = pulse_v_behavior(panelps,forward,angular,sideway,ttime,thisSpeed,trialTypes)
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
    % run ordering function to find motion pulses in order
    try
        [thisOrderR,thisOrderL] = order_motion_pulse(thisTrial,nSpeeds,thisSpeed,bufferWindow);
        % pull and store data based on ordered motion pulse indices
        % rows = data, columns = trials, z = sweeps
        for p = 1:size(thisOrderR,2)
            pulsePanelpsR(:,t,p) = panelps(thisOrderR(:,p),t);
            pulsePanelpsL(:,t,p) = panelps(thisOrderL(:,p),t);

            pulseForwardR(:,t,p) = forward(thisOrderR(:,p),t);
            pulseForwardL(:,t,p) = forward(thisOrderL(:,p),t);
            pulseAngularR(:,t,p) = angular(thisOrderR(:,p),t);
            pulseAngularL(:,t,p) = angular(thisOrderL(:,p),t);
            pulseSidewayR(:,t,p) = sideway(thisOrderR(:,p),t);
            pulseSidewayL(:,t,p) = sideway(thisOrderL(:,p),t);
        end
    catch
        pulsePanelpsR(:,t,p) = nan;
        pulsePanelpsL(:,t,p) = nan;

        pulseForwardR(:,t,p) = nan;
        pulseForwardL(:,t,p) = nan;
        pulseAngularR(:,t,p) = nan;
        pulseAngularL(:,t,p) = nan;
        pulseSidewayR(:,t,p) = nan;
        pulseSidewayL(:,t,p) = nan;
    end
end

% calculate means for each
meanPanelps(:,1,:) = mean(pulsePanelpsR,2,"omitnan");
meanPanelps(:,2,:) = mean(pulsePanelpsL,2,"omitnan");
meanForward(:,1,:) = mean(pulseForwardR,2,"omitnan");
meanForward(:,2,:) = mean(pulseForwardL,2,"omitnan");
meanAngular(:,1,:) = mean(pulseAngularR,2,"omitnan");
meanAngular(:,2,:) = mean(pulseAngularL,2,"omitnan");
meanSideway(:,1,:) = mean(pulseSidewayR,2,"omitnan");
meanSideway(:,2,:) = mean(pulseSidewayL,2,"omitnan");


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
figure; set(gcf,'Position',[100 100 1800 800])
colorlabels = {'#77AC30';'#0072BD';'#7E2F8E';'#D95319'};
lw = 1.5;
nSweep = size(pulsePanelpsR,3);
subx = nSweep*reSweep; % left and right OR stationary only
suby = 4; % panels, angular, sideway, forward
c = 1; %counter

tpulse = ttime(1:size(pulsePanelpsR,1))*1000;
trange = [0 max(tpulse)];
panelmax = round(max(meanPanelps,[],'all'),-1);
panelrange = [-panelmax panelmax];
angmax = 100;
%angmax = round(max(abs([meanAngularR meanAngularL]),[],'all'),-1);
angrange = [-angmax angmax];
sidmax = ceil(max(abs(meanSideway),[],'all'));
sidrange = [-sidmax sidmax];
fwdmax = ceil(max(meanForward,[],'all'));
fwdrange = [0 fwdmax];

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

% plot rightward angular
for s = 1:nSweep
    subplot(suby,subx,c)
    plot(tpulse,meanAngular(:,1,s),'Color',colorlabels{2},'LineWidth',lw)
    xlim(trange)
    ylim(angrange)
    yline(0)
    c = c+1; %update counter
end
if motionCheck
    % plot leftward angular
    for s = nSweep:-1:1
        subplot(suby,subx,c)
        plot(tpulse,meanAngular(:,2,s),'Color',colorlabels{2},'LineWidth',lw)
        xlim(trange)
        ylim(angrange)
        yline(0)
        c = c+1; %update counter
    end
end
xlabel('Time (msec)')
subplot(suby,subx,(nSweep*reSweep)+1); ylabel('Angular Vel (deg/sec)')

% plot rightward sideway
for s = 1:nSweep
    subplot(suby,subx,c)
    plot(tpulse,meanSideway(:,1,s),'Color',colorlabels{3},'LineWidth',lw)
    xlim(trange)
    ylim(sidrange)
    yline(0)
    c = c+1; %update counter
end
if motionCheck
    % plot leftward sideway
    for s = nSweep:-1:1
        subplot(suby,subx,c)
        plot(tpulse,meanSideway(:,2,s),'Color',colorlabels{3},'LineWidth',lw)
        xlim(trange)
        ylim(sidrange)
        yline(0)
        c = c+1; %update counter
    end
end
xlabel('Time (msec)')
subplot(suby,subx,(nSweep*2*reSweep)+1); ylabel('Sideway Vel (mm/sec)')

% plot rightward forward
for s = 1:nSweep
    subplot(suby,subx,c)
    plot(tpulse,meanForward(:,1,s),'Color',colorlabels{4},'LineWidth',lw)
    xlim(trange)
    ylim(fwdrange)
    yline(0)
    c = c+1; %update counter
end
if motionCheck
    % plot leftward sideway
    for s = nSweep:-1:1
        subplot(suby,subx,c)
        plot(tpulse,meanForward(:,2,s),'Color',colorlabels{4},'LineWidth',lw)
        xlim(trange)
        ylim(fwdrange)
        yline(0)
        c = c+1; %update counter
    end
end
xlabel('Time (msec)')
subplot(suby,subx,(nSweep*3*reSweep)+1); ylabel('Forward Vel (mm/sec)')


%% store variables for output

% store trials
pulseTrials.forwardR = pulseForwardR;
pulseTrials.forwardL = pulseForwardL;
pulseTrials.angularR = pulseAngularR;
pulseTrials.angularL = pulseAngularL;
pulseTrials.sidewayR = pulseSidewayR;
pulseTrials.sidewaydL = pulseSidewayL;

% store means
pulseMeans.panelps = meanPanelps;
pulseMeans.forward = meanForward;
pulseMeans.angular = meanAngular;
pulseMeans.sideway = meanSideway;


end