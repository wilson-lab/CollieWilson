% iinj_v_behavior_thresh
%
% analysis function that pulls the behavior responses from repeated current
% injections AND threshold according to the fly's behavior before, during,
% and/or after the pulse.
%
% INPUT
% iinject - current stim
% spikert - spikes/sec
% forward - velocity
% angular - velocity (can be normalized)
% sideway - velocity (can be normalized)
% ttime - time
% threshOpt - select threshold
% optPlot - 0 no plot, 1 to plot
%
% OUTPUT
% avgResponse - mean responses for each directional velocity
%
% CREATED: 12/12/23 MC adapted from iinj_v_behavior
%

function [avgResponse] = iinj_v_behavior_thresh(iinject,spikert,forward,angular,sideway,ttime,threshOpt, optPlot)
%% initialize
stepBuffer = 2000;
nTrials = size(iinject,2); %number of trials

tmsec = ttime*1000; %convert to msec

sr_limit = [0 14];
% for all data
f_limit = [0 6];
a_limit = [-100 100];
s_limit = [-1 2];
% for zoom data
f_limitz = [-.25 0.5];
a_limitz = [-20 20];
s_limitz = [-1 1];


%% find iinject pulses

% initialize
pulseSpikeRt = [];
pulseForward = [];
pulseAngular = [];
pulseSideway = [];
dp = 1;
hp = 1;

% for each trial
for nt = 1:nTrials
    % find all start/stops based on iinj output changes
    startstopidx = find(ischange(iinject(:,nt)));
    nSteps = length(startstopidx);

    % for each start step
    for ns = 1:2:nSteps
        thisStepVal = round(iinject(startstopidx(ns)+1,nt),-1); %step value
        thisStepIdx = startstopidx(ns)-stepBuffer:startstopidx(ns+1)+stepBuffer; %step idx

        % for + pulse
        if thisStepVal>0
            pulseIinject(:,dp,1) = iinject(thisStepIdx,nt);
            pulseSpikeRt(:,dp,1) = spikert(thisStepIdx,nt);
            pulseForward(:,dp,1) = forward(thisStepIdx,nt);
            pulseAngular(:,dp,1) = angular(thisStepIdx,nt);
            pulseSideway(:,dp,1) = sideway(thisStepIdx,nt);
            dp = dp+1;
            % for - pulse
        else
            pulseIinject(:,hp,2) = iinject(thisStepIdx,nt);
            pulseSpikeRt(:,hp,2) = spikert(thisStepIdx,nt);
            pulseForward(:,hp,2) = forward(thisStepIdx,nt);
            pulseAngular(:,hp,2) = angular(thisStepIdx,nt);
            pulseSideway(:,hp,2) = sideway(thisStepIdx,nt);
            hp = hp+1;
        end
    end
end
pulseDur = size(thisStepIdx,2);
pulseIdx = [tmsec(stepBuffer), tmsec(pulseDur-stepBuffer)];
pulseTime = tmsec(1:pulseDur);


%% threshold according to fly behavior during each pulse
% set pulse parameters
pulseStartStop = find(ischange(pulseIinject(:,1,1)));
% threshold parameters
fwdThresh = 2; %mm/sec, min fwd speed
runTime = 1; %points, min time spend running
minTrials = 5; %minimum number of trials the fly must have run to be considered

% determine which trials pass run threshold (both magnitude and duration)
select_full = sum(pulseForward>fwdThresh,1)>runTime; %full
select_bef = sum(pulseForward(1:pulseStartStop(1),:,:)>fwdThresh,1)>runTime; %before pulse
select_dur = sum(pulseForward(pulseStartStop(1):pulseStartStop(2),:,:)>fwdThresh,1)>runTime; %during pulse
select_aft = sum(pulseForward(pulseStartStop(2):pulseDur,:,:)>fwdThresh,1)>runTime; %after pulse


%% plot
if optPlot
    figure; set(gcf,'Position',[100 100 1500 800])
    colorlabels = {'#77AC30';'#D95319';'#0072BD';'#7E2F8E';[0.7 0.7 0.7]};
    lw = 1.5;

    % for each threshold window option
    for th = 1:4
        % select data according to threshold window
        switch th
            case 1
                this_select = select_full;
            case 2
                this_select = select_bef;
            case 3
                this_select = select_dur;
            case 4
                this_select = select_aft;
        end

        if sum(this_select(:,:,1))>=minTrials
            % plot spikerate
            subplot(4,8,th)
            plot(pulseTime,pulseSpikeRt(:,this_select(1,:,1),1),':','Color',colorlabels{5}); hold on
            meanSpikeRt(:,1) = mean(pulseSpikeRt(:,this_select(1,:,1),1),2,'omitnan');
            plot(pulseTime,meanSpikeRt(:,1),'LineWidth',lw,'Color',colorlabels{1});
            ylabel('Spikerate')
            ylim(sr_limit)

            % plot forward
            subplot(4,8,th+8)
            plot(pulseTime,pulseForward(:,this_select(1,:,1),1),':','Color',colorlabels{5}); hold on
            meanForward(:,1) = mean(pulseForward(:,this_select(1,:,1),1),2,'omitnan');
            plot(pulseTime,meanForward(:,1),'LineWidth',lw,'Color',colorlabels{2});
            ylabel('Forward (mm/s)')
            ylim(f_limit)

            % plot angular
            subplot(4,8,th+16)
            plot(pulseTime,pulseAngular(:,this_select(1,:,1),1),':','Color',colorlabels{5}); hold on
            meanAngular(:,1) = mean(pulseAngular(:,this_select(1,:,1),1),2,'omitnan');
            plot(pulseTime,meanAngular(:,1),'LineWidth',lw,'Color',colorlabels{3});
            ylabel('Angular(d/s)')
            ylim(a_limit)
            yline(0)

            % plot sideway
            subplot(4,8,th+24)
            plot(pulseTime,pulseSideway(:,this_select(1,:,1),1),':','Color',colorlabels{5}); hold on
            meanSideway(:,1) = mean(pulseSideway(:,this_select(1,:,1),1),2,'omitnan');
            plot(pulseTime,meanSideway(:,1),'LineWidth',lw,'Color',colorlabels{4});
            ylabel('Sideway(mm/s)')
            ylim(s_limit)
            yline(0)
        else
            meanSpikeRt(:,1) = nan(pulseDur,1);
            meanForward(:,1) = nan(pulseDur,1);
            meanAngular(:,1) = nan(pulseDur,1);
            meanSideway(:,1) = nan(pulseDur,1);
        end

        if sum(this_select(:,:,2))>=minTrials
            % plot spikerate
            subplot(4,8,th+4)
            plot(pulseTime,pulseSpikeRt(:,this_select(1,:,2),2),':','Color',colorlabels{5}); hold on
            meanSpikeRt(:,2) = mean(pulseSpikeRt(:,this_select(1,:,2),2),2,'omitnan');
            plot(pulseTime,meanSpikeRt(:,2),'LineWidth',lw,'Color',colorlabels{1});
            ylabel('Spikerate')
            ylim(sr_limit)

            % plot forward
            subplot(4,8,th+12)
            plot(pulseTime,pulseForward(:,this_select(1,:,2),2),':','Color',colorlabels{5}); hold on
            meanForward(:,2) = mean(pulseForward(:,this_select(1,:,2),2),2,'omitnan');
            plot(pulseTime,meanForward(:,2),'LineWidth',lw,'Color',colorlabels{2});
            ylabel('Forward(mm/s)')
            ylim(f_limit)

            % plot angular
            subplot(4,8,th+20)
            plot(pulseTime,pulseAngular(:,this_select(1,:,2),2),':','Color',colorlabels{5}); hold on
            meanAngular(:,2) = mean(pulseAngular(:,this_select(1,:,2),2),2,'omitnan');
            plot(pulseTime,meanAngular(:,2),'LineWidth',lw,'Color',colorlabels{3});
            ylabel('Angular(d/s)')
            ylim(a_limit)
            yline(0)

            % plot sideway
            subplot(4,8,th+28)
            plot(pulseTime,pulseSideway(:,this_select(1,:,2),2),':','Color',colorlabels{5}); hold on
            meanSideway(:,2) = mean(pulseSideway(:,this_select(1,:,2),2),2,'omitnan');
            plot(pulseTime,meanSideway(:,2),'LineWidth',lw,'Color',colorlabels{4});
            ylabel('Sideway(mm/s)')
            ylim(s_limit)
            yline(0)
        else
            meanSpikeRt(:,2) = nan(pulseDur,1);
            meanForward(:,2) = nan(pulseDur,1);
            meanAngular(:,2) = nan(pulseDur,1);
            meanSideway(:,2) = nan(pulseDur,1);
        end

        % store means for output based on input request
        if th == threshOpt
            avgResponse.iinject = reshape(mean(pulseIinject,2,'omitnan'),[],2);
            avgResponse.spikert = meanSpikeRt;
            avgResponse.forward = meanForward;
            avgResponse.angular = meanAngular;
            avgResponse.sideway = meanSideway;
        end
    end


end

end