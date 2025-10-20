% analysis function
%
% Pulls bar jumps from a dataset and plots the fly's corresponding angular
% velocity and changes in cell spikerate within a window around said jump.
%
% NOTE: does not threshold according to corrections, looks at ALL jumps
%
% INPUTS
% allPanelps - panel positions
% allJumpTrig - voltage indicators of jumps
% allForward - fly forward velocity
% allAngular - fly angular velocity
% allSpikeRt - cell spike rate
% exptTime - time array
% maxTimeStopped - max amount of time the fly can be stopped before and
% after a jump before the bout is disqualified (sec). To include all, set
% to 15 (max window size)
% optPlot - 0 for no plot, 1 for plot
%
% OUTPUTS
% jumpData - summary
%
% CREATED: 03/07/2024 MC modified from barjump_analysis to not threshold
% for corrections
%

function jumpData = barjump_nocorrect_analysis(allPanelPs,allJumpTrig,allForward,allAngular,allSpikeRt,exptTime, maxTimeStopped, optPlot)
%% settings

lw1 = 0.5; %trials
lw2 = 2; %means
plotWin = 10; %sec, time around jump to plot
ephysLog = length(allSpikeRt(:,1))>1; % check if spikert data present or empty

%% set jump parameters
% set jump window size for pulling data before/after
preWin = 15; %sec, window before jump
pstWin = 10; %sec, window after jump
% convert jump window to time
f = 2000;
pre_jumpWindow= preWin*f;
pst_jumpWindow= pstWin*f;

% frame lag
% there is a small frame lag between when the trigger is sent
% and when the bar actually jumps according to the data... unclear why...
% this corresponds to ~15msec
frameLag = 41;

% NOTE: the following parameters are hardcoded in python socket-client code
% set voltages that were used to denote jump triggers
voltages_inuse = 1:6;

% set jump sizes corresponding to said jump triggers
jumps_inuse = [30, -30, 60, -60, 180, -180]; %deg
turn_norm = [1, 1, 1, 1, 1, 1]; %do not flip turn vel

%% process jump trigger

% initialize
trialDur = size(allJumpTrig,1);
nTrials = size(allJumpTrig,2);

allJumpTrig_starts = zeros(trialDur,nTrials);
allJumpTrig_deg = zeros(trialDur,nTrials);

% pull only the first trigger value
% in most cases, a single jump trigger was sampled repeatedly making
% indexing harder
for t = 1:nTrials
    isolateTrigChanges = diff(allJumpTrig(:,t));
    trigStartIdx = find(isolateTrigChanges<0)+1;
    allJumpTrig_starts(trigStartIdx,t) = allJumpTrig(trigStartIdx-1,t);
end

% convert trigger voltages to degrees
for v = voltages_inuse
    allJumpTrig_deg(allJumpTrig_starts==v) = jumps_inuse(v);
end
nJumps = length(unique(jumps_inuse));

%% pull data according to jump size

% initialize
jumpPanelps = [];
jumpForward = [];
jumpAngular = [];
jumpSpikert = [];

% for each jump size
for v = 1:nJumps
    thisJump = jumps_inuse(v);
    c = 1; %initialize counter

    % for each trial
    for t = 1:nTrials
        % pull index for this jump size
        jumpIdx = find(allJumpTrig_deg(:,t)==thisJump);
        for j = 1:length(jumpIdx)
            % pull before, full, and after jump windows
            thisPreWindow = jumpIdx(j)-pre_jumpWindow+frameLag:jumpIdx(j);
            thisFullWindow = jumpIdx(j)-pre_jumpWindow+frameLag:jumpIdx(j)+pst_jumpWindow+frameLag;
            thisPstWindow = jumpIdx(j):jumpIdx(j)+pst_jumpWindow+frameLag;

            % confirm that this jump was identified correctly
            if length(unique(allJumpTrig(thisFullWindow,t)))<=2
                % include this jump ONLY IF the fly was walking prior
                % find trials when the fly was running above threshold
                fwdThresh = 1;
                timeNotMovingPre = exptTime(sum(allForward(thisPreWindow,t)<fwdThresh)+1);
                timeNotMovingPst = exptTime(sum(allForward(thisPstWindow,t)<fwdThresh)+1);

                if timeNotMovingPre<=maxTimeStopped && timeNotMovingPst<=maxTimeStopped
                    % pull panel and behavior data
                    jumpPanelps(:,c,v) = allPanelPs(thisFullWindow,t);
                    jumpForward(:,c,v) = allForward(thisFullWindow,t);
                    jumpAngular(:,c,v) = allAngular(thisFullWindow,t) * turn_norm(v);
                    % pull ephys data, if collected
                    if ephysLog
                        jumpSpikert(:,c,v) = allSpikeRt(thisFullWindow,t);
                    end

                else
                    % else, this jump should be excluded (nans)
                    jumpPanelps(:,c,v) = nan(length(thisFullWindow),1);
                    jumpForward(:,c,v) = nan(length(thisFullWindow),1);
                    jumpAngular(:,c,v) = nan(length(thisFullWindow),1);
                    if ephysLog
                        jumpSpikert(:,c,v) = nan(length(thisFullWindow),1);
                    end
                end
                % update counter
                c = c+1;
            end
        end
    end
end
% final counter
cFinal = c-1;

% pull window time, with zero indicating the jump trigger
windowTime = [-flip(exptTime(2:pre_jumpWindow+1)) exptTime(1:pst_jumpWindow+1)]; %sec


% transform panel data
jumpPanelps = unwrap(deg2rad(jumpPanelps)); %unwrap
jumpPanelps = rad2deg(jumpPanelps); %convert back to degrees
panelStarts = mean(jumpPanelps(1:pre_jumpWindow,:,:),1,'omitnan');
jumpPanelps = jumpPanelps - panelStarts;
% 180 jumps are a special case where the direction of turn is irrelevant
for v = find(abs(jumps_inuse)==180)
    jumpPanelps(:,:,v) = wrapTo180(jumpPanelps(:,:,v));
    for c = 1:cFinal
        assessWindow = 6 * f; %sec x f
        if mean(jumpPanelps(pre_jumpWindow:pre_jumpWindow+assessWindow,c,v))<0
            jumpPanelps(:,c,v) = jumpPanelps(:,c,v) *-1;
            jumpAngular(:,c,v) = jumpAngular(:,c,v) *-1;
        end
    end
     jumpPanelps(:,:,v) = abs(jumpPanelps(:,:,v));
     jumps_inuse(end) = 180;
    % jumpAngular(:,:,v) = abs(jumpAngular(:,:,v));
end


%% calculate mean and SD

% initialize
panelpsMean = [];
panelpsMean_P = [];
angularMean = [];
spikertMean = [];
panelpsSTD = [];
angularSTD = [];
spikertSTD = [];

% calculate mean/SD each experimental variable
for v = 1:nJumps
    
    panelpsMean(:,v) = mean(jumpPanelps(:,:,v),2,'omitnan');
    angularMean(:,v) = mean(jumpAngular(:,:,v),2,'omitnan');

    panelpsSTD(:,v) = std(jumpPanelps(:,:,v),0,2,'omitnan');
    angularSTD(:,v) = std(jumpAngular(:,:,v),0,2,'omitnan');

    if ephysLog
        spikertMean(:,v) = mean(jumpSpikert(:,:,v),2,'omitnan');
        spikertSTD(:,v) = std(jumpSpikert(:,:,v),0,2,'omitnan');
        jumpData.spikert = spikertMean;
        jumpData.spikertstd = spikertSTD;
        spikertMax = ceil(max(jumpSpikert,[],'all'));
    end

    % pull number of trials with running
    trialsWithRunning(v) = sum(~isnan(jumpPanelps(1,:,v)));
end

% store for output
jumpData.panelps = panelpsMean;
jumpData.panelpsPool = panelpsMean_P;
jumpData.angular = angularMean;
jumpData.panelpsstd = panelpsSTD;
jumpData.angularstd = angularSTD;
jumpData.time = windowTime;
jumpData.fwdMax = reshape(max(max(jumpForward)),[],nJumps);
jumpData.angMax = reshape(max(max(abs(jumpAngular))),[],nJumps);


%% plot data
% plot constraints
panelpsMax = 180;
angularMax = 100;

if optPlot
    % initialize
    nrow = 2+ephysLog;
    figure; set(gcf,'Position',[100 100 1800 900])
    tiledlayout(nrow,nJumps,'TileSpacing','compact')

    % for each jump size
    % plot panelps
    for v = 1:nJumps
        nexttile
        plot(windowTime,jumpPanelps(:,:,v), 'Color',[.8 .8 .8], 'LineWidth', lw1)
        hold on
        plot(windowTime, panelpsMean(:,v),'Color', '#77AC30','LineWidth',lw2)
        title([num2str(jumps_inuse(v)) ' jump'])
        xlabel('Time (msec)')
        if v==1
            ylabel('HD (deg)')
        end
        axis tight
        xlim([-2 plotWin])
        ylim([-panelpsMax panelpsMax])
        yline(0)
        xline(0)
    end
    % plot spikerate (if present)
    if ephysLog
        for v = 1:nJumps
            nexttile
            plot(windowTime,jumpSpikert(:,:,v), 'Color',[.8 .8 .8], 'LineWidth', lw1)
            hold on
            plot(windowTime, spikertMean(:,v),'Color', '#77AC30','LineWidth',lw2)
            title([num2str(jumps_inuse(v)) ' jump'])
            axis tight
            xlim([-2 plotWin])
            ylim([0 spikertMax])
            xline(0)
        end
    end
    % plot angular velocity
    for v = 1:nJumps
        nexttile
        plot(windowTime,jumpAngular(:,:,v), 'Color',[.8 .8 .8], 'LineWidth', lw1)
        hold on
        plot(windowTime, angularMean(:,v),'Color', '#0072BD','LineWidth',lw2)
        title([num2str(jumps_inuse(v)) ' jump'])
        xlabel('Time (msec)')
        if v==1
            ylabel('Angular Vel (deg/sec)')
        end
        axis tight
        xlim([-2 plotWin])
        ylim([-angularMax angularMax])
        yline(0)
        xline(0)
    end
end


end