% barjump_analysis
%
% Analyzes bar jumps from a dataset and generates corresponding plots of the fly's 
% angular velocity and changes in cell spike rate within a specified time window 
% around each jump event.
%
% INPUTS:
% allPanelPs   - Array of panel positions (degrees).
% allJumpTrig  - Voltage indicators for jump events (trigger signals).
% allForward    - Array of forward velocity of the fly.
% allAngular    - Array of angular velocity of the fly.
% allSpikeRt    - Array of cell spike rates.
% ttime         - Time array (in seconds).
% optPlot       - Binary option; set to 0 to disable plotting, or 1 to enable plotting.
%
% SETTINGS:
% The time window before and after each jump can be adjusted via settings.
% Default maximum window size for including jumps is 15 seconds.
%
% OUTPUTS:
% jumpData     - Struct containing summary data such as mean and standard deviation
%                of panel positions, angular velocity, and spike rates.
% restore       - Struct containing details on jump correction, including timing 
%                and variance metrics.
%
% CREATED: 02/27/2024 by MC
% UPDATED: 03/07/2024 by MC 
% - Added variance calculations.
% - Introduced exclusion criteria for jumps where the fly had already turned 
%   past the setpoint at the time of the jump.
% - Added option to specify maximum time the fly can be stationary before and 
%   after a jump.
%
function [jumpData,restore] = barjump_analysis(allPanelPs,allJumpTrig,allForward,allAngular,allSpikeRt,ttime, optPlot)
%% settings
settings = processSettings();
plotWin = 10; %sec, time around jump to plot
ephysLog = length(allSpikeRt(:,1))>1; % check if spikert data present or empty

%% set jump parameters
% convert jump window to time
[pre_jumpWindow] = fetchTimeIdx(ttime,settings.mtaxPreWin);
[pst_jumpWindow] = fetchTimeIdx(ttime,settings.mtaxPstWin);
[assessWindow] = fetchTimeIdx(ttime,6);

% frame lag
% there is a small frame lag between when the trigger is sent
% and when the bar actually jumps according to the data... unclear why...
% this corresponds to ~15msec
frameLag = 0; %11

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
for v = settings.voltages_inuse
    allJumpTrig_deg(allJumpTrig_starts==v) = settings.jumps_inuse(v);
end
nJumps = length(unique(settings.jumps_inuse));

%% pull data according to jump size

% initialize
jumpPanelps = [];
jumpForward = [];
jumpAngular = [];
jumpSpikert = [];

% for each jump size
for v = 1:nJumps
    thisJump = settings.jumps_inuse(v);
    c = 1; %initialize counter

    % for each trial
    for t = 1:nTrials
        % pull index for this jump size
        jumpIdx = find(allJumpTrig_deg(:,t)==thisJump);
        for j = 1:length(jumpIdx)
            % pull before, full, and after jump windows
            mtaxPreWindow = jumpIdx(j)-pre_jumpWindow+frameLag:jumpIdx(j);
            thisFullWindow = jumpIdx(j)-pre_jumpWindow+frameLag:jumpIdx(j)+pst_jumpWindow+frameLag;
            mtaxPstWindow = jumpIdx(j):jumpIdx(j)+pst_jumpWindow+frameLag;

            % confirm that this jump was identified correctly
            if length(unique(allJumpTrig(thisFullWindow,t)))<=2
                % include this jump ONLY IF the fly was walking prior
                % find trials when the fly was running above threshold
                fwdThresh = 0.1;
                timeNotMovingPre = ttime(sum(allForward(mtaxPreWindow,t)<fwdThresh)+1);
                timeNotMovingPst = ttime(sum(allForward(mtaxPstWindow,t)<fwdThresh)+1);

                if timeNotMovingPre<=settings.maxTimeStop && timeNotMovingPst<=settings.maxTimeStop
                    % pull panel and behavior data
                    jumpPanelps(:,c,v) = allPanelPs(thisFullWindow,t);
                    jumpForward(:,c,v) = allForward(thisFullWindow,t);
                    jumpAngular(:,c,v) = allAngular(thisFullWindow,t);
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
windowTime = [-flip(ttime(2:pre_jumpWindow+1)) ttime(1:pst_jumpWindow+1)]; %sec


% transform panel data
jumpPanelps = unwrap(deg2rad(jumpPanelps)); %unwrap
jumpPanelps = rad2deg(jumpPanelps); %convert back to degrees
panelStarts = mean(jumpPanelps(1:pre_jumpWindow,:,:),1,'omitnan');
jumpPanelps = jumpPanelps - panelStarts;
% 180 jumps are a special case where the direction of turn is irrelevant
for v = find(abs(settings.jumps_inuse)==180)
    jumpPanelps(:,:,v) = wrapTo180(jumpPanelps(:,:,v));
    for c = 1:cFinal
        if mean(jumpPanelps(pre_jumpWindow:pre_jumpWindow+assessWindow,c,v))<0
            jumpPanelps(:,c,v) = jumpPanelps(:,c,v) *-1;
            jumpAngular(:,c,v) = jumpAngular(:,c,v) *-1;
        end
    end
     jumpPanelps(:,:,v) = abs(jumpPanelps(:,:,v));
     settings.jumps_inuse(end) = 180;
    % jumpAngular(:,:,v) = abs(jumpAngular(:,:,v));
end

%% determine if jump was corrected for

% initialize
correct_restore = [];
correct_jumpvar_pre = [];
correct_jumpvar_pst = [];

% for each jump
for v = 1:nJumps
    for c = 1:cFinal
        % for jumps to the right
        if settings.jumps_inuse(v)>0
            HD_restore = find(jumpPanelps(pre_jumpWindow:end,c,v)<1,1);
        % for jumps to the left
        else
            HD_restore = find(jumpPanelps(pre_jumpWindow:end,c,v)>1,1);
        end
        
        % if jump was corrected, convert to time (s)
        if ~isempty(HD_restore)&HD_restore>1
            HD_restoreT = ttime(HD_restore+1);
            % calculate and store HD variance
            HD_var_pre = var(jumpPanelps(1:pre_jumpWindow,c,v));
            HD_var_pst = var(jumpPanelps(HD_restore:end,c,v));
        % else exclude
        else
            HD_restoreT= nan;
            HD_var_pre = nan;
            HD_var_pst = nan;
        end
        % store data for this correction
        correct_restore(c,v) = HD_restoreT;
        correct_jumpvar_pre(c,v) = HD_var_pre;
        correct_jumpvar_pst(c,v) = HD_var_pst;

        
    end
end
% count the number of corrected jumps
correctedJumps = ~isnan(correct_restore);
trialsWithCorrections = sum(~isnan(correct_restore));

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
    thisCorrect = correctedJumps(:,v);
    panelpsMean(:,v) = mean(jumpPanelps(:,thisCorrect,v),2,'omitnan');
    angularMean(:,v) = mean(jumpAngular(:,thisCorrect,v),2,'omitnan');

    panelpsSTD(:,v) = std(jumpPanelps(:,thisCorrect,v),0,2,'omitnan');
    angularSTD(:,v) = std(jumpAngular(:,thisCorrect,v),0,2,'omitnan');

    if ephysLog
        spikertMean(:,v) = mean(jumpSpikert(:,thisCorrect,v),2,'omitnan');
        spikertSTD(:,v) = std(jumpSpikert(:,thisCorrect,v),0,2,'omitnan');
        jumpData.spikert = spikertMean;
        jumpData.spikertstd = spikertSTD;
        spikertMax = ceil(max(jumpSpikert,[],'all'));
    end

    % pull number of trials with running
    trialsWithRunning(v) = sum(~isnan(jumpPanelps(1,:,v)));
end
% calculate mean/SD for time to correct
time2restoreMean = mean(correct_restore,'omitnan');
time2restoreSTD = std(correct_restore,'omitnan');
mean_pre_jumpvar = median(correct_jumpvar_pre,'omitnan');
mean_jumpvar = median(correct_jumpvar_pst,'omitnan');

% repeat but combine +/- jumps
c = 1; %counter
for v = 1:2:6
    if v<5
        flp = -1;
    else
        flp = 1;
    end
    thisCorrectR = jumpPanelps(:,correctedJumps(:,v),v);
    thisCorrectL = jumpPanelps(:,correctedJumps(:,v+1),v+1)*flp;
    panelpsMean_P(:,c) = median([thisCorrectR thisCorrectL],2,'omitnan');

    % calculate mean/SD for time to correct
    time2restoreMean_P(c) = median(correct_restore(:,v:v+1),'all','omitnan');
    time2restoreSTD_P(c) = std(correct_restore(:,v:v+1),1,'all','omitnan');
    mean_pre_jumpvar_P(c) = median(correct_jumpvar_pre(:,v:v+1),'all','omitnan');
    mean_jumpvar_P(c) = median(correct_jumpvar_pst(:,v:v+1),'all','omitnan');

    % pull number of trials with corrections and running
    trialsWithCorrections_P(c) = trialsWithCorrections(v) + trialsWithCorrections(v+1);
    trialsWithRunning_P(c) = trialsWithRunning(v) + trialsWithRunning(v+1);

    c = c+1; %update counter
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

restore.all = correct_restore;
restore.mean = time2restoreMean;
restore.std = time2restoreSTD;
restore.prevar = mean_pre_jumpvar;
restore.var = mean_jumpvar;

restore.meanPool = time2restoreMean_P;
restore.stdPool = time2restoreSTD_P;
restore.prevarPool = mean_pre_jumpvar_P;
restore.varPool = mean_jumpvar_P;

restore.trialsWithRun = trialsWithRunning;
restore.trialsWithCorrect = trialsWithCorrections;
restore.trialsWithoutCorrect = trialsWithRunning - trialsWithCorrections;
restore.percent = trialsWithCorrections./trialsWithRunning;
restore.percentPool = trialsWithCorrections_P./trialsWithRunning_P;

%% plot data
% plot constraints
panelpsMax = 180;
angularMax = 300;

if optPlot
    % initialize
    nrow = 2+ephysLog;
    figure; set(gcf,'Position',[100 100 1800 900])
    tiledlayout(nrow,nJumps,'TileSpacing','compact')

    % for each jump size
    % plot panelps
    for v = 1:nJumps
        nexttile
        if sum(correctedJumps(:,v))>0
            plot(windowTime,jumpPanelps(:,correctedJumps(:,v),v), 'Color',[.8 .8 .8], 'LineWidth', settings.lwTri)
            hold on
        end
        plot(windowTime, panelpsMean(:,v),'Color', '#77AC30','LineWidth',settings.lwAvg)
        title([num2str(settings.jumps_inuse(v)) ' jump'])
        xlabel('Time (msec)')
        if v==1
            ylabel('object pos (deg)')
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
            if sum(correctedJumps(:,v))>0
                plot(windowTime,jumpSpikert(:,correctedJumps(:,v),v), 'Color',[.8 .8 .8], 'LineWidth', settings.lwTri)
                hold on
            end
            plot(windowTime, spikertMean(:,v),'Color', '#77AC30','LineWidth',settings.lwAvg)
            title([num2str(settings.jumps_inuse(v)) ' jump'])
            axis tight
            xlim([-2 plotWin])
            ylim([0 spikertMax])
            xline(0)
        end
    end
    % plot angular velocity
    for v = 1:nJumps
        nexttile
        if sum(correctedJumps(:,v))>0
            plot(windowTime,jumpAngular(:,correctedJumps(:,v),v), 'Color',[.8 .8 .8], 'LineWidth', settings.lwTri)
            hold on
        end
        plot(windowTime, angularMean(:,v),'Color', '#0072BD','LineWidth',settings.lwAvg)
        title([num2str(settings.jumps_inuse(v)) ' jump'])
        xlabel('Time (msec)')
        if v==1
            ylabel('Angular Velocity (deg/sec)')
        end
        axis tight
        xlim([-2 plotWin])
        ylim([-angularMax angularMax])
        yline(0)
        xline(0)
    end
end


end