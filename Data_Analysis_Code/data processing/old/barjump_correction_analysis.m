% analysis function
%
% Pulls bar jumps from a dataset to determine which jumps the fly
% successfully corrected for within the specified window
%
% INPUTS
% allPanelps - panel positions of bar
% allJumpTrig - voltage indicators of jumps
% allAngular - angular velocity of fly
% exptTime - time array
% optPlot - 0 for no plot, 1 for plot
%
% OUTPUTS
% jumpData - summary
%
% CREATED: 02/13/2024 MC
%

function jumpData = barjump_correction_analysis(allPanelPs,allJumpTrig,allForward,allAngular,exptTime,optPlot)
%% settings
% based on Westeinde Wilson 2024 classifications for bar jumps

% set jump window size for pulling data before/after
preWin = 15; %sec, window before jump
prePlotWin = preWin-2; %sec, for plotting window before jump
pstWin = 10; %sec, window after jump

% frame lag
% there is a small frame lag between when the trigger is sent
% and when the bar actually jumps according to the data... unclear why...
% this corresponds to ~15msec
frameLag = 41;

% NOTE: the following parameters are hardcoded in python socket-client code
% set voltages that were used to denote jump triggers
voltages_inuse = 1:6;
% set jump sizes corresponding to said jump triggers
% note, 180 L vs R is the same, so v = 5 and 6 will be combined
jumps_inuse = [30, -30, 60, -60, 180, 180]; %deg
% set range for a jump to be considered "corrected for"
% values represent post jump heading must be within X degrees of original value
correct_min = abs(jumps_inuse)/3; %deg


%% process jump trigger

% initialize
trialDur = size(allJumpTrig,1);
nTrials = size(allJumpTrig,2);

allJumpTrig_starts = zeros(trialDur,nTrials);
allJumpTrig_deg = zeros(trialDur,nTrials);

% convert jump window to time
f = 2000;
preWint = preWin*f; %sec
prePlotWint = prePlotWin*f; %sec
pstWint = pstWin*f; %sec

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

%% pull jump window data for each jump

% initialize
panelpsData = [];
forwardData = [];
angularData = [];
correct_restore = [];
correct_partial = [];

% for each jump size
for v = 1:nJumps
    thisJump = jumps_inuse(v);
    c = 1; %initialize counter

    % for each trial
    for t = 1:nTrials
        % pull index for this jump size
        jumpIdx = find(allJumpTrig_deg(:,t)==thisJump);
        for j = 1:length(jumpIdx)
            % pull indices for this window
            thisWindow = jumpIdx(j)-preWint+frameLag:jumpIdx(j)+pstWint+frameLag;

            % store data for this jump
            panelpsData(:,c,v) = allPanelPs(thisWindow,t);
            forwardData(:,c,v) = allForward(thisWindow,t);
            angularData(:,c,v) = allAngular(thisWindow,t);
            c = c+1; %update counter
        end

    end
end
% number of times each jump was presented in total
nPresent = c-1;

% pull window time, with zero indicating the jump trigger
windowTime = [-flip(exptTime(2:preWint+1)) exptTime(1:pstWint+1)]; %sec

% transform panel data
panelpsData = unwrap(deg2rad(panelpsData)); %unwrap
panelpsData = rad2deg(panelpsData); %convert back to degrees
panelStarts = mean(panelpsData(1:preWin,:,:),1,'omitnan');
panelpsData = panelpsData - panelStarts;


% threshold whether jumps are included based on fly's behavior prior
% for a jump to be included, fly must have been running above threshold
% prior to the jump being made
minRunSpeed = 3; %mm/s
minRunTime = preWin*.75; %sec

timeSpentRunning = exptTime(sum(forwardData(1:preWint,:,:)>minRunSpeed)+1);
trialsWithoutRunning = timeSpentRunning<minRunTime;
% blank out trials with insufficient running
for v = 1:nJumps
    theseBadTrials = trialsWithoutRunning(:,:,v);
    panelpsData(1:end,theseBadTrials,v) = nan;
end

%% determine if jump was corrected for

% initialize
correct_restore = [];
partial_restore = [];

for v = 1:nJumps
    % min partial correction required (deg)
    thisMinCorrect = correct_min(v);
    for c = 1:nPresent
        HD_restore = find(abs(panelpsData(preWint:end,c,v))<1,1);
        HD_partial = find(abs(panelpsData(preWint:end,c,v))<thisMinCorrect,1);
        if isempty(HD_restore)
            HD_restore=nan;
        end
        if isempty(HD_partial)
            HD_partial=nan;
        end
        % store data for this correction
        correct_restore(c,v) = HD_restore;
        partial_restore(c,v) = HD_partial;
    end
end

%% pull max values for plotting
panelpsMax = 180;
angularMax = ceil(max(abs(angularData),[],'all'));
ptime = windowTime(prePlotWint:end);

if optPlot
    % initialize
    nrow = 2;
    figure; set(gcf,'Position',[100 100 1200 400])
    tiledlayout(nrow,nJumps,'TileSpacing','compact')

    % for each jump size
    % plot FULL corrections
    for v = 1:nJumps
        nexttile
        plot(ptime,panelpsData(prePlotWint:end,isnan(correct_restore(:,v)),v),'Color', [.8 .8 .8])
        hold on
        if sum(~isnan(correct_restore(:,v)))
            plot(ptime,panelpsData(prePlotWint:end,~isnan(correct_restore(:,v)),v),'Color', '#77AC30')
        end
        title([num2str(jumps_inuse(v)) ' jump'])
        xlabel('Time (msec)')
        if v==1
            ylabel('HD (deg)')
        end
        axis tight
        ylim([-panelpsMax panelpsMax])
        yline(0)
        xline(0)
    end

    
    % plot PARTIAL corrections
    for v = 1:nJumps
        nexttile
        plot(ptime,panelpsData(prePlotWint:end,isnan(partial_restore(:,v)),v),'Color', [.8 .8 .8])
        hold on
        if sum(~isnan(correct_restore(:,v)))
            plot(ptime,panelpsData(prePlotWint:end,~isnan(partial_restore(:,v)),v),'Color', '#77AC30')
        end
        title([num2str(jumps_inuse(v)) ' jump'])
        xlabel('Time (msec)')
        if v==1
            ylabel('HD (deg)')
        end
        axis tight
        ylim([-panelpsMax panelpsMax])
        yline(0)
        xline(0)
    end
end

jumpData = 1;

end