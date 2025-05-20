% optomotor_v_velocity
% This analysis function generates a summary plot of panel position versus 
% rotational sweeps during an optomotor assay. It plots leftward and rightward 
% sweeps separately, calculating the mean angular velocity for each direction 
% and adjusting for any baseline biases.
%
% INPUTS:
%   panelps   - Downsampled panel positions (degrees)
%   forward   - Downsampled forward velocities (mm/s)
%   angular   - Downsampled angular velocities (degrees/second)
%   ttime     - Time vector (in seconds)
%   runSelect - Behavioral parameters (0 for quiescent, >0 for running, -1 for all)
%   optPlot   - 1 to generate plots, 0 to skip plotting
%
% OUTPUTS:
%   mean_angularR - Mean angular velocities for rightward sweeps
%   mean_angularL - Mean angular velocities for leftward sweeps
%   panelWin      - Panel sweep positions for visualization
%
% ORIGINAL: 04/11/2024 - MC
%           08/01/2024 - MC (adjusted buffer window determination, changed behavior classifications)
%
function [mean_angularR,mean_angularL,panelWin] = optomotor_v_velocity(panelps,forward,angular,ttime,runSelect,optPlot)
%% initialize
nTrial = size(panelps,2);

% set analysis window
preT = 1; % s, pre stim
stmT = 0.5; %s, stim duration
pstT = 3; % s, post stim
[~,winPre] = min(abs(ttime - preT)); %find nearest index
[~,winStm] = min(abs(ttime - stmT)); %find nearest index
[~,winPst] = min(abs(ttime - pstT)); %find nearest index
sweepLength = winPre + winStm + winPst;

% for each pulse trial, set the percentage of time behavior must have been
% sufficient for said trial to be included in the average
minInclusion = 0.75;
% set number of pulse trials for average to be taken
minTrials = 5;

color_trials = [0.6 0.6 0.6];
color_mean = "#0072BD";

%% fetch behavioral index depending on run select
if runSelect==0 %quiescent only
    moveThresh = 0.25; %anything above/below considered moving
    runIdx = ~(forward<-moveThresh | forward>moveThresh);
elseif runSelect>0 %running only
    runIdx = schmittTrigger(forward,runSelect,0.1);
else %all
    runIdx = ones(size(forward));
end

%% fetch sweep data
% initialize
data_angularR=[];
data_angularL=[];
rc = 1; %counter
lc = 1; %counter

% for each trial
for t = 1:nTrial
    % find sweep starts
    sweepStartStop = diff(~isnan(panelps(:,t)));
    % set sweep start/stop indices based on analysis window
    sweepStart = find(sweepStartStop>0)-winPre;
    sweepStops = sweepStart+sweepLength;
    nSweep = size(sweepStart,1);

    % for each rightward sweep (odd), fetch data
    for s = 1:2:nSweep
        data_angularR(:,rc) = angular(sweepStart(s):sweepStops(s),t);
        runIdxR(:,rc) = runIdx(sweepStart(s):sweepStops(s),t);
        rc = rc+1; %update counter
    end
    % for each leftward sweep (even), fetch data
    for s = 2:2:nSweep
        data_angularL(:,lc) = angular(sweepStart(s):sweepStops(s),t);
        runIdxL(:,lc) = runIdx(sweepStart(s):sweepStops(s),t);
        lc=lc+1; %update counter
    end
end
optoDur = size(data_angularL,1);
t_opto = ttime(1:optoDur)*1000;
tmax=max(t_opto);

% fetch panel sweep examples
panelWin=[];
panelWin(:,1) = panelps(sweepStart(1):sweepStart(1)+sweepLength,t);
panelWin(:,2) = panelps(sweepStart(2):sweepStart(2)+sweepLength,t);

%% determine which trials met behavioral requirements (if any) for this run
minIdx = round(optoDur*minInclusion);
goodTrialsR = sum(runIdxR,1)>minIdx;
goodTrialsL = sum(runIdxL,1)>minIdx;

%% calculate mean

% if sufficient, calculate mean right
if sum(goodTrialsR)>=minTrials
    mean_angularR = mean(data_angularR,2,'omitnan');
else
    mean_angularR = nan(optoDur,1);
end
% if sufficient, calculate mean left
if sum(goodTrialsL)>=minTrials
    mean_angularL = mean(data_angularL,2,'omitnan');
else
    mean_angularL = nan(optoDur,1);
end
% adjust means for fly bias
adj_r = mean(mean_angularR(1:floor(winPre/2)));
adj_l = mean(mean_angularL(1:floor(winPre/2)));

mean_angularR = mean_angularR-adj_r;
mean_angularL = mean_angularL-adj_l;

%% optional: plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 600 700])
    tiledlayout(2,2,'TileSpacing','compact')
    ang_range = [-150 150];

    % plot sweep position
    nexttile
    plot(t_opto,panelWin(:,1),'k','Linewidth',1.5)
    xlim([0 tmax])
    nexttile
    plot(t_opto,panelWin(:,2),'k','Linewidth',1.5)
    xlim([0 tmax])
    
    % plot angular velocity
    nexttile; hold on
    plot(t_opto,data_angularR,'-','Color',color_trials)
    plot(t_opto,mean_angularR,'-','Linewidth',1.5,'Color',color_mean)
    axis tight
    yline(0)
    xlim([0 tmax])
    ylim(ang_range)
    ylabel('Angular Velocity (deg/s)')
    title('Rightward')

    nexttile; hold on
    plot(t_opto,data_angularL,'-','Color',color_trials)
    plot(t_opto,mean_angularL,'-','Linewidth',1.5,'Color',color_mean)
    axis tight
    yline(0)
    xlim([0 tmax])
    ylim(ang_range)
    xlabel('Time (msec.)')
    title('Leftward')
end

end

