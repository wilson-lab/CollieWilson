% fixationFinder.m
% Function used to pull only pursuit epochs by first thresholding with a
% Schmitt trigger and then by thresholding based on epoch duration
%
% INPUTS:
% panelps - panel position array containing all data
% forward - forward velocity array containing all data
% expttime - time array for 1 trial
% runThresh - minimum forward velocity for "running"
% 
% OUTPUTS:
% fixation - only fixation panelps and indices
%
% CREATED: 08/30/2024 MC
%

function fixation = fixationFinder(panelps,forward,expttime)
%% initialize
durTrial = size(panelps,1);
nTrial = size(panelps,2);
nCond = size(panelps,3);

% set min fixation time
minFixateTime = 4; %sec
minFixationIdx = fetchTimeIdx(expttime,minFixateTime);

% fetch settings
settings = processSettings();

% set threshold settings
fixateHighThresh = 10; % deg
fixateLowThresh = 130; % deg

%% threshold for when the fly was running

% find run indices using schmitt trigger
run_idx = zeros(durTrial,nTrial,nCond);
for c = 1:nCond
    run_idx(:,:,c) = schmittTrigger(forward(:,:,c),settings.runThreshB,0.1);
end

%% threshold for when the fly was fixating the target
% for each condition
fix_idx = zeros(durTrial,nTrial,nCond);
for c = 1:nCond
    fix_idx(:,:,c) = schmittTrigger(-abs(panelps(:,:,c)),-fixateHighThresh,-fixateLowThresh);
end
panelps_f = panelps;
panelps_f(~fix_idx) = nan;
forward_f = forward;
forward_f(~fix_idx) = nan;

% test
% clf
% for t = 1:nTrial
%     subplot(nTrial,1,t); hold on;
%     plot(expttime,panelps(:,t,7),'k')
%     plot(expttime,panelps_f(:,t,7),'r')
% end

%% threshold for pursuit epochs based on run duration
% for each condition
fix_idx2 = fix_idx;
panelps_fix = panelps_f;
forward_fix = forward_f;
for c = 1:nCond
    % for each trial
    for t = 1:nTrial
        % find pursuit start/stops for this trial
        fixation_startstop = diff(fix_idx(:,t,c));
        fixation_start = find(fixation_startstop==1);
        fixation_stop = find(fixation_startstop==-1);

        % if fly starts trial with pursuit (no first start)
        if fix_idx(1,t,c)
            fixation_start = [1 ; fixation_start];
        end
        % if fly ends trial with pursuit (no last stop)
        if fix_idx(end,t,c)
            fixation_stop = [fixation_stop ; durTrial];
        end

        % for each stop/start pair
        pursuit_duration = fixation_stop-fixation_start;
        for p = 1:length(pursuit_duration)
            % if this epoch is shorter than minimum
            if pursuit_duration(p)<minFixationIdx
                % blank out epochs that are too short
                fix_idx2(fixation_start(p):fixation_stop(p),t,c) = 0;
                panelps_fix(fixation_start(p):fixation_stop(p),t,c) = NaN;
                forward_fix(fixation_start(p):fixation_stop(p),t,c) = NaN;
            end
        end
    end
end
% test
clf
for t = 1:nTrial
    c = 4;
    subplot(nTrial,1,t); hold on;
    plot(expttime,panelps(:,t,c),'k')
    plot(expttime,panelps_f(:,t,c),'b')
    plot(expttime,panelps_fix(:,t,c),'r')
end

%% store for output
% store ALL data where fly was fixating
fixation.panelps_all = panelps_fix;
fixation.forward_all = forward_fix;
fixation.idx_all = logical(fix_idx2);

% store data where fly was fixating AND running
% fetch running timepoints
panelps_fix_r = panelps_fix;
panelps_fix_r(~run_idx) = nan;
forward_fix_r = forward_fix;
forward_fix_r(~run_idx) = nan;
fix_idxR = fix_idx2;
fix_idxR(~run_idx) = 0;

fixation.panelps_run = panelps_fix_r;
fixation.forward_run = forward_fix_r;
fixation.idx_run = logical(fix_idxR);

end