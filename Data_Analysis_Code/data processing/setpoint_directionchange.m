% setpoint_directionchange
% This function analyzes the intervals between direction changes in angular velocity.
% It identifies zero-crossings (sign changes) in the angular velocity input, calculates
% the time intervals between these changes, and filters out intervals exceeding a 
% specified maximum duration. Optionally, it plots histograms of the intervals for 
% each condition.
%
% INPUTS:
% angular_vel  - angular velocity input (matrix: time x trials x conditions)
% ttime        - time vector (in seconds)
% optPlot      - 1 to plot the results, 0 to skip plotting
%
% OUTPUT:
% interval_out - structure containing the intervals between direction changes for each condition
%
% CREATED: 10/09/24 - MC
%
function [interval_out] = setpoint_directionchange(angular_vel, ttime, optPlot)
%% initialize
% dataset info
nCond = size(angular_vel,3);
nTrial = size(angular_vel,2);

% set maximum interval duration (s)
maxIntervalDuration = 5; 

% store intervals between direction changes
storeIntervals = [];

%% analyze angular velocity
for c = 1:nCond
    tempIntervals = [];
    
    % for each trial
    for t = 1:nTrial
        % fetch data
        thisAngularVel = angular_vel(:,t,c);
        
        % remove nans
        thisAngularVel(isnan(thisAngularVel)) = 0;
        
        % find direction changes (zero-crossings in angular velocity)
        dir_changes = find(diff(sign(thisAngularVel)) ~= 0);
        
        % calculate time intervals between direction changes
        change_times = ttime(dir_changes);
        intervals = diff(change_times);
        
        % filter out intervals greater than maxIntervalDuration
        valid_idx = intervals <= maxIntervalDuration;
        valid_intervals = intervals(valid_idx);
        valid_changes = dir_changes(valid_idx); % Corresponding valid change points
        
        % Ensure each valid interval corresponds to a real change in direction
        for i = 1:length(valid_intervals)
            if valid_intervals(i) <= maxIntervalDuration
                tempIntervals{t}(i) = valid_intervals(i); 
            end
        end
    end
    
    % store intervals for each condition
    storeIntervals{c} = tempIntervals;
end

% store for output
interval_out.intervals = storeIntervals;

%% (optional) plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 800 600])
    tiledlayout(nCond,1,'TileSpacing','compact') % 1 column per condition
    
    for c = 1:nCond
        % Plot histogram of intervals for each condition
        nexttile
        all_intervals = [storeIntervals{c}{:}]; % concatenate all trials for this condition
        histogram(all_intervals, 'BinEdges', 0:0.1:max(all_intervals)); % adjust bins as needed
        xlabel('Interval (s)'); ylabel('Count'); title(['Condition ' num2str(c) ' Direction Change Intervals'])
    end
end
end