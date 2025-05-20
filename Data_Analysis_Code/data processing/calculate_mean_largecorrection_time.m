% calculate_mean_largecorrection_time
% This function analyzes the time it takes for the fly to change direction 
% in response to a motion pulse by calculating the average correction time 
% for both angular velocity and object position. It accounts for rightward 
% and leftward movements separately and provides optional plotting of the results.
%
% INPUTS:
%   panelps      - Panel positions (degrees) representing the target position
%   angular      - Angular velocity (degrees/second)
%   jumptrg      - Trigger for jumps or transitions in behavior
%   ttime        - Time vector (seconds)
%   optPlot      - Flag to determine whether to plot the results (1 for plotting, 0 to skip)
%
% OUTPUTS:
%   dc_output_vel_to_dir_change - Average time from maximum angular velocity to direction change
%   dc_output_pos_to_dir_change  - Average time from maximum object position to direction change
%   dc_t                          - Time vector for the windows used in analysis (in milliseconds)
%
% CREATED: [Date] MC
%
function [dc_output_vel_to_dir_change, dc_output_pos_to_dir_change, dc_t] = calculate_mean_largecorrection_time(panelps, angular, jumptrg, ttime, optPlot)
%% initialize
% fetch experiment conditions
nTrial = size(panelps,2);
nCond = size(panelps,3);

% set analysis parameters
pre_window = 0.25; % Pre-transition window in seconds
post_window = 0.75; % Post-transition window in seconds
pre_idx_window = fetchTimeIdx(ttime, pre_window);
post_idx_window = fetchTimeIdx(ttime, post_window);

% Preallocate time and output variables
dc_output_vel_to_dir_change = nan(nCond,1);  % Stores avg correction time (max angular vel to angular dir change)
dc_output_pos_to_dir_change = nan(nCond,1);  % Stores avg correction time (max obj position to angular dir change)
window_size = length(-pre_idx_window:post_idx_window);
dc_t = (ttime(1:window_size) - pre_window) * 1000;  % convert to ms

%% optional pre-process panel data by removing bar jumps
% omit timepoints right after a bar jump
tomit = 2; % s
iomit = fetchTimeIdx(ttime, tomit); % idx

for c = 1:nCond
    for t = 1:nTrial
        jumpidx = find(diff(jumptrg(:,t,c)) > 0);
        for j = 1:length(jumpidx)
            panelps(jumpidx(j):jumpidx(j) + iomit, t, c) = nan;
        end
    end
end

%% find crossings (triggering on transitions from large to small)
if optPlot
    % initialize plot
    figure; set(gcf,'Position',[100 100 1000 900])
    tiledlayout(nCond,3,'TileSpacing','compact') % Updated to 3 columns for plotting
end

% for each condition
for c = 1:nCond
    % initialize
    winpanelps_r = [];
    winpanelps_l = [];
    winangular_r = [];
    winangular_l = [];
    time_to_dir_change_r_vel = [];
    time_to_dir_change_l_vel = [];
    time_to_dir_change_r_pos = [];
    time_to_dir_change_l_pos = [];
    max_ang_vel_r = [];
    max_ang_vel_l = [];
    max_obj_pos_r = [];
    max_obj_pos_l = [];

    % estimate HD bias
    biasHD = mean(panelps(:,:,c), 'all', 'omitnan');

    % fetch data
    thispanelps = reshape(panelps(:,:,c) - biasHD, [], 1);
    thisangular = reshape(angular(:,:,c), [], 1);

    % Right crossings (from >40 to <40)
    transitions_right = diff(thispanelps > 40); % Find the points where it transitions from >40 to <40
    cross_right = find(transitions_right == -1); % Panel position goes from >40 to <40
    if ~isempty(cross_right)
        cross_right = cross_right([true; diff(cross_right) >= 200]); % Ensure indices are at least 200 apart
    end

    % Right crossings
    ri = 1;
    for i = 1:length(cross_right)
        idx = cross_right(i);
        % Fetch panel and angular data before and after the crossing
        thiswin_panelps = thispanelps(idx-pre_idx_window:idx+post_idx_window);
        thiswin_angular = thisangular(idx-pre_idx_window:idx+post_idx_window);

        % Omit if there are NaNs in the pre-transition window or if a large circular jump occurs
        if any(isnan(thiswin_panelps(1:pre_idx_window))) || any(abs(diff(thiswin_panelps)) > 100)
            continue;
        end

        % Find the point where the angular velocity crosses zero (direction change)
        zero_cross_idx = find(thiswin_angular(pre_idx_window+1:end) <= 0, 1, 'first');
        
        if ~isempty(zero_cross_idx)
            % Only consider the angular velocity up to this point
            max_vel_window = thiswin_angular(1:(pre_idx_window + zero_cross_idx - 1));

            % Find the maximum angular velocity before the crossing
            [max_vel, max_vel_idx] = max(max_vel_window);
            % Find the maximum panel position before the crossing
            [max_obj_pos, max_obj_pos_idx] = max(thiswin_panelps(1:pre_idx_window));

            winpanelps_r(:,ri) = thiswin_panelps;
            winangular_r(:,ri) = thiswin_angular;
            max_ang_vel_r(ri) = max_vel;
            max_obj_pos_r(ri) = max_obj_pos;
            time_to_dir_change_r_vel(ri) = dc_t(pre_idx_window + zero_cross_idx) - dc_t(max_vel_idx);  % Time from max angular velocity to direction change
            time_to_dir_change_r_pos(ri) = dc_t(pre_idx_window + zero_cross_idx) - dc_t(max_obj_pos_idx);  % Time from max object position to direction change
            ri = ri+1;
        end
    end

    % Left crossings (from <-40 to >-40)
    transitions_left = diff(thispanelps < -40); % Find the points where it transitions from <-40 to >-40
    cross_left = find(transitions_left == -1); % Panel position goes from <-40 to >-40
    if ~isempty(cross_left)
        cross_left = cross_left([true; diff(cross_left) >= 200]); % Ensure indices are at least 200 apart
    end

    % Left crossings
    li = 1;
    for i = 1:length(cross_left)
        idx = cross_left(i);
        % Fetch panel and angular data before and after the crossing
        thiswin_panelps = thispanelps(idx-pre_idx_window:idx+post_idx_window);
        thiswin_angular = thisangular(idx-pre_idx_window:idx+post_idx_window);

        % Omit if there are NaNs in the pre-transition window or if a large circular jump occurs
        if any(isnan(thiswin_panelps(1:pre_idx_window))) || any(abs(diff(thiswin_panelps)) > 100)
            continue;
        end

        % Find the point where the angular velocity crosses zero (direction change)
        zero_cross_idx = find(thiswin_angular(pre_idx_window+1:end) >= 0, 1, 'first');
        
        if ~isempty(zero_cross_idx)
            % Only consider the angular velocity up to this point
            max_vel_window = thiswin_angular(1:(pre_idx_window + zero_cross_idx - 1));

            % Find the maximum angular velocity before the crossing
            [max_vel, max_vel_idx] = min(max_vel_window);
            % Find the maximum panel position before the crossing
            [max_obj_pos, max_obj_pos_idx] = min(thiswin_panelps(1:pre_idx_window));

            winpanelps_l(:,li) = thiswin_panelps;
            winangular_l(:,li) = thiswin_angular;
            max_ang_vel_l(li) = -max_vel;
            max_obj_pos_l(li) = -max_obj_pos;
            time_to_dir_change_l_vel(li) = dc_t(pre_idx_window + zero_cross_idx) - dc_t(max_vel_idx);  % Time from max angular velocity to direction change
            time_to_dir_change_l_pos(li) = dc_t(pre_idx_window + zero_cross_idx) - dc_t(max_obj_pos_idx);  % Time from max object position to direction change
            li = li+1;
        end
    end

    % Combine right and left crossings
    combined_time_to_dir_change_vel = [time_to_dir_change_r_vel, time_to_dir_change_l_vel]; % combined times for velocity
    combined_time_to_dir_change_pos = [time_to_dir_change_r_pos, time_to_dir_change_l_pos]; % combined times for position

    %% Calculate average correction times for this condition
    if ~isempty(combined_time_to_dir_change_vel)
        dc_output_vel_to_dir_change(c) = mean(combined_time_to_dir_change_vel, 'omitnan'); % Store average correction time for velocity
    else
        dc_output_vel_to_dir_change(c) = nan; % No valid data for this condition
    end
    
    if ~isempty(combined_time_to_dir_change_pos)
        dc_output_pos_to_dir_change(c) = mean(combined_time_to_dir_change_pos, 'omitnan'); % Store average correction time for position
    else
        dc_output_pos_to_dir_change(c) = nan; % No valid data for this condition
    end

    %% Plotting
    if optPlot
        % Plot object position traces
        nexttile; hold on
        if ~isempty(winpanelps_r)
            plot(dc_t, winpanelps_r, 'Color', "#0072BD");
        end
        if ~isempty(winpanelps_l)
            plot(dc_t, winpanelps_l, 'Color', [0.2 0.2 0.2]);
        end
        axis tight; ylim([-180 180]); yline(0); xline(0); ylabel('Obj (deg)'); xlabel('Time (ms)');

        % Plot angular velocity traces
        nexttile; hold on
        if ~isempty(winangular_r)
            plot(dc_t, winangular_r, 'Color', "#0072BD");
        end
        if ~isempty(winangular_l)
            plot(dc_t, winangular_l, 'Color', [0.2 0.2 0.2]);
        end
        axis tight; ylim([-500 500]); yline(0); xline(0); ylabel('Ang (deg/s)'); xlabel('Time (ms)');

        % Scatter plot of max angular velocity vs time-to-direction-change
        nexttile; hold on
        if ~isempty(max_ang_vel_r) && ~isempty(time_to_dir_change_r_vel)
            scatter(max_ang_vel_r, time_to_dir_change_r_vel, '.', 'MarkerEdgeColor', "#0072BD");
        end
        if ~isempty(max_ang_vel_l) && ~isempty(time_to_dir_change_l_vel)
            scatter(max_ang_vel_l, time_to_dir_change_l_vel, '.', 'MarkerEdgeColor', [0.2 0.2 0.2]);
        end
        ylim([0 (pre_window+post_window)*1000]); xlim([00 500]); ylabel('Time 2 Dir Change (ms)'); xlabel('Max Angular Vel (deg/s)');
    end
end

end
