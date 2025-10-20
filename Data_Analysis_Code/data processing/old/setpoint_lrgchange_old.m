function [dc_output, dc_t, objbin_means_per_condition, objbin_center] = setpoint_lrgchange(panelps, angular, jumptrg, ttime, optPlot)
%% initialize
% fetch experiment conditions
nTrial = size(panelps,2);
nCond = size(panelps,3);

% set analysis parameters
pre_window = 0.5; % Pre-transition window in seconds
post_window = 0.75; % Post-transition window in seconds
pre_idx_window = fetchTimeIdx(ttime, pre_window);
post_idx_window = fetchTimeIdx(ttime, post_window);

% define object position bins
objbin_edges = [40:20:80 120 180];
objbin_center = (objbin_edges(1:end-1) + objbin_edges(2:end)) / 2; % Centers of the object position bins

minBinMembers = 5;

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
% initialize
dc_output = [];
objbin_means_per_condition = []; % To store bin means for each condition
window_size = length(-pre_idx_window:post_idx_window);
dc_t = (ttime(1:window_size) - pre_window) * 1000;  % convert to ms
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
    time_to_zero_r = [];
    time_to_zero_l = [];
    objmax_r = [];
    objmax_l = [];

    % estimate HD bias
    biasHD = mean(panelps(:,:,c), 'all', 'omitnan');

    % fetch data
    thispanelps = reshape(panelps(:,:,c) - biasHD, [], 1);
    thisangular = reshape(angular(:,:,c), [], 1);

    % Right crossings (from >40 to <40)
    transitions_right = diff(thispanelps > 40); % Find the points where it transitions from >40 to <40
    cross_right = find(transitions_right == -1); % Panel position goes from >40 to <40

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

        % Find the max object value before the transition
        [max_val, max_idx] = max(thiswin_panelps(1:pre_idx_window));
        
        % Find the first time the object crosses 0 after the max
        zero_cross_idx = find(thiswin_panelps(pre_idx_window+1:end) <= 0, 1, 'first');
        if ~isempty(zero_cross_idx)
            winpanelps_r(:,ri) = thiswin_panelps;
            winangular_r(:,ri) = thiswin_angular;
            objmax_r(ri) = max_val;
            time_to_zero_r(ri) = dc_t(pre_idx_window + zero_cross_idx) - dc_t(max_idx);
            ri = ri+1;
        end
    end

    % Left crossings (from <-40 to >-40)
    transitions_left = diff(thispanelps < -40); % Find the points where it transitions from <-40 to >-40
    cross_left = find(transitions_left == -1); % Panel position goes from <-40 to >-40

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

        % Find the max object value before the transition
        [max_val, max_idx] = min(thiswin_panelps(1:pre_idx_window));

        % Find the first time the object crosses 0 after the max
        zero_cross_idx = find(thiswin_panelps(pre_idx_window+1:end) >= 0, 1, 'first');
        if ~isempty(zero_cross_idx)
            winpanelps_l(:,li) = thiswin_panelps;
            winangular_l(:,li) = thiswin_angular;
            objmax_l(li) = -max_val;
            time_to_zero_l(li) = dc_t(pre_idx_window + zero_cross_idx) - dc_t(max_idx);
            li = li+1;
        end
    end

    % Combine right and left crossings
    combined_objmax = [objmax_r, objmax_l]; % combined max object values
    combined_time_to_zero = [time_to_zero_r, time_to_zero_l]; % combined times

    % Discretize the max object before crossing into bins
    [objbin_indices, ~] = discretize(combined_objmax, objbin_edges);
    objbin_means = nan(length(objbin_edges)-1, 1); % Store the mean time-to-zero for each object position bin
    for bin = 1:length(objbin_means)
        bin_members = combined_time_to_zero(objbin_indices == bin);
        if length(bin_members) >= minBinMembers % Only include bins with at least 3 members
            objbin_means(bin) = mean(bin_members, 'omitnan');
        end
    end
    % check for gaps
    minSequenceLength = 2;
    objbin_means_clean = removeIsolatedPoints(objbin_means, minSequenceLength);

    objbin_means_per_condition(:,c) = objbin_means_clean;

    % Store the results in the output
    time_to_zero_rl = [time_to_zero_r time_to_zero_l];
    if ~isempty(time_to_zero_rl)
        dc_output(:,c) = mean([time_to_zero_r time_to_zero_l], 2, 'omitnan');
    else
        dc_output(:,c) = nan;
    end

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

    % Scatter plot of object max vs time-to-zero (right)
    nexttile; hold on
    if ~isempty(objmax_r) && ~isempty(time_to_zero_r)
        scatter(objmax_r, time_to_zero_r, '.', 'MarkerEdgeColor', "#0072BD");
    end
    if ~isempty(objmax_l) && ~isempty(time_to_zero_l)
        scatter(objmax_l, time_to_zero_l, '.', 'MarkerEdgeColor', [0.2 0.2 0.2]);
    end
    if ~isempty(objbin_means)
        plot(objbin_center, objbin_means, 'Color', "r", 'LineWidth', 1.5);
    end
    ylim([0 (pre_window+post_window)*1000]); xlim([0 max(objbin_edges)]); ylabel('Time 2 Cross Zero (ms)'); xlabel('Obj Max (deg)');
end

end

end
