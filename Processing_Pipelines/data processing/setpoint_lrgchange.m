% SETPOINT_LRGCHANGE - This function analyzes transitions in object position from large to small
% deviations and calculates the time between the max angular velocity and the next direction change.
% It bins angular velocity data and calculates binned means for time-to-direction-change.

% INPUTS:
%   panelps   - 3D array of panel position data where each slice represents a condition, and each
%               column represents a trial.
%   angular   - 3D array of angular velocity data corresponding to the panel position.
%   jumptrg   - 3D array indicating bar jump triggers.
%   ttime     - Time vector corresponding to the temporal resolution of the data.
%   optPlot   - Optional flag (1/0) to plot the results (1 for yes, 0 for no).

% OUTPUTS:
%   binned_means_obj - Binned means for time-to-direction-change, based on angular velocity.
%   bins             - Structure containing angular velocity bin centers for plotting.

% Created: N/A by MC
% Updated: N/A

% The function:
% - Identifies transitions where the object position moves from large (>40) to small (<40) deviations.
% - Extracts and bins angular velocity values before the direction change and calculates the time 
%   between the max angular velocity and the zero-crossing.
% - Binned means are calculated for time-to-direction-change for each condition.
% - Optionally plots the panel position, angular velocity, and binned means for each condition.
%
function [binned_means_obj, bins] = setpoint_lrgchange(panelps, angular, jumptrg, ttime, optPlot)
%% initialize
% fetch experiment conditions
nTrial = size(panelps,2);
nCond = size(panelps,3);

% set analysis parameters
pre_window = 0.25; % Pre-transition window in seconds
post_window = 0.75; % Post-transition window in seconds
pre_idx_window = fetchTimeIdx(ttime, pre_window);
post_idx_window = fetchTimeIdx(ttime, post_window);

% Define bin edges for angular velocity and calculate bin centers
ang_vel_bin_edges = [50:100:450, 500]; % Modify this to suitable bin edges for your data
bins.ang_vel_center = ang_vel_bin_edges(1:end-1) + diff(ang_vel_bin_edges)/2; % Store bin centers

% Preallocate binned means for each condition
binned_means_obj = nan(length(bins.ang_vel_center), nCond); % Preallocate binned means

min_valid_idx = 10; % Minimum number of valid points required
min_bin_size = 3;   % Minimum number of points in each bin

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
    angmax_r = [];
    angmax_l = [];

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

        % Find the first time the angular velocity crosses zero
        zero_cross_idx_ang = find(thiswin_angular(pre_idx_window+1:end) .* thiswin_angular(pre_idx_window:end-1) <= 0, 1, 'first');
        if isempty(zero_cross_idx_ang)
            continue;
        end

        % Limit the angular velocity data to before the zero crossing
        angular_before_zero = thiswin_angular(1:(pre_idx_window + zero_cross_idx_ang));

        % Find the max angular velocity before the zero crossing
        [max_ang_val, max_ang_idx] = max(abs(angular_before_zero));

        if ~isempty(max_ang_val)
            winpanelps_r(:,ri) = thiswin_panelps;
            winangular_r(:,ri) = thiswin_angular;
            angmax_r(ri) = max_ang_val;
            time_to_zero_r(ri) = dc_t(pre_idx_window + zero_cross_idx_ang) - dc_t(max_ang_idx);
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

        % Find the first time the angular velocity crosses zero
        zero_cross_idx_ang = find(thiswin_angular(pre_idx_window+1:end) .* thiswin_angular(pre_idx_window:end-1) <= 0, 1, 'first');
        if isempty(zero_cross_idx_ang)
            continue;
        end

        % Limit the angular velocity data to before the zero crossing
        angular_before_zero = thiswin_angular(1:(pre_idx_window + zero_cross_idx_ang));

        % Find the max angular velocity before the zero crossing
        [max_ang_val, max_ang_idx] = max(abs(angular_before_zero));

        if ~isempty(max_ang_val)
            winpanelps_l(:,li) = thiswin_panelps;
            winangular_l(:,li) = thiswin_angular;
            angmax_l(li) = max_ang_val;
            time_to_zero_l(li) = dc_t(pre_idx_window + zero_cross_idx_ang) - dc_t(max_ang_idx);
            li = li+1;
        end
    end

    % Combine right and left crossings
    combined_angmax = [angmax_r, angmax_l]; % combined max angular velocity values
    combined_time_to_zero = [time_to_zero_r, time_to_zero_l]; % combined times

    %% Bin data by angular velocity and calculate binned means (if enough points)
    valid_ang_idx = combined_angmax >= 50; % Only use angular velocity max values above 50 (adjust threshold if necessary)
    nvalid = sum(valid_ang_idx);

    if nvalid > min_valid_idx
        % Bin data using the specified angular velocity bin edges
        [bin_counts, ~, bin_idx] = histcounts(combined_angmax, ang_vel_bin_edges);
        bin_means = nan(1, length(ang_vel_bin_edges)-1);

        for b = 1:length(ang_vel_bin_edges)-1
            bin_data = combined_time_to_zero(bin_idx == b);
            if length(bin_data) >= min_bin_size
                bin_means(b) = mean(bin_data); % Calculate mean for valid bins
            end
        end

        % Store the binned means for output
        binned_means_obj(:,c) = bin_means;
    else
        % If not enough valid points, set binned means to NaN
        binned_means_obj(:,c) = NaN;
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

        % Scatter plot of max angular velocity vs time-to-direction-change with binned means
        nexttile; hold on
        if ~isempty(angmax_r) && ~isempty(time_to_zero_r)
            scatter(angmax_r, time_to_zero_r, '.', 'MarkerEdgeColor', "#0072BD");
        end
        if ~isempty(angmax_l) && ~isempty(time_to_zero_l)
            scatter(angmax_l, time_to_zero_l, '.', 'MarkerEdgeColor', [0.2 0.2 0.2]);
        end
        if nvalid > min_valid_idx
            plot(bins.ang_vel_center, bin_means, 'o-', 'Color','k','MarkerEdgeColor', 'k', 'LineWidth', 1.5); % Plot binned means
        end
        legend(num2str(nvalid))
        set(gca, 'XScale', 'log')
        ylim([0 (pre_window+post_window)*1000]); xlim([50 500]); ylabel('Time 2 Change Dir (ms)'); xlabel('Max Angular (deg/s)');
    end
end

end
