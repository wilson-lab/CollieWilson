function [dc_output, dc_t, angbin_means_per_condition, angbin_center, posbin_means_per_condition, posbin_center] = setpoint_dirchange(panelps, angular, jumptrg, ttime, optPlot)
%% initialize
% fetch experiment conditions
nTrial = size(panelps,2);
nCond = size(panelps,3);

% set analysis parameters
cross_window = 0.3; % s
idx_window = fetchTimeIdx(ttime, cross_window);

% define angular velocity bins
angbin_edges = [0:40:200, 300];
angbin_center = (angbin_edges(1:end-1) + angbin_edges(2:end)) / 2; % Centers of the angular bins

% define object position bins
posbin_edges = [0:15:45 75];
posbin_center = (posbin_edges(1:end-1) + posbin_edges(2:end)) / 2; % Centers of the object position bins

minBinMembers = 10;

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

%% find crossings
% initialize
dc_output = [];
angbin_means_per_condition = []; % To store bin means for each condition
posbin_means_per_condition = []; % To store bin means for each condition
window_size = length(1-idx_window:1+idx_window);
dc_t = (ttime(1:window_size) - cross_window) * 1000;  % convert to ms
if optPlot
    % initialize plot
    figure; set(gcf,'Position',[100 100 1800 900])
    tiledlayout(nCond,6,'TileSpacing','compact')
end

% for each condition
for c = 1:nCond
    % initialize
    winpanelps_r = [];
    winpanelps_l = [];
    winangular_r = [];
    winangular_l = [];
    % estimate HD bias
    biasHD = mean(panelps(:,:,c), 'all', 'omitnan');

    % fetch data
    thispanelps = reshape(panelps(:,:,c) - biasHD, [], 1);
    thisangular = reshape(angular(:,:,c), [], 1);

    % calculate the difference in sign of the position values
    sign_changes = diff(sign(thispanelps));
    % find indices where sign changes occur (this indicates a zero-crossing)
    cross_right = find(sign_changes > 0);
    cross_left = find(sign_changes < 0);

    % Right crossings
    angzero_cross_times_r = [];
    angmax_before_cross_r = [];
    objmax_before_cross_r = [];
    ri = 1;
    for i = 1:length(cross_right)
        idx = cross_right(i);
        % Fetch angular data before and after the crossing
        thiswin_angular = thisangular(idx-idx_window:idx+idx_window);
        thiswin_panelps = thispanelps(idx-idx_window:idx+idx_window);

        % Check if the difference between any adjacent panel position values exceeds 10
        if any(abs(diff(thiswin_panelps)) > 10)
            continue; % Skip this window if any adjacent differences exceed 10
        end

        % Find zero-crossing time after the panel crossing
        sign_change_window = diff(sign(thiswin_angular(idx_window+1:end))); % search after the panel crossing
        zero_cross_idx = find(sign_change_window > 0, 1, 'first'); % first zero-crossing after the panel crossing
        if ~isempty(zero_cross_idx)
            % Store window
            winpanelps_r(:,ri) = thiswin_panelps;
            winangular_r(:,ri) = thiswin_angular;
            % Time of zero crossing
            angzero_cross_times_r(ri) = dc_t(idx_window+1+zero_cross_idx);
            % Find max angular value before the angular zero-crossing
            angmax_before_cross_r(ri) = -min(thiswin_angular(1:idx_window+zero_cross_idx), [], 'omitnan');
            % Find max object value before the angular zero-crossing
            objmax_before_cross_r(ri) = -min(thiswin_panelps(1:idx_window+zero_cross_idx), [], 'omitnan');
            ri = ri+1;
        end
    end


    % Left crossings
    angzero_cross_times_l = [];
    angmax_before_cross_l = [];
    objmax_before_cross_l = [];
    li = 1;
    for i = 1:length(cross_left)
        idx = cross_left(i);
        % Fetch angular data before and after the crossing
        thiswin_angular = thisangular(idx-idx_window:idx+idx_window);
        thiswin_panelps = thispanelps(idx-idx_window:idx+idx_window);

        % Check if the difference between any adjacent panel position values exceeds 10
        if any(abs(diff(thiswin_panelps)) > 10)
            continue; % Skip this window if any adjacent differences exceed 10
        end

        % Find zero-crossing time after the panel crossing
        sign_change_window = diff(sign(thiswin_angular(idx_window+1:end))); % search after the panel crossing
        zero_cross_idx = find(sign_change_window < 0, 1, 'first'); % first zero-crossing after the panel crossing
        if ~isempty(zero_cross_idx)
            % Store window
            winpanelps_l(:,li) = thiswin_panelps;
            winangular_l(:,li) = thiswin_angular;
            % Time of zero crossing
            angzero_cross_times_l(li) = dc_t(idx_window+1+zero_cross_idx);
            % Find max angular value before the angular zero-crossing
            angmax_before_cross_l(li) = max(thiswin_angular(1:idx_window+zero_cross_idx), [], 'omitnan');
            % Find max object value before the angular zero-crossing
            objmax_before_cross_l(li) = max(thiswin_panelps(1:idx_window+zero_cross_idx), [], 'omitnan');
            li = li+1;
        end
    end


    % Combine right and left crossings
    combined_zero_cross_times = [angzero_cross_times_r, angzero_cross_times_l]; % combined times
    combined_angmax_before_cross = [angmax_before_cross_r, angmax_before_cross_l]; % combined max values
    combined_posmax_before_cross = [objmax_before_cross_r, objmax_before_cross_l]; % combined max values

    % Discretize the max_before_cross into bins (angular)
    [angbin_indices, ~] = discretize(combined_angmax_before_cross, angbin_edges);
    angbin_means = nan(length(angbin_edges)-1, 1); % Store the mean zero-crossing times for each angular bin
    for bin = 1:length(angbin_means)
        bin_members = combined_zero_cross_times(angbin_indices == bin);
        if length(bin_members) >= minBinMembers % Only include bins with at least X members
            angbin_means(bin) = mean(bin_members, 'omitnan');
        end
    end
    % check for gaps
    minSequenceLength = 3;
    angbin_means_clean = removeIsolatedPoints(angbin_means, minSequenceLength);

    angbin_means_per_condition(:,c) = angbin_means_clean;

    % Discretize the max_before_cross into bins (object position)
    [posbin_indices, ~] = discretize(combined_posmax_before_cross, posbin_edges);
    posbin_means = nan(length(posbin_edges)-1, 1); % Store the mean zero-crossing times for each object position bin
    for bin = 1:length(posbin_means)
        bin_members = combined_zero_cross_times(posbin_indices == bin);
        if length(bin_members) >= minBinMembers % Only include bins with at least X members
            posbin_means(bin) = mean(bin_members, 'omitnan');
        end
    end
    % check for gaps
    minSequenceLength = 2;
    posbin_means_clean = removeIsolatedPoints(posbin_means, minSequenceLength);

    posbin_means_per_condition(:,c) = posbin_means_clean;

    % calculate mean angular velocity difference
    meanpangular = mean([winangular_r -winangular_l], 2, 'omitnan');

    % Store the results in the output
    dc_output(:,c) = meanpangular;


    if optPlot
        nexttile; hold on
        if ~isempty(winpanelps_l)
            plot(dc_t,winpanelps_l,'Color',[0.2 0.2 0.2])
        end
        if ~isempty(winpanelps_r)
            plot(dc_t,winpanelps_r,'Color',"#0072BD")
        end
        axis tight; ylim([-180 180]); yline(0); xline(0); ylabel('Obj (deg)'); xlabel('Time (ms)')

        nexttile; hold on
        if ~isempty(winangular_l)
            plot(dc_t,winangular_l,'Color',[0.2 0.2 0.2])
        end
        if ~isempty(winangular_r)
            plot(dc_t,winangular_r,'Color',"#0072BD")
        end
        axis tight; ylim([-400 400]); yline(0); xline(0); ylabel('Ang (deg/s)'); xlabel('Time (ms)')

        nexttile; hold on
        % Column 1: Plot binned means and linear fit on regular axis
        if ~isempty(angmax_before_cross_l)
            scatter(angmax_before_cross_l, angzero_cross_times_l, '.', 'MarkerEdgeColor', [0.2 0.2 0.2])
        end
        if ~isempty(angmax_before_cross_r)
            scatter(angmax_before_cross_r, angzero_cross_times_r, '.', 'MarkerEdgeColor', "#0072BD")
        end
        if ~isempty(angbin_center) && ~isempty(angbin_means)
            plot(angbin_center, angbin_means, 'Color', "r", 'LineWidth', 1.5)
        end

        % Add linear fit for full range
        if ~isempty(combined_angmax_before_cross) && ~isempty(combined_zero_cross_times)
            p_linear = polyfit(combined_angmax_before_cross, combined_zero_cross_times, 1); % Linear fit
            x_full_range = linspace(min(combined_angmax_before_cross), max(combined_angmax_before_cross), 100); % Full range for x-axis
            yfit_linear = polyval(p_linear, x_full_range);
            plot(x_full_range, yfit_linear, 'Color', 'b', 'LineWidth', 1.5) % Plot linear fit over full range
        end

        % Add hyperbolic fit for full range
        if ~isempty(combined_angmax_before_cross) && ~isempty(combined_zero_cross_times)
            % Define hyperbolic model
            hyperbolic_model = @(b, x) (b(1) * x) ./ (b(2) + x);

            % Initial parameter guess: [plateau, half-saturation point]
            initial_guess = [max(combined_zero_cross_times), mean(combined_angmax_before_cross)];

            % Fit the hyperbolic model to your data
            fit_params = nlinfit(combined_angmax_before_cross, combined_zero_cross_times, hyperbolic_model, initial_guess);

            % Generate fit values over full range for plotting
            yfit_hyperbolic = hyperbolic_model(fit_params, x_full_range);

            % Plot the hyperbolic fit in green over full range
            plot(x_full_range, yfit_hyperbolic, 'Color', 'g', 'LineWidth', 1.5) % Plot hyperbolic fit over full range
        end

        ylim([0 250]); xlim([0 max(angbin_edges)]); ylabel('Time 2 Change (ms)'); xlabel('Ang Max (deg/s)')

        nexttile; hold on
        % Column 2: Plot scatter with logx axis and log fit
        if ~isempty(angmax_before_cross_l)
            scatter(angmax_before_cross_l, angzero_cross_times_l, '.', 'MarkerEdgeColor', [0.2 0.2 0.2])
        end
        if ~isempty(angmax_before_cross_r)
            scatter(angmax_before_cross_r, angzero_cross_times_r, '.', 'MarkerEdgeColor', "#0072BD")
        end

        % Add logarithmic fit
        if ~isempty(combined_angmax_before_cross) && all(combined_angmax_before_cross > 0) % Logarithmic fits only work with positive values
            p_log = polyfit(log(combined_angmax_before_cross), combined_zero_cross_times, 1); % Log fit
            yfit_log = polyval(p_log, log(x_full_range));
            plot(x_full_range, yfit_log, 'Color', 'g', 'LineWidth', 1.5) % Plot logarithmic fit over full range
        end

        set(gca, 'XScale', 'log') % Set log scale for x-axis
        ylim([0 250]); xlim([min(angbin_edges(angbin_edges > 0)) max(angbin_edges)]); ylabel('Time 2 Change (ms)'); xlabel('Ang Max (deg/s) [log scale]')

        nexttile; hold on
        % Column 1: Plot binned means and linear fit for object max on regular axis
        if ~isempty(objmax_before_cross_l)
            scatter(objmax_before_cross_l, angzero_cross_times_l, '.', 'MarkerEdgeColor', [0.2 0.2 0.2])
        end
        if ~isempty(objmax_before_cross_r)
            scatter(objmax_before_cross_r, angzero_cross_times_r, '.', 'MarkerEdgeColor', "#0072BD")
        end
        if ~isempty(posbin_center) && ~isempty(posbin_means)
            plot(posbin_center, posbin_means, 'Color', "r", 'LineWidth', 1.5)
        end

        % Add linear fit for full range
        if ~isempty(combined_posmax_before_cross) && ~isempty(combined_zero_cross_times)
            p_linear_pos = polyfit(combined_posmax_before_cross, combined_zero_cross_times, 1); % Linear fit for object position
            x_full_range_pos = linspace(min(combined_posmax_before_cross), max(combined_posmax_before_cross), 100); % Full range for x-axis
            yfit_linear_pos = polyval(p_linear_pos, x_full_range_pos);
            plot(x_full_range_pos, yfit_linear_pos, 'Color', 'b', 'LineWidth', 1.5) % Plot linear fit for object position over full range
        end

        % Add hyperbolic fit for full range
        if ~isempty(combined_posmax_before_cross) && ~isempty(combined_zero_cross_times)
            % Define hyperbolic model for object position
            hyperbolic_model_pos = @(b, x) (b(1) * x) ./ (b(2) + x);

            % Initial parameter guess: [plateau, half-saturation point]
            initial_guess_pos = [max(combined_zero_cross_times), mean(combined_posmax_before_cross)];

            % Fit the hyperbolic model to object position data
            fit_params_pos = nlinfit(combined_posmax_before_cross, combined_zero_cross_times, hyperbolic_model_pos, initial_guess_pos);

            % Generate fit values for plotting over full range
            yfit_hyperbolic_pos = hyperbolic_model_pos(fit_params_pos, x_full_range_pos);

            % Plot the hyperbolic fit in green over full range
            plot(x_full_range_pos, yfit_hyperbolic_pos, 'Color', 'g', 'LineWidth', 1.5) % Plot hyperbolic fit for object position over full range
        end

        ylim([0 250]); xlim([0 max(posbin_edges)]); ylabel('Time 2 Change (ms)'); xlabel('Obj Max (deg)')

        nexttile; hold on
        % Column 2: Plot scatter with logx axis and log fit for object max
        if ~isempty(objmax_before_cross_l)
            scatter(objmax_before_cross_l, angzero_cross_times_l, '.', 'MarkerEdgeColor', [0.2 0.2 0.2])
        end
        if ~isempty(objmax_before_cross_r)
            scatter(objmax_before_cross_r, angzero_cross_times_r, '.', 'MarkerEdgeColor', "#0072BD")
        end

        % Add logarithmic fit for object position
        if ~isempty(combined_posmax_before_cross) && all(combined_posmax_before_cross > 0)
            p_log_pos = polyfit(log(combined_posmax_before_cross), combined_zero_cross_times, 1); % Log fit for object position
            yfit_log_pos = polyval(p_log_pos, log(x_full_range_pos));
            plot(x_full_range_pos, yfit_log_pos, 'Color', 'g', 'LineWidth', 1.5) % Plot logarithmic fit for object position over full range
        end

        set(gca, 'XScale', 'log') % Set log scale for x-axis
        ylim([0 250]); xlim([min(posbin_edges(posbin_edges > 0)) max(posbin_edges)]); ylabel('Time 2 Change (ms)'); xlabel('Obj Max (deg) [log scale]')

    end

end

end
