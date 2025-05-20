% SETPOINT_DIRCHANGE - This function analyzes direction changes by evaluating object position and 
% angular velocity data during setpoint fixation. It calculates the time between panel position 
% crossings and corresponding angular velocity zero-crossings and bins this data by both angular 
% velocity magnitude prior to and at the time of crossing.

% INPUTS:
%   panelps   - 3D array of panel position data, where each slice represents a condition, 
%               and each column represents a trial.
%   angular   - 3D array of angular velocity data corresponding to panel position data.
%   jumptrg   - 3D array indicating bar jump triggers.
%   ttime     - Time vector corresponding to the temporal resolution of the data.
%   optPlot   - Optional flag (1/0) for plotting results (1 for yes, 0 for no).

% OUTPUTS:
%   dc_t                - Time vector centered around direction changes (in milliseconds).
%   angbin_means        - Mean time to zero-crossing for each bin based on angular velocity prior to crossing.
%   cross_angvel_means  - Mean time to zero-crossing for each bin based on angular velocity at crossing.
%   bins                - Structure containing angular velocity bin edges and centers.

% Created: 10/??/2024 - MC
% Updated: 10/30/2024 - MC (Added binning by angular velocity at crossing)
%
function [dc_t, angbin_means, cross_angvel_means, bins] = setpoint_dirchange(panelps, angular, jumptrg, ttime, optPlot)
%% Initialize
% fetch experiment conditions
nTrial = size(panelps, 2);
nCond = size(panelps, 3);

% set analysis parameters
cross_window = 0.3; % s
idx_window = fetchTimeIdx(ttime, cross_window);

% Define fit ranges for angular velocity
ang_range = 10:400; % Avoid starting at 0 for log scale
bins.ang = ang_range;

% Define bin edges for angular velocity
angbin_edges = [0:40:160, 240];  % Bins for angular velocity
% define centers
angbin_centers = angbin_edges(1:end - 1) + diff(angbin_edges) / 2;
bins.ang2 = angbin_centers;

% Preallocate bins for angular means
angbin_means = nan(length(angbin_edges) - 1, nCond);
cross_angvel_means =  nan(length(angbin_edges) - 1, nCond);

min_valid_idx = 10;
min_bin_size = 2; % Minimum number of elements in bin for mean calculation

%% optional pre-process panel data by removing bar jumps
% omit timepoints right after a bar jump
tomit = 2; % s
iomit = fetchTimeIdx(ttime, tomit); % idx

for c = 1:nCond
    for t = 1:nTrial
        jumpidx = find(diff(jumptrg(:, t, c)) > 0);
        for j = 1:length(jumpidx)
            panelps(jumpidx(j):jumpidx(j) + iomit, t, c) = nan;
        end
    end
end

%% find crossings
% initialize
dc_output = [];
window_size = length(1 - idx_window:1 + idx_window);
dc_t = (ttime(1:window_size) - cross_window) * 1000;  % convert to ms
if optPlot
    % initialize plot
    figure; set(gcf, 'Position', [100 100 1500 900])
    tiledlayout(nCond, 4, 'TileSpacing', 'compact')
end

% for each condition
for c = 1:nCond
    % initialize window storage
    winpanelps_r = [];
    winpanelps_l = [];
    winangular_r = [];
    winangular_l = [];
    nang = 0; % Reset nang for each condition
    % estimate HD bias
    biasHD = mean(panelps(:, :, c), 'all', 'omitnan');

    % fetch data
    thispanelps = reshape(panelps(:, :, c) - biasHD, [], 1);
    thisangular = reshape(angular(:, :, c), [], 1);

    % calculate the difference in sign of the position values
    sign_changes = diff(sign(thispanelps));
    % find indices where sign changes occur (this indicates a zero-crossing)
    cross_right = find(sign_changes > 0);
    cross_left = find(sign_changes < 0);

    % Right crossings
    angzero_cross_times_r = [];
    angmax_before_cross_r = [];
    abs_angvel_at_cross_r = [];
    ri = 1;
    for i = 1:length(cross_right)
        idx = cross_right(i);
        % Fetch angular data before and after the crossing
        thiswin_angular = thisangular(idx - idx_window:idx + idx_window);
        thiswin_panelps = thispanelps(idx - idx_window:idx + idx_window);

        % Check if the difference between any adjacent panel position values exceeds 10
        if any(abs(diff(thiswin_panelps)) > 10)
            continue; % Skip this window if any adjacent differences exceed 10
        end

        % Find angular velocity zero-crossing time after the panel crossing
        sign_change_window = diff(sign(thiswin_angular(idx_window + 1:end))); % search after the panel crossing
        zero_cross_idx = find(sign_change_window > 0, 1, 'first'); % first angular velocity zero-crossing after the panel crossing
        if ~isempty(zero_cross_idx)
            % Store window data
            winpanelps_r(:, ri) = thiswin_panelps;
            winangular_r(:, ri) = thiswin_angular;
            % Store the time difference between panel crossing and angular velocity zero-crossing
            angzero_cross_times_r(ri) = dc_t(idx_window + 1 + zero_cross_idx) - dc_t(idx_window + 1);
            % Find max angular value before the angular zero-crossing
            angmax_before_cross_r(ri) = -min(thiswin_angular(1:idx_window + zero_cross_idx), [], 'omitnan');
            % Fetch absolute angular velocity at object crossing
            abs_angvel_at_cross_r(ri) = abs(thisangular(idx));
            ri = ri + 1;
            nang = nang + 1; % Increment nang for each valid crossing
        end
    end

    % Left crossings
    angzero_cross_times_l = [];
    angmax_before_cross_l = [];
    abs_angvel_at_cross_l = [];
    li = 1;
    for i = 1:length(cross_left)
        idx = cross_left(i);
        % Fetch angular data before and after the crossing
        thiswin_angular = thisangular(idx - idx_window:idx + idx_window);
        thiswin_panelps = thispanelps(idx - idx_window:idx + idx_window);

        % Check if the difference between any adjacent panel position values exceeds 10
        if any(abs(diff(thiswin_panelps)) > 10)
            continue; % Skip this window if any adjacent differences exceed 10
        end

        % Find angular velocity zero-crossing time after the panel crossing
        sign_change_window = diff(sign(thiswin_angular(idx_window + 1:end))); % search after the panel crossing
        zero_cross_idx = find(sign_change_window < 0, 1, 'first'); % first angular velocity zero-crossing after the panel crossing
        if ~isempty(zero_cross_idx)
            % Store window data
            winpanelps_l(:, li) = thiswin_panelps;
            winangular_l(:, li) = thiswin_angular;
            % Store the time difference between panel crossing and angular velocity zero-crossing
            angzero_cross_times_l(li) = dc_t(idx_window + 1 + zero_cross_idx) - dc_t(idx_window + 1);
            % Find max angular value before the angular zero-crossing
            angmax_before_cross_l(li) = max(thiswin_angular(1:idx_window + zero_cross_idx), [], 'omitnan');
            % Fetch absolute angular velocity at object crossing
            abs_angvel_at_cross_l(li) = abs(thisangular(idx));
            li = li + 1;
            nang = nang + 1; % Increment nang for each valid crossing
        end
    end

    % Combine right and left crossings
    combined_zero_cross_times = [angzero_cross_times_r, angzero_cross_times_l]; % combined times
    combined_angmax_before_cross = [angmax_before_cross_r, angmax_before_cross_l]; % combined max values before crossing
    combined_abs_angvel_at_cross = [abs_angvel_at_cross_r, abs_angvel_at_cross_l]; % combined absolute angular velocities at crossings

    %% Bin data
    % Bin data and calculate means for angular velocity prior to crossing
    [angbin_counts_prior, ~, angbin_idx_prior] = histcounts(combined_angmax_before_cross, angbin_edges);

    for b = 1:length(angbin_edges) - 1
        if angbin_counts_prior(b) >= min_bin_size && nang >= min_valid_idx
            % Calculate the mean time to zero-crossing for each angular velocity bin prior to crossing
            angbin_means(b, c) = mean(combined_zero_cross_times(angbin_idx_prior == b), 'omitnan');
        else
            angbin_means(b, c) = NaN;
        end
    end

    % Bin data and calculate means for angular velocity at crossing
    [angbin_counts_cross, ~, angbin_idx_cross] = histcounts(combined_abs_angvel_at_cross, angbin_edges);

    for b = 1:length(angbin_edges) - 1
        if angbin_counts_cross(b) >= min_bin_size && nang >= min_valid_idx
            % Calculate the mean time to zero-crossing for each angular velocity bin at crossing
            cross_angvel_means(b, c) = mean(combined_zero_cross_times(angbin_idx_cross == b), 'omitnan');
        else
            cross_angvel_means(b, c) = NaN;
        end
    end

    %% Plotting
    if optPlot

        % Plot the panel position data
        nexttile; hold on
        if ~isempty(winpanelps_l)
            plot(dc_t, winpanelps_l, 'Color', [0.2 0.2 0.2])
        end
        if ~isempty(winpanelps_r)
            plot(dc_t, winpanelps_r, 'Color', "#0072BD")
        end
        axis tight; ylim([-180 180]); yline(0); xline(0); ylabel('Obj (deg)'); xlabel('Time (ms)')

        % Plot the angular velocity data
        nexttile; hold on
        if ~isempty(winangular_l)
            plot(dc_t, winangular_l, 'Color', [0.2 0.2 0.2])
        end
        if ~isempty(winangular_r)
            plot(dc_t, winangular_r, 'Color', "#0072BD")
        end
        axis tight; ylim([-400 400]); yline(0); xline(0); ylabel('Ang (deg/s)'); xlabel('Time (ms)')

        % Plot angular max before crossing vs. time to zero-crossing
        nexttile; hold on
        if ~isempty(angmax_before_cross_l)
            scatter(angmax_before_cross_l, angzero_cross_times_l, '.', 'MarkerEdgeColor', [0.2 0.2 0.2])
        end
        if ~isempty(angmax_before_cross_r)
            scatter(angmax_before_cross_r, angzero_cross_times_r, '.', 'MarkerEdgeColor', "#0072BD")
        end
        plot(angbin_centers, angbin_means(:, c), 'o-', 'Color', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); % Plot binned means
        legend(num2str(nang))
        set(gca, 'XTick', [20, 50, 100, 200, 400]) % More ticks for angular max x-axis
        ylim([0 350]); xlim([0 400]); ylabel('Time 2 Change (ms)'); xlabel('Ang Max (deg/s)')

        % Plot angular velocity at crossing vs. time to zero-crossing
        nexttile; hold on
        if ~isempty(abs_angvel_at_cross_l)
            scatter(abs_angvel_at_cross_l, angzero_cross_times_l, '.', 'MarkerEdgeColor', [0.2 0.2 0.2])
        end
        if ~isempty(abs_angvel_at_cross_r)
            scatter(abs_angvel_at_cross_r, angzero_cross_times_r, '.', 'MarkerEdgeColor', "#0072BD")
        end
        plot(angbin_centers, cross_angvel_means(:, c), 'o-', 'Color', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); % Plot binned means for angular velocity at crossing
        ylim([0 350]); xlim([0 400]); ylabel('Time 2 Change (ms)'); xlabel('Ang at Crossing (deg/s)')
        set(gca, 'XTick', [20, 50, 100, 200, 400]) % More ticks for x-axis
    end

end
end
