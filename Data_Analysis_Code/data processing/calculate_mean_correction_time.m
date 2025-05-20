% calculate_mean_correction_time
% This function analyzes the time it takes for the fly to change direction 
% in response to motion pulses by calculating the average correction time 
% for small and large angular turns. It categorizes turns based on thresholds 
% and provides optional plotting of the results for each condition.
%
% INPUTS:
%   panelps      - Panel positions (degrees) representing the target position
%   angular      - Angular velocity (degrees/second)
%   jumptrg      - Trigger for jumps or transitions in behavior
%   ttime        - Time vector (in seconds)
%   optPlot      - Flag to determine whether to plot the results (1 for plotting, 0 to skip)
%
% OUTPUTS:
%   dc_output_small - Average time from angular velocity peaks to direction change for small turns
%   dc_output_large  - Average time from angular velocity peaks to direction change for large turns
%
% CREATED: [Date] MC
%
function [dc_output_small, dc_output_large] = calculate_mean_correction_time(panelps, angular, jumptrg, ttime, optPlot)
%% Initialize
% fetch experiment conditions
nTrial = size(panelps, 2);
nCond = size(panelps, 3);

% set analysis parameters
cross_window = 0.3; % s
idx_window = fetchTimeIdx(ttime, cross_window);

% Preallocate mean correction time matrices
dc_output_small = nan(nCond, 1);
dc_output_large = nan(nCond, 1);

turn_threshold = 100; % deg/s
turn_max = 300; %ddeg/s

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
window_size = length(1 - idx_window:1 + idx_window);
dc_t = (ttime(1:window_size) - cross_window) * 1000;  % convert to ms

if optPlot
    % initialize plot
    figure; set(gcf, 'Position', [100 100 1500 900])
    tiledlayout(nCond, 2, 'TileSpacing', 'compact')
end

% for each condition
for c = 1:nCond
    % initialize
    winangular_r = [];
    winangular_l = [];
    
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

        % Find zero-crossing time after the panel crossing
        sign_change_window = diff(sign(thiswin_angular(idx_window + 1:end))); % search after the panel crossing
        zero_cross_idx = find(sign_change_window > 0, 1, 'first'); % first zero-crossing after the panel crossing
        if ~isempty(zero_cross_idx)
            % Store window
            winangular_r(:, ri) = thiswin_angular;
            % Time of zero crossing
            angzero_cross_times_r(ri) = dc_t(idx_window + 1 + zero_cross_idx);
            % Find max angular value before the angular zero-crossing
            angmax_before_cross_r(ri) = -min(thiswin_angular(1:idx_window + zero_cross_idx), [], 'omitnan');
            ri = ri + 1;
        end
    end

    % Left crossings
    angzero_cross_times_l = [];
    angmax_before_cross_l = [];
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

        % Find zero-crossing time after the panel crossing
        sign_change_window = diff(sign(thiswin_angular(idx_window + 1:end))); % search after the panel crossing
        zero_cross_idx = find(sign_change_window < 0, 1, 'first'); % first zero-crossing after the panel crossing
        if ~isempty(zero_cross_idx)
            % Store window
            winangular_l(:, li) = thiswin_angular;
            % Time of zero crossing
            angzero_cross_times_l(li) = dc_t(idx_window + 1 + zero_cross_idx);
            % Find max angular value before the angular zero-crossing
            angmax_before_cross_l(li) = max(thiswin_angular(1:idx_window + zero_cross_idx), [], 'omitnan');
            li = li + 1;
        end
    end

    % Combine right and left crossings
    combined_zero_cross_times = [angzero_cross_times_r, angzero_cross_times_l]; % combined times
    combined_angmax_before_cross = [angmax_before_cross_r, angmax_before_cross_l]; % combined max angular velocities

    %% Calculate mean correction time for small and large turns
    small_turns_idx = combined_angmax_before_cross < turn_threshold;
    large_turns_idx = (combined_angmax_before_cross >= turn_threshold) & (combined_angmax_before_cross <= turn_max);

    if sum(small_turns_idx) > 0
        dc_output_small(c) = mean(combined_zero_cross_times(small_turns_idx), 'omitnan'); % mean correction time for small turns
    else
        dc_output_small(c) = NaN;
    end

    if sum(large_turns_idx) > 0
        dc_output_large(c) = mean(combined_zero_cross_times(large_turns_idx), 'omitnan'); % mean correction time for large turns
    else
        dc_output_large(c) = NaN;
    end

    %% Plotting
    if optPlot
        % Plot the angular velocity data
        nexttile; hold on
        if ~isempty(winangular_l)
            plot(dc_t, winangular_l, 'Color', [0.2 0.2 0.2])
        end
        if ~isempty(winangular_r)
            plot(dc_t, winangular_r, 'Color', "#0072BD")
        end
        axis tight; ylim([-400 400]); yline(0); xline(0); ylabel('Ang (deg/s)'); xlabel('Time (ms)')
        title(['Condition ' num2str(c)])
    end
end

end
