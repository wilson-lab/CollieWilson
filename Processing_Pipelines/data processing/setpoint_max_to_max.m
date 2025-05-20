%% setpoint_max_to_max
% This function calculates the time difference between the maximum object position
% (target inflection point) and the subsequent maximum angular velocity for
% leftward and rightward zero-crossings of the object position. The time difference
% is then binned by the angular velocity at the object maximum.
%
% INPUTS:
% - panelps: 3D array representing the object position over time.
%   Dimensions are (time, trials, conditions).
% - angular: 3D array representing the angular velocity over time.
%   Dimensions are (time, trials, conditions).
% - jumptrg: 3D array indicating the time of object jumps.
%   Dimensions are (time, trials, conditions).
% - ttime: 1D array representing the time vector for the data.
% - optPlot: Boolean flag indicating whether to plot the data or not.
%
% OUTPUTS:
% - angbin_means: 2D array containing the average time difference (in ms)
%   between the object maximum and subsequent angular velocity maximum,
%   binned according to the angular velocity at the time of the object maximum.
%   Dimensions are (bins, conditions).
% - bins: A structure containing the bin edges and centers for angular velocity.
%
% CREATED: 11/02/2024 - MC
%
function [angbin_means, bins] = setpoint_max_to_max(panelps, angular, jumptrg, ttime, optPlot)
%% Initialize
nTrial = size(panelps, 2);
nCond = size(panelps, 3);

% Set analysis parameters
cross_window = 0.3; % s
idx_window = fetchTimeIdx(ttime, cross_window);

% Define bin edges for angular velocity at the time of the object max
angbin_edges = [0:40:160, 240];
angbin_centers = angbin_edges(1:end - 1) + diff(angbin_edges) / 2;
bins.ang = angbin_centers;

% Preallocate bins for angular means
angbin_means = nan(length(angbin_edges) - 1, nCond);

min_valid_idx = 10;
min_bin_size = 2;

% Optional pre-process panel data by removing bar jumps
tomit = 2; % s
iomit = fetchTimeIdx(ttime, tomit);

for c = 1:nCond
    for t = 1:nTrial
        jumpidx = find(diff(jumptrg(:, t, c)) > 0);
        for j = 1:length(jumpidx)
            panelps(jumpidx(j):jumpidx(j) + iomit, t, c) = nan;
        end
    end
end

%% Find crossings and fetch angular velocity at object max
for c = 1:nCond
    thispanelps = reshape(panelps(:, :, c), [], 1);
    thisangular = reshape(angular(:, :, c), [], 1);

    % Find sign changes indicating zero-crossings
    sign_changes = diff(sign(thispanelps));
    cross_right = find(sign_changes > 0);
    cross_left = find(sign_changes < 0);

    ang_velocities_at_obj_max = [];
    timediff_pos2ang = [];

    for i = [cross_right; cross_left]'
        % Ensure index is within bounds for window extraction
        if i > idx_window && i + idx_window <= length(thispanelps)
            % Fetch windows around the crossing
            thiswin_panelps = thispanelps(i - idx_window:i + idx_window);
            thiswin_angular = thisangular(i - idx_window:i + idx_window);

            % Find object max index within the window
            [~, obj_max_idx] = max(abs(thiswin_panelps));
            angular_velocity_at_obj_max = thiswin_angular(obj_max_idx);

            % Find subsequent maximum angular velocity after object max
            % Start the search immediately after obj_max_idx
            [~, ang_max_relative_idx] = max(abs(thiswin_angular(obj_max_idx+1:end)));

            % Calculate absolute index of ang_max in the entire array
            ang_max_idx = obj_max_idx + ang_max_relative_idx;

            % Ensure ang_max_idx is within bounds of ttime
            if i + ang_max_idx - 1 <= length(ttime)
                timediff = ttime(i + ang_max_idx - 1) - ttime(i + obj_max_idx - 1);
                timediff_pos2ang = [timediff_pos2ang, timediff * 1000]; % Convert to ms
                ang_velocities_at_obj_max = [ang_velocities_at_obj_max, abs(angular_velocity_at_obj_max)];
            end
        end
    end


    %% Bin data by angular velocity at object max and calculate means
    [angbin_counts, ~, angbin_idx] = histcounts(ang_velocities_at_obj_max, angbin_edges);

    for b = 1:length(angbin_edges) - 1
        if angbin_counts(b) >= min_bin_size && length(timediff_pos2ang) >= min_valid_idx
            angbin_means(b, c) = mean(timediff_pos2ang(angbin_idx == b), 'omitnan');
        else
            angbin_means(b, c) = NaN;
        end
    end

    %% Plotting
    if optPlot
        % Plot Angular Velocity at Object Max vs. Time Difference to Angular Max
        figure;
        scatter(ang_velocities_at_obj_max, timediff_pos2ang, 'k.');
        hold on;
        plot(angbin_centers, angbin_means(:, c), 'o-', 'Color', 'k', 'LineWidth', 1.5);
        xlabel('Angular Velocity at Object Max (deg/s)');
        ylabel('Time Difference (ms)');
        title(['Condition ' num2str(c)]);
        set(gca, 'XScale', 'log');
    end
end
end
