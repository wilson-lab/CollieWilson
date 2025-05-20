% mp_p1difference.m
% CREATED: 11/18/2024 - MC

% This script compares the firing rate responses of AOTU019 and AOTU025 
% neurons during P1 activation versus no P1 activation. The analysis involves:

%% Initialize
clear
close all

% Define paths
plotPath = 'E:\Compare Motion Pulse';
dataPath = 'E:\Compare Motion Pulse\data';

%% Fetch data
% Change directory to data path
cd(dataPath)

% Define file names
AOTU019_P1_file = fullfile(dataPath, 'AOTU019_Motion_Pulse_sweeppeaks_25_dps_Quiescent.mat');
AOTU019_noP1_file = fullfile(dataPath, 'AOTU019_Motion_Pulse_No_P1_sweeppeaks_25_dps_Quiescent.mat');
AOTU025_P1_file = fullfile(dataPath, 'AOTU025_Motion_Pulse_sweeppeaks_25_dps_Quiescent.mat');
AOTU025_noP1_file = fullfile(dataPath, 'AOTU025_Motion_Pulse_No_P1_sweeppeaks_25_dps_Quiescent.mat');

% Load each file
AOTU019_P1_data = load(AOTU019_P1_file);
AOTU019_noP1_data = load(AOTU019_noP1_file);
AOTU025_P1_data = load(AOTU025_P1_file);
AOTU025_noP1_data = load(AOTU025_noP1_file);

% Extract fly names and firing rate data
AOTU019_P1_fly_names = AOTU019_P1_data.flyShortNames;
AOTU019_noP1_fly_names = AOTU019_noP1_data.flyShortNames;
AOTU025_P1_fly_names = AOTU025_P1_data.flyShortNames;
AOTU025_noP1_fly_names = AOTU025_noP1_data.flyShortNames;

% Find common animals between P1 and noP1 conditions for each AOTU
common_AOTU019_animals = intersect(AOTU019_P1_fly_names, AOTU019_noP1_fly_names);
common_AOTU025_animals = intersect(AOTU025_P1_fly_names, AOTU025_noP1_fly_names);

% Filter firing rate data (pulse_srR) for the common animals
AOTU019_P1_filtered = [];
AOTU019_noP1_filtered = [];
for i = 1:length(common_AOTU019_animals)
    idx_P1 = find(strcmp(AOTU019_P1_fly_names, common_AOTU019_animals{i}));
    idx_noP1 = find(strcmp(AOTU019_noP1_fly_names, common_AOTU019_animals{i}));
    AOTU019_P1_filtered = cat(2, AOTU019_P1_filtered, AOTU019_P1_data.pulse_srR(:, idx_P1, :));
    AOTU019_noP1_filtered = cat(2, AOTU019_noP1_filtered, AOTU019_noP1_data.pulse_srR(:, idx_noP1, :));
end

AOTU025_P1_filtered = [];
AOTU025_noP1_filtered = [];
for i = 1:length(common_AOTU025_animals)
    idx_P1 = find(strcmp(AOTU025_P1_fly_names, common_AOTU025_animals{i}));
    idx_noP1 = find(strcmp(AOTU025_noP1_fly_names, common_AOTU025_animals{i}));
    AOTU025_P1_filtered = cat(2, AOTU025_P1_filtered, AOTU025_P1_data.pulse_srR(:, idx_P1, :));
    AOTU025_noP1_filtered = cat(2, AOTU025_noP1_filtered, AOTU025_noP1_data.pulse_srR(:, idx_noP1, :));
end

%% Find peak responses for AOTU019 and AOTU025

% Find indices for AOTU019 where sweepPosR equals 11 and 34
AOTU019_sweepPos = AOTU019_P1_data.sweepPosR;
idx_AOTU019_peaks = find(AOTU019_sweepPos == 11 | AOTU019_sweepPos == 34);
%idx_AOTU019_peaks = find(AOTU019_sweepPos == 34);

% Fetch only the two z arrays corresponding to these indices across time and animals for AOTU019
AOTU019_P1_peaks = AOTU019_P1_filtered(:, :, idx_AOTU019_peaks);
AOTU019_noP1_peaks = AOTU019_noP1_filtered(:, :, idx_AOTU019_peaks);

% Find indices for AOTU025 where sweepPosR equals 56 and 79
AOTU025_sweepPos = AOTU025_P1_data.sweepPosR;
idx_AOTU025_peaks = find(AOTU025_sweepPos == 56 | AOTU025_sweepPos == 79);
%idx_AOTU025_peaks = find(AOTU025_sweepPos == 56);

% Fetch only the two z arrays corresponding to these indices across time and animals for AOTU025
AOTU025_P1_peaks = AOTU025_P1_filtered(:, :, idx_AOTU025_peaks);
AOTU025_noP1_peaks = AOTU025_noP1_filtered(:, :, idx_AOTU025_peaks);

% Display a message confirming the data extraction
disp('Peak sweep data extracted for AOTU019 and AOTU025.');

%% Calculate Difference in Activity (P1 - no P1) for Peak Sweeps

% Initialize containers for average differences
AOTU019_diff_avg = [];
AOTU025_diff_avg = [];

% Define middle half indices for rows (time dimension)
total_time_points = size(AOTU019_P1_peaks, 1); % Assuming rows correspond to time
start_idx = floor(total_time_points / 4) + 1; % Start from 1/4 in
end_idx = floor(3 * total_time_points / 4);   % End at 3/4 in

% Process AOTU019 data
for i = 1:size(AOTU019_P1_peaks, 2) % Loop through each animal (columns)
    % Extract only the middle half of the time data
    P1_peaks_mid = AOTU019_P1_peaks(start_idx:end_idx, i, :);
    noP1_peaks_mid = AOTU019_noP1_peaks(start_idx:end_idx, i, :);

    % Calculate the difference (P1 - no P1) and average across time (rows)
    diff_019 = median(mean(P1_peaks_mid - noP1_peaks_mid, 1, 'omitnan'),'omitnan');
    AOTU019_diff_avg = [AOTU019_diff_avg; diff_019(:)];
end

% Process AOTU025 data
for i = 1:size(AOTU025_P1_peaks, 2) % Loop through each animal (columns)
    % Extract only the middle half of the time data
    P1_peaks_mid = AOTU025_P1_peaks(start_idx:end_idx, i, :);
    noP1_peaks_mid = AOTU025_noP1_peaks(start_idx:end_idx, i, :);
    % Calculate the difference (P1 - no P1) and average across time (rows)
    diff_025 = median(mean(P1_peaks_mid - noP1_peaks_mid, 1, 'omitnan'),'omitnan');
    AOTU025_diff_avg = [AOTU025_diff_avg; diff_025(:)];
end

% Calculate the median differences for plotting
AOTU019_median_diff = median(AOTU019_diff_avg, 'omitnan');
AOTU025_median_diff = median(AOTU025_diff_avg, 'omitnan');

% Perform Independent t-test
[h, p] = ttest2(AOTU019_diff_avg, AOTU025_diff_avg);

%% Plot Scatter of Average Differences

figure;
hold on;

% Set jitter amount
jitterAmount = 0.005;

% Scatter plot of average differences with jitter
scatter(1 + randn(size(AOTU019_diff_avg)) * jitterAmount, AOTU019_diff_avg, ...
    '.', 'MarkerEdgeColor', [0.5, 0.5, 0.5]); % Grey points for AOTU019
scatter(2 + randn(size(AOTU025_diff_avg)) * jitterAmount, AOTU025_diff_avg, ...
    '.', 'MarkerEdgeColor', [0.5, 0.5, 0.5]); % Grey points for AOTU025

% Plot the median as a dash (no jitter)
plot(1, AOTU019_median_diff, '_', 'Color', 'b', 'MarkerSize', 12); % Blue dash for AOTU019
plot(2, AOTU025_median_diff, '_', 'Color', [0.5, 0, 0.5], 'MarkerSize', 12); % Purple dash for AOTU025

% Display the p-value on the plot
text(2.5, max([AOTU019_diff_avg; AOTU025_diff_avg]), sprintf('p = %.3f', p), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);

% Formatting the plot
xlim([0.5, 2.5]);
xticks([1, 2]);
xticklabels({'AOTU019', 'AOTU025'});
ylabel('Average Difference in Firing Rate (P1 - no P1)');
title('Average Difference in Firing Rate for AOTU019 and AOTU025');
hold off;

% Save plots as images
cd(plotPath)
saveas(gcf, 'p1_compare_visual.png');
set(gcf, 'renderer', 'Painters'); % Save vectorized version
saveas(gcf, 'p1_compare_visual.svg');


%% Calculate Means for P1 and noP1 Peak Responses

% Initialize containers for means
AOTU019_P1_means = [];
AOTU019_noP1_means = [];
AOTU025_P1_means = [];
AOTU025_noP1_means = [];

% Define middle half indices for rows (time dimension)
total_time_points = size(AOTU019_P1_peaks, 1); % Assuming rows correspond to time
start_idx = floor(total_time_points / 4) + 1; % Start from 1/4 in
end_idx = floor(3 * total_time_points / 4);   % End at 3/4 in

% Process AOTU019 data
for i = 1:size(AOTU019_P1_peaks, 2) % Loop through each animal (columns)
    % Extract only the middle half of the time data
    P1_peaks_mid = AOTU019_P1_peaks(start_idx:end_idx, i, :);
    noP1_peaks_mid = AOTU019_noP1_peaks(start_idx:end_idx, i, :);

    % Check for NaNs in the middle half of the data
    if all(~isnan(P1_peaks_mid), 'all') && all(~isnan(noP1_peaks_mid), 'all')
        % Calculate the mean across time (rows)
        AOTU019_P1_means = [AOTU019_P1_means; mean(P1_peaks_mid(:), 'omitnan')];
        AOTU019_noP1_means = [AOTU019_noP1_means; mean(noP1_peaks_mid(:), 'omitnan')];
    end
end

% Process AOTU025 data
for i = 1:size(AOTU025_P1_peaks, 2) % Loop through each animal (columns)
    % Extract only the middle half of the time data
    P1_peaks_mid = AOTU025_P1_peaks(start_idx:end_idx, i, :);
    noP1_peaks_mid = AOTU025_noP1_peaks(start_idx:end_idx, i, :);

    % Check for NaNs in the middle half of the data
    if all(~isnan(P1_peaks_mid), 'all') && all(~isnan(noP1_peaks_mid), 'all')
        % Calculate the mean across time (rows)
        AOTU025_P1_means = [AOTU025_P1_means; mean(P1_peaks_mid(:), 'omitnan')];
        AOTU025_noP1_means = [AOTU025_noP1_means; mean(noP1_peaks_mid(:), 'omitnan')];
    end
end

% Display a message confirming the mean calculation
disp('Means for P1 and noP1 peak responses calculated.');

%% Combine Data for ANOVA
% Create response vector
anova_data = [AOTU019_P1_means; AOTU019_noP1_means; AOTU025_P1_means; AOTU025_noP1_means];

% Create group labels
anova_group_celltype = [repmat({'AOTU019'}, length(AOTU019_P1_means), 1); ...
                        repmat({'AOTU019'}, length(AOTU019_noP1_means), 1); ...
                        repmat({'AOTU025'}, length(AOTU025_P1_means), 1); ...
                        repmat({'AOTU025'}, length(AOTU025_noP1_means), 1)];
anova_group_condition = [repmat({'P1'}, length(AOTU019_P1_means), 1); ...
                         repmat({'noP1'}, length(AOTU019_noP1_means), 1); ...
                         repmat({'P1'}, length(AOTU025_P1_means), 1); ...
                         repmat({'noP1'}, length(AOTU025_noP1_means), 1)];

% Create table for ANOVA
anova_tbl = table(anova_data, anova_group_celltype, anova_group_condition, ...
    'VariableNames', {'Response', 'CellType', 'Condition'});

%% Perform Two-Way ANOVA on Means
[p, tbl, stats] = anovan(anova_tbl.Response, ...
    {anova_tbl.CellType, anova_tbl.Condition}, ...
    'model', 'interaction', ...
    'varnames', {'Cell Type', 'Condition'});

% Display ANOVA table
disp('ANOVA Results:');
disp(tbl);

% Perform post-hoc tests if interaction is significant
if p(3) < 0.05 % Interaction term is the 3rd p-value
    disp('Significant interaction detected. Performing post-hoc tests...');
    [c, m] = multcompare(stats, 'Dimension', [1, 2]);
    disp('Post-hoc comparison results:');
    disp(c);
else
    disp('No significant interaction detected.');
end

