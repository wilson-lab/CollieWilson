% Comparison of Cross-Correlation Peak Values and Lag Times for AOTU019 and AOTU025
% 
% This script processes cross-correlation data for two Drosophila neurons,
% AOTU019 and AOTU025, under two experimental conditions: motion pulse and
% background pulse. The goal is to compare the peak correlation values and lag times
% of angular velocity cross-correlations for each neuron during dark and pursuit conditions.
%
% Inputs:
% - Cross-correlation data files in each specified folder. Each file should contain:
%   - r_val: Cross-correlation values for angular velocity.
%   - lag_t: Lag time values corresponding to r_val.
%
% Created: 11/08/2024 MC
%
%% Initialize parameters
% Load settings
close all; clear
settings = processSettings();

% Specify folders for AOTU019 and AOTU025 data
folder_motionpulse_019 = 'E:\AOTU019 Motion Pulse\interpolated\xcorr';
folder_backgroundpulse_019 = 'E:\AOTU019 Background P1\interpolated\xcorr';
folder_motionpulse_025 = 'E:\AOTU025 Motion Pulse\interpolated\xcorr';
folder_backgroundpulse_025 = 'E:\AOTU025 Background P1\interpolated\xcorr';

% Get list of files and find shared names for AOTU019
files_motionpulse_019 = dir(fullfile(folder_motionpulse_019, '*_xc.mat'));
files_backgroundpulse_019 = dir(fullfile(folder_backgroundpulse_019, '*_1_xc.mat'));

motionpulse_names_019 = cellfun(@(x) regexp(x, '^(.*)_xc', 'tokens', 'once'), {files_motionpulse_019.name}, 'UniformOutput', false);
backgroundpulse_names_019 = cellfun(@(x) regexp(x, '(?<=2023_)(.*)_1_xc', 'tokens', 'once'), {files_backgroundpulse_019.name}, 'UniformOutput', false);

motionpulse_names_019 = cellfun(@(x) x{1}, motionpulse_names_019, 'UniformOutput', false);
backgroundpulse_names_019 = cellfun(@(x) x{1}, backgroundpulse_names_019, 'UniformOutput', false);
shared_names_019 = intersect(motionpulse_names_019, backgroundpulse_names_019);

% Get list of files and find shared names for AOTU025
files_motionpulse_025 = dir(fullfile(folder_motionpulse_025, '*_xc.mat'));
files_backgroundpulse_025 = dir(fullfile(folder_backgroundpulse_025, '*_1_xc.mat'));

motionpulse_names_025 = cellfun(@(x) regexp(x, '^(.*)_xc', 'tokens', 'once'), {files_motionpulse_025.name}, 'UniformOutput', false);
backgroundpulse_names_025 = cellfun(@(x) regexp(x, '(?<=2024_)(.*)_1_xc', 'tokens', 'once'), {files_backgroundpulse_025.name}, 'UniformOutput', false);

motionpulse_names_025 = cellfun(@(x) x{1}, motionpulse_names_025, 'UniformOutput', false);
backgroundpulse_names_025 = cellfun(@(x) x{1}, backgroundpulse_names_025, 'UniformOutput', false);
shared_names_025 = intersect(motionpulse_names_025, backgroundpulse_names_025);

% Initialize storage arrays for AOTU019
r_val_ang_motion_019 = [];
pk_lag_ang_motion_019 = [];
pk_rval_ang_motion_019 = [];
r_val_ang_background_019 = [];
pk_lag_ang_background_019 = [];
pk_rval_ang_background_019 = [];

% Initialize storage arrays for AOTU025
r_val_ang_motion_025 = [];
pk_lag_ang_motion_025 = [];
pk_rval_ang_motion_025 = [];
r_val_ang_background_025 = [];
pk_lag_ang_background_025 = [];
pk_rval_ang_background_025 = [];

% Load data for AOTU019
for n = 1:length(shared_names_019)
    % Load motion pulse file
    motionFile = fullfile(folder_motionpulse_019, [shared_names_019{n} '_xc_v.mat']);
    load(motionFile, 'r_val', 'lag_t');
    
    % Find and store peak lag and peak r-values for angular velocity in motion pulse
    [peak_lag_motion, peak_rval_motion, r_val] = find_peak_lag_rval(r_val, lag_t, settings.minXCorrPromVm);
    r_val_ang_motion_019(:, n) = r_val.ang;
    pk_lag_ang_motion_019(n) = peak_lag_motion.ang;
    pk_rval_ang_motion_019(n) = peak_rval_motion.ang;

    % Load background pulse file
    backgroundFile = fullfile(folder_backgroundpulse_019, ['2023_' shared_names_019{n} '_1_xc.mat']);
    load(backgroundFile, 'r_val', 'lag_t');

    % Find and store peak lag and peak r-values for angular velocity in background pulse
    [peak_lag_background, peak_rval_background, r_val] = find_peak_lag_rval(r_val, lag_t, settings.minXCorrPromVm);
    r_val_ang_background_019(:, n) = r_val.ang;
    pk_lag_ang_background_019(n) = peak_lag_background.ang;
    pk_rval_ang_background_019(n) = peak_rval_background.ang;
end

% Load data for AOTU025
for n = 1:length(shared_names_025)
    % Load motion pulse file
    motionFile = fullfile(folder_motionpulse_025, [shared_names_025{n} '_xc_v.mat']);
    load(motionFile, 'r_val', 'lag_t');
    
    % Find and store peak lag and peak r-values for angular velocity in motion pulse
    [peak_lag_motion, peak_rval_motion, r_val] = find_peak_lag_rval(r_val, lag_t, settings.minXCorrPromVm);
    r_val_ang_motion_025(:, n) = r_val.ang;
    pk_lag_ang_motion_025(n) = peak_lag_motion.ang;
    pk_rval_ang_motion_025(n) = peak_rval_motion.ang;

    % Load background pulse file
    backgroundFile = fullfile(folder_backgroundpulse_025, ['2024_' shared_names_025{n} '_1_xc.mat']);
    load(backgroundFile, 'r_val', 'lag_t');

    % Find and store peak lag and peak r-values for angular velocity in background pulse
    [peak_lag_background, peak_rval_background, r_val] = find_peak_lag_rval(r_val, lag_t, settings.minXCorrPromVm);
    r_val_ang_background_025(:, n) = r_val.ang;
    pk_lag_ang_background_025(n) = peak_lag_background.ang;
    pk_rval_ang_background_025(n) = peak_rval_background.ang;
end

% Prepare combined data for plotting and analysis for AOTU019
combined_rval_ang_019 = [pk_rval_ang_background_019; pk_rval_ang_motion_019];
combined_lag_ang_019 = [pk_lag_ang_background_019; pk_lag_ang_motion_019];
valid_rval_idx_019 = ~any(isnan(combined_rval_ang_019), 1);
valid_lag_idx_019 = ~any(isnan(combined_lag_ang_019), 1);
filtered_rval_ang_019 = combined_rval_ang_019(:, valid_rval_idx_019);
filtered_lag_ang_019 = combined_lag_ang_019(:, valid_lag_idx_019);

% Prepare combined data for plotting and analysis for AOTU025
combined_rval_ang_025 = [pk_rval_ang_background_025; pk_rval_ang_motion_025];
combined_lag_ang_025 = [pk_lag_ang_background_025; pk_lag_ang_motion_025];
valid_rval_idx_025 = ~any(isnan(combined_rval_ang_025), 1);
valid_lag_idx_025 = ~any(isnan(combined_lag_ang_025), 1);
filtered_rval_ang_025 = combined_rval_ang_025(:, valid_rval_idx_025);
filtered_lag_ang_025 = combined_lag_ang_025(:, valid_lag_idx_025);

% Paired t-tests for AOTU019
[~, p_rval_019] = ttest(filtered_rval_ang_019(1, :), filtered_rval_ang_019(2, :));
[~, p_lag_019] = ttest(filtered_lag_ang_019(1, :), filtered_lag_ang_019(2, :));

% Paired t-tests for AOTU025
[~, p_rval_025] = ttest(filtered_rval_ang_025(1, :), filtered_rval_ang_025(2, :));
[~, p_lag_025] = ttest(filtered_lag_ang_025(1, :), filtered_lag_ang_025(2, :));

% Create a tiled layout with 1 row and 4 columns for both neurons
figure; set(gcf, 'Position', [100 100 400 800]);
tiledlayout(2,2, 'TileSpacing', 'compact');
marker_size = 10;

% Plot for AOTU019 - Peak Correlation
nexttile;
hold on;
for i = 1:size(filtered_rval_ang_019, 2)
    plot([1, 2], filtered_rval_ang_019(:, i), '-k');
end
scatter(ones(1, size(filtered_rval_ang_019, 2)), filtered_rval_ang_019(1, :), marker_size, 'k', 'filled');
scatter(2 * ones(1, size(filtered_rval_ang_019, 2)), filtered_rval_ang_019(2, :), marker_size, 'k', 'filled');
title('AOTU019 Peak Correlation');
xticks([1 2]); xticklabels({'Dark', 'Pursuit'});
ylabel('Correlation Value');
yline(0); ylim([-0.1 0.6]); xlim([0, 3]);
text(2.5, 0.55, sprintf('p(Cond) = %.3f', p_rval_019), 'FontSize', 8, 'HorizontalAlignment', 'right');

% Plot for AOTU019 - Lag Values
nexttile;
hold on;
for i = 1:size(filtered_lag_ang_019, 2)
    plot([1, 2], filtered_lag_ang_019(:, i), '-k');
end
scatter(ones(1, size(filtered_lag_ang_019, 2)), filtered_lag_ang_019(1, :), marker_size, 'k', 'filled');
scatter(2 * ones(1, size(filtered_lag_ang_019, 2)), filtered_lag_ang_019(2, :), marker_size, 'k', 'filled');
title('AOTU019 Lag Values');
xticks([1 2]); xticklabels({'Dark', 'Pursuit'});
ylabel('Lag (ms)');
xlim([0, 3]); ylim([-500 0]);
text(2.5, -30, sprintf('p(Cond) = %.3f', p_lag_019), 'FontSize', 8, 'HorizontalAlignment', 'right');

% Plot for AOTU025 - Peak Correlation
nexttile;
hold on;
for i = 1:size(filtered_rval_ang_025, 2)
    plot([1, 2], filtered_rval_ang_025(:, i), '-k');
end
scatter(ones(1, size(filtered_rval_ang_025, 2)), filtered_rval_ang_025(1, :), marker_size, 'k', 'filled');
scatter(2 * ones(1, size(filtered_rval_ang_025, 2)), filtered_rval_ang_025(2, :), marker_size, 'k', 'filled');
title('AOTU025 Peak Correlation');
xticks([1 2]); xticklabels({'Dark', 'Pursuit'});
ylabel('Correlation Value');
yline(0); ylim([-0.1 0.6]); xlim([0, 3]);
text(2.5, 0.55, sprintf('p(Cond) = %.3f', p_rval_025), 'FontSize', 8, 'HorizontalAlignment', 'right');

% Plot for AOTU025 - Lag Values
nexttile;
hold on;
for i = 1:size(filtered_lag_ang_025, 2)
    plot([1, 2], filtered_lag_ang_025(:, i), '-k');
end
scatter(ones(1, size(filtered_lag_ang_025, 2)), filtered_lag_ang_025(1, :), marker_size, 'k', 'filled');
scatter(2 * ones(1, size(filtered_lag_ang_025, 2)), filtered_lag_ang_025(2, :), marker_size, 'k', 'filled');
title('AOTU025 Lag Values');
xticks([1 2]); xticklabels({'Dark', 'Pursuit'});
ylabel('Lag (ms)');
xlim([0, 3]); ylim([-500 0]);
text(2.5, -30, sprintf('p(Cond) = %.3f', p_lag_025), 'FontSize', 8, 'HorizontalAlignment', 'right');

% Save the plots
saveas(gcf, 'AOTU_XCORRacross.png');
set(gcf, 'renderer', 'Painters');
saveas(gcf, 'AOTU_XCORRacross.svg');

disp('Comparison plot for AOTU019 and AOTU025 complete.');
