% xcorr_summary_motionpulse
% This script performs a cross-correlation analysis on the firing rate and 
% voltage data of AOTU019 and AOTU025 neurons in response to motion pulse stimuli. 
% Cross-correlation peak values (r_pk) and lag times (lag_pk) for rotational 
% (angular) velocity are loaded from preprocessed .mat files and compared between 
% the neurons. Plots are generated to display peak correlations and lag times for 
% both firing rate and voltage data, with statistical comparisons between AOTU019 
% and AOTU025 for each measure.
%
% CREATED: 11/08/2024 - MC
%
%% Load cross correlation data
clear; close all
% Define folder containing the .mat files
plotPath = 'E:\Compare Motion Pulse';
dataPath = 'E:\Compare Motion Pulse\data';

%% load cross-correlation data
disp('Analyzing xcorr for motion pulse.');
cd(dataPath);

% Load cross-correlation data for firing rate
AOTU019_fr_file = fullfile(dataPath, 'AOTU019_Motion_Pulse_xcorr.mat');
AOTU025_fr_file = fullfile(dataPath, 'AOTU025_Motion_Pulse_xcorr.mat');
AOTU019_fr_xc = load(AOTU019_fr_file);
AOTU025_fr_xc = load(AOTU025_fr_file);

% Load cross-correlation data for voltage
AOTU019_vm_file = fullfile(dataPath, 'AOTU019_Motion_Pulse_xcorr_vm.mat');
AOTU025_vm_file = fullfile(dataPath, 'AOTU025_Motion_Pulse_xcorr_vm.mat');
AOTU019_vm_xc = load(AOTU019_vm_file);
AOTU025_vm_xc = load(AOTU025_vm_file);

% Access specific variables within each loaded structure for firing rate
AOTU019_r_pk_ang = AOTU019_fr_xc.rotational_peaks.r_pk;  % Rotational peak r-value
AOTU019_lag_pk_ang = AOTU019_fr_xc.rotational_peaks.lag_pk;  % Rotational peak lag
AOTU025_r_pk_ang = AOTU025_fr_xc.rotational_peaks.r_pk;
AOTU025_lag_pk_ang = AOTU025_fr_xc.rotational_peaks.lag_pk;

% Access specific variables within each loaded structure for voltage
AOTU019_r_pk_angv = AOTU019_vm_xc.rotational_peaks.r_pk;  % Rotational peak r-value for voltage
AOTU019_lag_pk_angv = AOTU019_vm_xc.rotational_peaks.lag_pk;  % Rotational peak lag for voltage
AOTU025_r_pk_angv = AOTU025_vm_xc.rotational_peaks.r_pk;
AOTU025_lag_pk_angv = AOTU025_vm_xc.rotational_peaks.lag_pk;

% Forward data
AOTU019_r_fwd_pk_ang = AOTU019_fr_xc.forward_peaks.r_pk;  % Forward peak r-value for AOTU019
AOTU019_lag_fwd_pk_ang = AOTU019_fr_xc.forward_peaks.lag_pk;  % Forward peak lag for AOTU019
AOTU025_r_fwd_pk_ang = AOTU025_fr_xc.forward_peaks.r_pk;  % Forward peak r-value for AOTU025
AOTU025_lag_fwd_pk_ang = AOTU025_fr_xc.forward_peaks.lag_pk;  % Forward peak lag for AOTU025

%% Compare Cross-Correlation Data
% Define colors for each neuron
mean_colors = {[0 0 1], [0.5 0 0.5]}; % Blue for AOTU019, Purple for AOTU025

% Create a figure
figure; set(gcf, 'Position', [100 100 600 800])
tiledlayout(2, 2, 'TileSpacing', 'compact')

% Jitter amount for scatter plots
jitter_amount = 0.1;

% First row: Firing Rate Data
% Column 1: Peak Correlation (firing rate)
nexttile;
scatter(ones(size(AOTU019_r_pk_ang)) + jitter_amount * (rand(size(AOTU019_r_pk_ang)) - 0.5), AOTU019_r_pk_ang, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 * ones(size(AOTU025_r_pk_ang)) + jitter_amount * (rand(size(AOTU025_r_pk_ang)) - 0.5), AOTU025_r_pk_ang, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
median_AOTU019 = median(AOTU019_r_pk_ang, 'omitnan');
median_AOTU025 = median(AOTU025_r_pk_ang, 'omitnan');
plot(1, median_AOTU019, '_', 'MarkerSize', 15, 'Color', mean_colors{1}, 'LineWidth', 2);
plot(2, median_AOTU025, '_', 'MarkerSize', 15, 'Color', mean_colors{2}, 'LineWidth', 2);
title('Peak Correlation (Firing Rate)');
xticks([1 2]);
xlim([0 3]);
yline(0);
xticklabels({'AOTU019', 'AOTU025'});
ylabel('Correlation Value');
ylim([-0.1, 0.8]);

% Perform two-sample t-test for peak correlation (firing rate)
[~, p_corr_fr] = ttest2(AOTU019_r_pk_ang, AOTU025_r_pk_ang, 'Vartype', 'unequal');
text(2.8, -0.05, sprintf('p = %.3f', p_corr_fr), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Column 2: Peak Lag (firing rate)
nexttile;
scatter(ones(size(AOTU019_lag_pk_ang)) + jitter_amount * (rand(size(AOTU019_lag_pk_ang)) - 0.5), AOTU019_lag_pk_ang, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 * ones(size(AOTU025_lag_pk_ang)) + jitter_amount * (rand(size(AOTU025_lag_pk_ang)) - 0.5), AOTU025_lag_pk_ang, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
median_AOTU019 = median(AOTU019_lag_pk_ang, 'omitnan');
median_AOTU025 = median(AOTU025_lag_pk_ang, 'omitnan');
plot(1, median_AOTU019, '_', 'MarkerSize', 15, 'Color', mean_colors{1}, 'LineWidth', 2);
plot(2, median_AOTU025, '_', 'MarkerSize', 15, 'Color', mean_colors{2}, 'LineWidth', 2);
title('Lag Estimate (Firing Rate)');
xticks([1 2]);
xlim([0 3]);
yline(0);
xticklabels({'AOTU019', 'AOTU025'});
ylabel('Lag (ms)');
ylim([-400, 400]);

% Perform two-sample t-test for lag estimate (firing rate)
[~, p_lag_fr] = ttest2(AOTU019_lag_pk_ang, AOTU025_lag_pk_ang, 'Vartype', 'unequal');
text(2.8, -350, sprintf('p = %.3f', p_lag_fr), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Second row: Voltage Data
% Column 1: Peak Correlation (voltage)
nexttile;
scatter(ones(size(AOTU019_r_pk_angv)) + jitter_amount * (rand(size(AOTU019_r_pk_angv)) - 0.5), AOTU019_r_pk_angv, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 * ones(size(AOTU025_r_pk_angv)) + jitter_amount * (rand(size(AOTU025_r_pk_angv)) - 0.5), AOTU025_r_pk_angv, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
median_AOTU019 = median(AOTU019_r_pk_angv, 'omitnan');
median_AOTU025 = median(AOTU025_r_pk_angv, 'omitnan');
plot(1, median_AOTU019, '_', 'MarkerSize', 15, 'Color', mean_colors{1}, 'LineWidth', 2);
plot(2, median_AOTU025, '_', 'MarkerSize', 15, 'Color', mean_colors{2}, 'LineWidth', 2);
title('Peak Correlation (Voltage)');
xticks([1 2]);
xlim([0 3]);
yline(0);
xticklabels({'AOTU019', 'AOTU025'});
ylabel('Correlation Value');
ylim([-0.1, 0.8]);

% Perform two-sample t-test for peak correlation (voltage)
[~, p_corr_vm] = ttest2(AOTU019_r_pk_angv, AOTU025_r_pk_angv, 'Vartype', 'unequal');
text(2.8, -0.05, sprintf('p = %.3f', p_corr_vm), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Column 2: Peak Lag (voltage)
nexttile;
scatter(ones(size(AOTU019_lag_pk_angv)) + jitter_amount * (rand(size(AOTU019_lag_pk_angv)) - 0.5), AOTU019_lag_pk_angv, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 * ones(size(AOTU025_lag_pk_angv)) + jitter_amount * (rand(size(AOTU025_lag_pk_angv)) - 0.5), AOTU025_lag_pk_angv, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
median_AOTU019 = median(AOTU019_lag_pk_angv, 'omitnan');
median_AOTU025 = median(AOTU025_lag_pk_angv, 'omitnan');
plot(1, median_AOTU019, '_', 'MarkerSize', 15, 'Color', mean_colors{1}, 'LineWidth', 2);
plot(2, median_AOTU025, '_', 'MarkerSize', 15, 'Color', mean_colors{2}, 'LineWidth', 2);
title('Lag Estimate (Voltage)');
xticks([1 2]);
xlim([0 3]);
yline(0);
xticklabels({'AOTU019', 'AOTU025'});
ylabel('Lag (ms)');
ylim([-400, 400]);

% Perform two-sample t-test for lag estimate (voltage)
[~, p_lag_vm] = ttest2(AOTU019_lag_pk_angv, AOTU025_lag_pk_angv, 'Vartype', 'unequal');
text(2.8, -350, sprintf('p = %.3f', p_lag_vm), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Set global title and save plots
cd(plotPath);
sgtitle('Rotational XCORR: Firing Rate and Voltage');
saveas(gcf, 'xcorr_p1_motionpulse.png');
set(gcf, 'renderer', 'Painters'); % Save vectorized version
saveas(gcf, 'xcorr_p1_motionpulse.svg');

%% Compare Cross-Correlation Data for Forward Peaks
% Define colors for each neuron
mean_colors = {[0 0 1], [0.5 0 0.5]}; % Blue for AOTU019, Purple for AOTU025

% Create a figure
figure; set(gcf, 'Position', [100 100 600 800])
tiledlayout(1, 2, 'TileSpacing', 'compact')

% Jitter amount for scatter plots
jitter_amount = 0.1;

% Column 1: Peak Correlation (firing rate)
nexttile;
scatter(ones(size(AOTU019_r_fwd_pk_ang)) + jitter_amount * (rand(size(AOTU019_r_fwd_pk_ang)) - 0.5), AOTU019_r_fwd_pk_ang, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 * ones(size(AOTU025_r_fwd_pk_ang)) + jitter_amount * (rand(size(AOTU025_r_fwd_pk_ang)) - 0.5), AOTU025_r_fwd_pk_ang, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
median_AOTU019 = median(AOTU019_r_fwd_pk_ang, 'omitnan');
median_AOTU025 = median(AOTU025_r_fwd_pk_ang, 'omitnan');
plot(1, median_AOTU019, '_', 'MarkerSize', 15, 'Color', mean_colors{1}, 'LineWidth', 2);
plot(2, median_AOTU025, '_', 'MarkerSize', 15, 'Color', mean_colors{2}, 'LineWidth', 2);
title('Peak Correlation (Firing Rate)');
xticks([1 2]);
xlim([0 3]);
yline(0);
xticklabels({'AOTU019', 'AOTU025'});
ylabel('Correlation Value');
ylim([-0.1, 0.8]);

% Perform two-sample t-test for peak correlation (firing rate)
[~, p_corr_fr] = ttest2(AOTU019_r_fwd_pk_ang, AOTU025_r_fwd_pk_ang, 'Vartype', 'unequal');
text(2.8, -0.05, sprintf('p = %.3f', p_corr_fr), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Column 2: Peak Lag (firing rate)
nexttile;
scatter(ones(size(AOTU019_lag_fwd_pk_ang)) + jitter_amount * (rand(size(AOTU019_lag_fwd_pk_ang)) - 0.5), AOTU019_lag_fwd_pk_ang, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 * ones(size(AOTU025_lag_fwd_pk_ang)) + jitter_amount * (rand(size(AOTU025_lag_fwd_pk_ang)) - 0.5), AOTU025_lag_fwd_pk_ang, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
median_AOTU019 = median(AOTU019_lag_fwd_pk_ang, 'omitnan');
median_AOTU025 = median(AOTU025_lag_fwd_pk_ang, 'omitnan');
plot(1, median_AOTU019, '_', 'MarkerSize', 15, 'Color', mean_colors{1}, 'LineWidth', 2);
plot(2, median_AOTU025, '_', 'MarkerSize', 15, 'Color', mean_colors{2}, 'LineWidth', 2);
title('Lag Estimate (Firing Rate)');
xticks([1 2]);
xlim([0 3]);
yline(0);
xticklabels({'AOTU019', 'AOTU025'});
ylabel('Lag (ms)');
ylim([-150, 150]);

% Perform two-sample t-test for lag estimate (firing rate)
[~, p_lag_fr] = ttest2(AOTU019_lag_fwd_pk_ang, AOTU025_lag_fwd_pk_ang, 'Vartype', 'unequal');
text(2.8, -350, sprintf('p = %.3f', p_lag_fr), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Set global title and save plots
cd(plotPath);
sgtitle('Forward XCORR: Firing Rate');
saveas(gcf, 'xcorr_p1_motionpulse_forward.png');
set(gcf, 'renderer', 'Painters'); % Save vectorized version
saveas(gcf, 'xcorr_p1_motionpulse_forward.svg');
