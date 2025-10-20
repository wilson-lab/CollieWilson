% dark_velocityxcorr
% This script performs a cross-correlation analysis on the firing rate and 
% voltage data of AOTU019 and AOTU025 neurons under two conditions: with P1 
% stimulation (+P1) and without P1 stimulation (No P1). Cross-correlation 
% peak values and lag times are extracted from preprocessed .mat files and 
% plotted for comparison.
% 
% CREATED: 11/07/2024 - MC
%
%% Load cross correlation data
clear; close all

% Define folder containing the .mat files
plotPath = 'E:\Compare Motion Pulse';
dataPath = 'E:\Compare Motion Pulse\data';

%% load cross-correlation data for P1
disp('Analyzing analyzing xcorr in darkness +P1.')
% Load cross-correlation data for firing rate
cd(dataPath)
AOTU019_fr_file = fullfile(dataPath, 'AOTU019_Background_P1_xcorr.mat');
AOTU025_fr_file = fullfile(dataPath, 'AOTU025_Background_P1_xcorr.mat');
AOTU019_fr_xc = load(AOTU019_fr_file);
AOTU025_fr_xc = load(AOTU025_fr_file);

% Load cross-correlation data for voltage
AOTU019_vm_file = fullfile(dataPath, 'AOTU019_Background_P1_xcorrvm.mat');
AOTU025_vm_file = fullfile(dataPath, 'AOTU025_Background_P1_xcorrvm.mat');
AOTU019_vm_xc = load(AOTU019_vm_file);
AOTU025_vm_xc = load(AOTU025_vm_file);

% Access specific variables within each loaded structure for firing rate
AOTU019_r_pk_fwd = AOTU019_fr_xc.combined_data.r_pk_fwd;
AOTU019_r_pk_ang = AOTU019_fr_xc.combined_data.r_pk_ang;
AOTU019_lag_pk_fwd = AOTU019_fr_xc.combined_data.lag_pk_fwd;
AOTU019_lag_pk_ang = AOTU019_fr_xc.combined_data.lag_pk_ang;

AOTU025_r_pk_fwd = AOTU025_fr_xc.combined_data.r_pk_fwd;
AOTU025_r_pk_ang = AOTU025_fr_xc.combined_data.r_pk_ang;
AOTU025_lag_pk_fwd = AOTU025_fr_xc.combined_data.lag_pk_fwd;
AOTU025_lag_pk_ang = AOTU025_fr_xc.combined_data.lag_pk_ang;

% Access specific variables within each loaded structure for voltage
AOTU019_r_pk_fwdv = AOTU019_vm_xc.combined_data.r_pk_fwdv;
AOTU019_r_pk_angv = AOTU019_vm_xc.combined_data.r_pk_angv;
AOTU019_lag_pk_fwdv = AOTU019_vm_xc.combined_data.lag_pk_fwdv;
AOTU019_lag_pk_angv = AOTU019_vm_xc.combined_data.lag_pk_angv;

AOTU025_r_pk_fwdv = AOTU025_vm_xc.combined_data.r_pk_fwdv;
AOTU025_r_pk_angv = AOTU025_vm_xc.combined_data.r_pk_angv;
AOTU025_lag_pk_fwdv = AOTU025_vm_xc.combined_data.lag_pk_fwdv;
AOTU025_lag_pk_angv = AOTU025_vm_xc.combined_data.lag_pk_angv;

%% Compare cross correlation data
% Define colors for each neuron
mean_colors = {[0 0 1], [0.5 0 0.5]}; % Blue for AOTU019, Purple for AOTU025

% Create a figure
figure; set(gcf, 'Position', [100 100 600 800])
tiledlayout(2, 2, 'TileSpacing', 'compact')

% Jitter amount
jitter_amount = 0.1;

% First row: Angular Peak Correlation - Firing Rate
% First tile: Angular peak correlation (firing rate)
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
ylim([-0.0, 0.5]);

% Perform two-sample t-test for peak correlation (firing rate)
[~, p_corr_fr] = ttest2(AOTU019_r_pk_ang, AOTU025_r_pk_ang, 'Vartype', 'unequal');
text(2.8, 0.45, sprintf('p = %.3f', p_corr_fr), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Second tile: Angular lag estimate (firing rate)
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
ylim([-450, 000]);

% Perform two-sample t-test for lag estimate (firing rate)
[~, p_lag_fr] = ttest2(AOTU019_lag_pk_ang, AOTU025_lag_pk_ang, 'Vartype', 'unequal');
text(2.8, -50, sprintf('p = %.3f', p_lag_fr), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Second row: Angular Peak Correlation - Voltage
% Third tile: Angular peak correlation (voltage)
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
ylim([0, 0.7]);

% Perform two-sample t-test for peak correlation (voltage)
[~, p_corr_vm] = ttest2(AOTU019_r_pk_angv, AOTU025_r_pk_angv, 'Vartype', 'unequal');
text(2.8, -0.05, sprintf('p = %.3f', p_corr_vm), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Fourth tile: Angular lag estimate (voltage)
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
ylim([-300, 000]);

% Perform two-sample t-test for lag estimate (voltage)
[~, p_lag_vm] = ttest2(AOTU019_lag_pk_angv, AOTU025_lag_pk_angv, 'Vartype', 'unequal');
text(2.8, -350, sprintf('p = %.3f', p_lag_vm), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Adjust layout
sgtitle('Rotational XCORR (+P1)');

% Save plots as images
cd(plotPath)
saveas(gcf, 'xcorr_p1_darkness.png');
set(gcf, 'renderer', 'Painters'); % Save vectorized version
saveas(gcf, 'xcorr_p1_darkness.svg');

%% Load Cross-Correlation Data for no P1
disp('Analyzing analyzing xcorr in darkness -P1.')

% Access specific variables within each loaded structure for firing rate
AOTU019_r_pk_fwd = AOTU019_fr_xc.combined_data.r_pk_fwd_nop1;
AOTU019_r_pk_ang = AOTU019_fr_xc.combined_data.r_pk_ang_nop1;
AOTU019_lag_pk_fwd = AOTU019_fr_xc.combined_data.lag_pk_fwd_nop1;
AOTU019_lag_pk_ang = AOTU019_fr_xc.combined_data.lag_pk_ang_nop1;

AOTU025_r_pk_fwd = AOTU025_fr_xc.combined_data.r_pk_fwd_nop1;
AOTU025_r_pk_ang = AOTU025_fr_xc.combined_data.r_pk_ang_nop1;
AOTU025_lag_pk_fwd = AOTU025_fr_xc.combined_data.lag_pk_fwd_nop1;
AOTU025_lag_pk_ang = AOTU025_fr_xc.combined_data.lag_pk_ang_nop1;

% Access specific variables within each loaded structure for voltage
AOTU019_r_pk_fwdv = AOTU019_vm_xc.combined_data.r_pk_fwdv_nop1;
AOTU019_r_pk_angv = AOTU019_vm_xc.combined_data.r_pk_angv_nop1;
AOTU019_lag_pk_fwdv = AOTU019_vm_xc.combined_data.lag_pk_fwdv_nop1;
AOTU019_lag_pk_angv = AOTU019_vm_xc.combined_data.lag_pk_angv_nop1;

AOTU025_r_pk_fwdv = AOTU025_vm_xc.combined_data.r_pk_fwdv_nop1;
AOTU025_r_pk_angv = AOTU025_vm_xc.combined_data.r_pk_angv_nop1;
AOTU025_lag_pk_fwdv = AOTU025_vm_xc.combined_data.lag_pk_fwdv_nop1;
AOTU025_lag_pk_angv = AOTU025_vm_xc.combined_data.lag_pk_angv_nop1;

%% Compare Cross-Correlation Data
% Define colors for each neuron
mean_colors = {[0 0 1], [0.5 0 0.5]}; % Blue for AOTU019, Purple for AOTU025

% Create a 2x2 tiled layout
figure; set(gcf, 'Position', [100 100 600 800])
tiledlayout(2, 2, 'TileSpacing', 'compact')

% Jitter amount
jitter_amount = 0.1;

% First row: Angular Peak Correlation - Firing Rate
% First tile: Angular peak correlation (firing rate)
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
ylim([0, 0.7]);

% Perform two-sample t-test for peak correlation (firing rate)
[~, p_corr_fr] = ttest2(AOTU019_r_pk_ang, AOTU025_r_pk_ang, 'Vartype', 'unequal');
text(2.8, -0.05, sprintf('p = %.3f', p_corr_fr), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Second tile: Angular lag estimate (firing rate)
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
ylim([-300, 0]);

% Perform two-sample t-test for lag estimate (firing rate)
[~, p_lag_fr] = ttest2(AOTU019_lag_pk_ang, AOTU025_lag_pk_ang, 'Vartype', 'unequal');
text(2.8, -350, sprintf('p = %.3f', p_lag_fr), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Second row: Angular Peak Correlation - Voltage
% Third tile: Angular peak correlation (voltage)
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
ylim([0, 0.7]);

% Perform two-sample t-test for peak correlation (voltage)
[~, p_corr_vm] = ttest2(AOTU019_r_pk_angv, AOTU025_r_pk_angv, 'Vartype', 'unequal');
text(2.8, -0.05, sprintf('p = %.3f', p_corr_vm), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Fourth tile: Angular lag estimate (voltage)
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
ylim([-300, 0]);

% Perform two-sample t-test for lag estimate (voltage)
[~, p_lag_vm] = ttest2(AOTU019_lag_pk_angv, AOTU025_lag_pk_angv, 'Vartype', 'unequal');
text(2.8, -350, sprintf('p = %.3f', p_lag_vm), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Adjust layout
sgtitle('Rotational XCORR (-P1)');

% Save plots as images
cd(plotPath)
saveas(gcf, 'xcorr_nop1_darkness.png');
set(gcf, 'renderer', 'Painters'); % Save vectorized version
saveas(gcf, 'xcorr_nop1_darkness.svg');

