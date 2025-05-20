% dark_velocityxcorr (Firing Rate Only)
% This script performs a cross-correlation analysis on the firing rate of AOTU019 and 
% AOTU025 neurons under two conditions: with P1 stimulation (+P1) and without P1 (-P1).
% Voltage-related data and plots have been removed.
% 
% CREATED: 11/07/2024 - MC

%% Load Binned Data
clear; close all

% Define folder containing the .mat files
plotPath = 'E:\\Compare Motion Pulse';
dataPath = 'E:\\Compare Motion Pulse\\data';

disp('Analyzing binned data.')
cd(dataPath)

% Load binned data
AOTU019_binned = load(fullfile(dataPath, 'AOTU019_Motion_Pulse_velocity_binned.mat'));
AOTU025_binned = load(fullfile(dataPath, 'AOTU025_Motion_Pulse_velocity_binned.mat'));
fwdBins = AOTU019_binned.thisVelL.fwdBin;

AOTU019_srvfwd = AOTU019_binned.velL_srvfwd;
AOTU025_srvfwd = AOTU025_binned.velL_srvfwd;

% Prepare variables
nBins = size(AOTU019_srvfwd, 1);
nFlies_019 = size(AOTU019_srvfwd, 2);
nFlies_025 = size(AOTU025_srvfwd, 2);

% Total number of observations
totalRows = nBins * (nFlies_019 + nFlies_025);

% Initialize table variables
spikerate = [AOTU019_srvfwd(:); AOTU025_srvfwd(:)];
fwdBin = repmat(fwdBins(:), nFlies_019 + nFlies_025, 1);
cellType = [repmat({'AOTU019'}, nBins * nFlies_019, 1); repmat({'AOTU025'}, nBins * nFlies_025, 1)];
animalID = [repelem((1:nFlies_019)', nBins); repelem((1:nFlies_025)', nBins) + nFlies_019];

% Create table for LME
T = table(spikerate, fwdBin, categorical(cellType), categorical(animalID), ...
    'VariableNames', {'Spikerate', 'FwdBin', 'CellType', 'AnimalID'});

% Fit linear mixed-effects model
lme = fitlme(T, 'Spikerate ~ FwdBin * CellType + (1|AnimalID)');

% Display results
disp(anova(lme));

%% Load angular velocity binned data
cd(dataPath)

% Load binned data
AOTU019_binned = load(fullfile(dataPath, 'AOTU019_Motion_Pulse_velocity_binnedbyfr.mat'));
AOTU025_binned = load(fullfile(dataPath, 'AOTU025_Motion_Pulse_velocity_binnedbyfr.mat'));
angBins = AOTU019_binned.thisSpikeBin;

AOTU019_srvang = AOTU019_binned.sr_velang;
AOTU025_srvang = AOTU025_binned.sr_velang;

% Prepare variables
nBins = size(AOTU019_srvang, 1);
nFlies_019 = size(AOTU019_srvang, 2);
nFlies_025 = size(AOTU025_srvang, 2);

% Total number of observations
totalRows = nBins * (nFlies_019 + nFlies_025);

% Initialize table variables
spikerate = [AOTU019_srvang(:); AOTU025_srvang(:)];
frBin = repmat(angBins(:), nFlies_019 + nFlies_025, 1);
cellType = [repmat({'AOTU019'}, nBins * nFlies_019, 1); repmat({'AOTU025'}, nBins * nFlies_025, 1)];
animalID = [repelem((1:nFlies_019)', nBins); repelem((1:nFlies_025)', nBins) + nFlies_019];

% Create table for LME
T = table(spikerate, frBin, categorical(cellType), categorical(animalID), ...
    'VariableNames', {'Spikerate', 'FRBin', 'CellType', 'AnimalID'});

% Fit linear mixed-effects model
lme = fitlme(T, 'Spikerate ~ FRBin * CellType + (1|AnimalID)');

% Display results
disp(anova(lme));

%% Load cross correlation data
clear; close all

% Define folder containing the .mat files
plotPath = 'E:\\Compare Motion Pulse';
dataPath = 'E:\\Compare Motion Pulse\\data';

%% Load Cross-Correlation Data
disp('Analyzing xcorr in darkness.')
cd(dataPath)

% Load xcorr data
AOTU019_fr = load(fullfile(dataPath, 'AOTU019_Motion_Pulse_xcorr.mat'));
AOTU025_fr = load(fullfile(dataPath, 'AOTU025_Motion_Pulse_xcorr.mat'));
AOTU019_r_pk_fwd_p1 = AOTU019_fr.forward_peaks.r_pk;
AOTU019_lag_pk_fwd_p1 = AOTU019_fr.forward_peaks.lag_pk;
AOTU025_r_pk_fwd_p1 = AOTU025_fr.forward_peaks.r_pk;
AOTU025_lag_pk_fwd_p1 = AOTU025_fr.forward_peaks.lag_pk;

% Load fitted slope data
AOTU019_slope = load(fullfile(dataPath, 'AOTU019_Motion_Pulse_velocity_slopes.mat'));
AOTU025_slope = load(fullfile(dataPath, 'AOTU025_Motion_Pulse_velocity_slopes.mat'));
AOTU019_fwd_fit = AOTU019_slope.fits.fwd;
AOTU025_fwd_fit = AOTU025_slope.fits.fwd;
AOTU019_fwd_r2 = AOTU019_slope.fits.fwdr2;
AOTU025_fwd_r2 = AOTU025_slope.fits.fwdr2;

%% Filter by r_pk > 0.15
min_rpk = 0.15;
keep_019 = AOTU019_r_pk_fwd_p1 > min_rpk;
keep_025 = AOTU025_r_pk_fwd_p1 > min_rpk;

AOTU019_r_pk_fwd_p1 = AOTU019_r_pk_fwd_p1(keep_019);
AOTU019_lag_pk_fwd_p1 = AOTU019_lag_pk_fwd_p1(keep_019);
AOTU025_r_pk_fwd_p1 = AOTU025_r_pk_fwd_p1(keep_025);
AOTU025_lag_pk_fwd_p1 = AOTU025_lag_pk_fwd_p1(keep_025);

fprintf('Total included AOTU019: %d\n', sum(keep_019));
fprintf('Total included AOTU025: %d\n', sum(keep_025));

%% Plot Combined xcorr data
mean_colors = {[0 0 1], [0.5 0 0.5]};
jitter_amount = 0.1;

figure; set(gcf, 'Position', [100 100 600 400])
tiledlayout(1,2, 'TileSpacing', 'compact')

% Tile 1: +P1 Peak Correlation
nexttile;
scatter(1 + jitter_amount * (rand(size(AOTU019_r_pk_fwd_p1)) - 0.5), AOTU019_r_pk_fwd_p1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 + jitter_amount * (rand(size(AOTU025_r_pk_fwd_p1)) - 0.5), AOTU025_r_pk_fwd_p1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
plot(1, median(AOTU019_r_pk_fwd_p1, 'omitnan'), '_', 'Color', mean_colors{1}, 'MarkerSize', 15, 'LineWidth', 2);
plot(2, median(AOTU025_r_pk_fwd_p1, 'omitnan'), '_', 'Color', mean_colors{2}, 'MarkerSize', 15, 'LineWidth', 2);
title('Peak Correlation (+P1)');
xticks([1 2]); xticklabels({'AOTU019','AOTU025'});
ylabel('Correlation'); xlim([0 3]); ylim([0 0.7]); yline(0);
[~, p1] = ttest2(AOTU019_r_pk_fwd_p1, AOTU025_r_pk_fwd_p1);
text(2.8, 0.65, sprintf('p = %.3f', p1), 'HorizontalAlignment', 'right');

% Tile 2: +P1 Lag
nexttile;
scatter(1 + jitter_amount * (rand(size(AOTU019_lag_pk_fwd_p1)) - 0.5), AOTU019_lag_pk_fwd_p1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 + jitter_amount * (rand(size(AOTU025_lag_pk_fwd_p1)) - 0.5), AOTU025_lag_pk_fwd_p1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
plot(1, median(AOTU019_lag_pk_fwd_p1, 'omitnan'), '_', 'Color', mean_colors{1}, 'MarkerSize', 15, 'LineWidth', 2);
plot(2, median(AOTU025_lag_pk_fwd_p1, 'omitnan'), '_', 'Color', mean_colors{2}, 'MarkerSize', 15, 'LineWidth', 2);
title('Lag Estimate (+P1)');
xticks([1 2]); xticklabels({'AOTU019','AOTU025'});
ylabel('Lag (ms)'); xlim([0 3]); ylim([-175 175]); yline(0);
[~, p2] = ttest2(AOTU019_lag_pk_fwd_p1, AOTU025_lag_pk_fwd_p1);
text(2.8, -50, sprintf('p = %.3f', p2), 'HorizontalAlignment', 'right');

sgtitle('FWD XCORR +P1');
cd(plotPath)
saveas(gcf, 'xcorr_mp_fwd.png');
saveas(gcf, 'xcorr_mp_fwd.svg');

%% Plot Combined Slope Data
mean_colors = {[0 0 1], [0.5 0 0.5]};
jitter_amount = 0.05;

% Plot forward slope distribution for AOTU019 and AOTU025
figure; set(gcf, 'Position', [100 100 600 400])
tiledlayout(1,1, 'TileSpacing', 'compact');

nexttile;
hold on;

% Jittered x-values for distribution
x1 = 1 + jitter_amount * (rand(size(AOTU019_fwd_fit)) - 0.5);
x2 = 2 + jitter_amount * (rand(size(AOTU025_fwd_fit)) - 0.5);

% Plot AOTU019 forward slopes
scatter(x1, AOTU019_fwd_fit, '.', 'MarkerEdgeColor', mean_colors{1});
% Plot AOTU025 forward slopes
scatter(x2, AOTU025_fwd_fit, '.', 'MarkerEdgeColor', mean_colors{2});

% Plot medians
plot(1, median(AOTU019_fwd_fit, 'omitnan'), '_', 'Color', mean_colors{1}, 'MarkerSize', 15, 'LineWidth', 2);
plot(2, median(AOTU025_fwd_fit, 'omitnan'), '_', 'Color', mean_colors{2}, 'MarkerSize', 15, 'LineWidth', 2);

% Format axes
xticks([1 2]);
xticklabels({'AOTU019','AOTU025'});
ylabel('Fitted slope');
xlim([0 3]);
ylim([-1 4]);
yline(0, '--', 'Color', [0.5 0.5 0.5]);

% Statistical test (t-test)
[~, p1] = ttest2(AOTU019_fwd_fit, AOTU025_fwd_fit);
text(2.8, 3.75, sprintf('p = %.3f', p1), 'HorizontalAlignment', 'right');

sgtitle('Forward Slope Distribution (+P1)');

% Save the figure
cd(plotPath)
saveas(gcf, 'slopes_mp_fwd_distribution.png');
saveas(gcf, 'slopes_mp_fwd_distribution.svg');

%% Plot R² values for fitted forward slope
figure; set(gcf, 'Position', [100 100 600 400])
tiledlayout(1,1, 'TileSpacing', 'compact');

nexttile;
hold on;

% Jittered x-values
x1 = 1 + jitter_amount * (rand(size(AOTU019_fwd_r2)) - 0.5);
x2 = 2 + jitter_amount * (rand(size(AOTU025_fwd_r2)) - 0.5);

% Plot R² values
scatter(x1, AOTU019_fwd_r2, '.', 'MarkerEdgeColor', mean_colors{1});
scatter(x2, AOTU025_fwd_r2, '.', 'MarkerEdgeColor', mean_colors{2});

% Plot medians
plot(1, median(AOTU019_fwd_r2, 'omitnan'), '_', 'Color', mean_colors{1}, 'MarkerSize', 15, 'LineWidth', 2);
plot(2, median(AOTU025_fwd_r2, 'omitnan'), '_', 'Color', mean_colors{2}, 'MarkerSize', 15, 'LineWidth', 2);

% Format axes
xticks([1 2]);
xticklabels({'AOTU019','AOTU025'});
ylabel('R^2 (fit to velocity)');
xlim([0 3]);
ylim([0 1]);
yline(0.5, '--', 'Color', [0.5 0.5 0.5]);

% Statistical comparison (t-test)
[~, p_r2] = ttest2(AOTU019_fwd_r2, AOTU025_fwd_r2);
text(2.8, 0.95, sprintf('p = %.3f', p_r2), 'HorizontalAlignment', 'right');

sgtitle('Fit Quality (R^2) for Forward Slope (+P1)');

% Save the figure
cd(plotPath)
saveas(gcf, 'r2_forward_slope.png');
saveas(gcf, 'r2_forward_slope.svg');
