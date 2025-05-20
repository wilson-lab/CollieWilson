% dark_velocityxcorr (Firing Rate Only)
% This script performs a cross-correlation analysis on the firing rate of AOTU019 and 
% AOTU025 neurons under two conditions: with P1 stimulation (+P1) and without P1 (-P1).
% Voltage-related data and plots have been removed.
% 
% CREATED: 11/07/2024 - MC

%% Load cross correlation data
clear; close all

% Define folder containing the .mat files
plotPath = 'E:\\Compare Motion Pulse';
dataPath = 'E:\\Compare Motion Pulse\\data';

%% Load Cross-Correlation Data
disp('Analyzing xcorr in darkness.')
cd(dataPath)

AOTU019_fr = load(fullfile(dataPath, 'AOTU019_Background_P1_xcorr.mat'));
AOTU025_fr = load(fullfile(dataPath, 'AOTU025_Background_P1_xcorr.mat'));

% +P1 data
AOTU019_r_pk_fwd_p1 = AOTU019_fr.combined_data.r_pk_fwd;
AOTU019_lag_pk_fwd_p1 = AOTU019_fr.combined_data.lag_pk_fwd;
AOTU025_r_pk_fwd_p1 = AOTU025_fr.combined_data.r_pk_fwd;
AOTU025_lag_pk_fwd_p1 = AOTU025_fr.combined_data.lag_pk_fwd;

% -P1 data
AOTU019_r_pk_fwd_nop1 = AOTU019_fr.combined_data.r_pk_fwd_nop1;
AOTU019_lag_pk_fwd_nop1 = AOTU019_fr.combined_data.lag_pk_fwd_nop1;
AOTU025_r_pk_fwd_nop1 = AOTU025_fr.combined_data.r_pk_fwd_nop1;
AOTU025_lag_pk_fwd_nop1 = AOTU025_fr.combined_data.lag_pk_fwd_nop1;

%% Clean lag values and remove rows with NaNs

% Set lag values beyond threshold to NaN
AOTU019_lag_pk_fwd_p1(abs(AOTU019_lag_pk_fwd_p1) > 400) = NaN;
AOTU025_lag_pk_fwd_p1(abs(AOTU025_lag_pk_fwd_p1) > 400) = NaN;
AOTU019_lag_pk_fwd_nop1(abs(AOTU019_lag_pk_fwd_nop1) > 400) = NaN;
AOTU025_lag_pk_fwd_nop1(abs(AOTU025_lag_pk_fwd_nop1) > 400) = NaN;

% Remove flies with NaNs in either correlation or lag
valid_019_p1 = ~isnan(AOTU019_r_pk_fwd_p1) & ~isnan(AOTU019_lag_pk_fwd_p1);
valid_025_p1 = ~isnan(AOTU025_r_pk_fwd_p1) & ~isnan(AOTU025_lag_pk_fwd_p1);
valid_019_nop1 = ~isnan(AOTU019_r_pk_fwd_nop1) & ~isnan(AOTU019_lag_pk_fwd_nop1);
valid_025_nop1 = ~isnan(AOTU025_r_pk_fwd_nop1) & ~isnan(AOTU025_lag_pk_fwd_nop1);

% Filter the data
AOTU019_r_pk_fwd_p1 = AOTU019_r_pk_fwd_p1(valid_019_p1);
AOTU019_lag_pk_fwd_p1 = AOTU019_lag_pk_fwd_p1(valid_019_p1);
AOTU025_r_pk_fwd_p1 = AOTU025_r_pk_fwd_p1(valid_025_p1);
AOTU025_lag_pk_fwd_p1 = AOTU025_lag_pk_fwd_p1(valid_025_p1);

AOTU019_r_pk_fwd_nop1 = AOTU019_r_pk_fwd_nop1(valid_019_nop1);
AOTU019_lag_pk_fwd_nop1 = AOTU019_lag_pk_fwd_nop1(valid_019_nop1);
AOTU025_r_pk_fwd_nop1 = AOTU025_r_pk_fwd_nop1(valid_025_nop1);
AOTU025_lag_pk_fwd_nop1 = AOTU025_lag_pk_fwd_nop1(valid_025_nop1);

% Print retained n
fprintf('n AOTU019 -P1: %d\n', numel(AOTU019_r_pk_fwd_nop1));
fprintf('n AOTU019 +P1: %d\n', numel(AOTU019_r_pk_fwd_p1));
fprintf('n AOTU025 -P1: %d\n', numel(AOTU025_r_pk_fwd_nop1));
fprintf('n AOTU025 +P1: %d\n', numel(AOTU025_r_pk_fwd_p1));

%% Plot Combined
mean_colors = {[0 0 1], [0.5 0 0.5]};
jitter_amount = 0.1;

figure; set(gcf, 'Position', [100 100 600 800])
tiledlayout(2, 2, 'TileSpacing', 'compact')

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
ylabel('Lag (ms)'); xlim([0 3]); ylim([-400 400]); yline(0);
[~, p2] = ttest2(AOTU019_lag_pk_fwd_p1, AOTU025_lag_pk_fwd_p1);
text(2.8, -50, sprintf('p = %.3f', p2), 'HorizontalAlignment', 'right');

% Tile 3: -P1 Peak Correlation
nexttile;
scatter(1 + jitter_amount * (rand(size(AOTU019_r_pk_fwd_nop1)) - 0.5), AOTU019_r_pk_fwd_nop1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 + jitter_amount * (rand(size(AOTU025_r_pk_fwd_nop1)) - 0.5), AOTU025_r_pk_fwd_nop1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
plot(1, median(AOTU019_r_pk_fwd_nop1, 'omitnan'), '_', 'Color', mean_colors{1}, 'MarkerSize', 15, 'LineWidth', 2);
plot(2, median(AOTU025_r_pk_fwd_nop1, 'omitnan'), '_', 'Color', mean_colors{2}, 'MarkerSize', 15, 'LineWidth', 2);
title('Peak Correlation (-P1)');
xticks([1 2]); xticklabels({'AOTU019','AOTU025'});
ylabel('Correlation'); xlim([0 3]); ylim([0 0.7]); yline(0);
[~, p3] = ttest2(AOTU019_r_pk_fwd_nop1, AOTU025_r_pk_fwd_nop1);
text(2.8, 0.65, sprintf('p = %.3f', p3), 'HorizontalAlignment', 'right');

% Tile 4: -P1 Lag
nexttile;
scatter(1 + jitter_amount * (rand(size(AOTU019_lag_pk_fwd_nop1)) - 0.5), AOTU019_lag_pk_fwd_nop1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 + jitter_amount * (rand(size(AOTU025_lag_pk_fwd_nop1)) - 0.5), AOTU025_lag_pk_fwd_nop1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
plot(1, median(AOTU019_lag_pk_fwd_nop1, 'omitnan'), '_', 'Color', mean_colors{1}, 'MarkerSize', 15, 'LineWidth', 2);
plot(2, median(AOTU025_lag_pk_fwd_nop1, 'omitnan'), '_', 'Color', mean_colors{2}, 'MarkerSize', 15, 'LineWidth', 2);
title('Lag Estimate (-P1)');
xticks([1 2]); xticklabels({'AOTU019','AOTU025'});
ylabel('Lag (ms)'); xlim([0 3]); ylim([-400 400]); yline(0);
[~, p4] = ttest2(AOTU019_lag_pk_fwd_nop1, AOTU025_lag_pk_fwd_nop1);
text(2.8, -50, sprintf('p = %.3f', p4), 'HorizontalAlignment', 'right');

sgtitle('FWD XCORR: +P1 (top) vs -P1 (bottom)');
cd(plotPath)
saveas(gcf, 'xcorr_darkness_fwd.png');
saveas(gcf, 'xcorr_darkness_fwd.svg');

%% Full Comparison Plot and ANOVA

% Prepare data vectors and metadata for ANOVA
peak_corr_values = [AOTU019_r_pk_fwd_nop1(:); AOTU019_r_pk_fwd_p1(:); AOTU025_r_pk_fwd_nop1(:); AOTU025_r_pk_fwd_p1(:)];
lag_values = [AOTU019_lag_pk_fwd_nop1(:); AOTU019_lag_pk_fwd_p1(:); AOTU025_lag_pk_fwd_nop1(:); AOTU025_lag_pk_fwd_p1(:)];

n_019 = size(AOTU019_r_pk_fwd_nop1, 1);
n_025 = size(AOTU025_r_pk_fwd_nop1, 1);

cell_type = [repmat("AOTU019", 2*n_019, 1); repmat("AOTU025", 2*n_025, 1)];
arousal = [repmat("-P1", n_019, 1); repmat("+P1", n_019, 1); repmat("-P1", n_025, 1); repmat("+P1", n_025, 1)];
fly_id = [ (1:n_019)'; (1:n_019)'; (1:n_025)'+n_019; (1:n_025)'+n_019 ];

% Remove NaNs
valid_idx_peak = ~isnan(peak_corr_values);
valid_idx_lag = ~isnan(lag_values);

% Subset for peak
peak_corr_values = peak_corr_values(valid_idx_peak);
cell_type_peak = cell_type(valid_idx_peak);
arousal_peak = arousal(valid_idx_peak);
fly_id_peak = fly_id(valid_idx_peak);

% Subset for lag
lag_values = lag_values(valid_idx_lag);
cell_type_lag = cell_type(valid_idx_lag);
arousal_lag = arousal(valid_idx_lag);
fly_id_lag = fly_id(valid_idx_lag);

% Group labels for plotting
group_peak = strcat(cell_type_peak, "_", arousal_peak);
group_lag = strcat(cell_type_lag, "_", arousal_lag);

% Run ANOVA using fitlme (Linear Mixed Effects Model)
tbl_peak = table(peak_corr_values, cell_type_peak, arousal_peak, fly_id_peak, ...
    'VariableNames', {'Peak', 'Cell', 'Arousal', 'Fly'});
tbl_lag = table(lag_values, cell_type_lag, arousal_lag, fly_id_lag, ...
    'VariableNames', {'Lag', 'Cell', 'Arousal', 'Fly'});

lme_peak = fitlme(tbl_peak, 'Peak ~ Cell*Arousal + (1|Fly)');
peak_stats = anova(lme_peak);

lme_lag = fitlme(tbl_lag, 'Lag ~ Cell*Arousal + (1|Fly)');
lag_stats = anova(lme_lag);

% Plot with scatter + jitter + median lines
figure;
tiledlayout(1, 2);
jitter = 0.1;

% PEAK
nexttile;
groups = unique(group_peak, 'stable');
hold on;
for i = 1:numel(groups)
    idx = strcmp(group_peak, groups{i});
    x = i + jitter*(rand(sum(idx),1)-0.5);
    scatter(x, peak_corr_values(idx), '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
    plot(i, median(peak_corr_values(idx), 'omitnan'), '_', 'MarkerSize', 15, 'Color', 'k', 'LineWidth', 2);
end
xlim([0 numel(groups)+1]);
ylim([0 0.7]);
title('Peak Correlation');
ylabel('Correlation');
xticks(1:numel(groups));
xticklabels(groups);
ptext = sprintf('p(Cell) = %.3f\np(Arousal) = %.3f\np(Interact) = %.3f', ...
    peak_stats.pValue(2), peak_stats.pValue(3), peak_stats.pValue(4));
text(numel(groups)+0.9, 0.65, ptext, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

% LAG
nexttile;
groups = unique(group_lag, 'stable');
hold on;
for i = 1:numel(groups)
    idx = strcmp(group_lag, groups{i});
    x = i + jitter*(rand(sum(idx),1)-0.5);
    scatter(x, lag_values(idx), '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
    plot(i, median(lag_values(idx), 'omitnan'), '_', 'MarkerSize', 15, 'Color', 'k', 'LineWidth', 2);
end
xlim([0 numel(groups)+1]);
ylim([-450 450]);
title('Lag Estimate');
ylabel('Lag (ms)');
xticks(1:numel(groups));
xticklabels(groups);
ptext = sprintf('p(Cell) = %.3f\np(Arousal) = %.3f\np(Interact) = %.3f', ...
    lag_stats.pValue(2), lag_stats.pValue(3), lag_stats.pValue(4));
text(numel(groups)+0.9, 400, ptext, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

cd(plotPath)
saveas(gcf, 'xcorr_darkness_fwdcombined.png');
saveas(gcf, 'xcorr_darkness_fwdcombined.svg');

% -P1 data
AOTU019_r_pk_fwd_nop1 = AOTU019_fr.combined_data.r_pk_fwd_nop1;
...
AOTU025_lag_pk_fwd_nop1 = AOTU025_fr.combined_data.lag_pk_fwd_nop1;
