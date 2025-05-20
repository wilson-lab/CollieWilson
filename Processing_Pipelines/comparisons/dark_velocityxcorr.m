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
AOTU019_r_pk_ang_p1 = AOTU019_fr.combined_data.r_pk_ang;
AOTU019_lag_pk_ang_p1 = AOTU019_fr.combined_data.lag_pk_ang;
AOTU025_r_pk_ang_p1 = AOTU025_fr.combined_data.r_pk_ang;
AOTU025_lag_pk_ang_p1 = AOTU025_fr.combined_data.lag_pk_ang;

% -P1 data
AOTU019_r_pk_ang_nop1 = AOTU019_fr.combined_data.r_pk_ang_nop1;
AOTU019_lag_pk_ang_nop1 = AOTU019_fr.combined_data.lag_pk_ang_nop1;
AOTU025_r_pk_ang_nop1 = AOTU025_fr.combined_data.r_pk_ang_nop1;
AOTU025_lag_pk_ang_nop1 = AOTU025_fr.combined_data.lag_pk_ang_nop1;

%% Plot Combined
mean_colors = {[0 0 1], [0.5 0 0.5]};
jitter_amount = 0.1;

figure; set(gcf, 'Position', [100 100 600 800])
tiledlayout(2, 2, 'TileSpacing', 'compact')

% Tile 1: +P1 Peak Correlation
nexttile;
scatter(1 + jitter_amount * (rand(size(AOTU019_r_pk_ang_p1)) - 0.5), AOTU019_r_pk_ang_p1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 + jitter_amount * (rand(size(AOTU025_r_pk_ang_p1)) - 0.5), AOTU025_r_pk_ang_p1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
plot(1, median(AOTU019_r_pk_ang_p1, 'omitnan'), '_', 'Color', mean_colors{1}, 'MarkerSize', 15, 'LineWidth', 2);
plot(2, median(AOTU025_r_pk_ang_p1, 'omitnan'), '_', 'Color', mean_colors{2}, 'MarkerSize', 15, 'LineWidth', 2);
title('Peak Correlation (+P1)');
xticks([1 2]); xticklabels({'AOTU019','AOTU025'});
ylabel('Correlation'); xlim([0 3]); ylim([0 0.7]); yline(0);
[~, p1] = ttest2(AOTU019_r_pk_ang_p1, AOTU025_r_pk_ang_p1);
text(2.8, 0.65, sprintf('p = %.3f', p1), 'HorizontalAlignment', 'right');

% Tile 2: +P1 Lag
nexttile;
scatter(1 + jitter_amount * (rand(size(AOTU019_lag_pk_ang_p1)) - 0.5), AOTU019_lag_pk_ang_p1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 + jitter_amount * (rand(size(AOTU025_lag_pk_ang_p1)) - 0.5), AOTU025_lag_pk_ang_p1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
plot(1, median(AOTU019_lag_pk_ang_p1, 'omitnan'), '_', 'Color', mean_colors{1}, 'MarkerSize', 15, 'LineWidth', 2);
plot(2, median(AOTU025_lag_pk_ang_p1, 'omitnan'), '_', 'Color', mean_colors{2}, 'MarkerSize', 15, 'LineWidth', 2);
title('Lag Estimate (+P1)');
xticks([1 2]); xticklabels({'AOTU019','AOTU025'});
ylabel('Lag (ms)'); xlim([0 3]); ylim([-300 0]); yline(0);
[~, p2] = ttest2(AOTU019_lag_pk_ang_p1, AOTU025_lag_pk_ang_p1);
text(2.8, -50, sprintf('p = %.3f', p2), 'HorizontalAlignment', 'right');

% Tile 3: -P1 Peak Correlation
nexttile;
scatter(1 + jitter_amount * (rand(size(AOTU019_r_pk_ang_nop1)) - 0.5), AOTU019_r_pk_ang_nop1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 + jitter_amount * (rand(size(AOTU025_r_pk_ang_nop1)) - 0.5), AOTU025_r_pk_ang_nop1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
plot(1, median(AOTU019_r_pk_ang_nop1, 'omitnan'), '_', 'Color', mean_colors{1}, 'MarkerSize', 15, 'LineWidth', 2);
plot(2, median(AOTU025_r_pk_ang_nop1, 'omitnan'), '_', 'Color', mean_colors{2}, 'MarkerSize', 15, 'LineWidth', 2);
title('Peak Correlation (-P1)');
xticks([1 2]); xticklabels({'AOTU019','AOTU025'});
ylabel('Correlation'); xlim([0 3]); ylim([0 0.7]); yline(0);
[~, p3] = ttest2(AOTU019_r_pk_ang_nop1, AOTU025_r_pk_ang_nop1);
text(2.8, 0.65, sprintf('p = %.3f', p3), 'HorizontalAlignment', 'right');

% Tile 4: -P1 Lag
nexttile;
scatter(1 + jitter_amount * (rand(size(AOTU019_lag_pk_ang_nop1)) - 0.5), AOTU019_lag_pk_ang_nop1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
scatter(2 + jitter_amount * (rand(size(AOTU025_lag_pk_ang_nop1)) - 0.5), AOTU025_lag_pk_ang_nop1, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
plot(1, median(AOTU019_lag_pk_ang_nop1, 'omitnan'), '_', 'Color', mean_colors{1}, 'MarkerSize', 15, 'LineWidth', 2);
plot(2, median(AOTU025_lag_pk_ang_nop1, 'omitnan'), '_', 'Color', mean_colors{2}, 'MarkerSize', 15, 'LineWidth', 2);
title('Lag Estimate (-P1)');
xticks([1 2]); xticklabels({'AOTU019','AOTU025'});
ylabel('Lag (ms)'); xlim([0 3]); ylim([-300 0]); yline(0);
[~, p4] = ttest2(AOTU019_lag_pk_ang_nop1, AOTU025_lag_pk_ang_nop1);
text(2.8, -50, sprintf('p = %.3f', p4), 'HorizontalAlignment', 'right');

sgtitle('XCORR: +P1 (top) vs -P1 (bottom)');
cd(plotPath)
saveas(gcf, 'xcorr_darkness.png');
saveas(gcf, 'xcorr_darkness.svg');

%% Full Comparison Plot and ANOVA

% Prepare data vectors and metadata for ANOVA
peak_corr_values = [AOTU019_r_pk_ang_nop1(:); AOTU019_r_pk_ang_p1(:); AOTU025_r_pk_ang_nop1(:); AOTU025_r_pk_ang_p1(:)];
lag_values = [AOTU019_lag_pk_ang_nop1(:); AOTU019_lag_pk_ang_p1(:); AOTU025_lag_pk_ang_nop1(:); AOTU025_lag_pk_ang_p1(:)];

n_019 = size(AOTU019_r_pk_ang_nop1, 1);
n_025 = size(AOTU025_r_pk_ang_nop1, 1);

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
saveas(gcf, 'xcorr_darkness_combined.png');
saveas(gcf, 'xcorr_darkness_combined.svg');

%% Posthoc
% Get coefficient names
disp(lme_peak.CoefficientNames);

% Contrast vector: +P1 vs -P1 across both cell types
% This assumes the arousal effect is in 'Arousal_+P1'
C_arousal = zeros(1, length(lme_peak.CoefficientNames));
C_arousal(strcmp(lme_peak.CoefficientNames, 'Arousal_+P1')) = 1;

[p_arousal, Fval, df1, df2] = coefTest(lme_peak, C_arousal);
fprintf('Post hoc: +P1 vs -P1 (Peak Correlation): p = %.4e\n', p_arousal);

% Get coefficient names
disp(lme_lag.CoefficientNames);

% Contrast vector: AOTU025 vs AOTU019
% This assumes the cell type effect is in 'Cell_AOTU025'
C_cell = zeros(1, length(lme_lag.CoefficientNames));
C_cell(strcmp(lme_lag.CoefficientNames, 'Cell_AOTU025')) = 1;

[p_cell, Fval, df1, df2] = coefTest(lme_lag, C_cell);
fprintf('Post hoc: AOTU025 vs AOTU019 (Lag): p = %.4f\n', p_cell);

%% Display number of valid flies per cell type
valid_019_ids = unique(fly_id_peak(strcmp(cell_type_peak, "AOTU019")));
valid_025_ids = unique(fly_id_peak(strcmp(cell_type_peak, "AOTU025")));

fprintf('Number of valid flies:\n');
fprintf('  AOTU019: %d\n', numel(valid_019_ids));
fprintf('  AOTU025: %d\n', numel(valid_025_ids));
