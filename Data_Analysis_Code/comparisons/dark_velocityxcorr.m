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
ptext = sprintf('p(Cell) = %.6f\np(Arousal) = %.6f\np(Interact) = %.6f', ...
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
ptext = sprintf('p(Cell) = %.6f\np(Arousal) = %.6f\np(Interact) = %.6f', ...
    lag_stats.pValue(2), lag_stats.pValue(3), lag_stats.pValue(4));
text(numel(groups)+0.9, 400, ptext, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

cd(plotPath)
saveas(gcf, 'xcorr_darkness_combined.png');
saveas(gcf, 'xcorr_darkness_combined.svg');

%% Post hoc: paired t-tests (effect of arousal within each cell type) for PEAK CORR
% Bonferroni across two tests (AOTU019 and AOTU025)
alpha_raw  = 0.05;
alpha_bonf = alpha_raw / 2;

% --- Collapse to a per-fly vector (if matrix, average across columns) ---
x019_no = AOTU019_r_pk_ang_nop1(:);
x019_p1 = AOTU019_r_pk_ang_p1(:);
x025_no = AOTU025_r_pk_ang_nop1(:);
x025_p1 = AOTU025_r_pk_ang_p1(:);

% Pairwise non-NaN masks to keep matched pairs only
m019 = ~isnan(x019_no) & ~isnan(x019_p1);
m025 = ~isnan(x025_no) & ~isnan(x025_p1);

% --- AOTU019 paired t-test (+P1 vs −P1) ---
[~, p019_raw, ci019, stats019] = ttest(x019_p1(m019), x019_no(m019), 'Alpha', alpha_bonf);
p019_bonf = min(p019_raw*2, 1);  % Bonferroni adjust for 2 tests
d019 = x019_p1(m019) - x019_no(m019);
dz019 = mean(d019,'omitnan') / std(d019,'omitnan');  % Cohen's dz for paired

% --- AOTU025 paired t-test (+P1 vs −P1) ---
[~, p025_raw, ci025, stats025] = ttest(x025_p1(m025), x025_no(m025), 'Alpha', alpha_bonf);
p025_bonf = min(p025_raw*2, 1);  % Bonferroni adjust
d025 = x025_p1(m025) - x025_no(m025);
dz025 = mean(d025,'omitnan') / std(d025,'omitnan');

% Print concise results
fprintf('AOTU019 (+P1 vs -P1): t(%d)=%.3f, p_raw=%.4g, p_bonf=%.4g, dz=%.3f, CI[%.3f %.3f], n=%d\n', ...
    stats019.df, stats019.tstat, p019_raw, p019_bonf, dz019, ci019(1), ci019(2), sum(m019));

fprintf('AOTU025 (+P1 vs -P1): t(%d)=%.3f, p_raw=%.4g, p_bonf=%.4g, dz=%.3f, CI[%.3f %.3f], n=%d\n', ...
    stats025.df, stats025.tstat, p025_raw, p025_bonf, dz025, ci025(1), ci025(2), sum(m025));

%% Optional confirmatory check: per-fly means, Welch t-test (AOTU019 vs AOTU025)

% Collapse each fly to its mean across arousal
% AOTU019
lag019_no = AOTU019_lag_pk_ang_nop1;  % [flies x trials] or [flies x 1]
lag019_p1 = AOTU019_lag_pk_ang_p1;
m019 = mean([lag019_no, lag019_p1], 2, 'omitnan');  % per-fly mean across arousal

% AOTU025
lag025_no = AOTU025_lag_pk_ang_nop1;
lag025_p1 = AOTU025_lag_pk_ang_p1;
m025 = mean([lag025_no, lag025_p1], 2, 'omitnan');

% Remove NaNs
m019 = m019(~isnan(m019));
m025 = m025(~isnan(m025));

% Welch two-sample t-test
[~, p_welch, ~, stats_w] = ttest2(m025, m019, 'Vartype','unequal');

% Hedges' g for unequal n, unbiased
s019 = var(m019, 'omitnan'); s025 = var(m025, 'omitnan');
n019 = numel(m019);          n025 = numel(m025);
sp   = sqrt(((n019-1)*s019 + (n025-1)*s025) / (n019 + n025 - 2));
g    = (mean(m025)-mean(m019)) / sp;
J    = 1 - (3/(4*(n019+n025)-9));  % small-sample correction
g_unbiased = g * J;

fprintf('Welch t-test on per-fly mean lag (AOTU025 - AOTU019): t(%0.1f)=%.3f, p=%.4g, g=%.3f\n', ...
    stats_w.df, stats_w.tstat, p_welch, g_unbiased);

%% XCORR (dark): Paired no P1 vs +P1 per fly — r and lag (filter lag>0)
close all

% Colors
blue   = [0 0.35 0.90];
purple = [0.55 0 0.65];

% Axis ranges
r_range   = [-0.05 0.65];
lag_range = [-400 0];

figure; set(gcf,'Position',[100 100 900 900])
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% ---------------- AOTU019 ----------------
% Filter (valid + lag<=0 both conditions)
noP1_r = AOTU019_r_pk_ang_nop1(:);
p1_r   = AOTU019_r_pk_ang_p1(:);
noP1_l = AOTU019_lag_pk_ang_nop1(:);
p1_l   = AOTU019_lag_pk_ang_p1(:);
valid019 = isfinite(noP1_r) & isfinite(p1_r) & ...
           isfinite(noP1_l) & isfinite(p1_l) & ...
           (noP1_l <= 0) & (p1_l <= 0);
noP1_r = noP1_r(valid019); p1_r = p1_r(valid019);
noP1_l = noP1_l(valid019); p1_l = p1_l(valid019);

% r panel
ax = nexttile; hold on
for k = 1:numel(noP1_r)
    plot([1 2],[noP1_r(k) p1_r(k)],'.-','Color',[0 0 0],'MarkerSize',10);
end
plot(1,median(noP1_r,'omitnan'),'_','Color',blue,'MarkerSize',16,'LineWidth',1.8);
plot(2,median(p1_r,'omitnan'),'_','Color',blue,'MarkerSize',16,'LineWidth',1.8);
xlim([0.7 2.3]); ylim(r_range); yline(0,'-');
xticks([1 2]); xticklabels({'-P1','+P1'});
ylabel('Peak correlation (r)'); title('AOTU019'); box off

% lag panel
ax = nexttile; hold on
for k = 1:numel(noP1_l)
    plot([1 2],[noP1_l(k) p1_l(k)],'.-','Color',[0 0 0],'MarkerSize',10);
end
plot(1,median(noP1_l,'omitnan'),'_','Color',blue,'MarkerSize',16,'LineWidth',1.8);
plot(2,median(p1_l,'omitnan'),'_','Color',blue,'MarkerSize',16,'LineWidth',1.8);
xlim([0.7 2.3]); ylim(lag_range); yline(0,'-');
xticks([1 2]); xticklabels({'-P1','+P1'});
ylabel('Lag (ms)'); title('AOTU019'); box off

% ---------------- AOTU025 ----------------
noP1_r = AOTU025_r_pk_ang_nop1(:);
p1_r   = AOTU025_r_pk_ang_p1(:);
noP1_l = AOTU025_lag_pk_ang_nop1(:);
p1_l   = AOTU025_lag_pk_ang_p1(:);
valid025 = isfinite(noP1_r) & isfinite(p1_r) & ...
           isfinite(noP1_l) & isfinite(p1_l) & ...
           (noP1_l <= 0) & (p1_l <= 0);
noP1_r = noP1_r(valid025); p1_r = p1_r(valid025);
noP1_l = noP1_l(valid025); p1_l = p1_l(valid025);

% r panel
ax = nexttile; hold on
for k = 1:numel(noP1_r)
    plot([1 2],[noP1_r(k) p1_r(k)],'.-','Color',[0 0 0],'MarkerSize',10);
end
plot(1,median(noP1_r,'omitnan'),'_','Color',purple,'MarkerSize',16,'LineWidth',1.8);
plot(2,median(p1_r,'omitnan'),'_','Color',purple,'MarkerSize',16,'LineWidth',1.8);
xlim([0.7 2.3]); ylim(r_range); yline(0,'-');
xticks([1 2]); xticklabels({'-P1','+P1'});
ylabel('Peak correlation (r)'); title('AOTU025'); box off

% lag panel
ax = nexttile; hold on
for k = 1:numel(noP1_l)
    plot([1 2],[noP1_l(k) p1_l(k)],'.-','Color',[0 0 0],'MarkerSize',10);
end
plot(1,median(noP1_l,'omitnan'),'_','Color',purple,'MarkerSize',16,'LineWidth',1.8);
plot(2,median(p1_l,'omitnan'),'_','Color',purple,'MarkerSize',16,'LineWidth',1.8);
xlim([0.7 2.3]); ylim(lag_range); yline(0,'-');
xticks([1 2]); xticklabels({'-P1','+P1'});
ylabel('Lag (ms)'); title('AOTU025'); box off

sgtitle('XCORR (dark): Paired no P1 vs +P1 per fly (lag>0 excluded)','FontWeight','bold');

% Save (if plotPath exists)
if exist('plotPath','var') && isfolder(plotPath)
    cd(plotPath)
    saveas(gcf, 'xcorr_darkness_paired_filtered.png');
    set(gcf, 'renderer', 'Painters');
    saveas(gcf, 'xcorr_darkness_paired_filtered.svg');
end

% ---------------- Display valid flies ----------------
fprintf('Valid flies included (lag <= 0 in both conditions):\n');
fprintf('  AOTU019: %d\n', sum(valid019));
fprintf('  AOTU025: %d\n', sum(valid025));
