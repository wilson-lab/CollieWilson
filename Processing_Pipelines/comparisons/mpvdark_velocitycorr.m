% mpvdark_velocitycorr (updated version without -P1 darkness)

%% Load cross correlation data
clear; close all

plotPath = 'E:\\Compare Motion Pulse';
dataPath = 'E:\\Compare Motion Pulse\\data';

disp('Analyzing xcorr in darkness and visual pursuit.');
cd(dataPath)

% Load files
AOTU019_dark = load(fullfile(dataPath, 'AOTU019_Background_P1_xcorr.mat'));
AOTU025_dark = load(fullfile(dataPath, 'AOTU025_Background_P1_xcorr.mat'));
AOTU019_visual = load(fullfile(dataPath, 'AOTU019_Motion_Pulse_xcorr.mat'));
AOTU025_visual = load(fullfile(dataPath, 'AOTU025_Motion_Pulse_xcorr.mat'));

% Data
AOTU019_angr_dark_p1 = AOTU019_dark.combined_data.r_pk_ang;
AOTU025_angr_dark_p1 = AOTU025_dark.combined_data.r_pk_ang;
AOTU019_angr_visual_p1 = AOTU019_visual.rotational_peaks.r_pk;
AOTU025_angr_visual_p1 = AOTU025_visual.rotational_peaks.r_pk;

AOTU019_anglag_dark_p1 = AOTU019_dark.combined_data.lag_pk_ang;
AOTU025_anglag_dark_p1 = AOTU025_dark.combined_data.lag_pk_ang;
AOTU019_anglag_visual_p1 = AOTU019_visual.rotational_peaks.lag_pk;
AOTU025_anglag_visual_p1 = AOTU025_visual.rotational_peaks.lag_pk;

% Names
AOTU019_dark_names = AOTU019_dark.storeNames;
AOTU025_dark_names = AOTU025_dark.storeNames;
AOTU019_visual_names = AOTU019_visual.storeNames;
AOTU025_visual_names = AOTU025_visual.storeNames;

% Build metadata
peak_corr_values = [AOTU019_angr_dark_p1(:);
                    AOTU025_angr_dark_p1(:);
                    AOTU019_angr_visual_p1(:);
                    AOTU025_angr_visual_p1(:)];

lag_values = [AOTU019_anglag_dark_p1(:);
              AOTU025_anglag_dark_p1(:);
              AOTU019_anglag_visual_p1(:);
              AOTU025_anglag_visual_p1(:)];

% Sizes
n019d = size(AOTU019_angr_dark_p1,1);
n025d = size(AOTU025_angr_dark_p1,1);
n019v = size(AOTU019_angr_visual_p1,1);
n025v = size(AOTU025_angr_visual_p1,1);

rep019d = size(AOTU019_angr_dark_p1,2);
rep025d = size(AOTU025_angr_dark_p1,2);
rep019v = size(AOTU019_angr_visual_p1,2);
rep025v = size(AOTU025_angr_visual_p1,2);

cell_type = [repmat("AOTU019", n019d*rep019d,1);
             repmat("AOTU025", n025d*rep025d,1);
             repmat("AOTU019", n019v*rep019v,1);
             repmat("AOTU025", n025v*rep025v,1)];

condition = [repmat("+P1_Dark", n019d*rep019d,1);
             repmat("+P1_Dark", n025d*rep025d,1);
             repmat("+P1_Visual", n019v*rep019v,1);
             repmat("+P1_Visual", n025v*rep025v,1)];

fly_id = [repmat(AOTU019_dark_names(:), rep019d, 1);
          repmat(AOTU025_dark_names(:), rep025d, 1);
          repmat(AOTU019_visual_names(:), rep019v, 1);
          repmat(AOTU025_visual_names(:), rep025v, 1)];

% Remove NaNs
valid_peak = ~isnan(peak_corr_values);
valid_lag = ~isnan(lag_values);

peak_corr_values = peak_corr_values(valid_peak);
cell_type_peak = cell_type(valid_peak);
condition_peak = condition(valid_peak);
fly_id_peak = fly_id(valid_peak);

lag_values = lag_values(valid_lag);
cell_type_lag = cell_type(valid_lag);
condition_lag = condition(valid_lag);
fly_id_lag = fly_id(valid_lag);

group_peak = strcat(cell_type_peak, "_", condition_peak);
group_lag = strcat(cell_type_lag, "_", condition_lag);

% ---------------- Run LME ----------------
tbl_peak = table(peak_corr_values, cell_type_peak, condition_peak, fly_id_peak, ...
    'VariableNames', {'Peak', 'Cell', 'Condition', 'Fly'});
lme_peak = fitlme(tbl_peak, 'Peak ~ Cell*Condition + (1|Fly)');
peak_stats = anova(lme_peak);

tbl_lag = table(lag_values, cell_type_lag, condition_lag, fly_id_lag, ...
    'VariableNames', {'Lag', 'Cell', 'Condition', 'Fly'});
lme_lag = fitlme(tbl_lag, 'Lag ~ Cell*Condition + (1|Fly)');
lag_stats = anova(lme_lag);

% ---------------- Post-hoc (condition effect only) ----------------
[~,~,stats] = anovan(peak_corr_values, {condition_peak}, 'display', 'off');
figure;
[c,m,h,gnames] = multcompare(stats, 'ctype', 'tukey-kramer');
title('Post-hoc comparison for condition (Peak Corr)');
disp('All pairwise differences (Tukey):');
disp(array2table(c, ...
    'VariableNames', {'Group1','Group2','LowerCI','Estimate','UpperCI','pValue'}));

% ---------------- Plot ----------------
figure;
tiledlayout(1, 2);
jitter = 0.1;

colors = containers.Map({'AOTU019','AOTU025'}, {[0 0.4470 0.7410],[0.4940 0.1840 0.5560]});
order = ["AOTU019_+P1_Dark", "AOTU019_+P1_Visual", "AOTU025_+P1_Dark", "AOTU025_+P1_Visual"];
order_clean = strrep(order, '_', ' ');

% PEAK
nexttile;
hold on;
for i = 1:numel(order)
    idx = strcmp(group_peak, order(i));
    ctype = extractBefore(order(i), '_');
    x = i + jitter*(rand(sum(idx),1)-0.5);
    scatter(x, peak_corr_values(idx), '.', 'MarkerEdgeColor', colors(ctype));
    plot(i, median(peak_corr_values(idx),'omitnan'), '_', 'MarkerSize', 15, 'Color', colors(ctype), 'LineWidth', 2);
end
xlim([0 numel(order)+1]); ylim([0 0.7]);
xticks(1:numel(order)); xticklabels(order_clean); xtickangle(30);
ylabel('Correlation'); title('Peak Correlation');
ptext = sprintf(['p(Cell) = %.3f\np(Cond) = %.3f\np(Intxn) = %.3f'], ...
    peak_stats.pValue(2:4));
text(numel(order)+0.9, 0.65, ptext, 'HorizontalAlignment','right','VerticalAlignment','top');

% LAG
nexttile;
hold on;
for i = 1:numel(order)
    idx = strcmp(group_lag, order(i));
    ctype = extractBefore(order(i), '_');
    x = i + jitter*(rand(sum(idx),1)-0.5);
    scatter(x, lag_values(idx), '.', 'MarkerEdgeColor', colors(ctype));
    plot(i, median(lag_values(idx),'omitnan'), '_', 'MarkerSize', 15, 'Color', colors(ctype), 'LineWidth', 2);
end
xlim([0 numel(order)+1]); ylim([-350 0]);
xticks(1:numel(order)); xticklabels(order_clean); xtickangle(30);
ylabel('Lag (ms)'); title('Lag Estimate');
ptext = sprintf(['p(Cell) = %.3f\np(Cond) = %.3f\np(Intxn) = %.3f'], ...
    lag_stats.pValue(2:4));
text(numel(order)+0.9, -50, ptext, 'HorizontalAlignment','right','VerticalAlignment','top');

cd(plotPath)
saveas(gcf, 'xcorr_darkness_v_mp_combined.png');
saveas(gcf, 'xcorr_darkness_v_mp_combined.svg');

% ---------------- Count non-NaN entries per cell type and condition (Peak Correlation) ----------------
fprintf('\n--- Non-NaN entries per Cell Type and Condition (Peak Correlation) ---\n');
celltypes = unique(cell_type_peak);
conditions = unique(condition_peak);

for i = 1:numel(celltypes)
    for j = 1:numel(conditions)
        idx = strcmp(cell_type_peak, celltypes(i)) & strcmp(condition_peak, conditions(j));
        count = sum(idx);
        fprintf('%s, %s: %d\n', celltypes(i), conditions(j), count);
    end
end


