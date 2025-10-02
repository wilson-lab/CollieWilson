% mp_aipgdifference.m
% CREATED: 09/09/2025 - MC

% This script compares the firing rate responses of AOTU019 and AOTU025 
% neurons during aIPg activation versus no aIPg activation.

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
AOTU019_P1_file = fullfile(dataPath, 'AOTU019_aIPg_Motion_Pulse_sweeppeaks_25_dps_Quiescent.mat');
AOTU019_noP1_file = fullfile(dataPath, 'AOTU019_aIPg_Motion_Pulse_No_P1_sweeppeaks_25_dps_Quiescent.mat');
AOTU025_P1_file = fullfile(dataPath, 'AOTU025_aIPg_Motion_Pulse_sweeppeaks_25_dps_Quiescent.mat');
AOTU025_noP1_file = fullfile(dataPath, 'AOTU025_aIPg_Motion_Pulse_No_P1_sweeppeaks_25_dps_Quiescent.mat');

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
% use pulse_srRightward for raw
% use pulse_srRL_Rightward for R-L
AOTU019_P1_filtered = [];
AOTU019_noP1_filtered = [];
for i = 1:length(common_AOTU019_animals)
    idx_P1 = find(strcmp(AOTU019_P1_fly_names, common_AOTU019_animals{i}));
    idx_noP1 = find(strcmp(AOTU019_noP1_fly_names, common_AOTU019_animals{i}));
    AOTU019_P1_filtered = cat(2, AOTU019_P1_filtered, AOTU019_P1_data.pulse_srRightward(:, idx_P1, :));
    AOTU019_noP1_filtered = cat(2, AOTU019_noP1_filtered, AOTU019_noP1_data.pulse_srRightward(:, idx_noP1, :));
end

AOTU025_P1_filtered = [];
AOTU025_noP1_filtered = [];
for i = 1:length(common_AOTU025_animals)
    idx_P1 = find(strcmp(AOTU025_P1_fly_names, common_AOTU025_animals{i}));
    idx_noP1 = find(strcmp(AOTU025_noP1_fly_names, common_AOTU025_animals{i}));
    AOTU025_P1_filtered = cat(2, AOTU025_P1_filtered, AOTU025_P1_data.pulse_srRightward(:, idx_P1, :));
    AOTU025_noP1_filtered = cat(2, AOTU025_noP1_filtered, AOTU025_noP1_data.pulse_srRightward(:, idx_noP1, :));
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

% Fetch n
n019 = size(AOTU019_diff_avg,1);
n025 = size(AOTU025_diff_avg,1);

%% Plot Average Differences
% Initialize plot
figure; set(gcf,'Position',[100 100 200 600])
blue = [0 0.447 0.741]; purple = [0.494 0.184 0.556]; %colors
hold on;

% Perform Independent t-test
[h, p] = ttest2(AOTU019_diff_avg, AOTU025_diff_avg);

% Set jitter amount
jitterAmount = 0.005;

% Scatter plot of average differences with jitter
scatter(1 + randn(size(AOTU019_diff_avg)) * jitterAmount, AOTU019_diff_avg, '.', 'MarkerEdgeColor', 'k');
scatter(2 + randn(size(AOTU025_diff_avg)) * jitterAmount, AOTU025_diff_avg, '.', 'MarkerEdgeColor', 'k');

% Plot the median as a dash (no jitter)
plot(1, AOTU019_median_diff, '_', 'Color', blue, 'MarkerSize', 12);
plot(2, AOTU025_median_diff, '_', 'Color', purple, 'MarkerSize', 12);

% Formatting
xlim([0.5, 2.5]);
xticks([1, 2]);
xticklabels({'AOTU019', 'AOTU025'});
ylabel('Avg Difference in Firing Rate (aIPg - no aIPg)');

% Display the p-value on the plot
text(2.5, max([AOTU019_diff_avg; AOTU025_diff_avg]+3), sprintf('p = %.5f', p), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8);

% Save
sgtitle(['MP Compare (n = ' num2str(n019) ',' num2str(n025) ')']);
cd(plotPath)
saveas(gcf, 'aipg_compare_visual_diff.png');
set(gcf, 'renderer', 'Painters'); % Save vectorized version
saveas(gcf, 'aipg_compare_visual_diff.svg');

%% Plot Average before and after P1 stimulation

% Initialize plot
figure; set(gcf,'Position',[100 100 400 600]); tiledlayout(1,2,'TileSpacing','compact');
fr_range = [0 80];

% Initilize metrics for isolating sweep
total_time_points = size(AOTU019_P1_peaks, 1);
start_idx = floor(total_time_points/4) + 1;
end_idx   = floor(3*total_time_points/4);
mid_idx   = start_idx:end_idx;

% AOTU019 - Fetch averages during sweeps
id019 = strings(0,1); no019 = []; p1019 = [];
for i = 1:size(AOTU019_P1_peaks,2)
    % Middle segment
    P1_mid   = AOTU019_P1_peaks(mid_idx,i,:);
    noP1_mid = AOTU019_noP1_peaks(mid_idx,i,:);
    % Averages
    mP1   = mean(P1_mid(:),  'omitnan');
    mNoP1 = mean(noP1_mid(:),'omitnan');

    if isfinite(mP1) && isfinite(mNoP1)
        id019(end+1,1) = sprintf('019_%02d',i);
        no019(end+1)   = mNoP1;
        p1019(end+1)   = mP1;
    end
end

% AOTU019 - Fetch averages during sweeps
id025 = strings(0,1); no025 = []; p1025 = [];
for i = 1:size(AOTU025_P1_peaks,2)
    % Middle segment
    P1_mid   = AOTU025_P1_peaks(mid_idx,i,:);
    noP1_mid = AOTU025_noP1_peaks(mid_idx,i,:);
    % Averages
    mP1   = mean(P1_mid(:),  'omitnan');
    mNoP1 = mean(noP1_mid(:),'omitnan');

    if isfinite(mP1) && isfinite(mNoP1)
        id025(end+1,1) = sprintf('025_%02d',i);
        no025(end+1)   = mNoP1;
        p1025(end+1)   = mP1;
    end
end

% AOTU019
x = [1 2];
nexttile; hold on;
for k = 1:numel(no019)
    plot(x,[no019(k) p1019(k)],'.-','Color','k','MarkerSize',8);
end
plot(1,median(no019,'omitnan'),'_','Color',blue,'MarkerSize',14,'LineWidth',1.5);
plot(2,median(p1019,'omitnan'),'_','Color',blue,'MarkerSize',14,'LineWidth',1.5);
% Formatting
xlim([0.7 2.3]); xticks([1 2]); xticklabels({'no aIPg','aIPg'}); ylim(fr_range); yline(0)
ylabel('Mean firing rate'); title('AOTU019'); box off

% AOTU025
nexttile; hold on;
for k = 1:numel(no025)
    plot(x,[no025(k) p1025(k)],'.-','Color','k','MarkerSize',8);
end
plot(1,median(no025,'omitnan'),'_','Color',purple,'MarkerSize',14,'LineWidth',1.5);
plot(2,median(p1025,'omitnan'),'_','Color',purple,'MarkerSize',14,'LineWidth',1.5);
% Formatting
xlim([0.7 2.3]); xticks([1 2]); xticklabels({'no aIPg','aIPg'}); ylim(fr_range); yline(0)
ylabel('Mean R-L firing rate'); title('AOTU025'); box off


% Repeated-measures ANOVA using mixed effects model
Resp      = [no019(:); p1019(:); no025(:); p1025(:)];
Condition = [repmat("noP1",numel(no019),1); repmat("P1",numel(p1019),1); ...
             repmat("noP1",numel(no025),1); repmat("P1",numel(p1025),1)];
CellType  = [repmat("AOTU019",numel(no019)+numel(p1019),1); ...
             repmat("AOTU025",numel(no025)+numel(p1025),1)];
FlyID     = [id019; id019; id025; id025];

tbl = table(Resp, categorical(CellType), categorical(Condition), categorical(FlyID), ...
            'VariableNames', {'Response','CellType','Condition','FlyID'});

mdl = fitlme(tbl,'Response ~ CellType*Condition + (1|FlyID)');
disp('--- Mixed-effects ANOVA (no baseline subtraction) ---');
stats = anova(mdl,'DFMethod','Satterthwaite');
disp(stats)

% Display the p-value on the plot
text(2, fr_range(2)-1, sprintf('p(cell) = %.5f\np(aIPg) = %.5f\np(x) = %.5f', stats.pValue(2), stats.pValue(3), stats.pValue(4)), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8);

% Save
sgtitle(['MP Compare (n = ' num2str(n019) ',' num2str(n025) ')']);
cd(plotPath)
saveas(gcf, 'aipg_compare_visual_raw.png');
set(gcf, 'renderer', 'Painters');
saveas(gcf, 'aipg_compare_visual_raw.svg');

%% Posthoc: paired t-tests within each cell type (Bonferroni across 2 tests)
alpha = 0.05;

% AOTU019: paired by fly (P1 vs noP1)
[~, p019, ci019, stats019] = ttest(p1019(:), no019(:));   % two-sided, paired
dz019 = (mean(p1019(:) - no019(:), 'omitnan')) / std(p1019(:) - no019(:), 'omitnan'); % Cohen's dz

% AOTU025: paired by fly (P1 vs noP1)
[~, p025, ci025, stats025] = ttest(p1025(:), no025(:));   % two-sided, paired
dz025 = (mean(p1025(:) - no025(:), 'omitnan')) / std(p1025(:) - no025(:), 'omitnan'); % Cohen's dz

% Bonferroni correction across the two planned tests
p_bonf_019 = min(p019 * 2, 1);
p_bonf_025 = min(p025 * 2, 1);

fprintf('\n--- Posthoc paired t-tests (Bonferroni across 2 tests) ---\n');
fprintf('AOTU019: t(%d) = %.3f, p = %.4g, p_bonf = %.4g, dz = %.3f, CI = [%.3f, %.3f]\n', ...
    stats019.df, stats019.tstat, p019, p_bonf_019, dz019, ci019(1), ci019(2));
fprintf('AOTU025: t(%d) = %.3f, p = %.4g, p_bonf = %.4g, dz = %.3f, CI = [%.3f, %.3f]\n', ...
    stats025.df, stats025.tstat, p025, p_bonf_025, dz025, ci025(1), ci025(2));