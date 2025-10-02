% AOTU Spike and Voltage Difference Analysis
% This script loads firing rate and voltage differences from AOTU019 and 
% AOTU025 data and performs comparisons across all and quiescent-only 
% timepoints, running t-tests to check significance between the two cell types.
%
% CREATED: 09/03/2025 - MC adapted from p1 analysis

% Load AOTU019 and AOTU025 spike and voltage difference data
clear; close all
plotPath = 'E:\Compare Motion Pulse';  % Define folder for saving plots
dataPath = 'E:\Compare Motion Pulse\data';  % Define folder containing the .mat files
cd(dataPath);

% Load spike and voltage difference data for AOTU019 and AOTU025
AOTU019_data = load(fullfile(dataPath, 'AOTU019_Background_aIPg_spike_voltage_diffs.mat'));
AOTU025_data = load(fullfile(dataPath, 'AOTU025_Background_aIPg_spike_voltage_diffs.mat'));

% Access difference data
AOTU019_diff_sr = AOTU019_data.diff_sr;
AOTU019_diff_vm = AOTU019_data.diff_vm;
AOTU025_diff_sr = AOTU025_data.diff_sr;
AOTU025_diff_vm = AOTU025_data.diff_vm;

% Fetch n
n019 = size(AOTU019_diff_sr,1);
n025 = size(AOTU025_diff_sr,1);

%% Plot
% Create tiled layout for plotting
figure; set(gcf, 'Position', [100 100 300 800]);
tiledlayout(2, 2, 'TileSpacing', 'compact');

% Define jitter and color settings
jitter_amount = 0.025;  % Reduced jitter
scatterColor = [0.5 0.5 0.5];  % Grey for scatter points
color019 = [0 0 1];    % Blue for AOTU019
color025 = [0.5 0 0.5]; % Purple for AOTU025
text_size = 8;  % Slightly smaller text size for p-values

vm_range = [-1.5 5.5];
fr_range = [-6 58];

% Plot 1: Firing Rate Difference - All Timepoints
nexttile;
scatter(ones(size(AOTU019_diff_sr(:,1))) + jitter_amount * randn(size(AOTU019_diff_sr(:,1))), AOTU019_diff_sr(:,1), 'MarkerEdgeColor', scatterColor, 'Marker', '.');
hold on;
scatter(2 * ones(size(AOTU025_diff_sr(:,1))) + jitter_amount * randn(size(AOTU025_diff_sr(:,1))), AOTU025_diff_sr(:,1), 'MarkerEdgeColor', scatterColor, 'Marker', '.');
plot(1, median(AOTU019_diff_sr(:,1), 'omitnan'), 'LineStyle', 'none', 'Marker', '_', 'MarkerSize', 10, 'Color', color019, 'LineWidth', 2);
plot(2, median(AOTU025_diff_sr(:,1), 'omitnan'), 'LineStyle', 'none', 'Marker', '_', 'MarkerSize', 10, 'Color', color025, 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 1); % Add solid black y-line at 0
title('(All)');
xticks([1 2]); xticklabels({'AOTU019', 'AOTU025'});
ylabel('Firing Rate Difference');
xlim([0 3]); % Extend x-axis padding
ylim(fr_range); % Set y-axis limits for firing rate difference
[~, p_fr_all] = ttest2(AOTU019_diff_sr(:,1), AOTU025_diff_sr(:,1), 'Vartype', 'unequal');
text(0.5, -15, sprintf('p = %.5e', p_fr_all), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', text_size);

% Plot 2: Firing Rate Difference - Quiescent Only
nexttile;
scatter(ones(size(AOTU019_diff_sr(:,2))) + jitter_amount * randn(size(AOTU019_diff_sr(:,2))), AOTU019_diff_sr(:,2), 'MarkerEdgeColor', scatterColor, 'Marker', '.');
hold on;
scatter(2 * ones(size(AOTU025_diff_sr(:,2))) + jitter_amount * randn(size(AOTU025_diff_sr(:,2))), AOTU025_diff_sr(:,2), 'MarkerEdgeColor', scatterColor, 'Marker', '.');
plot(1, median(AOTU019_diff_sr(:,2), 'omitnan'), 'LineStyle', 'none', 'Marker', '_', 'MarkerSize', 10, 'Color', color019, 'LineWidth', 2);
plot(2, median(AOTU025_diff_sr(:,2), 'omitnan'), 'LineStyle', 'none', 'Marker', '_', 'MarkerSize', 10, 'Color', color025, 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 1); % Add solid black y-line at 0
title('(Quiet)');
xticks([1 2]); xticklabels({'AOTU019', 'AOTU025'});
ylabel('Firing Rate Difference');
xlim([0 3]); % Extend x-axis padding
ylim(fr_range); % Set y-axis limits for firing rate difference
[~, p_fr_quiet] = ttest2(AOTU019_diff_sr(:,2), AOTU025_diff_sr(:,2), 'Vartype', 'unequal');
text(0.5, -15, sprintf('p = %.5e', p_fr_quiet), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', text_size);

% Plot 3: Voltage Difference - All Timepoints
nexttile;
scatter(ones(size(AOTU019_diff_vm(:,1))) + jitter_amount * randn(size(AOTU019_diff_vm(:,1))), AOTU019_diff_vm(:,1), 'MarkerEdgeColor', scatterColor, 'Marker', '.');
hold on;
scatter(2 * ones(size(AOTU025_diff_vm(:,1))) + jitter_amount * randn(size(AOTU025_diff_vm(:,1))), AOTU025_diff_vm(:,1), 'MarkerEdgeColor', scatterColor, 'Marker', '.');
plot(1, median(AOTU019_diff_vm(:,1), 'omitnan'), 'LineStyle', 'none', 'Marker', '_', 'MarkerSize', 10, 'Color', color019, 'LineWidth', 2);
plot(2, median(AOTU025_diff_vm(:,1), 'omitnan'), 'LineStyle', 'none', 'Marker', '_', 'MarkerSize', 10, 'Color', color025, 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 1); % Add solid black y-line at 0
title('(All)');
xticks([1 2]); xticklabels({'AOTU019', 'AOTU025'});
ylabel('Voltage Difference');
xlim([0 3]); % Extend x-axis padding
ylim(vm_range); % Set y-axis limits for voltage difference
[~, p_vm_all] = ttest2(AOTU019_diff_vm(:,1), AOTU025_diff_vm(:,1), 'Vartype', 'unequal');
text(0.5, -1.2, sprintf('p = %.5e', p_vm_all), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', text_size);

% Plot 4: Voltage Difference - Quiescent Only
nexttile;
scatter(ones(size(AOTU019_diff_vm(:,2))) + jitter_amount * randn(size(AOTU019_diff_vm(:,2))), AOTU019_diff_vm(:,2), 'MarkerEdgeColor', scatterColor, 'Marker', '.');
hold on;
scatter(2 * ones(size(AOTU025_diff_vm(:,2))) + jitter_amount * randn(size(AOTU025_diff_vm(:,2))), AOTU025_diff_vm(:,2), 'MarkerEdgeColor', scatterColor, 'Marker', '.');
plot(1, median(AOTU019_diff_vm(:,2), 'omitnan'), 'LineStyle', 'none', 'Marker', '_', 'MarkerSize', 10, 'Color', color019, 'LineWidth', 2);
plot(2, median(AOTU025_diff_vm(:,2), 'omitnan'), 'LineStyle', 'none', 'Marker', '_', 'MarkerSize', 10, 'Color', color025, 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 1); % Add solid black y-line at 0
title('(Quiet)');
xticks([1 2]); xticklabels({'AOTU019', 'AOTU025'});
ylabel('Voltage Difference');
xlim([0 3]); % Extend x-axis padding
ylim(vm_range); % Set y-axis limits for voltage difference
[~, p_vm_quiet] = ttest2(AOTU019_diff_vm(:,2), AOTU025_diff_vm(:,2), 'Vartype', 'unequal');
text(0.5, -1.2, sprintf('p = %.5e', p_vm_quiet), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', text_size);

% Add main title
sgtitle(['Dark Compare (n = ' num2str(n019) ',' num2str(n025) ')']);
% Save plots as images
cd(plotPath)
saveas(gcf, 'aipg_compare_darkness.png');
set(gcf, 'renderer', 'Painters'); % Save vectorized version
saveas(gcf, 'aipg_compare_darkness.svg');


%% Percent Time Running – Arousal x Cell Type ANOVA
% Load percent run time data
cd(dataPath);
data019 = load('AOTU019_Background_aIPg_runpercent_data.mat');
data025 = load('AOTU025_Background_aIPg_runpercent_data.mat');

run019_off = data019.runpercent_off(:);
run019_on  = data019.runpercent_on(:);
run025_off = data025.runpercent_off(:);
run025_on  = data025.runpercent_on(:);

% Stack data for ANOVA
run_all = [run019_off; run019_on; run025_off; run025_on];
nFlies_019 = numel(run019_off);
nFlies_025 = numel(run025_off);

celltype = [repmat({'AOTU019'}, nFlies_019 * 2, 1); repmat({'AOTU025'}, nFlies_025 * 2, 1)];
arousal  = [repmat({'off'}, nFlies_019, 1); repmat({'on'}, nFlies_019, 1); ...
            repmat({'off'}, nFlies_025, 1); repmat({'on'}, nFlies_025, 1)];
T = table(run_all, celltype, arousal, 'VariableNames', {'RunTime', 'CellType', 'Arousal'});

% Run 2-way ANOVA
[p, tbl, stats] = anovan(T.RunTime, {T.CellType, T.Arousal}, ...
    'model', 2, 'varnames', {'Cell Type', 'Arousal'}, 'display', 'off');

% Plot run time with arousal state and cell type
figure; set(gcf, 'Position', [100 100 300 800]);
tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

% X positions
x019 = [1, 2];
x025 = [3, 4];

% Plot individual data with connecting lines
nexttile; hold on;
for i = 1:nFlies_019
    plot(x019, [run019_off(i), run019_on(i)], '-', 'Color', [0.7 0.7 0.7])
end
for i = 1:nFlies_025
    plot(x025, [run025_off(i), run025_on(i)], '-', 'Color', [0.7 0.7 0.7])
end

% Scatter points
scatter(repmat(1, nFlies_019, 1), run019_off, 20, 'k', '.');
scatter(repmat(2, nFlies_019, 1), run019_on, 20, 'k', '.');
scatter(repmat(3, nFlies_025, 1), run025_off, 20, 'k', '.');
scatter(repmat(4, nFlies_025, 1), run025_on, 20, 'k', '.');

% Medians
plot(1, median(run019_off, 'omitnan'), 'k_', 'MarkerSize', 12, 'LineWidth', 1.5);
plot(2, median(run019_on, 'omitnan'), 'r_', 'MarkerSize', 12, 'LineWidth', 1.5);
plot(3, median(run025_off, 'omitnan'), 'k_', 'MarkerSize', 12, 'LineWidth', 1.5);
plot(4, median(run025_on, 'omitnan'), 'r_', 'MarkerSize', 12, 'LineWidth', 1.5);

% Formatting
xlim([0.5 4.5]);
xticks([1 2 3 4]);
xticklabels({'019 off', '019 on', '025 off', '025 on'});
ylabel('% Time Running');
title('Arousal x Cell Type');
yl = ylim;
xl = xlim;
text(xl(2) - 0.1, yl(2) - 1, ...
    {sprintf('p(cell type) = %.3e', p(1)), sprintf('p(arousal) = %.3g', p(2))}, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8);

% Second tile: Change in percent run time (P1 on - off)
delta019 = run019_on - run019_off;
delta025 = run025_on - run025_off;

nexttile; hold on;

% Jittered x-values
x_jitter_019 = 1 + jitter_amount * randn(nFlies_019, 1);
x_jitter_025 = 2 + jitter_amount * randn(nFlies_025, 1);

% Plot individual differences with jitter
scatter(x_jitter_019, delta019, 20, [0 0 1], '.');         % AOTU019 - Blue
scatter(x_jitter_025, delta025, 20, [0.5 0 0.5], '.');     % AOTU025 - Purple

% Median bars
plot(1, median(delta019, 'omitnan'), '_', 'Color', [0 0 1], 'MarkerSize', 12, 'LineWidth', 2);
plot(2, median(delta025, 'omitnan'), '_', 'Color', [0.5 0 0.5], 'MarkerSize', 12, 'LineWidth', 2);

% Stats
[~, p_delta] = ttest2(delta019, delta025);

% Formatting
xlim([0.5 2.5]);
ylim([0 70]);
xticks([1 2]);
xticklabels({'AOTU019', 'AOTU025'});
ylabel('\Delta % Time Running (aIPg on - off)');
title('Change in Run Time with Arousal');
yl = ylim;
xl = xlim;
text(xl(2) - 0.05, yl(2) - 1, sprintf('p = %.3g', p_delta), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8);


% Save plot
cd(plotPath);
saveas(gcf, 'aipg_run_arousal_celltype.png');
set(gcf, 'renderer', 'Painters');
saveas(gcf, 'aipg_run_arousal_celltype.svg');
