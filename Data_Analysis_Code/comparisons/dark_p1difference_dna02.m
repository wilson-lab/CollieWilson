clear; close all
plotPath = 'E:\Compare Motion Pulse';  % Define folder for saving plots
dataPath = 'E:\Compare Motion Pulse\data';  % Define folder containing the .mat files
cd(dataPath);

% Load spike and voltage difference data for AOTU019 and AOTU025
DNa02_data = load(fullfile(dataPath, 'DNa02_Background_P1_spike_voltage_diffs.mat'));

% Access difference data
DNa02_diff_sr = DNa02_data.diff_sr;
DNa02_diff_vm = DNa02_data.diff_vm;

nFly = size(DNa02_diff_sr,1);

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
fr_range = [-6 58]-20;

% Plot 1: Firing Rate Difference - All Timepoints
nexttile;
scatter(ones(size(DNa02_diff_sr(:,1))) + jitter_amount * randn(size(DNa02_diff_sr(:,1))), DNa02_diff_sr(:,1), 'MarkerEdgeColor', scatterColor, 'Marker', '.');
hold on;
plot(1, median(DNa02_diff_sr(:,1), 'omitnan'), 'LineStyle', 'none', 'Marker', '_', 'MarkerSize', 10, 'Color', color019, 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 1); % Add solid black y-line at 0
title('(All)');
xticks(1); xticklabels({'DNa02'});
ylabel('Firing Rate Difference');
xlim([0 3]); % Extend x-axis padding
ylim(fr_range); % Set y-axis limits for firing rate difference
[~, p_fr_all] = ttest(DNa02_diff_sr(:,1));
text(0.5, -15, sprintf('p = %.5f', p_fr_all), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', text_size);

% Plot 2: Firing Rate Difference - Quiescent Only
nexttile;
scatter(ones(size(DNa02_diff_sr(:,2))) + jitter_amount * randn(size(DNa02_diff_sr(:,2))), DNa02_diff_sr(:,2), 'MarkerEdgeColor', scatterColor, 'Marker', '.');
hold on;
plot(1, median(DNa02_diff_sr(:,2), 'omitnan'), 'LineStyle', 'none', 'Marker', '_', 'MarkerSize', 10, 'Color', color019, 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 1); % Add solid black y-line at 0
title('(Quiet)');
xticks(1); xticklabels({'DNa02'});
ylabel('Firing Rate Difference');
xlim([0 3]); % Extend x-axis padding
ylim(fr_range); % Set y-axis limits for firing rate difference
[~, p_fr_quiet] = ttest(DNa02_diff_sr(:,2));
text(0.5, -15, sprintf('p = %.5f', p_fr_quiet), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', text_size);

% Plot 3: Voltage Difference - All Timepoints
nexttile;
scatter(ones(size(DNa02_diff_vm(:,1))) + jitter_amount * randn(size(DNa02_diff_vm(:,1))), DNa02_diff_vm(:,1), 'MarkerEdgeColor', scatterColor, 'Marker', '.');
hold on;
plot(1, median(DNa02_diff_vm(:,1), 'omitnan'), 'LineStyle', 'none', 'Marker', '_', 'MarkerSize', 10, 'Color', color019, 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 1); % Add solid black y-line at 0
title('(All)');
xticks(1); xticklabels({'DNa02'});
ylabel('Voltage Difference');
xlim([0 3]); % Extend x-axis padding
ylim(vm_range); % Set y-axis limits for voltage difference
[~, p_vm_all] = ttest(DNa02_diff_vm(:,1));
text(0.5, -1.2, sprintf('p = %.5f', p_vm_all), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', text_size);

% Plot 4: Voltage Difference - Quiescent Only
nexttile;
scatter(ones(size(DNa02_diff_vm(:,2))) + jitter_amount * randn(size(DNa02_diff_vm(:,2))), DNa02_diff_vm(:,2), 'MarkerEdgeColor', scatterColor, 'Marker', '.');
hold on;
plot(1, median(DNa02_diff_vm(:,2), 'omitnan'), 'LineStyle', 'none', 'Marker', '_', 'MarkerSize', 10, 'Color', color019, 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 1); % Add solid black y-line at 0
title('(Quiet)');
xticks(1); xticklabels({'DNa02'});
ylabel('Voltage Difference');
xlim([0 3]); % Extend x-axis padding
ylim(vm_range); % Set y-axis limits for voltage difference
[~, p_vm_quiet] = ttest(DNa02_diff_vm(:,2));
text(0.5, -1.2, sprintf('p = %.5f', p_vm_quiet), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', text_size);

% Add main title
sgtitle(['DNa01 +P1 (n = ' num2str(nFly) ')']);
% Save plots as images
cd(plotPath)
saveas(gcf, 'p1_compare_darkness_dna02.png');
set(gcf, 'renderer', 'Painters'); % Save vectorized version
saveas(gcf, 'p1_compare_darkness_dna02.svg');