% CREATED: 10/05/2025 - MC
%
%% Load direction selectivity data
clear; close all

plotPath = 'E:\Compare Motion Pulse';  % Define folder containing the .mat files
dataPath = 'E:\Compare Motion Pulse\data';  % Define folder containing the .mat files
cd(dataPath)

AOTU019_file = fullfile(dataPath, 'snr_AOTU019_aIPg_Motion_Pulse.mat');
AOTU025_file = fullfile(dataPath, 'snr_AOTU025_aIPg_Motion_Pulse.mat');

mean_colors = {"#0072BD","#7E2F8E"};

%% Load in dataset

% average between best sweep positions for each cell type
AOTU019 = load(AOTU019_file);
AOTU019_var = mean(AOTU019.storeVar(:,5:6),2);
AOTU019_avg = mean(AOTU019.storeMean(:,5:6),2);
AOTU019_snr = mean(AOTU019.storeSNR(:,5:6),2);

AOTU025 = load(AOTU025_file);
AOTU025_var = mean(AOTU025.storeVar(:,6:7),2);
AOTU025_avg = mean(AOTU025.storeMean(:,6:7),2);
AOTU025_snr = mean(AOTU025.storeSNR(:,6:7),2);

% Define number of flies for each cell type
n019 = sum(~isnan(AOTU019_snr));
n025 = sum(~isnan(AOTU025_snr));

%% Plot SNR comparison

% Data (vectors; rows = animals)
snr019 = AOTU019_snr(:);
snr025 = AOTU025_snr(:);

% Drop NaNs
snr019 = snr019(~isnan(snr019));
snr25  = snr025(~isnan(snr025));

% Colors and jitter
scatterColor = 'k';
jitter_amount = 0.015;
dashHalf = 0.15;                % half-length of median dash

% Figure
figure; set(gcf,'Position',[100 100 500 500]);
ax = axes; hold(ax,'on');

% Scatter (with jitter)
n1 = numel(snr019); n2 = numel(snr25);
scatter(1 + jitter_amount*randn(n1,1), snr019, '.', 'MarkerEdgeColor', scatterColor);
scatter(2 + jitter_amount*randn(n2,1), snr25,  '.', 'MarkerEdgeColor', scatterColor);

% Medians (horizontal dashes)
med019 = median(snr019,'omitnan');
med025 = median(snr25 ,'omitnan');
plot([1-dashHalf 1+dashHalf],[med019 med019],'Color',mean_colors{1},'LineWidth',2);
plot([2-dashHalf 2+dashHalf],[med025 med025],'Color',mean_colors{2},'LineWidth',2);

% Axes & labels
xlim([0.5 2.5]); 
xticks([1 2]); 
xticklabels({'AOTU019','AOTU025'});
ylabel('SNR (mean/var of per-trial means)');
grid on; box off;

% Two-sample t-test (AOTU019 vs AOTU025)
[~, pval] = ttest2(snr019, snr25);
% Add text annotation with p-value (top-right corner)
pstr = sprintf('p = %.5g', pval);
annotation('textbox', [0.65 0.85 0.3 0.1], 'String', pstr, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none', 'FontSize', 10);


% Save plots
cd(plotPath)
saveas(gcf, 'mp_snr_aipg.png');
set(gcf, 'renderer', 'Painters');
saveas(gcf, 'mp_snr_aipg.svg');
