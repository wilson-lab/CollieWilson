% CREATED: 10/05/2025 - MC
%
%% Load direction selectivity data
clear; close all

plotPath = 'E:\Compare Motion Pulse';  % Define folder containing the .mat files
dataPath = 'E:\Compare Motion Pulse\data';  % Define folder containing the .mat files
cd(dataPath)

AOTU019_file = fullfile(dataPath, 'snr_AOTU019_Motion_Pulse.mat');
AOTU025_file = fullfile(dataPath, 'snr_AOTU025_Motion_Pulse.mat');
DNa02_file = fullfile(dataPath,'snr_DNa02_Motion_Pulse.mat');

target_speeds = {'25 dps','75 dps'};
mean_colors = {"#0072BD","#7E2F8E"};

%% Load in dataset

% Load files
AOTU019 = load(AOTU019_file);
AOTU025 = load(AOTU025_file);

% Compute per-fly SNR averaged across best sweep bins
% AOTU019 uses sweep bins 5:6; AOTU025 uses 10:11

% AOTU019
snr019_slow = mean(AOTU019.storeSNR(:,5:6,2), 2);  % z=2
snr019_fast = mean(AOTU019.storeSNR(:,5:6,1), 2);  % z=1
pair019 = [snr019_slow snr019_fast];
pair019 = pair019(all(~isnan(pair019),2),:);

% AOTU025
snr025_slow = mean(AOTU025.storeSNR(:,10:11,2), 2); % z=2
snr025_fast = mean(AOTU025.storeSNR(:,10:11,1), 2); % z=1
pair025 = [snr025_slow snr025_fast];
pair025 = pair025(all(~isnan(pair025),2),:);

% Aesthetics
flyColor  = [0.65 0.65 0.65];  % lines/points per fly
dashHalf  = 0.15;              % half-length of median dash in x units

%% Plot (two tiles, one per cell type)
figure; set(gcf,'Position',[100 100 500 700]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for panel = 1:2
    nexttile; hold on;

    % Choose data by cell type
    switch panel
        case 1
            pairData = pair019;
            cellTitle = 'AOTU019';
        case 2
            pairData = pair025;
            cellTitle = 'AOTU025';
    end

    % Connect each fly (line) and plot dot at each speed (no jitter)
    for i = 1:size(pairData,1)
        plot([1 2], pairData(i,:), '-', 'Color', flyColor, 'LineWidth', 1);
        plot([1 2], pairData(i,:), '.', 'Color', flyColor, 'MarkerSize', 8);
    end

    % Medians (horizontal dashes) for each speed
    dashColor = mean_colors{panel};
    meds = median(pairData, 1, 'omitnan'); % [slow fast]
    plot([1-dashHalf 1+dashHalf], [meds(1) meds(1)], '-', 'Color', dashColor, 'LineWidth', 3);
    plot([2-dashHalf 2+dashHalf], [meds(2) meds(2)], '-', 'Color', dashColor, 'LineWidth', 3);

    xlim([0.5 2.5]); ylim([0 1.5]);
    xticks([1 2]); xticklabels(target_speeds);
    ylabel('SNR (mean/var of per-trial means)');
    title(cellTitle);
    grid on; box off;
end

sgtitle('SNR by Target Speed (paired within fly)');

%% Repeated-measures LME: effects of CellType and Target Speed
% (assumes pair019 and pair025 already exist as [slow fast] per-fly matrices)

% Build long-format table
n019 = size(pair019,1);
n025 = size(pair025,1);

SNR      = [pair019(:);                pair025(:)];
Speed    = [repmat({'slow';'fast'}, n019, 1); repmat({'slow';'fast'}, n025, 1)];
CellType = [repmat({'AOTU019'}, 2*n019, 1);  repmat({'AOTU025'}, 2*n025, 1)];
Animal   = [strcat("019_", string(repelem((1:n019)',2))); strcat("025_", string(repelem((1:n025)',2)))];

tbl = table(SNR, categorical(CellType), categorical(Speed), categorical(Animal), ...
            'VariableNames', {'SNR','CellType','Speed','Animal'});

% Fit linear mixed-effects model with random intercept per animal
lme = fitlme(tbl, 'SNR ~ CellType*Speed + (1|Animal)', 'FitMethod','REML');

% ANOVA table with Satterthwaite DF
stats = anova(lme, 'DFMethod', 'Satterthwaite');

% Extract p-values
p_cell  = stats.pValue(strcmp(stats.Term,'CellType'));
p_speed = stats.pValue(strcmp(stats.Term,'Speed'));
p_int   = stats.pValue(strcmp(stats.Term,'CellType:Speed'));

% Add p-values to figure (top-right) BEFORE saving
pstr = sprintf(['Cell Type: p = %.5g\n' ...
                'Speed: p = %.5g\n' ...
                'x: p = %.5g'], p_cell, p_speed, p_int);

annotation('textbox',[0.60 0.80 0.38 0.18], 'String', pstr, ...
           'HorizontalAlignment','right','VerticalAlignment','top', ...
           'EdgeColor','none', 'Interpreter','none', 'FontSize',10);

% Save
cd(plotPath)
saveas(gcf, 'mp_snr.png');
set(gcf, 'renderer', 'Painters');
saveas(gcf, 'mp_snr.svg');

