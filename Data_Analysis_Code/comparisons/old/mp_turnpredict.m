% mp_fr_to_turn_prediction
% Script to load AOTU019 and AOTU025 firing rate (fr) and turn data one at a time,
% find overlapping position bins, and generate a plot comparing normalized 
% firing rate and turn data with error bars.
%
% The script loads individual .mat files for AOTU019 and AOTU025, extracts 
% firing rate and turn data, and computes combined and normalized responses.
% It generates three plots: 
%   1. Comparing AOTU019 and AOTU025 firing rates with turning rate
%   2. Linear and nonlinear combined firing rate differences
%   3. Headclosed turn data, showing peakWT turn with linear and nonlinear combined responses
%
% The final plots are saved as PNG and SVG.
% CREATED: 11/10/2024 - MC

clear
close all

%% Load AOTU019 and AOTU025 Firing Rate and Turning Data
plotPath = 'E:\Compare Motion Pulse';  % Folder containing plots
dataPath = 'E:\Compare Motion Pulse\data';  % Folder containing the .mat files
cd(dataPath);

% Define file paths for AOTU019 and AOTU025 data
AOTU019_turn_file = fullfile(dataPath, 'AOTU019_Motion_Pulse_onlyturn_25_dps_Pursuit.mat');
AOTU019_fr_file = fullfile(dataPath, 'AOTU019_Motion_Pulse_fr_25_dps_Pursuit.mat');
AOTU025_turn_file = fullfile(dataPath, 'AOTU025_Motion_Pulse_onlyturn_25_dps_Pursuit.mat');
AOTU025_fr_file = fullfile(dataPath, 'AOTU025_Motion_Pulse_fr_25_dps_Pursuit.mat');

% Load data files
AOTU019_turn = load(AOTU019_turn_file);
AOTU019_fr = load(AOTU019_fr_file);
AOTU025_turn = load(AOTU025_turn_file);
AOTU025_fr = load(AOTU025_fr_file);

%% Extract and Normalize Firing Rates and Turning Data
% Extract position bins, firing rates, and SEM for both AOTU019 and AOTU025
positionBins_AOTU019 = AOTU019_fr.combinedData(:, 1);
normalizedFR_AOTU019 = AOTU019_fr.combinedData(:, 2);
sem_AOTU019_FR = AOTU019_fr.combinedData(:, 3);

positionBins_AOTU025 = AOTU025_fr.combinedData(:, 1);
normalizedFR_AOTU025 = AOTU025_fr.combinedData(:, 2);
sem_AOTU025_FR = AOTU025_fr.combinedData(:, 3);

% Turning rates and SEM for AOTU019 and AOTU025
normalizedTurn_AOTU019 = AOTU019_turn.combinedData(:, 2);
sem_AOTU019_Turn = AOTU019_turn.combinedData(:, 3);

normalizedTurn_AOTU025 = AOTU025_turn.combinedData(:, 2);
sem_AOTU025_Turn = AOTU025_turn.combinedData(:, 3);

% Find overlapping position bins
[commonBins, idxAOTU019, idxAOTU025] = intersect(positionBins_AOTU019, positionBins_AOTU025);

% Filter to include only positive object positions
positiveIdx = commonBins > 0;
commonBins = commonBins(positiveIdx);
overlapFR_AOTU019 = normalizedFR_AOTU019(idxAOTU019(positiveIdx));
overlapFR_AOTU025 = normalizedFR_AOTU025(idxAOTU025(positiveIdx));
sem_FR_AOTU019 = sem_AOTU019_FR(idxAOTU019(positiveIdx));
sem_FR_AOTU025 = sem_AOTU025_FR(idxAOTU025(positiveIdx));
overlapTurn_AOTU019 = normalizedTurn_AOTU019(idxAOTU019(positiveIdx));
overlapTurn_AOTU025 = normalizedTurn_AOTU025(idxAOTU025(positiveIdx));
sem_Turn_AOTU019 = sem_AOTU019_Turn(idxAOTU019(positiveIdx));
sem_Turn_AOTU025 = sem_AOTU025_Turn(idxAOTU025(positiveIdx));

% Normalize overlap firing rates to range [0, 1]
overlapFR_AOTU019 = (overlapFR_AOTU019 - min(overlapFR_AOTU019)) / (max(overlapFR_AOTU019) - min(overlapFR_AOTU019));
overlapFR_AOTU025 = (overlapFR_AOTU025 - min(overlapFR_AOTU025)) / (max(overlapFR_AOTU025) - min(overlapFR_AOTU025));

%% Calculate Combined Firing Rate and Turning Response with SEM
% Combined normalized turning rate
turn_headopen = (overlapTurn_AOTU019 + overlapTurn_AOTU025) / 2;
% Pooled SEM for combined rates
turn_headopenSEM = sqrt((sem_Turn_AOTU019.^2 + sem_Turn_AOTU025.^2) / 2);

%% Calculate Linear and Nonlinear Combined Differences
% Linear contra - ipsi firing rate difference, weighted
combinedContraFR = -overlapFR_AOTU019;    % AOTU019 as contralateral (inhibitory)
combinedIpsiFR = overlapFR_AOTU025; % AOTU025 as ipsilateral (excitatory)
contraMinusIpsiFR = combinedIpsiFR - combinedContraFR;
contraMinusIpsiFR_unweighted = contraMinusIpsiFR / max(abs(contraMinusIpsiFR));

% Linear contra - ipsi firing rate difference, weighted
combinedContraFR = -overlapFR_AOTU019;    % AOTU019 as contralateral (inhibitory)
combinedIpsiFR = overlapFR_AOTU025; % AOTU025 as ipsilateral (excitatory)
contraMinusIpsiFR = combinedIpsiFR - combinedContraFR;
contraMinusIpsiFR_weighted = contraMinusIpsiFR / max(abs(contraMinusIpsiFR));

% SEM for the combined contra - ipsi difference, normalized
combinedSEM_Difference = sqrt((sem_FR_AOTU019.^2 + sem_FR_AOTU025.^2) / 2);
combinedSEM_Difference_normalized = combinedSEM_Difference / max(abs(contraMinusIpsiFR));

% Define DN input range and output transformation using adjELU
DNa02input = linspace(-1, 1, 1000);  % Input range
alpha = 0.9;
DNa02output = linspace(-1, 1, 1000);  % Output range
%DNa02output = adjELU(DNa02input, alpha, -0.1);  % Apply ELU with alpha = 0.5 and shift = 0

% Apply nonlinearity to each data point in combinedContraFR and combinedIpsiFR
combinedContraFR = -overlapFR_AOTU019;    % AOTU019 as contralateral (inhibitory)
combinedIpsiFR = overlapFR_AOTU025; % AOTU025 as ipsilateral (excitatory)
contraNonlinearFR = arrayfun(@(x) interp1(DNa02input, DNa02output, x, 'nearest', 'extrap'), combinedContraFR);
ipsiNonlinearFR = arrayfun(@(x) interp1(DNa02input, DNa02output, x, 'nearest', 'extrap'), combinedIpsiFR);
% Nonlinear contra - ipsi difference, normalized
nonlinearContraMinusIpsiFR = ipsiNonlinearFR - contraNonlinearFR;
nonlinearContraMinusIpsiFR_weighted = nonlinearContraMinusIpsiFR / max(abs(nonlinearContraMinusIpsiFR));

% Unweighted
combinedContraFR = -overlapFR_AOTU019;    % AOTU019 as contralateral (inhibitory)
combinedIpsiFR = overlapFR_AOTU025; % AOTU025 as ipsilateral (excitatory)
contraNonlinearFR = arrayfun(@(x) interp1(DNa02input, DNa02output, x, 'nearest', 'extrap'), combinedContraFR);
ipsiNonlinearFR = arrayfun(@(x) interp1(DNa02input, DNa02output, x, 'nearest', 'extrap'), combinedIpsiFR);
% Nonlinear contra - ipsi difference, normalized
nonlinearContraMinusIpsiFR = ipsiNonlinearFR - contraNonlinearFR;
nonlinearContraMinusIpsiFR_unweighted = nonlinearContraMinusIpsiFR / max(abs(nonlinearContraMinusIpsiFR));

figure; set(gcf, 'Position', [100 100 500 500]); 
plot(DNa02input,DNa02output)
axis tight
% Save the combined plot
cd(plotPath)
sgtitle('DNa02 Nonlinearity');
saveas(gcf, 'DNa02nonlinearity.png');
saveas(gcf, 'DNa02nonlinearity.svg');

% Load Headclosed Behavior Data and Calculate Mean and SEM
cd(dataPath)
load('HeadClosedBehaviorThresh.mat', 'peakWT');
peakWT_selected = peakWT(7:11, :);  % Select rows
turn_headclosed = mean(peakWT_selected, 2, 'omitnan'); % Mean across animals
turnSEM_headclosed = std(peakWT_selected, 0, 2, 'omitnan') / sqrt(size(peakWT_selected, 2)); % SEM across animals

% Normalize mean and SEM to max mean of 1
maxMeanTurnWT = max(turn_headclosed);
turn_headclosed = turn_headclosed / maxMeanTurnWT;
turnSEM_headclosed = turnSEM_headclosed / maxMeanTurnWT;

% Calculate Combined Turn Response (Head Open + Head Closed) with SEM
combinedTurn_with_HeadClosed = turn_headopen;  % Combine head open and closed turn rates
combinedSEM_with_HeadClosed = turn_headopenSEM;  % Pooled SEM

% Plot Linear and Nonlinear Combined Difference (Contra - Ipsi) Normalized to Peak Turning Rate

figure; set(gcf, 'Position', [100 100 300 900]); 
tiledlayout(3, 1, "TileSpacing", "compact");

% First tile: Plot the peak normalized firing rate data for AOTU019 and AOTU025 with error bars
nexttile;
errorbar(commonBins, overlapFR_AOTU019, sem_FR_AOTU019, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0, 'Color', "#77AC30", 'LineWidth', 1, 'DisplayName', 'AOTU019 FR');
hold on;
errorbar(commonBins, overlapFR_AOTU025, sem_FR_AOTU025, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0, 'Color', "#7E2F8E", 'LineWidth', 1, 'DisplayName', 'AOTU025 FR');
xlabel('Position Bins');
lgd = legend('show', 'Location', 'southeast'); lgd.FontSize = 5;
xticks(0:30:150); % Set x-axis ticks every 30 degrees
xlim([0 120])
ylim([-0.2 1.3]);
grid off;

% Second tile: Linear and nonlinear combined firing rate differences
nexttile;
errorbar(commonBins, combinedTurn_with_HeadClosed, combinedSEM_with_HeadClosed, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0, 'Color', '#4DBEEE', 'LineWidth', 1, 'DisplayName', 'Headopen + Headclosed Turn');
hold on;
errorbar(commonBins, contraMinusIpsiFR_unweighted, combinedSEM_Difference_normalized, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName', 'Linear FRuw');
xlabel('Position Bins');
lgd = legend('show', 'Location', 'southeast'); lgd.FontSize = 5;
xticks(0:30:150); % Set x-axis ticks every 30 degrees
xlim([0 120])
ylim([0 1.3]);
grid off;

% Third tile: Headclosed turn data (peakWT) with linear and nonlinear responses
nexttile;
errorbar(commonBins, combinedTurn_with_HeadClosed, combinedSEM_with_HeadClosed, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0, 'Color', '#4DBEEE', 'LineWidth', 1, 'DisplayName', 'Headopen + Headclosed Turn');
hold on;
errorbar(commonBins, nonlinearContraMinusIpsiFR_unweighted, combinedSEM_Difference_normalized, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName', 'Nonlinear FR');
xlabel('Position Bins');
lgd = legend('show', 'Location', 'southeast'); lgd.FontSize = 5;
xticks(0:30:150); % Set x-axis ticks every 30 degrees
xlim([0 120])
ylim([0 1.3]);
grid off;

% Save the combined plot
cd(plotPath)
sgtitle('FR and Turn Response Comparison');
saveas(gcf, 'mp_fr_to_turn_prediction.png');
saveas(gcf, 'mp_fr_to_turn_prediction.svg');
