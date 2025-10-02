% SCRIPT: Compare AOTU019 and AOTU025 Firing Rates Across Sweep Positions
% CREATED: 12/21/2024 - MC
%
% This script loads peak firing rate data for AOTU019 and AOTU025 neurons,
% identifies common sweep positions presented to both, and performs the
% following analyses:
%
%
% INPUTS:
% - AOTU019_Motion_Pulse_sweeppeaks_25_dps_All.mat: Data for AOTU019 neurons.
% - AOTU025_Motion_Pulse_sweeppeaks_25_dps_All.mat: Data for AOTU025 neurons.
%
clear
close all

%% Load data
plotPath = 'E:\Compare Motion Pulse';  % Define folder for saving plots and results
dataPath = 'E:\Compare Motion Pulse\data';  % Define folder containing the .mat files
cd(dataPath)

% Define file paths
AOTU019_file = fullfile(dataPath, 'AOTU019_aIPg_Motion_Pulse_sweeppeaks_25_dps_All.mat');
AOTU025_file = fullfile(dataPath, 'AOTU025_aIPg_Motion_Pulse_sweeppeaks_25_dps_All.mat');

mean_colors = {"#0072BD", "#7E2F8E"};  % Colors for AOTU019 and AOTU025

%% Load datasets
AOTU019_data = load(AOTU019_file);  % Load AOTU019 dataset
AOTU025_data = load(AOTU025_file);  % Load AOTU025 dataset

% Extract relevant data
avg_srR_AOTU019 = AOTU019_data.avg_srRightward;  % Peak firing rates (rows: animals, columns: positions)
sweepPosR_AOTU019 = AOTU019_data.sweepPosR;  % Sweep positions

avg_srR_AOTU025 = AOTU025_data.avg_srRightward;  % Peak firing rates (rows: animals, columns: positions)
sweepPosR_AOTU025 = AOTU025_data.sweepPosR;  % Sweep positions

%% Restrict to overlapping sweep positions
common_sweepPos = intersect(sweepPosR_AOTU019, sweepPosR_AOTU025);  % Find common positions
common_sweepPos = common_sweepPos(4:end);

% Find indices of common positions
[~, idx_AOTU019] = ismember(common_sweepPos, sweepPosR_AOTU019);
[~, idx_AOTU025] = ismember(common_sweepPos, sweepPosR_AOTU025);

% Restrict data to common positions
restricted_srR_AOTU019 = avg_srR_AOTU019(:, idx_AOTU019);  % Columns are common positions
restricted_srR_AOTU025 = avg_srR_AOTU025(:, idx_AOTU025);

%% Calculate mean and SEM
mean_AOTU019 = mean(restricted_srR_AOTU019, 1, 'omitnan');  % Mean across animals for each position
sem_AOTU019 = std(restricted_srR_AOTU019, 0, 1, 'omitnan') ./ sqrt(size(restricted_srR_AOTU019, 1));

mean_AOTU025 = mean(restricted_srR_AOTU025, 1, 'omitnan');  % Mean across animals for each position
sem_AOTU025 = std(restricted_srR_AOTU025, 0, 1, 'omitnan') ./ sqrt(size(restricted_srR_AOTU025, 1));

%% Plot mean and SEM
figure;
hold on;
errorbar(common_sweepPos, mean_AOTU019, sem_AOTU019, 'Color', mean_colors{1}, 'LineWidth', 1.5, 'DisplayName', 'AOTU019');
errorbar(common_sweepPos, mean_AOTU025, sem_AOTU025, 'Color', mean_colors{2}, 'LineWidth', 1.5, 'DisplayName', 'AOTU025');

% Add plot details
xlabel('Sweep Position');
ylabel('Peak Firing Rate');
title('Comparison of AOTU019 and AOTU025 Peak Firing Rates');
legend('Location', 'Best');
grid on;

%% Prepare table for Linear Mixed-Effects Model

% Combine data
data_combined = [restricted_srR_AOTU019(:); restricted_srR_AOTU025(:)];

% Factors
cell_type_label = [repmat("AOTU019", numel(restricted_srR_AOTU019), 1); ...
                   repmat("AOTU025", numel(restricted_srR_AOTU025), 1)];

sweep_pos_label = [repmat(common_sweepPos, size(restricted_srR_AOTU019, 1), 1); ...
                   repmat(common_sweepPos, size(restricted_srR_AOTU025, 1), 1)];

fly_ids = [repelem("Fly019_" + (1:size(restricted_srR_AOTU019,1))', numel(common_sweepPos)); ...
           repelem("Fly025_" + (1:size(restricted_srR_AOTU025,1))', numel(common_sweepPos))];

% Create data table
T = table(data_combined, cell_type_label, sweep_pos_label, fly_ids, ...
    'VariableNames', {'FiringRate', 'CellType', 'SweepPos', 'FlyID'});

% Convert to categorical
T.CellType = categorical(T.CellType);
T.SweepPos = categorical(T.SweepPos);
T.FlyID = categorical(T.FlyID);

%% Run linear mixed-effects model with interaction

lme = fitlme(T, 'FiringRate ~ CellType*SweepPos + (1|FlyID)');

% Display results
disp(anova(lme,'DFMethod','satterthwaite'));
