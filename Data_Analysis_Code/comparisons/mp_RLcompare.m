%% Load RL data
clear; close all

% Define folder containing the .mat files
plotPath = 'E:\\Compare Motion Pulse';
dataPath = 'E:\\Compare Motion Pulse\\data';

% Load data
AOTU019_data = load(fullfile(dataPath, 'AOTU019_Motion_Pulse_RL_All.mat'));
AOTU025_data = load(fullfile(dataPath, 'AOTU025_Motion_Pulse_RL_All.mat'));
behavior_data = load(fullfile(dataPath, 'HeadClosedBehaviorThresh.mat'));

pulse_time = AOTU019_data.time_pulse;
pulse_posR = AOTU019_data.pulse_posR;
AOTU019_meanRL = AOTU019_data.mean_pulse_srRL_station;
AOTU019_semRL = AOTU019_data.sem_pulse_srRL_station;
AOTU025_meanRL_full = AOTU025_data.mean_pulse_srRL_station;
AOTU025_semRL_full = AOTU025_data.sem_pulse_srRL_station;

% Sweep alignment
nSweep_019 = size(AOTU019_meanRL,3);
nSweep_025 = size(AOTU025_meanRL_full,3);
start025 = floor((nSweep_025 - nSweep_019)/2) + 1;
end025 = start025 + nSweep_019 - 1;

AOTU025_meanRL = AOTU025_meanRL_full(:,:,start025:end025);
AOTU025_semRL = AOTU025_semRL_full(:,:,start025:end025);

% Extract behavior data: [time x animals x sweep]
beh_wt = behavior_data.wt_mopRLthresh(:,:,start025:end025); % [time x animals x sweep]
beh_na = behavior_data.na_mopRLthresh(:,:,start025:end025); % [time x animals x sweep]

% Combine animals
beh_all = cat(2, beh_wt, beh_na); % [time x total_animals x sweep]

% Compute mean and SEM across animals for each sweep
turn_meanRL = mean(beh_all, 2, 'omitnan'); % [time x 1 x sweep]
turn_semRL = std(beh_all, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(beh_all), 2)); % [time x 1 x sweep]

% Squeeze for easier plotting
turn_meanRL = squeeze(turn_meanRL); % [time x sweep]
turn_semRL = squeeze(turn_semRL);   % [time x sweep]

% Plot parameters
time_limit = [0 2000];
pos_limit = [-120 120];
fr_limit = [-55 55];
rot_limit = [-75 200];
semAlpha = 0.2;
color_AOTU019 = [0 0 1]; % blue
color_AOTU025 = [0.5 0 0.5]; % purple
color_SUM = [0.2 0.2 0.2]; % dark gray

% Propagate SEM for summed signal
sum_mean = AOTU019_meanRL + AOTU025_meanRL;
sum_sem = sqrt(AOTU019_semRL.^2 + AOTU025_semRL.^2);

% Create figure
figure; set(gcf, 'Position', [100 100 1500 1100])
tiledlayout(4, nSweep_019, 'TileSpacing', 'compact');

for s = 1:nSweep_019
    %% Row 1: Object Position
    nexttile(s)
    plot(pulse_time, pulse_posR(:,2,s), 'k', 'LineWidth', 1.5)
    xlim(time_limit); ylim(pos_limit); yline(0)
    if s == 1
        ylabel('Obj Pos (deg)')
    end

    %% Row 2: AOTU019 + AOTU025
    nexttile(nSweep_019 + s); hold on
    s019 = patch([pulse_time'; flipud(pulse_time')], ...
        [AOTU019_meanRL(:,2,s)-AOTU019_semRL(:,2,s); ...
         flipud(AOTU019_meanRL(:,2,s)+AOTU019_semRL(:,2,s))], ...
         'r', 'FaceAlpha', semAlpha, 'EdgeColor', 'none');
    s019.FaceColor = color_AOTU019;
    plot(pulse_time, AOTU019_meanRL(:,2,s), 'Color', color_AOTU019, 'LineWidth', 1.5)

    s025 = patch([pulse_time'; flipud(pulse_time')], ...
        [AOTU025_meanRL(:,2,s)-AOTU025_semRL(:,2,s); ...
         flipud(AOTU025_meanRL(:,2,s)+AOTU025_semRL(:,2,s))], ...
         'r', 'FaceAlpha', semAlpha, 'EdgeColor', 'none');
    s025.FaceColor = color_AOTU025;
    plot(pulse_time, AOTU025_meanRL(:,2,s), 'Color', color_AOTU025, 'LineWidth', 1.5)

    xlim(time_limit); ylim(fr_limit); yline(0)
    if s == 1
        ylabel('Mean RL (Hz)')
    end

    %% Row 3: Summed RL ± SEM
    nexttile(2*nSweep_019 + s); hold on
    sSUM = patch([pulse_time'; flipud(pulse_time')], ...
        [sum_mean(:,2,s)-sum_sem(:,2,s); ...
         flipud(sum_mean(:,2,s)+sum_sem(:,2,s))], ...
         'r', 'FaceAlpha', semAlpha, 'EdgeColor', 'none');
    sSUM.FaceColor = color_SUM;
    plot(pulse_time, sum_mean(:,2,s), 'Color', color_SUM, 'LineWidth', 1.5)
    xlim(time_limit); ylim(fr_limit); yline(0)
    if s == 1
        ylabel('019 + 025')
    end

    %% Row 4: Behavior
    nexttile(3*nSweep_019 + s); hold on
    bp = patch([pulse_time'; flipud(pulse_time')], ...
        [turn_meanRL(:,s) - turn_semRL(:,s); flipud(turn_meanRL(:,s) + turn_semRL(:,s))], ...
        'r', 'FaceAlpha', semAlpha, 'EdgeColor', 'none');
    bp.FaceColor = color_SUM;
    plot(pulse_time, turn_meanRL(:,s), 'Color', color_SUM, 'LineWidth', 1.5)
    xlim(time_limit); ylim(rot_limit); yline(0)
    if s == 1
        ylabel('Turn RL (Hz)')
    end
end

% Save the figure
cd(plotPath)
saveas(gcf, 'RL_comparisonovertime.png');
set(gcf, 'renderer', 'Painters')
saveas(gcf, 'RL_comparisonovertime.svg');


%% Compare raw peak firing rates across sweep positions (positive only)

% Load behavior data
behavior = load(fullfile(dataPath, 'HeadClosedBehaviorThresh.mat'));
peakWT = behavior.peakWT(3:11,:);   % Middle 9 sweeps
peakGFP = behavior.peakNA(3:11,:); % Middle 9 sweeps

% Define middle 9 sweep positions
sweepPos = [-78.75, -56.25, -33.75, -11.25, 11.25, 33.75, 56.25, 78.75, 101.25];
posIdx = sweepPos > 0; % Index for positive sweep positions
sweepPos_pos = sweepPos(posIdx);

% Combine and clean behavior data
all_behavior = [peakWT; peakGFP]; % 9 x n
validIdx = all(~isnan(all_behavior),1);
behavior_use = all_behavior(:,validIdx); % 9 x valid_animals
behavior_use = behavior_use ./ max(behavior_use, [], 1); % normalize per animal

% Compute behavior mean ± SEM
mean_behavior = mean(behavior_use,2);
sem_behavior = std(behavior_use,0,2) ./ sqrt(size(behavior_use,2));

% Load AOTU019 and AOTU025 raw peak spike rates (already 9 values)
AOTU019_peak_raw = AOTU019_data.mean_peak_srRL_station(1,:)';
AOTU019_sem_raw  = AOTU019_data.sem_peak_srRL_station(1,:)';
AOTU025_peak_raw = AOTU025_data.mean_peak_srRL_station(1,3:11)';
AOTU025_sem_raw  = AOTU025_data.sem_peak_srRL_station(1,3:11)';

% Normalize individual traces for tile 1
AOTU019_peak_norm = AOTU019_peak_raw ./ max(AOTU019_peak_raw);
AOTU025_peak_norm = AOTU025_peak_raw ./ max(AOTU025_peak_raw);
AOTU019_sem_norm  = AOTU019_sem_raw ./ max(AOTU019_peak_raw);
AOTU025_sem_norm  = AOTU025_sem_raw ./ max(AOTU025_peak_raw);

% Sum raw traces first, then normalize result for tile 2
sum_raw = AOTU019_peak_raw + AOTU025_peak_raw;
sum_sem = sqrt(AOTU019_sem_raw.^2 + AOTU025_sem_raw.^2);
sum_peak_norm = sum_raw ./ max(sum_raw);
sum_sem_norm  = sum_sem ./ max(sum_raw);

% Filter all signals by positive sweep positions
AOTU019_peak_pos = AOTU019_peak_norm(posIdx);
AOTU019_sem_pos  = AOTU019_sem_norm(posIdx);
AOTU025_peak_pos = AOTU025_peak_norm(posIdx);
AOTU025_sem_pos  = AOTU025_sem_norm(posIdx);
sum_peak_pos     = sum_peak_norm(posIdx);
sum_sem_pos      = sum_sem_norm(posIdx);
mean_behavior_pos = mean_behavior(posIdx);
sem_behavior_pos  = sem_behavior(posIdx);

% Plot normalized raw peak responses
figure;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

% Tile 1: Normalized AOTU019 and AOTU025 vs Behavior
nexttile; hold on
errorbar(sweepPos_pos, AOTU019_peak_pos, AOTU019_sem_pos, ...
    'Color', [0 0 1], 'LineWidth', 1.5, 'CapSize', 0) % blue
errorbar(sweepPos_pos, AOTU025_peak_pos, AOTU025_sem_pos, ...
    'Color', [0.5 0 0.5], 'LineWidth', 1.5, 'CapSize', 0) % purple
errorbar(sweepPos_pos, mean_behavior_pos, sem_behavior_pos, ...
    'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'CapSize', 0) % dark gray
ylim([0 1.2])
xlabel('Sweep Position (deg)')
ylabel('Norm. Peak Resp.')
title('AOTU019 vs AOTU025 vs Behavior')
legend({'AOTU019','AOTU025','Behavior'}, 'Location', 'northwest')

% Tile 2: Summed raw → normalized vs Behavior
nexttile; hold on
errorbar(sweepPos_pos, sum_peak_pos, sum_sem_pos, ...
    'Color', [1 0 0], 'LineWidth', 1.5, 'CapSize', 0) % red
errorbar(sweepPos_pos, mean_behavior_pos, sem_behavior_pos, ...
    'Color', 'k', 'LineWidth', 1.5, 'CapSize', 0) % red dashed
ylim([0 1.2])
xlabel('Sweep Position (deg)')
ylabel('Norm. Peak Resp.')
title('Summed Visual vs Behavior')
legend({'019+025','Behavior'}, 'Location', 'northwest')

% Save the figure
cd(plotPath)
saveas(gcf, 'RL_comparisonpeak.png');
set(gcf, 'renderer', 'Painters')
saveas(gcf, 'RL_comparisonpeak.svg');