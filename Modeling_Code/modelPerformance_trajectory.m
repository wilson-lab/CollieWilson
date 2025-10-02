% modelPerformance_trajectory
%
% This function runs and evaluates a fly steering model across different
% noise levels to assess how internal variability impacts trajectory control.
% It compares model performance across specified conditions (e.g., synaptic 
% strength or direction selectivity), calculates behavioral metrics, and 
% visualizes example trajectories at selected noise levels.
%
% INPUTS:
%   predicted_RF      - Struct containing the predicted receptive fields for model neurons.
%   comparisonType    - Comparison mode: 'synapse', 'strength', or 'dirselective'.
%
% CREATED: 04/06/2025 - MC
% UPDATED: 07/31/2025 - MC added repeated runs for DSI values

function modelPerformance_trajectory(predicted_RF, comparisonType)
%% Initialize

[folder, plotSettings, runSettings] = modelSettings();
close all;

% Define range of noise levels to test
noiseValues = 0:0.1:1.2;
nNoise = length(noiseValues);
kVal = 1; % fixed steering gain

stability_startPos = 0;
stability_simDuration = 60;
conditionColors = {"#0072BD"; "#D95319";"#7E2F8E"};

% Variables depending on comparison type
thisSynapse = "inhibitory";
if strcmp(comparisonType, 'strength')
    comparisonLabel = {'Strength 1', 'Strength 0'};
    nComp = 2;
elseif strcmp(comparisonType, 'dirselective')
    comparisonLabel = {'neither', '019 DS','025 DS'};
    nComp = 3;
end

% Preallocate arrays
metrics_prob = zeros(nNoise, nComp);
metrics_var = zeros(nNoise, nComp);
metrics_ISE = zeros(nNoise, nComp);
metrics_IAE = zeros(nNoise, nComp);
evt_avg = [];

% Select up to X evenly spaced noise values for plotting
maxPlotN = min(6, nNoise);
selected_noise_indices = round(linspace(1, nNoise, maxPlotN));

% Initialize figure for example plots
figure; set(gcf, 'Position', [100 100 1500 900]);
tiledlayout(maxPlotN, nComp, 'TileSpacing', 'compact');

% Loop over noise levels
for nIdx = 1:nNoise
    thisNoise = noiseValues(nIdx);
    runSettings.nNoise = nIdx;
    disp(['Running simulations for noise = ', num2str(thisNoise)]);

    for idx = 1:nComp
        % load fresh copy of RF
        thisTuning = predicted_RF;
        if strcmp(comparisonType, 'strength')
            switch idx
                case 2
                    thisTuning.AOTU019 = thisTuning.AOTU019 .* 0;
            end
        elseif strcmp(comparisonType, 'dirselective')
            runSettings.dirselective = idx - 1;
        end

        % Run model (stability only)
        [timebase, visobj_history, ~, rotvel_history] = aotu_steering_model(thisTuning, thisNoise, stability_startPos, thisSynapse, stability_simDuration, kVal, 2, runSettings);

        % Performance metrics
        metrics_results = calculatePerformanceMetrics(visobj_history, rotvel_history, timebase, plotSettings);
        metrics_prob(nIdx, idx) = metrics_results.prob;
        metrics_var(nIdx, idx) = metrics_results.var;
        metrics_ISE(nIdx, idx) = metrics_results.ISE;
        metrics_IAE(nIdx, idx) = metrics_results.IAE;

        % Analyze this run set
        [posvang, posBins] = analyzeErrorVsTurn(visobj_history, rotvel_history, runSettings);
        % Store results
        evt_avg(nIdx, :, idx) = posvang;

        % Plot example run for selected noise levels
        if any(nIdx == selected_noise_indices)
            visobj_plot = remove_large_jumps(visobj_history(1, :), 180);

            nexttile
            plot(timebase, visobj_plot, 'Color', conditionColors{idx});
            xlabel('Time (s)');
            ylabel('Position (deg)');
            title(['Noise = ' num2str(thisNoise)]);
            yline(0);
            axis tight
            ylim([-50 50]);

            if nIdx == 1
                legend(comparisonLabel{idx});
            end
        end
    end
end

% Save the example plot figure
sgtitle(['Example Runs for ' comparisonType ' across noise levels']);
saveas(gcf, fullfile(folder.final, [comparisonType '_NoiseExampleRuns.png']));
set(gcf,'renderer','Painters');
saveas(gcf, fullfile(folder.vectors, [comparisonType '_NoiseExampleRuns.svg']));

% Summary Plots
generateSummaryMetrics(noiseValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, comparisonLabel, comparisonType, folder);
generateSummaryEVT(nNoise, noiseValues, evt_avg, posBins, comparisonLabel, comparisonType, folder);


%% Run AGAIN for different DSI values

[folder, plotSettings, runSettings] = modelSettings();
close all;

% Define range of penalty levels to test
penaltyValues = 0:0.1:1;
dsiValues = (1 - penaltyValues) ./ (1 + penaltyValues); % convert to DSI
nPV = length(penaltyValues);

thisNoise = 0.7;

stability_startPos = 0;
stability_simDuration = 60;
conditionColors = {"#0072BD"; "#D95319";"#7E2F8E"};

% Preallocate arrays
metrics_prob = zeros(nPV, nComp);
metrics_var = zeros(nPV, nComp);
metrics_ISE = zeros(nPV, nComp);
metrics_IAE = zeros(nPV, nComp);
evt_avg = [];

% Select up to X evenly spaced noise values for plotting
maxPlotN = min(6, nPV);
selected_noise_indices = round(linspace(1, nPV, maxPlotN));

% Initialize figure for example plots
figure; set(gcf, 'Position', [100 100 1500 900]);
tiledlayout(maxPlotN, nComp, 'TileSpacing', 'compact');

% Loop over penalty levels
for nIdx = 1:nPV
    thisPenalty = penaltyValues(nIdx);
    thisDSI = dsiValues(nIdx);
    runSettings.AOTU019dsi_penalty = thisPenalty;
    disp(['Running simulations for dsi = ', num2str(thisDSI)]);

    for idx = 1:nComp
        % load fresh copy of RF
        thisTuning = predicted_RF;
        % select condition dependent features
        if strcmp(comparisonType, 'strength')
            switch idx
                case 2
                    thisTuning.AOTU019 = thisTuning.AOTU019 .* 0;
            end
        elseif strcmp(comparisonType, 'dirselective')
            % set no (0), 019 DS (1), or 025 DS (2) call
            runSettings.dirselective = idx - 1;    
        end

        % Run model (stability only)
        [timebase, visobj_history, ~, rotvel_history] = aotu_steering_model(thisTuning, thisNoise, stability_startPos, thisSynapse, stability_simDuration, kVal, 2, runSettings);

        % Performance metrics
        metrics_results = calculatePerformanceMetrics(visobj_history, rotvel_history, timebase, plotSettings);
        metrics_prob(nIdx, idx) = metrics_results.prob;
        metrics_var(nIdx, idx) = metrics_results.var;
        metrics_ISE(nIdx, idx) = metrics_results.ISE;
        metrics_IAE(nIdx, idx) = metrics_results.IAE;

        % Analyze this run set
        [posvang, posBins] = analyzeErrorVsTurn(visobj_history, rotvel_history, runSettings);
        % Store results
        evt_avg(nIdx, :, idx) = posvang;

        % Plot example run for selected noise levels
        if any(nIdx == selected_noise_indices)
            visobj_plot = remove_large_jumps(visobj_history(1, :), 180);

            nexttile
            plot(timebase, visobj_plot, 'Color', conditionColors{idx});
            xlabel('Time (s)');
            ylabel('Position (deg)');
            title(['DSI = ' num2str(thisDSI)]);
            yline(0);
            axis tight
            ylim([-50 50]);

            if nIdx == 1
                legend(comparisonLabel{idx});
            end
        end
    end
end

% Save the example plot figure
sgtitle(['Example Runs for ' comparisonType ' across dsi levels']);
saveas(gcf, fullfile(folder.final, [comparisonType '_DSIExampleRuns.png']));
set(gcf,'renderer','Painters');
saveas(gcf, fullfile(folder.vectors, [comparisonType '_DSIExampleRuns.svg']));

% Summary Plots
comparisonType2 = [comparisonType '_dsi'];
generateSummaryMetrics(dsiValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, comparisonLabel, comparisonType2, folder);
generateSummaryEVT(nPV, dsiValues, evt_avg, posBins, comparisonLabel, comparisonType2, folder);
end
