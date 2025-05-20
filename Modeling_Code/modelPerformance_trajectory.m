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

function modelPerformance_trajectory(predicted_RF, comparisonType)
%% Initialize

[folder, plotSettings, runSettings] = modelSettings();
close all;

% Define range of noise levels to test
noiseValues = 0:0.1:1.2;
nNoise = length(noiseValues);
kVal = 1.2; % fixed steering gain

stability_startPos = 0;
stability_simDuration = 60;
conditionColors = {"#0072BD"; "#D95319";"#7E2F8E"};

% Variables depending on comparison type
if strcmp(comparisonType, 'strength')
    comparisonLabel = {'Strength 1', 'Strength 0'};
    nComp = 2;
    thisSynapse = "inhibitory";
    strengthValues = [1, 0];
elseif strcmp(comparisonType, 'dirselective')
    comparisonLabel = {'none', 'selective','flipped'};
    nComp = 3;
    thisSynapse = "inhibitory";
    strengthValues = 1;
end

% Preallocate performance metric arrays
metrics_prob = zeros(nNoise, nComp);
metrics_var = zeros(nNoise, nComp);
metrics_ISE = zeros(nNoise, nComp);
metrics_IAE = zeros(nNoise, nComp);
metrics_dirChangeTime = zeros(nNoise, nComp);

indist_avg = [];
rotvel_binned_avg = [];
objvel_binned_avg = [];
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
        if strcmp(comparisonType, 'synapse')
            thisType = thisSynapse{idx};
            thisStrength = strengthValues;
        elseif strcmp(comparisonType, 'strength')
            thisStrength = strengthValues(idx);
            thisType = thisSynapse;
        elseif strcmp(comparisonType, 'dirselective')
            thisStrength = strengthValues;
            thisType = thisSynapse;
            runSettings.dirselective = idx - 1;
        end

        thisTuning = predicted_RF;
        thisTuning.AOTU019 = thisTuning.AOTU019 .* thisStrength;

        % Run model (stability only)
        [timebase, visobj_history, input_history, rotvel_history] = ...
            aotu_steering_model(thisTuning, thisNoise, stability_startPos, ...
            thisType, stability_simDuration, kVal, 2, runSettings);

        % Performance metrics
        metrics_results = calculatePerformanceMetrics(visobj_history, rotvel_history, timebase, plotSettings);
        metrics_prob(nIdx, idx) = metrics_results.prob;
        metrics_var(nIdx, idx) = metrics_results.var;
        metrics_ISE(nIdx, idx) = metrics_results.ISE;
        metrics_IAE(nIdx, idx) = metrics_results.IAE;

        % Timing and distribution analyses
        nTest = runSettings.numRuns;
        [rotvel_cross_times, rotvel_bins, rotvel_binned_avgs, objvel_binned_avgs] = ...
            analyzeObjectAndVelocityCrossings(visobj_history, rotvel_history, timebase, nTest);
        metrics_dirChangeTime(nIdx, idx) = mean(rotvel_cross_times, 'omitnan');

        [posvang, posBins] = analyzeErrorVsTurn(visobj_history, rotvel_history, runSettings);
        [input_distribution, inBins] = compute_normalized_distribution(input_history, runSettings);

        % Store results
        rotvel_binned_avg(nIdx, :, idx) = rotvel_binned_avgs;
        objvel_binned_avg(nIdx, :, idx) = objvel_binned_avgs;
        evt_avg(nIdx, :, idx) = posvang;
        indist_avg(nIdx,:, idx) = input_distribution;

        % Plot example run for selected noise levels
        if any(nIdx == selected_noise_indices)
            visobj_plot = remove_large_jumps(visobj_history(1, :), 180);

            nexttile
            plot(timebase, visobj_plot, 'Color', conditionColors{idx});
            xlabel('Time (s)');
            ylabel('Position (deg)');
            title(['Noise = ' num2str(thisNoise)]);
            yline(0);
            ylim([-180 180]);
            axis tight;

            if nIdx == 1
                legend(comparisonLabel{idx});
            end
        end
    end
end

%% Save the example plot figure
sgtitle(['Example Runs for ' comparisonLabel{1} ' vs ' comparisonLabel{2} ' across noise levels']);
saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_NoiseExampleRuns.png']));
set(gcf,'renderer','Painters');
saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_NoiseExampleRuns.svg']));

%% Summary Plots
generateSummaryMetrics(noiseValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, comparisonLabel, folder);
generateSummaryEVT(nNoise, noiseValues, evt_avg, posBins, comparisonLabel, folder);

end
