function modelPerformance_trajectory(predicted_RF, comparisonType)
% MODELPERFORMANCE_TRAJECTORY
% CREATED: 10/20/2025 - MC
%
% Runs and evaluates the fly steering model across noise or DSI levels to assess
% how variability and direction selectivity affect trajectory stability.
% Computes performance metrics, example trajectories, and summary plots for
% different comparison modes (e.g., 'strength' or 'dirselective').

%% INITIALIZE PARAMETERS AND SETTINGS
[folder, plotSettings, runSettings] = modelSettings();
close all;

switch lower(comparisonType)
    case 'strength'
        compLabels = {'Strength 1','Strength 0'};
        nComp = 2;
    case 'dirselective'
        compLabels = {'neither','019 DS','025 DS'};
        nComp = 3;
    otherwise
        compLabels = {'Cond A','Cond B'};
        nComp = 2;
end

thisSynapse = "inhibitory";
kVal        = runSettings.k;
condColors  = {"#0072BD"; "#D95319"; "#7E2F8E"};

%% RUN SIMULATIONS ACROSS NOISE LEVELS
noiseValues = 0:0.05:0.7;
nNoise      = numel(noiseValues);
maxShow     = min(6, nNoise);
showIdx     = round(linspace(1, nNoise, maxShow));

metrics_prob = zeros(nNoise, nComp);
metrics_var  = zeros(nNoise, nComp);
metrics_ISE  = zeros(nNoise, nComp);
metrics_IAE  = zeros(nNoise, nComp);
avg_obj_speed = zeros(nNoise, nComp);
evt_avg = [];
posBins = [];

startPos   = 0;
simDur     = 60;

close all
figure; set(gcf, 'Position', [100 100 1500 900]);
tiledlayout(maxShow, nComp, 'TileSpacing', 'compact');

% ---- Loop through noise levels ----
for nIdx = 1:nNoise
    thisNoise = noiseValues(nIdx);
    runSettings.nNoise = nIdx;

    for cIdx = 1:nComp
        % Load comparison-specific tuning
        thisTuning = predicted_RF;
        switch lower(comparisonType)
            case 'strength'
                if cIdx == 2
                    thisTuning.AOTU019 = thisTuning.AOTU019 .* 0;
                end
            case 'dirselective'
                runSettings.dirselective = cIdx - 1;
        end

        % Run steering model
        [t, vispos, ~, rotvel] = aotu_steering_model(thisTuning, thisNoise, ...
            startPos, thisSynapse, simDur, kVal, 2, runSettings);

        % Compute performance metrics
        mr = calculatePerformanceMetrics(vispos, rotvel, t, plotSettings);
        metrics_prob(nIdx, cIdx) = mr.prob;
        metrics_var(nIdx, cIdx)  = mr.var;
        metrics_ISE(nIdx, cIdx)  = mr.ISE;
        metrics_IAE(nIdx, cIdx)  = mr.IAE;
        avg_obj_speed(nIdx, cIdx) = compute_avg_obj_speed(t, vispos);

        % Error–turn analysis
        [posvang, posBins] = analyzeErrorVsTurn(vispos, rotvel, runSettings);
        evt_avg(nIdx, :, cIdx) = posvang;

        % Plot representative trajectories
        if any(nIdx == showIdx)
            trialPlot = remove_large_jumps(vispos(1,:), 180);
            nexttile
            plot(t, trialPlot, 'Color', condColors{cIdx});
            xlabel('Time (s)'); ylabel('Position (deg)');
            title(sprintf('Noise = %.2f', thisNoise));
            yline(0); axis tight; ylim([-50 50]);
            if nIdx == showIdx(1), legend(compLabels{cIdx}); end
        end
    end
end

% ---- Save and summarize noise-level results ----
sgtitle(['Example Runs for ' comparisonType ' across noise levels']);
set(gcf, 'Renderer', 'painters');
saveas(gcf, fullfile(folder.final,  [comparisonType '_NoiseExampleRuns.png']));
saveas(gcf, fullfile(folder.vectors,[comparisonType '_NoiseExampleRuns.svg']));

avg_obj_speed_across = mean(avg_obj_speed,2);
generateSummaryMetrics(avg_obj_speed_across, metrics_prob, compLabels, comparisonType, 'Object speed (deg/s)', folder)
generateSummaryEVT(nNoise, noiseValues, evt_avg, posBins, compLabels, comparisonType, folder);

%% RUN SIMULATIONS ACROSS DSI LEVELS
penaltyValues = 0:0.1:1;
dsiValues     = (1 - penaltyValues) ./ (1 + penaltyValues);
nDSI          = numel(penaltyValues);
maxShow2      = min(6, nDSI);
showIdx2      = round(linspace(1, nDSI, maxShow2));

metrics_prob = zeros(nDSI, nComp);
metrics_var  = zeros(nDSI, nComp);
metrics_ISE  = zeros(nDSI, nComp);
metrics_IAE  = zeros(nDSI, nComp);
evt_avg      = [];
posBins      = [];

thisNoise = 0.15;  % fixed for DSI sweep
close all
figure; set(gcf, 'Position', [100 100 1500 900]);
tiledlayout(maxShow2, nComp, 'TileSpacing', 'compact');

% ---- Loop through DSI values ----
for dIdx = 1:nDSI
    runSettings.AOTU019dsi_penalty = penaltyValues(dIdx);
    thisDSI = dsiValues(dIdx);

    for cIdx = 1:nComp
        % Load comparison-specific tuning
        thisTuning = predicted_RF;
        switch lower(comparisonType)
            case 'strength'
                if cIdx == 2
                    thisTuning.AOTU019 = thisTuning.AOTU019 .* 0;
                end
            case 'dirselective'
                runSettings.dirselective = cIdx - 1;
        end

        % Run steering model
        [t, vispos, ~, rotvel] = aotu_steering_model(thisTuning, thisNoise, ...
            startPos, thisSynapse, simDur, kVal, 2, runSettings);

        % Compute performance metrics
        mr = calculatePerformanceMetrics(vispos, rotvel, t, plotSettings);
        metrics_prob(dIdx, cIdx) = mr.prob;
        metrics_var(dIdx, cIdx)  = mr.var;
        metrics_ISE(dIdx, cIdx)  = mr.ISE;
        metrics_IAE(dIdx, cIdx)  = mr.IAE;

        % Error–turn analysis
        [posvang, posBins] = analyzeErrorVsTurn(vispos, rotvel, runSettings);
        evt_avg(dIdx, :, cIdx) = posvang;

        % Plot representative trajectories
        if any(dIdx == showIdx2)
            trialPlot = remove_large_jumps(vispos(1,:), 180);
            nexttile
            plot(t, trialPlot, 'Color', condColors{cIdx});
            xlabel('Time (s)'); ylabel('Position (deg)');
            title(sprintf('DSI = %.2f', thisDSI));
            yline(0); axis tight; ylim([-50 50]);
            if dIdx == showIdx2(1), legend(compLabels{cIdx}); end
        end
    end
end

% ---- Save and summarize DSI-level results ----
sgtitle(['Example Runs for ' comparisonType ' across DSI levels']);
set(gcf, 'Renderer', 'painters');
saveas(gcf, fullfile(folder.final,  [comparisonType '_DSIExampleRuns.png']));
saveas(gcf, fullfile(folder.vectors,[comparisonType '_DSIExampleRuns.svg']));

comparisonType2 = [comparisonType '_dsi'];
generateSummaryMetrics(dsiValues, metrics_prob, compLabels, comparisonType, 'DSI', folder)
generateSummaryEVT(nDSI, dsiValues, evt_avg, posBins, compLabels, comparisonType2, folder);
end