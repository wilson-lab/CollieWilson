function modelPerformance_step(predicted_RF, comparisonType)
% MODELPERFORMANCE_STEP
% CREATED: 10/20/2025 - MC
%
% Runs a step input test across initial positions to assess settling performance
% under different comparison modes ('synapse', 'strength', 'dirselective').
% Computes settling times, plots representative responses, and generates summaries.

%% INITIALIZE PARAMETERS AND SETTINGS
[folder, plotSettings, runSettings] = modelSettings();
close all;

% Comparison setup
switch lower(comparisonType)
    case 'synapse'
        compLabels = {'Inhibitory','Excitatory'};
        nComp      = 2;
        synapseOpt = {"inhibitory"; "excitatory"};
        strengthOpt = 1;
    case 'strength'
        compLabels  = {'Strength 1','Strength 0'};
        nComp       = 2;
        synapseOpt  = "inhibitory";
        strengthOpt = [1, 0];
    case 'dirselective'
        compLabels  = {'neither','019 DS','025 DS'};
        nComp       = 3;
        synapseOpt  = "inhibitory";
        strengthOpt = 1;
    otherwise
        compLabels  = {'Cond A','Cond B'};
        nComp       = 2;
        synapseOpt  = "inhibitory";
        strengthOpt = 1;
end

% Fixed model parameters
kVal          = runSettings.k;
noiseLevel    = 0.01;
simDuration   = 16;
startIndex    = 100;     % index in timebase at which step occurs
settleTolDeg  = 2;       % degrees from 0 for settling
condColors    = {"#0072BD"; "#D95319"; "#7E2F8E"};

% Sweep over initial step positions
startPosValues = 0:10:140;
nStart         = numel(startPosValues);

% Storage
avgSettlingTime = zeros(nStart, nComp);

% Representative plotting choices
maxShow   = min(6, nStart);
showIdx   = round(linspace(1, nStart, maxShow));

figure; set(gcf, 'Position', [100 100 1500 900]);
tiledlayout(maxShow, nComp, 'TileSpacing', 'compact');

% RUN STEP SIMULATIONS ACROSS START POSITIONS
for sIdx = 1:nStart
    thisStartPos = startPosValues(sIdx);
    runSettings.nStart = sIdx;

    for cIdx = 1:nComp
        % --- Configure synapse type and strength for this condition ---
        if iscell(synapseOpt)
            thisSynapse  = synapseOpt{cIdx};
        else
            thisSynapse  = synapseOpt;
        end
        if numel(strengthOpt) > 1
            thisStrength = strengthOpt(cIdx);
        else
            thisStrength = strengthOpt;
        end
        if strcmpi(comparisonType, 'dirselective')
            runSettings.dirselective = cIdx - 1; % 0,1,2
        end

        % --- Adjust tuning for condition ---
        thisTuning = predicted_RF;
        thisTuning.AOTU019 = thisTuning.AOTU019 .* thisStrength;

        % --- Run step simulation ---
        [t, vispos, ~, ~] = aotu_steering_model( ...
            thisTuning, noiseLevel, thisStartPos, thisSynapse, ...
            simDuration, kVal, startIndex, runSettings);

        % --- Settling time (relative to step onset time) ---
        tSettle = calculateAverageSettlingTime(t, vispos, settleTolDeg, 1);
        avgSettlingTime(sIdx, cIdx) = tSettle - t(startIndex);

        % --- Plot representative responses ---
        if any(sIdx == showIdx)
            trialPlot = remove_large_jumps(vispos(1,:), 180);
            nexttile
            plot(t, trialPlot, 'Color', condColors{cIdx});
            xlabel('Time (s)'); ylabel('Position (deg)');
            title(sprintf('Start Pos = %d°', thisStartPos));
            yline(0); xline(t(startIndex));
            axis tight; ylim([-180 180]);
            if sIdx == showIdx(1), legend(compLabels{cIdx}); end
        end
    end
end

% SAVE FIGURE AND SUMMARY (START POSITION SWEEP)
sgtitle(['Example Step Responses: ' comparisonType]);
set(gcf, 'Renderer', 'painters');
saveas(gcf, fullfile(folder.final,  [comparisonType '_StepExampleRuns.png']));
saveas(gcf, fullfile(folder.vectors,[comparisonType '_StepExampleRuns.svg']));

generateSettlingSummary(startPosValues, avgSettlingTime, compLabels, comparisonType, folder);

%% RUN DSI SWEEP (PENALTY → DSI) AT FIXED START POSITION
[folder, plotSettings, runSettings] = modelSettings(); %#ok<ASGLU>  % refresh settings for clean state
close all;

penaltyValues = 0:0.1:1;
dsiValues     = (1 - penaltyValues) ./ (1 + penaltyValues);
nDSI          = numel(penaltyValues);

fixedStartPos = 120;     % degrees
thisNoise     = noiseLevel;

% Reuse comparison labeling/size from above
avgSettlingTime = zeros(nDSI, nComp);

maxShow2 = min(6, nDSI);
showIdx2 = round(linspace(1, nDSI, maxShow2));

figure; set(gcf, 'Position', [100 100 1500 900]);
tiledlayout(maxShow2, nComp, 'TileSpacing', 'compact');

for dIdx = 1:nDSI
    runSettings.AOTU019dsi_penalty = penaltyValues(dIdx);
    thisDSI = dsiValues(dIdx);

    for cIdx = 1:nComp
        % --- Configure synapse/strength and DS mode per condition ---
        if iscell(synapseOpt)
            thisSynapse  = synapseOpt{cIdx};
        else
            thisSynapse  = synapseOpt;
        end
        if numel(strengthOpt) > 1
            thisStrength = strengthOpt(cIdx);
        else
            thisStrength = strengthOpt;
        end
        if strcmpi(comparisonType, 'dirselective')
            runSettings.dirselective = cIdx - 1; % 0,1,2
        end

        % --- Adjust tuning for condition ---
        thisTuning = predicted_RF;
        thisTuning.AOTU019 = thisTuning.AOTU019 .* thisStrength;

        % --- Run step simulation at fixed start position ---
        [t, vispos, ~, ~] = aotu_steering_model( ...
            thisTuning, thisNoise, fixedStartPos, thisSynapse, ...
            simDuration, kVal, startIndex, runSettings);

        % --- Settling time relative to step onset ---
        tSettle = calculateAverageSettlingTime(t, vispos, settleTolDeg, 1);
        avgSettlingTime(dIdx, cIdx) = tSettle - t(startIndex);

        % --- Plot representative responses ---
        if any(dIdx == showIdx2)
            trialPlot = remove_large_jumps(vispos(1,:), 180);
            nexttile
            plot(t, trialPlot, 'Color', condColors{cIdx});
            xlabel('Time (s)'); ylabel('Position (deg)');
            title(sprintf('DSI = %.2f', thisDSI));
            yline(0); xline(t(startIndex));
            axis tight; ylim([-180 180]);
            if dIdx == showIdx2(1), legend(compLabels{cIdx}); end
        end
    end
end

% SAVE FIGURE AND SUMMARY (DSI SWEEP)
sgtitle(['Example Step Responses: ' comparisonType ' (DSI sweep)']);
set(gcf, 'Renderer', 'painters');
saveas(gcf, fullfile(folder.final,  [comparisonType '_dsiStepExampleRuns.png']));
saveas(gcf, fullfile(folder.vectors,[comparisonType '_dsiStepExampleRuns.svg']));

comparisonType2 = [comparisonType '_dsi'];
generateSettlingSummary(dsiValues, avgSettlingTime, compLabels, comparisonType2, folder);
end
