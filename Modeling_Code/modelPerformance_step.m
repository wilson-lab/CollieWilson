% modelPerformance_step
%
% This function runs a step input test to evaluate the settling performance
% of a fly steering model across a range of initial stimulus positions.
% It compares model behavior under different experimental conditions
% (e.g., synapse type, strength, or direction selectivity), computes
% settling times for each condition, and generates example response plots
% at selected start positions.
%
% INPUTS:
%   predicted_RF      - Struct containing the predicted receptive fields for model neurons.
%   comparisonType    - Comparison mode: 'synapse', 'strength', or 'dirselective'.
%
% CREATED: 04/06/2025 - MC

function modelPerformance_step(predicted_RF, comparisonType)
%% Initialize

[folder, plotSettings, runSettings] = modelSettings();
close all;

% Define the range of step start positions to test
startPosValues = 0:10:130;
nStart = length(startPosValues);
kVal = 1.2; % fixed steering gain
step_noiseLevel = 0.1;
step_simDuration = 11;
startTime = 100;
settling_tolerance = 2; % degrees from 0
conditionColors = {"#0072BD"; "#D95319";"#7E2F8E"};

% Variables depending on comparison type
if strcmp(comparisonType, 'synapse')
    comparisonLabel = {'Inhibitory', 'Excitatory'};
    nComp = 2;
    thisSynapse = {"inhibitory"; "excitatory"};
    strengthValues = 1;
elseif strcmp(comparisonType, 'strength')
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

% Preallocate settling time array
avgSettlingTime = zeros(nStart, nComp);

% Select subset of start positions to plot
maxPlotS = min(6, nStart);
selected_start_indices = round(linspace(1, nStart, maxPlotS));

% Initialize figure for example plots
figure; set(gcf, 'Position', [100 100 1500 900]);
tiledlayout(maxPlotS, nComp, 'TileSpacing', 'compact');

% Loop over start positions
for sIdx = 1:nStart
    thisStartPos = startPosValues(sIdx);
    runSettings.nStart = sIdx;
    disp(['Running step test for startPos = ', num2str(thisStartPos)]);

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

        % Adjust tuning
        thisTuning = predicted_RF;
        thisTuning.AOTU019 = thisTuning.AOTU019 .* thisStrength;

        % Run step simulation
        [step_timebase, step_visobj_history, ~, ~] = aotu_steering_model(...
            thisTuning, step_noiseLevel, thisStartPos, thisType, ...
            step_simDuration, kVal, startTime, runSettings);

        % Compute settling time
        avgSettlingTime(sIdx, idx) = calculateAverageSettlingTime(...
            step_timebase, step_visobj_history, settling_tolerance, 1) ...
            - step_timebase(startTime);

        % Example plots
        if any(sIdx == selected_start_indices)
            step_visobj_plot = remove_large_jumps(step_visobj_history(1, :), 180);
            nexttile
            plot(step_timebase, step_visobj_plot, 'Color', conditionColors{idx});
            xlabel('Time (s)');
            ylabel('Position (deg)');
            title(['Start Pos = ' num2str(thisStartPos)]);
            yline(0);
            xline(step_timebase(startTime));
            ylim([-180 180]);
            axis tight;

            if sIdx == 1
                legend(comparisonLabel{idx});
            end
        end
    end
end

%% Save the figure of example runs
sgtitle(['Example Step Responses: ' comparisonLabel{1} ' vs ' comparisonLabel{2}]);
saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_StepExampleRuns.png']));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_StepExampleRuns.svg']));

%% Generate summary plot
generateSettlingSummary(startPosValues, avgSettlingTime, comparisonLabel, folder);

end
