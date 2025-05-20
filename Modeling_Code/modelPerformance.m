% modelPerformance
%
% This function runs and compares the performance of a fly steering model
% across different comparison types. It simulates fly behavior across a
% range of steering gain (k) values and computes performance metrics for each 
% comparison condition. The results are stored and visualized to compare how
% the model behaves under different conditions.
%
% INPUTS:
%   predicted_RF      - Struct containing the predicted receptive fields (RF) for neurons in the model.
%                       These RFs define the neural responses to visual stimuli in the azimuthal space.
%   comparisonType    - Type of comparison to perform: 'synapse', 'strength', or 'speed'.
%                       - 'synapse' compares inhibitory vs excitatory synapse types.
%                       - 'strength' compares two different strength levels for the model.
%                       - 'speed' compares fast and slow visuomotor delays.
%
% DESCRIPTION:
% The function initializes settings, runs the steering model for each condition (synapse type, strength, or speed), 
% and calculates several performance metrics (probability, variance, ISE, IAE, direction change times) across a range 
% of steering gain values (`k`). The results are plotted and saved for further comparison.
%
% CREATED: 10/30/2024 - MC
% UPDATED: 11/16/2024 - MC added step function for measuring settling time
%                          added direction selectivity model
%

function modelPerformance(predicted_RF, comparisonType)
%% Initialize

% Refresh settings
[folder, plotSettings, runSettings] = modelSettings();
close all;

% Define the range to test
kValues = [0,1.2];  % Steering gain (k) values to iterate over
nK = length(kValues);  % Number of k values
stability_noiseLevel = 0.7;
step_noiseLevel = 0.1;
stability_startPos = 0;
step_startPos = 100;
stability_simDuration = 60;
step_simDuration = 11;
settling_tolerance = 2; % degrees from 0
conditionColors = {"#0072BD"; "#D95319";"#7E2F8E"};

% Variables depending on comparison type
if strcmp(comparisonType, 'synapse')
    comparisonLabel = {'Inhibitory', 'Excitatory'};
    nComp = 2;
    thisSynapse = {"inhibitory"; "excitatory"};
    strengthValues = 1;  % Fixed
elseif strcmp(comparisonType, 'strength')
    comparisonLabel = {'Strength 1', 'Strength 0'};
    nComp = 2;
    thisSynapse = "inhibitory";  % Fixed
    strengthValues = [1, 0];
elseif strcmp(comparisonType, 'dirselective')
    comparisonLabel = {'none', 'selective','flipped'};
    nComp = 3;
    thisSynapse = "inhibitory";  % Fixed
    strengthValues = 1;  % Fixed
end

% Run model
% Preallocate to store performance metrics and binned averages
metrics_prob = zeros(nK, nComp);
metrics_var = zeros(nK, nComp);
metrics_ISE = zeros(nK, nComp);
metrics_IAE = zeros(nK, nComp);
metrics_dirChangeTime = zeros(nK, nComp);

% Preallocate to store binned averages
indist_avg = [];
rotvel_binned_avg = [];
objvel_binned_avg = [];
evt_avg = [];
avgSettlingTime = [];

% Select up to X evenly spaced k values for plotting
maxPlotK = min(6, nK);  % Limit to X `k` runs for plotting
selected_k_indices = round(linspace(1, nK, maxPlotK));  % Evenly spaced indices from `kValues`

% Initialize figure for example plots
figure; set(gcf, 'Position', [100 100 1500 900]);  % Set figure size
tiledlayout(maxPlotK, 4*nComp, 'TileSpacing', 'compact'); 

% Loop over k values (steering gain) for analysis
for kIdx = 1:nK
    thisK = kValues(kIdx);
    runSettings.nK = kIdx;
    % Output progress
    disp(['Running simulations for k = ', num2str(thisK)]);

    % Loop over comparison values
    for idx = 1:nComp
        if strcmp(comparisonType, 'synapse')
            thisType = thisSynapse{idx};
            thisStrength = strengthValues;  % Fixed strength for synapse comparison
        elseif strcmp(comparisonType, 'strength')
            thisStrength = strengthValues(idx);
            thisType = thisSynapse;  % Fixed synapse type for strength comparison
        elseif strcmp(comparisonType, 'dirselective')
            thisStrength = strengthValues; % Fixed
            thisType = thisSynapse;  % Fixed
            runSettings.dirselective = idx-1; %vary whether dsi included
        end

        % Adjust tuning for strength if applicable
        thisTuning = predicted_RF;
        thisTuning.AOTU019 = thisTuning.AOTU019 .* thisStrength;  % Adjust tuning by strength

        % Run the AOTU steering model - stability
        [timebase, visobj_history, input_history, rotvel_history] = aotu_steering_model(thisTuning, stability_noiseLevel, stability_startPos, thisType, stability_simDuration, thisK, 2, runSettings);
        % Run the AOTU steering model - step
        startTime = 100;
        [step_timebase, step_visobj_history, ~, ~] = aotu_steering_model(thisTuning, step_noiseLevel, step_startPos, thisType, step_simDuration, thisK, startTime, runSettings);

        % Calculate performance metrics
        metrics_results = calculatePerformanceMetrics(visobj_history, rotvel_history, timebase, plotSettings);
        % Store metrics for this comparison and k value
        metrics_prob(kIdx, idx) = metrics_results.prob;
        metrics_var(kIdx, idx) = metrics_results.var;
        metrics_ISE(kIdx, idx) = metrics_results.ISE;
        metrics_IAE(kIdx, idx) = metrics_results.IAE;

        % Measure performance at direction changes and store binned averages
        nTest = runSettings.numRuns;
        [rotvel_cross_times, rotvel_bins, rotvel_binned_avgs,objvel_binned_avgs] = analyzeObjectAndVelocityCrossings(visobj_history, rotvel_history, timebase, nTest);
        metrics_dirChangeTime(kIdx, idx) = mean(rotvel_cross_times, 'omitnan');

        % Analyze object position vs velocity
        [posvang, posBins] = analyzeErrorVsTurn(visobj_history, rotvel_history, runSettings); 
        
        % Compute normalized distribution of input history
        [input_distribution, inBins] = compute_normalized_distribution(input_history, runSettings);

        % Calculate settling time
        avgSettlingTime(kIdx, idx) = calculateAverageSettlingTime(step_timebase, step_visobj_history, settling_tolerance,1)-step_timebase(startTime);

        % Store the binned averages separately for each comparison
        rotvel_binned_avg(kIdx, :, idx) = rotvel_binned_avgs;
        objvel_binned_avg(kIdx, :, idx) = objvel_binned_avgs;
        evt_avg(kIdx, :, idx) = posvang;
        indist_avg(kIdx,:, idx) = input_distribution;

        % Plot example run
        if any(kIdx == selected_k_indices)
            % Call the function to remove large jumps
            visobj_plot = remove_large_jumps(visobj_history(1, :), 180);
            step_visobj_plot = remove_large_jumps(step_visobj_history(1, :), 180);
            
            nexttile([1 3])
            plot(timebase, visobj_plot, 'Color', conditionColors{idx});
            xlabel('Time (s)');
            ylabel('Pos');
            title(['k = ' num2str(thisK)]);
            axis tight;
            if kIdx == 1
                legend(comparisonLabel{idx});
            end
            yline(0)
            ylim([-180 180]);

            nexttile
            plot(step_timebase, step_visobj_plot, 'Color', conditionColors{idx});
            xlabel('Time (s)');
            ylabel('Pos');
            title(['k = ' num2str(thisK)]);
            axis tight;
            yline(0)
            xline(step_timebase(startTime))
            ylim([-180 180]);
        end
    end
end

%% Save the example plot comparing the two conditions across selected k values
sgtitle([' Example Runs for ' comparisonLabel{1} ' vs ' comparisonLabel{2}]);
saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_ExampleRuns' '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_ExampleRuns' '.svg']));

%% Generate summary plots
% Generate and save summary metrics plot
generateSummaryMetrics(kValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, comparisonLabel, folder);

% Generate and save summary plot of binned crossing times
%generateSummaryCrossing(nK, kValues, rotvel_binned_avg, objvel_binned_avg, rotvel_bins, comparisonLabel, folder);

% Generate and save summary plot of Angular Velocity vs Object Position (EVT)
generateSummaryEVT(nK, kValues, evt_avg, posBins, comparisonLabel,folder);

% Generate and save summary plot for steering drive
%generateSteeringPlot(steering_drive, comparisonLabel, runSettings,folder)

% Generate and save summary plot for inputs to DNa02
%generateSummaryInputDistributions(nK, kValues, indist_avg, inBins, comparisonLabel, runSettings, folder)

% Generate and save summary plot for settling time
generateSettlingSummary(kValues, avgSettlingTime, comparisonLabel, folder)

end
