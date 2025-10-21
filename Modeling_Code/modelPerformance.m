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
%   comparisonType    - Type of comparison to perform: 'synapse', 'aotu019/025silence', or 'speed'.
%                       - 'synapse' compares inhibitory vs excitatory synapse types.
%                       - 'aotu019/025silence' compares two different strength levels for the model.
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
%          07/30/2025 - MC expanded silencing model
%

function modelPerformance(predicted_RF, comparisonType)
%% Initialize

% Refresh settings
[folder, plotSettings, runSettings] = modelSettings();
close all;

% Define the range to test
kValues = runSettings.k;  % Steering gain (k) values to iterate over
nK = length(kValues);  % Number of k values
stability_noiseLevel = 0.7;
step_noiseLevel = 0.1;
stability_startPos = 0;
step_startPos = 100;
stability_simDuration = 60;
step_simDuration = 11;
settling_tolerance = 2; % degrees from 0
conditionColors = {"k"; "#6ca2e3";"#963977";"#e88598";"#ffc800"};

thisSynapse = "inhibitory";  % Fixed

% Variables depending on comparison type
if strcmp(comparisonType, 'dirselective')
    % compare no DS, 019 DS, and 025 DS
    comparisonLabel = {'none', 'selective','flipped'};
    nComp = 3;
elseif strcmp(comparisonType, 'silence')
    % compare normal, 019 silenced, 025 silenced, 019/025 silenced, minor silenced
    comparisonLabel = {'full', 'AOTU019 0'};
    nComp = 2;
end

% Run model
% Preallocate to store performance metrics and binned averages
metrics_prob = zeros(nK, nComp);
metrics_var = zeros(nK, nComp);
metrics_ISE = zeros(nK, nComp);
metrics_IAE = zeros(nK, nComp);
metrics_dirChangeTime = zeros(nK, nComp);

% Preallocate to store binned averages
all_binned_hist = [];
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
        % Fetch clean copy of RFs
        thisTuning = predicted_RF;

        if strcmp(comparisonType, 'dirselective') %direction selective model
            runSettings.dirselective = idx-1;
        elseif strcmp(comparisonType, 'silence')
            switch idx
                case 2 %silence 019
                    thisTuning.AOTU019 = thisTuning.AOTU019 .* 0;
                case 3 %silence 025
                    thisTuning.AOTU025 = thisTuning.AOTU025 .* 0;
                case 4 %silence 019/025
                    thisTuning.AOTU019 = thisTuning.AOTU019 .* 0;
                    thisTuning.AOTU025 = thisTuning.AOTU025 .* 0;
                case 5 %silence minor
                    thisTuning.sum = thisTuning.sum .* 0;
            end
        else
            close all
            error('Comparison not recognized! Check input to modelPerformance.')
        end

        % Run the AOTU steering model - stability
        [timebase, visobj_history, input_history, rotvel_history] = aotu_steering_model(thisTuning, stability_noiseLevel, stability_startPos, thisSynapse, stability_simDuration, thisK, 2, runSettings);
        % Run the AOTU steering model - step
        startTime = 100;
        [step_timebase, step_visobj_history, ~, ~] = aotu_steering_model(thisTuning, step_noiseLevel, step_startPos, thisSynapse, step_simDuration, thisK, startTime, runSettings);

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

        % Calculate settling time
        avgSettlingTime(kIdx, idx) = calculateAverageSettlingTime(step_timebase, step_visobj_history, settling_tolerance,1)-step_timebase(startTime);

        % Calculate input histogram
        bins = 36;
        binned_hist = input_histogram(input_history, bins);

        % Store the binned averages separately for each comparison
        rotvel_binned_avg(kIdx, :, idx) = rotvel_binned_avgs;
        objvel_binned_avg(kIdx, :, idx) = objvel_binned_avgs;
        evt_avg(kIdx, :, idx) = posvang;
        all_binned_hist(:,:,idx) = binned_hist;

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
sgtitle([' Example Runs for ' comparisonType ' comparison...']);
saveas(gcf, fullfile(folder.final, [comparisonType '_ExampleRuns' '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, [comparisonType '_ExampleRuns' '.svg']));

%% Generate summary plots
% Generate and save summary metrics plot
generateSummaryMetrics(kValues, metrics_prob, comparisonLabel, comparisonType, 'Gain', folder);

% Generate and save summary plot of Angular Velocity vs Object Position (EVT)
generateSummaryEVT(nK, kValues, evt_avg, posBins, comparisonLabel, comparisonType, folder);

% Generate and save summary plot for settling time
generateSettlingSummary(kValues, avgSettlingTime, comparisonLabel, comparisonType, folder)

%% Generate and save summary plot for activity
plot_binned_histograms(all_binned_hist, comparisonLabel)

end
