% SETPOINT_CROSSINGS - This function analyzes the time between direction changes and object positions
% during fixation and nonfixation periods across different conditions in a behavioral experiment. 
% It calculates the median time between direction changes and the median object position at each 
% crossing, separately for fixation and nonfixation periods. An optional swarm plot can be generated 
% to visualize the data.

% INPUTS:
%   panelps      - 3D array of panel position data, where each slice represents a condition.
%   angular      - 3D array of angular velocity data representing turn behavior.
%   forward      - 3D array of forward velocity data.
%   fixationIdx  - 3D logical array indicating whether the fly was fixating (1) or not (0).
%   ttime        - Time vector representing the temporal resolution of the data.
%   settings     - Structure containing parameters like run thresholds and gain conditions.

% OUTPUTS:
%   medianFixationTimes     - Median time between direction changes during fixation for each condition.
%   medianNonFixationTimes  - Median time between direction changes during nonfixation for each condition.
%   medianFixationObjPos    - Median object position during direction changes in fixation periods.
%   medianNonFixationObjPos - Median object position during direction changes in nonfixation periods.

% Created: N/A by MC
% Updated: N/A

% The function:
% - Identifies direction changes (zero-crossings) in angular velocity during fixation and nonfixation periods.
% - Calculates the time between direction changes and the corresponding object position for each condition.
% - Returns median values for both fixation and nonfixation periods.
% - Optionally generates a swarm plot showing the distribution of times and object positions across gain conditions.
%
function [medianFixationTimes, medianNonFixationTimes, medianFixationObjPos, medianNonFixationObjPos] = setpoint_crossings(panelps, angular, forward, fixationIdx, ttime, settings)
%% Initialization
% Get the number of conditions (z-dimension) and the number of trials (y-dimension)
[durTrial, numTrials, numConditions] = size(angular);

% Initialize arrays to store the median times and object positions for each condition
medianFixationTimes = zeros(1, numConditions);
medianNonFixationTimes = zeros(1, numConditions);
medianFixationObjPos = zeros(1, numConditions);
medianNonFixationObjPos = zeros(1, numConditions);

% Initialize arrays to store swarm data for plotting
swarmDataFixationTimes = [];
swarmDataFixationObjPos = [];
swarmDataNonFixationTimes = [];
swarmDataNonFixationObjPos = [];
gainLabelsFixation = [];
gainLabelsNonFixation = [];

% Retrieve the gains from settings
gains = settings.pursuitGain;
plotSwarm = 0;

%% Threshold for when the fly was running (nonfixation)
runIdx = zeros(durTrial, numTrials, numConditions);
for cond = 1:numConditions
    runIdx(:, :, cond) = schmittTrigger(forward(:, :, cond), settings.runThreshB, 0.1); % Threshold for running
end

%% Running the analysis
% Loop through each condition
for cond = 1:numConditions
    % Extract the data for the current condition across all trials
    angularCond = angular(:,:,cond);
    fixationCond = fixationIdx(:,:,cond);
    nonFixationCond = runIdx(:,:,cond);
    panelpsCond = panelps(:,:,cond);

    % Initialize arrays to store times and object positions for all trials in this condition
    allFixationTimes = [];
    allNonFixationTimes = [];
    allFixationObjPos = [];
    allNonFixationObjPos = [];

    % Loop through each trial (y-dimension)
    for trial = 1:numTrials
        % Extract the angular, fixation, nonfixation, and panel position data for this trial
        angularTrial = angularCond(:, trial);
        fixationTrial = fixationCond(:, trial);
        nonFixationTrial = nonFixationCond(:, trial);
        panelpsTrial = panelpsCond(:, trial);

        % Fetch angular and panel data for fixation and nonfixation periods
        angularFixation = angularTrial;
        angularNonFixation = angularTrial;

        % Set data for non-fixation periods to NaN in the fixation dataset
        angularFixation(fixationTrial == 0) = NaN;

        % Set data for non-running periods to NaN in the nonfixation dataset
        angularNonFixation(nonFixationTrial == 0 | fixationTrial == 1) = NaN; % Exclude fixation periods in nonfixation

        % Identify where direction changes occur for fixation periods
        changePointsFixation = find(diff(sign(angularFixation)) ~= 0);

        % Loop through direction changes for fixation to calculate time differences and fetch object position
        for i = 1:length(changePointsFixation)-1
            % Ensure no NaNs between adjacent direction changes
            if all(~isnan(angularFixation(changePointsFixation(i):changePointsFixation(i+1))))
                % Calculate time between current and next direction change (convert to msec)
                deltaTime = (ttime(changePointsFixation(i+1)) - ttime(changePointsFixation(i))) * 1000; % convert to msec
                allFixationTimes = [allFixationTimes, deltaTime];

                % Fetch the absolute object position at the time of the direction change
                objPos = abs(panelpsTrial(changePointsFixation(i)));
                allFixationObjPos = [allFixationObjPos, objPos];
            end
        end

        % Identify where direction changes occur for nonfixation periods
        changePointsNonFixation = find(diff(sign(angularNonFixation)) ~= 0);

        % Loop through direction changes for nonfixation to calculate time differences and fetch object position
        for i = 1:length(changePointsNonFixation)-1
            % Ensure no NaNs between adjacent direction changes
            if all(~isnan(angularNonFixation(changePointsNonFixation(i):changePointsNonFixation(i+1))))
                % Calculate time between current and next direction change (convert to msec)
                deltaTime = (ttime(changePointsNonFixation(i+1)) - ttime(changePointsNonFixation(i))) * 1000; % convert to msec
                allNonFixationTimes = [allNonFixationTimes, deltaTime];

                % Fetch the absolute object position at the time of the direction change
                objPos = abs(panelpsTrial(changePointsNonFixation(i)));
                allNonFixationObjPos = [allNonFixationObjPos, objPos];
            end
        end

    end

    % Calculate the median times and object positions across trials for this condition
    medianFixationTimes(cond) = median(allFixationTimes);
    medianNonFixationTimes(cond) = median(allNonFixationTimes);
    medianFixationObjPos(cond) = median(allFixationObjPos);
    medianNonFixationObjPos(cond) = median(allNonFixationObjPos);

    % Prepare swarm data for plotting (for both object positions and times)
    if plotSwarm
        % Store the times and object positions for swarm chart (fixation)
        swarmDataFixationTimes = [swarmDataFixationTimes, allFixationTimes];
        swarmDataFixationObjPos = [swarmDataFixationObjPos, allFixationObjPos];
        gainLabelsFixation = [gainLabelsFixation, repmat(gains(cond), 1, length(allFixationTimes))];

        % Store the times and object positions for swarm chart (nonfixation)
        swarmDataNonFixationTimes = [swarmDataNonFixationTimes, allNonFixationTimes];
        swarmDataNonFixationObjPos = [swarmDataNonFixationObjPos, allNonFixationObjPos];
        gainLabelsNonFixation = [gainLabelsNonFixation, repmat(gains(cond), 1, length(allNonFixationTimes))];
    end
end

% Plot the swarm chart (if requested)
if plotSwarm
    figure; set(gcf,'Position',[100 100 1500 900]);
    tiledlayout(2, 2);  % Two rows, two columns

    % Determine shared y-limits for both time plots and both object position plots
    allTimeData = [swarmDataFixationTimes, swarmDataNonFixationTimes]; % Combine time datasets to find common y-limits
    yLimitsTime = [0 max(allTimeData)];

    allObjPosData = [swarmDataFixationObjPos, swarmDataNonFixationObjPos]; % Combine object position datasets
    yLimitsObjPos = [0 max(allObjPosData)];

    % Define bar width for the median
    barWidth = 2; % You can adjust this to control the width of the bars

    % Plot fixation times
    nexttile;
    swarmchart(gainLabelsFixation, swarmDataFixationTimes, 10);  % Create swarm chart for fixation times
    hold on;
    for cond = 1:length(gains)
        bar(gains(cond), medianFixationTimes(cond), barWidth, 'FaceColor', 'none', 'EdgeColor', 'r'); % Plot bar at median
    end
    title('Fixation Times Between Direction Changes');
    xlabel('Gain Condition');
    ylabel('Time Between Direction Changes (ms)');
    xticks(gains);
    xticklabels(gains);
    ylim(yLimitsTime);  % Apply shared y-limits for times
    hold off;

    % Plot nonfixation times
    nexttile;
    swarmchart(gainLabelsNonFixation, swarmDataNonFixationTimes, 10);  % Create swarm chart for nonfixation times
    hold on;
    for cond = 1:length(gains)
        bar(gains(cond), medianNonFixationTimes(cond), barWidth, 'FaceColor', 'none', 'EdgeColor', 'r'); % Plot bar at median
    end
    title('Nonfixation Times Between Direction Changes');
    xlabel('Gain Condition');
    ylabel('Time Between Direction Changes (ms)');
    xticks(gains);
    xticklabels(gains);
    ylim(yLimitsTime);  % Apply shared y-limits for times
    hold off;

    % Plot object positions at time of crossing for fixation
    nexttile;
    swarmchart(gainLabelsFixation, swarmDataFixationObjPos, 10);  % Create swarm chart for fixation object positions
    hold on;
    for cond = 1:length(gains)
        bar(gains(cond), medianFixationObjPos(cond), barWidth, 'FaceColor', 'none', 'EdgeColor', 'r'); % Plot bar at median
    end
    title('Object Positions at Fixation Crossings');
    xlabel('Gain Condition');
    ylabel('Absolute Object Position');
    xticks(gains);
    xticklabels(gains);
    ylim(yLimitsObjPos);  % Apply shared y-limits for object positions
    hold off;

    % Plot object positions at time of crossing for nonfixation
    nexttile;
    swarmchart(gainLabelsNonFixation, swarmDataNonFixationObjPos, 10);  % Create swarm chart for nonfixation object positions
    hold on;
    for cond = 1:length(gains)
        bar(gains(cond), medianNonFixationObjPos(cond), barWidth, 'FaceColor', 'none', 'EdgeColor', 'r'); % Plot bar at median
    end
    title('Object Positions at Nonfixation Crossings');
    xlabel('Gain Condition');
    ylabel('Absolute Object Position');
    xticks(gains);
    xticklabels(gains);
    ylim(yLimitsObjPos);  % Apply shared y-limits for object positions
    hold off;
end


end
