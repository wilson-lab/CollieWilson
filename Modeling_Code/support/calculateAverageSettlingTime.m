% CALCULATEAVERAGESETTLINGTIME Compute the average settling time across multiple runs.
%
% This function calculates the time it takes for the object position in
% each run to settle within a specified degree range of the target value (0)
% and remain there for a minimum duration. It then computes the average
% settling time across all runs.
%
% INPUTS:
%   timebase        - 1D array of time points corresponding to columns in visobj_history.
%   visobj_history  - 2D array where rows represent runs, and columns represent time points.
%   toleranceDeg    - Scalar degree tolerance (e.g., 10 for ±10 degrees of 0).
%   minDuration     - Minimum duration (in seconds) the signal must stay within the tolerance.
%
% OUTPUT:
%   avgSettlingTime - Average settling time across all runs. Returns NaN if no settling
%                     time is found for any run.
%
% CREATED: 11/16/2024 - MC

function avgSettlingTime = calculateAverageSettlingTime(timebase, visobj_history, toleranceDeg, minDuration)
    % Calculate the average settling time across multiple runs of a model.
    
    % Number of runs
    numRuns = size(visobj_history, 1);
    % Preallocate array for settling times
    settlingTimes = nan(numRuns, 1);

    % Define the settling range based on the tolerance in degrees
    settlingRange = [-toleranceDeg, toleranceDeg];

    % Calculate the time step from the timebase
    dt = timebase(2) - timebase(1); % Assume uniform time steps
    minPoints = ceil(minDuration / dt); % Convert minDuration to number of points

    % Loop through each run
    for runIdx = 1:numRuns
        % Extract the time-series for this run
        runData = visobj_history(runIdx, :);

        % Check if the data is within the settling range
        withinRange = runData >= settlingRange(1) & runData <= settlingRange(2);

        % Find segments where the signal stays in the range for at least minDuration
        inRangeIndices = find(withinRange);
        if ~isempty(inRangeIndices)
            % Identify the start of continuous segments
            diffIndices = [true, diff(inRangeIndices) == 1]; % Continuity check
            segmentLengths = cumsum(diffIndices);
            segmentLengths(~diffIndices) = 0;
            validSegmentStart = find(segmentLengths >= minPoints, 1, 'first');
            if ~isempty(validSegmentStart)
                % Record the time of the first valid segment
                settlingTimes(runIdx) = timebase(inRangeIndices(validSegmentStart));
            end
        end
    end

    % Calculate the average settling time (ignoring NaNs)
    avgSettlingTime = mean(settlingTimes, 'omitnan');
    
    % If no runs have a valid settling time, return NaN
    if all(isnan(settlingTimes))
        avgSettlingTime = NaN;
    end
end
