% removeIsolatedPoints
% This function removes isolated points or small sequences of non-NaN values from the input data.
% It identifies continuous sequences of non-NaN values and retains only the largest sequence 
% in each column that meets a minimum length requirement. All other values are set to NaN.
%
% INPUTS:
%   data               - Matrix of data where isolated points or small sequences should be removed
%   minSequenceLength   - Minimum length of a sequence to be kept (default is 2 if not provided)
%
% OUTPUTS:
%   cleanedData         - Matrix with only the largest sequence of non-NaN values retained per column
%
% CREATED: [Date] MC
%
function cleanedData = removeIsolatedPoints(data, minSequenceLength)
% Check if minimum sequence length is provided, otherwise set default
if nargin < 2
    minSequenceLength = 2;  % Default value if not provided
end

% Get the size of the input data
[numRows, numCols] = size(data);

% Initialize the output as a copy of the input data
cleanedData = data;

% Loop through each column (condition)
for col = 1:numCols
    % Get the current column's data
    currentData = data(:, col);

    % Find where the data is not NaN
    nonNan = ~isnan(currentData);

    % Find the difference between consecutive elements in the logical array
    diffNonNan = diff([0; nonNan; 0]);

    % Find the start and end of sequences of non-NaN values
    startIdx = find(diffNonNan == 1);
    endIdx = find(diffNonNan == -1) - 1;

    % If no non-NaN values, continue to next column
    if isempty(startIdx)
        cleanedData(:, col) = nan;
        continue;
    end

    % Calculate the lengths of all sequences
    chunkLengths = endIdx - startIdx + 1;

    % Find the largest chunk
    [maxChunkLength, largestChunkIdx] = max(chunkLengths);

    % If the largest chunk is greater than or equal to the minimum sequence length
    if maxChunkLength >= minSequenceLength
        % Keep only the largest chunk, set everything else to NaN
        cleanedData(:, col) = nan; % Set the entire column to NaN first
        cleanedData(startIdx(largestChunkIdx):endIdx(largestChunkIdx), col) = currentData(startIdx(largestChunkIdx):endIdx(largestChunkIdx));
    else
        % If no chunk is greater than or equal to the minimum sequence length, set everything to NaN
        cleanedData(:, col) = nan;
    end
end
end
