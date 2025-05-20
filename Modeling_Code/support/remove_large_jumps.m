function data_processed = remove_large_jumps(data, threshold)
% remove_large_jumps
% This function sets values to NaN where the difference between adjacent points
% in the input `data` exceeds the specified `threshold`, across multiple trials.
%
% INPUTS:
%   data      - 2D array of data points (time × trials)
%   threshold - Threshold for detecting large jumps (e.g., 180 degrees)
%
% OUTPUT:
%   data_processed - Modified data where large jumps are replaced with NaN.

% Calculate the difference along time (rows)
diff_data = abs(diff(data, 1, 1));

% Initialize output
data_processed = data;

% Identify where the jump exceeds threshold
jump_mask = diff_data > threshold;

% Pad with a row of false to align with data (diff reduces length by 1)
jump_mask = [false(1, size(data, 2)); jump_mask];

% Set the offending point to NaN
data_processed(jump_mask) = NaN;
end
