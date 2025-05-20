% find_rise_time
% This function calculates the rise time of a signal represented by the input trial.
% The rise time is defined as the number of indices between the first occurrence of 
% the signal reaching 10% of its peak value and the first occurrence of the signal 
% reaching 90% of its peak value. If either threshold is not reached, the function returns 
% infinity to indicate an invalid rise time.
%
% INPUTS:
%   trial - Array representing the signal over time (e.g., voltage, position, etc.)
%
% OUTPUTS:
%   rise_time - The rise time in terms of the number of indices; returns 
%               infinity (inf) if the rise time cannot be determined
%
% CREATED: [Date] MC
%
function rise_time = find_rise_time(trial)
    min_val = min(trial);
    max_val = max(trial);
    ten_percent = min_val + 0.1 * (max_val - min_val);
    ninety_percent = min_val + 0.9 * (max_val - min_val);
    
    % Find indices where the trial reaches 10% and 90% of the peak
    start_idx = find(trial >= ten_percent, 1);
    end_idx = find(trial >= ninety_percent, 1);
    
    if isempty(start_idx) || isempty(end_idx)
        rise_time = inf;  % If no valid rise time is found
    else
        rise_time = end_idx - start_idx;
    end
end