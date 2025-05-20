% process_crossings
% This function processes the crossings of zero in angular velocity data, identifying when the fly's angular velocity
% changes direction within a specified time window around each provided crossing index.
%
% INPUTS:
%   cross_indices - Indices where initial crossings are detected (e.g., panel crossing or target position)
%   angular_data  - Angular velocity data (vector)
%   panelps_data  - Panel position data (vector)
%   idx_window    - Size of the time window (in indices) to search for zero crossings
%   ttime         - Time vector (seconds)
%
% OUTPUTS:
%   zero_cross_times - Array of times (in seconds) when the angular velocity crosses zero within the specified window
%
% CREATED: [Date] MC
%
function zero_cross_times = process_crossings(cross_indices, angular_data, panelps_data, idx_window, ttime)
    zero_cross_times = [];
    for i = 1:length(cross_indices)
        idx = cross_indices(i);
        if idx > idx_window && idx + idx_window <= length(angular_data)
            thiswin_angular = angular_data(idx - idx_window:idx + idx_window);
            thiswin_panelps = panelps_data(idx - idx_window:idx + idx_window);

            sign_change_window = diff(sign(thiswin_angular(idx_window + 1:end)));
            zero_cross_idx = find(sign_change_window > 0, 1, 'first');
            if ~isempty(zero_cross_idx)
                zero_cross_times(end+1) = ttime(idx_window + 1 + zero_cross_idx);
            end
        end
    end
end