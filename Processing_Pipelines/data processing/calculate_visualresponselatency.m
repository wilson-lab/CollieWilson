% CALCULATE_VISUALRESPONSELATENCY - This function calculates the response latency 
% between a visual motion pulse and a neuron firing response. The latency is defined
% as the time at which the first peak in the derivative of the firing rate occurs
% following the onset of the stimulus.

% INPUTS:
%   pulse_srR   - 3D array of firing rate data over time (rows) for each animal (columns)
%                 at each visual pulse position (3rd dimension).
%   pulse_posR  - 3D array of visual object positions over time (rows) at each visual 
%                 pulse position (3rd dimension). Columns represent the same positions 
%                 across all animals.
%   time_pulse  - Time vector (in ms) corresponding to the temporal resolution of the data.

% OUTPUT:
%   latency_times - Array of response latencies (in ms) for each animal. Each value
%                   represents the time relative to stimulus onset when the steepest 
%                   change in firing rate occurred.
%
% Created: 12/18/2024 by MC
%
function latency_times = calculate_visualresponselatency(pulse_srR, pulse_posR, time_pulse)
%% Initialize
% Fetch number of animals
nFlies = size(pulse_srR, 2);

% Find the x range for NaNs in pulse_posR for the first sweep position
nan_indices = isnan(pulse_posR(:, 1, 1));
start_motion = find(~nan_indices, 1, 'first'); % First non-NaN (start of motion)
end_motion = find(~nan_indices, 1, 'last');   % Last non-NaN (end of motion)

% Determine the duration of the motion pulse
pulse_duration = end_motion - start_motion + 1;
selected_time = time_pulse(1:pulse_duration);

% Set data outside the motion range to []
pulse_srR_visible = pulse_srR;
pulse_srR_visible(end_motion+1:end, :, :) = [];
pulse_srR_visible(1:start_motion-1, :, :) = [];

% Calculate average firing rate across animals for each pulse position
avg_firing_rate = squeeze(mean(pulse_srR_visible, 2, 'omitnan'));

% Find the pulse position with the highest average firing rate across animals
[~, max_pulse_idx] = max(mean(avg_firing_rate, 1, 'omitnan'));
% Extract firing rate data for the pulse with the highest average firing rate
selected_firing_rate = squeeze(pulse_srR_visible(:, :, max_pulse_idx));

%% Analyze firing latency

% Initialize output array
latency_times = nan(1, nFlies);

% Calculate the first derivative and find the peak latency for each animal
for animal_idx = 1:nFlies
    % Extract firing rate for the current animal
    firing_rate = selected_firing_rate(:, animal_idx);
    
    % Calculate the first derivative (rate of change)
    dFiringRate = diff(firing_rate);
    
    % Find the index of the first peak in the derivative
    [~, peak_idx] = max(dFiringRate);
    
    % Record the peak time relative to stimulus onset
    latency_times(animal_idx) = selected_time(peak_idx);
end

%% Filter out latencies > 150 ms
latency_threshold = 150; % Define the threshold in ms
valid_indices = latency_times <= latency_threshold; % Find valid points
latency_times = latency_times(valid_indices); % Retain only valid latencies

end