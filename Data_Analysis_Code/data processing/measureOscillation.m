% measureOscillation
% This function analyzes the oscillations in panel positions by identifying peaks 
% and troughs in the data. It calculates the amplitude of oscillations as the 
% difference between peaks and troughs and generates a plot showing the amplitude 
% over time.
%
% INPUTS:
%   panelps - Panel positions (1D array or vector) representing the oscillation data
%
% OUTPUTS:
%   osc_out - A structure containing the amplitudes of the oscillations (not currently used in the function)
%
% CREATED: [Date] MC
%
function osc_out = measureOscillation(panelps)
% Flatten the data if it's not a 1D array
panelps = panelps(:);

% Find peaks and troughs
[peaks, peak_locs] = findpeaks(panelps);
[troughs, trough_locs] = findpeaks(-panelps);

% Convert troughs back to positive values
troughs = -troughs;

% Ensure that peaks and troughs are paired correctly
min_length = min(length(peaks), length(troughs));
peaks = peaks(1:min_length);
troughs = troughs(1:min_length);
peak_locs = peak_locs(1:min_length);
trough_locs = trough_locs(1:min_length);

% Calculate amplitude over time (peak-to-trough)
amplitudes = peaks - troughs;

% Plot amplitude over time
figure;
plot(peak_locs, amplitudes, '-o');
xlabel('Time (index)');
ylabel('Oscillation Amplitude');
title('Oscillation Amplitude Over Time');
grid on;
end