% spikeFilter
% This function applies a median filter to a voltage trace to remove high-frequency 
% components caused by spikes, smoothing the signal over a specified time window.
%
% INPUTS:
% voltage - the raw voltage trace (vector)
% ttime   - time vector (in seconds) corresponding to the voltage trace
%
% OUTPUT:
% voltage_mf - median filtered voltage trace
%
% CREATED:
% 08/02/2024 - MC
%
function [voltage_mf] = spikeFilter(voltage,ttime)
% set filter parameters
mfSize = 0.015; %s, increase to increase smoothing
[~,mfBin] = min(abs(ttime - mfSize)); %find nearest time index

% median filter voltage signal to remove spikes
voltage_mf = medfilt1(voltage,mfBin,'truncate');
end