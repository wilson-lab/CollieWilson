% spikert_v_voltage
%
% Analyzes the relationship between spike rate (firing rate) and membrane 
% voltage by binning the voltage data and calculating the average firing 
% rate for each voltage bin.
%
% INPUTS:
% - spikert : Array of spike times or firing rates (e.g., spikes per second).
% - voltage : Array of corresponding membrane voltage readings (in mV).
%
% OUTPUT:
% - binOut : Structure containing:
%   - bin  : Array of voltage bin centers.
%   - fr   : Array of average firing rates corresponding to each voltage bin.
%   - prob : Probability distribution of voltage occurrences across bins.
%
% UPDATED: 10/14/2025 - MC added probability of voltage distribution
%
function binOut = spikert_v_voltage(spikert, voltage)
%% Initialize
settings = processSettings();

% Define binning parameters for voltage
vmMax = -30;   % Maximum voltage (mV)
vmMin = -80;   % Minimum voltage (mV)
vb = 2;        % Voltage bin width (mV)
minPoints = 0; % Minimum number of points required in each bin to calculate the mean

%% Reshape datasets
spikert_r = reshape(spikert, [], 1);
voltage_r = reshape(voltage, [], 1);

%% Discretize voltage data
vm_edge = vmMin - vb/2 : vb : vmMax + vb/2; % Edges
vm_bins = vmMin : vb : vmMax;               % Bin centers
voltage_disc = discretize(voltage_r, vm_edge, vm_bins);

%% Calculate mean firing rate per voltage bin
spikert_bin = nan(size(vm_bins));
volt_counts  = zeros(size(vm_bins));

for v = 1:length(vm_bins)
    thisVmIdx = find(voltage_disc == vm_bins(v));
    volt_counts(v) = numel(thisVmIdx); % Count samples in this bin
    
    if numel(thisVmIdx) > minPoints
        spikert_bin(v) = mean(spikert_r(thisVmIdx), 'omitnan');
    end
end

%% Compute probability of each voltage bin
prob_voltage = volt_counts / sum(volt_counts, 'omitnan');

%% Output
binOut.bin  = vm_bins;
binOut.fr   = spikert_bin;
binOut.prob = prob_voltage;

end
