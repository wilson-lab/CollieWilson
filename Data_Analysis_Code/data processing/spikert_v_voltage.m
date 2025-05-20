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
%   - bin : Array of voltage bin centers.
%   - fr  : Array of average firing rates corresponding to each voltage bin.
%
% UPDATED: 11/12/2024 - MC added bin min requirements
%
function binOut = spikert_v_voltage(spikert, voltage)
%% Initialize
% Fetch settings (any relevant settings used by this function are retrieved here)
settings = processSettings();

% Define binning parameters for voltage
vmMax = -30;   % Maximum voltage (mV)
vmMin = -80;   % Minimum voltage (mV)
vb = 2;        % Voltage bin width (mV)
minPoints = 0; % Minimum number of points required in each bin to calculate the mean

%% Reshape datasets
% Reshape spikert and voltage data into column vectors for easier indexing
spikert_r = reshape(spikert, [], 1);
voltage_r = reshape(voltage, [], 1);

%% Discretize voltage data
% Define bin edges and labels for voltage discretization
vm_edge = vmMin - vb / 2 : vb : vmMax + vb / 2; % Define edges of each bin
vm_bins = vmMin : vb : vmMax;                   % Define bin centers as labels
% Discretize voltage data, assigning each value to a bin based on bin edges
voltage_disc = discretize(voltage_r, vm_edge, vm_bins);

%% Calculate mean firing rate for each voltage bin
% Initialize array to store mean firing rates for each bin
spikert_bin = [];

% Loop through each voltage bin
for v = 1:length(vm_bins)
    thisBin = vm_bins(v);                    % Current bin center (voltage level)
    thisVmIdx = find(voltage_disc == thisBin); % Indices of voltage data in this bin
    
    % Check if bin contains at least minPoints before calculating mean
    if thisVmIdx > minPoints
        spikert_bin(v) = mean(spikert_r(thisVmIdx), 'omitnan'); % Calculate mean firing rate, ignoring NaNs
    else
        spikert_bin(v) = nan; % Assign NaN if bin has fewer than minPoints
    end
end

%% Output results
% Store bin centers and calculated firing rates in output structure
binOut.bin = vm_bins;    % Bin centers
binOut.fr = spikert_bin; % Firing rate per bin

end
