% processRFdata
% Processes the receptive field (RF) data from various neuron types 
% in the Pursuit_RFs dataset. Specifically filters and sums the RFs 
% of all neurons except AOTU019 and AOTU025, and prepares the final 
% predicted RF data.
%
% INPUTS:
%   Pursuit_RFs - Table containing RF data for different neuron types.
%   runSettings - Settings for the run, currently unused in this function.
%
% OUTPUT:
%   predicted_RF - Table containing the summed RF of all other neurons 
%                  (excluding AOTU019 and AOTU025) along with the RFs 
%                  of AOTU019 and AOTU025.
%
function predicted_RF = processRFdata(Pursuit_RFs)
    % Fetch all neuron cell types from the dataset
    Pursuit_RFs(:,1) = [];
    celltypes = Pursuit_RFs.Properties.VariableNames;

    % Filter out AOTU019 and AOTU025 (specific neurons of interest)
    minor_idx = ~ismember(celltypes, {'AOTU019','AOTU025'});
    minor_RFs = Pursuit_RFs(:, minor_idx); % RF data for all other neuron types

    % some neurons projects contralaterally, flip its RF data to match visual input alignment
    minor_contra_idx = ismember(minor_RFs.Properties.VariableNames, {'AOTU002_a','AOTU002_c'});
    minor_RFs(:, minor_contra_idx) = flip(minor_RFs(:, minor_contra_idx), 1); % Flip RF data for contralateral neurons

    % Sum RF data across all neurons except AOTU019 and AOTU025
    minor_RFs = sum(minor_RFs, 2);

    % (Optional) Apply a light Gaussian smoothing to the summed RF data
    % gwin = 15; % Smoothing window size
    % rfOther = smoothdata(rfOther, "gaussian", gwin);

    % Create the final predicted_RF table that includes:
    % 1. The summed RF of all other neurons (excluding AOTU019, AOTU025)
    % 2. The individual RFs of AOTU019 and AOTU025
    predicted_RF = [minor_RFs Pursuit_RFs(:, ismember(celltypes, {'AOTU019', 'AOTU025'}))];
end