% analyzeSpikingAndVoltageDifference
%
% Calculates the average difference in firing rate and median-filtered
% voltage between P1 on and off conditions for all data and quiescent-only data.
%
% INPUTS:
%   spikerate - 3D array of firing rate data with dimensions [time x trials x condition]
%               where condition 1 corresponds to P1 on and condition 2 to P1 off.
%   voltage   - 3D array of voltage data with dimensions [time x trials x condition],
%               structured similarly to spikerate.
%   ttime     - Time vector corresponding to rows of spikerate and voltage data.
%   forward   - 2D array of forward movement data with dimensions [time x trials],
%               used to determine when the fly was stationary or moving.
%   settings  - Structure containing the run threshold setting (runThreshE).
%
% OUTPUTS:
%   diffSR - 2-element vector with average firing rate difference:
%            diffSR(1) for all data and diffSR(2) for quiescent-only data.
%   diffVm - 2-element vector with average voltage difference:
%            diffVm(1) for all data and diffVm(2) for quiescent-only data.
%
% CREATED: 11/09/2024 - MC
%
function [diffSR, diffVm] = analyzeSpikingAndVoltageDifference(spikerate, voltage, ttime, forward, settings)
    % Apply median filter to voltage data for both conditions
    voltageFilteredP1 = spikeFilter(voltage(:,:,1), ttime);
    voltageFilteredP2 = spikeFilter(voltage(:,:,2), ttime);

    % Initialize arrays to store differences for all data and quiescent-only data
    allDiffSpikeRate = zeros(1, size(spikerate, 2));
    allDiffVoltage = zeros(1, size(voltage, 2));
    quiescentDiffSpikeRate = zeros(1, size(spikerate, 2));
    quiescentDiffVoltage = zeros(1, size(voltage, 2));

    % Loop through each trial
    for trial = 1:size(spikerate, 2)
        % Determine if the animal was stationary or moving (quiescent when runIdx == 0)
        runIdx = schmittTrigger(forward(:, trial), settings.runThreshE, 0.1);

        % Calculate average firing rate for all data (P1 on and P2 off)
        avgSpikeRateP1_all = mean(spikerate(:, trial, 1), 'omitnan');
        avgSpikeRateP2_all = mean(spikerate(:, trial, 2), 'omitnan');
        
        % Calculate average voltage for all data (P1 on and P2 off) (filtered)
        avgVoltageP1_all = mean(voltageFilteredP1(:, trial), 'omitnan');
        avgVoltageP2_all = mean(voltageFilteredP2(:, trial), 'omitnan');

        % Calculate average firing rate for quiescent-only data
        avgSpikeRateP1_quiescent = mean(spikerate(~runIdx, trial, 1), 'omitnan');
        avgSpikeRateP2_quiescent = mean(spikerate(~runIdx, trial, 2), 'omitnan');
        
        % Calculate average voltage for quiescent-only data (filtered)
        avgVoltageP1_quiescent = mean(voltageFilteredP1(~runIdx, trial), 'omitnan');
        avgVoltageP2_quiescent = mean(voltageFilteredP2(~runIdx, trial), 'omitnan');

        % Compute differences for this trial (P1 on - P2 off)
        allDiffSpikeRate(trial) = avgSpikeRateP1_all - avgSpikeRateP2_all;
        allDiffVoltage(trial) = avgVoltageP1_all - avgVoltageP2_all;
        quiescentDiffSpikeRate(trial) = avgSpikeRateP1_quiescent - avgSpikeRateP2_quiescent;
        quiescentDiffVoltage(trial) = avgVoltageP1_quiescent - avgVoltageP2_quiescent;
    end

    % Calculate the average difference across trials for all data and quiescent-only data
    diffSR = [mean(allDiffSpikeRate, 'omitnan'), mean(quiescentDiffSpikeRate, 'omitnan')];
    diffVm = [mean(allDiffVoltage, 'omitnan'), mean(quiescentDiffVoltage, 'omitnan')];
end
