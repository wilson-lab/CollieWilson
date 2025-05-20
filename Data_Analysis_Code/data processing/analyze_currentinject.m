% analyze_currentinject
% This function analyzes behavioral responses to repeated current injections. It calculates the 
% frequency of positive turns during depolarizing pulses and negative turns during hyperpolarizing 
% pulses while the animal is running forward. Additionally, it computes kernel density estimates 
% (KDEs) normalized to a max of 1 for the distribution of average turn magnitudes.
%
% INPUTS:
%   currentInj  - Current injection data (amperes) over time for each trial
%   angularVel  - Angular velocity data (degrees/second) over time for each trial
%   forwardVel  - Forward velocity data (mm/s) over time for each trial
%   settings    - Structure containing experimental parameters (e.g., run threshold)
%
% OUTPUTS:
%   depolarize_freq    - Frequency of positive turns during depolarizing pulses
%   hyperpolarize_freq  - Frequency of negative turns during hyperpolarizing pulses
%   depol_kde           - Kernel density estimate (KDE) for depolarizing turns, normalized to max=1
%   hyperpol_kde        - Kernel density estimate (KDE) for hyperpolarizing turns, normalized to max=1
%
% CREATED: 01/07/2025 - MC
% UPDATED: 01/10/2024 - MC added kde
%
function [depolarize_freq, hyperpolarize_freq, depol_kde, hyperpol_kde] = analyze_currentinject(currentInj, angularVel, forwardVel, settings)
    % Get indices where the animal was running forward
    runIdx = schmittTrigger(forwardVel, settings.runThreshE, 0.1);

    % Initialize variables
    numTrials = size(currentInj, 2);
    depolarize_turns = [];
    hyperpolarize_turns = [];
    turnMin = 15; % Minimum value to be considered a turn

    % Loop through each trial
    for trial = 1:numTrials
        % Fetch data for the current trial
        currentTrial = currentInj(:, trial);
        angularTrial = angularVel(:, trial);
        runningIndices = runIdx(:, trial);

        % Identify depolarize and hyperpolarize steps
        depolarizeIdx = find(currentTrial > 50);
        hyperpolarizeIdx = find(currentTrial < -50);

        % Analyze depolarizing steps
        for step = depolarizeIdx'
            if runningIndices(step) % Check if the animal is running
                % Calculate mean angular velocity during the depolarizing step
                depolarize_turns = [depolarize_turns; mean(angularTrial(step))];
            end
        end

        % Analyze hyperpolarizing steps
        for step = hyperpolarizeIdx'
            if runningIndices(step) % Check if the animal is running
                % Calculate mean angular velocity during the hyperpolarizing step
                hyperpolarize_turns = [hyperpolarize_turns; mean(angularTrial(step))];
            end
        end
    end

    % Calculate frequencies for turning
    depolarize_freq = sum(depolarize_turns > turnMin) / length(depolarize_turns);
    hyperpolarize_freq = sum(hyperpolarize_turns < -turnMin) / length(hyperpolarize_turns);

    % Generate kernel density estimates (KDEs) for turn magnitudes
    bins = 0:10:300; % Bin edges
    binCenters = bins(1:end-1) + diff(bins) / 2; % Calculate bin centers

    if ~isempty(depolarize_turns)
        [depol_density, depol_x] = ksdensity(depolarize_turns, binCenters, 'Function', 'pdf');
        depol_kde = depol_density / max(depol_density); % Normalize to a max of 1
    else
        depol_kde = zeros(size(binCenters)); % Empty case
    end

    if ~isempty(hyperpolarize_turns)
        [hyperpol_density, hyperpol_x] = ksdensity(abs(hyperpolarize_turns), binCenters, 'Function', 'pdf');
        hyperpol_kde = hyperpol_density / max(hyperpol_density); % Normalize to a max of 1
    else
        hyperpol_kde = zeros(size(binCenters)); % Empty case
    end

end
