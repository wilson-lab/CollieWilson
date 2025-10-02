function binned_hist = input_histogram(input_history, bins)
% CREATED: 07/30/2025 - MC

% Set max value
maxBin = 1.8;

% Combine across trials and time for each neuron
[nTrials, nTime, nNeurons] = size(input_history);
binned_hist = zeros(bins, nNeurons);
edges = linspace(0, maxBin, bins+1);  % bin edges from 0 to max

for z = 1:nNeurons
    % Flatten all values for this neuron
    data = reshape(input_history(:,:,z), [], 1);
    
    % Compute histogram counts
    counts = histcounts(data, edges);
    
    % Normalize to sum to 1
    binned_hist(:,z) = counts;
end
