function [bins, norm_counts] = compute_firingrate_distribution(int_spikert)
% CREATED: 10/20/2025 - MC
%
% Compute normalized firing rate distribution across all trials.
% Data are organized as time (x) by trials (y).

% Flatten across time and trials
all_spikes = int_spikert(:);

% Define bins
bins = 10:5:150;

% Compute histogram counts
counts = histcounts(all_spikes, bins);

% Normalize to probability (sum = 1)
norm_counts = counts / sum(counts, 'omitnan');

end
