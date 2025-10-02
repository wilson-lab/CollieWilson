function plot_binned_histograms(all_binned_hist, comparisonLabel)
% CREATED: 07/30/2025 - MC

% Check input size
[nBins, nNeurons, nConds] = size(all_binned_hist);
if nNeurons ~= 6
    error('Expected 6 neurons (3 L/R pairs).');
end

% Pair indices: [1 2], [3 4], [5 6]
pairIdx = [1 2; 3 4; 5 6];

% Omit first bin
omit_bin = 1;
valid_bins = (omit_bin+1):nBins;
nBins_used = numel(valid_bins);

% Define bin centers assuming linspace from 0 to 1.8
bin_centers = linspace(0, 1.8, nBins);
bin_centers = bin_centers(valid_bins);

% Average L/R pairs
binned_avg = zeros(nBins_used, 3, nConds);  % bins × neurons × conditions
for p = 1:3
    temp = mean(all_binned_hist(valid_bins, pairIdx(p,:), :), 2);  % omit bin 1
    binned_avg(:,p,:) = reshape(temp, [nBins_used, 1, nConds]);
end

% Normalize to total activity in first condition
for p = 1:3
    baseline_max = max(binned_avg(:,p,1));  % sum of condition 1 after omitting bin 1
    for c = 1:nConds
        binned_avg(:,p,c) = binned_avg(:,p,c) / baseline_max;
    end
end

% Plot
figure;
tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');

for p = 1:3
    nexttile;
    hold on;
    for c = 1:nConds
        plot(bin_centers, binned_avg(:,p,c), 'LineWidth', 1.5);
    end
    title(['Neuron ' num2str(p)]);
    xlabel('Normalized Activity');
    ylabel('Normalized Counts (Rel. to Condition 1)');
    xlim([0 1.8]);
    ylim([0 1.5])
    legend(comparisonLabel, 'Location', 'northeast');
end
