% check_nans_per_bin
% This function calculates the percentage of NaN values for each bin in three 
% 3D arrays corresponding to different datasets (KIR, WT, NA). It identifies 
% bins where the mean percentage of NaNs exceeds 70% and removes those bins 
% from the datasets.
%
% INPUTS:
%   kirData - 3D array of KIR data (bins x conditions x trials)
%   wtData  - 3D array of WT data (bins x conditions x trials)
%   naData  - 3D array of NA data (bins x conditions x trials)
%
% OUTPUTS:
%   processed_kirData - KIR data with bins having high NaN percentages removed
%   processed_wtData  - WT data with bins having high NaN percentages removed
%   processed_naData  - NA data with bins having high NaN percentages removed
%
% CREATED: [Date] MC
%
function [processed_kirData, processed_wtData, processed_naData] = check_nans_per_bin(kirData, wtData, naData)
    % Function to calculate NaN percentage for each bin in a 3D array
    function nan_percentage_per_bin = calculate_nan_percentage(array)
        num_bins = size(array, 1);  % Number of bins (x-axis)
        nan_percentage_per_bin = zeros(num_bins, 1);  % Initialize array to store NaN percentages per bin
        for x = 1:num_bins
            total_elements = numel(array(x, :, :));  % Total number of elements in the bin
            nan_count = sum(isnan(array(x, :, :)), 'all');  % Count NaNs in the bin
            nan_percentage_per_bin(x) = (nan_count / total_elements) * 100;  % Calculate NaN percentage
        end
    end

    % Calculate NaN percentages for each dataset separately
    kir_nan_percentages = calculate_nan_percentage(kirData);
    wt_nan_percentages = calculate_nan_percentage(wtData);
    na_nan_percentages = calculate_nan_percentage(naData);
    
    % Take the mean of NaN percentages across the datasets
    all_nan_percentages = [kir_nan_percentages, wt_nan_percentages, na_nan_percentages];
    mean_nan_percentages = mean(all_nan_percentages, 2, 'omitnan');  % Mean across datasets (columns), ignoring NaNs

    % Find bins with mean NaN percentage > 70%
    bins_to_remove = mean_nan_percentages > 70;

    % Remove the identified bins from each dataset
    processed_kirData = kirData(~bins_to_remove, :, :);
    processed_wtData = wtData(~bins_to_remove, :, :);
    processed_naData = naData(~bins_to_remove, :, :);
    
end
