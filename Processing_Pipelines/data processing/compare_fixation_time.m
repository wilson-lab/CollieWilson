% compare_fixation_time
%
% This function calculates the percentage of time a fly spent fixating a target 
% with and without P1 activation. It uses fixation data across multiple trials 
% and, if available, compares it to pre-experiment acclimatization data.
%
% INPUTS:
%   fix_idx   - Logical matrix indicating fixation status (1 = fixating, 0 = not)
%               across multiple trials (each column represents a trial)
%   filename  - Base filename used to locate the acclimation file
%   folder    - Structure containing folder paths, specifically folder.accl for
%               locating acclimatization data
%
% OUTPUTS:
%   fixation_percentage_P1   - Percentage of time spent fixating with P1 activation
%   fixation_percentage_noP1 - Percentage of time spent fixating without P1 activation
%
% CREATED: 11/10/2024 - MC
%
function [fixation_percentage_P1, fixation_percentage_noP1] = compare_fixation_time(fix_idx, filename, folder)
    %% Initialize outputs
    fixation_percentage_P1 = NaN;
    fixation_percentage_noP1 = NaN;

    %% Change to the acclimation folder
    cd(folder.accl);
    % Generate the acclimation file name by replacing '_int' with '_accl'
    accl_filename = strrep(filename, '_int', '_acc');

    %% Check if acclimatization data is available
    if exist(accl_filename, 'file')
        % Load acclimatization data
        load(accl_filename);

        % Define the last 5 minutes in seconds
        last_5_min = 5 * 60;  % 5 minutes in seconds

        % Find the start index for the last 5 minutes
        end_index = length(int_accl_time);
        start_index = find(int_accl_time >= int_accl_time(end_index) - last_5_min, 1, 'first');

        % Pull the last 5 minutes of data
        accl_panelPs_5min = int_accl_panelPs(start_index:end_index);
        accl_forwardVelocity_5min = int_accl_forwardVelocity(start_index:end_index);

        % Run fixationFinder on the last 5 minutes of acclimatization data
        thisFixation = fixationFinder(accl_panelPs_5min, accl_forwardVelocity_5min, int_accl_time(start_index:end_index), 0);

        % Calculate fixation percentage with P1 activation across all trials
        fixation_time_P1 = sum(fix_idx(:));  % Sum across all trials
        total_points_P1 = numel(fix_idx);
        fixation_percentage_P1 = (fixation_time_P1 / total_points_P1) * 100;

        % Calculate fixation percentage without P1 activation across all trials
        fixation_time_noP1 = sum(thisFixation.idx_run(:));  % Sum across all points
        total_points_noP1 = numel(thisFixation.idx_run);
        fixation_percentage_noP1 = (fixation_time_noP1 / total_points_noP1) * 100;
    end
end