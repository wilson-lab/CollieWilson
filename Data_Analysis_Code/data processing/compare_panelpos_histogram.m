% compare_panelpos_histogram
%
% This function compares histograms of panel position during the last
% 5 minutes of acclimatization (no P1) and during fixation trials with P1 activation.
%
% INPUTS:
%   panelpos_fix - Panel position data (time x trials) during P1 activation
%   filename     - Base filename used to locate the acclimatization file
%   folder       - Structure containing folder paths, specifically folder.accl
%
% OUTPUTS:
%   hist_noP1    - Histogram of panel positions during last 5 min of acclimation
%   hist_P1      - Histogram of panel positions during P1 activation
%
% CREATED: 04/10/2025 - MC
%
function [hist_noP1, hist_P1] = compare_panelpos_histogram(panelpos_fix, filename, folder)

    % Initialize outputs
    hist_noP1 = [];
    hist_P1 = [];

    %% Load acclimatization data
    cd(folder.accl);
    accl_filename = strrep(filename, '_int', '_acc');

    if exist(accl_filename, 'file')
        load(accl_filename);

        % Last 5 minutes
        last_5_min = 5 * 60; % seconds
        end_index = length(int_accl_time);
        start_index = find(int_accl_time >= int_accl_time(end_index) - last_5_min, 1, 'first');

        % Get panel position data
        accl_panelPs_5min = int_accl_panelPs(start_index:end_index);

        % Generate histogram for noP1
        [hist_noP1, ~] = panel_histogram(accl_panelPs_5min, [], 1); % normalize

        % Generate histogram for P1
        [hist_P1, ~] = panel_histogram(panelpos_fix, [], 1); % normalize
    else
        warning('Acclimatization file not found: %s', accl_filename);
    end
end
