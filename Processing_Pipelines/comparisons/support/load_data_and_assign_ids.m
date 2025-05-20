% LOAD_DATA_AND_ASSIGN_IDS
%
% This function loads cross-correlation data for each animal, assigns unique IDs,
% and extracts peak lag and correlation values for both angular (ang) and forward (fwd) motion.
% Data is extracted for two conditions: background (dark) and motion.
%
% INPUTS:
%   folder_motion      - Path to the folder containing motion condition data.
%   folder_background  - Path to the folder containing background (dark) condition data.
%   all_names          - Cell array of animal names.
%   settings           - Structure containing settings, including:
%                        settings.minXCorrProm - Minimum prominence for peak detection.
%
% OUTPUTS:
%   animal_ids         - Array of unique animal IDs.
%   conditions         - Array indicating condition type (1 = Dark, 2 = Motion).
%   rval_data          - Nx2 array of peak correlation values [ang, fwd].
%   lag_data           - Nx2 array of peak lag times [ang, fwd].
%
% CREATED: 03/18/2025 - MC

function [animal_ids, conditions, rval_data, lag_data] = load_data_and_assign_ids(folder_motion, folder_background, all_names, settings)
    % Load cross-correlation data, assigning animal IDs across conditions
    animal_ids = [];
    conditions = [];
    rval_data = [];
    lag_data = [];
    
    % Determine the year prefix based on folder name
    if contains(folder_background, 'AOTU019')
        year_prefix = '2023';
    elseif contains(folder_background, 'AOTU025')
        year_prefix = '2024';
    else
        error('Folder name must contain either "AOTU019" or "AOTU025" to set the correct year prefix.');
    end
    
    unique_id = 1; % Initialize unique ID counter
    animal_id_map = containers.Map; % Map to hold animal name-ID pairs
    
    for n = 1:length(all_names)
        name = all_names{n};
        
        % Assign an ID to each animal, reusing if animal is present in both conditions
        if ~isKey(animal_id_map, name)
            animal_id_map(name) = unique_id;
            unique_id = unique_id + 1;
        end
        current_id = animal_id_map(name);
        
        % Load background (dark) pulse data if available
        backgroundFile = fullfile(folder_background, [year_prefix '_' name '_1_xc.mat']);
        if isfile(backgroundFile)
            load(backgroundFile, 'r_val', 'lag_t');
            [peak_lag_background, peak_rval_background, r_val] = find_peak_lag_rval(r_val, lag_t, settings.minXCorrProm);
            rval_data = [rval_data; peak_rval_background.ang, peak_rval_background.fwd];
            lag_data = [lag_data; peak_lag_background.ang, peak_lag_background.fwd];
            animal_ids = [animal_ids; current_id];
            conditions = [conditions; 1]; % 1 = Dark condition
        end
        
        % Load motion pulse data if available
        motionFile = fullfile(folder_motion, [name '_xc.mat']);
        if isfile(motionFile)
            load(motionFile, 'r_val', 'lag_t');
            [peak_lag_motion, peak_rval_motion, r_val] = find_peak_lag_rval(r_val, lag_t, settings.minXCorrProm);
            rval_data = [rval_data; peak_rval_motion.ang, peak_rval_motion.fwd];
            lag_data = [lag_data; peak_lag_motion.ang, peak_lag_motion.fwd];
            animal_ids = [animal_ids; current_id];
            conditions = [conditions; 2]; % 2 = Motion condition
        end
    end
end
