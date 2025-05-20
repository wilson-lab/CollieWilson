function [avg_total_offset, avg_avg_offset, angular_distribution] = calculate_angular_offset(male_path, female_path)
    % Calculate the angular offset of the female relative to the male's travel vector.
    %
    % INPUTS:
    %   male_path    - [runs x time x 2] Array for male trajectories (x and y positions)
    %   female_path  - [runs x time x 2] Array for female trajectories (x and y positions)
    %
    % OUTPUTS:
    %   avg_total_offset    - Average total angular offset across all runs
    %   avg_avg_offset      - Average of the average angular offsets across all runs
    %   angular_distribution - Angular offset distribution (in degrees) normalized to 0-360 degrees

    % Validate inputs
    if size(male_path, 1) ~= size(female_path, 1) || ...
       size(male_path, 2) ~= size(female_path, 2) || ...
       size(male_path, 3) ~= 2 || size(female_path, 3) ~= 2
        error('Input arrays must have matching dimensions and a z-dimension of size 2 for spatial coordinates.');
    end

    % Get the number of runs and time points
    num_runs = size(male_path, 1);
    num_time = size(male_path, 2);

    % Initialize angular offset arrays
    total_offset = zeros(num_runs, 1); % Total angular offset for each run
    avg_offset = zeros(num_runs, 1);   % Average angular offset for each run
    all_offsets = [];                  % Collect angular offsets across all runs

    % Loop through each run
    for runIdx = 1:num_runs
        % Extract male and female paths for the current run
        male_positions = squeeze(male_path(runIdx, :, :));   % [time x 2]
        female_positions = squeeze(female_path(runIdx, :, :)); % [time x 2]

        % Calculate the male's travel vector at each time point
        male_vectors = diff(male_positions, 1, 1); % [time-1 x 2]
        male_vectors = [male_vectors; male_vectors(end, :)]; % Duplicate last vector to match dimensions

        % Calculate vectors from male to female
        male_to_female_vectors = female_positions - male_positions; % [time x 2]

        % Normalize vectors
        male_unit_vectors = male_vectors ./ vecnorm(male_vectors, 2, 2);
        male_to_female_unit_vectors = male_to_female_vectors ./ vecnorm(male_to_female_vectors, 2, 2);

        % Calculate angular offsets (in degrees)
        dot_products = dot(male_unit_vectors, male_to_female_unit_vectors, 2); % [time x 1]
        angles = acosd(max(-1, min(1, dot_products))); % Clamp values to avoid numerical issues

        % Collect angular offsets for this run
        all_offsets = [all_offsets; angles]; %#ok<AGROW>

        % Calculate total and average angular offset for the current run
        total_offset(runIdx) = sum(angles, 'omitnan');
        avg_offset(runIdx) = mean(angles, 'omitnan');
    end

    % Compute averages across runs
    avg_total_offset = mean(total_offset); % Average total angular offset across runs
    avg_avg_offset = mean(avg_offset);     % Average of the average angular offsets across runs

    % Compute the distribution of angular offsets (0-360 degrees)
    angular_distribution = histcounts(all_offsets, 0:5:360, 'Normalization', 'probability');
end
