function [avg_total_displacement, avg_avg_displacement, displacement_distribution] = calculate_displacement(male_path, female_path)
    % Calculate the displacement between male and female trajectories.
    %
    % INPUTS:
    %   male_path    - [runs x time x 2] Array for male trajectories (x and y positions)
    %   female_path  - [runs x time x 2] Array for female trajectories (x and y positions)
    %
    % OUTPUTS:
    %   avg_total_displacement   - Average total displacement across all runs
    %   avg_avg_displacement     - Average of the average displacement across all runs
    %   displacement_distribution - Distribution of displacements within a range of 0-20
    
    % Validate inputs
    if size(male_path, 1) ~= size(female_path, 1) || ...
       size(male_path, 2) ~= size(female_path, 2) || ...
       size(male_path, 3) ~= 2 || size(female_path, 3) ~= 2
        error('Input arrays must have matching dimensions and a z-dimension of size 2 for spatial coordinates.');
    end

    % Get the number of runs
    num_runs = size(male_path, 1);
    num_time = size(male_path, 2);

    % Initialize displacement arrays
    total_displacement = zeros(num_runs, 1); % Total displacement for each run
    avg_displacement = zeros(num_runs, 1);   % Average displacement for each run
    all_displacements = [];                  % Collect displacements within 0-20

    % Loop through each run
    for runIdx = 1:num_runs
        % Extract male and female paths for the current run
        male_positions = squeeze(male_path(runIdx, :, :));   % [time x 2]
        female_positions = squeeze(female_path(runIdx, :, :)); % [time x 2]

        % Calculate displacements at each time point
        displacements = sqrt(sum((male_positions - female_positions).^2, 2)); % [time x 1]

        % Filter displacements within the range of 0-20
        valid_displacements = displacements(displacements <= 20);

        % Append to the distribution array
        all_displacements = [all_displacements; valid_displacements]; %#ok<AGROW>

        % Calculate total and average displacement for the current run
        total_displacement(runIdx) = sum(displacements);
        avg_displacement(runIdx) = mean(displacements);
    end

    % Compute averages across runs
    avg_total_displacement = mean(total_displacement); % Average total displacement across runs
    avg_avg_displacement = mean(avg_displacement);     % Average of the average displacements across runs

    % Compute the distribution of displacements within the range 0-20
    displacement_distribution = histcounts(all_displacements, 0:0.1:20, 'Normalization', 'probability');
end
