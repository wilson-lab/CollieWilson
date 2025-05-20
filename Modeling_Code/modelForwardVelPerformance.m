function modelForwardVelPerformance(predicted_RF)
%% modelForwardVelPerformance
% Plots the trajectories of a female fly and a modeled male fly (with and without forward scaling)
% for a range of forward velocities.
%
% INPUTS:
%   predicted_RF - Struct containing the visual receptive fields for the model.
%
% OUTPUTS:
%   None. The function generates a tiled plot of the first run for each forward velocity.

%% Initialize
[folder, ~, runSettings] = modelSettings();
close all;
runSettings.numRuns = 100; % Number of times to run the simulation (ideal > 20k)

runSettings.minIn = -4.25;      % Minimum input value for downstream process
runSettings.maxIn = 6;          % Maximum input value for downstream process
runSettings.alpha = 0.5;        % alpha, specifies slope of negative non-linear portion
runSettings.shift = -.25;       % shift, determines at what point linear/non-linear portions begin
runSettings.DNa02input = linspace(runSettings.minIn, runSettings.maxIn, 1000); % Input range for DNa02
runSettings.DNa02output = adjELU(runSettings.DNa02input, runSettings.alpha, runSettings.shift); % Output curve based on ELU nonlinearity

% Define the range to test
kValues = [6, 10]; % K values to iterate over
nK = length(kValues);
fwdVel = 5:1:20; % Forward velocities to iterate over
nFwd = length(fwdVel);
simDuration = 30; % Simulation duration in seconds
conditionColors = {"#0072BD", "#D95319"}; % Define colors for conditions
fwdScaleOptions = [0,30]; % Forward scaling factor
nScale = length(fwdScaleOptions);
noiseLevel = 5;

% Initialize storage for average displacement results for each k
avg_displacement_per_k = zeros(nK, 1);

% Initialize a cell array to store results for plotting
all_displacement_results = cell(nK, 1);

% Create a single figure for all k values
figure; set(gcf, 'Position', [100 100 1500 900]);  % Set figure size
tiledlayout(nK, 2, 'TileSpacing', 'compact'); % One tile for each k value
sgtitle('Trajectories at Forward Speed = 10 Across Different k Values'); % Title for the figure

% Loop over k values
for kIdx = 1:nK
    % Assign k value
    runSettings.k = kValues(kIdx);

    % Initialize storage for current k value
    avg_displacement_results = zeros(nFwd, nScale);

    % Loop over forward velocities
    for fIdx = 1:nFwd
        % Assign forward velocity
        thisFwd = fwdVel(fIdx);

        % Loop over scaling conditions
        for sIdx = 1:nScale
            % Assign scaling factor
            thisScale = fwdScaleOptions(sIdx);

            % Display current iteration details for debugging
            fprintf('Running simulation: k = %d, FwdVel = %d, FwdScale = %d\n', runSettings.k, thisFwd, thisScale);

            % Run the AOTU steering model
            [female, male] = aotu_steering_fwdmodel(predicted_RF, simDuration, thisFwd, thisScale, noiseLevel, runSettings);

            % Call calculate_displacement to get metrics
            [~, avg_avg_disp, ~] = calculate_displacement(male, female);

            % Store the total displacement
            avg_displacement_results(fIdx, sIdx) = avg_avg_disp;

            % Extract the first run for plotting
            female_x = squeeze(female(1, :, 1)); % Female X trajectory
            female_y = squeeze(female(1, :, 2)); % Female Y trajectory
            male_x = squeeze(male(1, :, 1));     % Male X trajectory
            male_y = squeeze(male(1, :, 2));     % Male Y trajectory

            % Plot only if forward velocity is 10
            if thisFwd == 5 && mod(kIdx,1)==0
                if sIdx == 1
                    nexttile;
                    hold on;
                    % Plot female and male (without scaling)
                    plot(female_x, female_y, 'k', 'LineWidth', 1, 'DisplayName', 'Female');
                    plot(male_x, male_y, 'Color', conditionColors{1}, 'LineWidth', 1, 'DisplayName', 'Male');
                    xlabel('X Position');
                    ylabel('Y Position');
                    title(sprintf('k = %d, FwdVel = %d', kValues(kIdx), thisFwd));
                    if fIdx == 1 && kIdx == 1
                        legend('show', 'Location', 'northeast');
                    end
                    grid on;
                else
                    % Plot male (with scaling)
                    plot(male_x, male_y, 'Color', conditionColors{2}, 'LineWidth', 1, 'DisplayName', 'Male with Forward');
                end
            elseif thisFwd == 15 && mod(kIdx,1)==0
                if sIdx == 1
                    nexttile;
                    hold on;
                    % Plot female and male (without scaling)
                    plot(female_x, female_y, 'k', 'LineWidth', 1, 'DisplayName', 'Female');
                    plot(male_x, male_y, 'Color', conditionColors{1}, 'LineWidth', 1, 'DisplayName', 'Male');
                    xlabel('X Position');
                    ylabel('Y Position');
                    title(sprintf('k = %d, FwdVel = %d', kValues(kIdx), thisFwd));
                    if fIdx == 1 && kIdx == 1
                        legend('show', 'Location', 'northeast');
                    end
                    grid on;
                else
                    % Plot male (with scaling)
                    plot(male_x, male_y, 'Color', conditionColors{2}, 'LineWidth', 1, 'DisplayName', 'Male with Forward');
                end
            end
        end
    end

    % Calculate the average displacement across all velocities and scaling conditions
    avg_displacement_per_k(kIdx) = mean(avg_displacement_results(:));

    % Store results for plotting
    all_displacement_results{kIdx} = avg_displacement_results;
end

sgtitle('Example Runs for forward velocity model');
saveas(gcf, fullfile(folder.final, ['forwardvelocity_ExampleRuns' '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, ['forwardvelocity_ExampleRuns' '.svg']));

%% Plot the displacement results
figure; set(gcf, 'Position', [100 100 1200 600]);  % Set figure size
tiledlayout(1,nK,"TileSpacing","compact")

for kIdx = 1:nK
    nexttile;
    results = all_displacement_results{kIdx};
    hold on;
    % Plot results for each scaling condition
    plot(fwdVel, results(:, 1), '-', 'Color', conditionColors{1}, 'DisplayName', 'Scale = 0');
    plot(fwdVel, results(:, 2), '-', 'Color', conditionColors{2}, 'DisplayName', ['Scale = ' num2str(fwdScaleOptions(2))]);
    xlabel('Forward velocity');
    ylabel('Average displacement');
    title(sprintf('k = %d', kValues(kIdx)));
    legend('show', 'Location', 'best');
    ylim([0 15])
    grid on;
end
saveas(gcf, fullfile(folder.final, ['forwardvelocity_displacement' '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, ['forwardvelocity_displacement' '.svg']));

end
