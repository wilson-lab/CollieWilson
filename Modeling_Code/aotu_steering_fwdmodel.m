% aotu_steering_fwdmodel
%
% Simulates the steering behavior of a modeled male fly based on visual and neural inputs.
% The model computes the trajectories of a female and male fly, incorporating visual 
% object position, rotational velocity, and AOTU019 and AOTU025 neural activity. The male's
% behavior is determined by its response to the female's position relative to its current heading.
%
% INPUTS:
%   noiTuning        - Struct containing visual tuning curves for different neurons (AOTU019, AOTU025, sum).
%                      Specifies how these neurons respond to object positions in the azimuthal space.
%   simDuration      - Duration of the simulation in seconds.
%   forwardVelocity  - Forward velocity of the fly (units per second).
%   fwdScale         - Scale factor for normalizing forward velocity.
%   noiseLevel       - Scale factor for adjusting trajectory noise of
%   female
%   runSettings      - Struct containing additional model parameters.
%
% OUTPUTS:
%   female           - [numRuns x nTime x 2] Trajectory of the female fly, containing x and y positions.
%   male             - [numRuns x nTime x 2] Trajectory of the male fly, containing x and y positions.
%
% CREATED: 11/17/24 - MC
%
function [female, male] = aotu_steering_fwdmodel(noiTuning, simDuration, forwardVelocity, fwdScale, noiseLevel, runSettings)
%% Initialize
% Extract model parameters from runSettings
visobj_position = runSettings.visObjPosition; % Object position (azimuthal space)

numRuns = runSettings.numRuns;   % Number of simulation runs
fs = runSettings.fs;             % Simulation update rate (Hz)

% Compute the number of time steps based on the simulation duration and fs
nTime = round(simDuration * fs); % Total number of time steps based on duration and update rate

%% Generate neural parameters
% Generate delays for each neuron in the model
delay_AOTU019 = round(runSettings.AOTU019_delay / 1000 * fs); % Delay for AOTU019 in samples
delay_Others = round(runSettings.Others_delay / 1000 * fs);   % Delay for AOTU025 and noiSum in samples

% Generate input-output for each neuron
% Visual tuning curves for each neuron type over 360 deg azimuthal space
AOTU019R = abs(noiTuning.AOTU019);  % Right side for AOTU019
AOTU019L = flip(AOTU019R);     % Left side (flipped)
AOTU025R = abs(noiTuning.AOTU025);  % Right side for AOTU025
AOTU025L = flip(AOTU025R);     % Left side (flipped)
groupedR = abs(noiTuning.sum);       % Right side for other neurons
groupedL = flip(groupedR);       % Left side (flipped)

AOTU019input = linspace(0, max(AOTU019R), 1000); % Input range
AOTU019output = linspace(0, 1, 1000); % Output range

AOTU025input = linspace(0, max(AOTU025R), 1000); % Input range
AOTU025output = linspace(0, 1, 1000); % Output range

% Motor input-output for DNa02
DNa02input = runSettings.DNa02input;  % Input range for downstream neurons
DNa02output = runSettings.DNa02output; % Output range for downstream neurons

%% Run the simulation
% Initialize parameters
femaleForwardVelocity = forwardVelocity;

if fwdScale==0
    forwardActivity = 1;
else
    forwardActivity = forwardVelocity/fwdScale;
end

k = runSettings.k; % Steering gain
dt = 1/fs; %change in time

% Generate random noise for the female's motion
rng(13); % Set random seed for reproducibility
inputTrajectory = randn(numRuns, nTime); % Random Gaussian noise
inputTrajectory = noiseLevel * zscore(inputTrajectory); % Scale noise amplitude
inputTrajectory = lowpass(inputTrajectory', 0.1, fs)'; % Apply low-pass filter to smooth noise

% Initialize angular trajectories
female_angle = zeros(numRuns, nTime); % Angular trajectory for the female
male_angle = zeros(numRuns, nTime);   % Angular trajectory for the male

% Initialize 2D positions
female_x = zeros(numRuns, nTime); % X position for the female
female_y = zeros(numRuns, nTime); % Y position for the female
male_x = zeros(numRuns, nTime);   % X position for the male
male_y = zeros(numRuns, nTime);   % Y position for the male

% Set initial positions
female_x(:, 1) = 0; 
female_y(:, 1) = 0; % Female starts at origin
male_x(:, 1) = 0; % Male starts to the left of the female
male_y(:, 1) = 0; % Male starts below the female

% Calculate initial direction of the female relative to the male
delta_x_initial = female_x(:, 1) - male_x(:, 1); % Relative X distance at t=1
delta_y_initial = female_y(:, 1) - male_y(:, 1); % Relative Y distance at t=1
initial_direction = atan2d(delta_y_initial, delta_x_initial); % Direction to the female in global coordinates

% Set male's initial heading directly toward the female
male_angle(:, 1) = wrapTo180(initial_direction); % Initialize male's angle

% Initialize neural history and outputs
input_history = zeros(numRuns, nTime, 2); % Inputs to downstream neurons (right/left)
DNa02R_history = zeros(numRuns, nTime);   % DNa02 right-side activity
DNa02L_history = zeros(numRuns, nTime);   % DNa02 left-side activity
DNa02RLdiff_history = zeros(numRuns, nTime); % Difference in right and left DNa02 activity
rotvel_history = zeros(numRuns, nTime);  % Male's rotational velocity

% Simulation loop
for t = 2:nTime
    % Update the female's angular trajectory by adding noise
    female_angle(:, t) = wrapTo180(female_angle(:, t-1) + inputTrajectory(:, t)); % Add noise and wrap angle

    % Compute the forward velocity components based on the female's heading
    fwd_x_female = femaleForwardVelocity * dt .* cosd(female_angle(:, t)); % Forward velocity in the x direction
    fwd_y_female = femaleForwardVelocity * dt .* sind(female_angle(:, t)); % Forward velocity in the y direction

    % Update the female's 2D position based on transformed forward velocity
    female_x(:, t) = female_x(:, t-1) + fwd_x_female;
    female_y(:, t) = female_y(:, t-1) + fwd_y_female;

    % Calculate the relative angular offset of the female from the male
    delta_x = female_x(:, t) - male_x(:, t-1); % Relative X distance
    delta_y = female_y(:, t) - male_y(:, t-1); % Relative Y distance
    direction_to_female = atan2d(delta_y, delta_x); % Direction to the female in global coordinates
    angular_offset = wrapTo180(direction_to_female - male_angle(:, t-1)); % Offset in male's coordinate frame

    % Clamp the angular offset to the male's visual receptive field
    angular_offset = max(min(angular_offset, 120), -120);

    % Store the angular offset for neural computations
    visobj_history(:, t) = angular_offset;

    % Apply delays for AOTU019 and AOTU025 neurons
    t_AOTU019 = max(t - delay_AOTU019, 1); % Delay for AOTU019
    t_Others = max(t - delay_Others, 1);   % Delay for AOTU025

    % Map delayed angular offsets to visual positions
    [~, p1_AOTU019] = ismember(round(visobj_history(:, t_AOTU019)), visobj_position);
    [~, p2_Others] = ismember(round(visobj_history(:, t_Others)), visobj_position);

    % Compute neural inputs
    r025_input = AOTU025R(p2_Others);
    l025_input = AOTU025L(p2_Others);
    r019_input = AOTU019R(p1_AOTU019)*forwardActivity;
    l019_input = AOTU019L(p1_AOTU019)*forwardActivity;

    % Handle cases where AOTU019 input exceeds available activity
    if any(r019_input>max(AOTU019input))
        r019_input(r019_input>max(AOTU019input)) = max(AOTU019input); % Set to the index of the maximum value
    end
    if any(l019_input>max(AOTU019input))
        l019_input(l019_input>max(AOTU019input)) = max(AOTU019input); % Set to the index of the maximum value
    end

    % Compute neural outputs for AOTU neurons
    r025o = AOTU025output(interp1(AOTU025input, 1:length(AOTU025input), r025_input, 'nearest'))';
    l025o = AOTU025output(interp1(AOTU025input, 1:length(AOTU025input), l025_input, 'nearest'))';
    r019o = AOTU019output(interp1(AOTU019input, 1:length(AOTU019input), r019_input, 'nearest'))';
    l019o = AOTU019output(interp1(AOTU019input, 1:length(AOTU019input), l019_input, 'nearest'))';

    % Compute steering drive
    current_inputR = r025o - l019o + groupedR(p2_Others);
    current_inputL = l025o - r019o + groupedL(p2_Others);
    input_history(:, t, 1) = current_inputR; % Right input
    input_history(:, t, 2) = current_inputL; % Left input

    % Compute DNa02 activity
    dr = interp1(DNa02input, 1:length(DNa02input), current_inputR, 'nearest');
    dl = interp1(DNa02input, 1:length(DNa02input), current_inputL, 'nearest');

    % Handle NaN cases by setting them to the maximum value in DNa02output
    if any(isnan(dr))
        dr(isnan(dr)) = length(DNa02output); % Set to the index of the maximum value
    end
    if any(isnan(dl))
        dl(isnan(dl)) = length(DNa02output); % Set to the index of the maximum value
    end

    DNa02R_history(:, t) = DNa02output(dr); % Right-side DNa02 activity
    DNa02L_history(:, t) = DNa02output(dl); % Left-side DNa02 activity

    % Calculate rotational velocity (difference in DNa02 activity scaled by gain)
    DNa02RLdiff_history(:, t) = DNa02R_history(:, t) - DNa02L_history(:, t);
    rotvel_history(:, t) = k * DNa02RLdiff_history(:, t);

    % Update the male's angular trajectory
    male_angle(:, t) = wrapTo180(male_angle(:, t-1) + rotvel_history(:, t));

    % Compute the forward velocity components based on the male's heading
    fwd_x = forwardVelocity * dt .* cosd(male_angle(:, t)); % Forward velocity in the x direction
    fwd_y = forwardVelocity * dt .* sind(male_angle(:, t)); % Forward velocity in the y direction

    % Compute the rotational velocity components (optional, if you want to add them as well)
    rot_x = rotvel_history(:, t) * dt .* cosd(male_angle(:, t)); % Rotational velocity in the x direction
    rot_y = rotvel_history(:, t) * dt .* sind(male_angle(:, t)); % Rotational velocity in the y direction

    % Update the male's 2D position by combining forward and rotational velocity components
    male_x(:, t) = male_x(:, t-1) + fwd_x + rot_x;
    male_y(:, t) = male_y(:, t-1) + fwd_y + rot_y;

end

% Combine the male and female positions into a single output array
female = cat(3, female_x, female_y); % Female trajectory as [x, y]
male = cat(3, male_x, male_y);       % Male trajectory as [x, y]

end