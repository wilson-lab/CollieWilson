% aotu_steering_model
%
% Simulates the steering behavior of a fly based on AOTU019 and AOTU025 neural input without feedforward control.
% The model generates visual object position, rotational velocity, and inputs from various neurons, 
% and computes the corresponding output based on the AOTU019 synapse type (inhibitory or excitatory).
%
% INPUTS:
%   noiTuning        - Struct containing visual tuning curves for different neurons (AOTU019, AOTU025, sum).
%                      It specifies how these neurons respond to object positions in the azimuthal space.
%   noiseLevel       - Noise level added to the rotational velocity signal.
%   startPos         - Initial object position (azimuthal angle) in degrees.
%   AOTU019synapse   - Type of synapse for AOTU019 (e.g., 'inhibitory', 'excitatory').
%   simDuration      - Duration of the simulation in seconds.
%   k                - Gain factor for computing rotational velocity from neural inputs.
%   runSettings      - Struct containing additional model parameters.
%
% OUTPUTS:
%   timebase         - Time vector for the simulation, from 0 to simDuration, based on the simulation rate (fs).
%   visobj_history   - History of the visual object positions over time for each simulation run.
%   input_history    - History of input values to the DNa02 neurons for each simulation run, for right and left sides.
%   rotvel_history   - History of the fly's rotational velocity (in degrees per second) for each simulation run.
%
% CREATED: 10/30/24 - MC
% UPDATED: 11/06/24 - MC added nonlinearity to AOTU019 and AOTU025
%          11/16/24 - MC added direction selectivity
%
function [timebase, visobj_history, input_history, rotvel_history] = aotu_steering_model(...
    noiTuning, noiseLevel, startPos, AOTU019synapse, simDuration, k, startTime, runSettings)

%% Initialize
% Extract model parameters from runSettings
visobj_position = runSettings.visObjPosition; % Object position (azimuthal space)

numRuns = runSettings.numRuns;   % Number of simulation runs
fs = runSettings.fs;             % Simulation update rate (Hz)
fpass = runSettings.fpass;       % Lowpass filter cutoff (Hz)

% Compute the number of time steps based on the simulation duration and fs
nTime = round(simDuration * fs); % Total number of time steps based on duration and update rate
jumptime=20;
jump_interval = jumptime * fs;

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

% Option - include silencing
if max(AOTU019R)==0
    AOTU019input = linspace(0, 1, 1000); % Input range
else
    AOTU019input = linspace(0, max(AOTU019R), 1000); % Input range
end
AOTU019output = ELU(AOTU019input); % Output range
if max(AOTU025R)==0
    AOTU025input = linspace(0, 1, 1000); % Input range
else
    AOTU025input = linspace(0, max(AOTU025R), 1000); % Input range
end
AOTU025output = ELU(AOTU025input); % Output range

% Optional - include direction selectivity
if runSettings.dirselective == 1
    dsi019 = runSettings.AOTU019dsi_penalty;
    dsi025 = runSettings.AOTU025dsi_penalty;
elseif runSettings.dirselective == 2
    dsi019 = runSettings.AOTU025dsi_penalty;
    dsi025 = runSettings.AOTU019dsi_penalty;
else
    dsi019 = 1;
    dsi025 = 1;
end
% Motor input-output for DNa02
DNa02input = runSettings.DNa02input;  % Input range for downstream neurons
DNa02output = runSettings.DNa02output; % Output range for downstream neurons

%% Run the simulation
%Initialize simulation outputs
input_history = zeros(numRuns, nTime, 6); % Inputs to each left/right neuron
DNa02R_history = zeros(numRuns, nTime);   % DNa02 right side activity
DNa02L_history = zeros(numRuns, nTime);   % DNa02 left side activity
DNa02RLdiff_history = zeros(numRuns, nTime);   % DNa02 right-left difference activity
rotvel_history = zeros(numRuns, nTime);   % Fly's rotational velocity over time

% Random noise component for rotational velocity (modeling the influence of noise)
rng(13); % Set seed for reproducibility
inputTrajectory = randn(numRuns, nTime);       % Random noise drawn from Gaussian distribution
inputTrajectory = lowpass(inputTrajectory', fpass, fs)'; % Low-pass filter the noise
inputTrajectory = noiseLevel * zscore(inputTrajectory);  % Scale the noise

% Initialize visual object position and history
startSide = ones(numRuns, 1); 
startSide(2:2:end) = -1;  % Alternate which side the object starts on
visobj_history = zeros(numRuns, nTime);  % Visual object position history
visobj_history(:, 1:startTime) = repmat((startPos * startSide'), startTime, 1)'; % Set initial position

% For each time step
for t = startTime+1:nTime
    % Compute delayed indices for AOTU019, AOTU025, and noiSum
    t_AOTU019 = max(t - delay_AOTU019, 1); % Apply delay for AOTU019
    t_Others = max(t - delay_Others, 1);   % Apply delay for AOTU025 and noiSum

    % Get previous object positions at the delayed times
    [~, p1_AOTU019] = ismember(wrapTo180(round(visobj_history(:, t_AOTU019))), visobj_position);
    [~, p2_Others] = ismember(wrapTo180(round(visobj_history(:, t_Others))), visobj_position);

    % Determine object motion direction based on delayed positions
    % Compute rightward and leftward motion for AOTU019
    if t_AOTU019 - 1 < 1
        right_motion_019 = zeros(size(visobj_history(:, t_AOTU019))); % Set motion to 0
        left_motion_019 = zeros(size(visobj_history(:, t_AOTU019)));  % Set motion to 0
    else
        right_motion_019 = visobj_history(:, t_AOTU019) - visobj_history(:, t_AOTU019 - 1);
        left_motion_019 = -right_motion_019; % Left motion is the negative of right motion
    end

    % Compute rightward and leftward motion for AOTU025
    if t_Others - 1 < 1
        right_motion_025 = zeros(size(visobj_history(:, t_Others))); % Set motion to 0
        left_motion_025 = zeros(size(visobj_history(:, t_Others)));  % Set motion to 0
    else
        right_motion_025 = visobj_history(:, t_Others) - visobj_history(:, t_Others - 1);
        left_motion_025 = -right_motion_025; % Left motion is the negative of right motion
    end

    % Initialize direction selectivity multipliers as ones
    ds_multiplier_right_019 = ones(size(right_motion_019));
    ds_multiplier_left_019 = ones(size(left_motion_019));
    ds_multiplier_right_025 = ones(size(right_motion_025));
    ds_multiplier_left_025 = ones(size(left_motion_025));

    % Apply direction selectivity for AOTU019
    ds_multiplier_right_019(right_motion_019 < 0) = dsi019; % Moving left affects right
    ds_multiplier_left_019(left_motion_019 < 0) = dsi019;  % Moving right affects left

    % Apply direction selectivity for AOTU025
    ds_multiplier_right_025(right_motion_025 < 0) = dsi025; % Moving left affects right
    ds_multiplier_left_025(left_motion_025 < 0) = dsi025;  % Moving right affects left

    % Apply direction selectivity to inputs for AOTU025
    r025_input = AOTU025R(p2_Others) .* ds_multiplier_right_025; % Adjust for rightward motion
    l025_input = AOTU025L(p2_Others) .* ds_multiplier_left_025;  % Adjust for leftward motion

    % Apply direction selectivity to inputs for AOTU019
    r019_input = AOTU019R(p1_AOTU019) .* ds_multiplier_right_019; % Adjust for rightward motion
    l019_input = AOTU019L(p1_AOTU019) .* ds_multiplier_left_019;  % Adjust for leftward motion

    % Compute input-output for AOTU025
    r025i = interp1(AOTU025input, 1:length(AOTU025input), r025_input, 'nearest');
    l025i = interp1(AOTU025input, 1:length(AOTU025input), l025_input, 'nearest');
    r025o = AOTU025output(r025i)'; % Generate 025 output
    l025o = AOTU025output(l025i)';

    % Compute input-output for AOTU019
    r019i = interp1(AOTU019input, 1:length(AOTU019input), r019_input, 'nearest');
    l019i = interp1(AOTU019input, 1:length(AOTU019input), l019_input, 'nearest');
    r019o = AOTU019output(r019i)'; % Generate 019 output
    l019o = AOTU019output(l019i)';

    % Compute steering drive based on AOTU019 synapse type
    switch AOTU019synapse
        case "inhibitory"
            % AOTU019 inhibitory and contralateral
            current_inputR = r025o - l019o + groupedR(p2_Others);
            current_inputL = l025o - r019o + groupedL(p2_Others);
        case "excitatory"
            % AOTU019 excitatory and ipsilateral
            current_inputR = r025o + r019o + groupedR(p2_Others);
            current_inputL = l025o + l019o + groupedL(p2_Others);
    end

    % Store the input history
    input_history(:, t, 1) = r019o;
    input_history(:, t, 2) = l019o;
    input_history(:, t, 3) = r025o;
    input_history(:, t, 4) = l025o;
    input_history(:, t, 5) = groupedR(p2_Others);
    input_history(:, t, 6) = groupedL(p2_Others);

    %% Apply delay for DNa02 output
    % Compute DNa02 activity without delay
    dr = interp1(DNa02input, 1:length(DNa02input), current_inputR, 'nearest');
    DNa02R_history(:, t) = DNa02output(dr); % Right side DNa02 activity at time t
    dl = interp1(DNa02input, 1:length(DNa02input), current_inputL, 'nearest');
    DNa02L_history(:, t) = DNa02output(dl); % Left side DNa02 activity at time t

    % Calculate difference in DNa02 activity
    DNa02RLdiff_history(:,t) = DNa02R_history(:, t) - DNa02L_history(:, t);
    % Compute the rotational velocity at time t based on scaled DNa02 activity
    rotvel_history(:, t) = k * DNa02RLdiff_history(:,t);

    % Update the visual object position with the delayed rotational velocity influence
    visobj_history(:, t) = visobj_history(:, t) + visobj_history(:, t-1) - rotvel_history(:, t);

    % Add noise to the delayed visual object trajectory
    visobj_history(:, t) = visobj_history(:, t) + inputTrajectory(:, t);

    % Ensure that the visual object position stays within the range [-180, 180]
    visobj_history(:, t) = wrapTo180(visobj_history(:, t));

    % Random jump every jumptime seconds
    if mod(t, round(jumptime * fs)) == 0
        % Randomly add or subtract 100° to each run
        jumpDirection = randi([0 1], numRuns, 1) * 2 - 1; % Returns -1 or +1
        visobj_history(:, t) = visobj_history(:, t) + 100 * jumpDirection;
        visobj_history(:, t) = wrapTo180(visobj_history(:, t)); % Keep within [-180, 180]
    end

end

% Generate timebase for the simulation based on actual sampling rate
timebase = linspace(0, simDuration, nTime); % Simulate from 0 to total simulation duration

% Calculate time step based on sampling rate (fs)
dt = 2 / fs;
rotvel_history = rotvel_history / dt;

