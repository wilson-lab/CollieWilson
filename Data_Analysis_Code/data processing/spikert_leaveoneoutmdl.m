% spikert_leaveoneoutmdl
% This function fits a combined Gaussian and linear model to estimate the contribution of 
% forward velocity, angular velocity, visual object position, and object direction to the 
% firing rate. It uses a leave-one-out approach to assess the contribution of each variable 
% by training models that omit one predictor at a time and performing cross-validation 
% to compare their performance.
%
% INPUT
% spikert   - array of spike rates (time x trials)
% panelps   - array of visual object positions (time x trials)
% forward   - array of forward velocities (time x trials)
% angular   - array of angular velocities (time x trials)
% ttime     - time vector (in seconds) for lag shifting the velocities
%
% OUTPUT
% R2_fold   - performance (R-squared) across folds for models with and without each variable
% mdl_full  - the full Gaussian and linear model trained on the entire dataset
%
% CREATED: [Date] MC
%
function [R2_fold, mdl_full] = spikert_leaveoneoutmdl(spikert, panelps, forward, angular, ttime)
%% Initialize

% Fetch processing settings
settings = processSettings;

% Adjust behavior for lag estimate
forward_lag = lagShift(forward, ttime, settings.fwdLag);
angular_lag = lagShift(angular, ttime, settings.angLag);

% Determine object direction and visual object position for each pulse
paneldir = fetchPanelDir(panelps);

% Reshape data
firing_rate = reshape(spikert, [], 1);
forward_velocity = reshape(forward_lag, [], 1);
angular_velocity = reshape(angular_lag, [], 1);
visual_object_position = reshape(panelps, [], 1);
object_direction = reshape(paneldir, [], 1);  % Reintroducing object direction

% Combine variables in a single data table
T = table(visual_object_position, object_direction, forward_velocity, angular_velocity, firing_rate);

%% Test model against data folds

% Define predictors and response
X1 = T.forward_velocity;  % Forward velocity (linear term)
X2 = T.angular_velocity;  % Angular velocity (linear term)
X3 = T.visual_object_position;  % Visual object position (Gaussian term)
X4 = T.object_direction;  % Object direction (linear term)
y = T.firing_rate;  % Firing rate

% Find rows with no NaNs in the predictors
valid_idx = ~any(isnan([X1, X2, X3, X4]), 2);

% Filter the data to only include valid rows
X1_valid = X1(valid_idx);
X2_valid = X2(valid_idx);
X3_valid = X3(valid_idx);
X4_valid = X4(valid_idx);
y_valid = y(valid_idx);

% Set the number of folds for cross-validation (e.g., 5)
numFolds = 5;

% Generate the cross-validation indices
cv = cvpartition(length(X1_valid), 'KFold', numFolds);

% Initialize arrays to store performance metrics
R2_full_model = zeros(numFolds,1);        % R-squared for full model
R2_wo_X1 = zeros(numFolds,1);             % R-squared without forward velocity (X1)
R2_wo_X2 = zeros(numFolds,1);             % R-squared without angular velocity (X2)
R2_wo_X3 = zeros(numFolds,1);             % R-squared without visual object position (X3)
R2_wo_X4 = zeros(numFolds,1);             % R-squared without object direction (X4)

% Define the combined Gaussian + Linear model with object direction
gaussian_linear_model = @(b, X) b(1) + b(2) * X(:,1) + b(3) * X(:,2) + ...
                               b(4) * exp(-(X(:,3) - b(5)).^2 / (2 * b(6)^2)) + b(7) * X(:,4);

% Loop through each fold
parfor i = 1:numFolds
    % Get the training and testing indices
    trainIdx = training(cv, i);
    testIdx = test(cv, i);
    
    % Split the data into training and testing sets (valid rows only)
    X1_train = X1_valid(trainIdx);
    X2_train = X2_valid(trainIdx);
    X3_train = X3_valid(trainIdx);
    X4_train = X4_valid(trainIdx);
    y_train = y_valid(trainIdx);
    
    X1_test = X1_valid(testIdx);
    X2_test = X2_valid(testIdx);
    X3_test = X3_valid(testIdx);
    X4_test = X4_valid(testIdx);
    y_test = y_valid(testIdx);

    % Initial guesses for the combined model parameters
    initial_guess = [0, 1, 1, 1, median(X3_train), iqr(X3_train), 1];
    
    % Fit the combined model (Gaussian for X3, linear for X1, X2, and X4)
    mdl_full = fitnlm([X1_train, X2_train, X3_train, X4_train], y_train, gaussian_linear_model, initial_guess);
    
    % Predict on the test set using the full model
    y_pred_full = predict(mdl_full, [X1_test, X2_test, X3_test, X4_test]);
    R2_full_model(i) = calculate_r2(y_test, y_pred_full);
    
    % Without Forward Velocity (X1)
    mdl_wo_X1 = fitnlm([X2_train, X3_train, X4_train], y_train, ...
                       @(b, X) b(1) + b(2) * X(:,1) + b(3) * exp(-(X(:,2) - b(4)).^2 / (2 * b(5)^2)) + b(6) * X(:,3), ...
                       initial_guess([1, 3:end]));
    y_pred_wo_X1 = predict(mdl_wo_X1, [X2_test, X3_test, X4_test]);
    R2_wo_X1(i) = calculate_r2(y_test, y_pred_wo_X1);
    
    % Without Angular Velocity (X2)
    mdl_wo_X2 = fitnlm([X1_train, X3_train, X4_train], y_train, ...
                       @(b, X) b(1) + b(2) * X(:,1) + b(3) * exp(-(X(:,2) - b(4)).^2 / (2 * b(5)^2)) + b(6) * X(:,3), ...
                       initial_guess([1, 2, 4:end]));
    y_pred_wo_X2 = predict(mdl_wo_X2, [X1_test, X3_test, X4_test]);
    R2_wo_X2(i) = calculate_r2(y_test, y_pred_wo_X2);
    
    % Without Visual Object Position (X3)
    mdl_wo_X3 = fitlm([X1_train, X2_train, X4_train], y_train);  % No Gaussian term for X3
    y_pred_wo_X3 = predict(mdl_wo_X3, [X1_test, X2_test, X4_test]);
    R2_wo_X3(i) = calculate_r2(y_test, y_pred_wo_X3);
    
    % Without Object Direction (X4)
    mdl_wo_X4 = fitnlm([X1_train, X2_train, X3_train], y_train, ...
                       @(b, X) b(1) + b(2) * X(:,1) + b(3) * X(:,2) + b(4) * exp(-(X(:,3) - b(5)).^2 / (2 * b(6)^2)), ...
                       initial_guess([1, 2, 3, 5:end]));
    y_pred_wo_X4 = predict(mdl_wo_X4, [X1_test, X2_test, X3_test]);
    R2_wo_X4(i) = calculate_r2(y_test, y_pred_wo_X4);
end

% Store fold performance at each step for output
R2_fold(1) = mean(R2_full_model); % full model
R2_fold(2) = mean(R2_wo_X1);      % w/o forward velocity
R2_fold(3) = mean(R2_wo_X2);      % w/o angular velocity
R2_fold(4) = mean(R2_wo_X3);      % w/o visual object position
R2_fold(5) = mean(R2_wo_X4);      % w/o object direction


%% Retrain on full dataset

% Initial guesses for the combined model parameters
initial_guess = [0, 1, 1, 1, median(X3_valid), iqr(X3_valid), 1];

% Fit the combined Gaussian + Linear model on the full dataset
mdl_full = fitnlm([X1_valid, X2_valid, X3_valid, X4_valid], y_valid, gaussian_linear_model, initial_guess);

% Optionally plot predictions based on the full model (if enabled)
if 0
    % Define a range of visual object positions for plotting
    position_range = linspace(min(X3_valid), max(X3_valid), 1000);  % Increased resolution with 1000 points

    % Use fixed mean velocities for forward and angular velocity (example values)
    fixed_X1 = mean(X1_valid);  % Forward velocity
    fixed_X2 = mean(X2_valid);  % Angular velocity
    fixed_X4 = mean(X4_valid);  % Object direction (fixed or can vary)

    % Create input data for prediction
    predict_input = [repmat(fixed_X1, length(position_range), 1), ...  % Repeat fixed forward velocity
                     repmat(fixed_X2, length(position_range), 1), ...  % Repeat fixed angular velocity
                     position_range', ...                              % Vary visual object position
                     repmat(fixed_X4, length(position_range), 1)];     % Repeat fixed object direction

    % Predict using the Gaussian + Linear model trained on the full dataset
    full_predictions = predict(mdl_full, predict_input);

    % Plot the results
    figure;
    plot(position_range, full_predictions, 'r-', 'LineWidth', 2);
    xlabel('Visual Object Position');
    ylabel('Predicted Firing Rate');
    title('Combined Gaussian + Linear Model Predictions vs Visual Object Position');
    grid on;
end
end
