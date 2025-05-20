% spikert_velocitymdl
% This function fits and evaluates multiple linear models using cross-validation 
% to determine the relationship between spike rate and two behavioral velocities 
% (forward and angular). It calculates R-squared values for models with added 
% terms and interaction terms between the velocities.
%
% INPUT
% spikert - array of spike rates (time x trials)
% forward - array of forward velocities (time x trials)
% angular - array of angular velocities (time x trials)
% ttime   - time vector (in seconds) for lag shifting the velocities
%
% OUTPUT
% R2_out - 3-element array containing mean R-squared values for the following:
%          [1] R-squared for added X1 and X2 (linear terms)
%          [2] R-squared for interaction X1*X2
%          [3] R-squared for added X1, X2 and interaction term
%
% CREATED: [Date] MC
%
function R2_out = spikert_velocitymdl(spikert, forward, angular, ttime)
%% Initialize

% Fetch processing settings
settings = processSettings;

% Adjust behavior for lag estimate
forward_lag = lagShift(forward, ttime, settings.fwdLag);
angular_lag = lagShift(angular, ttime, settings.angLag);

% Reshape data
firing_rate = reshape(spikert, [], 1);
forward_velocity = reshape(forward_lag, [], 1);
angular_velocity = reshape(angular_lag, [], 1);

% Combine variables into a table
T = table(forward_velocity, angular_velocity, firing_rate);

%% Cross-validation setup

% Define predictors and response
X1 = T.forward_velocity;  % Forward velocity (linear term)
X2 = T.angular_velocity;  % Angular velocity (linear term)
y = T.firing_rate;  % Firing rate

% Find rows with no NaNs in the predictors
valid_idx = ~any(isnan([X1, X2]), 2);

% Filter the data to only include valid rows
X1_valid = X1(valid_idx);
X2_valid = X2(valid_idx);
y_valid = y(valid_idx);

% Set the number of folds for cross-validation (e.g., 5)
numFolds = 5;

% Generate the cross-validation indices
cv = cvpartition(length(X1_valid), 'KFold', numFolds);

% Initialize arrays to store performance metrics
R2_added = zeros(numFolds,1);               % R-squared for added X1 and X2
R2_interaction = zeros(numFolds,1);          % R-squared for interaction X1*X2
R2_added_and_interaction = zeros(numFolds,1); % R-squared for added + interaction

% Define the models
linear_model = @(b, X) b(1) + b(2) * X(:,1) + b(3) * X(:,2);  % Linear model for added X1 and X2
interaction_model = @(b, X) b(1) + b(2) * X(:,1) .* X(:,2);   % Interaction-only model
added_and_interaction_model = @(b, X) b(1) + b(2) * X(:,1) + b(3) * X(:,2) + b(4) * X(:,1) .* X(:,2); % Added + interaction

% Loop through each fold
parfor i = 1:numFolds
    % Get the training and testing indices
    trainIdx = training(cv, i);
    testIdx = test(cv, i);
    
    % Split the data into training and testing sets (valid rows only)
    X1_train = X1_valid(trainIdx);
    X2_train = X2_valid(trainIdx);
    y_train = y_valid(trainIdx);
    
    X1_test = X1_valid(testIdx);
    X2_test = X2_valid(testIdx);
    y_test = y_valid(testIdx);
    
    % Initial guesses for model parameters
    initial_guess_linear = [0, 1, 1];  % 3 parameters for linear model
    initial_guess_interaction = [0, 1]; % 2 parameters for interaction-only model
    initial_guess_added_and_interaction = [0, 1, 1, 1]; % 4 parameters for added + interaction model
    
    % Fit the models
    % Added X1 and X2 (linear model)
    mdl_added = fitnlm([X1_train, X2_train], y_train, linear_model, initial_guess_linear);
    
    % Interaction X1*X2
    mdl_interaction = fitnlm([X1_train, X2_train], y_train, interaction_model, initial_guess_interaction);
    
    % Added + Interaction (linear + interaction model)
    mdl_added_and_interaction = fitnlm([X1_train, X2_train], y_train, added_and_interaction_model, initial_guess_added_and_interaction);
    
    % Predict on the test set using all models
    y_pred_added = predict(mdl_added, [X1_test, X2_test]);
    y_pred_interaction = predict(mdl_interaction, [X1_test, X2_test]);
    y_pred_added_and_interaction = predict(mdl_added_and_interaction, [X1_test, X2_test]);
    
    % Calculate R-squared for all models
    R2_added(i) = calculate_r2(y_test, y_pred_added);
    R2_interaction(i) = calculate_r2(y_test, y_pred_interaction);
    R2_added_and_interaction(i) = calculate_r2(y_test, y_pred_added_and_interaction);
end

% Store fold performance for output
R2_out(1) = mean(R2_added);  % Mean R-squared for added X1 and X2
R2_out(2) = mean(R2_interaction);  % Mean R-squared for interaction X1*X2
R2_out(3) = mean(R2_added_and_interaction);  % Mean R-squared for added + interaction

end
