% spikert_objmdl
% This function fits and evaluates multiple models relating spike rate to
% forward and angular velocities, visual object position, and object direction
% using cross-validation. It computes R-squared values for different combinations 
% of predictors, including Gaussian modeling of object position.
%
% INPUT
% spikert   - array of spike rates (time x trials)
% panelps   - array of visual object positions (time x trials)
% forward   - array of forward velocities (time x trials)
% angular   - array of angular velocities (time x trials)
% ttime     - time vector (in seconds) for lag shifting the velocities
%
% OUTPUT
% R2_obj - 8-element array containing mean R-squared values for the following:
%          [1] R-squared for Gaussian object position only
%          [2] R-squared for forward + Gaussian object position
%          [3] R-squared for forward * Gaussian object position (interaction)
%          [4] R-squared for angular + Gaussian object position
%          [5] R-squared for angular * Gaussian object position (interaction)
%          [6] R-squared for forward + angular + Gaussian object position
%          [7] R-squared for Gaussian object position + object direction
%          [8] R-squared for Gaussian object position * object direction (interaction)
%
% CREATED: [Date] MC
%
function R2_obj = spikert_objmdl(spikert, panelps, forward, angular, ttime)
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
visual_object_position = reshape(panelps, [], 1);
object_direction = reshape(fetchPanelDir(panelps), [], 1);  % Object direction

% Combine variables into a table
T = table(forward_velocity, angular_velocity, visual_object_position, object_direction, firing_rate);

%% Cross-validation setup

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
R2_gaussian_only = zeros(numFolds,1);               % R-squared for Gaussian object position only
R2_fwd_gauss = zeros(numFolds,1);                   % R-squared for forward + Gaussian object position
R2_fwdTimesGauss = zeros(numFolds,1);               % R-squared for forward * Gaussian object position
R2_ang_gauss = zeros(numFolds,1);                   % R-squared for angular + Gaussian object position
R2_angTimesGauss = zeros(numFolds,1);               % R-squared for angular * Gaussian object position
R2_fwd_ang_gauss = zeros(numFolds,1);               % R-squared for forward + angular + Gaussian object position
R2_gauss_obj_dir = zeros(numFolds,1);               % R-squared for Gaussian object position + object direction
R2_gaussTimes_obj_dir = zeros(numFolds,1);          % R-squared for Gaussian object position * object direction

% Define the Gaussian model (consistent across all models)
gaussian_model = @(b, X) b(1) + b(2) * exp(-(X(:,1) - b(3)).^2 / (2 * b(4)^2));  % Gaussian for visual object position only
fwd_gauss_model = @(b, X) b(1) + b(2) * X(:,1) + b(3) * exp(-(X(:,2) - b(4)).^2 / (2 * b(5)^2));  % Forward + Gaussian
fwdTimesGauss_model = @(b, X) b(1) + b(2) * (X(:,1) .* exp(-(X(:,2) - b(3)).^2 / (2 * b(4)^2)));  % Forward * Gaussian
ang_gauss_model = @(b, X) b(1) + b(2) * X(:,1) + b(3) * exp(-(X(:,2) - b(4)).^2 / (2 * b(5)^2));  % Angular + Gaussian
angTimesGauss_model = @(b, X) b(1) + b(2) * (X(:,1) .* exp(-(X(:,2) - b(3)).^2 / (2 * b(4)^2)));  % Angular * Gaussian
fwd_ang_gauss_model = @(b, X) b(1) + b(2) * X(:,1) + b(3) * X(:,2) + b(4) * exp(-(X(:,3) - b(5)).^2 / (2 * b(6)^2));  % Forward + Angular + Gaussian
gauss_obj_dir_model = @(b, X) b(1) + b(2) * exp(-(X(:,1) - b(3)).^2 / (2 * b(4)^2)) + b(5) * X(:,2);  % Gaussian object position + object direction
gaussTimes_obj_dir_model = @(b, X) b(1) + b(2) * (exp(-(X(:,1) - b(3)).^2 / (2 * b(4)^2)) .* X(:,2));  % Gaussian object position * object direction

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
    
    % Initial guesses for model parameters
    initial_guess_gaussian = [0, 1, median(X3_train), iqr(X3_train)];  % 4 parameters for Gaussian model
    initial_guess_fwd_gauss = [0, 1, 1, median(X3_train), iqr(X3_train)];  % 5 parameters for forward + Gaussian
    initial_guess_timesGauss = [0, 1, median(X3_train), iqr(X3_train)];  % 4 parameters for interaction models
    initial_guess_fwd_ang_gauss = [0, 1, 1, 1, median(X3_train), iqr(X3_train)];  % 6 parameters for forward + angular + Gaussian
    initial_guess_gauss_obj_dir = [0, 1, median(X3_train), iqr(X3_train), 1];  % 5 parameters for Gaussian + object direction
    initial_guess_gaussTimes_obj_dir = [0, 1, median(X3_train), iqr(X3_train)];  % 4 parameters for Gaussian * object direction
    
    % Fit the models
    % Gaussian object position only
    mdl_gaussian_only = fitnlm(X3_train, y_train, gaussian_model, initial_guess_gaussian);
    
    % Forward + Gaussian object position
    mdl_fwd_gauss = fitnlm([X1_train, X3_train], y_train, fwd_gauss_model, initial_guess_fwd_gauss);
    
    % Forward * Gaussian object position (interaction model)
    mdl_fwdTimesGauss = fitnlm([X1_train, X3_train], y_train, fwdTimesGauss_model, initial_guess_timesGauss);
    
    % Angular + Gaussian object position
    mdl_ang_gauss = fitnlm([X2_train, X3_train], y_train, ang_gauss_model, initial_guess_fwd_gauss);
    
    % Angular * Gaussian object position (interaction model)
    mdl_angTimesGauss = fitnlm([X2_train, X3_train], y_train, angTimesGauss_model, initial_guess_timesGauss);
    
    % Forward + Angular + Gaussian object position
    mdl_fwd_ang_gauss = fitnlm([X1_train, X2_train, X3_train], y_train, fwd_ang_gauss_model, initial_guess_fwd_ang_gauss);
    
    % Gaussian object position + object direction
    mdl_gauss_obj_dir = fitnlm([X3_train, X4_train], y_train, gauss_obj_dir_model, initial_guess_gauss_obj_dir);
    
    % Gaussian object position * object direction (interaction model)
    mdl_gaussTimes_obj_dir = fitnlm([X3_train, X4_train], y_train, gaussTimes_obj_dir_model, initial_guess_gaussTimes_obj_dir);
    
    % Predict on the test set using all models
    y_pred_gaussian_only = predict(mdl_gaussian_only, X3_test);
    y_pred_fwd_gauss = predict(mdl_fwd_gauss, [X1_test, X3_test]);
    y_pred_fwdTimesGauss = predict(mdl_fwdTimesGauss, [X1_test, X3_test]);
    y_pred_ang_gauss = predict(mdl_ang_gauss, [X2_test, X3_test]);
    y_pred_angTimesGauss = predict(mdl_angTimesGauss, [X2_test, X3_test]);
    y_pred_fwd_ang_gauss = predict(mdl_fwd_ang_gauss, [X1_test, X2_test, X3_test]);
    y_pred_gauss_obj_dir = predict(mdl_gauss_obj_dir, [X3_test, X4_test]);
    y_pred_gaussTimes_obj_dir = predict(mdl_gaussTimes_obj_dir, [X3_test, X4_test]);
    
    % Calculate R-squared for all models
    R2_gaussian_only(i) = calculate_r2(y_test, y_pred_gaussian_only);
    R2_fwd_gauss(i) = calculate_r2(y_test, y_pred_fwd_gauss);
    R2_fwdTimesGauss(i) = calculate_r2(y_test, y_pred_fwdTimesGauss);
    R2_ang_gauss(i) = calculate_r2(y_test, y_pred_ang_gauss);
    R2_angTimesGauss(i) = calculate_r2(y_test, y_pred_angTimesGauss);
    R2_fwd_ang_gauss(i) = calculate_r2(y_test, y_pred_fwd_ang_gauss);
    R2_gauss_obj_dir(i) = calculate_r2(y_test, y_pred_gauss_obj_dir);
    R2_gaussTimes_obj_dir(i) = calculate_r2(y_test, y_pred_gaussTimes_obj_dir);
end

% Store fold performance for output
R2_obj(1) = mean(R2_gaussian_only);  % Mean R-squared for Gaussian object position only
R2_obj(2) = mean(R2_fwd_gauss);  % Mean R-squared for forward + Gaussian object position
R2_obj(3) = mean(R2_fwdTimesGauss);  % Mean R-squared for forward * Gaussian object position
R2_obj(4) = mean(R2_ang_gauss);  % Mean R-squared for angular + Gaussian object position
R2_obj(5) = mean(R2_angTimesGauss);  % Mean R-squared for angular * Gaussian object position
R2_obj(6) = mean(R2_fwd_ang_gauss);  % Mean R-squared for forward + angular + Gaussian object position
R2_obj(7) = mean(R2_gauss_obj_dir);  % Mean R-squared for Gaussian object position + object direction
R2_obj(8) = mean(R2_gaussTimes_obj_dir);  % Mean R-squared for Gaussian object position * object direction

end
