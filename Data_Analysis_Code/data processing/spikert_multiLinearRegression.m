% spikert_multiLinearRegression
% This function performs multi-linear regression to estimate how much each
% independent variable (visual object position, forward velocity, and angular velocity) 
% contributes to the dependent variable (firing rate). The model uses cross-validation 
% to evaluate performance across multiple folds.
%
% INPUT
% spikert   - array of spike rates (time x trials)
% panelps   - array of visual object positions (time x trials)
% forward   - array of forward velocities (time x trials)
% angular   - array of angular velocities (time x trials)
% ttime     - time vector (in seconds) for lag shifting the velocities
%
% OUTPUT
% model_out - structure containing model performance metrics (e.g., R-squared, MSE)
%
% CREATED 09/10/2024 - MC
%
function model_out = spikert_multiLinearRegression(spikert,panelps,forward,angular,ttime)
%% intiailize

% fetch processing settings
settings = processSettings;

% adjust behavior for lag estimate
forward_lag = lagShift(forward,ttime,settings.fwdLag);
angular_lag = lagShift(angular,ttime,settings.angLag);

% reshape
firing_rate = reshape(spikert,[],1);
visual_object_position = reshape(panelps,[],1);
forward_velocity = reshape(forward_lag,[],1);
angular_velocity = reshape(angular_lag,[],1);

% combine variables in a single data table
T = table(visual_object_position,forward_velocity,angular_velocity,firing_rate);
% define predictors and response
X = T{:,{'forward_velocity','angular_velocity','visual_object_position'}};
y = T.firing_rate;

% set number of folkds (e.g., 5)
numFolds = 5;
% generate the cross-validation indices
cv = cvpartition(size(X,1), 'KFold', numFolds);

%% fit data

% initialize
R2_values = zeros(numFolds,1);  % r-squared values for each fold
MSE_values = zeros(numFolds,1); % Mean Squared Error for each fold

% Loop through each fold
for i = 1:numFolds
    % Get the training and testing indices
    trainIdx = training(cv, i);  % Logical indices for training data
    testIdx = test(cv, i);       % Logical indices for testing data
    
    % Split the data into training and testing sets
    X_train = X(trainIdx, :);  % Training data (predictors)
    y_train = y(trainIdx);     % Training data (response)
    
    X_test = X(testIdx, :);    % Testing data (predictors)
    y_test = y(testIdx);       % Testing data (response)

    % For motion pulse experiments, object position is largely nans which
    % are incompatible with fitlm and predict
    % Train data with velocity predictors (all time points)
    mdl_velonly = fitlm(X_train(:,1:2),y_train,'interactions');
    % Train data with all predictors (only motion pulse time points)
    mdl_all = fitlm(X_train,y_train,'interactions');
    
    % fetch test data where predictors are available
    idx_velonly = ~any(isnan(X_test(:,1:2)),2);
    idx_all = ~any(isnan(X_test),2);
    % predict for the test data
    y_predict_velonly = predict(mdl_velonly,X_test(idx_velonly,1:2));
    y_predict_all = predict(mdl_all,X_test(idx_all,:));

    % Calculate performance metrics for velocity only data
    % R-squared
    SS_residuals = sum((y_test(idx_velonly) - y_predict_velonly).^2); % Sum of squared residuals
    SS_total = sum((y_test(idx_velonly) - mean(y_test(idx_velonly))).^2); % Total sum of squares
    R2_values(i,1) = 1 - (SS_residuals / SS_total);
    % Mean Squared Error (MSE)
    MSE_values(i,1) = mean((y_test(idx_velonly) - y_predict_velonly).^2,'omitnan');
    % Calculate performance metrics for all predictor data
    % R-squared
    SS_residuals = sum((y_test(idx_all) - y_predict_all).^2); % Sum of squared residuals
    SS_total = sum((y_test(idx_all) - mean(y_test(idx_all))).^2); % Total sum of squares
    R2_values(i,2) = 1 - (SS_residuals / SS_total);
    % Mean Squared Error (MSE)
    MSE_values(i,2) = mean((y_test(idx_all) - y_predict_all).^2,'omitnan');
    
end

% store for output



end