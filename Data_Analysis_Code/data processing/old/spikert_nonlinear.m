function [R2_linear, R2_truncatedGaussian, mdl_linear, mdl_truncatedGaussian] = spikert_nonlinear(spikert, panelps, forward, angular, ttime, peak_position, spread)
    %% Initialize

    % Fetch processing settings
    settings = processSettings;

    % Adjust behavior for lag estimate
    forward_lag = lagShift(forward, ttime, settings.fwdLag);
    angular_lag = lagShift(angular, ttime, settings.angLag);

    % Replace NaNs in visual_object_position with a placeholder value
    placeholder_value = -999;
    panelps(isnan(panelps)) = placeholder_value;

    % Reshape variables
    firing_rate = reshape(spikert, [], 1);
    forward_velocity = reshape(forward_lag, [], 1);
    angular_velocity = reshape(angular_lag, [], 1);
    visual_object_position = reshape(panelps, [], 1);

    % Define valid indices (non-NaN) for linear and Gaussian models
    valid_idx_linear = ~any(isnan([forward_velocity, angular_velocity]), 2);
    valid_idx_gaussian = visual_object_position ~= placeholder_value;

    %% Linear Model (Forward and Angular Velocity Only)
    % Separate valid data for the linear model
    X_linear = [forward_velocity(valid_idx_linear), angular_velocity(valid_idx_linear)];
    y_linear = firing_rate(valid_idx_linear);

    % Cross-validation for linear model
    numFolds = 5;
    cv = cvpartition(length(y_linear), 'KFold', numFolds);
    R2_folds_linear = zeros(numFolds, 1);

    parfor i = 1:numFolds
        trainIdx = training(cv, i);
        testIdx = test(cv, i);

        % Fit linear model on forward and angular velocity
        mdl_linear = fitlm(X_linear(trainIdx, :), y_linear(trainIdx));
        y_pred = predict(mdl_linear, X_linear(testIdx, :));
        R2_folds_linear(i) = calculate_r2(y_linear(testIdx), y_pred);
    end

    % Mean R-squared across folds for the linear model
    R2_linear = mean(R2_folds_linear);

    % Fit the final linear model on the full dataset
    mdl_linear = fitlm(X_linear, y_linear);

    %% Truncated Gaussian Model (Visual Object Position Only)
    X_gaussian = visual_object_position(valid_idx_gaussian);
    y_gaussian = firing_rate(valid_idx_gaussian);

    % Calculate a fixed baseline as the mean of spikert
    baseline_estimate = mean(spikert(:), 'omitnan');

    % Define Truncated Gaussian model with fixed baseline, peak, and spread, truncated within [-15, 180]
    gaussian_model = @(b, X) baseline_estimate + b * exp(-((X - peak_position).^2) / (2 * spread^2)) .* (X >= -15 & X <= 180);
    % baseline_estimate is fixed, b is the amplitude, peak is fixed at peak_position, spread is provided as input, active within -15 to 180

    % Initial guess for the amplitude parameter
    initial_guess_gaussian = 1;

    % Cross-validation for Gaussian model
    cv_gaussian = cvpartition(length(y_gaussian), 'KFold', numFolds);
    R2_folds_gaussian = zeros(numFolds, 1);

    parfor i = 1:numFolds
        trainIdx = training(cv_gaussian, i);
        testIdx = test(cv_gaussian, i);

        % Fit Gaussian model with fixed baseline, peak, and spread
        mdl_gaussian = fitnlm(X_gaussian(trainIdx), y_gaussian(trainIdx), gaussian_model, initial_guess_gaussian);
        y_pred_gaussian = predict(mdl_gaussian, X_gaussian(testIdx));
        R2_folds_gaussian(i) = calculate_r2(y_gaussian(testIdx), y_pred_gaussian);
    end

    % Mean R-squared across folds for the Gaussian model
    R2_truncatedGaussian = mean(R2_folds_gaussian);

    % Fit the final Gaussian model on the full dataset
    mdl_truncatedGaussian = fitnlm(X_gaussian, y_gaussian, gaussian_model, initial_guess_gaussian);

    %% Optional Plotting for Truncated Gaussian Model
    position_range = linspace(-180, 180, 100);
    truncated_gaussian_predictions = predict(mdl_truncatedGaussian, position_range');

    figure;
    plot(position_range, truncated_gaussian_predictions, 'r-', 'LineWidth', 2);
    xlabel('Visual Object Position (degrees)');
    ylabel('Predicted Firing Rate');
    title('Truncated Gaussian Model Predictions vs Visual Object Position');
    grid on;
end
