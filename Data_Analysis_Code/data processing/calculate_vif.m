% calculate_vif
% This function calculates the Variance Inflation Factor (VIF) for each predictor variable 
% in a matrix of predictors. VIF quantifies how much the variance of a regression coefficient 
% is inflated due to multicollinearity with other predictors.
%
% INPUTS:
%   X - Matrix of predictor variables (observations x predictors), where each column 
%       represents a different predictor variable.
%
% OUTPUTS:
%   vif_values - Array of VIF values for each predictor in the input matrix X.
%
% CREATED: [Date] MC
%
function vif_values = calculate_vif(X)
    % X is the matrix of predictors, where each column is a predictor variable
    % Initialize VIF values
    vif_values = zeros(1, size(X, 2));
    
    for j = 1:size(X, 2)
        % Select the j-th predictor
        X_j = X(:, j);
        
        % Select all other predictors except the j-th one
        X_others = X(:, [1:j-1, j+1:end]);
        
        % Fit a linear model to predict X_j using the other predictors
        mdl = fitlm(X_others, X_j);
        
        % Get the R-squared value from the model
        R2_j = mdl.Rsquared.Ordinary;
        
        % Calculate the VIF for the j-th predictor
        vif_values(j) = 1 / (1 - R2_j);
    end
end
