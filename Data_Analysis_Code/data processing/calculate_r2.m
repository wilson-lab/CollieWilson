% calculate_r2
% This function calculates the R-squared (coefficient of determination) value
% between observed and predicted values.
%
% INPUT
% y_true    - array of observed values
% y_pred    - array of predicted values
%
% OUTPUT
% R2        - R-squared value
%
% CREATED: 11/07/2024 - MC
%
function R2 = calculate_r2(y_true, y_pred)
    % Remove NaNs from both observed and predicted values
    valid_idx = ~isnan(y_true) & ~isnan(y_pred);
    y_true = y_true(valid_idx);
    y_pred = y_pred(valid_idx);

    % Calculate the total sum of squares
    SS_tot = sum((y_true - mean(y_true)).^2);

    % Calculate the residual sum of squares
    SS_res = sum((y_true - y_pred).^2);

    % Calculate R-squared
    R2 = 1 - (SS_res / SS_tot);
end
