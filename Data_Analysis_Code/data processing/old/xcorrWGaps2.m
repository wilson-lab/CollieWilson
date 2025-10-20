% xcorrWGaps.m
%
% Function to compute cross-correlation when there are non-random gaps in
%  the data.
% Operates on 2 input vectors. Gaps in vectors should be NaNs. Input
%  vectors must have same sampling rate.
% Uses corr function to calculate Pierson's correlation for each lag, after
%  removing missing data points.
%
% INPUTS:
%   x - first input vector
%   y - second input vector
%   numLags - number of lags to calculate, should be odd number (0 plus
%       equal number of lags on each side of 0)
%
% OUTPUTS:
%   allCorr - all correlations, of length numLags  
%   lags - all lags used, in samples, matched to allCorr. Negative lags are
%       later time points in x correlated with earlier time points in y.
%
% CREATED: 8/31/23 - HHY
%
% UPDATED:
%   8/31/23 - HHY
%
function [allCorr, lags] = xcorrWGaps2(x, y, numLags)

    % convert x and y to columns
    if ~iscolumn(x)
        x = x';
    end
    if ~iscolumn(y)
        y = y';
    end

    % get number of lags to each side of 0
    numLags1Side = ceil((numLags-1)/2);
    lags = (-numLags1Side:numLags1Side)';

    % precompute NaN logical arrays
    xNaNLog = isnan(x);
    yNaNLog = isnan(y);
    allNaNLog = xNaNLog | yNaNLog;

    % Remove NaNs
    x(allNaNLog) = NaN;
    y(allNaNLog) = NaN;

    % Create lagged versions of x and y
    n = length(x);
    laggedX = NaN(n, numLags1Side*2 + 1);
    laggedY = NaN(n, numLags1Side*2 + 1);

    for i = 1:numLags1Side
        laggedX((i+1):end, numLags1Side + 1 - i) = x(1:end-i);
        laggedY(1:end-i, numLags1Side + 1 - i) = y((i+1):end);
        laggedX(1:end-i, numLags1Side + 1 + i) = x((i+1):end);
        laggedY((i+1):end, numLags1Side + 1 + i) = y(1:end-i);
    end

    laggedX(:, numLags1Side + 1) = x;
    laggedY(:, numLags1Side + 1) = y;

    % Compute correlations
    allCorr = NaN(numLags1Side * 2 + 1, 1);
    for i = 1:(numLags1Side*2 + 1)
        validMask = ~isnan(laggedX(:, i)) & ~isnan(laggedY(:, i));
        if any(validMask)
            allCorr(i) = corr(laggedX(validMask, i), laggedY(validMask, i));
        end
    end
end
