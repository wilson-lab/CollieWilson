function [width_deg, halfwidth_deg, bounds_deg, N] = centeredArcWidth(theta_deg)
% CREATED: 10/17/2025 - MC
% Returns the median absolute angular deviation from 0°.

% ---- flatten and remove NaNs
x = theta_deg(:);
x = x(~isnan(x));
N = numel(x);
if N == 0
    width_deg = NaN; halfwidth_deg = NaN; bounds_deg = [NaN NaN];
    return
end

% ---- compute absolute circular distance from 0° (wrapped to [-180, 180])
dist_deg = abs(mod(x + 180, 360) - 180);

% ---- median deviation
halfwidth_deg = median(dist_deg);
width_deg = 2 * halfwidth_deg;
bounds_deg = [-halfwidth_deg, +halfwidth_deg];
end