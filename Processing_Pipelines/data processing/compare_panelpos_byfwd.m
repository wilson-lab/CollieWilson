% compare_panelpos_by_speed
%
% This function compares histograms of panel position during slow versus fast walking periods
% for a given animal. Forward velocity is divided into tertiles, and only the bottom and top thirds
% are used to define "slow" and "fast" walking, respectively.
%
% INPUTS:
%   panelpos - Panel position data (time x trials)
%   forward  - Forward velocity data (time x trials)
%
% OUTPUTS:
%   hist_slow - Histogram of panel positions during slow walking (lowest tertile)
%   hist_fast - Histogram of panel positions during fast walking (highest tertile)
%
% CREATED: 04/16/2025 - MC
%
function [hist_slow, hist_fast] = compare_panelpos_byfwd(panelpos, forward)

    % Reshape data into column vectors
    panelpos_flat = reshape(panelpos, [], 1);
    forward_flat = reshape(forward, [], 1);

    % Only include time points where forward velocity is > 0
    validIdx = forward_flat > 0 & ~isnan(forward_flat) & ~isnan(panelpos_flat);
    forward_valid = forward_flat(validIdx);
    panelpos_valid = panelpos_flat(validIdx);

    % Return empty if there's no valid data
    if isempty(forward_valid)
        hist_slow = [];
        hist_fast = [];
        return
    end

    % Define tertiles
    cutoff1 = quantile(forward_valid, 1/3);
    cutoff2 = quantile(forward_valid, 2/3);

    % Get indices for slow and fast walking
    slowIdx = forward_valid <= cutoff1;
    fastIdx = forward_valid > cutoff2;

    % Compute histograms for slow and fast walking
    [hist_slow, ~] = panel_histogram(panelpos_valid(slowIdx), [], 1); % normalize
    [hist_fast, ~] = panel_histogram(panelpos_valid(fastIdx), [], 1); % normalize
end
