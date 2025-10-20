% This function compares angular arc width of panel position during slow versus fast walking periods
% for a given animal. Forward velocity is divided into tertiles, and only the bottom and top thirds
% are used to define "slow" and "fast" walking, respectively.
%
% INPUTS:
%   panelpos - Panel position data (time x trials)
%   forward  - Forward velocity data (time x trials)
%
% CREATED: 09/29/25 - MC
%
function [arc_slow, arc_fast] = arc_byfwd(panelpos, forward)

    % Reshape data into column vectors
    panelpos_flat = reshape(panelpos, [], 1);
    forward_flat = reshape(forward, [], 1);

    % Only include time points where forward velocity is > 0
    validIdx = forward_flat > 0 & ~isnan(forward_flat) & ~isnan(panelpos_flat);
    forward_valid = forward_flat(validIdx);
    panelpos_valid = panelpos_flat(validIdx);

    % Return empty if there's no valid data
    if isempty(forward_valid)
        arc_slow = [];
        arc_fast = [];
        return
    end

    % Define tertiles
    cutoff1 = quantile(forward_valid, 1/3);
    cutoff2 = quantile(forward_valid, 2/3);

    % Get indices for slow and fast walking
    slowIdx = forward_valid <= cutoff1;
    fastIdx = forward_valid > cutoff2;

    % Convert panel positions to radians (wrapped to 0–360°)
    ang_slow = panelpos_valid(slowIdx);
    ang_fast = panelpos_valid(fastIdx);

    % Compute 50% arc width
    [~, arc_slow, ~, ~] = centeredArcWidth(ang_slow);
    [~, arc_fast, ~, ~] = centeredArcWidth(ang_fast);
end
