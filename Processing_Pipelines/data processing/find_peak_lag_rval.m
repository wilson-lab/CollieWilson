% Function to find the peak lag, peak r-value, and handle r_val output
% Inputs:
%   - r_val: Structure containing cross-correlation values for forward (fwd),
%            angular (ang), and sideway (sid) motion.
%   - lag_t: Array of lag times corresponding to r_val.
%   - minProm: Minimum prominence threshold for peak detection.
% Outputs:
%   - peak_lag: Structure with peak lag values or NaN if no peak is prominent.
%   - peak_rval: Structure with peak r-values or NaN if no peak is prominent.
%   - r_val_out: Structure with same fields as r_val, with NaN for fields
%                where no prominent peak was detected.
%
% Created: 11/08/2024 - MC
%
function [peak_lag, peak_rval, r_val_out] = find_peak_lag_rval(r_val, lag_t, minProm)
    % Initialize output structures
    peak_lag = struct('fwd', nan, 'ang', nan, 'sid', nan);
    peak_rval = struct('fwd', nan, 'ang', nan, 'sid', nan);
    r_val_out = struct('fwd', nan(size(r_val.fwd)), 'ang', nan(size(r_val.ang)), 'sid', nan(size(r_val.sid)));

    % Detect peaks in forward correlation
    [f_rpk, f_locs] = findpeaks(r_val.fwd, 'MinPeakProminence', minProm, 'SortStr', 'descend');
    if ~isempty(f_locs)
        peak_lag.fwd = lag_t(f_locs(1));
        peak_rval.fwd = f_rpk(1);
        r_val_out.fwd = r_val.fwd;  % Keep original r_val if peak is found
    end

    % Detect peaks in angular correlation
    [a_rpk, a_locs] = findpeaks(r_val.ang, 'MinPeakProminence', minProm, 'SortStr', 'descend');
    if ~isempty(a_locs)
        peak_lag.ang = lag_t(a_locs(1));
        peak_rval.ang = a_rpk(1);
        r_val_out.ang = r_val.ang;
    end

    % Detect peaks in sideway correlation
    [s_rpk, s_locs] = findpeaks(r_val.sid, 'MinPeakProminence', minProm, 'SortStr', 'descend');
    if ~isempty(s_locs)
        peak_lag.sid = lag_t(s_locs(1));
        peak_rval.sid = s_rpk(1);
        r_val_out.sid = r_val.sid;
    end
end
