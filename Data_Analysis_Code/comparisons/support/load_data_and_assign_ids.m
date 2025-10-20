function [names, ids, rval_ang, lag_ang] = load_xc_overlaps(folder_motion, folder_bg, all_names, settings)
% Rows = overlapping animals; Cols = [bg1 (dark+arousal), bg2 (dark), motion].
% Outputs:
%   rval_ang : N x 3 matrix of peak correlation values
%   lag_ang  : N x 3 matrix of peak lag values
%
% CREATED: 10/19/2025 - MC

    % Determine year prefix
    if contains(folder_bg,'AOTU019')
        yp = '2023';
    elseif contains(folder_bg,'AOTU025')
        yp = '2024';
    else
        error('folder_bg must contain AOTU019 or AOTU025');
    end

    % Check which files exist
    has_m  = false(size(all_names));
    has_b1 = false(size(all_names));
    has_b2 = false(size(all_names));
    for i = 1:numel(all_names)
        nm = all_names{i};
        has_m(i)  = isfile(fullfile(folder_motion, [nm '_xc.mat']));
        has_b1(i) = isfile(fullfile(folder_bg, [yp '_' nm '_1_xc.mat']));
        has_b2(i) = isfile(fullfile(folder_bg, [yp '_' nm '_2_xc.mat']));
    end

    % Keep overlapping animals (motion + at least one background)
    keep = has_m & (has_b1 | has_b2);
    names = all_names(keep);
    has_b1 = has_b1(keep);
    has_b2 = has_b2(keep);
    N = numel(names);
    ids = (1:N)';

    % Initialize output matrices
    rval_ang = nan(N,3);
    lag_ang  = nan(N,3);

    % Load data
    for r = 1:N
        nm = names{r};

        if has_b1(r)
            S = load(fullfile(folder_bg, [yp '_' nm '_1_xc.mat']), 'r_val','lag_t');
            [pl, pr] = find_peak_lag_rval(S.r_val, S.lag_t, settings.minXCorrProm);
            rval_ang(r,1) = pr.ang;
            lag_ang(r,1)  = pl.ang;
        end

        if has_b2(r)
            S = load(fullfile(folder_bg, [yp '_' nm '_2_xc.mat']), 'r_val','lag_t');
            [pl, pr] = find_peak_lag_rval(S.r_val, S.lag_t, settings.minXCorrProm);
            rval_ang(r,2) = pr.ang;
            lag_ang(r,2)  = pl.ang;
        end

        S = load(fullfile(folder_motion, [nm '_xc.mat']), 'r_val','lag_t');
        [pl, pr] = find_peak_lag_rval(S.r_val, S.lag_t, settings.minXCorrProm);
        rval_ang(r,3) = pr.ang;
        lag_ang(r,3)  = pl.ang;
    end
end
