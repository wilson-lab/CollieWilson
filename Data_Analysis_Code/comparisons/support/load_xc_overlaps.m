function [names, ids, rval_ang, lag_ang] = load_xc_overlaps(folder_motion, folder_bg, all_names, settings)
% Rows = overlapping animals; Cols = [bg1 (dark+arousal), bg2 (dark), motion].
% Per-animal year prefixes are auto-detected from filenames in folder_bg.

    % Pre-scan motion & background availability, capturing per-animal year prefixes
    Nall = numel(all_names);
    has_m  = false(Nall,1);
    has_b1 = false(Nall,1);
    has_b2 = false(Nall,1);
    yp_b1  = repmat({''}, Nall, 1);
    yp_b2  = repmat({''}, Nall, 1);

    for i = 1:Nall
        nm = all_names{i};
        has_m(i) = isfile(fullfile(folder_motion, [nm '_xc.mat']));

        % Find bg1 and extract its year prefix for this animal if present
        L1 = dir(fullfile(folder_bg, ['*_' nm '_1_xc.mat']));
        if ~isempty(L1)
            has_b1(i) = true;
            yp_b1{i} = regexp(L1(1).name, '^(\d{4})_', 'tokens','once');
            yp_b1{i} = yp_b1{i}{1}; % 4-digit year
        end

        % Find bg2 and extract its year prefix for this animal if present
        L2 = dir(fullfile(folder_bg, ['*_' nm '_2_xc.mat']));
        if ~isempty(L2)
            has_b2(i) = true;
            yp_b2{i} = regexp(L2(1).name, '^(\d{4})_', 'tokens','once');
            yp_b2{i} = yp_b2{i}{1}; % 4-digit year
        end
    end

    % Keep animals that overlap: motion + at least one background
    keep = has_m & (has_b1 | has_b2);
    names = all_names(keep);
    yp_b1 = yp_b1(keep);
    yp_b2 = yp_b2(keep);
    has_b1 = has_b1(keep);
    has_b2 = has_b2(keep);

    N = numel(names);
    ids = (1:N)';

    % Preallocate outputs (N x 3). Cols: [bg1, bg2, motion]
    rval_ang = nan(N,3);
    lag_ang  = nan(N,3);

    % Load each animal
    for r = 1:N
        nm = names{r};

        % bg1: dark + arousal
        if has_b1(r)
            f1 = fullfile(folder_bg, [yp_b1{r} '_' nm '_1_xc.mat']);
            S = load(f1, 'r_val','lag_t');
            [pl, pr] = find_peak_lag_rval(S.r_val, S.lag_t, settings.minXCorrProm);
            rval_ang(r,1) = pr.ang;  lag_ang(r,1) = pl.ang;
        end

        % bg2: dark (no arousal)
        if has_b2(r)
            f2 = fullfile(folder_bg, [yp_b2{r} '_' nm '_2_xc.mat']);
            S = load(f2, 'r_val','lag_t');
            [pl, pr] = find_peak_lag_rval(S.r_val, S.lag_t, settings.minXCorrProm);
            rval_ang(r,2) = pr.ang;  lag_ang(r,2) = pl.ang;
        end

        % motion: with arousal (no year prefix)
        fm = fullfile(folder_motion, [nm '_xc.mat']);
        S = load(fm, 'r_val','lag_t');
        [pl, pr] = find_peak_lag_rval(S.r_val, S.lag_t, settings.minXCorrProm);
        rval_ang(r,3) = pr.ang;    lag_ang(r,3) = pl.ang;
    end
end
