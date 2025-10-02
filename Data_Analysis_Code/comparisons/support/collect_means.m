function [flyID, noP1_mean, P1_mean] = collect_means(P1_arr, noP1_arr, tag)
    nFlies = size(P1_arr, 2);
    flyID = strings(0,1); noP1_mean = []; P1_mean = [];
    for i = 1:nFlies
        P1_mid   = P1_arr(start_idx:end_idx, i, :);
        noP1_mid = noP1_arr(start_idx:end_idx, i, :);

        mP1   = mean(P1_mid(:),   'omitnan');
        mNoP1 = mean(noP1_mid(:), 'omitnan');

        if isfinite(mP1) && isfinite(mNoP1)
            flyID(end+1,1)   = sprintf('%s_%02d', tag, i); %#ok<AGROW>
            noP1_mean(end+1) = mNoP1;                      %#ok<AGROW>
            P1_mean(end+1)   = mP1;                        %#ok<AGROW>
        end
    end
end