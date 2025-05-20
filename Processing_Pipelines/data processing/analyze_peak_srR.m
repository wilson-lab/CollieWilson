function analyze_peak_srR(folder, peak_srR_store)
    % Change to summary folder
    cd(folder.summary);

    % Remove flies (rows) containing NaN values across any bin or speed
    nan_mask = any(isnan(peak_srR_store), [2, 3]); % Check for NaNs in any bin or speed
    peak_srR_filtered = peak_srR_store(~nan_mask, :, :); % Keep only rows without NaNs

    % Reshape filtered data into a long format for ANOVA
    [num_flies, num_bins, num_speeds] = size(peak_srR_filtered);
    data_long = reshape(peak_srR_filtered, [], 1); % Flatten firing rate data
    fly_ids = repmat((1:num_flies)', [num_bins * num_speeds, 1]); % Fly identifiers
    bin_ids = repmat(repelem((1:num_bins)', num_flies), num_speeds, 1); % Bin identifiers
    speed_ids = repelem((1:num_speeds)', num_bins * num_flies); % Speed identifiers

    % Create a table for the ANOVA
    anova_table = table(fly_ids, bin_ids, speed_ids, data_long, ...
                        'VariableNames', {'Fly', 'Bin', 'Speed', 'FiringRate'});
    if ~isempty(anova_table)
    % Run three-way ANOVA
    [p, tbl, stats] = anovan(anova_table.FiringRate, ...
        {anova_table.Bin, anova_table.Speed}, ...
        'model', 'interaction', ...
        'varnames', {'Bin', 'Speed'}, ...
        'alpha', 0.05);

    % Perform multiple comparisons for each factor
    results_bin = multcompare(stats, 'Dimension', 1); % Bin
    results_speed = multcompare(stats, 'Dimension', 2); % Speed

    % Save results to Excel
    output_file = 'ANOVA_Results.xlsx';
    writetable(cell2table(tbl), output_file, 'Sheet', 'ANOVA_Table');
    writematrix(results_bin, output_file, 'Sheet', 'Bin_Comparisons');
    writematrix(results_speed, output_file, 'Sheet', 'Speed_Comparisons');

    disp(['ANOVA results and comparisons saved to ' fullfile(folder.summary, output_file)]);
    end
end
