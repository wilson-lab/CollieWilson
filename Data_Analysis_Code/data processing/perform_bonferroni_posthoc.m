function comparisons = perform_bonferroni_posthoc(lme, levels, n_comparisons, factor_name)
    comparisons = {};
    predictions = zeros(length(levels), 1);

    % Predict response for each level of the factor
    for i = 1:length(levels)
        temp_data = lme.Variables(1, :);  % template row
        temp_data.(factor_name)(:) = categorical(string(levels(i)));
        predictions(i) = predict(lme, temp_data);
    end

    % Pairwise comparisons
    for i = 1:length(levels)
        for j = i+1:length(levels)
            diff_value = predictions(i) - predictions(j);
            se_diff = sqrt(lme.CoefficientCovariance(i, i) + ...
                           lme.CoefficientCovariance(j, j) - ...
                           2 * lme.CoefficientCovariance(i, j));
            t_stat = diff_value / se_diff;
            p_value = 2 * (1 - tcdf(abs(t_stat), lme.DFE));
            p_adj = min(p_value * n_comparisons, 1);  % Bonferroni adjustment

            comparisons = [comparisons; {char(levels(i)), char(levels(j)), ...
                diff_value, se_diff, t_stat, p_adj}];
        end
    end
end
