function result = compare_conditions(lme)
    % Extract fixed effects
    fe = fixedEffects(lme);
    covFE = lme.CoefficientCovariance;

    % Get condition names
    levels = categories(lme.Variables.Condition);

    % Build contrast matrix for all pairwise comparisons
    nLevels = numel(levels);
    result = table();
    idx = 1;

    for i = 1:nLevels
        for j = i+1:nLevels
            % Create contrast vector
            c = zeros(size(fe));
            c(strcmp(lme.CoefficientNames, ['Condition_' levels{j}])) = 1;
            c(strcmp(lme.CoefficientNames, ['Condition_' levels{i}])) = -1;

            % Estimate difference and SE
            diff = c' * fe;
            se = sqrt(c' * covFE * c);
            tval = diff / se;
            df = lme.DFE;
            pval = 2 * tcdf(-abs(tval), df);

            % Store
            result(idx,:) = table(levels{i}, levels{j}, diff, se, tval, df, pval, ...
                'VariableNames', {'Group1', 'Group2', 'Estimate', 'SE', 'tStat', 'DF', 'pValue'});
            idx = idx + 1;
        end
    end
end
