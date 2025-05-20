function run_forward_velocity_anova(panelps, kir, wt, na)
% CREATED: 04/17/2025 - MC
% Bins forward velocity by panel position using 6° bins (-36 to +36),
% subtracts per-fly baseline, and runs repeated-measures ANOVA with
% Genotype × PanelPos interaction and FlyID as random effect.

% Trim to middle 50% of data
nTime = length(panelps);
q1 = round(nTime * 0.25);
q3 = round(nTime * 0.75);
midIdx = q1:q3;

panelps = panelps(midIdx);
kir = kir(midIdx, :);
wt = wt(midIdx, :);
na = na(midIdx, :);

allData = {kir, wt, na};
genotypes = {'KIR', 'WT', 'NA'};
nGroups = numel(allData);

% Define bin edges and centers (6° bins)
edges = -37.5:1.75:37.5;
binCenters = edges(1:end-1) + 3;
nBins = numel(binCenters);

% Preallocate long format data
ForwardVel = [];
PanelPos = [];
Genotype = {};
FlyID = {};

flyCounter = 1;

for g = 1:nGroups
    data = allData{g};  % time x flies
    nFlies = size(data, 2);

    for f = 1:nFlies
        flyData = data(:, f);
        flyMean = mean(flyData);  % baseline per fly

        for b = 1:nBins
            idx = panelps >= edges(b) & panelps < edges(b+1);
            if any(idx)
                thisBinVal = mean(flyData(idx)) - flyMean;
                ForwardVel(end+1,1) = thisBinVal;
                PanelPos(end+1,1) = binCenters(b);
                Genotype{end+1,1} = genotypes{g};
                FlyID{end+1,1} = sprintf('fly%d', flyCounter);
            end
        end
        flyCounter = flyCounter + 1;
    end
end

% Build table for fitlme
T = table(ForwardVel, categorical(Genotype), ...
          categorical(FlyID), PanelPos, ...
          'VariableNames', {'ForwardVel','Genotype','FlyID','PanelPos'});

% Fit linear mixed-effects model
T.PanelPos = categorical(T.PanelPos);
lme = fitlme(T, 'ForwardVel ~ Genotype * PanelPos + (1|FlyID)');

disp('Repeated Measures ANOVA via Linear Mixed-Effects Model:')
disp(anova(lme))

% Post-hoc Tukey (collapsed)
disp('Tukey post-hoc test across Genotypes (collapsed across bins):')
T_avg = groupsummary(T, {'FlyID','Genotype'}, 'mean', 'ForwardVel');
[p, tbl, stats] = anova1(T_avg.mean_ForwardVel, T_avg.Genotype, 'off');
posthoc = multcompare(stats, 'CType', 'tukey-kramer');
disp(array2table(posthoc, ...
    'VariableNames', {'Group1','Group2','Lower','Estimate','Upper','pValue'}));

% --------- Plotting ---------
colors = lines(nGroups);
panelList = categories(T.PanelPos);
panelVals = str2double(panelList);

figure;
hold on
for g = 1:nGroups
    thisGeno = genotypes{g};
    mask = T.Genotype == thisGeno;

    means = nan(size(panelVals));
    sems = nan(size(panelVals));

    for i = 1:numel(panelVals)
        vals = T.ForwardVel(mask & T.PanelPos == panelList{i});
        means(i) = mean(vals);
        sems(i) = std(vals) / sqrt(numel(vals));
    end

    % Plot line
    plot(panelVals, means, '-', 'Color', colors(g,:), 'LineWidth', 1.5, ...
         'DisplayName', thisGeno)

    % Plot patch for SEM (like your example)
    x_patch = [panelVals; flipud(panelVals)];
    y_patch = [means - sems; flipud(means + sems)];
    sp = patch(x_patch, y_patch, 'r', ...
               'FaceAlpha', 0.1, 'EdgeColor', 'none');
    sp.FaceColor = colors(g,:);
end

xlabel('Panel Position (deg)')
ylabel('Baseline-Subtracted Forward Velocity (mm/s)')
legend('Location', 'best')
title('Forward Velocity vs Panel Position (Baseline-Subtracted)')
xline(0, 'k--', 'LineWidth', 1)
yline(0, 'k--', 'LineWidth', 1)

% --------- Pairwise Post-hoc Comparisons ---------
fprintf('\nPost-hoc pairwise comparisons of tuning curves (Genotype × PanelPos):\n')
genoList = categories(T.Genotype);

for i = 1:numel(genoList)
    for j = i+1:numel(genoList)
        g1 = genoList{i};
        g2 = genoList{j};

        mask = T.Genotype == g1 | T.Genotype == g2;
        T_sub = T(mask, :);

        T_sub.Genotype = removecats(categorical(T_sub.Genotype));
        T_sub.PanelPos = removecats(categorical(T_sub.PanelPos));

        try
            lme_pair = fitlme(T_sub, 'ForwardVel ~ Genotype * PanelPos + (1|FlyID)');
            anova_tbl = anova(lme_pair);
            pval = anova_tbl.pValue(strcmp(anova_tbl.Term, 'Genotype:PanelPos'));
            fprintf('%s vs %s: interaction p = %.4e\n', g1, g2, pval)
        catch ME
            fprintf('%s vs %s: model fitting failed (%s)\n', g1, g2, ME.message)
        end
    end
end
