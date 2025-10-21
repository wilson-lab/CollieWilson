function generateSummaryMetrics(xvalues, metrics_prob, comparisonLabel, comparisonType, xvariable, folder)
% GENERATESUMMARYPROB
% CREATED: 10/20/2025 - MC
%
% Plots probability near setpoint across xvalues for multiple conditions.
% X-axis label is provided via xvariable.

%% SETUP
numConds = size(metrics_prob, 2);
colors   = lines(numConds);

%% PLOT PROBABILITY
figure; set(gcf, 'Position', [100 100 700 450]);
hold on;
for c = 1:numConds
    plot(xvalues, metrics_prob(:, c), '-', ...
        'Color', colors(c, :), 'LineWidth', 1);
end
xlabel(xvariable);
ylabel('Probability Near Setpoint');
legend(comparisonLabel, 'Location', 'best');
title(['Probability vs ' xvariable ' — ' comparisonType]);
ylim([0 1]);

%% SAVE
set(gcf, 'Renderer', 'painters');
saveas(gcf, fullfile(folder.final,  [comparisonType '_SummaryProb.png']));
saveas(gcf, fullfile(folder.vectors,[comparisonType '_SummaryProb.svg']));
end
