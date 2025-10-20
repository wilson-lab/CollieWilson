function local_plot_by_condition(tbl, condLabels, ttl, ylab, ylim_use)
    vals = table2array(tbl);          % rows = animals, cols = conditions
    nCond = size(vals,2);
    hold on; box off;
    % jittered dots
    for c = 1:nCond
        y = vals(:,c);
        x = c + 0.1*(rand(size(y))-0.5);  % jitter width ~0.18
        plot(x, y, 'k.', 'MarkerSize', 5);
        % median as a solid dash (short horizontal line)
        medc = median(y,'omitnan');
        if ~isnan(medc)
            plot([c-0.25, c+0.25], [medc, medc], '-', 'LineWidth', 2);
        end
    end
    xlim([0.5, nCond+0.5]);
    xticks(1:nCond); xticklabels(condLabels); xtickangle(15);
    if nargin >= 5 && all(isfinite(ylim_use)), ylim(ylim_use); end
    ylabel(ylab);
    title(ttl);
    grid on;
end