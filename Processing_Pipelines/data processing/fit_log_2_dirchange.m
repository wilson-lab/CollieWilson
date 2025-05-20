function fit_log_2_dirchange(meanKIRDC, meanWTDC, meanNADC, settings)
    % Inputs:
    % - meanKIRDC, meanWTDC, meanNADC: datasets for each genotype
    % - settings: structure containing geneColor for plotting colors

    % Combine inputs into a structure for iteration
    genotypes = {'Kir', 'WT', 'Na'};
    data = {meanKIRDC, meanWTDC, meanNADC};
    colors = settings.geneColor; % Colors for each genotype

    slopes = cell(1, length(genotypes)); % Store slopes for each genotype
    r2_values = cell(1, length(genotypes)); % Store R^2 values for each genotype
    all_slopes = [];
    all_r2 = [];
    group_labels = [];
    animal_labels = [];

    % Fit slope for each animal using log base 2
    for g = 1:length(genotypes)
        genotype_data = data{g};
        [numRows, numAnimals] = size(genotype_data);
        slopes{g} = zeros(1, numAnimals);
        r2_values{g} = zeros(1, numAnimals);

        for a = 1:numAnimals
            x = (1:numRows)'; % Rotational velocities
            y = genotype_data(:, a); % Direction change times

            % Exclude invalid data
            validIdx = x > 0 & ~isnan(y);
            x = x(validIdx);
            y = y(validIdx);

            if ~isempty(x)
                % Fit the slope (a) with log base 2
                log_x = log2(x);
                p = polyfit(log_x, y, 1); % Linear fit
                fitted_y = polyval(p, log_x);
                
                % Calculate R^2
                ss_res = sum((y - fitted_y).^2); % Residual sum of squares
                ss_tot = sum((y - mean(y)).^2); % Total sum of squares
                r2 = 1 - (ss_res / ss_tot);

                slopes{g}(a) = p(1); % Store the slope
                r2_values{g}(a) = r2; % Store the R^2 value

                % Store for ANOVA
                all_slopes = [all_slopes; p(1)];
                all_r2 = [all_r2; r2];
                group_labels = [group_labels; g];
                animal_labels = [animal_labels; a + (g-1)*numAnimals];
            else
                slopes{g}(a) = NaN;
                r2_values{g}(a) = NaN;
            end
        end
    end

    % Create figure
    figure; set(gcf,'Position',[100 100 1200 400])
    tiledlayout(1, 3); % Single row, 3 tiles

    % First tile: Kernel density distribution of R^2 values
    nexttile;
    hold on;
    for g = 1:length(genotypes)
        r2_vals = r2_values{g};
        r2_vals = r2_vals(~isnan(r2_vals)); % Remove NaNs
        [density, xvals] = ksdensity(r2_vals, 'Bandwidth', 0.05);
        plot(xvals, density, 'Color', colors{g}, 'LineWidth', 2);
    end
    xlabel('R^2');
    ylabel('Density');
    legend(genotypes, 'Location', 'Best');
    title('Kernel Density of R^2 Values');
    hold off;

    % Second tile: Kernel density distribution of slopes
    nexttile;
    hold on;
    for g = 1:length(genotypes)
        slope_vals = slopes{g};
        slope_vals = slope_vals(~isnan(slope_vals)); % Remove NaNs
        [density, xvals] = ksdensity(slope_vals, 'Bandwidth', 1.5);
        plot(xvals, density, 'Color', colors{g}, 'LineWidth', 2);
    end
    xlabel('Slope (a)');
    ylabel('Density');
    legend(genotypes, 'Location', 'Best');
    title('Kernel Density of Slopes');
    hold off;

    % Third tile: Scatter plot with medians of slopes
    nexttile;
    hold on;
    x_padding = 0.5; % Padding for x-axis
    for g = 1:length(genotypes)
        scatter(repmat(g, size(slopes{g})), slopes{g}, '.', 'MarkerEdgeColor', colors{g});
        % Add median as a horizontal dash
        median_slope = median(slopes{g}, 'omitnan');
        plot(g, median_slope, '_', 'Color', 'k', 'LineWidth', 1);
    end
    set(gca, 'XTick', 1:length(genotypes), 'XTickLabel', genotypes);
    xlabel('Genotype');
    ylabel('Slope (a)');
    title('Slopes for Each Genotype');
    xlim([1-x_padding, length(genotypes)+x_padding]); % Pad x-axis
    hold off;

    % Perform ANOVA on slopes
    tbl = table(all_slopes, group_labels, animal_labels, ...
                'VariableNames', {'Slope', 'Genotype', 'Animal'});
    tbl.Genotype = categorical(tbl.Genotype);
    tbl.Animal = categorical(tbl.Animal);
    lme = fitlme(tbl, 'Slope ~ Genotype + (1|Animal)');
    p_gene = lme.anova.pValue(2); % p-value for genotype effect

    % Add p-value to the third plot
    nexttile(3);
    text(max(xlim)-0.5, max(ylim)-0.1, ...
         sprintf('p(gene) = %.3g', p_gene), ...
         'HorizontalAlignment', 'right', 'FontSize', 10);
end
