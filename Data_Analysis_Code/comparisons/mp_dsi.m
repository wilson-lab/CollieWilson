% dsi_summary
% This script analyzes direction selectivity data for neurons AOTU019 and 
% AOTU025 at two target speeds. It loads preprocessed datasets, extracts
% direction selectivity index (DSI) values, and organizes the data into a
% table format for a three-way ANOVA. The main effects of cell type, speed,
% and fly ID on DSI are tested.
% 
% CREATED: 11/07/2024 - MC
%
%% Load direction selectivity data
plotPath = 'E:\Compare Motion Pulse';  % Define folder containing the .mat files
dataPath = 'E:\Compare Motion Pulse\data';  % Define folder containing the .mat files
cd(dataPath)

AOTU019_file = fullfile(dataPath, 'ds_AOTU019_Motion_Pulse.mat');
AOTU025_file = fullfile(dataPath, 'ds_AOTU025_Motion_Pulse.mat');

target_speeds = {'75dps';'25ps'};
mean_colors = {"#0072BD","#7E2F8E"};

%% Load in dataset
AOTU019_DSI = load(AOTU019_file);
AOTU025_DSI = load(AOTU025_file);

% Extract the store_ds data for each neuron
AOTU019_store_ds = AOTU019_DSI.store_ds;  % size: [2 x flies x 3]
AOTU025_store_ds = AOTU025_DSI.store_ds;

% Define number of flies for each cell type
num_flies_019 = size(AOTU019_store_ds, 2);
num_flies_025 = size(AOTU025_store_ds, 2);

% Define z condition names for labeling
condition_names = {'All', 'Quiescent', 'Pursuit'};
target_speeds = {'75dps', '25dps'};
mean_colors = {[0 0 1], [0.5 0 0.5]}; % Blue and Purple

% Create figure for 3 conditions × 2 speeds
figure; set(gcf, 'Position', [100 100 500 900])
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact')

% Loop over z-dimension (conditions)
for z = 1:3
    % Create vectors for data, CellType, Speed, and FlyID
    DSI_data = [AOTU019_store_ds(1, :, z)'; AOTU019_store_ds(2, :, z)'; ...
                AOTU025_store_ds(1, :, z)'; AOTU025_store_ds(2, :, z)'];

    CellType = [repmat("AOTU019", num_flies_019 * 2, 1); repmat("AOTU025", num_flies_025 * 2, 1)];

    Speed = [repmat("75dps", num_flies_019, 1); repmat("25dps", num_flies_019, 1); ...
             repmat("75dps", num_flies_025, 1); repmat("25dps", num_flies_025, 1)];

    FlyID = [repelem("Fly019_" + (1:num_flies_019)', 2); ...
             repelem("Fly025_" + (1:num_flies_025)', 2)];

    % Create table for repeated-measures ANOVA (within each cell type)
    data_table = table(DSI_data, CellType, Speed, FlyID, ...
        'VariableNames', {'DSI', 'CellType', 'Speed', 'FlyID'});

    % Convert to categorical
    data_table.CellType = categorical(data_table.CellType);
    data_table.Speed = categorical(data_table.Speed);
    data_table.FlyID = categorical(data_table.FlyID);

    % Fit mixed-effects model with CellType (between-subjects), Speed (within-subject), and interaction
    % Random intercept per fly to model repeated measures
    lme = fitlme(data_table, 'DSI ~ CellType*Speed + (1|FlyID)');

    % Extract ANOVA table
    lme_anova = anova(lme);

    % Extract p-values for each fixed effect
    p_celltype = lme_anova.pValue(strcmp(lme_anova.Term, 'CellType'));
    p_speed = lme_anova.pValue(strcmp(lme_anova.Term, 'Speed'));
    p_interaction = lme_anova.pValue(strcmp(lme_anova.Term, 'CellType:Speed'));

    % --- Post hoc comparisons: CellType at each Speed ---
fprintf('\nPost hoc cell type comparisons for condition: %s\n', condition_names{z});

for s = 1:2
    this_speed = target_speeds{s};

    % Subset data by speed
    T_sub = data_table(data_table.Speed == this_speed, :);

    % Fit model for cell type comparison at this speed
    lme_posthoc = fitlme(T_sub, 'DSI ~ CellType + (1|FlyID)');
    a = anova(lme_posthoc);

    % Extract p-value for cell type
    pval = a.pValue(strcmp(a.Term, 'CellType'));

    % Print result with stars
    stars = '';
    if pval < 0.0001
        stars = '****';
    elseif pval < 0.001
        stars = '***';
    elseif pval < 0.01
        stars = '**';
    elseif pval < 0.05
        stars = '*';
    end

    fprintf('  Speed %s: p = %.4f %s\n', this_speed, pval, stars);
end



    % Plotting loop over the 2 speeds
    for speed_idx = 1:2
        nexttile((z-1)*2 + speed_idx);  % Tile index across rows
        jitter_amount = 0.1;

        % Extract data
        AOTU019_data = AOTU019_store_ds(speed_idx, :, z);
        AOTU025_data = AOTU025_store_ds(speed_idx, :, z);

        % Scatter individual points
        scatter(ones(size(AOTU019_data)) + jitter_amount * (rand(size(AOTU019_data)) - 0.5), AOTU019_data, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
        hold on;
        scatter(2 * ones(size(AOTU025_data)) + jitter_amount * (rand(size(AOTU025_data)) - 0.5), AOTU025_data, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);

        % Plot means
        plot(1, median(AOTU019_data, 'omitnan'), '_', 'MarkerSize', 15, 'Color', mean_colors{1}, 'LineWidth', 2);
        plot(2, median(AOTU025_data, 'omitnan'), '_', 'MarkerSize', 15, 'Color', mean_colors{2}, 'LineWidth', 2);

        % Labels and limits
        set(gca, 'XTick', [1 2], 'XTickLabel', {'AOTU019', 'AOTU025'});
        ylabel('DSI');
        title([condition_names{z}, ' | Speed: ', target_speeds{speed_idx}]);
        xlim([0, 3]);
        ylim([-0.5 1]);
        yline(0)
    end

    % Add p-values to last tile of the row
    axes_handle = gca;
    text(axes_handle, 2.5, 1, sprintf('p(celltype): %.3f', p_celltype), ...
         'FontSize', 8, 'HorizontalAlignment', 'right');
    text(axes_handle, 2.5, 0.9, sprintf('p(speed): %.3f', p_speed), ...
         'FontSize', 8, 'HorizontalAlignment', 'right');
    text(axes_handle, 2.5, 0.8, sprintf('p(interaction): %.3f', p_interaction), ...
         'FontSize', 8, 'HorizontalAlignment', 'right');
end


% Save plots
cd(plotPath)
saveas(gcf, 'mp_dsi_conditions.png');
set(gcf, 'renderer', 'Painters');
saveas(gcf, 'mp_dsi_conditions.svg');

%% Plot DSI across positions

% Extract data
data019 = AOTU019_DSI.store_all25_ds; % size: x (animals) x y (positions) x z (conditions)
data025_full = AOTU025_DSI.store_all25_ds;
% (Optional) Trim AOTU025 data
data025 = data025_full(:,3:end,:);

% Get the number of sweep positions in AOTU019 to trim AOTU025
nPos019 = size(data019, 2);
nPos025 = size(data025, 2);
% Define x-axis values (panel positions)
panel_positions = [-11.25, 11.25, 33.75, 56.25, 78.75, 101.25, 123.75];
panel_positions019 = panel_positions(1:nPos019);
panel_positions025 = panel_positions(3:end);

% Define colors
color019 = [0 0.447 0.741];     % blue
color025 = [0.494 0.184 0.556]; % purple

% Prepare figure
figure; set(gcf, 'Position', [100 100 300 900])
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact')
cond_names = {'All', 'Quiescent', 'Pursuit'};

for z = 1:3
    nexttile;
    
    % Get condition data
    d019 = squeeze(data019(:,:,z)); % x = animals, y = positions
    d025 = squeeze(data025(:,:,z));
    
    % Compute mean and SEM across animals
    m019 = mean(d019, 1, 'omitnan');
    sem019 = std(d019, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(d019), 1));
    
    m025 = mean(d025, 1, 'omitnan');
    sem025 = std(d025, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(d025), 1));
    
    % Plot error bars without markers or caps
    hold on;
    errorbar(panel_positions019, m019, sem019, '-', ...
        'Color', color019, 'LineWidth', 1.5, 'CapSize', 0);

    errorbar(panel_positions025, m025, sem025, '-', ...
        'Color', color025, 'LineWidth', 1.5, 'CapSize', 0);

    % Formatting
    title(cond_names{z});
    xlabel('Sweep Position');
    ylabel('DS Index');
    ylim([-0.5 1]);
    yline(0)
    xticks([0 30 60 90 120 150])
    xlim([min(panel_positions)-5, max(panel_positions)+5]);
    box off;
end
% Save plots
cd(plotPath)
saveas(gcf, 'mp_dsi__across_conditions.png');
set(gcf, 'renderer', 'Painters');
saveas(gcf, 'mp_dsi_across_conditions.svg');

%% Load DNa02 direction selectivity data
DNa02_file = fullfile(dataPath, 'ds_DNa02_Motion_Pulse.mat');
DNa02_DSI = load(DNa02_file);
DNa02_store_ds = DNa02_DSI.store_ds;  % size: [1 x flies x 3]
num_flies_DNa02 = size(DNa02_store_ds, 2);

% Define color for DNa02
color_DNa02 = [0.25 0.25 0.25]; % dark gray

% Create figure
figure; set(gcf, 'Position', [700 100 300 900])
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact')

for z = 1:3
    nexttile;
    jitter_amount = 0.025;

    % Extract data for condition z at speed 25dps (row 1)
    DNa02_data = squeeze(DNa02_store_ds(1,:,z));

    % Scatter individual points
    scatter(ones(size(DNa02_data)) + jitter_amount*(rand(size(DNa02_data))-0.5), DNa02_data, '.', ...
        'MarkerEdgeColor', [0.5 0.5 0.5]);
    hold on;

    % Plot mean
    plot(1, mean(DNa02_data, 'omitnan'), '_', 'MarkerSize', 15, 'Color', color_DNa02, 'LineWidth', 2);

    % Labels and limits
    set(gca, 'XTick', 1, 'XTickLabel', {'25dps'});
    ylabel('DSI');
    title(['DNa02 | ' condition_names{z}]);
    xlim([0.5 1.5]);
    ylim([0 2]);
    yline(0)
end

% Save plots
cd(plotPath)
saveas(gcf, 'mp_dsi_dna02_only.png');
set(gcf, 'renderer', 'Painters');
saveas(gcf, 'mp_dsi_dna02_only.svg');
