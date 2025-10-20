% Comparison of Cross-Correlation Peak Values and Lag Times for AOTU019 and AOTU025
% 
% This script processes cross-correlation data for two Drosophila neurons,
% AOTU019 and AOTU025, under two experimental conditions: motion pulse and
% background pulse. It performs a two-way ANOVA to evaluate differences in
% peak correlation values and lag times based on condition (dark vs. pursuit) 
% and animal ID, where condition and animal ID are main effects (no interaction).
%
% Created: 11/08/2024 MC
%
%% Initialize parameters
% Load settings
close all; clear
settings = processSettings();
folder_plots = 'E:\Compare Motion Pulse';

% Specify folders for AOTU019 and AOTU025 data
folder_motionpulse_019 = 'E:\AOTU019 Motion Pulse\interpolated\xcorr';
folder_backgroundpulse_019 = 'E:\AOTU019 Background P1\interpolated\xcorr';
folder_motionpulse_025 = 'E:\AOTU025 Motion Pulse\interpolated\xcorr';
folder_backgroundpulse_025 = 'E:\AOTU025 Background P1\interpolated\xcorr';

% Get list of files for AOTU019
files_motionpulse_019 = dir(fullfile(folder_motionpulse_019, '*_xc.mat'));
files_backgroundpulse_019 = dir(fullfile(folder_backgroundpulse_019, '*_1_xc.mat'));

% Get list of files for AOTU025
files_motionpulse_025 = dir(fullfile(folder_motionpulse_025, '*_xc.mat'));
files_backgroundpulse_025 = dir(fullfile(folder_backgroundpulse_025, '*_1_xc.mat'));

% Extract animal names from filenames for both conditions
motionpulse_names_019 = extract_animal_names(files_motionpulse_019, '^(.*)_xc');
backgroundpulse_names_019 = extract_animal_names(files_backgroundpulse_019, '(?<=2023_)(.*)_1_xc');
motionpulse_names_025 = extract_animal_names(files_motionpulse_025, '^(.*)_xc');
backgroundpulse_names_025 = extract_animal_names(files_backgroundpulse_025, '(?<=202[45]_)(.*)_1_xc');

% Get unique names across conditions for AOTU019 and AOTU025
all_names_019 = union(motionpulse_names_019, backgroundpulse_names_019);
all_names_025 = union(motionpulse_names_025, backgroundpulse_names_025);

% Assign unique animal IDs, shared across conditions where applicable
[animal_ids_019, animal_conditions_019, rval_data_019, lag_data_019] = load_data_and_assign_ids(folder_motionpulse_019, folder_backgroundpulse_019, all_names_019, settings);
[animal_ids_025, animal_conditions_025, rval_data_025, lag_data_025] = load_data_and_assign_ids(folder_motionpulse_025, folder_backgroundpulse_025, all_names_025, settings);

%% Perform two-way ANOVA on peak correlation and lag data for each neuron
% Analyze Angular Velocity
% AOTU019: Two-way ANOVA for Peak Correlation (Angular)
[p_rval_019_ang, tbl_rval_019_ang] = anovan(rval_data_019(:,1), {animal_conditions_019, animal_ids_019}, ...
                                    'model', 'linear', 'varnames', {'Condition', 'AnimalID'}, 'display', 'off');

% AOTU019: Two-way ANOVA for Lag (Angular)
[p_lag_019_ang, tbl_lag_019_ang] = anovan(lag_data_019(:,1), {animal_conditions_019, animal_ids_019}, ...
                                  'model', 'linear', 'varnames', {'Condition', 'AnimalID'}, 'display', 'off');

% AOTU025: Two-way ANOVA for Peak Correlation (Angular)
[p_rval_025_ang, tbl_rval_025_ang] = anovan(rval_data_025(:,1), {animal_conditions_025, animal_ids_025}, ...
                                    'model', 'linear', 'varnames', {'Condition', 'AnimalID'}, 'display', 'off');

% AOTU025: Two-way ANOVA for Lag (Angular)
[p_lag_025_ang, tbl_lag_025_ang] = anovan(lag_data_025(:,1), {animal_conditions_025, animal_ids_025}, ...
                                  'model', 'linear', 'varnames', {'Condition', 'AnimalID'}, 'display', 'off');

% Analyze Forward Velocity
% AOTU019: Two-way ANOVA for Peak Correlation (Forward)
[p_rval_019_fwd, tbl_rval_019_fwd] = anovan(rval_data_019(:,2), {animal_conditions_019, animal_ids_019}, ...
                                    'model', 'linear', 'varnames', {'Condition', 'AnimalID'}, 'display', 'off');

% AOTU019: Two-way ANOVA for Lag (Forward)
[p_lag_019_fwd, tbl_lag_019_fwd] = anovan(lag_data_019(:,2), {animal_conditions_019, animal_ids_019}, ...
                                  'model', 'linear', 'varnames', {'Condition', 'AnimalID'}, 'display', 'off');

% AOTU025: Two-way ANOVA for Peak Correlation (Forward)
[p_rval_025_fwd, tbl_rval_025_fwd] = anovan(rval_data_025(:,2), {animal_conditions_025, animal_ids_025}, ...
                                    'model', 'linear', 'varnames', {'Condition', 'AnimalID'}, 'display', 'off');

% AOTU025: Two-way ANOVA for Lag (Forward)
[p_lag_025_fwd, tbl_lag_025_fwd] = anovan(lag_data_025(:,2), {animal_conditions_025, animal_ids_025}, ...
                                  'model', 'linear', 'varnames', {'Condition', 'AnimalID'}, 'display', 'off');

%% Compare Cross-Correlation Data with Medians
% Define colors for each condition
dark_color = [0 0 0];          % Black for dark condition
pursuit_color_AOTU019 = [0 0 1]; % Blue for AOTU019 pursuit condition
pursuit_color_AOTU025 = [0.5 0 0.5]; % Purple for AOTU025 pursuit condition

% Create a figure
figure; set(gcf, 'Position', [100 100 600 800])
tiledlayout(2, 2, 'TileSpacing', 'compact')

% Jitter amount for scatter plots
jitter_amount = 0.1;
y_limits = [0 1]; % Y-axis limit for all plots
x_limits = [0 3]; % X-axis limit for all plots
p_text_fontsize = 8; % Font size for p-value text

% --- AOTU019 Angular Velocity Peak Correlation ---
nexttile;
scatter(animal_conditions_019 + jitter_amount * (rand(size(rval_data_019(:,1))) - 0.5), ...
        rval_data_019(:,1), '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
median_rval_019_dark = median(rval_data_019(animal_conditions_019 == 1,1), 'omitnan');
median_rval_019_pursuit = median(rval_data_019(animal_conditions_019 == 2,1), 'omitnan');

plot(1, median_rval_019_dark, '_', 'MarkerSize', 15, 'Color', dark_color, 'LineWidth', 2);
plot(2, median_rval_019_pursuit, '_', 'MarkerSize', 15, 'Color', pursuit_color_AOTU019, 'LineWidth', 2);

title('AOTU019 Angular Correlation');
xticks([1 2]); xticklabels({'Dark', 'Pursuit'});
xlim(x_limits);
ylim(y_limits);
ylabel('Correlation');
text(2.9, 0.95, sprintf('p(cond) = %.3f\np(fly) = %.3f', p_rval_019_ang(1), p_rval_019_ang(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', p_text_fontsize);

% --- AOTU019 Forward Velocity Peak Correlation ---
nexttile;
scatter(animal_conditions_019 + jitter_amount * (rand(size(rval_data_019(:,2))) - 0.5), ...
        rval_data_019(:,2), '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
median_rval_019_dark_fwd = median(rval_data_019(animal_conditions_019 == 1,2), 'omitnan');
median_rval_019_pursuit_fwd = median(rval_data_019(animal_conditions_019 == 2,2), 'omitnan');

plot(1, median_rval_019_dark_fwd, '_', 'MarkerSize', 15, 'Color', dark_color, 'LineWidth', 2);
plot(2, median_rval_019_pursuit_fwd, '_', 'MarkerSize', 15, 'Color', pursuit_color_AOTU019, 'LineWidth', 2);

title('AOTU019 Forward Correlation');
xticks([1 2]); xticklabels({'Dark', 'Pursuit'});
xlim(x_limits);
ylim(y_limits);
ylabel('Correlation');
text(2.9, 0.95, sprintf('p(cond) = %.3f\np(fly) = %.3f', p_rval_019_fwd(1), p_rval_019_fwd(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', p_text_fontsize);

% --- AOTU025 Angular Velocity Peak Correlation ---
nexttile;
scatter(animal_conditions_025 + jitter_amount * (rand(size(rval_data_025(:,1))) - 0.5), ...
        rval_data_025(:,1), '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
median_rval_025_dark = median(rval_data_025(animal_conditions_025 == 1,1), 'omitnan');
median_rval_025_pursuit = median(rval_data_025(animal_conditions_025 == 2,1), 'omitnan');

plot(1, median_rval_025_dark, '_', 'MarkerSize', 15, 'Color', dark_color, 'LineWidth', 2);
plot(2, median_rval_025_pursuit, '_', 'MarkerSize', 15, 'Color', pursuit_color_AOTU025, 'LineWidth', 2);

title('AOTU025 Angular Correlation');
xticks([1 2]); xticklabels({'Dark', 'Pursuit'});
xlim(x_limits);
ylim(y_limits);
ylabel('Correlation');
text(2.9, 0.95, sprintf('p(cond) = %.3f\np(fly) = %.3f', p_rval_025_ang(1), p_rval_025_ang(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', p_text_fontsize);

% --- AOTU025 Forward Velocity Peak Correlation ---
nexttile;
scatter(animal_conditions_025 + jitter_amount * (rand(size(rval_data_025(:,2))) - 0.5), ...
        rval_data_025(:,2), '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
median_rval_025_dark_fwd = median(rval_data_025(animal_conditions_025 == 1,2), 'omitnan');
median_rval_025_pursuit_fwd = median(rval_data_025(animal_conditions_025 == 2,2), 'omitnan');

plot(1, median_rval_025_dark_fwd, '_', 'MarkerSize', 15, 'Color', dark_color, 'LineWidth', 2);
plot(2, median_rval_025_pursuit_fwd, '_', 'MarkerSize', 15, 'Color', pursuit_color_AOTU025, 'LineWidth', 2);

title('AOTU025 Forward Correlation');
xticks([1 2]); xticklabels({'Dark', 'Pursuit'});
xlim(x_limits);
ylim(y_limits);
ylabel('Correlation');
text(2.9, 0.95, sprintf('p(cond) = %.3f\np(fly) = %.3f', p_rval_025_fwd(1), p_rval_025_fwd(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', p_text_fontsize);

% Save plots
cd(folder_plots)
saveas(gcf, 'AOTU_Comparison_dark_pursuit.png');
set(gcf, 'renderer', 'Painters'); % Save vectorized version
saveas(gcf, 'AOTU_Comparison_dark_pursuit.svg');

%% Compare Cross-Correlation Data with Medians
% Define colors for each condition
dark_color = [0 0 0];          % Black for dark condition
pursuit_color_AOTU019 = [0 0 1]; % Blue for AOTU019 pursuit condition
pursuit_color_AOTU025 = [0.5 0 0.5]; % Purple for AOTU025 pursuit condition

% Create a figure
figure; set(gcf, 'Position', [100 100 600 800])
tiledlayout(2, 2, 'TileSpacing', 'compact')

% Jitter amount for scatter plots
jitter_amount = 0.1;
y_limits = [-400 400]; % Y-axis limit for all plots
x_limits = [0 3]; % X-axis limit for all plots
p_text_fontsize = 8; % Font size for p-value text

% --- AOTU019 Angular Velocity Peak Lag ---
nexttile;
scatter(animal_conditions_019 + jitter_amount * (rand(size(lag_data_019(:,1))) - 0.5), ...
        lag_data_019(:,1), '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
median_lag_019_dark = median(lag_data_019(animal_conditions_019 == 1,1), 'omitnan');
median_lag_019_pursuit = median(lag_data_019(animal_conditions_019 == 2,1), 'omitnan');

plot(1, median_lag_019_dark, '_', 'MarkerSize', 15, 'Color', dark_color, 'LineWidth', 2);
plot(2, median_lag_019_pursuit, '_', 'MarkerSize', 15, 'Color', pursuit_color_AOTU019, 'LineWidth', 2);

title('AOTU019 Angular Lag');
xticks([1 2]); xticklabels({'Dark', 'Pursuit'}); yline(0);
xlim(x_limits);
ylim(y_limits);
ylabel('Lag');
text(2.9, 300, sprintf('p(cond) = %.3f\np(fly) = %.3f', p_lag_019_ang(1), p_lag_019_ang(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', p_text_fontsize);

% --- AOTU019 Forward Velocity Peak Lag ---
nexttile;
scatter(animal_conditions_019 + jitter_amount * (rand(size(lag_data_019(:,2))) - 0.5), ...
        lag_data_019(:,2), '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
median_lag_019_dark_fwd = median(lag_data_019(animal_conditions_019 == 1,2), 'omitnan');
median_lag_019_pursuit_fwd = median(lag_data_019(animal_conditions_019 == 2,2), 'omitnan');

plot(1, median_lag_019_dark_fwd, '_', 'MarkerSize', 15, 'Color', dark_color, 'LineWidth', 2);
plot(2, median_lag_019_pursuit_fwd, '_', 'MarkerSize', 15, 'Color', pursuit_color_AOTU019, 'LineWidth', 2);

title('AOTU019 Forward Lag');
xticks([1 2]); xticklabels({'Dark', 'Pursuit'}); yline(0);
xlim(x_limits);
ylim(y_limits);
ylabel('Lag');
text(2.9, 300, sprintf('p(cond) = %.3f\np(fly) = %.3f', p_lag_019_fwd(1), p_lag_019_fwd(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', p_text_fontsize);

% --- AOTU025 Angular Velocity Peak Lag ---
nexttile;
scatter(animal_conditions_025 + jitter_amount * (rand(size(lag_data_025(:,1))) - 0.5), ...
        lag_data_025(:,1), '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
median_lag_025_dark = median(lag_data_025(animal_conditions_025 == 1,1), 'omitnan');
median_lag_025_pursuit = median(lag_data_025(animal_conditions_025 == 2,1), 'omitnan');

plot(1, median_lag_025_dark, '_', 'MarkerSize', 15, 'Color', dark_color, 'LineWidth', 2);
plot(2, median_lag_025_pursuit, '_', 'MarkerSize', 15, 'Color', pursuit_color_AOTU025, 'LineWidth', 2);

title('AOTU025 Angular Lag');
xticks([1 2]); xticklabels({'Dark', 'Pursuit'}); yline(0);
xlim(x_limits);
ylim(y_limits);
ylabel('Lag');
text(2.9, 300, sprintf('p(cond) = %.3f\np(fly) = %.3f', p_lag_025_ang(1), p_lag_025_ang(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', p_text_fontsize);

% --- AOTU025 Forward Velocity Peak Lag ---
nexttile;
scatter(animal_conditions_025 + jitter_amount * (rand(size(lag_data_025(:,2))) - 0.5), ...
        lag_data_025(:,2), '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;
median_lag_025_dark_fwd = median(lag_data_025(animal_conditions_025 == 1,2), 'omitnan');
median_lag_025_pursuit_fwd = median(lag_data_025(animal_conditions_025 == 2,2), 'omitnan');

plot(1, median_lag_025_dark_fwd, '_', 'MarkerSize', 15, 'Color', dark_color, 'LineWidth', 2);
plot(2, median_lag_025_pursuit_fwd, '_', 'MarkerSize', 15, 'Color', pursuit_color_AOTU025, 'LineWidth', 2);

title('AOTU025 Forward Lag');
xticks([1 2]); xticklabels({'Dark', 'Pursuit'}); yline(0);
xlim(x_limits);
ylim(y_limits);
ylabel('Lag');
text(2.9, 300, sprintf('p(cond) = %.3f\np(fly) = %.3f', p_lag_025_fwd(1), p_lag_025_fwd(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', p_text_fontsize);

% Save plots
cd(folder_plots)
saveas(gcf, 'AOTU_Comparison_dark_pursuit_lag.png');
set(gcf, 'renderer', 'Painters'); % Save vectorized version
saveas(gcf, 'AOTU_Comparison_dark_pursuit_lag.svg');
