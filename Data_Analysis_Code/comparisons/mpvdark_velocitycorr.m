% mpvdark_velocitycorr (updated version without -P1 darkness)

%% Load cross correlation data
clear; close all

plotPath = 'E:\\Compare Motion Pulse';
dataPath = 'E:\\Compare Motion Pulse\\data';

disp('Analyzing xcorr in darkness and visual pursuit.');
cd(dataPath)

% Load files
AOTU019_dark = load(fullfile(dataPath, 'AOTU019_Background_P1_xcorr.mat'));
AOTU025_dark = load(fullfile(dataPath, 'AOTU025_Background_P1_xcorr.mat'));
AOTU019_visual = load(fullfile(dataPath, 'AOTU019_Motion_Pulse_xcorr.mat'));
AOTU025_visual = load(fullfile(dataPath, 'AOTU025_Motion_Pulse_xcorr.mat'));

% Data
AOTU019_angr_dark_p1 = AOTU019_dark.combined_data.r_pk_ang;
AOTU025_angr_dark_p1 = AOTU025_dark.combined_data.r_pk_ang;
AOTU019_angr_dark_nop1 = AOTU019_dark.combined_data.r_pk_ang_nop1;
AOTU025_angr_dark_nop1 = AOTU025_dark.combined_data.r_pk_ang_nop1;
AOTU019_angr_visual_p1 = AOTU019_visual.rotational_peaks.r_pk;
AOTU025_angr_visual_p1 = AOTU025_visual.rotational_peaks.r_pk;

AOTU019_anglag_dark_p1 = AOTU019_dark.combined_data.lag_pk_ang;
AOTU025_anglag_dark_p1 = AOTU025_dark.combined_data.lag_pk_ang;
AOTU019_anglag_dark_nop1 = AOTU019_dark.combined_data.lag_pk_ang_nop1;
AOTU025_anglag_dark_nop1 = AOTU025_dark.combined_data.lag_pk_ang_nop1;
AOTU019_anglag_visual_p1 = AOTU019_visual.rotational_peaks.lag_pk;
AOTU025_anglag_visual_p1 = AOTU025_visual.rotational_peaks.lag_pk;

%% Organize cross correlation data
% --- AOTU019 ---
names_dark_p1_019    = AOTU019_dark.storeNames(:);
names_dark_nop1_019  = AOTU019_dark.storeNames(:);
names_visual_p1_019  = AOTU019_visual.storeNames(:);

all_animals_019 = unique([names_dark_p1_019; names_dark_nop1_019; names_visual_p1_019]);

r_data_019   = nan(numel(all_animals_019), 3);
lag_data_019 = nan(numel(all_animals_019), 3);

% Fill correlation
r_data_019 = fill_column(r_data_019, all_animals_019, names_dark_p1_019,   AOTU019_angr_dark_p1,    1);
r_data_019 = fill_column(r_data_019, all_animals_019, names_dark_nop1_019, AOTU019_angr_dark_nop1,  2);
r_data_019 = fill_column(r_data_019, all_animals_019, names_visual_p1_019, AOTU019_angr_visual_p1,  3);

% Fill lag
lag_data_019 = fill_column(lag_data_019, all_animals_019, names_dark_p1_019,   AOTU019_anglag_dark_p1,    1);
lag_data_019 = fill_column(lag_data_019, all_animals_019, names_dark_nop1_019, AOTU019_anglag_dark_nop1,  2);
lag_data_019 = fill_column(lag_data_019, all_animals_019, names_visual_p1_019, AOTU019_anglag_visual_p1,  3);
r_data_019(lag_data_019>0) = nan; %bad xcorr
lag_data_019(lag_data_019>0) = nan; %bad xcorr

r_table_019   = array2table(r_data_019,   'VariableNames', {'Dark_P1','Dark_noP1','Visual_P1'}, 'RowNames', all_animals_019);
lag_table_019 = array2table(lag_data_019, 'VariableNames', {'Dark_P1','Dark_noP1','Visual_P1'}, 'RowNames', all_animals_019);

% --- AOTU025 ---
names_dark_p1_025    = AOTU025_dark.storeNames(:);
names_dark_nop1_025  = AOTU025_dark.storeNames(:);
names_visual_p1_025  = AOTU025_visual.storeNames(:);

all_animals_025 = unique([names_dark_p1_025; names_dark_nop1_025; names_visual_p1_025]);

r_data_025   = nan(numel(all_animals_025), 3);
lag_data_025 = nan(numel(all_animals_025), 3);

% Fill correlation
r_data_025 = fill_column(r_data_025, all_animals_025, names_dark_p1_025,   AOTU025_angr_dark_p1,    1);
r_data_025 = fill_column(r_data_025, all_animals_025, names_dark_nop1_025, AOTU025_angr_dark_nop1,  2);
r_data_025 = fill_column(r_data_025, all_animals_025, names_visual_p1_025, AOTU025_angr_visual_p1,  3);

% Fill lag
lag_data_025 = fill_column(lag_data_025, all_animals_025, names_dark_p1_025,   AOTU025_anglag_dark_p1,    1);
lag_data_025 = fill_column(lag_data_025, all_animals_025, names_dark_nop1_025, AOTU025_anglag_dark_nop1,  2);
lag_data_025 = fill_column(lag_data_025, all_animals_025, names_visual_p1_025, AOTU025_anglag_visual_p1,  3);
r_data_025(lag_data_025>0) = nan; %bad xcorr
lag_data_025(lag_data_025>0) = nan; %bad xcorr

r_table_025   = array2table(r_data_025,   'VariableNames', {'Dark_P1','Dark_noP1','Visual_P1'}, 'RowNames', all_animals_025);
lag_table_025 = array2table(lag_data_025, 'VariableNames', {'Dark_P1','Dark_noP1','Visual_P1'}, 'RowNames', all_animals_025);

% Number of animals (non-NaN entries) per condition
n_lag_019 = sum(~isnan(lag_data_019), 1);
n_lag_025 = sum(~isnan(lag_data_025), 1);

% Display results
disp('AOTU019 lag data counts per condition (Dark+P1, Dark−P1, Visual+P1):');
disp(n_lag_019);

disp('AOTU025 lag data counts per condition (Dark+P1, Dark−P1, Visual+P1):');
disp(n_lag_025);

%% Plot
close all
condLabels = {'Dark+arousal','Dark−arousal','Visual+arousal'};

% Compute consistent y-limits per metric across both cell types
r_all   = [table2array(r_table_019); table2array(r_table_025)];
lag_all = [table2array(lag_table_019); table2array(lag_table_025)];
r_ylim   = [min(r_all,[],'all','omitnan'),   max(r_all,[],'all','omitnan')];
lag_ylim = [min(lag_all,[],'all','omitnan'), 0];

figure('Color','w'); set(gcf, 'Position', [100 100 600 800])
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% Helper to plot one metric/table
plot_tile = @(tbl, ttl, ylab, ylim_use) ...
    local_plot_by_condition(tbl, condLabels, ttl, ylab, ylim_use);

% Row 1: r
nexttile; plot_tile(r_table_019, 'AOTU019 — correlation (r)', 'r', r_ylim);
nexttile; plot_tile(r_table_025, 'AOTU025 — correlation (r)', 'r', r_ylim);

% Row 2: lag
nexttile; plot_tile(lag_table_019, 'AOTU019 — lag (samples)', 'lag', lag_ylim);
nexttile; plot_tile(lag_table_025, 'AOTU025 — lag (samples)', 'lag', lag_ylim);

xlabel(tl, 'Condition');

% Save plot
cd(plotPath)
plotname = 'xcorr_dark_v_mp';
saveas(gcf, [plotname '.png']);

% Save vectorized plot
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
disp('Complete.')

%% Stats
% ---------- Helper: reshape wide 3-col table -> long table with factors ----------
% Returns columns: Animal (categorical), Value (double), Condition (categorical),
% Visual (0=Dark, 1=Motion), Arousal (0=NoP1, 1=P1)
make_long = @(T) local_make_long(T);

% Build long tables
r019   = make_long(r_table_019);
lag019 = make_long(lag_table_019);
r025   = make_long(r_table_025);
lag025 = make_long(lag_table_025);

% ---------- Fit LME and run tests ----------
run_lme_and_tests(r019,   'AOTU019 — r');
run_lme_and_tests(lag019, 'AOTU019 — lag');
run_lme_and_tests(r025,   'AOTU025 — r');
run_lme_and_tests(lag025, 'AOTU025 — lag');