% pipeline_iinj_hold
%
% Pipeline Function
% Pulls all processed files from ALL flies in a given experiment, performs
% necessary analyses and plots accordingly.
%
% INPUTS
% exptFolder - Overarching experiment folder
%
% Created: 09/05/2023 by MC
% Updated: 07/08/2024 by MC (cleaned up and simplified)
%
function pipeline_iinj_hold(exptFolder)
%% Initialize 
disp('STARTING ANALYSES FOR POOLED IHOLD STIM...')
close all  % Close any open figures

% Set filename info and create directories
filebase = strrep(exptFolder,' ','_');  % Replace spaces in folder name with underscores
folder = generateFolders(exptFolder);  % Create necessary directories

% Load processing settings
settings = processSettings();

%% Set plot variables
% Define voltage and spike rate ranges for plotting
voltage_range = [-100 -45];  % Voltage range for first plot
spikert_range = [-80 -25];   % Spike rate range for first plot
voltage_range2 = [-85 -70];  % Voltage range for second plot
spikert_range2 = [-60 -45];  % Spike rate range for second plot

% Define limits for forward, angular, and side movement
fwd_limit = [0 6];           % Forward movement limits
ang_limit = [-160 160];      % Angular movement limits
sid_limit = [-2.5 2.5];      % Side movement limits

%% Load and pool all trials from each experiment folder
disp('Loading in datasets...')

% Set up and load all files in the directory
cd(folder.int)
allFiles = dir('*int.mat');  % Get all .mat files
nFlies = length(allFiles);   % Number of flies

% Initialize data storage arrays
nFliesThresh = nFlies;  % Threshold for flies in analysis
nt_t = 0;               % Trial counter
allForward = [];
allSideway = [];
allAngular = [];
allSpikeRt = [];
allVoltage = [];

% Loop through each fly (trial)
for nt = 1:nFlies
    % Load trial data
    thisTrial = allFiles(nt).name;
    thisFly = thisTrial(6:16);  % Extract fly name
    flyShortNames{nt} = strrep(thisFly, '_', ' ');
    cd(folder.int)
    load(thisTrial)  % Load the trial data
    nTypes = size(int_voltage, 3);  % Number of trial types

    % Pool trial data across flies
    allForward = [allForward, int_forward];
    allSideway = [allSideway, int_sideway];
    allAngular = [allAngular, int_angular];
    allSpikeRt = [allSpikeRt, int_spikert];
    allVoltage = [allVoltage, int_voltage];

    % Calculate running time for each fly
    flyRunTime(nt, 1) = (sum(int_forward > settings.runThreshE, 'all') / length(int_time)) * 60;

    % Only include flies with sufficient running time
    if flyRunTime(nt, 1) > settings.minRunTime
        nt_t = nt_t + 1;  % Increment trial counter

        % Process each trial type for this fly
        for tt = 1:nTypes
            % Extract data for this trial type
            thisForward = int_forward(:, :, tt);
            thisSideway = int_sideway(:, :, tt);
            thisAngular = int_angular(:, :, tt);

            % Process voltage data (remove spikes for certain trial types)
            if tt == 2 || tt == 4
                % Remove spike events from voltage data
                int_voltage_xspk = int_voltage;
                for n = 1:size(int_voltage, 2)
                    unclamped_spike = find(int_spikert(:, n, tt) > 0);
                    int_voltage_xspk(unclamped_spike, n, tt) = nan;
                end
                thisActivity = int_voltage_xspk(:, :, tt);  % Store voltage without spikes
            else
                thisActivity = int_voltage(:, :, tt);  % Store unprocessed voltage
            end

            % Apply median filter to smooth voltage data and remove remaining spikes
            thisActivity_mf = spikeFilter(thisActivity, int_time);

            % Analyze behavior tuning for this trial type
            thisBehave = spikert_binvelocity(thisForward, thisAngular, thisSideway, thisActivity_mf, int_time, 0);
            thisBehaveL = spikert_binvelocity(thisForward, thisAngular, thisSideway, thisActivity_mf, int_time, 1);

            % Store behavior tuning results
            all_srvfwd{nt_t, tt} = thisBehave.fwdMean';
            all_srvang{nt_t, tt} = thisBehave.angMean';
            all_srvsid{nt_t, tt} = thisBehave.sidMean';
            allL_srvfwd{nt_t, tt} = thisBehaveL.fwdMean';
            allL_srvang{nt_t, tt} = thisBehaveL.angMean';
            allL_srvsid{nt_t, tt} = thisBehaveL.sidMean';

            % Run or load cross-correlation analysis for this trial type
            cd(folder.xcorr)
            thisXCorrFile = [thisFly '_' num2str(tt) '_xc.mat'];
            if exist(thisXCorrFile, 'file')
                load(thisXCorrFile)  % Load precomputed cross-correlation data
            else
                [r_val, lag_t] = spikert_xcorr(thisActivity_mf, thisForward, thisAngular, thisSideway, int_time);
                save(thisXCorrFile, 'r_val', 'lag_t', '-v7.3');  % Save new cross-correlation data
            end

            % Use the find_peak_lag_rval function to detect peaks in cross-correlation
            [peak_lag, peak_rval, r_val] = find_peak_lag_rval(r_val, lag_t, settings.minXCorrPromVm);
            % Store results in nt_t for each motion type
            r_val_fwd(:, nt_t,tt) = r_val.fwd;
            lag_pk_fwd(nt_t,tt) = peak_lag.fwd;
            r_pk_fwd(nt_t,tt) = peak_rval.fwd;
            r_val_ang(:, nt_t,tt) = r_val.ang;
            lag_pk_ang(nt_t,tt) = peak_lag.ang;
            r_pk_ang(nt_t,tt) = peak_rval.ang;
            r_val_sid(:, nt_t,tt) = r_val.sid;
            lag_pk_sid(nt_t,tt) = peak_lag.sid;
            r_pk_sid(nt_t,tt) = peak_rval.sid;
        end
        disp(['Fly ' num2str(nt) ' ready to go!'])
    else
        disp([thisFly ' omitted from behavior analyses.'])
        nFliesThresh = nFliesThresh - 1;
    end
end

% Store velocity bin information from the last trial type processed
fwdBin = thisBehave.fwdBin;
angBin = thisBehave.angBin;
sidBin = thisBehave.sidBin;

disp('Data loaded.')


%% Plot activity versus directional velocity by trial type
cd(folder.summary)
disp('By trial type: plotting spike rate vs behavior...')

% Loop through each trial type
for tt = 1:nTypes
    % Set voltage or spike rate limits based on trial type
    if tt == 2 || tt == 4
        cell_limit = voltage_range;
        cell_limit2 = voltage_range2;
    else
        cell_limit = spikert_range;
        cell_limit2 = spikert_range2;
    end

    % Initialize figure and layout
    figure; set(gcf,'Position',[100 100 500 900])
    tiledlayout(3,2,"TileSpacing","compact")  % 3 rows, 2 columns

    % Loop through each velocity type (forward, angular, sideways)
    for v = 1:3
        switch v
            case 1  % Forward velocity
                velData = cat(2, all_srvfwd{:,tt});
                velBin = fwdBin;
            case 2  % Angular velocity
                velData = cat(2, all_srvang{:,tt});
                velBin = angBin;
            case 3  % Sideways velocity
                velData = cat(2, all_srvsid{:,tt});
                velBin = sidBin;
        end
        
        % Calculate mean and SEM of velocity data
        velMean = mean(velData, 2, 'omitnan');
        velSEM = std(velData, 0, 2, 'omitnan') / sqrt(nFliesThresh);

        % Plot individual trials and mean
        nexttile; hold on
        plot(velBin, velData, 'Color', settings.trialColor, 'LineWidth', settings.lwTri)
        plot(velBin, velMean, 'Color', settings.velColor{v}, 'LineWidth', settings.lwAvg)
        axis tight; ylim(cell_limit);
        ylabel('Voltage (mV)'); xlabel(settings.velLabel{v})

        % Plot mean with SEM as shaded area
        nexttile; hold on
        sem = patch([velBin'; flipud(velBin')], [velMean - velSEM; flipud(velMean + velSEM)], ...
                    'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
        sem.FaceColor = settings.velColor{v};  % Set color of SEM patch
        plot(velBin, velMean, 'Color', settings.velColor{v}, 'LineWidth', settings.lwAvg)
        axis tight; ylim(cell_limit2);
        ylabel('Voltage (mV)'); xlabel(settings.velLabel{v})
    end
    
    % Add title with trial information
    thisletter = settings.letters(tt);
    sgtitle([strrep(filebase,'_',' ') ' ' settings.holdLabel{tt} ' (n = ' num2str(nFliesThresh) ')'])

    % Save plot as PNG and SVG
    plotname = ['ihold_' thisletter '_sp_v_vel'];
    cd(folder.summary)
    saveas(gcf, [plotname '.png']);  % Save PNG
    copyfile([plotname '.png'], folder.dropbox, 'f');  % Copy to Dropbox

    cd(folder.vector)
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg']);  % Save vector plot as SVG
    copyfile([plotname '.svg'], folder.dropbox, 'f');  % Copy to Dropbox
end


%% Plot activity versus directional velocity by trial type (with lag)
cd(folder.summary)
disp('By trial type: plotting spike rate vs behavior (with lag)...')

% Loop through each trial type
for tt = 1:nTypes
    % Set voltage or spike rate limits based on trial type
    if tt == 2 || tt == 4
        cell_limit = voltage_range;
        cell_limit2 = voltage_range2;
    else
        cell_limit = spikert_range;
        cell_limit2 = spikert_range2;
    end

    % Initialize figure and layout
    figure; set(gcf,'Position',[100 100 500 900])
    tiledlayout(3,2,"TileSpacing","compact")  % 3 rows, 2 columns

    % Loop through each velocity type (forward, angular, sideways)
    for v = 1:3
        switch v
            case 1  % Forward velocity
                velData = cat(2, allL_srvfwd{:,tt});
                velBin = fwdBin;
            case 2  % Angular velocity
                velData = cat(2, allL_srvang{:,tt});
                velBin = angBin;
            case 3  % Sideways velocity
                velData = cat(2, allL_srvsid{:,tt});
                velBin = sidBin;
        end
        
        % Calculate mean and SEM for velocity data
        velMean = mean(velData, 2, 'omitnan');
        velSEM = std(velData, 0, 2, 'omitnan') / sqrt(nFliesThresh);

        % Plot individual trials and mean
        nexttile; hold on
        plot(velBin, velData, 'Color', settings.trialColor, 'LineWidth', settings.lwTri)
        plot(velBin, velMean, 'Color', settings.velColor{v}, 'LineWidth', settings.lwAvg)
        axis tight; ylim(cell_limit);
        ylabel('Voltage (mV)'); xlabel(settings.velLabel{v})

        % Plot mean with SEM as shaded area
        nexttile; hold on
        sem = patch([velBin'; flipud(velBin')], [velMean - velSEM; flipud(velMean + velSEM)], ...
                    'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
        sem.FaceColor = settings.velColor{v};  % Set color of SEM patch
        plot(velBin, velMean, 'Color', settings.velColor{v}, 'LineWidth', settings.lwAvg)
        axis tight; ylim(cell_limit2);
        ylabel('Voltage (mV)'); xlabel(settings.velLabel{v})
    end
    
    % Add title with trial and lag information
    thisletter = settings.letters(tt);
    sgtitle([strrep(filebase, '_', ' ') ' ' settings.holdLabel{tt} ' w/lag (n = ' num2str(nFliesThresh) ')'])

    % Save plot as PNG and SVG
    plotname = ['ihold_' thisletter '_sp_v_velL'];
    cd(folder.summary)
    saveas(gcf, [plotname '.png']);  % Save PNG
    copyfile([plotname '.png'], folder.dropbox, 'f');  % Copy to Dropbox

    cd(folder.vector)
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg']);  % Save vector plot as SVG
    copyfile([plotname '.svg'], folder.dropbox, 'f');  % Copy to Dropbox
end

%% Plot lag and cross-correlation by trial type
disp('Performing cross-correlation for firing vs pursuit behavior across conditions...')
cd(folder.summary)

% Calculate median r-values and lags across flies, ignoring NaNs
r_val_fwd_mean = median(r_val_fwd, 2, 'omitnan');
r_val_ang_mean = median(r_val_ang, 2, 'omitnan');
r_val_sid_mean = median(r_val_sid, 2, 'omitnan');
lagpeak_fwd_mean = median(lag_pk_fwd, 1, 'omitnan');
lagpeak_ang_mean = median(lag_pk_ang, 1, 'omitnan');
lagpeak_sid_mean = median(lag_pk_sid, 1, 'omitnan');
rpeak_fwd_mean = median(r_pk_fwd, 1, 'omitnan');
rpeak_ang_mean = median(r_pk_ang, 1, 'omitnan');
rpeak_sid_mean = median(r_pk_sid, 1, 'omitnan');

% Initialize figure and plot settings
figure;
set(gcf, 'Position', [100, 100, 1200, 800])
xc_lim = [-400, 400];  % x-axis limits for lag plots
x = 1:nTypes;  % Condition indices for summary plots
suby = nTypes + 2;  % Number of columns for subplots
r_range = [-0.5, 1];  % y-axis limits for r-values

% Plot r-values over time lags for each condition
for tt = 1:nTypes
    % Forward velocity cross-correlation plot
    subplot(3, suby, tt)
    plot(lag_t, r_val_fwd(:,:,tt), '-.', 'LineWidth', settings.lwTri, 'Color', settings.trialColor); hold on
    plot(lag_t, r_val_fwd_mean(:, :, tt), 'LineWidth', settings.lwAvg, 'Color', settings.velColor{1})
    title(settings.holdLabel{tt})
    ylim(r_range)
    xline(0);
    if tt == 1
        ylabel('Forward r-value')
    end

    % Angular velocity cross-correlation plot
    subplot(3, suby, suby + tt)
    plot(lag_t, r_val_ang(:,:,tt), '-.', 'LineWidth', settings.lwTri, 'Color', settings.trialColor); hold on
    plot(lag_t, r_val_ang_mean(:, :, tt), 'LineWidth', settings.lwAvg, 'Color', settings.velColor{2})
    ylim(r_range)
    xline(0);
    if tt == 1
        ylabel('Angular r-value')
    end

    % Sideways velocity cross-correlation plot
    subplot(3, suby, 2 * suby + tt)
    plot(lag_t, r_val_sid(:,:,tt), '-.', 'LineWidth', settings.lwTri, 'Color', settings.trialColor); hold on
    plot(lag_t, r_val_sid_mean(:, :, tt), 'LineWidth', settings.lwAvg, 'Color', settings.velColor{3})
    ylim(r_range)
    xline(0);
    if tt == 1
        ylabel('Sideway r-value')
    end
    xlabel('Time (ms)')
end

% Summary plots of r-values and optimal lags for forward, angular, and sideway velocities

% Forward velocity summary plot
subplot(3, suby, suby - 1)
plot(x, r_pk_fwd, '.', 'Color', settings.trialColor); hold on
plot(x, rpeak_fwd_mean, 'Marker', '_', 'LineStyle', 'none', 'Color', settings.velColor{1});
xlim([0, nTypes + 1]);
xticks(x)
xticklabels([]);  % No x-tick labels
ylim(r_range);
ylabel('Peak r')
yline(0);

% Lag for forward velocity
subplot(3, suby, suby)
plot(x, lag_pk_fwd, '.', 'Color', settings.trialColor); hold on
plot(x, lagpeak_fwd_mean, 'Marker', '_', 'LineStyle', 'none', 'Color', settings.velColor{1});
xlim([0, nTypes + 1]);
xticks(x)
xticklabels([]);  % No x-tick labels
ylim(xc_lim);
ylabel('Lag (ms)')
yline(0);

% Label optimal lag points for forward velocity
for i = x
    text(x(i), 250, num2str(lagpeak_fwd_mean(i)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 8, 'Rotation', 90);
end

% Angular velocity summary plot
subplot(3, suby, 2 * suby - 1)
plot(x, r_pk_ang, '.', 'Color', settings.trialColor); hold on
plot(x, rpeak_ang_mean, 'Marker', '_', 'LineStyle', 'none', 'Color', settings.velColor{2});
xlim([0, nTypes + 1]);
xticks(x)
xticklabels([]);  % No x-tick labels
ylim(r_range);
ylabel('Peak r')
yline(0);

% Lag for angular velocity
subplot(3, suby, 2 * suby)
plot(x, lag_pk_ang, '.', 'Color', settings.trialColor); hold on
plot(x, lagpeak_ang_mean, 'Marker', '_', 'LineStyle', 'none', 'Color', settings.velColor{2});
xlim([0, nTypes + 1]);
xticks(x)
xticklabels([]);  % No x-tick labels
ylim(xc_lim);
ylabel('Lag (ms)')
yline(0);

% Label optimal lag points for angular velocity
for i = x
    text(x(i), 250, num2str(lagpeak_ang_mean(i)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 8, 'Rotation', 90);
end

% Sideway velocity summary plot
subplot(3, suby, 3 * suby - 1)
plot(x, r_pk_sid, '.', 'Color', settings.trialColor); hold on
plot(x, rpeak_sid_mean, 'Marker', '_', 'LineStyle', 'none', 'Color', settings.velColor{3});
xlim([0, nTypes + 1]);
xticks(x)
xticklabels(settings.holdLabel);
ylim(r_range);
ylabel('Peak r')
yline(0);

% Lag for sideway velocity
subplot(3, suby, 3 * suby)
plot(x, lag_pk_sid, '.', 'Color', settings.trialColor); hold on
plot(x, lagpeak_sid_mean, 'Marker', '_', 'LineStyle', 'none', 'Color', settings.velColor{3});
xlim([0, nTypes + 1]);
xticks(x)
xticklabels(settings.holdLabel);
ylim(xc_lim);
ylabel('Lag (ms)')
yline(0);

% Save plot as PNG and SVG
cd(folder.summary)
plotname = 'ihold_xcorr';
saveas(gcf, [plotname '.png']);  % Save as PNG
copyfile([plotname '.png'], folder.dropbox, 'f');  % Copy to Dropbox

cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg']);  % Save as SVG
copyfile([plotname '.svg'], folder.dropbox, 'f');  % Copy to Dropbox


%% End
disp('ALL ANALYSES COMPLETE.')
end

