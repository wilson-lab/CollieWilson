% pipeline_p1background
%
% Pipeline Function
% Pulls all processed files from ALL flies in a given experiment, performs
% necessary analyses and plots accordingly. Can be used for both behavior-only
% experiments and ephys experiments.
%
% INPUTS
% exptFolder - overarching experiment folder
%
% This function pools data across all flies, performs cross-correlation
% analyses, spike rate vs velocity analyses, and voltage vs firing rate
% analyses. The results are saved as plots, and both raster and summary
% figures are generated for P1 stimulation background experiments.
%
% 06/21/2022 - MC adapted from visual pursuit pipeline
% 07/08/2024 - MC cleaned up and simplified
% 11/07/2024 - MC expanded and improved xcorr code
%
function pipeline_p1background(exptFolder)
%% Initialize
disp('STARTING ANALYSES FOR POOLED P1 STIM BACKGROUND ONLY...')
close all  % Close any open figures to start fresh for new analyses

% Load processing settings
settings = processSettings();

% Set filename info and create necessary directories
filebase = strrep(exptFolder, ' ', '_');

% Generate folder structure
folder = generateFolders(exptFolder);

% Define condition info
trialTypes = {'P1 Stim'; 'No Stim'};
nTypes = 2;  % Specifies the number of conditions for clarity

%% Set plotting variables

% Define plotting ranges for various data types
srRange = 80;         % Range for spike rate (spikes/s)
fwdRange = 10;        % Range for forward velocity (mm/s)
angRange = 250;       % Range for angular velocity (deg/s)
sidRange = 2.5;       % Range for sideways velocity (mm/s)

% Range for correlation coefficient (r) in plots
r_range = [-0.2 0.8];

%% Load in and pool all pulse trials from each experiment
disp('Loading in pulse datasets...')

% Find all files in this directory
cd(folder.int)
pulseFiles = dir('*pulse.mat');  % Pull all file info for pulse datasets
nFlies = length(pulseFiles);     % Total number of flies based on file count
flylist = [];

% Initialize data storage arrays
nFliesThresh = nFlies;  % Initialize counter for flies in behavior analysis
nt_t = 0;               % Counter for trials
binFRVm_on = []; binFRVm_off = [];
binProbVm_on = []; binProbVm_off = [];
storeNames = {};

for nt = 1:nFlies
    % Load this trial file
    thisTrial = pulseFiles(nt).name;
    thisFly = extractBefore(thisTrial, '_pulse');  % Extract fly ID
    disp(['Processing: ' thisFly])
    cd(folder.int)
    load(thisTrial)

    % Calculate distribution of firing rates
    [fr_bins, fr_counts_on(nt,:)] = compute_firingrate_distribution(int_spikert(:,:,1));
    [fr_bins, fr_counts_off(nt,:)] = compute_firingrate_distribution(int_spikert(:,:,2));

    % Threshold for flies with sufficient running time
    flyRunTime(nt,1) = (sum(int_forward > settings.runThreshE, 'all') / length(int_time)) * 60;

    % Bin spike rate by voltage levels
    thisVoltage = spikeFilter(int_voltage(:,:,1), int_time);
    sr_binOn = spikert_v_voltage(int_spikert(:,:,1), thisVoltage);
    thisVoltage = spikeFilter(int_voltage(:,:,2), int_time);
    sr_binOff = spikert_v_voltage(int_spikert(:,:,2), thisVoltage);
    binFRVm_on(:,nt) = sr_binOn.fr';
    binProbVm_on(:,nt) = sr_binOn.prob';
    binFRVm_off(:,nt) = sr_binOff.fr';
    binProbVm_off(:,nt) = sr_binOff.prob';

    % Determine average change in activity with/without P1
    [diffSR, diffVm] = analyzeSpikingAndVoltageDifference(int_spikert, int_voltage, int_time, int_forward,settings);
    diff_sr(nt,:) = diffSR;
    diff_vm(nt,:) = diffVm;

    if flyRunTime(nt,1) > settings.minRunTime
        % Update trial counter
        nt_t = nt_t + 1;
        storeNames{nt_t} = thisFly;

        % Run cross-correlation analysis for P1 stim and no stim conditions
        for tt = 1:nTypes
            % Retrieve data for cross-correlation analysis
            thisSpikeRt = int_spikert(:,:,tt);
            thisVoltage = spikeFilter(int_voltage(:,:,tt), int_time);
            thisForward = int_forward(:,:,tt);
            thisAngular = int_angular(:,:,tt);
            thisSideway = int_sideway(:,:,tt);

            % Bin spike rate by velocity without and with lag estimate
            [binVel_noLag] = spikert_binvelocity(thisForward, thisAngular, thisSideway, thisSpikeRt, int_time, 0);
            [binVel_withLag] = spikert_binvelocity(thisForward, thisAngular, thisSideway, thisSpikeRt, int_time, 1);
            % Bin voltage by velocity without and with lag estimate
            [binVolt_noLag] = spikert_binvelocity(thisForward, thisAngular, thisSideway, thisVoltage, int_time, 0);
            [binVolt_withLag] = spikert_binvelocity(thisForward, thisAngular, thisSideway, thisVoltage, int_time, 1);
            % Bin spikerate by velocity with lag estimate
            [binSR_withLag] = velocity_binbyspikerate(thisForward, thisAngular, thisSideway, thisSpikeRt, int_time, 1);

            % Fit spiker rate to velocity with lag estimate
            [slope_fwd, slope_ang, slope_sid] = fit_velocity_data(thisForward, thisAngular, thisSideway, thisSpikeRt, int_time, 1);

            % Determine running stats
            runIdx = schmittTrigger(thisForward,settings.runThreshE,0.1);
            runPercent = median(sum(runIdx)/size(runIdx,1)*100);

            % Store spike rate data according to condition
            if tt == 1
                srfwd_on_nolag(:,nt_t) = binVel_noLag.fwdMean';
                srang_on_nolag(:,nt_t) = binVel_noLag.angMean';
                srsid_on_nolag(:,nt_t) = binVel_noLag.sidMean';
                srfwd_on_lag(:,nt_t) = binVel_withLag.fwdMean';
                srang_on_lag(:,nt_t) = binVel_withLag.angMean';
                srsid_on_lag(:,nt_t) = binVel_withLag.sidMean';

                fwdsr_on_lag(:,nt_t) = binSR_withLag.fwdMean';
                angsr_on_lag(:,nt_t) = binSR_withLag.angMean';
                sidsr_on_lag(:,nt_t) = binSR_withLag.sidMean';

                % Store voltage data for 'on' condition
                vfw_on_nolag(:,nt_t) = binVolt_noLag.fwdMean';
                vang_on_nolag(:,nt_t) = binVolt_noLag.angMean';
                vsid_on_nolag(:,nt_t) = binVolt_noLag.sidMean';
                vfw_on_lag(:,nt_t) = binVolt_withLag.fwdMean';
                vang_on_lag(:,nt_t) = binVolt_withLag.angMean';
                vsid_on_lag(:,nt_t) = binVolt_withLag.sidMean';

                % Store fit for 'on' condition
                fwd_on_fits(nt_t) = slope_fwd;
                ang_on_fits(nt_t) = slope_ang;
                sid_on_fits(nt_t) = slope_sid;

                % Store run time for 'on' condition
                runpercent_on(nt_t) = runPercent;

            else
                srfwd_off_nolag(:,nt_t) = binVel_noLag.fwdMean';
                srang_off_nolag(:,nt_t) = binVel_noLag.angMean';
                srsid_off_nolag(:,nt_t) = binVel_noLag.sidMean';
                srfwd_off_lag(:,nt_t) = binVel_withLag.fwdMean';
                srang_off_lag(:,nt_t) = binVel_withLag.angMean';
                srsid_off_lag(:,nt_t) = binVel_withLag.sidMean';

                fwdsr_off_lag(:,nt_t) = binSR_withLag.fwdMean';
                angsr_off_lag(:,nt_t) = binSR_withLag.angMean';
                sidsr_off_lag(:,nt_t) = binSR_withLag.sidMean';

                % Store voltage data for 'off' condition
                vfw_off_nolag(:,nt_t) = binVolt_noLag.fwdMean';
                vang_off_nolag(:,nt_t) = binVolt_noLag.angMean';
                vsid_off_nolag(:,nt_t) = binVolt_noLag.sidMean';
                vfw_off_lag(:,nt_t) = binVolt_withLag.fwdMean';
                vang_off_lag(:,nt_t) = binVolt_withLag.angMean';
                vsid_off_lag(:,nt_t) = binVolt_withLag.sidMean';

                % Store fit for 'off' condition
                fwd_off_fits(nt_t) = slope_fwd;
                ang_off_fits(nt_t) = slope_ang;
                sid_off_fits(nt_t) = slope_sid;

                % Store run time for 'off' condition
                runpercent_off(nt_t) = runPercent;
            end


            % Check if cross-correlation file exists
            cd(folder.xcorr)
            thisXCorrFile = [thisFly '_' num2str(tt) '_xc.mat'];
            if exist(thisXCorrFile, 'file')
                disp('Loading previous xcorr.')
                load(thisXCorrFile)
            else
                % Run cross-correlation if no file found
                [r_val, lag_t] = spikert_xcorr(thisSpikeRt, thisForward, thisAngular, thisSideway, int_time);
                save(thisXCorrFile, 'r_val', 'lag_t', '-v7.3');  % Save results
            end

            % Find peaks for each motion type using the function
            [peak_lag, peak_rval,r_val] = find_peak_lag_rval(r_val, lag_t, settings.minXCorrProm);

            % Store peak results for each trial type
            r_val_fwd(:, nt_t, tt) = r_val.fwd;
            lag_pk_fwd(nt_t, tt) = peak_lag.fwd;
            r_pk_fwd(nt_t, tt) = peak_rval.fwd;
            r_val_ang(:, nt_t, tt) = r_val.ang;
            lag_pk_ang(nt_t, tt) = peak_lag.ang;
            r_pk_ang(nt_t, tt) = peak_rval.ang;
            r_val_sid(:, nt_t, tt) = r_val.sid;
            lag_pk_sid(nt_t, tt) = peak_lag.sid;
            r_pk_sid(nt_t, tt) = peak_rval.sid;

        end
        flylist{nt_t} = thisFly;
    else
        disp([thisFly ' omitted from behavior analyses.'])
        nFliesThresh = nFliesThresh - 1;
    end
end

disp('Complete.')

% Store bins for future use
if nt_t>0
    fwdBins = binVel_noLag.fwdBin';
    angBins = binVel_noLag.angBin';
    sidBins = binVel_noLag.sidBin';
    vmBins = sr_binOn.bin';

    srBins = binSR_withLag.spikeRateBin';
end

%% Cross-Correlation Analysis for Firing Rate vs. Pursuit Behavior
if nt_t>1
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
    set(gcf, 'Position', [100, 100, 600, 800])
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
        title(trialTypes{tt})
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
    xticklabels([]);
    ylim(r_range);
    ylabel('Peak r')
    yline(0);

    % Lag for forward velocity
    subplot(3, suby, suby)
    plot(x, lag_pk_fwd, '.', 'Color', settings.trialColor); hold on
    plot(x, lagpeak_fwd_mean, 'Marker', '_', 'LineStyle', 'none', 'Color', settings.velColor{1});
    xlim([0, nTypes + 1]);
    xticks(x)
    xticklabels([]);
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
    xticklabels([]);
    ylim(r_range);
    ylabel('Peak r')
    yline(0);

    % Lag for angular velocity
    subplot(3, suby, 2 * suby)
    plot(x, lag_pk_ang, '.', 'Color', settings.trialColor); hold on
    plot(x, lagpeak_ang_mean, 'Marker', '_', 'LineStyle', 'none', 'Color', settings.velColor{2});
    xlim([0, nTypes + 1]);
    xticks(x)
    xticklabels([]);
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
    xticklabels(trialTypes);
    ylim(r_range);
    ylabel('Peak r')
    yline(0);

    % Lag for sideway velocity
    subplot(3, suby, 3 * suby)
    plot(x, lag_pk_sid, '.', 'Color', settings.trialColor); hold on
    plot(x, lagpeak_sid_mean, 'Marker', '_', 'LineStyle', 'none', 'Color', settings.velColor{3});
    xlim([0, nTypes + 1]);
    xticks(x)
    xticklabels(trialTypes);
    ylim(xc_lim);
    ylabel('Lag (ms)')
    yline(0);

    % Label optimal lag points for sideway velocity
    for i = x
        text(x(i), 250, num2str(lagpeak_sid_mean(i)), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', 8, 'Rotation', 90);
    end

    % Add title and save plot
    sgtitle([strrep(filebase, '_', ' ') ' (n = ' num2str(nFliesThresh) ') Xcorr FR'])

    % Save plot as PNG and SVG
    cd(folder.summary)
    plotname = 'p1background_xcorr';
    saveas(gcf, [plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox, 'f');

    cd(folder.vector)
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg']);
    copyfile([plotname '.svg'], folder.dropbox, 'f');

    % Save key cross-correlation variables in a .mat file
    cd(folder.compare);
    combined_data.r_pk_fwd = r_pk_fwd(:, 1);
    combined_data.r_pk_ang = r_pk_ang(:, 1);
    combined_data.lag_pk_fwd = lag_pk_fwd(:, 1);
    combined_data.lag_pk_ang = lag_pk_ang(:, 1);
    combined_data.r_pk_fwd_nop1 = r_pk_fwd(:, 2);
    combined_data.r_pk_ang_nop1 = r_pk_ang(:, 2);
    combined_data.lag_pk_fwd_nop1 = lag_pk_fwd(:, 2);
    combined_data.lag_pk_ang_nop1 = lag_pk_ang(:, 2);

    filename = [filebase '_xcorr.mat'];
    save(filename, 'combined_data', 'storeNames');

    disp('Cross-correlation analysis complete.')
end

%% plot spike rate versus directional velocity for w/ w/o P1 activation
if nt_t>1
    disp('Analyzing behavior tuning for P1 on vs off...')

    % replace missing bins with nans
    srfwd_on_nolag(srfwd_on_nolag==0) = nan;
    srang_on_nolag(srang_on_nolag==0) = nan;
    srsid_on_nolag(srsid_on_nolag==0) = nan;
    srfwd_off_nolag(srfwd_off_nolag==0) = nan;
    srang_off_nolag(srang_off_nolag==0) = nan;
    srsid_off_nolag(srsid_off_nolag==0) = nan;

    % calculate means
    mean_srfwd_on = mean(srfwd_on_nolag,2,'omitnan');
    mean_srang_on = mean(srang_on_nolag,2,'omitnan');
    mean_srsid_on = mean(srsid_on_nolag,2,'omitnan');
    mean_srfwd_off = mean(srfwd_off_nolag,2,'omitnan');
    mean_srang_off = mean(srang_off_nolag,2,'omitnan');
    mean_srsid_off = mean(srsid_off_nolag,2,'omitnan');
    % calculate SEMs
    sem_srfwd_on = std(srfwd_on_nolag,0,2,'omitnan')/sqrt(nFliesThresh);
    sem_srang_on = std(srang_on_nolag,0,2,'omitnan')/sqrt(nFliesThresh);
    sem_srsid_on = std(srsid_on_nolag,0,2,'omitnan')/sqrt(nFliesThresh);
    sem_srfwd_off = std(srfwd_off_nolag,0,2,'omitnan')/sqrt(nFliesThresh);
    sem_srang_off = std(srang_off_nolag,0,2,'omitnan')/sqrt(nFliesThresh);
    sem_srsid_off = std(srsid_off_nolag,0,2,'omitnan')/sqrt(nFliesThresh);


    % initialize
    figure; set(gcf,'Position',[100 100 1000 400])

    % plot forward velocity
    subplot(1,3,1)
    % plot SEM band
    r(1) = patch([fwdBins; flipud(fwdBins)],[(mean_srfwd_off-sem_srfwd_off); flipud((mean_srfwd_off+sem_srfwd_off))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.bckColor{1});
    hold on
    r(2) = patch([fwdBins; flipud(fwdBins)],[(mean_srfwd_on-sem_srfwd_on); flipud((mean_srfwd_on+sem_srfwd_on))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.bckColor{2});
    % plot average
    plot(fwdBins,mean_srfwd_off,'Color',settings.bckColor{1},'MarkerFaceColor','w','LineWidth',settings.lwAvg)
    plot(fwdBins,mean_srfwd_on,'Color',settings.bckColor{2},'MarkerFaceColor','w','LineWidth',settings.lwAvg)
    ylabel(settings.spkLabel)
    xlabel(settings.velLabel{1})
    ylim([0 srRange])
    xlim([0 fwdRange])

    % plot angular
    subplot(1,3,2)
    % plot SEM band
    r(1) = patch([angBins; flipud(angBins)],[(mean_srang_off-sem_srang_off); flipud((mean_srang_off+sem_srang_off))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.bckColor{1});
    hold on
    r(2) = patch([angBins; flipud(angBins)],[(mean_srang_on-sem_srang_on); flipud((mean_srang_on+sem_srang_on))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.bckColor{2});
    % plot average
    plot(angBins,mean_srang_off,'Color',settings.bckColor{1},'MarkerFaceColor','w','LineWidth',settings.lwAvg)
    plot(angBins,mean_srang_on,'Color',settings.bckColor{2},'MarkerFaceColor','w','LineWidth',settings.lwAvg)
    xlabel(settings.velLabel{2})
    ylim([0 srRange])
    xlim([-angRange angRange])
    xline(0)

    % plot sideway
    subplot(1,3,3)
    % plot SEM band
    r(1) = patch([sidBins; flipud(sidBins)],[(mean_srsid_off-sem_srsid_off); flipud((mean_srsid_off+sem_srsid_off))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.bckColor{1});
    hold on
    r(2) = patch([sidBins; flipud(sidBins)],[(mean_srsid_on-sem_srsid_on); flipud((mean_srsid_on+sem_srsid_on))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.bckColor{2});
    % plot average
    plot(sidBins,mean_srsid_off,'Color',settings.bckColor{1},'MarkerFaceColor','w','LineWidth',settings.lwAvg)
    plot(sidBins,mean_srsid_on,'Color',settings.bckColor{2},'MarkerFaceColor','w','LineWidth',settings.lwAvg)
    ylabel(settings.spkLabel)
    xlabel(settings.velLabel{3})
    ylim([0 srRange])
    xlim([-sidRange sidRange])
    xline(0)

    sgtitle([strrep(filebase,'_',' ') ' (n = ' num2str(nFliesThresh) ')'])
    % save plot
    cd(folder.summary)
    plotname = 'p1background_sr_v_vel';
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vector plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

    disp('Complete.')
end

%% plot spike rate versus directional velocity for w/ w/o P1 activation w/lag
if nt_t>1
    disp('Analyzing behavior tuning for P1 on vs off w/lag...')

    % replace missing bins with nans
    srfwd_on_lag(srfwd_on_lag==0) = nan;
    srang_on_lag(srang_on_lag==0) = nan;
    srsid_on_lag(srsid_on_lag==0) = nan;
    srfwd_off_lag(srfwd_off_lag==0) = nan;
    srang_off_lag(srang_off_lag==0) = nan;
    srsid_off_lag(srsid_off_lag==0) = nan;

    % calculate means
    mean_srfwd_on = mean(srfwd_on_lag,2,'omitnan');
    mean_srang_on = mean(srang_on_lag,2,'omitnan');
    mean_srsid_on = mean(srsid_on_lag,2,'omitnan');
    mean_srfwd_off = mean(srfwd_off_lag,2,'omitnan');
    mean_srang_off = mean(srang_off_lag,2,'omitnan');
    mean_srsid_off = mean(srsid_off_lag,2,'omitnan');
    % calculate SEMs
    sem_srfwd_on = std(srfwd_on_lag,0,2,'omitnan')/sqrt(nFliesThresh);
    sem_srang_on = std(srang_on_lag,0,2,'omitnan')/sqrt(nFliesThresh);
    sem_srsid_on = std(srsid_on_lag,0,2,'omitnan')/sqrt(nFliesThresh);
    sem_srfwd_off = std(srfwd_off_lag,0,2,'omitnan')/sqrt(nFliesThresh);
    sem_srang_off = std(srang_off_lag,0,2,'omitnan')/sqrt(nFliesThresh);
    sem_srsid_off = std(srsid_off_lag,0,2,'omitnan')/sqrt(nFliesThresh);

    mean_srang_on(isnan(mean_srang_on)) = 0;
    sem_srang_on(isnan(sem_srang_on)) = 0;


    % initialize
    figure; set(gcf,'Position',[100 100 1000 400])

    % plot forward velocity
    subplot(1,3,1)
    % plot SEM band
    r(1) = patch([fwdBins; flipud(fwdBins)],[(mean_srfwd_off-sem_srfwd_off); flipud((mean_srfwd_off+sem_srfwd_off))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.bckColor{1});
    hold on
    r(2) = patch([fwdBins; flipud(fwdBins)],[(mean_srfwd_on-sem_srfwd_on); flipud((mean_srfwd_on+sem_srfwd_on))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.bckColor{2});
    % plot average
    plot(fwdBins,mean_srfwd_off,'Color',settings.bckColor{1},'MarkerFaceColor','w','LineWidth',settings.lwAvg)
    plot(fwdBins,mean_srfwd_on,'Color',settings.bckColor{2},'MarkerFaceColor','w','LineWidth',settings.lwAvg)
    ylabel(settings.spkLabel)
    xlabel(settings.velLabel{1})
    ylim([0 srRange])
    xlim([0 fwdRange])

    % plot angular
    subplot(1,3,2)
    % plot SEM band
    r(1) = patch([angBins; flipud(angBins)],[(mean_srang_off-sem_srang_off); flipud((mean_srang_off+sem_srang_off))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.bckColor{1});
    hold on
    r(2) = patch([angBins; flipud(angBins)],[(mean_srang_on-sem_srang_on); flipud((mean_srang_on+sem_srang_on))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.bckColor{2});
    % plot average
    plot(angBins,mean_srang_off,'Color',settings.bckColor{1},'MarkerFaceColor','w','LineWidth',settings.lwAvg)
    plot(angBins,mean_srang_on,'Color',settings.bckColor{2},'MarkerFaceColor','w','LineWidth',settings.lwAvg)
    xlabel(settings.velLabel{2})
    ylim([0 srRange])
    xlim([-200 200])
    xline(0)

    % plot sideway
    subplot(1,3,3)
    % plot SEM band
    r(1) = patch([sidBins; flipud(sidBins)],[(mean_srsid_off-sem_srsid_off); flipud((mean_srsid_off+sem_srsid_off))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.bckColor{1});
    hold on
    r(2) = patch([sidBins; flipud(sidBins)],[(mean_srsid_on-sem_srsid_on); flipud((mean_srsid_on+sem_srsid_on))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.bckColor{2});
    % plot average
    plot(sidBins,mean_srsid_off,'Color',settings.bckColor{1},'MarkerFaceColor','w','LineWidth',settings.lwAvg)
    plot(sidBins,mean_srsid_on,'Color',settings.bckColor{2},'MarkerFaceColor','w','LineWidth',settings.lwAvg)
    ylabel(settings.spkLabel)
    xlabel(settings.velLabel{3})
    ylim([0 srRange])
    xlim([-sidRange sidRange])
    xline(0)

    sgtitle([strrep(filebase,'_',' ') ' w/lag (n = ' num2str(nFliesThresh) ')'])
    % save plot
    cd(folder.summary)
    plotname = 'p1background_sr_v_vel_lag';
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vector plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

    disp('Complete.')
end

%% Plot slopes for directional velocity data
if nt_t>1
    % Define colors from settings
    scatter_color = settings.trialColor;      % Color for individual points
    p1_off_color = settings.bckColor{1};      % Gray for P1 off median
    p1_on_color = settings.bckColor{2};       % Red for P1 on median
    jitter_amount = 0.1;                      % Jitter for scatter points

    % Create figure and tiled layout
    figure; set(gcf, 'Position', [100 100 400 500])
    tiledlayout(1, 3, 'TileSpacing', 'compact')

    % Tile 1: Forward Velocity Slopes
    nexttile;
    hold on;
    scatter(ones(size(fwd_on_fits)) + jitter_amount * (rand(size(fwd_on_fits)) - 0.5), fwd_on_fits, ...
        '.', 'MarkerEdgeColor', scatter_color, 'MarkerFaceColor', scatter_color);
    scatter(2 * ones(size(fwd_off_fits)) + jitter_amount * (rand(size(fwd_off_fits)) - 0.5), fwd_off_fits, ...
        '.', 'MarkerEdgeColor', scatter_color, 'MarkerFaceColor', scatter_color);
    median_fwd_on = median(fwd_on_fits, 'omitnan');
    median_fwd_off = median(fwd_off_fits, 'omitnan');
    plot(1, median_fwd_on, '_', 'MarkerSize', 15, 'Color', p1_on_color, 'LineWidth', 2); % P1 on median
    plot(2, median_fwd_off, '_', 'MarkerSize', 15, 'Color', p1_off_color, 'LineWidth', 2); % P1 off median
    xticks([1 2]);
    xticklabels({'P1 On', 'P1 Off'});
    xlim([0 3]);
    yline(0);
    ylabel('Forward Velocity Slope');
    hold off;

    % Tile 2: Angular Velocity Slopes
    nexttile;
    hold on;
    scatter(ones(size(ang_on_fits)) + jitter_amount * (rand(size(ang_on_fits)) - 0.5), ang_on_fits, ...
        '.', 'MarkerEdgeColor', scatter_color, 'MarkerFaceColor', scatter_color);
    scatter(2 * ones(size(ang_off_fits)) + jitter_amount * (rand(size(ang_off_fits)) - 0.5), ang_off_fits, ...
        '.', 'MarkerEdgeColor', scatter_color, 'MarkerFaceColor', scatter_color);
    median_ang_on = median(ang_on_fits, 'omitnan');
    median_ang_off = median(ang_off_fits, 'omitnan');
    plot(1, median_ang_on, '_', 'MarkerSize', 15, 'Color', p1_on_color, 'LineWidth', 2); % P1 on median
    plot(2, median_ang_off, '_', 'MarkerSize', 15, 'Color', p1_off_color, 'LineWidth', 2); % P1 off median
    xticks([1 2]);
    xticklabels({'P1 On', 'P1 Off'});
    xlim([0 3]);
    yline(0);
    ylabel('Angular Velocity Slope');
    hold off;

    % Tile 3: Sideways Velocity Slopes
    nexttile;
    hold on;
    scatter(ones(size(sid_on_fits)) + jitter_amount * (rand(size(sid_on_fits)) - 0.5), sid_on_fits, ...
        '.', 'MarkerEdgeColor', scatter_color, 'MarkerFaceColor', scatter_color);
    scatter(2 * ones(size(sid_off_fits)) + jitter_amount * (rand(size(sid_off_fits)) - 0.5), sid_off_fits, ...
        '.', 'MarkerEdgeColor', scatter_color, 'MarkerFaceColor', scatter_color);
    median_sid_on = median(sid_on_fits, 'omitnan');
    median_sid_off = median(sid_off_fits, 'omitnan');
    plot(1, median_sid_on, '_', 'MarkerSize', 15, 'Color', p1_on_color, 'LineWidth', 2); % P1 on median
    plot(2, median_sid_off, '_', 'MarkerSize', 15, 'Color', p1_off_color, 'LineWidth', 2); % P1 off median
    xticks([1 2]);
    yline(0);
    xticklabels({'P1 On', 'P1 Off'});
    xlim([0 3]);
    ylabel('Sideways Velocity Slope');
    hold off;

    % Global title
    sgtitle('Velocity Slope Distributions');
    % save plot
    cd(folder.summary)
    plotname = 'sr_v_vel_fitted';
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

    % Define the structure to hold the slopes for P1 on and P1 off
    fits_on.fwd = fwd_on_fits;
    fits_on.ang = ang_on_fits;
    fits_on.sid = sid_on_fits;
    fits_off.fwd = fwd_off_fits;
    fits_off.ang = ang_off_fits;
    fits_off.sid = sid_off_fits;

    % Save the slopes as a .mat file in the specified folder
    cd(folder.compare);
    savename = [filebase '_vel_slopes.mat'];
    save(savename, 'fits_on', 'fits_off','flylist');
end

%% Plot difference in firing rate and voltage with and without P1

% Create figure and tiled layout
figure; set(gcf, 'Position', [100 100 400 500])
tiledlayout(1, 2, 'TileSpacing', 'compact')
sr_max = round(ceil(max(diff_sr,[],'all'))+5,-1);

% Define x-axis positions with reduced jitter for all data (1) and quiescent only (2)
x_all = 1 + 0.02 * randn(size(diff_sr, 1), 1); % Reduced jitter for all data
x_quiet = 2 + 0.02 * randn(size(diff_sr, 1), 1); % Reduced jitter for quiet-only data

% Calculate medians for all data and quiet-only data
median_diff_sr_all = median(diff_sr(:,1), 'omitnan');
median_diff_sr_quiet = median(diff_sr(:,2), 'omitnan');
median_diff_vm_all = median(diff_vm(:,1), 'omitnan');
median_diff_vm_quiet = median(diff_vm(:,2), 'omitnan');

% Plot for Firing Rate Difference
nexttile;
hold on;
scatter(x_all, diff_sr(:,1), [], settings.trialColor, '.');
scatter(x_quiet, diff_sr(:,2), [], settings.trialColor, '.');
plot(1, median_diff_sr_all, 'k--', 'Marker', '_', 'MarkerSize', 15); % Dash marker for median (all)
plot(2, median_diff_sr_quiet, 'k--', 'Marker', '_', 'MarkerSize', 15); % Dash marker for median (quiet)
title('Firing Rate Difference');
ylabel('Difference in Firing Rate');
xticks([1 2]);
xlim([0 3])
ylim([-sr_max sr_max])
yline(0)
xticklabels({'all', 'rest'});


% Plot for Voltage Difference
nexttile;
hold on;
scatter(x_all, diff_vm(:,1), [], settings.trialColor, '.');
scatter(x_quiet, diff_vm(:,2), [], settings.trialColor, '.');
plot(1, median_diff_vm_all, 'k--', 'Marker', '_', 'MarkerSize', 15); % Dash marker for median (all)
plot(2, median_diff_vm_quiet, 'k--', 'Marker', '_', 'MarkerSize', 15); % Dash marker for median (quiet)
title('Voltage Difference');
ylabel('Difference in Voltage');
xticks([1 2]);
xlim([0 3])
yline(0)
xticklabels({'all', 'rest'});

sgtitle([strrep(filebase,'_',' ') ' (n = ' num2str(nFlies) ')'])
% save plot
cd(folder.summary)
plotname = 'p1background_change';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vector plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

% Save the differences as a .mat file in the specified folder
cd(folder.compare);
savename = [filebase '_spike_voltage_diffs.mat'];
save(savename, 'diff_sr', 'diff_vm');

disp('Complete.')

%% Plot percent run time for P1 on v off
if nt_t>1
    % Plot percent time running with and without P1
    figure; set(gcf, 'Position', [600 100 300 400])
    hold on;

    % X positions for off and on
    xvals = [1 2];

    % Plot individual data and lines connecting points per fly
    for i = 1:length(runpercent_off)
        plot(xvals, [runpercent_off(i), runpercent_on(i)], '-', 'Color', [0.7 0.7 0.7]) % light grey lines
    end

    % Scatter individual points
    scatter(repmat(1, size(runpercent_off)), runpercent_off, 20, 'k', '.');
    scatter(repmat(2, size(runpercent_on)), runpercent_on, 20, 'k', '.');

    % Plot medians
    med_off = median(runpercent_off, 'omitnan');
    med_on = median(runpercent_on, 'omitnan');
    plot(1, med_off, 'k_', 'MarkerSize', 15, 'LineWidth', 1.5);
    plot(2, med_on, 'r_', 'MarkerSize', 15, 'LineWidth', 1.5); % red median for P1 on

    % Axis labels and formatting
    xlim([0.5 2.5]);
    ylim([0 100])
    xticks([1 2]);
    xticklabels({'P1 off', 'P1 on'});
    ylabel('% Time Running');
    title('Running with vs. without P1');

    % Save figure
    cd(folder.summary);
    saveas(gcf, 'p1_runpercent.png');
    copyfile('p1_runpercent.png', folder.dropbox, 'f');
    cd(folder.vector);
    set(gcf,'renderer','Painters');
    saveas(gcf, 'p1_runpercent.svg');
    copyfile('p1_runpercent.svg', folder.dropbox, 'f');

    % Save run percent data
    cd(folder.compare);
    save([filebase '_runpercent_data.mat'], 'runpercent_off', 'runpercent_on');
end

%% analyze directional velocity v spike rate w/lag
if nt_t>1
    disp('Analyzing directional velocity v spike rate w/lag...')

    thisSpikeBin = binSR_withLag.spikeRateBin;
    sr_lim = [0 80];

    figure; set(gcf,'Position',[100 100 800 800])
    for v = 1:3
        switch v
            case 1 % fwd
                thisData = fwdsr_on_lag;
                thisYlim = [0 10];
            case 2 % ang
                thisData = angsr_on_lag;
                thisYlim = [-200 200];
            case 3 % sid
                thisData = sidsr_on_lag;
                thisYlim = [-2 2];
        end

        % count number of animals with data per bin
        validN = sum(~isnan(thisData), 2);

        % only use bins with at least 3 animals
        validIdx = validN >= 5;
        xVals = thisSpikeBin(validIdx);
        yData = thisData(validIdx, :);
        validN_used = validN(validIdx); % same size as xVals

        % compute mean and SEM across valid data
        thisMean = mean(yData, 2, 'omitnan');
        thisSEM  = std(yData, [], 2, 'omitnan') ./ sqrt(validN_used);

        % ensure xVals is a column vector
        xVals = xVals(:);

        % plot mean + individual fly traces
        subplot(2,3,v); hold on
        plot(xVals, yData, 'Color', settings.trialColor)
        plot(xVals, thisMean, 'Color', settings.velColor{v}, 'LineWidth', settings.lwAvg)
        axis tight; ylim(thisYlim); xlim(sr_lim); yline(0);
        ylabel(settings.velLabel{v}); xlabel('Spike rate (spikes/s)')

        % plot mean ± SEM
        subplot(2,3,v+3); hold on
        plot(xVals, thisMean, 'Color', settings.velColor{v}, 'LineWidth', settings.lwAvg)
        er = patch([xVals; flipud(xVals)], ...
            [(thisMean - thisSEM); flipud(thisMean + thisSEM)], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
        if ~isempty(er)
            er.FaceColor = settings.velColor{v};
        end
        axis tight; ylim(thisYlim); xlim(sr_lim); yline(0)
        ylabel(settings.velLabel{v}); xlabel('Spike rate (spikes/s)')
    end

    sgtitle([strrep(filebase,'_','/') ' Velocity by Spike Rate (n = ' num2str(nFliesThresh) ')'])

    % save plot
    cd(folder.summary)
    plotname = 'pulse_vel_v_sr_lag';
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');

    % save vectorized plot
    cd(folder.vector)
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

end

%% Plot relationship between change in membrane voltage and firing rate
figure('Color','w');
set(gcf, 'Position', [100, 100, 850, 800]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% ---------- (1) TOP-LEFT: Measured FR vs Vm (OFF black, ON red) ----------
nexttile; hold on;

% OFF (black)
for nt = 1:size(binFRVm_off,2)
    plot(vmBins, binFRVm_off(:,nt), 'k', 'LineWidth', 0.5);
end

% ON (red)
for nt = 1:size(binFRVm_on,2)
    plot(vmBins, binFRVm_on(:,nt), 'r', 'LineWidth', 0.5);
end

% Limits
allY = [binFRVm_off(:); binFRVm_on(:)];
ymin = min(allY(~isnan(allY))); if isempty(ymin), ymin = 0; end
ymax = max(allY(~isnan(allY))); if isempty(ymax), ymax = 1; end
xlim([min(vmBins) max(vmBins)]);
ylim([ymin ymax]);
grid on;
xlabel('Membrane voltage (mV)');
ylabel('Firing rate (Hz)');
title('Measured FR vs Vm');

% ---------- (2) TOP-RIGHT: ΔFR vs ΔVm (zero-shifted per fly) ----------
nexttile; hold on;

nFlies = max(size(binFRVm_off,2), size(binFRVm_on,2));
xmins = []; xmaxs = []; ymins = []; ymaxs = [];

for nt = 1:nFlies
    % Handle unequal columns safely
    y_off = nan(size(vmBins)); if nt <= size(binFRVm_off,2), y_off = binFRVm_off(:,nt); end
    y_on  = nan(size(vmBins)); if nt <= size(binFRVm_on,2),  y_on  = binFRVm_on(:,nt);  end

    validBins = ~isnan(y_off) | ~isnan(y_on);
    if ~any(validBins), continue; end

    minVm_nt = min(vmBins(validBins));
    minFR_nt = min([y_off(validBins); y_on(validBins)], [], 'omitnan');

    x_shift = vmBins - minVm_nt;
    y_off_shift = y_off - minFR_nt;
    y_on_shift  = y_on  - minFR_nt;

    plot(x_shift, y_off_shift, 'k', 'LineWidth', 0.5);
    plot(x_shift, y_on_shift,  'r', 'LineWidth', 0.5);

    xmins(end+1,1) = min(x_shift(validBins));
    xmaxs(end+1,1) = max(x_shift(validBins));
    ymins(end+1,1) = min([y_off_shift(validBins); y_on_shift(validBins)], [], 'omitnan');
    ymaxs(end+1,1) = max([y_off_shift(validBins); y_on_shift(validBins)], [], 'omitnan');
end

if ~isempty(xmins), xlim([min(xmins) max(xmaxs)]); else, xlim([0 1]); end
if ~isempty(ymins), ylim([min(ymins) max(ymaxs)]); else, ylim([0 1]); end

grid on;
xlabel('\Delta Membrane voltage (mV)');
ylabel('\Delta Firing rate (Hz)');
title('\DeltaFR vs \DeltaVm');

% ---------- (3) BOTTOM-LEFT: Vm distribution (probability) ----------
nexttile; hold on;

% OFF (black probabilities)
for nt = 1:size(binProbVm_off,2)
    plot(vmBins, binProbVm_off(:,nt), 'k', 'LineWidth', 0.5);
end

% ON (red probabilities)
for nt = 1:size(binProbVm_on,2)
    plot(vmBins, binProbVm_on(:,nt), 'r', 'LineWidth', 0.5);
end

% Limits for probabilities
allP = [binProbVm_off(:); binProbVm_on(:)];
pmax = max(allP(~isnan(allP))); if isempty(pmax), pmax = 1; end
xlim([min(vmBins) max(vmBins)]);
ylim([0 pmax]);
grid on;
xlabel('Membrane voltage (mV)');
ylabel('Probability');
title('Vm distribution (prob.)');

% ---------- (4) BOTTOM-RIGHT: \DeltaVm distribution (probability) ----------
nexttile; hold on;

xmins = []; xmaxs = [];  % reuse for shifted x-lims

nFliesP = max(size(binProbVm_off,2), size(binProbVm_on,2));
for nt = 1:nFliesP
    p_off = nan(size(vmBins)); if nt <= size(binProbVm_off,2), p_off = binProbVm_off(:,nt); end
    p_on  = nan(size(vmBins)); if nt <= size(binProbVm_on,2),  p_on  = binProbVm_on(:,nt);  end

    % Use same validBins definition as FR (bins with any data)
    validBins = (~isnan(p_off) & p_off>0) | (~isnan(p_on) & p_on>0);
    % If probability arrays may have zeros but FR had data, optionally fall back:
    if ~any(validBins)
        % fallback to any non-NaN (keeps alignment with vmBins if sparsity exists)
        validBins = ~isnan(p_off) | ~isnan(p_on);
    end
    if ~any(validBins), continue; end

    % min Vm per fly across bins with any probability mass (or fallback)
    minVm_nt = min(vmBins(validBins));
    x_shift = vmBins - minVm_nt;

    % Plot shifted probabilities (y stays as prob; x is shifted)
    plot(x_shift, p_off, 'k', 'LineWidth', 0.5);
    plot(x_shift, p_on,  'r', 'LineWidth', 0.5);

    xmins(end+1,1) = min(x_shift(validBins));
    xmaxs(end+1,1) = max(x_shift(validBins));
end

if ~isempty(xmins), xlim([min(xmins) max(xmaxs)]); else, xlim([0 1]); end
ylim([0 pmax]); % re-use global prob max for consistent scaling
grid on;
xlabel('\Delta Membrane voltage (mV)');
ylabel('Probability');
title('\DeltaVm distribution (prob.)');

% Overall title
title(tl, 'FR–Vm and Vm distributions (OFF=black, ON=red)', 'FontWeight','bold');

% save plot
cd(folder.summary)
plotname = 'pulse_nonlinear';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');

% save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

%% Plot firing rate distribution
% Bin centers for plotting
binCenters = fr_bins(1:end-1) + diff(fr_bins)/2;   % 0:20:200 -> centers at 10:20:190

% Mean and SEM across animals
meanCounts_on = mean(fr_counts_on, 1, 'omitnan');
semCounts_on  = std(fr_counts_on, 0, 1, 'omitnan') ./ sqrt(nFlies);
meanCounts_off = mean(fr_counts_off, 1, 'omitnan');
semCounts_off  = std(fr_counts_off, 0, 1, 'omitnan') ./ sqrt(nFlies);

% Figure and layout
figure; set(gcf, 'Position', [100 100 500 800]); hold on;

% SEM patch
x_patch = [binCenters, fliplr(binCenters)];
y_patch = [meanCounts_on - semCounts_on, fliplr(meanCounts_on + semCounts_on)];
sp = patch(x_patch(:), y_patch(:), 'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
sp.FaceColor = 'r';  % match line color
% Mean line
plot(binCenters, meanCounts_on, 'LineWidth', 2, 'Color', 'r');

% SEM patch
x_patch = [binCenters, fliplr(binCenters)];
y_patch = [meanCounts_off - semCounts_off, fliplr(meanCounts_off + semCounts_off)];
sp2 = patch(x_patch(:), y_patch(:), 'k', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
sp2.FaceColor = 'k';  % match line color
% Mean line
plot(binCenters, meanCounts_off, 'LineWidth', 2, 'Color', 'k');

% Axes/labels
xlim([fr_bins(1) fr_bins(end)]);
ylim([0, 0.5]);
xlabel('Firing rate (Hz)');
ylabel('Probability');
title('Firing rate distribution (mean \pm SEM across animals)');
box off; set(gca, 'Layer', 'top');

% save plot
cd(folder.summary)
plotname = 'fr_distribution';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

%% end
disp('ALL ANALYSES COMPLETE.')
end

