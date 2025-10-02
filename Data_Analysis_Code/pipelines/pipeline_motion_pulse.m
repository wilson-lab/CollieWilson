% pipeline_motion_pulse
%
% Pipeline Function
% Pulls all processed files from ALL flies in a given experiment, performs
% necessary analyses and plots accordingly. Can be used for both behavior
% only experiments and ephys experiments.
%
% INPUTS
% exptFolder  - overarching experiment folder
% pulseSpeeds - cell array of pulse speeds presented (e.g., 0, 25, 75 dps)
%
% This function processes motion pulse data, performing analyses on
% spikerate, directional velocity, voltage responses, cross-correlation,
% and directional selectivity. Results are saved as summary figures and plots.
%
% 06/14/2023 - MC adapted from battery pipeline
% 07/05/2024 - MC streamlined
% 11/07/2024 - MC added model fit to firing rate
% 12/18/2024 - MC added visual response latency calculation
%
function pipeline_motion_pulse(exptFolder,pulseSpeeds)
disp('ANALYZING POOLED TARGET PULSE EXPT...')
%% Initialize Folder Structure and Settings
% Set up filename based on the experiment folder by replacing spaces with underscores
filebase = strrep(exptFolder, ' ', '_');

% Generate directories required for file storage and organization
folder = generateFolders(exptFolder);

% Load processing settings that configure parameters for subsequent analyses
settings = processSettings();

%% Initialize Motion Pulse Experiment Variables
% Set variables based on the motion pulse stimulus parameters

% Determine the number of pulse speeds presented in the experiment
nPulses = length(pulseSpeeds);

% Check if the stimulus includes motion by confirming that the first pulse speed is non-zero
motionCheck = ~contains(pulseSpeeds{1}, '0dps');

% Configure sweep count, plot dimensions, and model fit parameters based on the experiment paradigm
if motionCheck && contains(filebase, 'AOTU019')
    % Settings specific to AOTU019 motion experiments
    nSweep = 9;         % Number of sweeps for AOTU019 experiments
    tileX = 3;          % Number of tiles in the x-direction for plots
    plotWidth = 800;    % Plot width in pixels

    % Define model fit parameters for position peaks
    peak_position = 35;         % Peak position (degrees)
    spread = 20;                % Expected spread of RF(degrees)

elseif motionCheck && contains(filebase, 'AOTU025')
    % Settings specific to AOTU025 motion experiments
    nSweep = 13;
    tileX = 3;
    plotWidth = 800;

    % Define model fit parameters for position peaks
    peak_position = 65;
    spread = 20;

elseif motionCheck && contains(filebase, 'DNa02')
    % Settings specific to DNa02 motion experiments
    nSweep = 9;
    tileX = 3;
    plotWidth = 800;

    % Define model fit parameters for position peaks
    peak_position = 65;
    spread = 20;

elseif (~motionCheck) && contains(filebase, 'AOTU019')
    % Default settings for stationary or other conditions
    nSweep = 17;
    tileX = 2;
    plotWidth = 600;

    % Define model fit parameters for position peaks
    peak_position = 35;
    spread = 20;

elseif (~motionCheck) && contains(filebase, 'AOTU025')
    % Default settings for stationary or other conditions
    nSweep = 25;
    tileX = 2;
    plotWidth = 600;

    % Define model fit parameters for position peaks
    peak_position = 65;
    spread = 20;
end

if contains(filebase, 'aIPg')
    % Settings specific to aIPG motion experiments
    nSweep = 9;
    tileX = 3;
    plotWidth = 800;
end
%% Initialize Plotting Variables

% Spike Rate Limits
sr_limitT = [0 80];     % Spike rate limits for trials
sr_limitA = [0 80];     % Spike rate limits for averaged data
sr_limitD = [-35 60];   % Spike rate limits for difference plots

% Voltage Limits
vm_diff = [-8 8];       % Voltage difference range

% Directional Velocity Limits
ang_limit = [-100 100]; % Angular velocity limits
ang_limit2 = [-150 150]; % Angular velocity limits

% General Settings
runSelect = [-1, 0, settings.runThreshE];  % Selection threshold for running state
n_limit = [-1.5 1.5];                      % Normalized average limits
n_limitD = [0 1.5];                        % Normalized difference limits
ds_limit = [-1 1];                         % Direction selectivity limits
ps_limit = [-150 150];                     % Position selectivity limits

% Sweep Colors for Plotting
nSweep_ipsi = ceil(nSweep / 2);        % Number of sweeps in the ipsilateral direction
nSweep_cntr = floor(nSweep / 2);       % Number of sweeps in the contralateral direction

% Colormap setup: bone for contralateral, hsv for ipsilateral sweeps
cm1 = flip(colormap(bone(nSweep_cntr + 3))); % Adjusted colormap for contralateral sweeps
cm1(end, :) = []; cm1(1, :) = [];            % Remove outer colors for better contrast
cm2 = flip(colormap(hsv(nSweep_ipsi)));      % Colormap for ipsilateral sweeps

% Set color scheme for rightward and leftward sweeps
color_rightward = [cm1(2:end, :); cm2];      % Colors for rightward sweeps
color_leftward = [cm1; cm2(1:end-1, :)];     % Colors for leftward sweeps

% Close all previous figures to ensure a clean slate for plotting
close all

%% Load and Pool All Trials from Each Experiment Folder

disp('Loading datasets...')
% Navigate to the folder containing intermediate data files
cd(folder.int)

% Identify all .mat files in the directory for analysis
allFiles = dir('*int.mat');
nFlies = length(allFiles);  % Total number of flies in the experiment
flylist = [];

% Initialize data storage arrays for each fly's data
nFliesThresh = nFlies;  % Thresholded count of flies meeting behavior criteria
nt_t = 0;               % Counter for valid trials meeting run criteria
allPanelPs = cell(1, nFlies);    % Storage for panel positions
allForward = cell(1, nFlies);    % Storage for forward velocities
allAngular = cell(1, nFlies);    % Storage for angular velocities
allSideway = cell(1, nFlies);    % Storage for sideways velocities
allSpikeRt = cell(1, nFlies);    % Storage for spike rates
storeNames = {};

% Loop through each trial, loading and processing data for each fly
for nt = 1:nFlies
    % Load trial data for the current fly
    cd(folder.int)
    thisTrial = allFiles(nt).name;
    thisFly = thisTrial(6:16);              % Extract the fly name from the filename
    flyShortNames{nt} = strrep(thisFly, '_', ' ');  % Convert fly name format for readability
    disp(['Loading ' num2str(nt) '/' num2str(nFlies) ': ' flyShortNames{nt}])
    load(thisTrial)

    % Store data arrays for this fly's trial
    allPanelPs{nt} = int_panelps;
    allForward{nt} = int_forward;
    allAngular{nt} = int_angular;
    allSideway{nt} = int_sideway;
    allSpikeRt{nt} = int_spikert;

    % Apply median filtering to voltage data to remove transient spikes
    mf_voltage = spikeFilter(int_voltage, int_time);
    allVoltage{nt} = mf_voltage;

    % Calculate the fly's total walking time and check if it meets run threshold
    flyRunTime(nt, 1) = (sum(int_forward > settings.runThreshE, 'all') / length(int_time)) * 60;

    % Only include flies that exhibited sufficient running behavior
    if flyRunTime(nt, 1) > settings.minRunTime
        nt_t = nt_t + 1;  % Increment counter for valid trials
        thisFlyLong = extractBefore(thisTrial, '_int');  % Extract fly ID
        storeNames{nt_t} = thisFlyLong;

        % Analyze behavior tuning for the current trial
        thisVel = spikert_binvelocity(int_forward, int_angular, int_sideway, int_spikert, int_time, 0);
        thisVelL = spikert_binvelocity(int_forward, int_angular, int_sideway, int_spikert, int_time, 1);
        thisAcc = spikert_binacceleration(int_forward, int_angular, int_sideway, int_spikert, int_time, 1);
        % Bin spikerate by velocity with lag estimate
        [binSR_withLag] = velocity_binbyspikerate(int_forward, int_angular, int_sideway, int_spikert, int_time, 1);
        % Fit relationship with velocity
        [slope_fwd, slope_ang, slope_sid, r2_fwd, ~, ~] = fit_velocity_data(int_forward, int_angular, int_sideway, int_spikert, int_time, 1);

        % Store the behavior tuning results for this trial
        vel_srvfwd(:, nt_t) = thisVel.fwdMean';
        vel_srvang(:, nt_t) = thisVel.angMean';
        vel_srvsid(:, nt_t) = thisVel.sidMean';
        velL_srvfwd(:, nt_t) = thisVelL.fwdMean';
        velL_srvang(:, nt_t) = thisVelL.angMean';
        velL_srvang_fastfwd(:, nt_t) = thisVelL.angMean_fastFwd';
        velL_srvang_slowfwd(:, nt_t) = thisVelL.angMean_slowFwd';
        velL_srvsid(:, nt_t) = thisVelL.sidMean';

        sr_velfwd(:,nt_t) = binSR_withLag.fwdMean';
        sr_velang(:,nt_t) = binSR_withLag.angMean';
        sr_velsid(:,nt_t) = binSR_withLag.sidMean';
        acc_srvfwd(:, nt_t) = thisAcc.fwdMean';
        acc_srvang(:, nt_t) = thisAcc.angMean';
        acc_srvsid(:, nt_t) = thisAcc.sidMean';

        fwd_fits(nt_t) = slope_fwd;
        fwd_r2(nt_t) = r2_fwd;
        ang_fits(nt_t) = slope_ang;
        sid_fits(nt_t) = slope_sid;

        % Model firing rate as a function of motor and visual parameters
        [R2_motor,~] = spikert_seqfit(int_spikert, int_panelps, int_forward, int_angular, int_time, peak_position, spread);
        R2_fwd = compareSpikeRateModels(int_spikert, int_panelps, int_forward);
        % Store sequential R2 values
        R2_seq(nt_t,:) = R2_motor;

        R2_fwdmodel(nt_t,1) = R2_fwd.additive;
        R2_fwdmodel(nt_t,2) = R2_fwd.multiplicative;
        R2_fwdmodel(nt_t,3) = R2_fwd.full;

        % Perform or load cross-correlation analysis for firing rate
        cd(folder.xcorr)
        thisXCorrFile = [thisFly '_xc.mat'];
        if exist(thisXCorrFile, 'file')
            % Load previously computed cross-correlation results
            disp('Loading previous cross-correlation data.')
            load(thisXCorrFile)
        else
            % Perform new cross-correlation analysis
            [r_val, lag_t] = spikert_xcorr(int_spikert, int_forward, int_angular, int_sideway, int_time);
            % Save results for future access
            save(thisXCorrFile, 'r_val', 'lag_t', '-v7.3');
        end
        % Find and store peak lag and peak r-values for each motion type
        [peak_lag, peak_rval,r_val] = find_peak_lag_rval(r_val, lag_t, 0.1);

        % Store results in nt_t for each motion type
        r_val_fwd(:, nt_t) = r_val.fwd;
        pk_lag_fwd(nt_t) = peak_lag.fwd;
        pk_rval_fwd(nt_t) = peak_rval.fwd;
        [peak_lag, peak_rval,r_val] = find_peak_lag_rval(r_val, lag_t, settings.minXCorrProm);
        r_val_ang(:, nt_t) = r_val.ang;
        pk_lag_ang(nt_t) = peak_lag.ang;
        pk_rval_ang(nt_t) = peak_rval.ang;
        r_val_sid(:, nt_t) = r_val.sid;
        pk_lag_sid(nt_t) = peak_lag.sid;
        pk_rval_sid(nt_t) = peak_rval.sid;

        flylist{nt_t} = thisFly;
    else
        % Exclude flies that did not meet the minimum running behavior requirement
        disp([thisFly ' omitted from behavior analyses.'])
        nFliesThresh = nFliesThresh - 1;
    end
end

%% Plot model fits
if nt_t>1
    % Transpose R2_seq so each row represents a model and each column represents an animal
    R2_seq_transposed = R2_seq';
    num_models = size(R2_seq_transposed, 1);
    num_animals = size(R2_seq_transposed, 2);

    % Calculate R² improvement (absolute increase) for each addition
    R2_improvement = zeros(3, num_animals); % Rows: [Forward, Angular, Gaussian]

    for animal = 1:num_animals
        % Initial R² for forward velocity component
        R2_improvement(1, animal) = R2_seq_transposed(1, animal); % R² for forward only

        % Improvement from adding Angular Velocity
        R2_improvement(2, animal) = R2_seq_transposed(2, animal) - R2_seq_transposed(1, animal);

        % Improvement from adding Gaussian Model (Object Position)
        R2_improvement(3, animal) = R2_seq_transposed(3, animal) - R2_seq_transposed(2, animal);
    end

    % Calculate median values across animals
    R2_median = median(R2_seq_transposed, 2); % Median R² values for each model
    R2_improvement_median = median(R2_improvement, 2); % Median R² improvement for each addition

    % Create a tiled layout with two plots
    figure; set(gcf, 'Position', [100 100 700 800]);
    tiledlayout(1, 2);

    % First plot: Absolute R² values across model runs
    nexttile;
    hold on;
    for animal = 1:num_animals
        plot(1:num_models, R2_seq_transposed(:, animal), '-', 'Color', settings.trialColor);
    end
    % Plot median R² values across animals
    plot(1:num_models, R2_median, '-', 'Color', 'k','LineWidth',1);
    xlabel('Model Runs');
    ylabel('R² Value');
    title('Model Fits');
    xticks(1:num_models);
    xticklabels({'Fwd', 'Fwd+Ang', 'Fwd+Ang+Obj'});
    axis padded; ylim([0 0.6]);
    grid on;
    hold off;

    % Second plot: Absolute R² improvement for each addition
    nexttile;
    hold on;
    for animal = 1:num_animals
        plot(1:3, R2_improvement(:, animal), '-', 'Color', settings.trialColor);
    end
    % Plot median R² improvements across animals
    plot(1:3, R2_improvement_median, '-', 'Color', 'k','LineWidth',1);
    xlabel('Predictor Added');
    ylabel('R² Improvement');
    title('Model Improvements');
    xticks([1 2 3]);
    xticklabels({'Fwd', '+Angular', '+Object'});
    axis padded; ylim([0 0.6]);
    grid on;
    hold off;

    % Save plot
    cd(folder.summary)
    plotname = 'incremental_model';
    saveas(gcf, [plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox, 'f');
    % Save vectorized plot
    cd(folder.vector)
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox, 'f');
end

%% analyze each motion pulse vs turn response
disp('Analying motion pulse vs turn response...')
cd(folder.summary)
if nt_t>1
    % run (1) all (3) running only
    for r = [1 3]
        % pull run info
        thisName = settings.behaviorGroup{r};
        thisRun = runSelect(r);

        % for each pulse speed
        for p = 1:nPulses
            % initialize
            thisSpeed = pulseSpeeds{p}(1:2);
            pulse_turnR = [];
            pulse_turnL = [];
            pulse_turnRL = [];
            c = 1; %fly counter

            % for each fly
            for nt = 1:nFlies
                if flyRunTime(nt)>settings.minRunTime
                    % select data
                    thisPanelPs = allPanelPs{nt};
                    thisForward = allForward{nt};
                    thisAngular = allAngular{nt};
                    % determine relationship between pulse and behavior
                    [~, thisMean] = pulse_v_output(thisPanelPs,thisForward,thisAngular,int_time,p,pulseSpeeds,nSweep,thisRun);
                    % if first run, store motion pulse positions
                    if c==1
                        % sweep positions
                        pulse_posR = thisMean.panelpsR;
                        pulse_posL = thisMean.panelpsL;
                        pulseDur = size(pulse_posR,1);
                        % sweep centers for each position
                        sweepPosR = round(pulse_posR(round(pulseDur/2),1,:));
                        sweepPosL = round(pulse_posL(round(pulseDur/2),1,:));
                        % sweep indices
                        sweepIdx = find(~isnan(pulse_posR(:,1)));
                        sweepIdx2End = sweepIdx(1):length(pulse_posR(:,1));
                    end
                    % store turn averages for each fly
                    pulse_turnR(:,c,:) = thisMean.varOutR;
                    pulse_turnL(:,c,:) = thisMean.varOutL;
                    pulse_turnRL(:,c,:) = thisMean.varOutRL;
                    % for each sweep, store peak turning for each fly
                    for s = 1:nSweep
                        if sweepPosR(s) > -settings.binoc % for rightward sweeps on right or binoc, take max right
                            peak_turnR(c,s) = max(thisMean.varOutR(sweepIdx2End,:,s));
                            peak_turnRL(c,s) = max(thisMean.varOutRL(sweepIdx2End,:,s));
                        else % else take max left
                            peak_turnR(c,s) = min(thisMean.varOutR(sweepIdx2End,:,s));
                            peak_turnRL(c,s) = min(thisMean.varOutRL(sweepIdx2End,:,s));
                        end
                        if sweepPosL(s) < settings.binoc % for leftward sweeps on left or binoc, take max left
                            peak_turnL(c,s) = min(thisMean.varOutL(sweepIdx2End,:,s));
                        else % else take max right
                            peak_turnL(c,s) = max(thisMean.varOutL(sweepIdx2End,:,s));
                        end
                    end
                    c = c+1; %update counter
                end
            end
            % calculate turn means
            mean_pulse_turnR = mean(pulse_turnR,2,'omitnan');
            mean_pulse_turnL = mean(pulse_turnL,2,'omitnan');
            mean_pulse_turnRL = mean(pulse_turnRL,2,'omitnan');
            % calculate turn SEMs
            sem_pulse_turnR = std(pulse_turnR,0,2,'omitnan')./sqrt(nFliesThresh);
            sem_pulse_turnL = std(pulse_turnL,0,2,'omitnan')./sqrt(nFliesThresh);
            sem_pulse_turnRL = std(pulse_turnRL,0,2,'omitnan')./sqrt(nFliesThresh);

            % calculate turn mean peaks
            for s = 1:nSweep
                if sweepPosR(s) > -settings.binoc
                    mean_peak_turnR(s) = max(mean_pulse_turnR(sweepIdx2End,:,s));
                    mean_peak_turnRL(s) = max(mean_pulse_turnRL(sweepIdx2End,:,s));
                else
                    mean_peak_turnR(s) = min(mean_pulse_turnR(sweepIdx2End,:,s));
                    mean_peak_turnRL(s) = min(mean_pulse_turnRL(sweepIdx2End,:,s));
                end
                if sweepPosL(s) < settings.binoc
                    mean_peak_turnL(s) = min(mean_pulse_turnL(sweepIdx2End,:,s));
                else
                    mean_peak_turnL(s) = max(mean_pulse_turnL(sweepIdx2End,:,s));
                end
            end
            % calculate turn mean peak SEMs
            sem_peak_turnR = std(peak_turnR,0,1,'omitnan')./sqrt(nFliesThresh);
            sem_peak_turnL = std(peak_turnL,0,1,'omitnan')./sqrt(nFliesThresh);
            sem_peak_turnRL = std(peak_turnRL,0,1,'omitnan')./sqrt(nFliesThresh);

            % plot averages
            % initialize
            figure; set(gcf,'Position',[100 100 plotWidth 800])
            tiledlayout(4,tileX,"TileSpacing","compact")
            % initialize remaining variables
            time_pulse = int_time(1:size(pulse_posR,1))*1000;
            time_limit = [0 max(time_pulse)];

            % plot target position
            % plot rightward target
            nexttile; hold on
            for s = 1:nSweep
                plot(time_pulse,pulse_posR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
            end
            title('Rightward Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(ps_limit); yline(0)
            if motionCheck
                % plot leftward target
                nexttile; hold on
                for s = nSweep:-1:1
                    plot(time_pulse,pulse_posL(:,s),'Color',color_leftward(s,:),'LineWidth',settings.lwAvg)
                end
                title('Leftward Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
                xlim(time_limit); ylim(ps_limit); yline(0)
            end
            % plot combined
            nexttile; hold on
            for s = nSweep:-1:1
                plot(time_pulse,pulse_posR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
            end
            title('R+L Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(ps_limit); yline(0)

            % plot behavior
            % plot rightward turning
            nexttile([3 1]); hold on
            for s = 1:nSweep
                sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_turnR(:,s)-sem_pulse_turnR(:,s); flipud(mean_pulse_turnR(:,s)+sem_pulse_turnR(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                sp.FaceColor = color_rightward(s,:);
                plot(time_pulse,mean_pulse_turnR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
            end
            ylabel('Angular Velocity (deg/sec)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(ang_limit); yline(0)
            if motionCheck
                % plot leftward turning
                nexttile([3 1]); hold on
                for s = 1:nSweep
                    sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_turnL(:,s)-sem_pulse_turnL(:,s); flipud(mean_pulse_turnL(:,s)+sem_pulse_turnL(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                    sp.FaceColor = color_leftward(s,:);
                    plot(time_pulse,mean_pulse_turnL(:,s),'Color',color_leftward(s,:),'LineWidth',settings.lwAvg)
                end
                ylabel('Angular Velocity (deg/sec)'); xlabel('Time (msec)')
                xlim(time_limit); ylim(ang_limit); yline(0)
            end
            % plot combined
            nexttile([3 1]); hold on
            for s = 1:nSweep
                sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_turnRL(:,s)-sem_pulse_turnRL(:,s); flipud(mean_pulse_turnRL(:,s)+sem_pulse_turnRL(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                sp.FaceColor = color_rightward(s,:);
                plot(time_pulse,mean_pulse_turnRL(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
            end
            ylabel('Angular Velocity (deg/sec)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(ang_limit); yline(0)

            sgtitle([strrep([filebase ' ' pulseSpeeds{p} ' ' thisName],'_','/') ' (n = ' num2str(nFliesThresh) ')'])
            % save plot
            cd(folder.summary)
            plotname = strjoin({'pulse_v_turn', thisSpeed, 'dps', thisName},'_');
            saveas(gcf,[plotname '.png']);
            copyfile([plotname '.png'], folder.dropbox,'f');
            % save vectorized plot
            cd(folder.vector)
            set(gcf,'renderer','Painters')
            saveas(gcf, [plotname '.svg'])
            copyfile([plotname '.svg'], folder.dropbox,'f');

            % estimate delay
            try
                thisDelay = pulse_estdelay(pulse_posR,mean_pulse_turnRL,time_pulse);
                estDelay(:,p,r) = thisDelay;
            catch
                estDelay(:,p,r) = nan(nSweep,1);
            end

            % plot normalized peak
            % initialize
            figure; set(gcf,'Position',[100 100 plotWidth 500])
            tiledlayout(1,tileX,"TileSpacing","compact")
            % normalize peak to 1
            normPeakR = max(abs(mean_peak_turnR));
            normPeakL = max(abs(mean_peak_turnL));
            normPeakRL = max(abs(mean_peak_turnRL));
            % plot peak turns for rightward sweeps
            nexttile;hold on
            sweepPosR = reshape(sweepPosR,[],1);
            norm_mean_peak_turnR = mean_peak_turnR ./ normPeakR;
            norm_sem_peak_turnR = sem_peak_turnR ./ normPeakR;
            errorbar(sweepPosR, norm_mean_peak_turnR, norm_sem_peak_turnR, 'LineStyle', '-', 'Marker', 'none', 'CapSize',0);
            title('Rightward Sweeps'); xlabel('Sweep Pos (deg)'); ylabel('Normalized Peak Turn')
            xlim(ps_limit); ylim(n_limit); xticks(sweepPosR); yline(0); xline(0)
            if motionCheck
                % plot peak turns for leftward sweeps
                nexttile;hold on
                sweepPosL = reshape(sweepPosL, [], 1);
                norm_mean_peak_turnL = mean_peak_turnL ./ normPeakL;
                norm_sem_peak_turnL = sem_peak_turnL ./ normPeakL;
                errorbar(sweepPosL, norm_mean_peak_turnL, norm_sem_peak_turnL, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
                title('Leftward Sweeps'); xlabel('Sweep Pos (deg)');
                xlim(ps_limit); ylim(n_limit); xticks(sweepPosL); yline(0); xline(0)
            end
            % plot peak turns for R+L sweeps
            nexttile;hold on
            sweepPosR = reshape(sweepPosR, [], 1);
            norm_mean_peak_turnRL = mean_peak_turnRL ./ normPeakRL;
            norm_sem_peak_turnRL = sem_peak_turnRL ./ normPeakRL;
            errorbar(sweepPosR, norm_mean_peak_turnRL, norm_sem_peak_turnRL, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
            title('R+L Sweeps'); xlabel('Sweep Pos (deg)');
            xlim(ps_limit); ylim(n_limit); xticks(sweepPosR); yline(0); xline(0)

            sgtitle([strrep([filebase ' ' pulseSpeeds{p} ' ' thisName],'_','/') ' (n = ' num2str(nFliesThresh) ')'])
            % save plot
            cd(folder.summary)
            plotname = strjoin({'pulse_v_turn', thisSpeed, 'dps', 'norm',thisName},'_');
            saveas(gcf,[plotname '.png']);
            copyfile([plotname '.png'], folder.dropbox,'f');
            % save vectorized plot
            cd(folder.vector)
            set(gcf,'renderer','Painters')
            saveas(gcf, [plotname '.svg'])
            copyfile([plotname '.svg'], folder.dropbox,'f');

            % save normalized RL data
            cd(folder.compare)
            dataname = strjoin({filebase, 'turn', num2str(thisSpeed), 'dps', thisName}, '_');
            combinedData = [sweepPosR(:), mean_peak_turnRL(:)./normPeakRL, sem_peak_turnRL(:)./normPeakRL];
            save([dataname '.mat'], 'combinedData');

        end
    end

    %% analyze each motion pulse vs turn response
    disp('Analying motion pulse vs turn response...')
    cd(folder.summary)

    % run (1) all (3) running only
    for r = 3
        % pull run info
        thisName = settings.behaviorGroup{r};
        thisRun = runSelect(r);

        % for each pulse speed
        for p = 1:nPulses
            % initialize
            thisSpeed = pulseSpeeds{p}(1:2);
            pulse_turnR = [];
            pulse_turnL = [];
            pulse_turnRL = [];
            c = 1; %fly counter

            % for each fly
            for nt = 1:nFlies
                if flyRunTime(nt)>settings.minRunTime
                    % select data
                    thisPanelPs = allPanelPs{nt};
                    thisForward = allForward{nt};
                    thisAngular = allAngular{nt};
                    % determine relationship between pulse and behavior
                    [~, thisMean] = pulse_v_onlyturns(thisPanelPs,thisForward,thisAngular,int_time,p,pulseSpeeds,nSweep,thisRun);
                    % if first run, store motion pulse positions
                    if c==1
                        % sweep positions
                        pulse_posR = thisMean.panelpsR;
                        pulse_posL = thisMean.panelpsL;
                        pulseDur = size(pulse_posR,1);
                        % sweep centers for each position
                        sweepPosR = round(pulse_posR(round(pulseDur/2),1,:));
                        sweepPosL = round(pulse_posL(round(pulseDur/2),1,:));
                        % sweep indices
                        sweepIdx = find(~isnan(pulse_posR(:,1)));
                        sweepIdx2End = sweepIdx(1):length(pulse_posR(:,1));
                    end
                    % store turn averages for each fly
                    pulse_turnR(:,c,:) = thisMean.varOutR;
                    pulse_turnL(:,c,:) = thisMean.varOutL;
                    pulse_turnRL(:,c,:) = thisMean.varOutRL;
                    % for each sweep, store peak turning for each fly
                    for s = 1:nSweep
                        if sweepPosR(s) > -settings.binoc % for rightward sweeps on right or binoc, take max right
                            peak_turnR(c,s) = max(thisMean.varOutR(sweepIdx2End,:,s));
                            peak_turnRL(c,s) = max(thisMean.varOutRL(sweepIdx2End,:,s));
                        else % else take max left
                            peak_turnR(c,s) = min(thisMean.varOutR(sweepIdx2End,:,s));
                            peak_turnRL(c,s) = min(thisMean.varOutRL(sweepIdx2End,:,s));
                        end
                        if sweepPosL(s) < settings.binoc % for leftward sweeps on left or binoc, take max left
                            peak_turnL(c,s) = min(thisMean.varOutL(sweepIdx2End,:,s));
                        else % else take max right
                            peak_turnL(c,s) = max(thisMean.varOutL(sweepIdx2End,:,s));
                        end
                    end
                    c = c+1; %update counter
                end
            end
            % calculate turn means
            mean_pulse_turnR = mean(pulse_turnR,2,'omitnan');
            mean_pulse_turnL = mean(pulse_turnL,2,'omitnan');
            mean_pulse_turnRL = mean(pulse_turnRL,2,'omitnan');
            % calculate turn SEMs
            sem_pulse_turnR = std(pulse_turnR,0,2,'omitnan')./sqrt(nFliesThresh);
            sem_pulse_turnL = std(pulse_turnL,0,2,'omitnan')./sqrt(nFliesThresh);
            sem_pulse_turnRL = std(pulse_turnRL,0,2,'omitnan')./sqrt(nFliesThresh);

            % calculate turn mean peaks
            for s = 1:nSweep
                if sweepPosR(s) > -settings.binoc
                    mean_peak_turnR(s) = max(mean_pulse_turnR(sweepIdx2End,:,s));
                    mean_peak_turnRL(s) = max(mean_pulse_turnRL(sweepIdx2End,:,s));
                else
                    mean_peak_turnR(s) = min(mean_pulse_turnR(sweepIdx2End,:,s));
                    mean_peak_turnRL(s) = min(mean_pulse_turnRL(sweepIdx2End,:,s));
                end
                if sweepPosL(s) < settings.binoc
                    mean_peak_turnL(s) = min(mean_pulse_turnL(sweepIdx2End,:,s));
                else
                    mean_peak_turnL(s) = max(mean_pulse_turnL(sweepIdx2End,:,s));
                end
            end
            % calculate turn mean peak SEMs
            sem_peak_turnR = std(peak_turnR,0,1,'omitnan')./sqrt(nFliesThresh);
            sem_peak_turnL = std(peak_turnL,0,1,'omitnan')./sqrt(nFliesThresh);
            sem_peak_turnRL = std(peak_turnRL,0,1,'omitnan')./sqrt(nFliesThresh);

            % plot averages
            % initialize
            figure; set(gcf,'Position',[100 100 plotWidth 800])
            tiledlayout(4,tileX,"TileSpacing","compact")
            % initialize remaining variables
            time_pulse = int_time(1:size(pulse_posR,1))*1000;
            time_limit = [0 max(time_pulse)];

            % plot target position
            % plot rightward target
            nexttile; hold on
            for s = 1:nSweep
                plot(time_pulse,pulse_posR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
            end
            title('Rightward Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(ps_limit); yline(0)
            if motionCheck
                % plot leftward target
                nexttile; hold on
                for s = nSweep:-1:1
                    plot(time_pulse,pulse_posL(:,s),'Color',color_leftward(s,:),'LineWidth',settings.lwAvg)
                end
                title('Leftward Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
                xlim(time_limit); ylim(ps_limit); yline(0)
            end
            % plot combined
            nexttile; hold on
            for s = nSweep:-1:1
                plot(time_pulse,pulse_posR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
            end
            title('R+L Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(ps_limit); yline(0)

            % plot behavior
            % plot rightward turning
            nexttile([3 1]); hold on
            for s = 1:nSweep
                sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_turnR(:,s)-sem_pulse_turnR(:,s); flipud(mean_pulse_turnR(:,s)+sem_pulse_turnR(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                sp.FaceColor = color_rightward(s,:);
                plot(time_pulse,mean_pulse_turnR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
            end
            ylabel('Angular Velocity (deg/sec)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(ang_limit2); yline(0)
            if motionCheck
                % plot leftward turning
                nexttile([3 1]); hold on
                for s = 1:nSweep
                    sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_turnL(:,s)-sem_pulse_turnL(:,s); flipud(mean_pulse_turnL(:,s)+sem_pulse_turnL(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                    sp.FaceColor = color_leftward(s,:);
                    plot(time_pulse,mean_pulse_turnL(:,s),'Color',color_leftward(s,:),'LineWidth',settings.lwAvg)
                end
                ylabel('Angular Velocity (deg/sec)'); xlabel('Time (msec)')
                xlim(time_limit); ylim(ang_limit2); yline(0)
            end
            % plot combined
            nexttile([3 1]); hold on
            for s = 1:nSweep
                sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_turnRL(:,s)-sem_pulse_turnRL(:,s); flipud(mean_pulse_turnRL(:,s)+sem_pulse_turnRL(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                sp.FaceColor = color_rightward(s,:);
                plot(time_pulse,mean_pulse_turnRL(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
            end
            ylabel('Angular Velocity (deg/sec)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(ang_limit2); yline(0)

            sgtitle([strrep([filebase ' ' pulseSpeeds{p} ' ' thisName],'_','/') ' (n = ' num2str(nFliesThresh) ')'])
            % save plot
            cd(folder.summary)
            plotname = strjoin({'pulse_v_onlyturn', thisSpeed, 'dps', thisName},'_');
            saveas(gcf,[plotname '.png']);
            copyfile([plotname '.png'], folder.dropbox,'f');
            % save vectorized plot
            cd(folder.vector)
            set(gcf,'renderer','Painters')
            saveas(gcf, [plotname '.svg'])
            copyfile([plotname '.svg'], folder.dropbox,'f');

            % estimate delay
            try
                thisDelay = pulse_estdelay(pulse_posR,mean_pulse_turnRL,time_pulse);
                estDelay(:,p,r) = thisDelay;
            catch
                estDelay(:,p,r) = nan(nSweep,1);
            end

            % plot normalized peak
            % initialize
            figure; set(gcf,'Position',[100 100 plotWidth 500])
            tiledlayout(1,tileX,"TileSpacing","compact")
            % normalize peak to 1
            normPeakR = max(abs(mean_peak_turnR));
            normPeakL = max(abs(mean_peak_turnL));
            normPeakRL = max(abs(mean_peak_turnRL));
            % plot peak turns for rightward sweeps
            nexttile;hold on
            sweepPosR = reshape(sweepPosR,[],1);
            norm_mean_peak_turnR = mean_peak_turnR ./ normPeakR;
            norm_sem_peak_turnR = sem_peak_turnR ./ normPeakR;
            errorbar(sweepPosR, norm_mean_peak_turnR, norm_sem_peak_turnR, 'LineStyle', '-', 'Marker', 'none', 'CapSize',0);
            title('Rightward Sweeps'); xlabel('Sweep Pos (deg)'); ylabel('Normalized Peak Turn')
            xlim(ps_limit); ylim(n_limit); xticks(sweepPosR); yline(0); xline(0)
            if motionCheck
                % plot peak turns for leftward sweeps
                nexttile;hold on
                sweepPosL = reshape(sweepPosL, [], 1);
                norm_mean_peak_turnL = mean_peak_turnL ./ normPeakL;
                norm_sem_peak_turnL = sem_peak_turnL ./ normPeakL;
                errorbar(sweepPosL, norm_mean_peak_turnL, norm_sem_peak_turnL, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
                title('Leftward Sweeps'); xlabel('Sweep Pos (deg)');
                xlim(ps_limit); ylim(n_limit); xticks(sweepPosL); yline(0); xline(0)
            end
            % plot peak turns for R+L sweeps
            nexttile;hold on
            sweepPosR = reshape(sweepPosR, [], 1);
            norm_mean_peak_turnRL = mean_peak_turnRL ./ normPeakRL;
            norm_sem_peak_turnRL = sem_peak_turnRL ./ normPeakRL;
            errorbar(sweepPosR, norm_mean_peak_turnRL, norm_sem_peak_turnRL, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
            title('R+L Sweeps'); xlabel('Sweep Pos (deg)');
            xlim(ps_limit); ylim(n_limit); xticks(sweepPosR); yline(0); xline(0)

            sgtitle([strrep([filebase ' ' pulseSpeeds{p} ' ' thisName],'_','/') ' (n = ' num2str(nFliesThresh) ')'])
            % save plot
            cd(folder.summary)
            plotname = strjoin({'pulse_v_onlyturn', thisSpeed, 'dps', 'norm',thisName},'_');
            saveas(gcf,[plotname '.png']);
            copyfile([plotname '.png'], folder.dropbox,'f');
            % save vectorized plot
            cd(folder.vector)
            set(gcf,'renderer','Painters')
            saveas(gcf, [plotname '.svg'])
            copyfile([plotname '.svg'], folder.dropbox,'f');

            % save normalized RL data
            cd(folder.compare)
            dataname = strjoin({filebase, 'onlyturn', num2str(thisSpeed), 'dps', thisName}, '_');
            combinedData = [sweepPosR(:), mean_peak_turnRL(:)./normPeakRL, sem_peak_turnRL(:)./normPeakRL];
            save([dataname '.mat'], 'combinedData');

        end
    end

    %% analyze each motion pulse vs forward response
    disp('Analying motion pulse vs forward response...')
    cd(folder.summary)
    fwd_limit = [0 10];

    % run (1) all (3) running only
    for r = [1 3]
        % pull run info
        thisName = settings.behaviorGroup{r};
        thisRun = runSelect(r);

        % for each pulse speed
        for p = 1:nPulses
            % initialize
            thisSpeed = pulseSpeeds{p}(1:2);
            pulse_forwardR = [];
            pulse_forwardL = [];
            pulse_forwardRL = [];
            c = 1; %fly counter

            % for each fly
            for nt = 1:nFlies
                if flyRunTime(nt)>settings.minRunTime
                    % select data
                    thisPanelPs = allPanelPs{nt};
                    thisForward = allForward{nt};
                    % determine relationship between pulse and behavior
                    [~, thisMean] = pulse_v_output(thisPanelPs,thisForward,thisForward,int_time,p,pulseSpeeds,nSweep,thisRun);
                    % if first run, store motion pulse positions
                    if c==1
                        % sweep positions
                        pulse_posR = thisMean.panelpsR;
                        pulse_posL = thisMean.panelpsL;
                        pulseDur = size(pulse_posR,1);
                        % sweep centers for each position
                        sweepPosR = round(pulse_posR(round(pulseDur/2),1,:));
                        sweepPosL = round(pulse_posL(round(pulseDur/2),1,:));
                        % sweep indices
                        sweepIdx = find(~isnan(pulse_posR(:,1)));
                        sweepIdx2End = sweepIdx(1):length(pulse_posR(:,1));
                    end
                    % store forward averages for each fly
                    pulse_forwardR(:,c,:) = thisMean.varOutR;
                    pulse_forwardL(:,c,:) = thisMean.varOutL;
                    pulse_forwardRL(:,c,:) = mean([thisMean.varOutR, flip(thisMean.varOutL,3)],2);
                    % for each sweep, store peak forwarding for each fly
                    for s = 1:nSweep
                        peak_forwardR(c,s) = max(thisMean.varOutR(sweepIdx2End,:,s));
                        peak_forwardRL(c,s) = max(thisMean.varOutRL(sweepIdx2End,:,s));
                        peak_forwardL(c,s) = max(thisMean.varOutL(sweepIdx2End,:,s));
                    end
                    c = c+1; %update counter
                end
            end
            % calculate forward means
            mean_pulse_forwardR = mean(pulse_forwardR,2,'omitnan');
            mean_pulse_forwardL = mean(pulse_forwardL,2,'omitnan');
            mean_pulse_forwardRL = mean(pulse_forwardRL,2,'omitnan');
            % calculate forward SEMs
            sem_pulse_forwardR = std(pulse_forwardR,0,2,'omitnan')./sqrt(nFliesThresh);
            sem_pulse_forwardL = std(pulse_forwardL,0,2,'omitnan')./sqrt(nFliesThresh);
            sem_pulse_forwardRL = std(pulse_forwardRL,0,2,'omitnan')./sqrt(nFliesThresh);

            % calculate forward mean peaks
            for s = 1:nSweep
                mean_peak_forwardR(s) = max(mean_pulse_forwardR(sweepIdx2End,:,s));
                mean_peak_forwardRL(s) = max(mean_pulse_forwardRL(sweepIdx2End,:,s));
                mean_peak_forwardL(s) = max(mean_pulse_forwardL(sweepIdx2End,:,s));
            end
            % calculate forward mean peak SEMs
            sem_peak_forwardR = std(peak_forwardR,0,1,'omitnan')./sqrt(nFliesThresh);
            sem_peak_forwardL = std(peak_forwardL,0,1,'omitnan')./sqrt(nFliesThresh);
            sem_peak_forwardRL = std(peak_forwardRL,0,1,'omitnan')./sqrt(nFliesThresh);

            % plot averages
            % initialize
            figure; set(gcf,'Position',[100 100 plotWidth 800])
            tiledlayout(4,tileX,"TileSpacing","compact")
            % initialize remaining variables
            time_pulse = int_time(1:size(pulse_posR,1))*1000;
            time_limit = [0 max(time_pulse)];

            % plot target position
            % plot rightward target
            nexttile; hold on
            for s = 1:nSweep
                plot(time_pulse,pulse_posR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
            end
            title('Rightward Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(ps_limit); yline(0)
            if motionCheck
                % plot leftward target
                nexttile; hold on
                for s = nSweep:-1:1
                    plot(time_pulse,pulse_posL(:,s),'Color',color_leftward(s,:),'LineWidth',settings.lwAvg)
                end
                title('Leftward Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
                xlim(time_limit); ylim(ps_limit); yline(0)
            end
            % plot combined
            nexttile; hold on
            for s = nSweep:-1:1
                plot(time_pulse,pulse_posR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
            end
            title('R+L Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(ps_limit); yline(0)

            % plot behavior
            % plot rightward forwarding
            nexttile([3 1]); hold on
            for s = 1:nSweep
                sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_forwardR(:,s)-sem_pulse_forwardR(:,s); flipud(mean_pulse_forwardR(:,s)+sem_pulse_forwardR(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                sp.FaceColor = color_rightward(s,:);
                plot(time_pulse,mean_pulse_forwardR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
            end
            ylabel('Forward Velocity (deg/sec)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(fwd_limit); yline(0)
            if motionCheck
                % plot leftward forwarding
                nexttile([3 1]); hold on
                for s = 1:nSweep
                    sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_forwardL(:,s)-sem_pulse_forwardL(:,s); flipud(mean_pulse_forwardL(:,s)+sem_pulse_forwardL(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                    sp.FaceColor = color_leftward(s,:);
                    plot(time_pulse,mean_pulse_forwardL(:,s),'Color',color_leftward(s,:),'LineWidth',settings.lwAvg)
                end
                ylabel('Forward Velocity (deg/sec)'); xlabel('Time (msec)')
                xlim(time_limit); ylim(fwd_limit); yline(0)
            end
            % plot combined
            nexttile([3 1]); hold on
            for s = 1:nSweep
                sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_forwardRL(:,s)-sem_pulse_forwardRL(:,s); flipud(mean_pulse_forwardRL(:,s)+sem_pulse_forwardRL(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                sp.FaceColor = color_rightward(s,:);
                plot(time_pulse,mean_pulse_forwardRL(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
            end
            ylabel('Forward Velocity (deg/sec)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(fwd_limit); yline(0)

            sgtitle([strrep([filebase ' ' pulseSpeeds{p} ' ' thisName],'_','/') ' (n = ' num2str(nFliesThresh) ')'])
            % save plot
            cd(folder.summary)
            plotname = strjoin({'pulse_v_forward', thisSpeed, 'dps', thisName},'_');
            saveas(gcf,[plotname '.png']);
            copyfile([plotname '.png'], folder.dropbox,'f');
            % save vectorized plot
            cd(folder.vector)
            set(gcf,'renderer','Painters')
            saveas(gcf, [plotname '.svg'])
            copyfile([plotname '.svg'], folder.dropbox,'f');

            % estimate delay
            try
                thisDelay = pulse_estdelay(pulse_posR,mean_pulse_forwardRL,time_pulse);
                estDelay(:,p,r) = thisDelay;
            catch
                estDelay(:,p,r) = nan(nSweep,1);
            end

            % plot normalized peak
            % initialize
            figure; set(gcf,'Position',[100 100 plotWidth 500])
            tiledlayout(1,tileX,"TileSpacing","compact")
            % normalize peak to 1
            normPeakR = max(abs(mean_peak_forwardR));
            normPeakL = max(abs(mean_peak_forwardL));
            normPeakRL = max(abs(mean_peak_forwardRL));
            % plot peak forwards for rightward sweeps
            nexttile;hold on
            sweepPosR = reshape(sweepPosR, [], 1);
            norm_mean_peak_forwardR = mean_peak_forwardR ./ normPeakR;
            norm_sem_peak_forwardR = sem_peak_forwardR ./ normPeakR;
            errorbar(sweepPosR, norm_mean_peak_forwardR, norm_sem_peak_forwardR, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
            title('Rightward Sweeps'); xlabel('Sweep Pos (deg)'); ylabel('Normalized Peak forward')
            xlim(ps_limit); ylim(n_limit); xticks(sweepPosR); yline(0); xline(0)
            if motionCheck
                % plot peak forwards for leftward sweeps
                nexttile;hold on
                sweepPosL = reshape(sweepPosL, [], 1);
                norm_mean_peak_forwardL = mean_peak_forwardL ./ normPeakL;
                norm_sem_peak_forwardL = sem_peak_forwardL ./ normPeakL;
                errorbar(sweepPosL, norm_mean_peak_forwardL, norm_sem_peak_forwardL, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
                title('Leftward Sweeps'); xlabel('Sweep Pos (deg)');
                xlim(ps_limit); ylim(n_limit); xticks(sweepPosL); yline(0); xline(0)
            end
            % plot peak forwards for R+L sweeps
            nexttile;hold on
            sweepPosR = reshape(sweepPosR, [], 1);
            norm_mean_peak_forwardRL = mean_peak_forwardRL ./ normPeakRL;
            norm_sem_peak_forwardRL = sem_peak_forwardRL ./ normPeakRL;
            errorbar(sweepPosR, norm_mean_peak_forwardRL, norm_sem_peak_forwardRL, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
            title('R+L Sweeps'); xlabel('Sweep Pos (deg)');
            xlim(ps_limit); ylim(n_limit); xticks(sweepPosR); yline(0); xline(0)

            sgtitle([strrep([filebase ' ' pulseSpeeds{p} ' ' thisName],'_','/') ' (n = ' num2str(nFliesThresh) ')'])
            % save plot
            cd(folder.summary)
            plotname = strjoin({'pulse_v_forward', thisSpeed, 'dps', 'norm',thisName},'_');
            saveas(gcf,[plotname '.png']);
            copyfile([plotname '.png'], folder.dropbox,'f');
            % save vectorized plot
            cd(folder.vector)
            set(gcf,'renderer','Painters')
            saveas(gcf, [plotname '.svg'])
            copyfile([plotname '.svg'], folder.dropbox,'f');

            % save normalized RL data
            cd(folder.compare)
            dataname = strjoin({filebase, 'forward', num2str(thisSpeed), 'dps', thisName}, '_');
            combinedData = [sweepPosR(:), mean_peak_forwardRL(:)./normPeakRL, sem_peak_forwardRL(:)./normPeakRL];
            actualData = [sweepPosR(:), mean_peak_turnRL(:), sem_peak_turnRL(:)];
            save([dataname '.mat'], 'combinedData', 'actualData');

        end
    end
end

%% analyze each motion pulse vs spikerate response
disp('Analying motion pulse vs spikerate response...')
cd(folder.summary)
store_ds = []; store_all75_ds = []; store_all25_ds = [];
close all
peak_srR_store = [];
ds_limit = [-1.2 1.2];

% run (1) all, (2) quiescent only, (3) running only
for r = 2
    % pull run info
    thisName = settings.behaviorGroup{r};
    thisRun = runSelect(r);
    % initialize
    mean_pulse_srR = [];
    mean_pulse_srL = [];
    sem_pulse_srR = [];
    sem_pulse_srL = [];

    % for each pulse speed
    for p = 1:nPulses
        % initialize
        thisSpeed = pulseSpeeds{p}(1:2);
        pulse_srRightward = [];
        pulse_srLeftward = [];

        % for each fly
        for nt = 1:nFlies
            % select data
            thisPanelPs = allPanelPs{nt};
            thisForward = allForward{nt};
            thisSpikert = allSpikeRt{nt};
            % determine relationship between pulse and spikerate
            [~, thisMean] = pulse_v_output(thisPanelPs,thisForward,thisSpikert,int_time,p,pulseSpeeds,nSweep,thisRun);
            % if first run, store motion pulse positions
            if nt == 1
                % sweep positions
                pulse_posR(:,p,:) = thisMean.panelpsR;
                pulse_posL(:,p,:) = thisMean.panelpsL;
                pulseDur = size(pulse_posR,1);
                time_pulse = int_time(1:size(pulse_posR,1))*1000;
                % sweep centers for each position
                sweepPosR = round(pulse_posR(round(pulseDur/2),1,:));
                sweepPosL = round(pulse_posL(round(pulseDur/2),1,:));
                % sweep indices
                sweepIdx = find(~isnan(pulse_posR(:,p,1)));
                sweepIdx2End = sweepIdx(1):length(pulse_posR(:,p,1));
            end
            % store spikerate averages for each fly
            pulse_srRightward(:,nt,:) = thisMean.varOutR;
            pulse_srLeftward(:,nt,:) = thisMean.varOutL;
            pulse_srRL_Rightward(:,nt,:) = thisMean.varOutR - flip(thisMean.varOutL,3);
            pulse_srRL_Leftward(:,nt,:) = thisMean.varOutL - flip(thisMean.varOutR,3);
            % for each sweep, store peak spikerate for each fly
            for s = 1:nSweep
                avg_srRightward(nt,s) = mean(pulse_srRightward(sweepIdx,nt,s),'omitnan');
                avg_srLeftward(nt,s) = mean(pulse_srLeftward(sweepIdx,nt,s),'omitnan');
                avg_srRL_Rightward(nt,s) = mean(pulse_srRL_Rightward(sweepIdx,nt,s),'omitnan');
                avg_srRL_Leftward(nt,s) = mean(pulse_srRL_Leftward(sweepIdx,nt,s),'omitnan');
            end
        end
        % fetch adjusted n
        nFlies_adj_station = sum(~isnan(pulse_srRL_Rightward(1,:,1)));

        % calculate response latency
        latency_times = calculate_visualresponselatency(pulse_srRightward, pulse_posR(:,p,:), time_pulse);
        mean_latency = mean(latency_times,'omitnan');
        sem_latency = std(latency_times,'omitnan')./sqrt(nFlies_adj_station);

        % calculate spikerate means and sem (per timepoint)
        mean_pulse_srR(:,p,:) = mean(pulse_srRightward, 2, 'omitnan');
        mean_pulse_srL(:,p,:) = mean(pulse_srLeftward, 2, 'omitnan');
        mean_pulse_srRL(:,p,:) = mean_pulse_srR(:,p,:) - flip(mean_pulse_srL(:,p,:), 3);

        % calculate SEM (per timepoint)
        sem_pulse_srR(:,p,:) = std(pulse_srRightward, 0, 2, 'omitnan') ./ sqrt(nFlies_adj_station);
        sem_pulse_srL(:,p,:) = std(pulse_srLeftward, 0, 2, 'omitnan') ./ sqrt(nFlies_adj_station);
        sem_pulse_srRL(:,p,:) = sqrt(sem_pulse_srR(:,p,:).^2 + flip(sem_pulse_srL(:,p,:).^2, 3));

        % (Optional) If you still want LR:
        mean_pulse_srLR(:,p,:) = mean_pulse_srL(:,p,:) - flip(mean_pulse_srR(:,p,:), 3);
        sem_pulse_srLR(:,p,:) = sqrt(sem_pulse_srL(:,p,:).^2 + flip(sem_pulse_srR(:,p,:).^2, 3));

        % calculate spikerate peak means (across time)
        mean_peak_srR = mean(avg_srRightward, 1, 'omitnan');
        mean_peak_srL = mean(avg_srLeftward, 1, 'omitnan');
        mean_peak_srRL_rightward = mean_peak_srR - mean(avg_srLeftward(:,nSweep:-1:1), 1, 'omitnan');

        % calculate SEM for peak values
        sem_peak_srR = std(avg_srRightward, 0, 1, 'omitnan') ./ sqrt(nFlies_adj_station);
        sem_peak_srL = std(avg_srLeftward, 0, 1, 'omitnan') ./ sqrt(nFlies_adj_station);
        sem_peak_srRL = sqrt(sem_peak_srR.^2 + std(avg_srLeftward(:,nSweep:-1:1), 0, 1, 'omitnan').^2 ./ nFlies_adj_station);

        % Optional: LR version
        mean_peak_srRL_leftward = mean_peak_srL - mean(avg_srRightward(:,nSweep:-1:1), 1, 'omitnan');
        sem_peak_srLR = sqrt(sem_peak_srL.^2 + std(avg_srRightward(:,nSweep:-1:1), 0, 1, 'omitnan').^2 ./ nFlies_adj_station);

        % Calculate baseline subtraction for pulse_srR and pulse_srL
        animal_baseline_srR = mean(pulse_srRightward(1:sweepIdx(1),:,:), 1, 'omitnan');  % Baseline across time per animal for pulse_srR
        animal_baseline_srL = mean(pulse_srLeftward(1:sweepIdx(1),:,:), 1, 'omitnan');  % Baseline across time per animal for pulse_srL

        pulse_srR_baseline_subtracted = pulse_srRightward - animal_baseline_srR;  % Subtract baseline from pulse_srR
        pulse_srL_baseline_subtracted = pulse_srLeftward - animal_baseline_srL;  % Subtract baseline from pulse_srL

        % Calculate baseline subtraction for peak_srR and peak_srL
        animal_baseline_peak_srR = mean(avg_srRightward, 1, 'omitnan');  % Baseline across animals for peak_srR
        animal_baseline_peak_srL = mean(avg_srLeftward, 1, 'omitnan');  % Baseline across animals for peak_srL

        peak_srR_baseline_subtracted = avg_srRightward - animal_baseline_peak_srR;  % Subtract baseline from peak_srR
        peak_srL_baseline_subtracted = avg_srLeftward - animal_baseline_peak_srL;  % Subtract baseline from peak_srL

        % Calculate direction selectivity using peak
        ds_sr_all = mean((pulse_srR_baseline_subtracted(sweepIdx,:,1:nSweep-1) - pulse_srL_baseline_subtracted(sweepIdx,:,2:nSweep)) ./ ...
            (pulse_srR_baseline_subtracted(sweepIdx,:,1:nSweep-1) + pulse_srL_baseline_subtracted(sweepIdx,:,2:nSweep)), 'omitnan');
        % Calculate direction selectivity using peak
        ds_sr_peak = (peak_srR_baseline_subtracted(:,1:nSweep-1) - peak_srL_baseline_subtracted(:,2:nSweep)) ./ ...
            (peak_srR_baseline_subtracted(:,1:nSweep-1) + peak_srL_baseline_subtracted(:,2:nSweep));

        mean_ds_sr_all = mean(ds_sr_all,2,'omitnan');
        sem_ds_sr_all = std(ds_sr_all,0,2,'omitnan')./sqrt(nFlies);
        mean_ds_sr_peak = mean(ds_sr_peak,1,'omitnan');
        sem_ds_sr_peak = std(ds_sr_peak,0,1,'omitnan')./sqrt(nFlies);

        % For ipsi sweeps, fetch DSI for saving later
        ipsiIdx = sweepPosR(1:end-1)>0;
        ipsiIdx2 = sweepPosR(1:end-1)>-15;
        if p == 1
            store_ds(1,:,r) = median(ds_sr_all(:,:,ipsiIdx),3);
            store_all75_ds(:,:,r) = squeeze(ds_sr_all(:,:,ipsiIdx2));
        else
            store_ds(2,:,r) = median(ds_sr_all(:,:,ipsiIdx),3);
            store_all25_ds(:,:,r) = squeeze(ds_sr_all(:,:,ipsiIdx2));
        end

        % plot averages for each speed
        % initialize
        figure; set(gcf,'Position',[100 100 plotWidth 800])
        tiledlayout(4,tileX,"TileSpacing","compact")
        % initialize remaining variables
        time_pulse = int_time(1:size(pulse_posR,1))*1000;
        time_limit = [0 max(time_pulse)];

        % plot target position
        % plot rightward target
        nexttile; hold on
        for s = 1:nSweep
            plot(time_pulse,pulse_posR(:,p,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
        end
        title('Rightward Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(ps_limit); yline(0)
        if motionCheck
            % plot leftward target
            nexttile; hold on
            for s = nSweep:-1:1
                plot(time_pulse,pulse_posL(:,p,s),'Color',color_leftward(s,:),'LineWidth',settings.lwAvg)
            end
            title('Leftward Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(ps_limit); yline(0)
        end
        % plot R-L
        nexttile; hold on
        for s = 1:nSweep
            plot(time_pulse,pulse_posR(:,p,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
        end
        title('R-L'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(ps_limit); yline(0)

        % plot spikerate
        % plot rightward turning
        nexttile([3 1]); hold on
        for s = 1:nSweep
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srR(:,p,s)-sem_pulse_srR(:,p,s); flipud(mean_pulse_srR(:,p,s)+sem_pulse_srR(:,p,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = color_rightward(s,:);
            plot(time_pulse,mean_pulse_srR(:,p,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
        end
        % Calculate the string to display
        latency_text = sprintf('%.2f \\pm %.2f', mean_latency, sem_latency);
        text(time_limit(2), sr_limitA(2), latency_text,'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 6);
        ylabel('Firing Rate (spikes/s)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(sr_limitA); yline(0)
        if motionCheck
            % plot leftward spikerate
            nexttile([3 1]); hold on
            for s = 1:nSweep
                sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srL(:,p,s)-sem_pulse_srL(:,p,s); flipud(mean_pulse_srL(:,p,s)+sem_pulse_srL(:,p,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                sp.FaceColor = color_leftward(s,:);
                plot(time_pulse,mean_pulse_srL(:,p,s),'Color',color_leftward(s,:),'LineWidth',settings.lwAvg)
            end
            ylabel('Firing Rate (spikes/s)'); xlabel('Time (msec)')
            xlim(time_limit); ylim(sr_limitA); yline(0)
        end
        % plot R-L
        nexttile([3 1]); hold on
        for s = 1:nSweep
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srRL(:,p,s)-sem_pulse_srRL(:,p,s); flipud(mean_pulse_srRL(:,p,s)+sem_pulse_srRL(:,p,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = color_rightward(s,:);
            plot(time_pulse,mean_pulse_srRL(:,p,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
        end
        ylabel('Expected Firing Rate Difference (spikes/s)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(sr_limitD); yline(0)

        sgtitle([strrep([filebase ' ' pulseSpeeds{p} ' ' thisName],'_','/') ' (n = ' num2str(nFlies_adj_station) ')'])
        % save plot
        cd(folder.summary)
        plotname = strjoin({'pulse_v_fr', thisSpeed, 'dps', thisName},'_');
        saveas(gcf,[plotname '.png']);
        copyfile([plotname '.png'], folder.dropbox,'f');
        % save vectorized plot
        cd(folder.vector)
        set(gcf,'renderer','Painters')
        saveas(gcf, [plotname '.svg'])
        copyfile([plotname '.svg'], folder.dropbox,'f');

        % initialize
        figure; set(gcf,'Position',[100 100 900 800])
        tiledlayout(2,4,"TileSpacing","compact")
        % plot raw peaks
        % plot peak fr for rightward sweeps
        nexttile;hold on
        sweepPosR = reshape(sweepPosR, [], 1);
        errorbar(sweepPosR, mean_peak_srR, sem_peak_srR, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
        title('Rightward Sweeps'); xlabel('Sweep Pos (deg)'); ylabel('Average FR')
        xlim(ps_limit); ylim(sr_limitA); xticks(sweepPosR); yline(0); xline(0)
        if motionCheck
            % plot peak fr for leftward sweeps
            nexttile;hold on
            sweepPosL = reshape(sweepPosL, [], 1);
            errorbar(sweepPosL, mean_peak_srL, sem_peak_srL, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
            title('Leftward Sweeps'); xlabel('Sweep Pos (deg)'); ylabel('Average FR')
            xlim(ps_limit); ylim(sr_limitA); xticks(sweepPosL); yline(0); xline(0)
        end
        % plot peak fr for R-L sweeps
        nexttile;hold on
        sweepPosR = reshape(sweepPosR, [], 1);
        errorbar(sweepPosR, mean_peak_srRL_rightward, sem_peak_srRL, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
        if motionCheck
            errorbar(sweepPosL, mean_peak_srRL_leftward, sem_peak_srLR, 'k', 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
        end
        title('R-L Rightward, Leftward'); xlabel('Sweep Pos (deg)'); ylabel('Average FR')
        xlim(ps_limit); ylim([-60 60]); xticks(sweepPosR); yline(0); xline(0)
        if motionCheck
            % plot direction selectivity index
            nexttile; hold on
            sweepPosR = reshape(sweepPosR, [], 1);
            mean_ds_sr_all = reshape(mean_ds_sr_all, [], 1);
            sem_ds_sr_all = reshape(sem_ds_sr_all, [], 1);
            errorbar(sweepPosR(1:end-1), mean_ds_sr_all, sem_ds_sr_all, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
            title('Sweep DS'); xlabel('Sweep Pos (deg)');
            xlim(ps_limit); ylim(ds_limit); xticks(sweepPosL); yline(0)
        end

        % plot normalized peak
        % normalize peak to 1
        normPeakR = max(abs(mean_peak_srR));
        normPeakL = max(abs(mean_peak_srL));
        normPeakRL = max(abs(mean_peak_srRL_rightward));
        % plot peak fr for rightward sweeps
        nexttile;hold on
        sweepPosR = reshape(sweepPosR, [], 1);
        norm_mean_peak_srR = mean_peak_srR ./ normPeakR;
        norm_sem_peak_srR = sem_peak_srR ./ normPeakR;
        errorbar(sweepPosR, norm_mean_peak_srR, norm_sem_peak_srR, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
        title('Rightward Sweeps'); xlabel('Sweep Pos (deg)'); ylabel('Normalized Peak FR')
        xlim(ps_limit); ylim(n_limitD); xticks(sweepPosR); yline(0); xline(0)
        if motionCheck
            % plot peak fr for leftward sweeps
            nexttile;hold on
            sweepPosL = reshape(sweepPosL, [], 1);
            norm_mean_peak_srL = mean_peak_srL ./ normPeakL;
            norm_sem_peak_srL = sem_peak_srL ./ normPeakL;
            errorbar(sweepPosL, norm_mean_peak_srL, norm_sem_peak_srL, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
            title('Leftward Sweeps'); xlabel('Sweep Pos (deg)'); ylabel('Normalized Peak FR')
            xlim(ps_limit); ylim(n_limitD); xticks(sweepPosL); yline(0); xline(0)
        end
        % plot peak fr for R-L sweeps
        nexttile;hold on
        sweepPosR = reshape(sweepPosR, [], 1);
        norm_mean_peak_srRL = mean_peak_srRL_rightward ./ normPeakRL;
        norm_sem_peak_srRL = sem_peak_srRL ./ normPeakRL;
        errorbar(sweepPosR, norm_mean_peak_srRL, norm_sem_peak_srRL, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
        title('R-L Sweeps'); xlabel('Sweep Pos (deg)'); ylabel('Normalized Peak FR')
        xlim(ps_limit); ylim(n_limit); xticks(sweepPosR); yline(0); xline(0)
        if motionCheck
            % plot direction selectivity index
            nexttile; hold on
            sweepPosR = reshape(sweepPosR, [], 1);
            errorbar(sweepPosR(1:end-1), mean_ds_sr_peak, sem_ds_sr_peak, 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
            title('Peak DS'); xlabel('Sweep Pos (deg)');
            xlim(ps_limit); ylim(ds_limit); xticks(sweepPosL); yline(0)
        end

        sgtitle(strrep([filebase ' ' pulseSpeeds{p} ' ' thisName],'_','/'))
        % save plot
        cd(folder.summary)
        plotname = strjoin({'pulse_v_fr', thisSpeed, 'dps', 'norm',thisName},'_');
        saveas(gcf,[plotname '.png']);
        copyfile([plotname '.png'], folder.dropbox,'f');
        % save vectorized plot
        cd(folder.vector)
        set(gcf,'renderer','Painters')
        saveas(gcf, [plotname '.svg'])
        copyfile([plotname '.svg'], folder.dropbox,'f');

        % store peak data
        if r==2
            peak_srR_store(:,:,p) = avg_srRightward;
        end
        norm_peakmean_RL(p,:) = norm_mean_peak_srR;
        norm_peaksem_RL(p,:) = norm_sem_peak_srR;
    end

    % plot pulse speeds together
    % initialize
    figure; set(gcf,'Position',[100 100 1500 600])
    tiledlayout(6,nSweep,"TileSpacing","compact")
    % plot rightward target positions overlayed
    for s = 1:nSweep
        nexttile; hold on
        for p = 1:nPulses
            plot(time_pulse,pulse_posR(:,p,s),'Color',settings.mopColor{p},'LineWidth',settings.lwAvg)
        end
        if s==1
            ylabel('Target Pos (deg)');
        end
        xlim(time_limit); ylim(ps_limit); yline(0)
    end
    % plot rightward firing rate averages overlayed
    for s = 1:nSweep
        nexttile([2 1]); hold on
        for p = 1:nPulses
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srR(:,p,s)-sem_pulse_srR(:,p,s); flipud(mean_pulse_srR(:,p,s)+sem_pulse_srR(:,p,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = settings.mopColor{p};
            plot(time_pulse,mean_pulse_srR(:,p,s),'Color',settings.mopColor{p},'LineWidth',settings.lwAvg)
        end
        if s==1
            ylabel('Firing Rate (spikes/s)');
        end
        xlim(time_limit); ylim(sr_limitA); yline(0)
    end
    if motionCheck
        % plot leftward target positions overlayed
        for s = 1:nSweep
            nexttile; hold on
            for p = 1:nPulses
                plot(time_pulse,pulse_posL(:,p,s),'Color',settings.mopColor{p},'LineWidth',settings.lwAvg)
            end
            if s==1
                ylabel('Target Pos (deg)');
            end
            xlim(time_limit); ylim(ps_limit); yline(0)
        end
        % plot leftward firing rate averages overlayed
        for s = 1:nSweep
            nexttile([2 1]); hold on
            for p = 1:nPulses
                sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srL(:,p,s)-sem_pulse_srL(:,p,s); flipud(mean_pulse_srL(:,p,s)+sem_pulse_srL(:,p,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                sp.FaceColor = settings.mopColor{p};
                plot(time_pulse,mean_pulse_srL(:,p,s),'Color',settings.mopColor{p},'LineWidth',settings.lwAvg)
            end
            if s==1
                ylabel('Firing Rate (spikes/s)');
            end
            xlim(time_limit); ylim(sr_limitA); yline(0)
        end
    end
    sgtitle(strrep([filebase ' ' thisName],'_','/'))
    % save plot
    cd(folder.summary)
    plotname = strjoin({'pulse_v_fr_overlay' thisName},'_');
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

    % plot R-L
    if motionCheck
        % initialize
        figure; set(gcf,'Position',[100 100 1500 600])
        tiledlayout(6,nSweep,"TileSpacing","compact")
        % plot rightward target positions overlayed
        for s = 1:nSweep
            nexttile; hold on
            for p = 1:nPulses
                plot(time_pulse,pulse_posR(:,p,s),'Color',settings.mopColor{p},'LineWidth',settings.lwAvg)
            end
            if s==1
                ylabel('Target Pos (deg)');
            end
            xlim(time_limit); ylim(ps_limit); yline(0)
        end
        % plot R-L firing rate averages overlayed
        for s = 1:nSweep
            nexttile([2 1]); hold on
            for p = 1:nPulses
                sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srRL(:,p,s)-sem_pulse_srRL(:,p,s); flipud(mean_pulse_srRL(:,p,s)+sem_pulse_srRL(:,p,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                sp.FaceColor = settings.mopColor{p};
                plot(time_pulse,mean_pulse_srRL(:,p,s),'Color',settings.mopColor{p},'LineWidth',settings.lwAvg)
            end
            if s==1
                ylabel('R-L Firing Rate (spikes/s)');
            end
            xlim(time_limit); ylim(sr_limitD); yline(0)
        end
        % plot leftward target positions overlayed
        for s = 1:nSweep
            nexttile; hold on
            for p = 1:nPulses
                plot(time_pulse,pulse_posL(:,p,s),'Color',settings.mopColor{p},'LineWidth',settings.lwAvg)
            end
            if s==1
                ylabel('Target Pos (deg)');
            end
            xlim(time_limit); ylim(ps_limit); yline(0)
        end
        % plot L-R firing rate averages overlayed
        for s = 1:nSweep
            nexttile([2 1]); hold on
            for p = 1:nPulses
                sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srLR(:,p,s)-sem_pulse_srLR(:,p,s); flipud(mean_pulse_srLR(:,p,s)+sem_pulse_srLR(:,p,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                sp.FaceColor = settings.mopColor{p};
                plot(time_pulse,mean_pulse_srLR(:,p,s),'Color',settings.mopColor{p},'LineWidth',settings.lwAvg)
            end
            if s==1
                ylabel('R-L Firing Rate (spikes/s)');
            end
            xlim(time_limit); ylim(sr_limitD); yline(0)
        end
        sgtitle(strrep([filebase ' ' thisName ' R-L'],'_','/'))
        % save plot
        cd(folder.summary)
        plotname = strjoin({'pulse_v_frRL_overlay' thisName},'_');
        saveas(gcf,[plotname '.png']);
        copyfile([plotname '.png'], folder.dropbox,'f');
        % save vectorized plot
        cd(folder.vector)
        set(gcf,'renderer','Painters')
        saveas(gcf, [plotname '.svg'])
        copyfile([plotname '.svg'], folder.dropbox,'f');

        % save RL data
        cd(folder.compare)
        dataname = strjoin({filebase, 'RL', thisName}, '_');
        save([dataname '.mat'], 'pulse_posR', 'mean_pulse_srRL','sem_pulse_srRL','mean_peak_srRL_rightward','sem_peak_srRL','time_pulse');
    end

    % save normalized RL data
    cd(folder.compare)
    dataname = strjoin({filebase, 'fr', num2str(thisSpeed), 'dps', thisName}, '_');
    combinedData = [sweepPosR(:), mean_peak_srRL_rightward(:)./normPeakRL, sem_peak_srRL(:)./normPeakRL];
    save([dataname '.mat'], 'combinedData');

    cd(folder.compare)
    dataname = strjoin({filebase, 'sweeppeaks', num2str(thisSpeed), 'dps', thisName}, '_');
    save([dataname '.mat'], 'flyShortNames', 'sweepPosR','pulse_srRightward','pulse_srRL_Rightward','avg_srRightward');

end

if motionCheck
    analyze_peak_srR(folder, peak_srR_store)
    [~, stats, ~] = run_dir_sweep_lme(avg_srRL_Rightward, avg_srRL_Leftward)
end

% save direction selectivity data
cd(folder.compare)
dataname = strjoin({'ds', filebase}, '_');
save([dataname '.mat'], 'store_ds','store_all75_ds','store_all25_ds');

%% Model FR to turning
if nt_t>1
    % Set model parameters
    k = 50;      % gain — scale to ~200 deg/s
    f = 40;       % friction offset
    w = 100;       % ~100 ms window at ~889 Hz sampling rate

    [r_model, local_auc] = model_rotvel_from_dna02(mean_pulse_srRL, time_pulse, k, f, w);


    % Setup
    [~, nPulses, nSweep] = size(mean_pulse_srRL);
    figure;set(gcf,'Position',[100 100 1500 600])
    tiledlayout(3, nSweep, 'TileSpacing', 'compact');

    % --- Row 1: Target position ---
    for s = 1:nSweep
        nexttile(s); hold on
        for p = 1:nPulses
            plot(time_pulse, pulse_posR(:,p,s), ...
                'Color', 'k', 'LineWidth', settings.lwAvg);
        end
        if s == 1
            ylabel('Target Pos (deg)');
        end
        xlim(time_limit); ylim(ps_limit); yline(0);
    end

    % --- Row 2: R-L firing rate ---
    for s = 1:nSweep
        nexttile(nSweep + s); hold on
        for p = 1:nPulses
            m = mean_pulse_srRL(:,p,s);
            e = sem_pulse_srRL(:,p,s);
            patch([time_pulse(:); flipud(time_pulse(:))], ...
                [m - e; flipud(m + e)], ...
                'k', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
            plot(time_pulse, m, 'Color', settings.mopColor{p}, 'LineWidth', settings.lwAvg);
        end
        if s == 1
            ylabel('R-L FR (sp/s)');
        end
        xlim(time_limit); ylim(sr_limitD); yline(0);
    end

    % --- Row 3: Modeled rotVel ---
    for s = 1:nSweep
        nexttile(2*nSweep + s); hold on
        plot(time_pulse, r_model(:,s), 'k', 'LineWidth', 1.5);
        if s == 1
            ylabel('Model rotVel');
        end
        xlim(time_limit); yline(0); ylim([-20 200])
    end

    % save plot
    cd(folder.summary)
    plotname = 'pulse_torque';
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

    %% analyze spike rate v directional velocity
    disp('Analyzing spike rate v directional velocity...')

    figure; set(gcf,'Position',[100 100 800 800])
    for v = 1:3
        switch v
            case 1 %fwd
                thisBin = thisVel.fwdBin;
                thisData = vel_srvfwd;
            case 2 %ang
                thisBin = thisVel.angBin;
                thisData = vel_srvang;
            case 3 %sid
                thisBin = thisVel.sidBin;
                thisData = vel_srvsid;
        end
        % calculate mean and sem
        thisMean = mean(thisData,2,'omitnan');
        thisSEM= std(thisData,[],2,'omitnan')./sqrt(nFliesThresh);

        % plot mean + trials
        subplot(2,3,v); hold on
        plot(thisBin,thisData,'Color',settings.trialColor)
        plot(thisBin,thisMean,'Color',settings.velColor{v},'Linewidth',settings.lwAvg)
        axis tight; ylim(sr_limitT); xline(0)
        ylabel('Firing rate (spikes/s)'); xlabel(settings.velLabel{v})

        % plot mean +/- sem
        subplot(2,3,v+3); hold on
        plot(thisBin,thisMean,'Color',settings.velColor{v},'Linewidth',settings.lwAvg)
        er = patch([thisBin'; flipud(thisBin')],[(thisMean-thisSEM); flipud((thisMean+thisSEM))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
        er.FaceColor = settings.velColor{v};
        axis tight; ylim(sr_limitA); xline(0)
        ylabel('Firing rate (spikes/s)'); xlabel(settings.velLabel{v})

    end
    sgtitle([strrep(filebase,'_','/') ' Velocity (n = ' num2str(nFliesThresh) ')'])
    % save plot
    cd(folder.summary)
    plotname = 'pulse_sr_v_vel';
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

    %% analyze spike rate v directional velocity w/lag
    disp('Analyzing spike rate v directional velocity w/lag...')

    figure; set(gcf,'Position',[100 100 800 800])
    for v = 1:3
        switch v
            case 1 %fwd
                thisBin = thisVelL.fwdBin;
                thisData = velL_srvfwd;
                thisXlim = [0 10];
            case 2 %ang
                thisBin = thisVelL.angBin;
                thisData = velL_srvang;
                thisXlim = [-200 200];
            case 3 %sid
                thisBin = thisVelL.sidBin;
                thisData = velL_srvsid;
                thisXlim = [-2 2];
        end
        % calculate mean and sem
        thisMean = mean(thisData,2,'omitnan');
        thisSEM= std(thisData,[],2,'omitnan')./sqrt(nFliesThresh);

        % plot mean + trials
        subplot(2,3,v); hold on
        plot(thisBin,thisData,'Color',settings.trialColor)
        plot(thisBin,thisMean,'Color',settings.velColor{v},'Linewidth',settings.lwAvg)
        axis tight; ylim(sr_limitT); xlim(thisXlim); xline(0)
        ylabel('Firing rate (spikes/s)'); xlabel(settings.velLabel{v})

        % plot mean +/- sem
        subplot(2,3,v+3); hold on
        plot(thisBin,thisMean,'Color',settings.velColor{v},'Linewidth',settings.lwAvg)
        er = patch([thisBin'; flipud(thisBin')],[(thisMean-thisSEM); flipud((thisMean+thisSEM))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
        er.FaceColor = settings.velColor{v};
        axis tight; ylim(sr_limitA); xlim(thisXlim); xline(0)
        ylabel('Firing rate (spikes/s)'); xlabel(settings.velLabel{v})

    end
    sgtitle([strrep(filebase,'_','/') ' Velocity w/lag (n = ' num2str(nFliesThresh) ')'])
    % save plot
    cd(folder.summary)
    plotname = 'pulse_sr_v_vel_lag';
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

    % Save the slopes as a .mat file in the specified folder
    cd(folder.compare);
    savename = [filebase '_velocity_binned.mat'];
    save(savename, 'thisVelL','velL_srvfwd','velL_srvang','velL_srvsid');

    %% analyze CHANGE in spike rate v directional velocity w/lag
    disp('Analyzing change in spike rate v directional velocity w/lag...')
    sr_limitA2 = [-25 25];

    figure; set(gcf,'Position',[100 100 800 800])
    for v = 1:3
        switch v
            case 1 %fwd
                thisBin = thisVelL.fwdBin;
                thisData = velL_srvfwd-velL_srvfwd(1,:);
            case 2 %ang
                thisBin = thisVelL.angBin;
                thisData = velL_srvang-velL_srvang(thisBin==0,:);
            case 3 %sid
                thisBin = thisVelL.sidBin;
                thisData = velL_srvsid-velL_srvsid(thisBin==0,:);
        end
        % calculate mean and sem
        thisMean = mean(thisData,2,'omitnan');
        thisSEM= std(thisData,[],2,'omitnan')./sqrt(nFliesThresh);

        % plot mean + trials
        subplot(2,3,v); hold on
        plot(thisBin,thisData,'Color',settings.trialColor)
        plot(thisBin,thisMean,'Color',settings.velColor{v},'Linewidth',settings.lwAvg)
        axis tight; ylim(sr_limitT-mean(sr_limitT)); xline(0)
        ylabel('Firing rate change (spikes/s)'); xlabel(settings.velLabel{v})

        % plot mean +/- sem
        subplot(2,3,v+3); hold on
        plot(thisBin,thisMean,'Color',settings.velColor{v},'Linewidth',settings.lwAvg)
        er = patch([thisBin'; flipud(thisBin')],[(thisMean-thisSEM); flipud((thisMean+thisSEM))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
        er.FaceColor = settings.velColor{v};
        axis tight; ylim(sr_limitA2); xline(0)
        ylabel('Firing rate change (spikes/s)'); xlabel(settings.velLabel{v})

    end
    sgtitle([strrep(filebase,'_','/') ' Velocity w/lag (n = ' num2str(nFliesThresh) ')'])
    % save plot
    cd(folder.summary)
    plotname = 'pulse_changesr_v_vel_lag';
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

    %% Plot angular velocity tuning split by forward velocity
    disp('Plotting angular tuning split by forward velocity...')

    sr_limitA2 = [-25 25]; % Set Y-axis limits
    v = 2; % Angular velocity index in settings

    thisBin = thisVelL.angBin;
    thisData_slow = velL_srvang_slowfwd;
    thisData_fast = velL_srvang_fastfwd;

    % Compute mean and SEM for each
    thisMean_slow = mean(thisData_slow, 2, 'omitnan');
    thisMean_slow(isnan(thisMean_slow)) = 0;
    thisSEM_slow = std(thisData_slow, [], 2, 'omitnan') ./ sqrt(nFliesThresh);
    thisSEM_slow(isnan(thisSEM_slow)) = 0;

    thisMean_fast = mean(thisData_fast, 2, 'omitnan');
    thisMean_fast(isnan(thisMean_fast)) = 0;
    thisSEM_fast = std(thisData_fast, [], 2, 'omitnan') ./ sqrt(nFliesThresh);
    thisSEM_fast(isnan(thisSEM_fast)) = 0;

    % Plotting
    figure; set(gcf,'Position',[100 100 550 900])

    % SLOW forward tile
    subplot(2,2,1); hold on
    plot(thisBin, thisData_slow, 'Color', settings.trialColor)
    plot(thisBin, thisMean_slow, 'Color', settings.velColor{v}, 'LineWidth', settings.lwAvg)
    axis tight; ylim(sr_limitA); xline(0)
    title('Angular tuning: slow forward')
    xlabel('Angular velocity (deg/s)')
    ylabel('Firing rate (spikes/s)')

    subplot(2,2,3); hold on
    plot(thisBin, thisMean_slow, 'Color', settings.velColor{v}, 'LineWidth', settings.lwAvg)
    er = patch([thisBin'; flipud(thisBin')], ...
        [(thisMean_slow - thisSEM_slow); flipud(thisMean_slow + thisSEM_slow)], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    er.FaceColor = settings.velColor{2};
    axis tight; ylim(sr_limitA); xline(0)
    xlabel('Angular velocity (deg/s)')
    ylabel('Firing rate (spikes/s)')

    % FAST forward tile
    subplot(2,2,2); hold on
    plot(thisBin, thisData_fast, 'Color', settings.trialColor)
    plot(thisBin, thisMean_fast, 'Color', settings.velColor{v}, 'LineWidth', settings.lwAvg)
    axis tight; ylim(sr_limitA); xline(0)
    title('Angular tuning: fast forward')
    xlabel('Angular velocity (deg/s)')
    ylabel('Firing rate (spikes/s)')

    subplot(2,2,4); hold on
    plot(thisBin, thisMean_fast, 'Color', settings.velColor{v}, 'LineWidth', settings.lwAvg)
    er = patch([thisBin'; flipud(thisBin')], ...
        [(thisMean_fast - thisSEM_fast); flipud(thisMean_fast + thisSEM_fast)], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    er.FaceColor = settings.velColor{2};
    axis tight; ylim(sr_limitA); xline(0)
    xlabel('Angular velocity (deg/s)')
    ylabel('Firing rate (spikes/s)')

    % Save plot
    sgtitle([strrep(filebase,'_','/') ' Angular split by forward (n = ' num2str(nFliesThresh) ')'])
    cd(folder.summary)
    plotname = 'split_fwd_ang_tuning';
    saveas(gcf, [plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox, 'f');

    % Save vector format
    cd(folder.vector)
    set(gcf,'Renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox, 'f');


    %% Plot slopes for directional velocity data
    % Define colors from settings
    scatter_color = settings.trialColor;      % Color for individual points
    median_color = 'k';      % Red for median marker
    jitter_amount = 0.1;                      % Jitter for scatter points

    % Create figure and tiled layout
    figure; set(gcf, 'Position', [100 100 300 500])
    tiledlayout(1, 3, 'TileSpacing', 'compact')

    % Tile 1: Forward Velocity Slopes
    nexttile;
    hold on;
    scatter(ones(size(fwd_fits)) + jitter_amount * (rand(size(fwd_fits)) - 0.5), fwd_fits, ...
        '.', 'MarkerEdgeColor', scatter_color, 'MarkerFaceColor', scatter_color);
    median_fwd = median(fwd_fits, 'omitnan');
    plot(1, median_fwd, '_', 'MarkerSize', 15, 'Color', median_color, 'LineWidth', 2); % Median marker
    xticks(1);
    xticklabels({'Condition'});
    xlim([0 2]);
    yline(0);
    ylabel('Forward Velocity Slope');
    hold off;

    % Tile 2: Angular Velocity Slopes
    nexttile;
    hold on;
    scatter(ones(size(ang_fits)) + jitter_amount * (rand(size(ang_fits)) - 0.5), ang_fits, ...
        '.', 'MarkerEdgeColor', scatter_color, 'MarkerFaceColor', scatter_color);
    median_ang = median(ang_fits, 'omitnan');
    plot(1, median_ang, '_', 'MarkerSize', 15, 'Color', median_color, 'LineWidth', 2); % Median marker
    xticks(1);
    xticklabels({'Condition'});
    xlim([0 2]);
    yline(0);
    ylabel('Angular Velocity Slope');
    hold off;

    % Tile 3: Sideways Velocity Slopes
    nexttile;
    hold on;
    scatter(ones(size(sid_fits)) + jitter_amount * (rand(size(sid_fits)) - 0.5), sid_fits, ...
        '.', 'MarkerEdgeColor', scatter_color, 'MarkerFaceColor', scatter_color);
    median_sid = median(sid_fits, 'omitnan');
    plot(1, median_sid, '_', 'MarkerSize', 15, 'Color', median_color, 'LineWidth', 2); % Median marker
    xticks(1);
    xticklabels({'Condition'});
    xlim([0 2]);
    yline(0);
    ylabel('Sideways Velocity Slope');
    hold off;

    % Global title
    sgtitle('Velocity Slope Distributions for Condition');

    % Save plot
    cd(folder.summary);
    plotname = 'velocity_slope_distribution';
    saveas(gcf, [plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox, 'f');
    % Save vectorized plot
    cd(folder.vector);
    set(gcf, 'renderer', 'Painters');
    saveas(gcf, [plotname '.svg']);
    copyfile([plotname '.svg'], folder.dropbox, 'f');

    % Define the structure to hold the slopes
    fits.fwd = fwd_fits;
    fits.fwdr2 = fwd_r2;
    fits.ang = ang_fits;
    fits.sid = sid_fits;

    % Save the slopes as a .mat file in the specified folder
    cd(folder.compare);
    savename = [filebase '_velocity_slopes.mat'];
    save(savename, 'fits','flylist');

    %% analyze spike rate v directional acceleration
    disp('Analyzing spike rate v directional acceleration...')

    figure; set(gcf,'Position',[100 100 800 800])
    for v = 1:3
        switch v
            case 1 %fwd
                thisBin = thisAcc.fwdBin;
                thisData = acc_srvfwd;
            case 2 %ang
                thisBin = thisAcc.angBin;
                thisData = acc_srvang;
            case 3 %sid
                thisBin = thisAcc.sidBin;
                thisData = acc_srvsid;
        end
        % calculate mean and sem
        thisMean = mean(thisData,2,'omitnan');
        thisSEM= std(thisData,[],2,'omitnan')./sqrt(nFliesThresh);

        % plot mean + trials
        subplot(2,3,v); hold on
        plot(thisBin,thisData,'Color',settings.trialColor)
        plot(thisBin,thisMean,'Color',settings.velColor{v},'Linewidth',settings.lwAvg)
        axis tight; ylim(sr_limitT); xline(0)
        ylabel('Firing rate (spikes/s)'); xlabel(settings.accLabel{v})

        % plot mean +/- sem
        subplot(2,3,v+3); hold on
        plot(thisBin,thisMean,'Color',settings.velColor{v},'Linewidth',settings.lwAvg)
        er = patch([thisBin'; flipud(thisBin')],[(thisMean-thisSEM); flipud((thisMean+thisSEM))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
        er.FaceColor = settings.velColor{v};
        axis tight; ylim(sr_limitA); xline(0)
        ylabel('Firing rate (spikes/s)'); xlabel(settings.accLabel{v})

    end
    sgtitle([strrep(filebase,'_','/') ' Acceleration (n = ' num2str(nFliesThresh) ')'])
    % save plot
    cd(folder.summary)
    plotname = 'pulse_sr_v_acc';
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

    %% Determine lag and cross-correlation for each trial type for firing rate
    disp('Performing cross-correlation for firing vs pursuit behavior across conditions...')
    cd(folder.summary)

    % Initialize figure and tiled layout
    figure; set(gcf, 'Position', [100 100 600 800])  % Adjust width to fit 3 columns
    tiledlayout(3, 3, 'TileSpacing', 'compact')
    xc_lim = [-400 400];
    r_range = [-0.1 1];

    % For each directional velocity
    for v = 1:3
        switch v
            case 1
                r_val = r_val_fwd;
                lag_pk = pk_lag_fwd;
                r_pk = pk_rval_fwd;
                thisV = 'Forward Rval';
            case 2
                r_val = r_val_ang;
                lag_pk = pk_lag_ang;
                r_pk = pk_rval_ang;
                thisV = 'Angular Rval';
            case 3
                r_val = r_val_sid;
                lag_pk = pk_lag_sid;
                r_pk = pk_rval_sid;
                thisV = 'Sideways Rval';
        end

        % Calculate median r-values across lags
        r_median = median(r_val, 2, 'omitnan');

        % Calculate median for peak lag and peak r values
        peak_lag_median = median(lag_pk, 'omitnan');
        peak_r_median = median(r_pk, 'omitnan');

        % Plot cross-correlation r-values over time lags
        nexttile; hold on
        plot(lag_t, r_val, ':', 'LineWidth', 0.5, 'Color', settings.trialColor)
        plot(lag_t, r_median, 'LineWidth', 1.5, 'Color', settings.velColor{v})
        xline(0); yline(0); ylim(r_range); ylabel(thisV); xlabel('Lag (msec)')

        % Plot peak correlation (r_pk) as dash marker for median
        nexttile; hold on
        plot(1, r_pk, '.', 'Color', settings.trialColor)
        plot(1, peak_r_median, '_', 'MarkerSize', 15, 'Color', settings.velColor{v}, 'LineWidth', 2)
        ylabel('Peak Correlation (Rval)')
        ylim(r_range)
        yline(0);
        set(gca, 'XTick', [])

        % Plot peak lags (lag_pk) as dash marker for median
        nexttile; hold on
        plot(1, lag_pk, '.', 'Color', settings.trialColor)
        plot(1, peak_lag_median, '_', 'MarkerSize', 15, 'Color', settings.velColor{v}, 'LineWidth', 2)
        ylabel('Peak Lag (msec)')
        ylim(xc_lim)
        yline(0)
        set(gca, 'XTick', [])
        % Add label above median point for peak lag
        text(1, 250, num2str(peak_lag_median), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', 8, 'Rotation', 90);
    end

    % Set global title and save plot
    sgtitle([strrep(filebase, '_', ' ') ' (n = ' num2str(nFliesThresh) ')'])
    cd(folder.summary)
    plotname = 'pulse_xcorr';
    saveas(gcf, [plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox, 'f');

    % Save vectorized version
    cd(folder.vector)
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox, 'f');

    % Store peak data for rotational (angular) velocity after the loop
    rotational_peaks.lag_pk = pk_lag_ang;
    rotational_peaks.r_pk = pk_rval_ang;
    forward_peaks.lag_pk = pk_lag_fwd;
    forward_peaks.r_pk = pk_rval_fwd;
    % Save rotational peak data as a .mat file in folder.compare
    cd(folder.compare)
    save([ filebase '_xcorr.mat'], 'rotational_peaks','forward_peaks','storeNames');




    %% analyze directional velocity v spike rate w/lag
    disp('Analyzing directional velocity v spike rate w/lag...')

    thisSpikeBin = binSR_withLag.spikeRateBin;
    sr_lim = [0 80];

    figure; set(gcf,'Position',[100 100 800 800])
    for v = 1:3
        switch v
            case 1 % fwd
                thisData = sr_velfwd;
                thisYlim = [0 10];
            case 2 % ang
                thisData = sr_velang;
                thisYlim = [-200 200];
            case 3 % sid
                thisData = sr_velsid;
                thisYlim = [-2 2];
        end

        % count number of animals with data per bin
        validN = sum(~isnan(thisData), 2);

        % only use bins with at least 3 animals
        validIdx = validN >= 5;
        xVals = thisSpikeBin(validIdx);
        yData = thisData(validIdx, :);
        validN_used = validN(validIdx);  % match yData size

        % compute mean and SEM across valid data
        thisMean = mean(yData, 2, 'omitnan');
        thisSEM  = std(yData, [], 2, 'omitnan') ./ sqrt(validN_used);

        % ensure xVals is column vector for patch compatibility
        xVals = xVals(:);

        % plot individual trials and mean
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
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

    % Save the slopes as a .mat file in the specified folder
    cd(folder.compare);
    savename = [filebase '_velocity_binnedbyfr.mat'];
    save(savename, 'thisSpikeBin','sr_velfwd','sr_velang','sr_velsid');

    %% Analyze motion pulse vs spikerate across forward velocity groups
    if motionCheck
        disp('Analyzing motion pulse vs spikerate by forward velocity group...')
        cd(folder.summary)
        sr_limitD = [-100 100];

        fwd_groups = {'slow fwd';'fast fwd'};
        minFwd = 4;

        % initialize
        mean_pulse_srR = [];
        mean_pulse_srL = [];
        sem_pulse_srR = [];
        sem_pulse_srL = [];
        peak_srR_store = [];

        thisName = 'byfwd';

        % for each pulse speed
        for p = 1:nPulses
            % initialize
            thisSpeed = pulseSpeeds{p}(1:2);
            pulse_srRightward = [];
            pulse_srLeftward = [];
            c = 1; %counter

            % for each fly
            for nt = 1:nFlies
                % select data
                thisPanelPs = allPanelPs{nt};
                thisForward = allForward{nt};
                thisSpikert = allSpikeRt{nt};
                % determine relationship between pulse and spikerate
                [~, thisMean] = pulse_v_output_by_fwd(thisPanelPs,thisForward,thisSpikert,int_time,p,pulseSpeeds,nSweep,minFwd);
                % if first run, store motion pulse positions
                if c == 1
                    % sweep positions
                    pulse_posR(:,p,:) = thisMean.panelpsR;
                    pulse_posL(:,p,:) = thisMean.panelpsL;
                    pulseDur = size(pulse_posR,1);
                    time_pulse = int_time(1:size(pulse_posR,1))*1000;
                    % sweep centers for each position
                    sweepPosR = round(pulse_posR(round(pulseDur/2),1,:));
                    sweepPosL = round(pulse_posL(round(pulseDur/2),1,:));
                    % sweep indices
                    sweepIdx = find(~isnan(pulse_posR(:,p,1)));
                    sweepIdx2End = sweepIdx(1):length(pulse_posR(:,p,1));
                end
                % store spikerate averages for each fly
                pulse_srRL_Rightward(:,c,:) = thisMean.stationary.varOutR - flip(thisMean.stationary.varOutL,3);
                pulse_srRL_slowfwd(:,c,:) = thisMean.slow.varOutR - flip(thisMean.slow.varOutL,3);
                pulse_srRL_fastfwd(:,c,:) = thisMean.fast.varOutR - flip(thisMean.fast.varOutL,3);
                % for each sweep, store peak spikerate for each fly
                for s = find(sweepPosR<0)
                    avg_srRL_Rightward(c,s) = mean(pulse_srRL_Rightward(sweepIdx,c,s),'omitnan');
                    peak_srRL_slowfwd(c,s) = mean(pulse_srRL_slowfwd(sweepIdx,c,s),'omitnan');
                    peak_srRL_fastfwd(c,s) = mean(pulse_srRL_fastfwd(sweepIdx,c,s),'omitnan');
                end
                for s = find(sweepPosR>0)
                    avg_srRL_Rightward(c,s) = mean(pulse_srRL_Rightward(sweepIdx,c,s),'omitnan');
                    peak_srRL_slowfwd(c,s) = mean(pulse_srRL_slowfwd(sweepIdx,c,s),'omitnan');
                    peak_srRL_fastfwd(c,s) = mean(pulse_srRL_fastfwd(sweepIdx,c,s),'omitnan');
                end
                c=c+1;
            end
            % Find flies with complete data in both forward velocity groups
            validFlies = all(~isnan(peak_srRL_slowfwd), 2) & all(~isnan(peak_srRL_fastfwd), 2);

            % Reduce all peak data to valid flies only
            avg_srRL_Rightward = avg_srRL_Rightward(validFlies, :);
            peak_srRL_slowfwd = peak_srRL_slowfwd(validFlies, :);
            peak_srRL_fastfwd = peak_srRL_fastfwd(validFlies, :);
            pulse_srRL_Rightward = pulse_srRL_Rightward(:, validFlies, :);
            pulse_srRL_slowfwd = pulse_srRL_slowfwd(:, validFlies, :);
            pulse_srRL_fastfwd = pulse_srRL_fastfwd(:, validFlies, :);


            % Count new nFlies
            nFlies_adj_station = sum(validFlies);
            nFlies_adj_fwd = nFlies_adj_station;

            % calculate spikerate means and sem
            mean_pulse_srRL_slowfwd(:,p,:) = mean(pulse_srRL_slowfwd,2,'omitnan');
            sem_pulse_srRL_slowfwd(:,p,:) = std(pulse_srRL_slowfwd,0,2,'omitnan')./sqrt(nFlies_adj_fwd);
            mean_pulse_srRL_fastfwd(:,p,:) = mean(pulse_srRL_fastfwd,2,'omitnan');
            sem_pulse_srRL_fastfwd(:,p,:) = std(pulse_srRL_fastfwd,0,2,'omitnan')./sqrt(nFlies_adj_fwd);
            % calculate spikerate peaks mean and sem
            mean_peak_srRL_slowfwd = mean(peak_srRL_slowfwd,1,'omitnan');
            sem_peak_srRL_slowfwd = std(peak_srRL_slowfwd,0,1,'omitnan')./sqrt(nFlies_adj_fwd);
            mean_peak_srRL_fastfwd = mean(peak_srRL_fastfwd,1,'omitnan');
            sem_peak_srRL_fastfwd = std(peak_srRL_fastfwd,0,1,'omitnan')./sqrt(nFlies_adj_fwd);

            % plot averages for each speed
            % initialize
            figure; set(gcf,'Position',[100 100 plotWidth 800])
            tiledlayout(4,2,"TileSpacing","compact")
            % initialize remaining variables
            time_pulse = int_time(1:size(pulse_posR,1))*1000;
            time_limit = [0 max(time_pulse)];

            % plot target position
            for g = 1:2
                nexttile; hold on
                for s = 1:nSweep
                    plot(time_pulse,pulse_posR(:,p,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
                end
                title(fwd_groups{g}); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
                xlim(time_limit); ylim(ps_limit); yline(0)
            end

            % plot R-L spikerate
            for g = 1:2
                switch g
                    case 1
                        meanRL = mean_pulse_srRL;
                        semRL = sem_pulse_srRL;
                    case 2
                        meanRL = mean_pulse_srRL_slowfwd;
                        semRL = sem_pulse_srRL_slowfwd;
                    case 3
                        meanRL = mean_pulse_srRL_fastfwd;
                        semRL = sem_pulse_srRL_fastfwd;
                end

                nexttile([3 1]); hold on
                for s = 1:nSweep
                    sp = patch([time_pulse'; flipud(time_pulse')],[meanRL(:,p,s)-semRL(:,p,s); flipud(meanRL(:,p,s)+semRL(:,p,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
                    sp.FaceColor = color_rightward(s,:);
                    plot(time_pulse,meanRL(:,p,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
                end
                ylabel('Expected Firing Rate Difference (spikes/s)'); xlabel('Time (msec)')
                xlim(time_limit); ylim(sr_limitD); yline(0)
            end

            sgtitle([strrep([filebase ' ' pulseSpeeds{p} ' ' thisName],'_','/') ' (n = ' num2str(nFlies_adj_station) ')'])
            % save plot
            cd(folder.summary)
            plotname = strjoin({'pulse_v_fr', thisSpeed, 'dps', thisName},'_');
            saveas(gcf,[plotname '.png']);
            copyfile([plotname '.png'], folder.dropbox,'f');
            % save vectorized plot
            cd(folder.vector)
            set(gcf,'renderer','Painters')
            saveas(gcf, [plotname '.svg'])
            copyfile([plotname '.svg'], folder.dropbox,'f');


            % initialize
            figure; set(gcf,'Position',[100 100 400 800])

            % plot peak fr for R-L sweeps
            hold on
            sweepPosR = reshape(sweepPosR, [], 1);
            %errorbar(sweepPosR, mean_peak_srRL_station, sem_peak_srRL_station, 'k', 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
            errorbar(sweepPosR, mean_peak_srRL_slowfwd, sem_peak_srRL_slowfwd, 'b', 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
            errorbar(sweepPosR, mean_peak_srRL_fastfwd, sem_peak_srRL_fastfwd, 'r', 'LineStyle', '-', 'Marker', 'none', 'CapSize', 0);
            xlabel('Sweep Pos (deg)'); ylabel('Peak FR')
            legend(fwd_groups)
            xlim([0 150]); ylim([-10 50]); xticks([0:30:150]); yline(0); xline(0)

            sgtitle(strrep([filebase ' ' pulseSpeeds{p} ' ' thisName],'_','/'))

            % Combine data
            [nAnimals, nSweep] = size(avg_srRL_Rightward);

            % Identify positive sweep positions
            keepSweepIdx = find(sweepPosR > 0);
            nKeep = numel(keepSweepIdx);

            % Extract just those sweep positions from each matrix
            peak_slow_keep = peak_srRL_slowfwd(:, keepSweepIdx);
            peak_fastfwd_keep = peak_srRL_fastfwd(:, keepSweepIdx);
            sweepPos_keep = sweepPosR(keepSweepIdx);

            % Reshape into long-form vectors
            peak_all = [peak_slow_keep(:); peak_fastfwd_keep(:)];
            sweepPos = repmat(sweepPos_keep, 1, nAnimals * 2)';
            sweepPos = sweepPos(:);
            animalID = repmat((1:nAnimals)', nKeep * 2, 1);
            behavior = [ones(nAnimals * nKeep, 1); 2 * ones(nAnimals * nKeep, 1)];

            T = table(peak_all, sweepPos, categorical(behavior), categorical(animalID), ...
                'VariableNames', {'PeakFR', 'Sweep', 'Behavior', 'Animal'});
            if ~isempty(T)
                lme = fitlme(T, 'PeakFR ~ Sweep * Behavior + (1|Animal)');

                % Display ANOVA table for fixed effects
                anovaTbl = anova(lme);  % <-- capture the output
                disp(anovaTbl);         % optional, to still display it

                % Extract p-values from ANOVA table
                p_sweep     = anovaTbl.pValue(strcmp(anovaTbl.Term, 'Sweep'));
                p_behavior  = anovaTbl.pValue(strcmp(anovaTbl.Term, 'Behavior'));
                p_interact  = anovaTbl.pValue(strcmp(anovaTbl.Term, 'Sweep:Behavior'));

                % Format text including interaction
                txt = sprintf('p(Sweep) = %.3g\np(Fwd) = %.3g\np(Intx) = %.3g', ...
                    p_sweep, p_behavior, p_interact);

                % Add to plot (left-aligned)
                xTxt = 0;    % left edge of x-axis
                yTxt = 40;   % top of y-axis

                text(xTxt, yTxt, txt, ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'top', ...
                    'FontSize', 8);

                % save plot
                cd(folder.summary)
                plotname = strjoin({'pulse_v_fr', thisSpeed, 'dps', 'norm',thisName},'_');
                saveas(gcf,[plotname '.png']);
                copyfile([plotname '.png'], folder.dropbox,'f');
                % save vectorized plot
                cd(folder.vector)
                set(gcf,'renderer','Painters')
                saveas(gcf, [plotname '.svg'])
                copyfile([plotname '.svg'], folder.dropbox,'f');
            end

        end
    end
end
%% end
disp('ALL ANALYSES COMPLETE.')
end

