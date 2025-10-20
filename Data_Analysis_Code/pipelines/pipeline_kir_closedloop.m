% pipeline_kir_closedloop
%
% This pipeline function processes and analyzes data from all flies
% in a given closed-loop experiment. It pulls all relevant processed files,
% performs the necessary analyses, and generates plots accordingly.
%
% INPUT:
%   exptFolder - String representing the path to the overarching experiment
%                folder containing processed data for all flies.
%
% Created: 08/23/2024 MC Adapted from the open loop pipeline
% Updated: 08/18/2025 MC made friendly to 1 gain conditions (aka 025)
%
function pipeline_kir_closedloop(exptFolder)
%% Initialize
disp('STARTING ANALYSES FOR POOLED KIR PURSUIT...')
close all  % Close all open figures

% Load processing settings
settings = processSettings();

% Replace spaces in the folder name with underscores
filebase = strrep(exptFolder, ' ', '_');

% Generate folder structure for saving outputs
folder = generateFolders(exptFolder);

% Change to the intermediate folder
cd(folder.int)

% Find all processed data files for the experiment
dataFiles = dir('*int.mat');

% Number of flies based on data files
nFlies = length(dataFiles);

% Number of pursuit gains from settings
if contains(exptFolder,'AOTU025 KIR')
    nGain = 1;
else
    nGain = length(settings.pursuitGain);
end

% Minimum fixation time to include a fly (in seconds)
minFixationTime = settings.minFixationTime;

%% Load in and pool pursuit data
disp('Loading in and analyzing pursuit datasets...')
nKIR = 0; nWT = 0; nNA = 0;  % Initialize counters for each group
oKIR = 0; oWT = 0; oNA = 0;  % Initialize counters for omitted trials
walkThresh = 1; %mm/s

for e = 1:nFlies
    disp(['Processing fly ' num2str(e) '/' num2str(nFlies) '...'])
    % Load current trial
    cd(folder.int)
    thisTrial = dataFiles(e).name;
    load(thisTrial)
    etime = int_time; etime(end+1) = 60;  % Append trial end time

    % Find fixation periods
    thisFixation = fixationFinder(int_panelps,int_forward,int_time,0);
    fix_panelps = thisFixation.panelps_run;
    fix_panelvel = int_panelvel;
    fix_panelvel(~thisFixation.idx_run) = nan;  % Exclude non-fixation data
    fix_forward = thisFixation.forward_run;
    fix_angular = int_angular;
    fix_angular(~thisFixation.idx_run) = nan;
    fix_sideway = int_sideway;
    fix_sideway(~thisFixation.idx_run) = nan;
    % Find walking periods
    walk_panelps = int_panelps;
    walkIdx = [];
    for g = 1:nGain
        walkIdx(:,:,g) = schmittTrigger(int_forward(:,:,g), walkThresh, 0.1);
    end
    walk_panelps(~walkIdx) = nan;

    % Calculate run parameters per trial condition
    thisFixationT = reshape(sum(etime(sum(~isnan(fix_panelps))+1)),1,nGain);  % Fixation time (s)
    thisFixationFwdV = reshape(mean(fix_forward,[1 2],'omitnan'),1,nGain);  % Avg forward speed (mm/s)
    thisFixationAngV = reshape(mean(abs(fix_angular),[1 2],'omitnan'),1,nGain);  % Avg angular speed (deg/s)

    % Calculate run parameters across all conditions
    thisFixationT(end+1) = mean(sum(etime(sum(~isnan(fix_panelps))+1)));  % Overall fixation time
    thisFixationFwdV(end+1) = mean(fix_forward,'all','omitnan');  % Overall forward speed
    thisFixationAngV(end+1) = mean(abs(fix_angular),'all','omitnan');  % Overall angular speed

    % If fly meets fixation criteria
    if all(thisFixationT(1:end-1) > minFixationTime)
        % Analyze pursuit performance
        [opt_lag, r_val] = behavior_v_panelps_lag(fix_panelps, fix_angular, int_time, settings);
        sp_out = setpoint_performance(fix_panelps, int_jumptrg, fix_forward, int_time, 0);
        freq_hd_out = setpoint_freq(fix_panelps, int_time, 0);
        freq_ang_out = setpoint_freq(fix_angular, int_time, 0);
        [smallJumps, largeJumps] = setpoint_jumps(fix_panelps, int_jumptrg, int_time); %binned
        [objAtJump, corrTimes] = setpoint_jumps2(fix_panelps, int_jumptrg, int_time); %separated
        [posvang, posvangRL, posBins] = setpoint_errorvturn(fix_panelps, fix_angular, int_time, settings, 1, 0);
        [posvfwd, ~, ~] = setpoint_errorvturn(fix_panelps, fix_forward, int_time, settings, 0, 0);
        [~, posvangRL_byfwd, ~] = setpoint_errorvturn_byFwd(fix_panelps, fix_angular, fix_forward, int_time, settings, 1, 0);
        [hist_slow, hist_fast] = compare_panelpos_byfwd(fix_panelps, fix_forward);
        [velvang, velvangFrontFOVRL, velBins] = setpoint_velvturn(fix_panelps,fix_panelvel,fix_angular,int_time,settings,1,0);
        [~, posvaccelRL, ~] = setpoint_errorvaccel(fix_panelps, fix_angular, int_time, settings, 1, 0);
        [dc_t, angbin_prior_means , angbin_cross_means, bins] = setpoint_dirchange(fix_panelps, fix_angular, int_jumptrg, int_time, 0);
        [ang_maxbins, bins2] = setpoint_max_to_max(fix_panelps, fix_angular, int_jumptrg, int_time, 0);
        [medianFixationTimes, medianNonFixationTimes, medianFixationObjPos, medianNonFixationObjPos] = setpoint_crossings(int_panelps, int_angular, int_forward, thisFixation.idx_run, int_time, settings);
        [bc_output, binsbig] = setpoint_lrgchange(fix_panelps, fix_angular, int_jumptrg, int_time, 0);

        fixidx_95x = thisFixation.idx_run(:,:,1);
        [fixation_percentage_P1, fixation_percentage_noP1] = compare_fixation_time(fixidx_95x, thisTrial, folder);
        [hist_noP1, hist_P1] = compare_panelpos_histogram(int_panelps(:,:,1), thisTrial, folder);

        % Analyze velocity distributions
        for c = 1:nGain
            [fwdHist, angHist, sidHist] = velocity_histogram(fix_forward(:,:,c), fix_angular(:,:,c), fix_sideway(:,:,c), 1);
            thisFwd(:,c) = fwdHist(:,2);
            thisAng(:,c) = angHist(:,2);
            thisSid(:,c) = sidHist(:,2);
            [posHist, velHist] = panel_histogram(fix_panelps(:,:,c), fix_panelvel(:,:,c), 1);
            thisPanelpos(:,c) = posHist(:,2);
            thisPanelvel(:,c) = velHist(:,2);
        end

        % Pool trial data by condition
        if contains(thisTrial,'KIR')
            nKIR = nKIR + 1;
            kirRunTime(nKIR,:) = thisFixationT;
            kirRunSpeed(nKIR,:) = thisFixationFwdV;
            kirTurnSpeed(nKIR,:) = thisFixationAngV;
            kirOptLag(nKIR,:) = opt_lag';
            kirRVal(:,:,nKIR) = r_val;
            kirSPprob(nKIR,:) = sp_out.prob;
            kirSPiae(nKIR,:) = sp_out.iae;
            kirSPise(nKIR,:) = sp_out.ise;
            kirSPcstd(nKIR,:) = sp_out.cstd;
            kirSPcvar(nKIR,:) = sp_out.cvar;
            kirFFT(:,:,nKIR) = freq_hd_out.fft;
            kirPSD(:,:,nKIR) = freq_hd_out.psd;
            kirFFTang(:,:,nKIR) = freq_ang_out.fft;
            kirPSDang(:,:,nKIR) = freq_ang_out.psd;
            kirEVT(:,:,nKIR) = posvang;
            kirEVTRL(:,:,nKIR) = posvangRL;
            kirEVTRL_lowfwd(:,:,nKIR) = posvangRL_byfwd(:,:,1);
            kirEVTRL_highfwd(:,:,nKIR) = posvangRL_byfwd(:,:,end);
            kirEVAccRL(:,:,nKIR) = posvaccelRL;
            kirVVT(:,:,nKIR) = velvang;
            kirVVTfront(:,:,nKIR) = velvangFrontFOVRL;
            kirEVFwd(:,:,nKIR) = posvfwd;
            kirBigChange(:,:,nKIR) = bc_output;
            kirDCAngPriorBins(:,:,nKIR) = angbin_prior_means;
            kirDCAngCrossBins(:,:,nKIR) = angbin_cross_means;
            kirMCAngBins(:,:,nKIR) = ang_maxbins;
            kirFwdHist(:,:,nKIR) = thisFwd;
            kirAngHist(:,:,nKIR) = thisAng;
            kirSidHist(:,:,nKIR) = thisSid;
            kirPosHist(:,:,nKIR) = thisPanelpos;
            kirVelHist(:,:,nKIR) = thisPanelvel;
            kirFixInt(nKIR,:) = medianFixationTimes;
            kirNotInt(nKIR,:) = medianNonFixationTimes;
            kirFixPos(nKIR,:) = medianFixationObjPos;
            kirNotPos(nKIR,:) = medianNonFixationObjPos;
            kirJumpSmall(nKIR,:) = smallJumps;
            kirJumpLarge(nKIR,:) = largeJumps;
            kirJumpAll_pos(nKIR,:) = objAtJump;
            kirJumpAll_time(nKIR,:) = corrTimes;

            kirFixRatio(nKIR,1) = fixation_percentage_noP1;
            kirFixRatio(nKIR,2) = fixation_percentage_P1;

            kirHDslow(:,nKIR) = hist_slow(:,2);
            kirHDfast(:,nKIR) = hist_fast(:,2);
            kirHDnoP1(:,nKIR) = hist_noP1(:,2);
            kirHDP1(:,nKIR) = hist_P1(:,2);

        elseif contains(thisTrial,'WT') % WT flies
            nWT = nWT + 1;
            wtRunTime(nWT,:) = thisFixationT;
            wtRunSpeed(nWT,:) = thisFixationFwdV;
            wtTurnSpeed(nWT,:) = thisFixationAngV;
            wtOptLag(nWT,:) = opt_lag';
            wtRVal(:,:,nWT) = r_val;
            wtSPprob(nWT,:) = sp_out.prob;
            wtSPiae(nWT,:) = sp_out.iae;
            wtSPise(nWT,:) = sp_out.ise;
            wtSPcstd(nWT,:) = sp_out.cstd;
            wtSPcvar(nWT,:) = sp_out.cvar;
            wtFFT(:,:,nWT) = freq_hd_out.fft;
            wtPSD(:,:,nWT) = freq_hd_out.psd;
            wtFFTang(:,:,nWT) = freq_ang_out.fft;
            wtPSDang(:,:,nWT) = freq_ang_out.psd;
            wtEVT(:,:,nWT) = posvang;
            wtEVTRL(:,:,nWT) = posvangRL;
            wtEVTRL_lowfwd(:,:,nWT) = posvangRL_byfwd(:,:,1);
            wtEVTRL_highfwd(:,:,nWT) = posvangRL_byfwd(:,:,end);
            wtEVAccRL(:,:,nWT) = posvaccelRL;
            wtVVT(:,:,nWT) = velvang;
            wtVVTfront(:,:,nWT) = velvangFrontFOVRL;
            wtEVFwd(:,:,nWT) = posvfwd;
            wtBigChange(:,:,nWT) = bc_output;
            wtDCAngPriorBins(:,:,nWT) = angbin_prior_means;
            wtDCAngCrossBins(:,:,nWT) = angbin_cross_means;
            wtMCAngBins(:,:,nWT) = ang_maxbins;
            wtFwdHist(:,:,nWT) = thisFwd;
            wtAngHist(:,:,nWT) = thisAng;
            wtSidHist(:,:,nWT) = thisSid;
            wtPosHist(:,:,nWT) = thisPanelpos;
            wtVelHist(:,:,nWT) = thisPanelvel;
            wtFixInt(nWT,:) = medianFixationTimes;
            wtNotInt(nWT,:) = medianNonFixationTimes;
            wtFixPos(nWT,:) = medianFixationObjPos;
            wtNotPos(nWT,:) = medianNonFixationObjPos;
            wtJumpSmall(nWT,:) = smallJumps;
            wtJumpLarge(nWT,:) = largeJumps;
            wtJumpAll_pos(nWT,:) = objAtJump;
            wtJumpAll_time(nWT,:) = corrTimes;

            wtFixRatio(nWT,1) = fixation_percentage_noP1;
            wtFixRatio(nWT,2) = fixation_percentage_P1;
            wtHDslow(:,nWT) = hist_slow(:,2);
            wtHDfast(:,nWT) = hist_fast(:,2);
            wtHDnoP1(:,nWT) = hist_noP1(:,2);
            wtHDP1(:,nWT) = hist_P1(:,2);

        elseif contains(thisTrial,'NA') % NA flies
            nNA = nNA + 1;
            naRunTime(nNA,:) = thisFixationT;
            naRunSpeed(nNA,:) = thisFixationFwdV;
            naTurnSpeed(nNA,:) = thisFixationAngV;
            naOptLag(nNA,:) = opt_lag';
            naRVal(:,:,nNA) = r_val;
            naSPprob(nNA,:) = sp_out.prob;
            naSPiae(nNA,:) = sp_out.iae;
            naSPise(nNA,:) = sp_out.ise;
            naSPcstd(nNA,:) = sp_out.cstd;
            naSPcvar(nNA,:) = sp_out.cvar;
            naFFT(:,:,nNA) = freq_hd_out.fft;
            naPSD(:,:,nNA) = freq_hd_out.psd;
            naFFTang(:,:,nNA) = freq_ang_out.fft;
            naPSDang(:,:,nNA) = freq_ang_out.psd;
            naEVT(:,:,nNA) = posvang;
            naEVTRL(:,:,nNA) = posvangRL;
            naEVTRL_lowfwd(:,:,nNA) = posvangRL_byfwd(:,:,1);
            naEVTRL_highfwd(:,:,nNA) = posvangRL_byfwd(:,:,end);
            naEVAccRL(:,:,nNA) = posvaccelRL;
            naVVT(:,:,nNA) = velvang;
            naVVTfront(:,:,nNA) = velvangFrontFOVRL;
            naEVFwd(:,:,nNA) = posvfwd;
            naBigChange(:,:,nNA) = bc_output;
            naDCAngPriorBins(:,:,nNA) = angbin_prior_means;
            naDCAngCrossBins(:,:,nNA) = angbin_cross_means;
            naMCAngBins(:,:,nNA) = ang_maxbins;
            naFwdHist(:,:,nNA) = thisFwd;
            naAngHist(:,:,nNA) = thisAng;
            naSidHist(:,:,nNA) = thisSid;
            naPosHist(:,:,nNA) = thisPanelpos;
            naVelHist(:,:,nNA) = thisPanelvel;
            naFixInt(nNA,:) = medianFixationTimes;
            naNotInt(nNA,:) = medianNonFixationTimes;
            naFixPos(nNA,:) = medianFixationObjPos;
            naNotPos(nNA,:) = medianNonFixationObjPos;
            naJumpSmall(nNA,:) = smallJumps;
            naJumpLarge(nNA,:) = largeJumps;
            naJumpAll_pos(nNA,:) = objAtJump;
            naJumpAll_time(nNA,:) = corrTimes;

            naFixRatio(nNA,1) = fixation_percentage_noP1;
            naFixRatio(nNA,2) = fixation_percentage_P1;
            naHDslow(:,nNA) = hist_slow(:,2);
            naHDfast(:,nNA) = hist_fast(:,2);
            naHDnoP1(:,nNA) = hist_noP1(:,2);
            naHDP1(:,nNA) = hist_P1(:,2);
        end
    else
        disp('Omitted.')  % If trial is omitted

        % Pool omitted trial data
        if contains(thisTrial,'KIR') % KIR flies omitted
            oKIR = oKIR + 1;
            okirRunTime(oKIR,:) = thisFixationT;
            okirRunSpeed(oKIR,:) = thisFixationFwdV;
            okirTurnSpeed(oKIR,:) = thisFixationAngV;
        elseif contains(thisTrial,'WT') % WT flies omitted
            oWT = oWT + 1;
            owtRunTime(oWT,:) = thisFixationT;
            owtRunSpeed(oWT,:) = thisFixationFwdV;
            owtTurnSpeed(oWT,:) = thisFixationAngV;
        elseif contains(thisTrial,'NA') % NA flies omitted
            oNA = oNA + 1;
            onaRunTime(oNA,:) = thisFixationT;
            onaRunSpeed(oNA,:) = thisFixationFwdV;
            onaTurnSpeed(oNA,:) = thisFixationAngV;
        end
    end
end

% Output frequency data
freq_fft = freq_hd_out.f_fft;
freq_psd = freq_hd_out.f_psd;
binsbig = binsbig.ang_vel_center;

disp('Complete.')


%% Basic parameters
disp('Comparing base parameters...')
close all

% Call the general ANOVA function to analyze the contribution of genotype and gain to general performance
[p_runTime, ~] = run_genotype_anova_repeated(kirRunTime(:, 1:nGain), wtRunTime(:, 1:nGain), naRunTime(:, 1:nGain), 'RunTime', folder);
[p_runSpeed, ~] = run_genotype_anova_repeated(kirRunSpeed(:, 1:nGain), wtRunSpeed(:, 1:nGain), naRunSpeed(:, 1:nGain), 'RunSpeed', folder);
[p_turnSpeed, ~] = run_genotype_anova_repeated(kirTurnSpeed(:, 1:nGain), wtTurnSpeed(:, 1:nGain), naRunSpeed(:, 1:nGain), 'TurnSpeed', folder);
anova_pvals = {p_runTime, p_runSpeed, p_turnSpeed};

% Call the general ANOVA function to analyze the contribution of genotype to general performance
[p_runTime2, ~] = run_genotype_anova1(kirRunTime(:, end), wtRunTime(:, end), naRunTime(:, end), 'RunTime2', folder);
[p_runSpeed2, ~] = run_genotype_anova1(kirRunSpeed(:, end), wtRunSpeed(:, end), naRunSpeed(:, end), 'RunSpeed2', folder);
[p_turnSpeed2, ~] = run_genotype_anova1(kirTurnSpeed(:, end), wtTurnSpeed(:, end), naRunSpeed(:, end), 'TurnSpeed2', folder);
anova_pvals2 = {p_runTime2, p_runSpeed2, p_turnSpeed2};

% Initialize figure
figure; set(gcf, 'Position', [100 100 600 800])
tiledlayout(3, 4, 'TileSpacing', 'compact')

for p = 1:3
    % Fetch parameters based on index
    switch p
        case 1
            kirData = kirRunTime; wtData = wtRunTime; naData = naRunTime;
            nameData = 'Fixation Time (s)'; yrange = [0 400];
            pval_text2 = ['p(gene) = ' num2str(anova_pvals2{p}(1))];
        case 2
            kirData = kirRunSpeed; wtData = wtRunSpeed; naData = naRunSpeed;
            nameData = 'Fwd Speed (mm/s)'; yrange = [0 20];
            pval_text2 = ['p(gene) = ' num2str(anova_pvals2{p}(1))];
        case 3
            kirData = kirTurnSpeed; wtData = wtTurnSpeed; naData = naTurnSpeed;
            nameData = 'Turn Speed (deg/s)'; yrange = [0 150];
            pval_text2 = ['p(gene) = ' num2str(anova_pvals2{p}(1))];
    end

    % Calculate median for each condition
    kirMedian = median(kirData, 1, 'omitnan');
    wtMedian = median(wtData, 1, 'omitnan');
    naMedian = median(naData, 1, 'omitnan');

    % Plot steering gain separately
    nexttile([1 3]); hold on
    plot(settings.pursuitGain(1:nGain) - 2, kirData(:, 1:nGain), '.', 'Color', settings.trialColor)
    plot(settings.pursuitGain(1:nGain) - 2, kirMedian(1:nGain), '_', 'Color', settings.geneColor{1}, 'MarkerSize', 8)

    plot(settings.pursuitGain(1:nGain), wtData(:, 1:nGain), '.', 'Color', settings.trialColor)
    plot(settings.pursuitGain(1:nGain), wtMedian(1:nGain), '_', 'Color', settings.geneColor{2}, 'MarkerSize', 8)

    plot(settings.pursuitGain(1:nGain) + 2, naData(:, 1:nGain), '.', 'Color', settings.trialColor)
    plot(settings.pursuitGain(1:nGain) + 2, naMedian(1:nGain), '_', 'Color', settings.geneColor{3}, 'MarkerSize', 8)

    axis padded; ylabel(nameData); xticks(settings.pursuitGain(1:nGain)); ylim(yrange); xlabel('Steering Gain (k)')
    if p == 1, yline(minFixationTime, ':'); end

    % Plot steering gain together
    nexttile; hold on

    % Define jitter magnitude
    jitterMagnitude = 0.2;

    % Plot Kir data with jitter
    jitteredX_KIR = 1 + jitterMagnitude * (rand(1, nKIR) - 0.5);
    plot(jitteredX_KIR, kirData(:, end), '.', 'Color', settings.trialColor);
    plot(1, kirMedian(end), '_', 'Color', settings.geneColor{1}, 'MarkerSize', 8);

    % Plot WT data with jitter
    jitteredX_WT = 2 + jitterMagnitude * (rand(1, nWT) - 0.5);
    plot(jitteredX_WT, wtData(:, end), '.', 'Color', settings.trialColor);
    plot(2, wtMedian(end), '_', 'Color', settings.geneColor{2}, 'MarkerSize', 8);

    % Plot NA data with jitter
    jitteredX_NA = 3 + jitterMagnitude * (rand(1, nNA) - 0.5);
    plot(jitteredX_NA, naData(:, end), '.', 'Color', settings.trialColor);
    plot(3, naMedian(end), '_', 'Color', settings.geneColor{3}, 'MarkerSize', 8);


    axis padded; ylabel(nameData); xlim([0 4]); xticks(1:3); xticklabels(settings.geneLabel); ylim(yrange); xlabel('Fly Average')
    text(0.95, 0.95, pval_text2, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    if p == 1, yline(minFixationTime, ':'); end
end

sgtitle([ 'Included Summary (n = ' num2str(nKIR) 'kir, ' num2str(nWT) 'wt, ' num2str(nNA) 'na)'])

% Save plot
cd(folder.summary)
plotname = 'basics_summary_included';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.')


%% Repeat for omitted flies
% % Call the general ANOVA function to analyze the contribution of genotype and gain to omitted data
% [p_oRunTime, ~] = run_genotype_anova_repeated(okirRunTime, owtRunTime, onaRunTime, 'Omitted_RunTime', folder);
% [p_oRunSpeed, ~] = run_genotype_anova_repeated(okirRunSpeed, owtRunSpeed, onaRunSpeed, 'Omitted_RunSpeed', folder);
% [p_oTurnSpeed, ~] = run_genotype_anova_repeated(okirTurnSpeed, owtTurnSpeed, onaTurnSpeed, 'Omitted_TurnSpeed', folder);
%
% % Store p-values from ANOVA for each omitted parameter
% anova_pvals_omitted = {p_oRunTime, p_oRunSpeed, p_oTurnSpeed};
%
% % Initialize omitted figure
% figure; set(gcf, 'Position', [100 100 600 800])
% tiledlayout(3, 4, 'TileSpacing', 'compact')
%
% for p = 1:3
%     % Fetch parameters for omitted data
%     switch p
%         case 1
%             kirData = okirRunTime; wtData = owtRunTime; naData = onaRunTime;
%             nameData = 'Fixation Time (s)'; yrange = [0 400];
%             pval_text = ['p(gene) = ' num2str(anova_pvals_omitted{p}(1)) ', p(k) = ' num2str(anova_pvals_omitted{p}(2))];
%         case 2
%             kirData = okirRunSpeed; wtData = owtRunSpeed; naData = onaRunSpeed;
%             nameData = 'Fwd Speed (mm/s)'; yrange = [0 30];
%             pval_text = ['p(gene) = ' num2str(anova_pvals_omitted{p}(1)) ', p(k) = ' num2str(anova_pvals_omitted{p}(2))];
%         case 3
%             kirData = okirTurnSpeed; wtData = owtTurnSpeed; naData = onaTurnSpeed;
%             nameData = 'Turn Speed (deg/s)'; yrange = [0 150];
%             pval_text = ['p(gene) = ' num2str(anova_pvals_omitted{p}(1)) ', p(k) = ' num2str(anova_pvals_omitted{p}(2))];
%     end
%
%     % Calculate median for each condition
%     kirMedian = median(kirData, 1, 'omitnan');
%     wtMedian = median(wtData, 1, 'omitnan');
%     naMedian = median(naData, 1, 'omitnan');
%
%     % Plot steering gain separately
%     nexttile([1 3]); hold on
%     plot(settings.pursuitGain - 2, kirData(:, 1:nGain), '.', 'Color', settings.trialColor)
%     plot(settings.pursuitGain - 2, kirMedian(1:nGain), '_', 'Color', settings.geneColor{1}, 'MarkerSize', 8)
%
%     plot(settings.pursuitGain, wtData(:, 1:nGain), '.', 'Color', settings.trialColor)
%     plot(settings.pursuitGain, wtMedian(1:nGain), '_', 'Color', settings.geneColor{2}, 'MarkerSize', 8)
%
%     plot(settings.pursuitGain + 2, naData(:, 1:nGain), '.', 'Color', settings.trialColor)
%     plot(settings.pursuitGain + 2, naMedian(1:nGain), '_', 'Color', settings.geneColor{3}, 'MarkerSize', 8)
%
%     axis padded; ylabel(nameData); xticks(settings.pursuitGain); ylim(yrange); xlabel('Steering Gain (k)')
%     if p == 1, yline(minFixationTime, ':'); end
%
%     % Add p-value annotation to the plot (top right corner)
%     text(0.95, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
%
%     % Plot steering gain together
%     nexttile; hold on
%
%     % Define jitter magnitude
%     jitterMagnitude = 0.1;
%
%     % Plot Kir data with jitter
%     jitteredX_KIR = 1 + jitterMagnitude * (rand(1, oKIR) - 0.5);
%     plot(jitteredX_KIR, kirData(:, end), '.', 'Color', settings.trialColor);
%     plot(1, kirMedian(end), '_', 'Color', settings.geneColor{1}, 'MarkerSize', 8);
%
%     % Plot WT data with jitter
%     jitteredX_WT = 2 + jitterMagnitude * (rand(1, oWT) - 0.5);
%     plot(jitteredX_WT, wtData(:, end), '.', 'Color', settings.trialColor);
%     plot(2, wtMedian(end), '_', 'Color', settings.geneColor{2}, 'MarkerSize', 8);
%
%     % Plot NA data with jitter
%     jitteredX_NA = 3 + jitterMagnitude * (rand(1, oNA) - 0.5);
%     plot(jitteredX_NA, naData(:, end), '.', 'Color', settings.trialColor);
%     plot(3, naMedian(end), '_', 'Color', settings.geneColor{3}, 'MarkerSize', 8);
%
%
%     axis padded; ylabel(nameData); xlim([0 4]); xticks(1:3); xticklabels(settings.geneLabel); ylim(yrange); xlabel('Fly Average')
%     if p == 1, yline(minFixationTime, ':'); end
% end
%
% sgtitle([ 'Omitted Summary (n = ' num2str(oKIR) 'kir, ' num2str(oWT) 'wt, ' num2str(oNA) 'na)'])
%
% % Save omitted plot
% cd(folder.summary)
% plotname = 'basics_summary_omitted';
% saveas(gcf, [plotname '.png']);
% copyfile([plotname '.png'], folder.dropbox, 'f');
%
% % Save omitted vectorized plot
% cd(folder.vector)
% set(gcf, 'renderer', 'Painters')
% saveas(gcf, [plotname '.svg'])
% copyfile([plotname '.svg'], folder.dropbox, 'f');
% disp('Complete.')

%% Compare fixation ratios
% Initialize figure and tiled layout
figure;
set(gcf, 'Position', [100 100 700 600]);
tiledlayout(1, 4, 'TileSpacing', 'compact');

% Define colors and marker properties
scatterColor = [0.5 0.5 0.5];  % Grey color for scatter points
noP1_medianColor = 'k';         % Black dash for no P1 median
P1_medianColor = 'r';           % Red dash for P1 median

% Plot each genotype in separate tiles
genotypes = {'KIR', 'WT', 'NA'};
fixRatios = {kirFixRatio, wtFixRatio, naFixRatio};

for i = 1:3
    nexttile;
    hold on;

    % Get fixation data for the current genotype
    fixRatio = fixRatios{i};

    % X-coordinates for no P1 and P1 points
    x_coords = [1, 2];

    % Use plot to connect points between columns 1 and 2 for each row
    for p = 1:size(fixRatio, 1)
        plot(x_coords, fixRatio(p, :), '-', 'Marker', '.', 'MarkerEdgeColor', scatterColor, 'LineWidth', 0.5);
        hold on;
    end
    % Plot median markers with dash style
    noP1_median = median(fixRatio(:, 1), 'omitnan');
    P1_median = median(fixRatio(:, 2), 'omitnan');
    scatter(1, noP1_median, '_', 'MarkerEdgeColor', noP1_medianColor, 'MarkerFaceColor', noP1_medianColor, 'LineWidth', 1);
    scatter(2, P1_median, '_', 'MarkerEdgeColor', P1_medianColor, 'MarkerFaceColor', P1_medianColor, 'LineWidth', 1);

    % Set axis properties
    ylim([0 100])
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'No P1', 'With P1'});
    ylabel('Fixation Percentage (%)');
    title(genotypes{i});

    hold off;
end

% Combine data for all genotypes
dependentVar_kir = kirFixRatio;
dependentVar_wt = wtFixRatio;
dependentVar_na = naFixRatio;

% Define the name of the dependent variable
dependent_var_name = 'Fixation_Percentage';

% Run linear mixed-effects model (Genotype × Condition, with Animal as random effect)
% Assign genotype labels
nKIR = size(dependentVar_kir,1);
nWT  = size(dependentVar_wt,1);
nNA  = size(dependentVar_na,1);

% Stack data
data_all = [
    reshape(dependentVar_kir.', [], 1);
    reshape(dependentVar_wt.',  [], 1);
    reshape(dependentVar_na.',  [], 1)
    ];

% Define categorical variables
condition = repmat({'NoP1'; 'P1'}, [nKIR + nWT + nNA, 1]);
condition = categorical(condition);

genotype = [
    repmat({'KIR'}, nKIR * 2, 1);
    repmat({'WT'},  nWT * 2, 1);
    repmat({'NA'},  nNA * 2, 1)
    ];
genotype = categorical(genotype);

% Assign animal ID (repeated per condition)
animalID = [
    repelem((1:nKIR)', 2);
    repelem((1:nWT)' + nKIR, 2);
    repelem((1:nNA)' + nKIR + nWT, 2)
    ];
animalID = categorical(animalID);

% Create table
T = table(data_all, genotype, condition, animalID, ...
    'VariableNames', {'Fixation', 'Genotype', 'Condition', 'AnimalID'});

% Run LME
lme = fitlme(T, 'Fixation ~ Genotype * Condition + (1|AnimalID)');

% Output ANOVA table
anovaTbl = anova(lme);
disp(anovaTbl);

% Extract and display p-values
p_genotype = anovaTbl.pValue(strcmp(anovaTbl.Term, 'Genotype'));
p_condition = anovaTbl.pValue(strcmp(anovaTbl.Term, 'Condition'));
p_interaction = anovaTbl.pValue(strcmp(anovaTbl.Term, 'Genotype:Condition'));

disp(['p(Genotype): ', num2str(p_genotype)]);
disp(['p(Condition): ', num2str(p_condition)]);
disp(['p(Interaction): ', num2str(p_interaction)]);

% Plot combined data in the fourth tile
allFixRatios = cat(1, kirFixRatio, wtFixRatio, naFixRatio);
nexttile;
hold on;

% Use plot to connect points between columns 1 and 2 for each row
for x = 1:size(allFixRatios, 1)
    plot([1, 2], allFixRatios(x, :), '-', 'Marker', '.', 'MarkerEdgeColor', scatterColor, 'LineWidth', 0.5);
end

% Plot median markers with dash style
noP1_median_all = median(allFixRatios(:, 1), 'omitnan');
P1_median_all = median(allFixRatios(:, 2), 'omitnan');
scatter(1, noP1_median_all, '_', 'MarkerEdgeColor', noP1_medianColor, 'MarkerFaceColor', noP1_medianColor, 'LineWidth', 1);
scatter(2, P1_median_all, '_', 'MarkerEdgeColor', P1_medianColor, 'MarkerFaceColor', P1_medianColor, 'LineWidth', 1);

% Set axis properties
ylim([0 100])
xlim([0 3]);
xticks([1 2]);
xticklabels({'No P1', 'With P1'});
ylabel('Fixation Percentage (%)');
title('All Genotypes');
hold off;

% Adjust overall figure properties
sgtitle('Fixation Percentage');

% Save the figure
cd(folder.summary)
plotname = 'basics_fixationpercent';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters');
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.');

[~, p_kir] = ttest(kirFixRatio(:,1), kirFixRatio(:,2));
[~, p_wt]  = ttest(wtFixRatio(:,1),  wtFixRatio(:,2));
[~, p_na]  = ttest(naFixRatio(:,1),  naFixRatio(:,2));
pvals = [p_kir, p_wt, p_na];
adj_p = min(pvals * 3, 1);  % Bonferroni correction


%% HD with and without P1
% Combine across genotypes
all_noP1 = [kirHDnoP1, wtHDnoP1, naHDnoP1];  % [bins x flies]
all_P1   = [kirHDP1, wtHDP1, naHDP1];

% Extract bin centers (assumes same binning across genotypes)
bin_centers = hist_noP1(:,1);

% Calculate mean and SEM
mean_noP1 = mean(all_noP1, 2, 'omitnan');
sem_noP1  = std(all_noP1, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(all_noP1), 2));
mean_P1   = mean(all_P1, 2, 'omitnan');
sem_P1    = std(all_P1, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(all_P1), 2));

% Define interpolation resolution
smooth_resolution = 0.5;  % finer spacing in degrees
bin_interp = min(bin_centers):smooth_resolution:max(bin_centers);

% Interpolate mean and SEM with smoothing splines
mean_noP1_smooth = interp1(bin_centers, mean_noP1, bin_interp, 'pchip');
sem_noP1_smooth  = interp1(bin_centers, sem_noP1, bin_interp, 'pchip');

mean_P1_smooth   = interp1(bin_centers, mean_P1, bin_interp, 'pchip');
sem_P1_smooth    = interp1(bin_centers, sem_P1, bin_interp, 'pchip');

% Create new figure with two rows
figure;
tiledlayout(2,2, 'TileSpacing', 'compact');

% ==== Tile 1: Line plot ====
nexttile(1); hold on

% No P1 (black)
patch([bin_interp, fliplr(bin_interp)], ...
    [mean_noP1_smooth - sem_noP1_smooth, fliplr(mean_noP1_smooth + sem_noP1_smooth)], ...
    'k', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
plot(bin_interp, mean_noP1_smooth, 'k-', 'LineWidth', 1.5);

% P1 (red)
patch([bin_interp, fliplr(bin_interp)], ...
    [mean_P1_smooth - sem_P1_smooth, fliplr(mean_P1_smooth + sem_P1_smooth)], ...
    'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
plot(bin_interp, mean_P1_smooth, 'r-', 'LineWidth', 1.5);

xlabel('Panel Position (°)');
ylabel('Probability');
legend({'No P1 ± SEM','No P1','P1 ± SEM','P1'}, 'Location', 'best');
title('Panel Position Histogram Across Genotypes');
xlim([-180 180]); yline(0); box off;
xticks(-150:30:150)

% --- Make bin EDGES from your bin CENTERS (uniform bins assumed) ---
bin_centers = bin_centers(:);                       % ensure column
bin_width   = median(diff(bin_centers),'omitnan');  % degrees
bin_edges_deg = [bin_centers - bin_width/2; bin_centers(end) + bin_width/2];
bin_edges_rad = deg2rad(bin_edges_deg);             % radians

% ==== Tile 2: Polar histogram - No P1 ====
rmax = 0.125;
nexttile(3);
polarhistogram('BinEdges', bin_edges_rad', ...
    'BinCounts', mean_noP1(:), ...
    'FaceColor', 'k', 'FaceAlpha', 0.6, ...
    'Normalization', 'probability');
rlim([0 rmax]);
title('Polar Histogram - No P1');

% ==== Tile 3: Polar histogram - P1 ====
nexttile(4);
polarhistogram('BinEdges', bin_edges_rad', ...
    'BinCounts', mean_P1(:), ...
    'FaceColor', 'r', 'FaceAlpha', 0.6, ...
    'Normalization', 'probability');
rlim([0 rmax]);
title('Polar Histogram - P1');

% Save the figure
cd(folder.summary)
plotname = 'basics_hdwithwithoutp1';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters');
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.');


%% Optimal visual-motor lag
disp('Comparing optimal lag time...')

% Call the general ANOVA function to analyze the contribution of genotype and gain to lag times
if exptFolder == 'AOTU019 KIR'
    [p_lag, ~] = run_genotype_anova_repeated(kirOptLag, wtOptLag, naOptLag, 'Visual-Motor Lag', folder);
    pval_text_lag = ['p(gene) = ' num2str(p_lag(1)) ', p(k) = ' num2str(p_lag(2))];
else
    [p_lag, ~] = run_genotype_anova1(kirOptLag, wtOptLag, naOptLag, 'Visual-Motor Lag', folder);
    pval_text_lag = ['p(gene) = ' num2str(p_lag(1))];
end

% Initialize figure layout for optimal lag comparison
figure; set(gcf, 'Position', [100 100 600 300])
tiledlayout(1, 5, 'TileSpacing', 'compact')

% Plot estimated lag for each genotype with median dash markers
nexttile([1, 5])
for g = 1:3
    switch g
        case 1
            this_lag = kirOptLag; thisN = nKIR; x = settings.pursuitGain(1:nGain) - 2;
        case 2
            this_lag = wtOptLag; thisN = nWT; x = settings.pursuitGain(1:nGain);
        case 3
            this_lag = naOptLag; thisN = nNA; x = settings.pursuitGain(1:nGain) + 2;
    end

    % Calculate median for each gain condition
    median_lag = median(this_lag, 'omitnan');

    % Plot individual lag data points and dash markers for the median
    hold on
    plot(x, this_lag, '.', 'Color', settings.trialColor)      % Individual data points
    plot(x, median_lag, '_', 'Color', settings.geneColor{g}, 'MarkerSize', 8)  % Median as dash marker
    axis padded; ylim([0 250]); xticks(settings.pursuitGain(1:nGain))
    ylabel('Est. Visual-Motor Turn Lag (ms)');
end
xlabel('Steering Gain (X)')

% Add p-value annotation to the plot (top right corner)
text(0.95, 0.95, pval_text_lag, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

% Save plot as PNG and SVG formats
cd(folder.summary)
plotname = 'basics_visualmotorlag';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.')

%% Setpoint performance for low gain only
close all

% --- Extract first column per genotype
kir = kirSPprob(:,1);
wt  = wtSPprob(:,1);
na  = naSPprob(:,1);

% --- X positions and labels
xpos   = [1 2 3];
labels = {'Kir','WT','Na'};
colKir = settings.geneColor{1};
colWT  = settings.geneColor{2};
colNa  = settings.geneColor{3};

% --- Plot
figure('Color','w'); set(gcf, 'Position', [100 100 200 400]); hold on;
yl = [0 1];
jit = 0.06;              % horizontal jitter
ms  = 5;                % marker size for dots

% helper to scatter with jitter in black, ignoring NaNs
plot_jitter = @(x,y) arrayfun(@(yy) ...
    plot(x + (rand*2-1)*jit, yy, '.', 'Color', [0 0 0], 'MarkerSize', ms), ...
    y(~isnan(y)));

% Plot per genotype
plot_jitter(xpos(1), kir);
plot_jitter(xpos(2), wt);
plot_jitter(xpos(3), na);

% Medians (horizontal dashed lines in genotype color)
medKir = median(kir,'omitnan');
medWT  = median(wt, 'omitnan');
medNa  = median(na, 'omitnan');

line([xpos(1)-0.25 xpos(1)+0.25], [medKir medKir], 'Color', colKir, 'LineStyle','-', 'LineWidth',2);
line([xpos(2)-0.25 xpos(2)+0.25], [medWT  medWT ], 'Color', colWT,  'LineStyle','-', 'LineWidth',2);
line([xpos(3)-0.25 xpos(3)+0.25], [medNa  medNa ], 'Color', colNa,  'LineStyle','-', 'LineWidth',2);

% Axes & cosmetics
xlim([0.5 3.5]);
ylim(yl);
xticks(xpos); xticklabels(labels);
ylabel('Setpoint Probability');
box on;

% --- Mixed-effects model: SPprob ~ Genotype + (1 | Animal)
% Build a long table with Genotype and Animal IDs
vals       = [kir; wt; na];
geno       = [repelem({'Kir'}, numel(kir))'; repelem({'WT'}, numel(wt))'; repelem({'Na'}, numel(na))'];
animal_id  = arrayfun(@(k) sprintf('Kir_%02d',k), (1:numel(kir))', 'uni',0);
animal_id  = [animal_id; arrayfun(@(k) sprintf('WT_%02d',k), (1:numel(wt))', 'uni',0)];
animal_id  = [animal_id; arrayfun(@(k) sprintf('Na_%02d',k), (1:numel(na))', 'uni',0)];

T = table(vals, categorical(geno), categorical(animal_id), ...
          'VariableNames', {'SPprob','Genotype','Animal'});

% Drop rows with NaN
T = T(~isnan(T.SPprob), :);

% Fit LME (random intercept per animal, fixed Genotype)
lme = fitlme(T, 'SPprob ~ 1 + Genotype + (1|Animal)');

% Display results
disp('--- Mixed-Effects Model (SPprob ~ Genotype + (1|Animal)) ---');
disp(lme);
disp('--- Fixed Effects ANOVA (effect of Genotype) ---');
result = anova(lme,'DFMethod','Satterthwaite');
disp(result);

text(0.5, 0.9, ['p = ' num2str(result.pValue(2))], 'Units', 'normalized', 'FontSize', 7);

cd(folder.summary)
plotname = 'basics_setpoint_stats95';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.')


%% Setpoint performance parameters
disp('Comparing setpoint performance...')

% Initialize figure layout for setpoint performance comparison
figure; set(gcf, 'Position', [100 100 1200 400])
tiledlayout(1, 5, 'TileSpacing', 'compact')

% Initialize p-values storage
anova_pvals_setpoint = cell(1, 5);

% Plot setpoint performance metrics
for p = 1:5
    % Select dataset and labels for each metric
    switch p
        case 1, kirData = kirSPprob; wtData = wtSPprob; naData = naSPprob; nameData = 'Setpoint Probability'; yrange = [0 0.7];
        case 2, kirData = kirSPiae; wtData = wtSPiae; naData = naSPiae; nameData = 'Integral of Absolute Error'; yrange = [0 3e3];
        case 3, kirData = kirSPise; wtData = wtSPise; naData = naSPise; nameData = 'Integral of Squared Error'; yrange = [0 2e5];
        case 4, kirData = kirSPcstd; wtData = wtSPcstd; naData = naSPcstd; nameData = 'Circular STD'; yrange = [0 1];
        case 5, kirData = kirSPcvar; wtData = wtSPcvar; naData = naSPcvar; nameData = '1-Circular Variance'; yrange = [0 1];
    end

    % Call the general ANOVA function to analyze the contribution of genotype and gain
    if contains(exptFolder,'AOTU019 KIR')
        [p_val, ~] = run_genotype_anova_repeated(kirData, wtData, naData, nameData, folder);
        anova_pvals_setpoint{p} = p_val;  % Store p-values
        pval_text = {['p(gene) = ' num2str(anova_pvals_setpoint{p}(1))] ; ['p(k) = ' num2str(anova_pvals_setpoint{p}(2))] ; ['p(g*k) = ' num2str(anova_pvals_setpoint{p}(3))]};
    else
        [p_val, ~] = run_genotype_anova1(kirData, wtData, naData, nameData, folder);
        anova_pvals_setpoint{p} = p_val;  % Store p-values
        pval_text = {['p(gene) = ' num2str(anova_pvals_setpoint{p}(1))]};
    end

    % Calculate median and SEM for each condition
    kirMedian = median(kirData, 1, 'omitnan'); kirSEM = std(kirData, 0, 1, 'omitnan') ./ sqrt(nKIR);
    wtMedian = median(wtData, 1, 'omitnan'); wtSEM = std(wtData, 0, 1, 'omitnan') ./ sqrt(nWT);
    naMedian = median(naData, 1, 'omitnan'); naSEM = std(naData, 0, 1, 'omitnan') ./ sqrt(nNA);

    % Plot data with error bars for KIR, WT, and NA without caps
    nexttile; hold on
    errorbar(settings.pursuitGain(1:nGain)-2, kirMedian, kirSEM, '-', 'Color', settings.geneColor{1}, 'CapSize', 0)
    errorbar(settings.pursuitGain(1:nGain), wtMedian, wtSEM, '-', 'Color', settings.geneColor{2}, 'CapSize', 0)
    errorbar(settings.pursuitGain(1:nGain)+2, naMedian, naSEM, '-', 'Color', settings.geneColor{3}, 'CapSize', 0)
    axis padded; ylim(yrange); xticks(settings.pursuitGain); ylabel(nameData); xlabel('k')

    % Add p-value annotation to the plot (top right corner)
    if p == 5
        text(0.5, 0.5, pval_text, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    else
        text(0.95, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    end
end

legend(settings.geneLabel, 'Location', 'southeast')

% Save setpoint performance plot as PNG and SVG formats
cd(folder.summary)
plotname = 'basics_setpoint_stats';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.')


%% Frequency of HD position oscillations
disp('Comparing frequency of setpoint oscillations for HD...')

% Initialize figure for FFT comparison
figure; set(gcf,'Position',[100 100 1800 800])
tiledlayout(5,nGain,'TileSpacing','compact')
fftLim = [0 3]; yrange = [0 1]; yrange2 = [0 3];

% Plot FFT trials and means for each genotype
for g = 1:3
    switch g
        case 1, thisFFT = kirFFT; thisN = nKIR;
        case 2, thisFFT = wtFFT; thisN = nWT;
        case 3, thisFFT = naFFT; thisN = nNA;
    end
    % Plot each condition's FFT
    for c = 1:nGain
        conFFT = reshape(thisFFT(:,c,:),[],thisN);  % Fetch data
        meanFFT = mean(conFFT,2,'omitnan');  % Calculate mean

        % Plot FFT data and mean
        nexttile; hold on
        plot(freq_fft,conFFT,'Color',settings.trialColor,'LineWidth',settings.lwTri)
        plot(freq_fft,meanFFT,'Color',settings.geneColor{g},'LineWidth',settings.lwTri)
        xlim(fftLim); ylim(yrange2); xlabel('Freq (Hz)')
        if g==1, title([num2str(settings.pursuitGain(c)) 'X']); end
        if c==1, ylabel('Heading Amp (deg)'); end
    end
end

% Plot FFT means together for each condition
for c = 1:nGain
    % Fetch data and calculate mean and SEM
    thisKIR = reshape(kirFFT(:,c,:),[],nKIR);
    thisWT = reshape(wtFFT(:,c,:),[],nWT);
    thisNA = reshape(naFFT(:,c,:),[],nNA);
    meanKIR = mean(thisKIR,2,'omitnan');
    meanWT = mean(thisWT,2,'omitnan');
    meanNA = mean(thisNA,2,'omitnan');
    semKIR = std(thisKIR,0,2,'omitnan')./sqrt(nKIR);
    semWT = std(thisWT,0,2,'omitnan')./sqrt(nWT);
    semNA = std(thisNA,0,2,'omitnan')./sqrt(nNA);

    % Plot mean FFT with SEM
    nexttile([2 1]); hold on
    plot(freq_fft,meanKIR,'Color',settings.geneColor{1},'LineWidth',settings.lwTri)
    sem1 = patch([freq_fft'; flipud(freq_fft')],[meanKIR-semKIR; flipud(meanKIR+semKIR)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sem1.FaceColor = settings.geneColor{1};
    plot(freq_fft,meanWT,'Color',settings.geneColor{2},'LineWidth',settings.lwTri)
    sem2 = patch([freq_fft'; flipud(freq_fft')],[meanWT-semWT; flipud(meanWT+semWT)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sem2.FaceColor = settings.geneColor{2};
    plot(freq_fft,meanNA,'Color',settings.geneColor{3},'LineWidth',settings.lwTri)
    sem3 = patch([freq_fft'; flipud(freq_fft')],[meanNA-semNA; flipud(meanNA+semNA)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sem3.FaceColor = settings.geneColor{3};
    xlim(fftLim); ylim(yrange); xlabel('Freq (Hz)');
    title([num2str(settings.pursuitGain(c)) 'X'])
    if c==1, ylabel('Heading Amplitude (deg)'); end
end
sgtitle('FFT')

% Save FFT plot as PNG and SVG
cd(folder.summary)
plotname = 'frequency_HD_FFT';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

% Initialize figure for PSD comparison
figure; set(gcf,'Position',[100 100 1800 800])
tiledlayout(5,nGain,'TileSpacing','compact')
psdLim = [0 3]; yrange = [0 40];

% Plot PSD trials and means for each genotype
for g = 1:3
    switch g
        case 1, thisPSD = kirPSD; thisN = nKIR;
        case 2, thisPSD = wtPSD; thisN = nWT;
        case 3, thisPSD = naPSD; thisN = nNA;
    end
    % Plot each condition's PSD
    for c = 1:nGain
        conPSD = reshape(thisPSD(:,c,:),[],thisN);  % Fetch data
        meanPSD = mean(conPSD,2,'omitnan');  % Calculate mean

        % Plot PSD data and mean
        nexttile; hold on
        plot(freq_psd,10*log10(conPSD),'Color',settings.trialColor,'LineWidth',settings.lwTri)
        plot(freq_psd,10*log10(meanPSD),'Color',settings.geneColor{g},'LineWidth',settings.lwTri)
        xlim(psdLim); ylim(yrange); xlabel('Freq (Hz)')
        if g==1, title([num2str(settings.pursuitGain(c)) 'X']); end
        if c==1, ylabel('Heading Amp (Deg^2/Hz)'); end
    end
end

% Plot PSD means together for each condition
for c = 1:nGain
    % Fetch data and calculate mean and SEM
    thisKIR = 10*log10(reshape(kirPSD(:,c,:),[],nKIR));
    thisWT = 10*log10(reshape(wtPSD(:,c,:),[],nWT));
    thisNA = 10*log10(reshape(naPSD(:,c,:),[],nNA));
    meanKIR = mean(thisKIR,2,'omitnan');
    meanWT = mean(thisWT,2,'omitnan');
    meanNA = mean(thisNA,2,'omitnan');
    semKIR = std(thisKIR,0,2,'omitnan')./sqrt(nKIR);
    semWT = std(thisWT,0,2,'omitnan')./sqrt(nWT);
    semNA = std(thisNA,0,2,'omitnan')./sqrt(nNA);

    % Plot mean PSD with SEM
    nexttile([2 1]); hold on
    plot(freq_psd,meanKIR,'Color',settings.geneColor{1},'LineWidth',settings.lwTri)
    sem1 = patch([freq_psd; flipud(freq_psd)],[meanKIR-semKIR; flipud(meanKIR+semKIR)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sem1.FaceColor = settings.geneColor{1};
    plot(freq_psd,meanWT,'Color',settings.geneColor{2},'LineWidth',settings.lwTri)
    sem2 = patch([freq_psd; flipud(freq_psd)],[meanWT-semWT; flipud(meanWT+semWT)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sem2.FaceColor = settings.geneColor{2};
    plot(freq_psd,meanNA,'Color',settings.geneColor{3},'LineWidth',settings.lwTri)
    sem3 = patch([freq_psd; flipud(freq_psd)],[meanNA-semNA; flipud(meanNA+semNA)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sem3.FaceColor = settings.geneColor{3};
    xlim(psdLim); ylim(yrange); xlabel('Freq (Hz)');
    title([num2str(settings.pursuitGain(c)) 'X'])
    if c==1, ylabel('Heading Amplitude (Deg^2/Hz)'); end
end
sgtitle('Welchs PSD')

% Save PSD plot as PNG and SVG
cd(folder.summary)
plotname = 'frequency_HD_PSD';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

disp('Complete.')


%% Frequency of angular velocity oscillations
disp('Comparing frequency of setpoint oscillations for angular velocity...')

% Initialize figure for FFT comparison
figure; set(gcf,'Position',[100 100 1800 800])
tiledlayout(5,nGain,'TileSpacing','compact')
fftLim = [0 3]; yrange = [0 4]; yrange2 = [0 10];

% Plot FFT trials and means for each genotype
for g = 1:3
    switch g
        case 1, thisFFT = kirFFTang; thisN = nKIR;
        case 2, thisFFT = wtFFTang; thisN = nWT;
        case 3, thisFFT = naFFTang; thisN = nNA;
    end
    % Plot each condition's FFT
    for c = 1:nGain
        conFFT = reshape(thisFFT(:,c,:),[],thisN);  % Fetch data
        meanFFT = mean(conFFT,2,'omitnan');  % Calculate mean

        % Plot FFT data and mean
        nexttile; hold on
        plot(freq_fft,conFFT,'Color',settings.trialColor,'LineWidth',settings.lwTri)
        plot(freq_fft,meanFFT,'Color',settings.geneColor{g},'LineWidth',settings.lwTri)
        xlim(fftLim); ylim(yrange2); xlabel('Freq (Hz)')
        if g==1, title([num2str(settings.pursuitGain(c)) 'X']); end
        if c==1, ylabel('Amp (deg/s)'); end
    end
end

% Plot FFT means together for each condition
for c = 1:nGain
    % Fetch data and calculate mean and SEM
    thisKIR = reshape(kirFFTang(:,c,:),[],nKIR);
    thisWT = reshape(wtFFTang(:,c,:),[],nWT);
    thisNA = reshape(naFFTang(:,c,:),[],nNA);
    meanKIR = mean(thisKIR,2,'omitnan');
    meanWT = mean(thisWT,2,'omitnan');
    meanNA = mean(thisNA,2,'omitnan');
    semKIR = std(thisKIR,0,2,'omitnan')./sqrt(nKIR);
    semWT = std(thisWT,0,2,'omitnan')./sqrt(nWT);
    semNA = std(thisNA,0,2,'omitnan')./sqrt(nNA);

    % Plot mean FFT with SEM
    nexttile([2 1]); hold on
    plot(freq_fft,meanKIR,'Color',settings.geneColor{1},'LineWidth',settings.lwTri)
    sem1 = patch([freq_fft'; flipud(freq_fft')],[meanKIR-semKIR; flipud(meanKIR+semKIR)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sem1.FaceColor = settings.geneColor{1};
    plot(freq_fft,meanWT,'Color',settings.geneColor{2},'LineWidth',settings.lwTri)
    sem2 = patch([freq_fft'; flipud(freq_fft')],[meanWT-semWT; flipud(meanWT+semWT)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sem2.FaceColor = settings.geneColor{2};
    plot(freq_fft,meanNA,'Color',settings.geneColor{3},'LineWidth',settings.lwTri)
    sem2 = patch([freq_fft'; flipud(freq_fft')],[meanNA-semNA; flipud(meanNA+semNA)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sem2.FaceColor = settings.geneColor{3};
    xlim(fftLim); ylim(yrange); xlabel('Freq (Hz)');
    title([num2str(settings.pursuitGain(c)) 'X'])
    if c==1, ylabel('Angular Velocity Amp (deg/s)'); end
end
sgtitle('FFT')

% Save FFT plot as PNG and SVG
cd(folder.summary)
plotname = 'frequency_ang_FFT';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

% Initialize figure for PSD comparison
figure; set(gcf,'Position',[100 100 1800 800])
tiledlayout(5,nGain,'TileSpacing','compact')
psdLim = [0 5]; yrange = [0 50];

% Plot PSD trials and means for each genotype
for g = 1:3
    switch g
        case 1, thisPSD = kirPSDang; thisN = nKIR;
        case 2, thisPSD = wtPSDang; thisN = nWT;
        case 3, thisPSD = naPSDang; thisN = nNA;
    end
    % Plot each condition's PSD
    for c = 1:nGain
        conPSD = reshape(thisPSD(:,c,:),[],thisN);  % Fetch data
        meanPSD = mean(conPSD,2,'omitnan');  % Calculate mean

        % Plot PSD data and mean
        nexttile; hold on
        plot(freq_psd,10*log10(conPSD),'Color',settings.trialColor,'LineWidth',settings.lwTri)
        plot(freq_psd,10*log10(meanPSD),'Color',settings.geneColor{g},'LineWidth',settings.lwTri)
        xlim(psdLim); ylim(yrange); xlabel('Freq (Hz)')
        if g==1, title([num2str(settings.pursuitGain(c)) 'X']); end
        if c==1, ylabel('Amp ((Deg/s)^2/Hz)'); end
    end
end

% Plot PSD means together for each condition
for c = 1:nGain
    % Fetch data and calculate mean and SEM
    thisKIR = 10*log10(reshape(kirPSDang(:,c,:),[],nKIR));
    thisWT = 10*log10(reshape(wtPSDang(:,c,:),[],nWT));
    thisNA = 10*log10(reshape(naPSDang(:,c,:),[],nNA));
    meanKIR = mean(thisKIR,2,'omitnan');
    meanWT = mean(thisWT,2,'omitnan');
    meanNA = mean(thisNA,2,'omitnan');
    semKIR = std(thisKIR,0,2,'omitnan')./sqrt(nKIR);
    semWT = std(thisWT,0,2,'omitnan')./sqrt(nWT);
    semNA = std(thisNA,0,2,'omitnan')./sqrt(nNA);

    % Plot mean PSD with SEM
    nexttile([2 1]); hold on
    plot(freq_psd,meanKIR,'Color',settings.geneColor{1},'LineWidth',settings.lwTri)
    sem1 = patch([freq_psd; flipud(freq_psd)],[meanKIR-semKIR; flipud(meanKIR+semKIR)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sem1.FaceColor = settings.geneColor{1};
    plot(freq_psd,meanWT,'Color',settings.geneColor{2},'LineWidth',settings.lwTri)
    sem2 = patch([freq_psd; flipud(freq_psd)],[meanWT-semWT; flipud(meanWT+semWT)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sem2.FaceColor = settings.geneColor{2};
    plot(freq_psd,meanNA,'Color',settings.geneColor{3},'LineWidth',settings.lwTri)
    sem3 = patch([freq_psd; flipud(freq_psd)],[meanNA-semNA; flipud(meanNA+semNA)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sem3.FaceColor = settings.geneColor{3};
    xlim(psdLim); ylim(yrange); xlabel('Freq (Hz)');
    title([num2str(settings.pursuitGain(c)) 'X'])
    if c==1, ylabel('Angular Velocity Amplitude ((Deg/s)^2/Hz)'); end
end
legend(settings.geneLabel)
sgtitle('Welchs PSD')

% Save PSD plot as PNG and SVG
cd(folder.summary)
plotname = 'frequency_ang_PSD';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

disp('Complete.')


%% Relationship between error and turning
disp('Comparing relationship between error and turning...')

% Initialize figure for separate plot
figure; set(gcf,'Position',[100 100 1200 800])
tiledlayout(3,nGain,'TileSpacing','compact')
angLim = [-400 400]; angLimA = [-300 300];

% Plot separate trials and means for each genotype
for g = 1:3
    switch g
        case 1, evt = kirEVTRL; thisN = nKIR;
        case 2, evt = wtEVTRL; thisN = nWT;
        case 3, evt = naEVTRL; thisN = nNA;
    end
    meanEVT = mean(evt,3,'omitnan');  % Calculate mean EVT
    % Plot each condition's EVT
    for c = 1:nGain
        nexttile; hold on
        thisTrialEVT = reshape(evt(:,c,:),[],thisN);  % Fetch data
        plot(posBins,thisTrialEVT,'Color',settings.trialColor,'LineWidth',settings.lwTri)
        plot(posBins,meanEVT(:,c),'Color',settings.geneColor{g},'LineWidth',settings.lwAvg)
        axis tight; xline(0); yline(0); xlim([-100 100]); ylim(angLim)
        if g==1, title([num2str(settings.pursuitGain(c)) 'X']); end
        if g==3, xlabel('Target Error (deg)'); end
        if c==1, ylabel('Angular Velocity (deg/s)'); end
    end
end

% Save plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_turn_sep';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

% Initialize figure for combined plot
figure; set(gcf,'Position',[100 100 1700 600])
tiledlayout(1,nGain,'TileSpacing','compact')

% Plot combined data with SEM for each condition
for c = 1:nGain
    % Plot for each genotype (KIR, WT, NA)
    nexttile; hold on
    for g = 1:3
        switch g
            case 1
                evt = kirEVTRL;
                thisN = nKIR;
                evt_color = settings.geneColor{1};  % Color for KIR
            case 2
                evt = wtEVTRL;
                thisN = nWT;
                evt_color = settings.geneColor{2};  % Color for WT
            case 3
                evt = naEVTRL;
                thisN = nNA;
                evt_color = settings.geneColor{3};  % Color for NA
        end

        % Calculate mean and SEM
        meanEVT = mean(evt, 3, 'omitnan');
        semEVT = std(evt, 0, 3, 'omitnan') ./ sqrt(thisN);  % SEM calculation
        validE = sum(~isnan(evt(:, c, :)), 3) > thisN / 3;

        % Plot EVT with SEM
        plot(posBins(validE), meanEVT(validE, c), 'Color', evt_color, 'LineWidth', settings.lwAvg)
        % Plot SEM using patch
        sem_patch = patch([posBins(validE)'; flipud(posBins(validE)')], ...
            [meanEVT(validE, c) - semEVT(validE, c); flipud(meanEVT(validE, c) + semEVT(validE, c))], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
        sem_patch.FaceColor = evt_color;

        axis tight; xline(0); yline(0); xlim([-80 80]); ylim(angLimA);
    end
    title([num2str(settings.pursuitGain(c)) 'X'])
    xlabel('Target Error (deg)')
    if c == 1
        ylabel('Angular Velocity (deg/s)')
    end
end

% Save combined plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_turn_combined';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox,'f');

%% save zoomed version
% Zoom into the plot for each condition
for c = 1:nGain
    nexttile(c); ylim([-130 130]); xlim([-40 40])
end

% Save zoomed plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_turn_combined_20x';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');
disp('Complete.')
%%
% Initialize figure for combined plot
figure; set(gcf,'Position',[100 100 300 600])
tiledlayout(1,1,'TileSpacing','compact')  % Adjusting to a single combined plot

% Combined plot with SEM for each genotype across all gain conditions
nexttile; hold on
for g = 1:3
    switch g
        case 1
            evt = kirEVTRL;
            thisN = nKIR;
            evt_color = settings.geneColor{1};  % Color for KIR
        case 2
            evt = wtEVTRL;
            thisN = nWT;
            evt_color = settings.geneColor{2};  % Color for WT
        case 3
            evt = naEVTRL;
            thisN = nNA;
            evt_color = settings.geneColor{3};  % Color for NA
    end

    % Calculate the mean across gain conditions per animal
    meanEVT_animal = mean(evt, 2, 'omitnan');  % Mean across gains for each animal
    meanEVT = mean(meanEVT_animal, 3, 'omitnan');  % Mean across animals
    semEVT = std(meanEVT_animal, 0, 3, 'omitnan') ./ sqrt(thisN);  % SEM calculation
    meanEVT(isnan(meanEVT)) = 0;
    semEVT(isnan(semEVT)) = 0;

    % Plot EVT with SEM
    plot(posBins, meanEVT, 'Color', evt_color, 'LineWidth', settings.lwAvg)

    % Plot SEM using patch
    sem_patch = patch([posBins'; flipud(posBins')], ...
        [meanEVT - semEVT; flipud(meanEVT + semEVT)], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    sem_patch.FaceColor = evt_color;

    axis tight; xline(0); yline(0); xlim([-80 80]); ylim(angLimA);
end

% Labels and title
title('Combined Gain Condition Averages Across Genotypes')
xlabel('Target Error (deg)')
ylabel('Angular Velocity (deg/s)')

% Save combined plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_turn_combined_across_gains';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox,'f');
% save zoomed version
% Zoom into the plot for each condition
ylim([-90 90]); xlim([-30 30])

% Save zoomed plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_turn_across_gains_20x';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');
disp('Complete.')

%% Fit relationship between turn velocity and object position with ANOVA
disp('Comparing linear fit between object position and turn velocity for each genotype...')

% Set rangeValue for plotting (e.g., -20 to 20)
rangeValue = 20;

% Restrict the posBins to the specified range
rangeMask = (posBins >= -rangeValue) & (posBins <= rangeValue);
posBinsRestricted = posBins(rangeMask);

% Number of gain conditions (nGain)
nGain = size(kirEVT, 2);

% Call the fitting function once to get the outputs for all conditions
[kirSlopes, wtSlopes, naSlopes, kirCI, wtCI, naCI, kirR, wtR, naR] = fitTurnVelocity(posBins, kirEVT, wtEVT, naEVT, rangeValue);
kirSlopes = cell2mat(kirSlopes);
wtSlopes = cell2mat(wtSlopes);
naSlopes = cell2mat(naSlopes);

% Call the general ANOVA function to analyze the contribution of genotype and gain to lag times
if exptFolder == 'AOTU019 KIR'
    [p_evt, ~] = run_genotype_anova_repeated(kirSlopes', wtSlopes', naSlopes', 'Error_V_Turn', folder);
    pval_text_evt = {['p(gene) = ' num2str(p_evt(1))];['p(k) = ' num2str(p_evt(2))]};
else
    [p_evt, ~] = run_genotype_anova1(kirSlopes', wtSlopes', naSlopes', 'Error_V_Turn', folder);
    pval_text_evt = {['p(gene) = ' num2str(p_evt(1))]};
end

% Initialize figure for separate plot
figure; set(gcf,'Position',[100 100 1200 900])  % Adjust width and height for 4 columns, nGain rows
tiledlayout(nGain, 3, 'TileSpacing', 'compact')  % Four columns: fits, slopes, r values, and bootstrapped differences
angLim = [-110 110];

% Colors for genotypes
kirColor = settings.geneColor{1};
wtColor = settings.geneColor{2};
naColor = settings.geneColor{3};

% Loop over each condition to plot linear fits, slope distributions, r values, and bootstrapped differences
for c = 1:nGain
    % Column 1: Linear fits for each condition
    nexttile; hold on

    % Kir group
    kirGroupSlopes = kirSlopes(c,:); % this condition
    kirSlope = mean(kirGroupSlopes);  % Mean slope for plotting fit
    kirCICond = kirCI{c}.slope;  % Confidence interval for the slope
    kirFitLine = kirSlope * posBinsRestricted;  % Fit line

    % WT group
    wtGroupSlopes = wtSlopes(c,:);
    wtSlope = mean(wtGroupSlopes);
    wtCICond = wtCI{c}.slope;
    wtFitLine = wtSlope * posBinsRestricted;

    % NA group
    naGroupSlopes = naSlopes(c,:);
    naSlope = mean(naGroupSlopes);
    naCICond = naCI{c}.slope;
    naFitLine = naSlope * posBinsRestricted;

    % Plot linear fits and confidence intervals using patch
    plot(posBinsRestricted, kirFitLine, 'Color', kirColor, 'LineWidth', settings.lwAvg)
    kir_patch = patch([posBinsRestricted'; flipud(posBinsRestricted')], ...
        [(kirFitLine + kirCICond(2))'; flipud((kirFitLine - kirCICond(1))')], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    kir_patch.FaceColor = kirColor;

    plot(posBinsRestricted, wtFitLine, 'Color', wtColor, 'LineWidth', settings.lwAvg)
    wt_patch = patch([posBinsRestricted'; flipud(posBinsRestricted')], ...
        [(wtFitLine + wtCICond(2))'; flipud((wtFitLine - wtCICond(1))')], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    wt_patch.FaceColor = wtColor;

    plot(posBinsRestricted, naFitLine, 'Color', naColor, 'LineWidth', settings.lwAvg)
    na_patch = patch([posBinsRestricted'; flipud(posBinsRestricted')], ...
        [(naFitLine + naCICond(2))'; flipud((naFitLine - naCICond(1))')], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    na_patch.FaceColor = naColor;

    % Set axis limits and labels
    xlim([-rangeValue rangeValue]);  % Show the specified range
    ylim(angLim)
    xline(0); yline(0);
    if c == nGain, xlabel('Object Position (deg)'); end
    ylabel('AngVel (deg/s)')
    title([num2str(settings.pursuitGain(c)) 'X - Linear Fits']);

    % Column 2: KDE of slopes for each genotype
    nexttile; hold on
    bandwidth = 0.25;  % Adjust this value as needed

    % Plot KDE for each genotype (distributions of slopes for each condition)
    [fKir, xiKir] = ksdensity(kirGroupSlopes, 'Bandwidth', bandwidth);  % KDE for KIR group slopes
    [fWt, xiWt] = ksdensity(wtGroupSlopes, 'Bandwidth', bandwidth);    % KDE for WT group slopes
    [fNa, xiNa] = ksdensity(naGroupSlopes, 'Bandwidth', bandwidth);    % KDE for NA group slopes

    % Plot the KDE curves for slopes
    plot(xiKir, fKir, 'Color', kirColor, 'LineWidth', settings.lwAvg)
    plot(xiWt, fWt, 'Color', wtColor, 'LineWidth', settings.lwAvg)
    plot(xiNa, fNa, 'Color', naColor, 'LineWidth', settings.lwAvg)
    xlim([0 7])

    % Add labels and title
    if c == nGain, xlabel('Slope (Angular Velocity / Object Position)'); end
    ylabel('Density')
    title([num2str(settings.pursuitGain(c)) 'X - Slope KDE']);

    % Column 3: KDE of r values for each genotype
    nexttile; hold on
    bandwidth = 0.0075;  % Adjust this value as needed

    % Plot KDE for each genotype (distributions of r values for each condition)
    [fKirR, xiKirR] = ksdensity(kirR{c}, 'Bandwidth', bandwidth);  % KDE for KIR group r values
    [fWtR, xiWtR] = ksdensity(wtR{c}, 'Bandwidth', bandwidth);    % KDE for WT group r values
    [fNaR, xiNaR] = ksdensity(naR{c}, 'Bandwidth', bandwidth);    % KDE for NA group r values

    % Plot the KDE curves for r values
    plot(xiKirR, fKirR, 'Color', kirColor, 'LineWidth', settings.lwAvg)
    plot(xiWtR, fWtR, 'Color', wtColor, 'LineWidth', settings.lwAvg)
    plot(xiNaR, fNaR, 'Color', naColor, 'LineWidth', settings.lwAvg)
    xlim([0.4 1.1]); xline(1,':')

    % Add labels and title
    if c == nGain, xlabel('R-Squared (r values)'); end
    ylabel('Density')
    title([num2str(settings.pursuitGain(c)) 'X - R Value KDE']);
end
% Add p-value annotation to the plot (top right corner)
text(0.05, 0.95, pval_text_evt, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

sgtitle('Linear Fit, Slope Distribution, R Value, and ANOVA Each Gain Condition')
% Save plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_turn_fit_sep';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');


%% Fit relationship between object position and turn velocity across gains with ANOVA
disp('Comparing linear fit between object position and turn velocity across conditions for all genotypes...')

% Set rangeValue for plotting (e.g., -20 to 20)
rangeValue = 20;

% Restrict the posBins to the specified range
rangeMask = (posBins >= -rangeValue) & (posBins <= rangeValue);
posBinsRestricted = posBins(rangeMask);

% Initialize figure for separate plot
figure; set(gcf,'Position',[100 100 1200 400])
tiledlayout(1, 4, 'TileSpacing', 'compact')
angLim = [-110 110];

% Call the fitting function once to get the outputs for all conditions
[kirSlopes, wtSlopes, naSlopes, ~, ~, ~, kirR, wtR, naR] = fitTurnVelocity(posBins, kirEVT, wtEVT, naEVT, rangeValue);
% Combine across gains
kirSlopes = cell2mat(kirSlopes)';
wtSlopes = cell2mat(wtSlopes)';
naSlopes = cell2mat(naSlopes)';
kirSlopeMeans = mean(kirSlopes,2);
wtSlopeMeans = mean(wtSlopes,2);
naSlopeMeans = mean(naSlopes,2);
kirR = mean(cell2mat(kirR)',2);
wtR = mean(cell2mat(wtR)',2);
naR = mean(cell2mat(naR)',2);
% Compute CI
confidence_level = 0.95;
[kirCI_lower, kirCI_upper] = calculate_confidence_intervals(kirSlopeMeans, confidence_level);
[wtCI_lower, wtCI_upper] = calculate_confidence_intervals(wtSlopeMeans, confidence_level);
[naCI_lower, naCI_upper] = calculate_confidence_intervals(naSlopeMeans, confidence_level);

% Call the general ANOVA function to analyze the contribution of genotype and gain to lag times
[p_evt, ~] = run_genotype_anova1(kirSlopeMeans, wtSlopeMeans, naSlopeMeans, 'Error_V_TurnAll', folder);
% Store p-values from ANOVA for the visual-motor lag plot
pval_text_evt = ['p(gene) = ' num2str(p_evt(1))];

% Colors for genotypes
kirColor = settings.geneColor{1};
wtColor = settings.geneColor{2};
naColor = settings.geneColor{3};

% Tile 1: Plot the results for all genotypes in one tile
nexttile; hold on

% Kir group
kirFitLine = mean(kirSlopeMeans) * posBinsRestricted;  % Use mean slope for KIR group
plot(posBinsRestricted, kirFitLine, 'Color', kirColor, 'LineWidth', settings.lwAvg)
kirCI_upper = kirFitLine + kirCI_upper;
kirCI_lower = kirFitLine - kirCI_lower;
% Plot the SEM patch for KIR group
kir_patch = patch([posBinsRestricted'; flipud(posBinsRestricted')], ...
    [kirCI_upper'; flipud(kirCI_lower')], ...
    'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
kir_patch.FaceColor = kirColor;
% WT group
wtFitLine = mean(wtSlopeMeans) * posBinsRestricted;  % Use mean slope for WT group
plot(posBinsRestricted, wtFitLine, 'Color', wtColor, 'LineWidth', settings.lwAvg)
wtCI_upper = wtFitLine + wtCI_upper;
wtCI_lower = wtFitLine - wtCI_lower;
% Plot the SEM patch for WT group
wt_patch = patch([posBinsRestricted'; flipud(posBinsRestricted')], ...
    [wtCI_upper'; flipud(wtCI_lower')], ...
    'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
wt_patch.FaceColor = wtColor;
% NA group
naFitLine = mean(naSlopeMeans) * posBinsRestricted;  % Use mean slope for NA group
plot(posBinsRestricted, naFitLine, 'Color', naColor, 'LineWidth', settings.lwAvg)
naCI_upper = naFitLine + naCI_upper;
naCI_lower = naFitLine - naCI_lower;
% Plot the SEM patch for NA group
na_patch = patch([posBinsRestricted'; flipud(posBinsRestricted')], ...
    [naCI_upper'; flipud(naCI_lower')], ...
    'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
na_patch.FaceColor = naColor;

% Set axis limits and labels
xlim([-rangeValue rangeValue]);  % Show the specified range
ylim(angLim)
xline(0); yline(0);
xlabel('Object Position (deg)')
ylabel('Angular Velocity (deg/s)')
title('Linear Fit Across Conditions with CI')

% Add p-value annotation to the plot (top right corner)
text(0.05, 0.95, pval_text_evt, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Tile 2: Plot the kernel density estimate (KDE) of slopes for each genotype
nexttile; hold on
bandwidth = 0.25;  % Adjust this value as needed
% Plot KDE for each genotype (slopes)
% Define genotypes and data
genotypes = {'KIR', 'WT', 'NA'};
slopeMeans = {kirSlopeMeans, wtSlopeMeans, naSlopeMeans};

% Define jitter and median marker properties
jitterAmount = 0.1; % Adjust as needed for visibility
medianColor = 'k';  % Black color for median marker
% Loop over each genotype to plot scatter with jitter and median
for i = 1:length(genotypes)
    % Add jitter to x-position for each data point
    xJitter = i + jitterAmount * (rand(size(slopeMeans{i})) - 0.5);

    % Scatter plot with grey '.' markers
    scatter(xJitter, slopeMeans{i}, 'Marker', '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);

    % Plot median with a black diamond marker
    medianValue = median(slopeMeans{i});
    plot(i, medianValue, '_', 'MarkerFaceColor', medianColor, 'MarkerEdgeColor', medianColor);
end

% Set x-ticks and labels for genotypes
set(gca, 'XTick', 1:length(genotypes), 'XTickLabel', genotypes);
xlim([0 4])
% Labels and formatting
xlabel('Genotype');
ylabel('Slope Means');
title('Scatter Plot of Slope Means with Jitter and Median');
hold off;


% Tile 3: Plot the kernel density estimate (KDE) of slopes for each genotype
nexttile; hold on
bandwidth = 0.25;  % Adjust this value as needed
% Plot KDE for each genotype (slopes)
[fKir, xiKir] = ksdensity(kirSlopeMeans, 'Bandwidth', bandwidth);  % KDE for KIR group
[fWt, xiWt] = ksdensity(wtSlopeMeans, 'Bandwidth', bandwidth);    % KDE for WT group
[fNa, xiNa] = ksdensity(naSlopeMeans, 'Bandwidth', bandwidth);    % KDE for NA group

% Plot the KDE curves for slopes
plot(xiKir, fKir, 'Color', kirColor, 'LineWidth', settings.lwAvg)
plot(xiWt, fWt, 'Color', wtColor, 'LineWidth', settings.lwAvg)
plot(xiNa, fNa, 'Color', naColor, 'LineWidth', settings.lwAvg)
% Add labels and title
xlabel('Slope (Angular Velocity / Object Position)')
ylabel('Density')
title('Slope Distribution (Kernel Density)')

% Tile 4: Plot the KDE of r values for each genotype
nexttile; hold on
bandwidth = 0.0075;  % Adjust this value as needed
% Plot KDE for each genotype (r values)
[fKirR, xiKirR] = ksdensity(kirR, 'Bandwidth', bandwidth);  % KDE for KIR group r values
[fWtR, xiWtR] = ksdensity(wtR, 'Bandwidth', bandwidth);    % KDE for WT group r values
[fNaR, xiNaR] = ksdensity(naR, 'Bandwidth', bandwidth);    % KDE for NA group r values
% Plot the KDE curves for r values
plot(xiKirR, fKirR, 'Color', kirColor, 'LineWidth', settings.lwAvg)
plot(xiWtR, fWtR, 'Color', wtColor, 'LineWidth', settings.lwAvg)
plot(xiNaR, fNaR, 'Color', naColor, 'LineWidth', settings.lwAvg)
% Add labels and title
xline(1,':')
xlabel('R-Squared (r values)')
ylabel('Density')
title('R Value Distribution (Kernel Density)')

% Save plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_turn_fit_across';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');



%% Relationship between error velocity and turning
disp('Comparing relationship between error velocity and turning...')

% initialize
figure; set(gcf,'Position',[100 100 1700 800])
tiledlayout(3,nGain,'TileSpacing','compact')
angLim = [-100 100];
angLimA = [-40 40];

% plot separate
% for each genotype
for g = 1:3
    switch g
        case 1
            evt = kirVVT;
            thisN = nKIR;
        case 2
            evt = wtVVT;
            thisN = nWT;
        case 3
            evt = naVVT;
            thisN = nNA;
    end
    meanEVT = mean(evt,3,'omitnan');
    % for each condition
    for c = 1:nGain
        nexttile; hold on
        thisTrialEVT = reshape(evt(:,c,:),[],thisN);
        plot(velBins,thisTrialEVT,'Color',settings.trialColor,'LineWidth',settings.lwTri)
        plot(velBins,meanEVT(:,c),'Color',settings.geneColor{g},'LineWidth',settings.lwAvg)
        axis tight; xline(0);yline(0);ylim(angLim); xlim([-400 400])
        if g==1
            title([num2str(settings.pursuitGain(c)) 'X'])
        elseif g==3
            xlabel('Target Velocity (deg/s)')
        end
        if c==1
            ylabel('Angular Velocity (deg/s)')
        end
    end
end
sgtitle('Target Velocity v Turning for Full FOV')
% save plot
cd(folder.summary)
plotname = 'velocity_v_turn_sep';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

% initialize
figure; set(gcf,'Position',[100 100 1700 800])
tiledlayout(3,nGain,'TileSpacing','compact')
angLim = [-100 100];
angLimA = [-40 40];

% plot separate
% for each genotype
for g = 1:3
    switch g
        case 1
            evt = kirVVTfront;
            thisN = nKIR;
        case 2
            evt = wtVVTfront;
            thisN = nWT;
        case 3
            evt = naVVTfront;
            thisN = nNA;
    end
    meanEVT = mean(evt,3,'omitnan');
    % for each condition
    for c = 1:nGain
        nexttile; hold on
        thisTrialEVT = reshape(evt(:,c,:),[],thisN);
        plot(velBins,thisTrialEVT,'Color',settings.trialColor,'LineWidth',settings.lwTri)
        plot(velBins,meanEVT(:,c),'Color',settings.geneColor{g},'LineWidth',settings.lwAvg)
        axis tight; xline(0);yline(0);ylim(angLim); xlim([-400 400])
        if g==1
            title([num2str(settings.pursuitGain(c)) 'X'])
        elseif g==3
            xlabel('Target Velocity (deg/s)')
        end
        if c==1
            ylabel('Angular Velocity (deg/s)')
        end
    end
end
sgtitle('Target Velocity v Turning for Front FOV')
% save plot
cd(folder.summary)
plotname = 'velocity_v_turn_sepfront';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

% initialize
figure; set(gcf,'Position',[100 100 1700 900])
tiledlayout(2,nGain,'TileSpacing','compact')

% plot together
% for each condition
for c = 1:nGain
    nexttile; hold on
    % for each genotype
    for g = 1:3
        switch g
            case 1
                evt = kirVVT;
                thisN = nKIR;
            case 2
                evt = wtVVT;
                thisN = nWT;
            case 3
                evt = naVVT;
                thisN = nNA;
        end
        meanEVT = mean(evt,3,'omitnan');
        semEVT = std(evt,0,3,'omitnan')./sqrt(thisN);
        % fetch valid points
        validE = sum(~isnan(evt(:,c,:)),3)>thisN/3;

        plot(velBins(validE),meanEVT(validE,c),'Color',settings.geneColor{g},'LineWidth',settings.lwAvg)
        sem = patch([velBins(validE)'; flipud(velBins(validE)')],[meanEVT(validE,c)-semEVT(validE,c); flipud(meanEVT(validE,c)+semEVT(validE,c))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
        sem.FaceColor = settings.geneColor{g};
        axis tight; xline(0);yline(0);ylim(angLimA); xlim([-300 300])
    end
    title([num2str(settings.pursuitGain(c)) 'X all'])
    xlabel('Target Velocity (deg/s)')
    if c==1
        ylabel('Angular Velocity (deg/s)')
    end
end
% plot together for front FOV
% for each condition
for c = 1:nGain
    nexttile; hold on
    % for each genotype
    for g = 1:3
        switch g
            case 1
                evt = kirVVTfront;
                thisN = nKIR;
            case 2
                evt = wtVVTfront;
                thisN = nWT;
            case 3
                evt = naVVTfront;
                thisN = nNA;
        end
        meanEVT = mean(evt,3,'omitnan');
        semEVT = std(evt,0,3,'omitnan')./sqrt(thisN);
        % fetch valid points
        validE = sum(~isnan(evt(:,c,:)),3)>thisN/3;

        plot(velBins(validE),meanEVT(validE,c),'Color',settings.geneColor{g},'LineWidth',settings.lwAvg)
        sem = patch([velBins(validE)'; flipud(velBins(validE)')],[meanEVT(validE,c)-semEVT(validE,c); flipud(meanEVT(validE,c)+semEVT(validE,c))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
        sem.FaceColor = settings.geneColor{g};
        axis tight; xline(0);yline(0);ylim(angLimA); xlim([-300 300])
    end
    title([num2str(settings.pursuitGain(c)) 'X frontFOV'])
    xlabel('Target Velocity (deg/s)')
    if c==1
        ylabel('Angular Velocity (deg/s)')
    end
end
% save plot
cd(folder.summary)
plotname = 'velocity_v_turn_combined';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

disp('Complete.')

%% Combined error velocity across gain conditions
% Initialize figure for combined plot
figure; set(gcf,'Position',[100 100 600 600])
tiledlayout(1,1,'TileSpacing','compact')  % Adjusting to a single combined plot
angLimA = [-60 60];

% Combined plot with SEM for each genotype across all gain conditions
nexttile; hold on
storeVVT = {};

for g = 1:3
    switch g
        case 1
            evt = kirVVT;
            thisN = nKIR;
            evt_color = settings.geneColor{1};  % Color for KIR
        case 2
            evt = wtVVT;
            thisN = nWT;
            evt_color = settings.geneColor{2};  % Color for WT
        case 3
            evt = naVVT;
            thisN = nNA;
            evt_color = settings.geneColor{3};  % Color for NA
    end

    % Calculate the mean across gain conditions per animal
    meanVVT_animal = mean(evt, 2, 'omitnan');  % Mean across gains for each animal
    meanVVT = mean(meanVVT_animal, 3, 'omitnan');  % Mean across animals
    semVVT = std(meanVVT_animal, 0, 3, 'omitnan') ./ sqrt(thisN);  % SEM calculation
    meanVVT(isnan(meanVVT)) = 0;
    semVVT(isnan(semVVT)) = 0;

    % Store
    sz = size(meanVVT_animal);
    idx = find(sz == 1, 1, 'first');  % find singleton dimension
    newOrder = [setdiff(1:numel(sz), idx, 'stable'), idx];
    storeVVT{g} = permute(meanVVT_animal, newOrder);

    % Plot EVT with SEM
    plot(velBins, meanVVT, 'Color', evt_color, 'LineWidth', settings.lwAvg)

    % Plot SEM using patch
    sem_patch = patch([velBins'; flipud(velBins')], ...
        [meanVVT - semVVT; flipud(meanVVT + semVVT)], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    sem_patch.FaceColor = evt_color;

    axis tight; xline(0); yline(0); xlim([-500 500]); ylim(angLimA);
end

% Labels and title
title('Combined Gains')
xlabel('Object Velocity(deg/s)')
ylabel('Rotational Velocity (deg/s)')

% Save combined plot as PNG and SVG
cd(folder.summary)
plotname = 'errorvel_v_turn_combined_across_gains';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox,'f');

%% stats
genotypes = {'KIR', 'WT', 'NA'};

allData = [];
for g = 1:numel(storeVVT)
    thisData = storeVVT{g};                 % nBins × nAnimals
    [nBins, nAnimals] = size(thisData);

    % Create long table for this genotype
    T = table;
    T.Bin       = repmat(velBins(:), nAnimals, 1);
    T.Animal    = categorical(repelem((1:nAnimals)', nBins));
    T.Genotype  = categorical(repmat(genotypes(g), nBins * nAnimals, 1));
    T.Response  = thisData(:);

    % Append
    allData = [allData; T];
end

% --- Fit linear mixed-effects model ---
lme = fitlme(allData, 'Response ~ Bin * Genotype + (1|Animal)');

% --- Display results ---
anova(lme)

% Inputs:
% storeVVT{1}=KIR, {2}=WT, {3}=NA  (nBins x nAnimals, NaNs allowed)
% velBins: nBins x 1 vector

genos = {'KIR','WT','NA'};
nBins = numel(velBins);

rows = {};
for b = 1:nBins
    y = []; g = {};
    for gi = 1:3
        v = storeVVT{gi}(b, :)';
        v = v(~isnan(v));
        y = [y; v]; %#ok<AGROW>
        g = [g; repmat(genos(gi), numel(v), 1)]; %#ok<AGROW>
    end
    if numel(unique(g)) < 2 || numel(y) < 3
        continue  % not enough data for this bin
    end

    % One-way ANOVA (no figure)
    [~,~,stats] = anova1(y, g, 'off');

    % Tukey-Kramer within this bin
    C = multcompare(stats, 'ctype', 'tukey-kramer', 'display', 'off'); 
    % C columns: [i j lower diff upper p]
    giNames = stats.gnames;  % order used by multcompare

    for k = 1:size(C,1)
        rows(end+1, :) = { ...
            velBins(b), giNames{C(k,1)}, giNames{C(k,2)}, ...   % Bin, Group1, Group2
            C(k,4), C(k,3), C(k,5), C(k,6)}; %#ok<AGROW>
    end
end

posthoc = cell2table(rows, 'VariableNames', ...
    {'Bin','Group1','Group2','Estimate','CI_Lower','CI_Upper','pValue_Tukey_withinBin'});

p = posthoc.pValue_Tukey_withinBin;
posthoc.p_FDR_BH = fdr_bh(p);
posthoc.sig_FDR_05 = posthoc.p_FDR_BH < 0.05;

%% Relationship between error and acceleration
disp('Comparing relationship between error and acceleration...')

% Initialize figure for separate plot
figure; set(gcf,'Position',[100 100 1200 800])
tiledlayout(3, nGain, 'TileSpacing', 'compact')
accelLim = [-400 400]; accelLimA = [-300 300];

% Plot separate trials and means for each genotype
for g = 1:3
    switch g
        case 1, evacc = kirEVAccRL; thisN = nKIR;
        case 2, evacc = wtEVAccRL; thisN = nWT;
        case 3, evacc = naEVAccRL; thisN = nNA;
    end
    meanEVAcc = mean(evacc, 3, 'omitnan');  % Calculate mean EVAcc
    % Plot each condition's EVAcc
    for c = 1:nGain
        nexttile; hold on
        thisTrialEVAcc = reshape(evacc(:,c,:),[],thisN);  % Fetch data
        plot(posBins, thisTrialEVAcc, 'Color', settings.trialColor, 'LineWidth', settings.lwTri)
        plot(posBins, meanEVAcc(:,c), 'Color', settings.geneColor{g}, 'LineWidth', settings.lwAvg)
        axis tight; xline(0); yline(0); xlim([-100 100]); ylim(accelLim)
        if g == 1, title([num2str(settings.pursuitGain(c)) 'X']); end
        if g == 3, xlabel('Target Error (deg)'); end
        if c == 1, ylabel('Angular Acceleration (deg/s²)'); end
    end
end

% Save plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_accel_sep';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox, 'f');

% Initialize figure for combined plot
figure; set(gcf,'Position',[100 100 1700 600])
tiledlayout(1, nGain, 'TileSpacing', 'compact')

% Plot combined data with SEM for each condition
for c = 1:nGain
    % Plot for each genotype (KIR, WT, NA)
    nexttile; hold on
    for g = 1:3
        switch g
            case 1
                evacc = kirEVAccRL;
                thisN = nKIR;
                evacc_color = settings.geneColor{1};  % Color for KIR
            case 2
                evacc = wtEVAccRL;
                thisN = nWT;
                evacc_color = settings.geneColor{2};  % Color for WT
            case 3
                evacc = naEVAccRL;
                thisN = nNA;
                evacc_color = settings.geneColor{3};  % Color for NA
        end

        % Calculate mean and SEM
        meanEVAcc = mean(evacc, 3, 'omitnan');
        semEVAcc = std(evacc, 0, 3, 'omitnan') ./ sqrt(thisN);  % SEM calculation
        validE = sum(~isnan(evacc(:, c, :)), 3) > thisN / 3;

        % Plot EVAcc with SEM
        plot(posBins(validE), meanEVAcc(validE, c), 'Color', evacc_color, 'LineWidth', settings.lwAvg)
        % Plot SEM using patch
        sem_patch = patch([posBins(validE)'; flipud(posBins(validE)')], ...
            [meanEVAcc(validE, c) - semEVAcc(validE, c); flipud(meanEVAcc(validE, c) + semEVAcc(validE, c))], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
        sem_patch.FaceColor = evacc_color;

        axis tight; xline(0); yline(0); xlim([-80 80]); ylim(accelLimA);
    end
    title([num2str(settings.pursuitGain(c)) 'X'])
    xlabel('Target Error (deg)')
    if c == 1
        ylabel('Angular Acceleration (deg/s²)')
    end
end

% Save combined plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_accel_combined';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox, 'f');

% Save zoomed version
% Zoom into the plot for each condition
for c = 1:nGain
    nexttile(c); ylim([-130 130]); xlim([-40 40])
end

% Save zoomed plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_accel_combined_20x';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.')


%% Direction change v object crossing binned by crossing
lagrange = [80 220];
angrange = [0 220];
dc_velbins = bins.ang2;

if exptFolder == 'AOTU019 KIR'
    % compare relationship between max turn and time of crossing binned
    % Initialize figure for combined plot
    figure; set(gcf,'Position',[100 100 1200 900])
    tiledlayout(5, nGain, 'TileSpacing', 'compact')

    % Plot separately for each group
    for g = 1:3
        % Fetch data for each group
        switch g
            case 1
                thisBinned = kirDCAngCrossBins;
                thisN = nKIR;
            case 2
                thisBinned = wtDCAngCrossBins;
                thisN = nWT;
            case 3
                thisBinned = naDCAngCrossBins;
                thisN = nNA;
        end

        % Calculate mean
        meanBinned = mean(thisBinned, 3, 'omitnan');

        % Plot for each condition
        for c = 1:nGain
            nexttile; hold on
            plot(dc_velbins, reshape(thisBinned(:,c,:), [], thisN), 'Color', settings.trialColor, 'Linewidth', settings.lwTri)
            plot(dc_velbins, meanBinned(:, c), 'Color', settings.geneColor{g}, 'Linewidth', settings.lwAvg)
            xlim(angrange); ylim(lagrange)
            if c == 1
                ylabel({'Direction Change (ms)'})
            end
            if g == 1
                title([num2str(settings.pursuitGain(c)) 'X'])
            end
        end
    end

    % Plot together
    % Calculate mean and SEM across gains for each genotype
    meanKIR = mean(kirDCAngCrossBins, 3, 'omitnan');
    meanWT = mean(wtDCAngCrossBins, 3, 'omitnan');
    meanNA = mean(naDCAngCrossBins, 3, 'omitnan');
    semKIR = std(kirDCAngCrossBins, 0, 3, 'omitnan') ./ sqrt(nKIR);
    semWT = std(wtDCAngCrossBins, 0, 3, 'omitnan') ./ sqrt(nWT);
    semNA = std(naDCAngCrossBins, 0, 3, 'omitnan') ./ sqrt(nNA);

    % Run three-way ANOVA on DCAngBins data (assuming each genotype is in the corresponding variable)
    [p, ~] = run_genotype_anova3_repeated(kirDCAngCrossBins, wtDCAngCrossBins, naDCAngCrossBins, dc_velbins, 'DirectionChangeSep', folder);
    pval_text = {['p(gene) = ' num2str(p(1)) ] [ 'p(gain) = ' num2str(p(2)) ] [ 'p(angbin) = ' num2str(p(3))]};

    % Plot setpoint performance metrics
    for c = 1:nGain
        nexttile([2, 1]); hold on
        % Plot KIR
        plot(dc_velbins(:), meanKIR(:, c), 'Color', settings.geneColor{1}, 'Linewidth', settings.lwAvg)
        sem1 = patch([dc_velbins(:); flipud(dc_velbins(:))], ...
            [meanKIR(:, c) - semKIR(:, c); flipud(meanKIR(:, c) + semKIR(:, c))], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
        sem1.FaceColor = settings.geneColor{1};

        % Plot WT
        plot(dc_velbins(:), meanWT(:, c), 'Color', settings.geneColor{2}, 'Linewidth', settings.lwAvg)
        sem2 = patch([dc_velbins(:); flipud(dc_velbins(:))], ...
            [meanWT(:, c) - semWT(:, c); flipud(meanWT(:, c) + semWT(:, c))], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
        sem2.FaceColor = settings.geneColor{2};

        % Plot NA
        plot(dc_velbins(:), meanNA(:, c), 'Color', settings.geneColor{3}, 'Linewidth', settings.lwAvg)
        sem3 = patch([dc_velbins(:); flipud(dc_velbins(:))], ...
            [meanNA(:, c) - semNA(:, c); flipud(meanNA(:, c) + semNA(:, c))], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
        sem3.FaceColor = settings.geneColor{3};
        xlim(angrange); ylim(lagrange)
        title([num2str(settings.pursuitGain(c)) 'X'])
        if c == 1
            ylabel({'Time to Angular Velocity', 'Direction Change (ms)'})
        end
        xlabel({'Rotational Velocity', 'at Obj Crossing (deg/s)'})
    end

    % Add p-value annotation to the last plot (top left corner)
    text(0.05, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 7, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    sgtitle({'Time to Angular Velocity Direction Change After Object Crossing Midline', 'Binned by Rotational Velocity at Object Crossing'})

    % Save plot
    cd(folder.summary)
    plotname = 'dirchange_timing_angcrossbin';
    saveas(gcf, [plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox, 'f');
    % Save vectorized plot
    cd(folder.vector)
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox, 'f');
end

% Initialize figure for combined
figure; set(gcf,'Position',[100 100 500 900])
tiledlayout(5, 2, 'TileSpacing', 'compact') % 3 columns for separate genotypes + combined
dc_velbins = bins.ang2;

% Average across gains for each genotype
meanKIRDC = reshape(mean(kirDCAngCrossBins, 2, 'omitnan'), [], nKIR); % Reshape for KIR
meanWTDC = reshape(mean(wtDCAngCrossBins, 2, 'omitnan'), [], nWT); % Reshape for WT
meanNADC = reshape(mean(naDCAngCrossBins, 2, 'omitnan'), [], nNA); % Reshape for NA

% Plot each genotype separately
for g = 1:3
    % Fetch data for each group
    switch g
        case 1
            thisBinned = meanKIRDC;
            thisN = nKIR;
            thisColor = settings.geneColor{1};
            titleName = 'KIR';
        case 2
            thisBinned = meanWTDC;
            thisN = nWT;
            thisColor = settings.geneColor{2};
            titleName = 'WT';
        case 3
            thisBinned = meanNADC;
            thisN = nNA;
            thisColor = settings.geneColor{3};
            titleName = 'NA';
    end
    % Calculate mean across gain conditions
    meanBinned = mean(thisBinned, 2, 'omitnan');  % Assuming allMeanBinned is a 3D array

    % Plot across gain conditions for the current genotype
    nexttile; hold on
    plot(dc_velbins, reshape(thisBinned, [], thisN), 'Color', settings.trialColor, 'LineWidth', settings.lwTri);
    plot(dc_velbins, meanBinned, 'Color', thisColor, 'LineWidth', settings.lwAvg);
    xlim(angrange); ylim(lagrange);
    set(gca, 'XTick', dc_velbins)
    ylabel('Direction Change (ms)');
    title(titleName);

    % Plot across gain conditions for the current genotype
    nexttile; hold on
    plot(dc_velbins, reshape(thisBinned, [], thisN), 'Color', settings.trialColor, 'LineWidth', settings.lwTri);
    plot(dc_velbins, meanBinned, 'Color', thisColor, 'LineWidth', settings.lwAvg);
    xlim(angrange); ylim(lagrange);
    set(gca, 'XTick', dc_velbins, 'XScale', 'log'); % Set x-axis to logarithmic scale
    title(titleName);
end
sgtitle('Direction Change by Genotype');  % Add a common title

[p, ~] = run_genotype_anova_repeated(meanKIRDC', meanWTDC', meanNADC', 'DirectionChangeAcross', folder);

% Plot all genotypes together
% Calculate mean and SEM across genotypes
meanKIR_across = mean(meanKIRDC, 2, 'omitnan');
meanWT_across = mean(meanWTDC, 2, 'omitnan');
meanNA_across = mean(meanNADC, 2, 'omitnan');
semKIR_across = std(meanKIRDC, 0, 2, 'omitnan') ./ sqrt(nKIR);
semWT_across = std(meanWTDC, 0, 2, 'omitnan') ./ sqrt(nWT);
semNA_across = std(meanNADC, 0, 2, 'omitnan') ./ sqrt(nNA);

% Plot all genotypes together in a new tile
for t = 1:2
    nexttile([2, 1]); hold on
    % Plot KIR
    plot(dc_velbins, meanKIR_across, 'Color', settings.geneColor{1}, 'Linewidth', settings.lwAvg)
    semPatchKIR = patch([dc_velbins'; flipud(dc_velbins')], ...
        [meanKIR_across - semKIR_across; flipud(meanKIR_across + semKIR_across)], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    semPatchKIR.FaceColor = settings.geneColor{1};

    % Plot WT
    plot(dc_velbins, meanWT_across, 'Color', settings.geneColor{2}, 'Linewidth', settings.lwAvg)
    semPatchWT = patch([dc_velbins'; flipud(dc_velbins')], ...
        [meanWT_across - semWT_across; flipud(meanWT_across + semWT_across)], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    semPatchWT.FaceColor = settings.geneColor{2};

    % Plot NA
    plot(dc_velbins, meanNA_across, 'Color', settings.geneColor{3}, 'Linewidth', settings.lwAvg)
    semPatchNA = patch([dc_velbins'; flipud(dc_velbins')], ...
        [meanNA_across - semNA_across; flipud(meanNA_across + semNA_across)], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    semPatchNA.FaceColor = settings.geneColor{3};

    % Add p-value annotation to the plot (top left corner)
    pval_text = {['p(gene) = ' num2str(p(1))]; ['p(angbin) = ' num2str(p(2))]};  % Create text for p-values
    text(0.05, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    xlim(angrange); ylim(lagrange)
    title('All Gains')
    if t==1
        ylabel('Direction Change (ms)')
    else
        set(gca, 'XTick', dc_velbins, 'XScale', 'log'); % Set x-axis to logarithmic scale
    end
    xlabel('Rotational Velocity at Cross (deg/s)')
end

% Save plot
cd(folder.summary)
plotname = 'dirchange_timing_angcrossbin_across';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');
% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');


% Generate log fit
fit_log_2_dirchange(meanKIRDC, meanWTDC, meanNADC, settings)
% Save plot
cd(folder.summary)
plotname = 'dirchange_timing_angcrossbin_fitacross';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');
% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');

%% Direction change v object crossing binned by prior
lagrange = [80 220];
angrange = [0 220];
dc_velbins = bins.ang2;

if exptFolder == 'AOTU019 KIR'
    % compare relationship between max turn and time of crossing binned
    % Initialize figure for combined plot
    figure; set(gcf,'Position',[100 100 1200 900])
    tiledlayout(5, nGain, 'TileSpacing', 'compact')

    % Plot separately for each group
    for g = 1:3
        % Fetch data for each group
        switch g
            case 1
                thisBinned = kirDCAngPriorBins;
                thisN = nKIR;
            case 2
                thisBinned = wtDCAngPriorBins;
                thisN = nWT;
            case 3
                thisBinned = naDCAngPriorBins;
                thisN = nNA;
        end

        % Calculate mean
        meanBinned = mean(thisBinned, 3, 'omitnan');

        % Plot for each condition
        for c = 1:nGain
            nexttile; hold on
            plot(dc_velbins, reshape(thisBinned(:,c,:), [], thisN), 'Color', settings.trialColor, 'Linewidth', settings.lwTri)
            plot(dc_velbins, meanBinned(:, c), 'Color', settings.geneColor{g}, 'Linewidth', settings.lwAvg)
            xlim(angrange); ylim(lagrange)
            set(gca, 'XScale', 'log')  % Set x-axis to log scale
            set(gca, 'XTick', dc_velbins)
            if c == 1
                ylabel({'Direction Change (ms)'})
            end
            if g == 1
                title([num2str(settings.pursuitGain(c)) 'X'])
            end
        end
    end

    % Plot together
    % Calculate mean and SEM across gains for each genotype
    meanKIR = mean(kirDCAngPriorBins, 3, 'omitnan');
    meanWT = mean(wtDCAngPriorBins, 3, 'omitnan');
    meanNA = mean(naDCAngPriorBins, 3, 'omitnan');
    semKIR = std(kirDCAngPriorBins, 0, 3, 'omitnan') ./ sqrt(nKIR);
    semWT = std(wtDCAngPriorBins, 0, 3, 'omitnan') ./ sqrt(nWT);
    semNA = std(naDCAngPriorBins, 0, 3, 'omitnan') ./ sqrt(nNA);

    % Run three-way ANOVA on DCAngBins data (assuming each genotype is in the corresponding variable)
    [p, ~] = run_genotype_anova3_repeated(kirDCAngPriorBins, wtDCAngPriorBins, naDCAngPriorBins, dc_velbins, 'DirectionChangePriorSep', folder);

    % Store p-values from ANOVA
    pval_text = {['p(gene) = ' num2str(p(1)) ] [ 'p(gain) = ' num2str(p(2)) ] [ 'p(angbin) = ' num2str(p(3))]};

    % Plot setpoint performance metrics
    for c = 1:nGain
        validkir = sum(~isnan(kirDCAngPriorBins(:,c,:)), 3) > nKIR / 3;
        validwt = sum(~isnan(wtDCAngPriorBins(:,c,:)), 3) > nWT / 3;
        validna = sum(~isnan(naDCAngPriorBins(:,c,:)), 3) > nNA / 3;

        nexttile([2, 1]); hold on
        % Plot KIR
        plot(dc_velbins(validkir), meanKIR(validkir, c), 'Color', settings.geneColor{1}, 'Linewidth', settings.lwAvg)
        sem1 = patch([dc_velbins(validkir)'; flipud(dc_velbins(validkir)')], ...
            [meanKIR(validkir, c) - semKIR(validkir, c); flipud(meanKIR(validkir, c) + semKIR(validkir, c))], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
        sem1.FaceColor = settings.geneColor{1};

        % Plot WT
        plot(dc_velbins(validwt), meanWT(validwt, c), 'Color', settings.geneColor{2}, 'Linewidth', settings.lwAvg)
        sem2 = patch([dc_velbins(validwt)'; flipud(dc_velbins(validwt)')], ...
            [meanWT(validwt, c) - semWT(validwt, c); flipud(meanWT(validwt, c) + semWT(validwt, c))], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
        sem2.FaceColor = settings.geneColor{2};

        % Plot NA
        plot(dc_velbins(validna), meanNA(validna, c), 'Color', settings.geneColor{3}, 'Linewidth', settings.lwAvg)
        sem3 = patch([dc_velbins(validna)'; flipud(dc_velbins(validna)')], ...
            [meanNA(validna, c) - semNA(validna, c); flipud(meanNA(validna, c) + semNA(validna, c))], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
        sem3.FaceColor = settings.geneColor{3};

        xlim(angrange); ylim(lagrange)
        set(gca, 'XScale', 'log')  % Set x-axis to log scale
        set(gca, 'XTick', dc_velbins)
        title([num2str(settings.pursuitGain(c)) 'X'])
        if c == 1
            ylabel({'Time to Angular Velocity', 'Direction Change (ms)'})
        end
        xlabel({'Rotational Velocity', 'Max Prior (deg/s)'})
    end

    % Add p-value annotation to the last plot (top left corner)
    text(0.05, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 7, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    sgtitle({'Time to Angular Velocity Direction Change After Object Crossing Midline', 'Binned by Rotational Velocity Max of Prior Turn'})

    % Save plot
    cd(folder.summary)
    plotname = 'dirchange_timing_angpriorbin';
    saveas(gcf, [plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox, 'f');
    % Save vectorized plot
    cd(folder.vector)
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox, 'f');
end

% Initialize figure for combined
figure; set(gcf,'Position',[100 100 300 900])
tiledlayout(5, 1, 'TileSpacing', 'compact') % 3 columns for separate genotypes + combined

lagrange = [100 200];
angrange = [0 350];
dc_velbins = bins.ang2;

% Average across gains for each genotype
meanKIRDC = reshape(mean(kirDCAngPriorBins, 2, 'omitnan'), [], nKIR); % Reshape for KIR
meanWTDC = reshape(mean(wtDCAngPriorBins, 2, 'omitnan'), [], nWT); % Reshape for WT
meanNADC = reshape(mean(naDCAngPriorBins, 2, 'omitnan'), [], nNA); % Reshape for NA

% Plot each genotype separately
for g = 1:3
    % Fetch data for each group
    switch g
        case 1
            thisBinned = meanKIRDC;
            thisN = nKIR;
            thisColor = settings.geneColor{1};
            titleName = 'KIR';
        case 2
            thisBinned = meanWTDC;
            thisN = nWT;
            thisColor = settings.geneColor{2};
            titleName = 'WT';
        case 3
            thisBinned = meanNADC;
            thisN = nNA;
            thisColor = settings.geneColor{3};
            titleName = 'NA';
    end
    % Calculate mean across gain conditions
    meanBinned = mean(thisBinned, 2, 'omitnan');  % Assuming allMeanBinned is a 3D array

    % Plot across gain conditions for the current genotype
    nexttile; hold on
    plot(dc_velbins, reshape(thisBinned, [], thisN), 'Color', settings.trialColor, 'LineWidth', settings.lwTri);
    plot(dc_velbins, meanBinned, 'Color', thisColor, 'LineWidth', settings.lwAvg);
    xlim(angrange); ylim(lagrange);
    set(gca, 'XScale', 'log')  % Set x-axis to log scale
    set(gca, 'XTick', dc_velbins)
    ylabel('Direction Change (ms)');
    title(titleName);
end
sgtitle('Direction Change by Genotype');  % Add a common title

[p, ~] = run_genotype_anova_repeated(meanKIRDC', meanWTDC', meanNADC', 'DirectionChangePriorAcross', folder);

% Plot all genotypes together
% Calculate mean and SEM across genotypes
meanKIR_across = mean(meanKIRDC, 2, 'omitnan');
meanWT_across = mean(meanWTDC, 2, 'omitnan');
meanNA_across = mean(meanNADC, 2, 'omitnan');
semKIR_across = std(meanKIRDC, 0, 2, 'omitnan') ./ sqrt(nKIR);
semWT_across = std(meanWTDC, 0, 2, 'omitnan') ./ sqrt(nWT);
semNA_across = std(meanNADC, 0, 2, 'omitnan') ./ sqrt(nNA);

% Plot all genotypes together in a new tile
nexttile([2, 1]); hold on
% Plot KIR
plot(dc_velbins, meanKIR_across, 'Color', settings.geneColor{1}, 'Linewidth', settings.lwAvg)
semPatchKIR = patch([dc_velbins'; flipud(dc_velbins')], ...
    [meanKIR_across - semKIR_across; flipud(meanKIR_across + semKIR_across)], ...
    'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
semPatchKIR.FaceColor = settings.geneColor{1};

% Plot WT
plot(dc_velbins, meanWT_across, 'Color', settings.geneColor{2}, 'Linewidth', settings.lwAvg)
semPatchWT = patch([dc_velbins'; flipud(dc_velbins')], ...
    [meanWT_across - semWT_across; flipud(meanWT_across + semWT_across)], ...
    'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
semPatchWT.FaceColor = settings.geneColor{2};

% Plot NA
plot(dc_velbins, meanNA_across, 'Color', settings.geneColor{3}, 'Linewidth', settings.lwAvg)
semPatchNA = patch([dc_velbins'; flipud(dc_velbins')], ...
    [meanNA_across - semNA_across; flipud(meanNA_across + semNA_across)], ...
    'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
semPatchNA.FaceColor = settings.geneColor{3};

% Add p-value annotation to the plot (top left corner)
pval_text = {['p(gene) = ' num2str(p(1))]; ['p(angbin) = ' num2str(p(2))]};  % Create text for p-values
text(0.05, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

xlim(angrange); ylim(lagrange)
set(gca, 'XScale', 'log')  % Set x-axis to log scale
set(gca, 'XTick', dc_velbins)
title('All Gains')
ylabel('Direction Change (ms)')
xlabel('Rotational Velocity Max Prior (deg/s)')

% Save plot
cd(folder.summary)
plotname = 'dirchange_timing_angpriorbin_across';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');
% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');

%% compare relationship between max-to-max time and max angular velocity binned
lagrange = [0 300]; % Time range in ms
angrange = [0 220];   % Max angular velocity range in deg/s
dc_velbins = bins.ang2;

if exptFolder == 'AOTU019 KIR'
    % Initialize figure for combined plot
    figure; set(gcf, 'Position', [100 100 1200 900])
    tiledlayout(5, nGain, 'TileSpacing', 'compact')

    % Plot separately for each group
    for g = 1:3
        % Fetch data for each group
        switch g
            case 1
                thisBinned = kirMCAngBins;
                thisN = nKIR;
            case 2
                thisBinned = wtMCAngBins;
                thisN = nWT;
            case 3
                thisBinned = naMCAngBins;
                thisN = nNA;
        end

        % Calculate mean across animals for each bin and gain condition
        meanBinned = mean(thisBinned, 3, 'omitnan');

        % Plot for each gain condition
        for c = 1:nGain
            nexttile; hold on
            plot(dc_velbins, reshape(thisBinned(:, c, :), [], thisN), 'Color', settings.trialColor, 'Linewidth', settings.lwTri)
            plot(dc_velbins, meanBinned(:, c), 'Color', settings.geneColor{g}, 'Linewidth', settings.lwAvg)
            xlim(angrange); ylim(lagrange)
            if c == 1
                ylabel({'Max-to-Max Time (ms)'})
            end
            if g == 1
                title([num2str(settings.pursuitGain(c)) 'X'])
            end
        end
    end

    % Plot all genotypes together
    % Calculate mean and SEM across gain conditions for each genotype
    meanKIR = mean(kirMCAngBins, 3, 'omitnan');
    meanWT = mean(wtMCAngBins, 3, 'omitnan');
    meanNA = mean(naMCAngBins, 3, 'omitnan');
    semKIR = std(kirMCAngBins, 0, 3, 'omitnan') ./ sqrt(nKIR);
    semWT = std(wtMCAngBins, 0, 3, 'omitnan') ./ sqrt(nWT);
    semNA = std(naMCAngBins, 0, 3, 'omitnan') ./ sqrt(nNA);

    % Run three-way ANOVA on MCAngBins data (assuming each genotype is in the corresponding variable)
    [p, ~] = run_genotype_anova3_repeated(kirMCAngBins, wtMCAngBins, naMCAngBins, dc_velbins, 'MaxToMaxSep', folder);

    % Store p-values from ANOVA
    pval_text = {['p(gene) = ' num2str(p(1))] ['p(gain) = ' num2str(p(2))] ['p(angbin) = ' num2str(p(3))]};

    % Plot combined results
    for c = 1:nGain
        validkir = sum(~isnan(kirMCAngBins(:, c, :)), 3) > nKIR / 3;
        validwt = sum(~isnan(wtMCAngBins(:, c, :)), 3) > nWT / 3;
        validna = sum(~isnan(naMCAngBins(:, c, :)), 3) > nNA / 3;

        nexttile([2, 1]); hold on
        % Plot KIR
        plot(dc_velbins(validkir), meanKIR(validkir, c), 'Color', settings.geneColor{1}, 'Linewidth', settings.lwAvg)
        patch([dc_velbins(validkir)'; flipud(dc_velbins(validkir)')], ...
            [meanKIR(validkir, c) - semKIR(validkir, c); flipud(meanKIR(validkir, c) + semKIR(validkir, c))], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none', 'FaceColor', settings.geneColor{1});

        % Plot WT
        plot(dc_velbins(validwt), meanWT(validwt, c), 'Color', settings.geneColor{2}, 'Linewidth', settings.lwAvg)
        patch([dc_velbins(validwt)'; flipud(dc_velbins(validwt)')], ...
            [meanWT(validwt, c) - semWT(validwt, c); flipud(meanWT(validwt, c) + semWT(validwt, c))], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none', 'FaceColor', settings.geneColor{2});

        % Plot NA
        plot(dc_velbins(validna), meanNA(validna, c), 'Color', settings.geneColor{3}, 'Linewidth', settings.lwAvg)
        patch([dc_velbins(validna)'; flipud(dc_velbins(validna)')], ...
            [meanNA(validna, c) - semNA(validna, c); flipud(meanNA(validna, c) + semNA(validna, c))], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none', 'FaceColor', settings.geneColor{3});

        xlim(angrange); ylim(lagrange)
        title([num2str(settings.pursuitGain(c)) 'X'])
        if c == 1
            ylabel({'Max-to-Max Time (ms)'})
        end
        xlabel({'Max Angular Velocity (deg/s)'})
    end

    % Add p-value annotation to the last plot
    text(0.05, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    sgtitle({'Max-to-Max Time Lag', 'Binned by Max Angular Velocity'})
    % Save plot
    cd(folder.summary)
    plotname = 'max_to_max_timing_angbin';
    saveas(gcf, [plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox, 'f');
    % Save vectorized plot
    cd(folder.vector)
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox, 'f');

    % Plot all genotypes together in a final tile
    figure
    nexttile([2, 1]); hold on
    % Calculate mean and SEM across genotypes
    meanKIR_across = mean(meanKIR, 2, 'omitnan');
    meanWT_across = mean(meanWT, 2, 'omitnan');
    meanNA_across = mean(meanNA, 2, 'omitnan');
    semKIR_across = std(meanKIR, 0, 2, 'omitnan') ./ sqrt(nKIR);
    semWT_across = std(meanWT, 0, 2, 'omitnan') ./ sqrt(nWT);
    semNA_across = std(meanNA, 0, 2, 'omitnan') ./ sqrt(nNA);

    [p, ~] = run_genotype_anova_repeated(meanKIRDC', meanWTDC', meanNADC', 'DirectionChangePriorAcross', folder);

    % Plot KIR
    plot(dc_velbins, meanKIR_across, 'Color', settings.geneColor{1}, 'Linewidth', settings.lwAvg)
    patch([dc_velbins'; flipud(dc_velbins')], ...
        [meanKIR_across - semKIR_across; flipud(meanKIR_across + semKIR_across)], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none', 'FaceColor', settings.geneColor{1});

    % Plot WT
    plot(dc_velbins, meanWT_across, 'Color', settings.geneColor{2}, 'Linewidth', settings.lwAvg)
    patch([dc_velbins'; flipud(dc_velbins')], ...
        [meanWT_across - semWT_across; flipud(meanWT_across + semWT_across)], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none', 'FaceColor', settings.geneColor{2});

    % Plot NA
    plot(dc_velbins, meanNA_across, 'Color', settings.geneColor{3}, 'Linewidth', settings.lwAvg)
    patch([dc_velbins'; flipud(dc_velbins')], ...
        [meanNA_across - semNA_across; flipud(meanNA_across + semNA_across)], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none', 'FaceColor', settings.geneColor{3});

    % Add p-value annotation to the plot (top left corner)
    text(0.05, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    xlim(angrange); ylim(lagrange)
    set(gca, 'XScale', 'log')  % Set x-axis to log scale
    title('All Gains')
    ylabel('Max-to-Max Time (ms)')
    xlabel('Max Angular Velocity (deg/s)')

    % Save final combined plot
    plotname = 'max_to_max_combined';
    saveas(gcf, [plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox, 'f');
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox, 'f');
end

%% Interval between direction changes
if exptFolder == 'AOTU019 KIR'
    % Initialize figure with a tiled layout
    figure;
    tiledlayout(2, 1);  % Two rows, one column
    set(gcf,'Position',[100 100 400 900]);

    % Call the general ANOVA function to analyze the contribution of genotype and gain to general performance
    [p_fixint, ~] = run_genotype_anova_repeated(kirFixInt, wtFixInt, naFixInt, 'IntervalFixation', folder);
    [p_notint, ~] = run_genotype_anova_repeated(kirNotInt, wtNotInt, naNotInt, 'IntervalOther', folder);

    % Format p-value text for display
    pval_text_fixation = ['p(gene) = ' num2str(p_fixint(1)) ', p(gain) = ' num2str(p_fixint(2))];
    pval_text_nonfixation = ['p(gene) = ' num2str(p_notint(1)) ', p(gain) = ' num2str(p_notint(2))];

    % Get the number of gain conditions
    gainConditions = settings.pursuitGain;

    % Loop for plotting fixation and non-fixation times
    for row = 1:2
        if row == 1
            % First row: Fixation times
            dataKir = kirFixInt;
            dataWt = wtFixInt;
            dataNa = naFixInt;
            plotTitle = 'Fixation';
            pval_text = pval_text_fixation;
        else
            % Second row: Non-fixation times
            dataKir = kirNotInt;
            dataWt = wtNotInt;
            dataNa = naNotInt;
            plotTitle = 'Non-Fixation';
            pval_text = pval_text_nonfixation;
        end

        % Select the current tile
        nexttile;
        hold on;

        % Plot KIR genotype
        xVals = gainConditions - 2;  % Offset for KIR
        color = settings.geneColor{1};
        for cond = 1:length(gainConditions)
            scatter(repmat(xVals(cond), size(dataKir, 1), 1), dataKir(:, cond), 40, settings.trialColor, '.');
            meanVal = mean(dataKir(:, cond));
            semVal = std(dataKir(:, cond)) / sqrt(size(dataKir, 1));
            errorbar(xVals(cond), meanVal, semVal, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'Color', color, 'CapSize', 10, 'LineWidth', 1);
        end

        % Plot WT genotype
        xVals = gainConditions;  % No offset for WT
        color = settings.geneColor{2};
        for cond = 1:length(gainConditions)
            scatter(repmat(xVals(cond), size(dataWt, 1), 1), dataWt(:, cond), 40, settings.trialColor, '.');
            meanVal = mean(dataWt(:, cond));
            semVal = std(dataWt(:, cond)) / sqrt(size(dataWt, 1));
            errorbar(xVals(cond), meanVal, semVal, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'Color', color, 'CapSize', 10, 'LineWidth', 1);
        end

        % Plot NA genotype
        xVals = gainConditions + 2;  % Offset for NA
        color = settings.geneColor{3};
        for cond = 1:length(gainConditions)
            scatter(repmat(xVals(cond), size(dataNa, 1), 1), dataNa(:, cond), 40, settings.trialColor, '.');
            meanVal = mean(dataNa(:, cond));
            semVal = std(dataNa(:, cond)) / sqrt(size(dataNa, 1));
            errorbar(xVals(cond), meanVal, semVal, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'Color', color, 'CapSize', 10, 'LineWidth', 1);
        end

        % Set x-axis labels and limits
        axis padded
        ylim([150 500])
        xticks(gainConditions);
        xticklabels(arrayfun(@num2str, gainConditions, 'UniformOutput', false));
        xlabel('Gain');
        ylabel({'Median Interval'; 'Between Direction Changes (ms)'});

        % Add title for each row
        title(plotTitle);

        % Add p-value annotation to the plot (top right corner)
        text(0.95, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

        hold off;
    end

    % Save plot
    cd(folder.summary)
    plotname = 'dirchange_interval';
    saveas(gcf, [plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox, 'f');
    % Save vectorized plot
    cd(folder.vector)
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox, 'f');
end

%% Position of the object during direction changes
if exptFolder == 'AOTU019 KIR'
    % Initialize figure with a tiled layout
    figure;
    tiledlayout(2, 1);  % Two rows, one column
    set(gcf,'Position',[100 100 400 900]);

    % Call the general ANOVA function to analyze the contribution of genotype and gain to object position during direction changes
    [p_fixpos, ~] = run_genotype_anova_repeated(kirFixPos, wtFixPos, naFixPos, 'PositionFixation', folder);
    [p_notpos, ~] = run_genotype_anova_repeated(kirNotPos, wtNotPos, naNotPos, 'PositionOther', folder);

    % Format p-value text for display
    pval_text_fixation = ['p(gene) = ' num2str(p_fixpos(1)) ', p(gain) = ' num2str(p_fixpos(2))];
    pval_text_nonfixation = ['p(gene) = ' num2str(p_notpos(1)) ', p(gain) = ' num2str(p_notpos(2))];

    % Get the number of gain conditions
    gainConditions = settings.pursuitGain;

    % Loop for plotting fixation and non-fixation object positions
    for row = 1:2
        if row == 1
            % First row: Fixation object positions
            dataKir = kirFixPos;
            dataWt = wtFixPos;
            dataNa = naFixPos;
            plotTitle = 'Fixation';
            pval_text = pval_text_fixation;
        else
            % Second row: Non-fixation object positions
            dataKir = kirNotPos;
            dataWt = wtNotPos;
            dataNa = naNotPos;
            plotTitle = 'Non-Fixation';
            pval_text = pval_text_nonfixation;
        end

        % Select the current tile
        nexttile;
        hold on;

        % Plot KIR genotype
        xVals = gainConditions - 2;  % Offset for KIR
        color = settings.geneColor{1};
        for cond = 1:length(gainConditions)
            scatter(repmat(xVals(cond), size(dataKir, 1), 1), dataKir(:, cond), 40, settings.trialColor, '.');
            meanVal = mean(dataKir(:, cond));
            semVal = std(dataKir(:, cond)) / sqrt(size(dataKir, 1));
            errorbar(xVals(cond), meanVal, semVal, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'Color', color, 'CapSize', 10, 'LineWidth', 1);
        end

        % Plot WT genotype
        xVals = gainConditions;  % No offset for WT
        color = settings.geneColor{2};
        for cond = 1:length(gainConditions)
            scatter(repmat(xVals(cond), size(dataWt, 1), 1), dataWt(:, cond), 40, settings.trialColor, '.');
            meanVal = mean(dataWt(:, cond));
            semVal = std(dataWt(:, cond)) / sqrt(size(dataWt, 1));
            errorbar(xVals(cond), meanVal, semVal, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'Color', color, 'CapSize', 10, 'LineWidth', 1);
        end

        % Plot NA genotype
        xVals = gainConditions + 2;  % Offset for NA
        color = settings.geneColor{3};
        for cond = 1:length(gainConditions)
            scatter(repmat(xVals(cond), size(dataNa, 1), 1), dataNa(:, cond), 40, settings.trialColor, '.');
            meanVal = mean(dataNa(:, cond));
            semVal = std(dataNa(:, cond)) / sqrt(size(dataNa, 1));
            errorbar(xVals(cond), meanVal, semVal, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'Color', color, 'CapSize', 10, 'LineWidth', 1);
        end

        % Set x-axis labels and limits
        axis padded
        ylim([0 120])
        xticks(gainConditions);
        xticklabels(arrayfun(@num2str, gainConditions, 'UniformOutput', false));
        xlabel('Gain');
        ylabel({'Median Object Position'; 'During Direction Changes (deg)'});

        % Add title for each row
        title(plotTitle);

        % Add p-value annotation to the plot (top right corner)
        text(0.95, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

        hold off;
    end

    % Save plot
    cd(folder.summary)
    plotname = 'dirchange_objectpos';
    saveas(gcf, [plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox, 'f');
    % Save vectorized plot
    cd(folder.vector)
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox, 'f');
end

%% Plot Jump Corrections with binned means
% Colors (kir, wt, na)
cols = {settings.geneColor{1}, settings.geneColor{2}, settings.geneColor{3}};

% Flatten (animals × conditions cell arrays) → paired x/y vectors
[x_kir, y_kir] = flattenJump(kirJumpAll_pos, kirJumpAll_time);
[x_wt,  y_wt ] = flattenJump(wtJumpAll_pos,  wtJumpAll_time);
[x_na,  y_na ] = flattenJump(naJumpAll_pos,  naJumpAll_time);

X = {x_kir, x_wt, x_na};
Y = {y_kir, y_wt, y_na};
labels = {'Kir','WT','NA'};

figure; hold on;

% Scatter for each genotype
for g = 1:3
    scatter(X{g}, Y{g}, 20, 'Marker', '.', ...
        'MarkerEdgeColor', cols{g});
end

% Global 10° bins
x_all = vertcat(X{:});
edges = 0:10:ceil(max(x_all,[],'omitnan')/10)*10;
centers = edges(1:end-1) + diff(edges)/2;
nb = numel(edges)-1;

% Overlay binned mean ± SEM
for g = 1:3
    xg = X{g}; yg = Y{g};
    m = isfinite(xg) & isfinite(yg);
    xg = xg(m); yg = yg(m);

    [~,~,bin] = histcounts(xg, edges);
    mu = nan(1,nb); se = nan(1,nb);
    for b = 1:nb
        yy = yg(bin==b);
        if ~isempty(yy)
            mu(b) = mean(yy,'omitnan');
            if numel(yy) > 1
                se(b) = std(yy,0,'omitnan')/sqrt(numel(yy));
            else
                se(b) = 0;
            end
        end
    end

    keep = isfinite(mu);
    fill([centers(keep) fliplr(centers(keep))], ...
         [mu(keep)-se(keep) fliplr(mu(keep)+se(keep))], ...
         'k', 'FaceAlpha',settings.semAlpha  , 'EdgeColor','none');
    plot(centers(keep), mu(keep), 'Color', cols{g}, 'LineWidth', 2);
end

xlabel('Object displacement at jump (deg)');
ylabel('Correction time to zero-crossing (s)');
title('Jump size vs. correction time');
grid on; box on;
legend(labels, 'Location','best');


% Save plot
cd(folder.summary)
plotname = 'jump_correction_times';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');

%% Time of Jump Corrections for Small and Large Deviations
close all;

% Initialize figure with a tiled layout
figure;
tiledlayout(2, 1);  % Two rows, one column
set(gcf,'Position',[100 100 400 900]);

% Call the general ANOVA function to analyze the contribution of genotype and gain to small and large jump corrections
if exptFolder == 'AOTU019 KIR'
    [p_smalljump, ~] = run_genotype_anova_repeated(kirJumpSmall, wtJumpSmall, naJumpSmall, 'SmallJump', folder);
    [p_largejump, ~] = run_genotype_anova_repeated(kirJumpLarge, wtJumpLarge, naJumpLarge, 'LargeJump', folder);
    pval_text_smalljump = ['p(gene) = ' num2str(p_smalljump(1)) ', p(gain) = ' num2str(p_smalljump(2))];
    pval_text_largejump = ['p(gene) = ' num2str(p_largejump(1)) ', p(gain) = ' num2str(p_largejump(2))];
else
    [p_smalljump, ~] = run_genotype_anova1(kirJumpSmall, wtJumpSmall, naJumpSmall, 'SmallJump', folder);
    [p_largejump, ~] = run_genotype_anova1(kirJumpLarge, wtJumpLarge, naJumpLarge, 'LargeJump', folder);
    pval_text_smalljump = ['p(gene) = ' num2str(p_smalljump(1))];
    pval_text_largejump = ['p(gene) = ' num2str(p_largejump(1))];
end

% Get the number of gain conditions
gainConditions = settings.pursuitGain(1:nGain);

% Loop for plotting small and large jump corrections
for row = 1:2
    if row == 1
        % First row: Small deviations (jump corrections for small deviations)
        dataKir = kirJumpSmall.*100;
        dataWt = wtJumpSmall.*100;
        dataNa = naJumpSmall.*100;
        plotTitle = 'Small Jump Corrections';
        pval_text = pval_text_smalljump;
    else
        % Second row: Large deviations (jump corrections for large deviations)
        dataKir = kirJumpLarge.*100;
        dataWt = wtJumpLarge.*100;
        dataNa = naJumpLarge.*100;
        plotTitle = 'Large Jump Corrections';
        pval_text = pval_text_largejump;
    end

    % Select the current tile
    nexttile;
    hold on;

    % Plot KIR genotype
    xVals = gainConditions - 2;  % Offset for KIR
    color = settings.geneColor{1};
    for cond = 1:length(gainConditions)
        scatter(repmat(xVals(cond), size(dataKir, 1), 1), dataKir(:, cond), 40, settings.trialColor, '.');
        meanVal = mean(dataKir(:, cond), 'omitnan');
        semVal = std(dataKir(:, cond), 'omitnan') / sqrt(sum(~isnan(dataKir(:, cond))));
        errorbar(xVals(cond), meanVal, semVal, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'Color', color, 'CapSize', 10, 'LineWidth', 1);
    end

    % Plot WT genotype
    xVals = gainConditions;  % No offset for WT
    color = settings.geneColor{2};
    for cond = 1:length(gainConditions)
        scatter(repmat(xVals(cond), size(dataWt, 1), 1), dataWt(:, cond), 40, settings.trialColor, '.');
        meanVal = mean(dataWt(:, cond), 'omitnan');
        semVal = std(dataWt(:, cond), 'omitnan') / sqrt(sum(~isnan(dataWt(:, cond))));
        errorbar(xVals(cond), meanVal, semVal, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'Color', color, 'CapSize', 10, 'LineWidth', 1);
    end

    % Plot NA genotype
    xVals = gainConditions + 2;  % Offset for NA
    color = settings.geneColor{3};
    for cond = 1:length(gainConditions)
        scatter(repmat(xVals(cond), size(dataNa, 1), 1), dataNa(:, cond), 40, settings.trialColor, '.');
        meanVal = mean(dataNa(:, cond), 'omitnan');
        semVal = std(dataNa(:, cond), 'omitnan') / sqrt(sum(~isnan(dataNa(:, cond))));
        errorbar(xVals(cond), meanVal, semVal, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'Color', color, 'CapSize', 10, 'LineWidth', 1);
    end

    % Set x-axis labels and limits
    axis padded
    ylim([0 500])
    xticks(gainConditions);
    xticklabels(arrayfun(@num2str, gainConditions, 'UniformOutput', false));
    xlabel('Gain');
    ylabel({'Median Time to Correct Jumps'; 'for Small and Large Deviations (ms)'});

    % Add title for each row
    title(plotTitle);

    % Add p-value annotation to the plot (top right corner)
    text(0.95, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

    hold off;
end

% Save plot
cd(folder.summary)
plotname = 'jump_correction_binned';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');


%% Fly velocity distributions
disp('Analyzing histogram distributions...')
hLabel = 'Norm Count';
for c = 1:nGain
    thisGain = num2str(settings.pursuitGain(c));
    figure; set(gcf,'Position',[100 100 1800 800])
    tiledlayout(3,4,'TileSpacing','compact')
    fBin = fwdHist(:,1);
    aBin = angHist(:,1);
    sBin = sidHist(:,1);

    for v = 1:3
        switch v
            case 1 %forward
                hBin = fBin;
                kirVel = reshape(kirFwdHist(:,c,:),[],nKIR);
                wtVel = reshape(wtFwdHist(:,c,:),[],nWT);
                naVel = reshape(naFwdHist(:,c,:),[],nNA);
                xvar = 'Binned Forward Velocity (mm/s)';
            case 2 %angular
                hBin = aBin;
                kirVel = reshape(kirAngHist(:,c,:),[],nKIR);
                wtVel = reshape(wtAngHist(:,c,:),[],nWT);
                naVel = reshape(naAngHist(:,c,:),[],nNA);
                xvar = 'Binned Angular Velocity (deg/s)';
            case 3 %sideways
                hBin = sBin;
                kirVel = reshape(kirSidHist(:,c,:),[],nKIR);
                wtVel = reshape(wtSidHist(:,c,:),[],nWT);
                naVel = reshape(naSidHist(:,c,:),[],nNA);
                xvar = 'Binned Sideways Velocity (mm/s)';
        end
        % calculate mean and sem
        kirVel_mean = mean(kirVel,2);
        wtVel_mean = mean(wtVel,2);
        naVel_mean = mean(naVel,2);
        kirVel_sem = std(kirVel,[],2)./sqrt(nKIR);
        wtVel_sem = std(wtVel,[],2)./sqrt(nWT);
        naVel_sem = std(naVel,[],2)./sqrt(nNA);

        ax1 = nexttile; hold on
        plot(hBin,kirVel,'Color',settings.trialColor)
        plot(hBin,kirVel_mean,'Color',settings.geneColor{1},'LineWidth',1.5)
        xlabel(xvar); ylabel(hLabel); xline(0); axis tight
        if v==1
            title(settings.geneLabel{1})
        end
        ax2 = nexttile; hold on
        plot(hBin,wtVel,'Color',settings.trialColor)
        plot(hBin,wtVel_mean,'Color',settings.geneColor{2},'LineWidth',1.5)
        xlabel(xvar); ylabel(hLabel); xline(0); axis tight
        if v==1
            title(settings.geneLabel{2})
        end
        ax3 = nexttile; hold on
        plot(hBin,naVel,'Color',settings.trialColor)
        plot(hBin,naVel_mean,'Color',settings.geneColor{3},'LineWidth',1.5)
        xlabel(xvar); ylabel(hLabel); xline(0); axis tight
        if v==1
            title(settings.geneLabel{3})
        end
        linkaxes([ax1 ax2 ax3], 'xy')

        nexttile, hold on
        plot(hBin,kirVel_mean,'Color',settings.geneColor{1},'LineWidth',1.5)
        sp1 = patch([hBin; flipud(hBin)],[kirVel_mean-kirVel_sem; flipud(kirVel_mean+kirVel_sem)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
        sp1.FaceColor = settings.geneColor{1};
        plot(hBin,wtVel_mean,'Color',settings.geneColor{2},'LineWidth',1.5)
        sp2 = patch([hBin; flipud(hBin)],[wtVel_mean-wtVel_sem; flipud(wtVel_mean+wtVel_sem)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
        sp2.FaceColor = settings.geneColor{2};
        plot(hBin,naVel_mean,'Color',settings.geneColor{3},'LineWidth',1.5)
        sp3 = patch([hBin; flipud(hBin)],[naVel_mean-naVel_sem; flipud(naVel_mean+naVel_sem)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
        sp3.FaceColor = settings.geneColor{3};
        xlabel(xvar); ylabel(hLabel); xline(0); axis tight
    end
    % save plot
    sgtitle([thisGain 'X Directional Velocity Distribution'])
    cd(folder.summary)
    plotname = ['hist_' thisGain 'x_gain'];
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');
end
disp('Complete.')

%% Object distributions
disp('Analyzing object histogram distributions...')
hLabel = 'Norm Count';
pBin = posHist(:,1);
vBin = velHist(:,1);

posLim = [-100 100];
velLim = [-400 400];
hLimPos = [0 0.3];
hLimVel = [0 0.3];

% initialize, plot object position
figure; set(gcf,'Position',[100 100 1800 800])
tiledlayout(nGain,4,'TileSpacing','compact')

% for each condition
for c = 1:nGain
    % fetch data
    thisGain = num2str(settings.pursuitGain(c));
    thiskir = reshape(kirPosHist(:,c,:),[],nKIR);
    thiswt = reshape(wtPosHist(:,c,:),[],nWT);
    thisna = reshape(naPosHist(:,c,:),[],nNA);

    % calculate mean and sem
    meankir = mean(thiskir,2);
    semkir = std(thiskir,[],2)./sqrt(nKIR);
    meanwt = mean(thiswt,2);
    semwt = std(thiswt,[],2)./sqrt(nWT);
    meanna = mean(thisna,2);
    semna = std(thisna,[],2)./sqrt(nNA);

    % plot separately
    nexttile; hold on
    plot(pBin,thiskir,'Color',settings.trialColor)
    plot(pBin,meankir,'Color',settings.geneColor{1},'Linewidth',1.5)
    xlim(posLim); ylim(hLimPos); xline(0); ylabel([thisGain 'X'])
    if c==1
        title(settings.geneLabel{1})
    elseif c==nGain
        xlabel('Object Position (deg)')
    end
    nexttile; hold on
    plot(pBin,thiswt,'Color',settings.trialColor)
    plot(pBin,meanwt,'Color',settings.geneColor{2},'Linewidth',1.5)
    xlim(posLim); ylim(hLimPos); xline(0);
    if c==1
        title(settings.geneLabel{2})
    elseif c==nGain
        xlabel('Object Position (deg)')
    end
    nexttile; hold on
    plot(pBin,thisna,'Color',settings.trialColor)
    plot(pBin,meanna,'Color',settings.geneColor{3},'Linewidth',1.5)
    xlim(posLim); ylim(hLimPos); xline(0);
    if c==1
        title(settings.geneLabel{3})
    elseif c==nGain
        xlabel('Object Position (deg)')
    end

    % plot together
    nexttile; hold on
    plot(pBin,meankir,'Color',settings.geneColor{1},'Linewidth',1.5)
    sp1 = patch([pBin; flipud(pBin)],[meankir-semkir; flipud(meankir+semkir)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{1};
    plot(pBin,meanwt,'Color',settings.geneColor{2},'Linewidth',1.5)
    sp1 = patch([pBin; flipud(pBin)],[meanwt-semwt; flipud(meanwt+semwt)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{2};
    plot(pBin,meanna,'Color',settings.geneColor{3},'Linewidth',1.5)
    sp1 = patch([pBin; flipud(pBin)],[meanna-semna; flipud(meanna+semna)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{3};
    xlim(posLim); ylim(hLimPos); xline(0);
    if c==nGain
        xlabel('Object Position (deg)')
    end

end
% save plot
sgtitle('Object Position Distribution')
cd(folder.summary)
plotname = 'hist_objpos';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');


% initialize, plot object velocity
figure; set(gcf,'Position',[100 100 1800 800])
tiledlayout(nGain,4,'TileSpacing','compact')

% for each condition
for c = 1:nGain
    % fetch data
    thisGain = num2str(settings.pursuitGain(c));
    thiskir = reshape(kirVelHist(:,c,:),[],nKIR);
    thiswt = reshape(wtVelHist(:,c,:),[],nWT);
    thisna = reshape(naVelHist(:,c,:),[],nNA);

    % calculate mean and sem
    meankir = mean(thiskir,2);
    semkir = std(thiskir,[],2)./sqrt(nKIR);
    meanwt = mean(thiswt,2);
    semwt = std(thiswt,[],2)./sqrt(nWT);
    meanna = mean(thisna,2);
    semna = std(thisna,[],2)./sqrt(nNA);

    % plot separately
    nexttile; hold on
    plot(vBin,thiskir,'Color',settings.trialColor)
    plot(vBin,meankir,'Color',settings.geneColor{1},'Linewidth',1.5)
    xlim(velLim); ylim(hLimVel); xline(0); ylabel([thisGain 'X'])
    nexttile; hold on
    if c==1
        title(settings.geneLabel{1})
    elseif c==nGain
        xlabel('Object Velocity (deg/s)')
    end
    plot(vBin,thiswt,'Color',settings.trialColor)
    plot(vBin,meanwt,'Color',settings.geneColor{2},'Linewidth',1.5)
    xlim(velLim); ylim(hLimVel); xline(0);
    if c==1
        title(settings.geneLabel{2})
    elseif c==nGain
        xlabel('Object Velocity (deg/s)')
    end
    nexttile; hold on
    plot(vBin,thisna,'Color',settings.trialColor)
    plot(vBin,meanna,'Color',settings.geneColor{3},'Linewidth',1.5)
    xlim(velLim); ylim(hLimVel); xline(0);
    if c==1
        title(settings.geneLabel{3})
    elseif c==nGain
        xlabel('Object Velocity (deg/s)')
    end

    % plot together
    nexttile; hold on
    plot(vBin,meankir,'Color',settings.geneColor{1},'Linewidth',1.5)
    sp1 = patch([vBin; flipud(vBin)],[meankir-semkir; flipud(meankir+semkir)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{1};
    plot(vBin,meanwt,'Color',settings.geneColor{2},'Linewidth',1.5)
    sp1 = patch([vBin; flipud(vBin)],[meanwt-semwt; flipud(meanwt+semwt)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{2};
    plot(vBin,meanna,'Color',settings.geneColor{3},'Linewidth',1.5)
    sp1 = patch([vBin; flipud(vBin)],[meanna-semna; flipud(meanna+semna)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{3};
    xlim(velLim); ylim(hLimVel); xline(0);
    if c==nGain
        xlabel('Object Velocity (deg/s)')
    end

end
% save plot
sgtitle('Object Velocity Distribution')
cd(folder.summary)
plotname = 'hist_objvel';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

%% Relationship between error and turning
disp('Comparing relationship between error and forward velocity...')

% Initialize figure for separate plot
figure; set(gcf,'Position',[100 100 1200 800])
tiledlayout(3,nGain,'TileSpacing','compact')
fwdLim = [0 30]; fwdLimA = [0 15];

% Plot separate trials and means for each genotype
for g = 1:3
    switch g
        case 1, evt = kirEVFwd; thisN = nKIR;
        case 2, evt = wtEVFwd; thisN = nWT;
        case 3, evt = naEVFwd; thisN = nNA;
    end
    meanEVT = mean(evt,3,'omitnan');  % Calculate mean EVT
    % Plot each condition's EVT
    for c = 1:nGain
        nexttile; hold on
        thisTrialEVT = reshape(evt(:,c,:),[],thisN);  % Fetch data
        plot(posBins,thisTrialEVT,'Color',settings.trialColor,'LineWidth',settings.lwTri)
        plot(posBins,meanEVT(:,c),'Color',settings.geneColor{g},'LineWidth',settings.lwAvg)
        axis tight; xline(0); yline(0); xlim([-100 100]); ylim(fwdLim)
        if g==1, title([num2str(settings.pursuitGain(c)) 'X']); end
        if g==3, xlabel('Target Error (deg)'); end
        if c==1, ylabel('Forward Velocity (mm/s)'); end
    end
end

% Save plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_fwd_sep';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

% Initialize figure for combined plot
figure; set(gcf,'Position',[100 100 1700 600])
tiledlayout(1,nGain,'TileSpacing','compact')

% Plot combined data with SEM for each condition
for c = 1:nGain
    % Plot for each genotype (KIR, WT, NA)
    nexttile; hold on
    for g = 1:3
        switch g
            case 1
                evt = kirEVFwd;
                thisN = nKIR;
                evt_color = settings.geneColor{1};  % Color for KIR
            case 2
                evt = wtEVFwd;
                thisN = nWT;
                evt_color = settings.geneColor{2};  % Color for WT
            case 3
                evt = naEVFwd;
                thisN = nNA;
                evt_color = settings.geneColor{3};  % Color for NA
        end

        % Calculate mean and SEM
        meanEVT = mean(evt, 3, 'omitnan');
        semEVT = std(evt, 0, 3, 'omitnan') ./ sqrt(thisN);  % SEM calculation
        validE = sum(~isnan(evt(:, c, :)), 3) > thisN / 3;

        % Plot EVT with SEM
        plot(posBins(validE), meanEVT(validE, c), 'Color', evt_color, 'LineWidth', settings.lwAvg)
        % Plot SEM using patch
        sem_patch = patch([posBins(validE)'; flipud(posBins(validE)')], ...
            [meanEVT(validE, c) - semEVT(validE, c); flipud(meanEVT(validE, c) + semEVT(validE, c))], ...
            'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
        sem_patch.FaceColor = evt_color;

        axis tight; xline(0); yline(0); xlim([-80 80]); ylim(fwdLimA);
    end
    title([num2str(settings.pursuitGain(c)) 'X'])
    xlabel('Target Error (deg)')
    if c == 1
        ylabel('Forward Velocity (mm/s)')
    end
end

% Save combined plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_fwd_combined';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox,'f');

% save zoomed version
% Zoom into the plot for each condition
for c = 1:nGain
    nexttile(c); ylim([0 15]); xlim([-40 40])
end

% Save zoomed plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_fwd_combined_20x';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');
disp('Complete.')
%%
% Initialize figure for combined plot
figure; set(gcf,'Position',[100 100 300 600])
tiledlayout(1,1,'TileSpacing','compact')  % Adjusting to a single combined plot

% Combined plot with SEM for each genotype across all gain conditions
nexttile; hold on
for g = 1:3
    switch g
        case 1
            evt = kirEVFwd;
            thisN = nKIR;
            evt_color = settings.geneColor{1};  % Color for KIR
        case 2
            evt = wtEVFwd;
            thisN = nWT;
            evt_color = settings.geneColor{2};  % Color for WT
        case 3
            evt = naEVFwd;
            thisN = nNA;
            evt_color = settings.geneColor{3};  % Color for NA
    end

    % Calculate the mean across gain conditions per animal
    meanEVT_animal = mean(evt, 2, 'omitnan');  % Mean across gains for each animal
    meanEVT = mean(meanEVT_animal, 3, 'omitnan');  % Mean across animals
    semEVT = std(meanEVT_animal, 0, 3, 'omitnan') ./ sqrt(thisN);  % SEM calculation
    meanEVT(isnan(meanEVT)) = 0;
    semEVT(isnan(semEVT)) = 0;

    % Plot EVT with SEM
    plot(posBins, meanEVT, 'Color', evt_color, 'LineWidth', settings.lwAvg)

    % Plot SEM using patch
    sem_patch = patch([posBins'; flipud(posBins')], ...
        [meanEVT - semEVT; flipud(meanEVT + semEVT)], ...
        'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    sem_patch.FaceColor = evt_color;

    axis tight; xline(0); yline(0); xlim([-80 80]); ylim(fwdLimA);
end

% Labels and title
title('Combined Gain Condition Averages Across Genotypes')
xlabel('Target Error (deg)')
ylabel('Forward Velocity (mm/s)')

% Save combined plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_fwd_combined_across_gains';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox,'f');

% save zoomed version
% Zoom into the plot for each condition
xlim([-30 30])

% Save zoomed plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_fwd_across_gains_20x';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');
disp('Complete.')

%% Compare error vs turn for low vs high forward velocity (Control only, avg across gain)
disp('Comparing low vs high forward velocity in controls (averaged across gain)...')

% Combine WT and NA into a control group
ctrl_low = cat(3, wtEVTRL_lowfwd, naEVTRL_lowfwd);     % [pos x gain x fly]
ctrl_high = cat(3, wtEVTRL_highfwd, naEVTRL_highfwd);
nCTRL = size(ctrl_low, 3);

% Average across gain conditions
meanLow = mean(ctrl_low, 2, 'omitnan');   % [pos x 1 x fly]
meanHigh = mean(ctrl_high, 2, 'omitnan'); % [pos x 1 x fly]
meanLow = squeeze(meanLow);   % [pos x fly]
meanHigh = squeeze(meanHigh); % [pos x fly]

% Compute group mean and SEM
avgLow = mean(meanLow, 2, 'omitnan');
semLow = std(meanLow, 0, 2, 'omitnan') ./ sqrt(nCTRL);

avgHigh = mean(meanHigh, 2, 'omitnan');
semHigh = std(meanHigh, 0, 2, 'omitnan') ./ sqrt(nCTRL);

validH = sum(~isnan(meanLow), 2) > nCTRL / 3;
%validH = sum(~isnan(meanHigh), 2) > nCTRL / 3;

% Plotting
figure; set(gcf,'Position',[200 200 600 500])
hold on

% Plot LOW (black)
plot(posBins(validH), avgLow(validH), 'k', 'LineWidth', settings.lwAvg, 'LineStyle','-')
sem_patch_low = patch([posBins(validH)'; flipud(posBins(validH)')], ...
    [avgLow(validH)-semLow(validH); flipud(avgLow(validH)+semLow(validH))], ...
    'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
sem_patch_low.FaceColor = [0 0 0];  % black

% Plot HIGH (dashed color)
plot(posBins(validH), avgHigh(validH), 'Color', settings.geneColor{2}, ...
    'LineWidth', settings.lwAvg, 'LineStyle','--')
sem_patch_high = patch([posBins(validH)'; flipud(posBins(validH)')], ...
    [avgHigh(validH)-semHigh(validH); flipud(avgHigh(validH)+semHigh(validH))], ...
    'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
sem_patch_high.FaceColor = settings.geneColor{2};

xline(0); yline(0)
xlim([-60 60]); ylim([-200 200])
xlabel('Target Error (deg)')
ylabel('Angular Velocity (deg/s)')
title('Control: Low vs High Forward Velocity (Avg. Across Gain)')

% Save
cd(folder.summary)
saveas(gcf,'error_v_turn_control_avgAcrossGain.png');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf,'error_v_turn_control_avgAcrossGain.svg');

%% Estimate linear turning gain: slope fit for low and high forward velocity
disp('Fitting linear slopes to turning vs. error curves...')

% Initialize
slopesLow = nan(nCTRL, 1);
slopesHigh = nan(nCTRL, 1);

for f = 1:nCTRL
    % --- Low forward velocity ---
    yLow = meanLow(:,f); % angular velocity
    x = posBins(:);
    valid = ~isnan(yLow);
    if sum(valid) > 5
        x_valid = x(valid);
        y_valid = yLow(valid);
        slope = (x_valid \ y_valid); % fit without intercept
        slopesLow(f) = slope;
    end

    % --- High forward velocity ---
    yHigh = meanHigh(:,f); % angular velocity
    valid = ~isnan(yHigh);
    if sum(valid) > 5
        x_valid = x(valid);
        y_valid = yHigh(valid);
        slope = (x_valid \ y_valid); % fit without intercept
        slopesHigh(f) = slope;
    end
end

% Plotting
figure; set(gcf,'Position',[300 300 400 400]); hold on

% Plot paired lines between low and high slopes
for f = 1:nCTRL
    if ~isnan(slopesLow(f)) && ~isnan(slopesHigh(f))
        plot([1 2], [slopesLow(f) slopesHigh(f)], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5)
    end
end

% Scatter points (no jitter)
scatter(ones(nCTRL,1), slopesLow, 10, 'k', 'filled')
scatter(2*ones(nCTRL,1), slopesHigh, 10, 'k', 'filled')

% Median lines
plot([0.85 1.15], repmat(median(slopesLow,'omitnan'),1,2), 'k', 'LineWidth', 2)
plot([1.85 2.15], repmat(median(slopesHigh,'omitnan'),1,2), 'r', 'Color', settings.geneColor{2}, 'LineWidth', 2)

% Axes and labels
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Low Fwd','High Fwd'})
ylabel('Linear Slope (Turn Gain)')
title('Turn Gain by Forward Velocity')

% T-test
[~,pval] = ttest(slopesLow, slopesHigh);

% Display p-value in top right corner with full precision
yl = ylim;
text(2.45, yl(2), sprintf('p(fwd) = %.3e', pval), ...
    'HorizontalAlignment','right', 'VerticalAlignment','top', 'FontSize', 10)

% Save
cd(folder.summary)
saveas(gcf,'error_v_turn_control_gain_slope.png');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf,'error_v_turn_control_gain_slope.svg');


%% HD during slow vs fast walking (WT + NA only)
% Combine across genotypes (omit Kir)
all_slow = [wtHDslow, naHDslow];  % [bins x flies]
all_fast = [wtHDfast, naHDfast];

% Extract histogram values
bin_centers = hist_slow(:,1);          % center of bins in degrees
prob_slow = mean(all_slow, 2, 'omitnan');
prob_fast = mean(all_fast, 2, 'omitnan');

% Compute SEM
sem_slow = std(all_slow, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(all_slow), 2));
sem_fast = std(all_fast, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(all_fast), 2));

% Interpolation for smooth line plot
smooth_resolution = 0.5;
bin_interp = min(bin_centers):smooth_resolution:max(bin_centers);
mean_slow_smooth = interp1(bin_centers, prob_slow, bin_interp, 'pchip');
sem_slow_smooth  = interp1(bin_centers, sem_slow, bin_interp, 'pchip');
mean_fast_smooth = interp1(bin_centers, prob_fast, bin_interp, 'pchip');
sem_fast_smooth  = interp1(bin_centers, sem_fast, bin_interp, 'pchip');

% Define bin edges for polar histogram (from centers)
bin_width = median(diff(bin_centers));
bin_edges = [bin_centers' - bin_width/2, bin_centers(end) + bin_width/2]; % row vector
bin_edges_rad = deg2rad(bin_edges);

% Create new figure
figure;
tiledlayout(2,2, 'TileSpacing', 'compact');

% ==== Tile 1: Line plot ====
nexttile(1); hold on

% Slow (black)
patch([bin_interp, fliplr(bin_interp)], ...
    [mean_slow_smooth - sem_slow_smooth, fliplr(mean_slow_smooth + sem_slow_smooth)], ...
    'k', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
plot(bin_interp, mean_slow_smooth, 'k-', 'LineWidth', 1.5);

% Fast (red)
patch([bin_interp, fliplr(bin_interp)], ...
    [mean_fast_smooth - sem_fast_smooth, fliplr(mean_fast_smooth + sem_fast_smooth)], ...
    'r', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
plot(bin_interp, mean_fast_smooth, 'r-', 'LineWidth', 1.5);

xlabel('Panel Position (°)');
ylabel('Probability');
legend({'Slow ± SEM','Slow','Fast ± SEM','Fast'}, 'Location', 'best');
title('Panel Position Histogram (WT + NA)');
xlim([-180 180]); yline(0); box off;
xticks(-150:30:150)

% ==== Tile 2: Polar histogram - Slow ====
nexttile(3);
polarhistogram('BinEdges', bin_edges_rad, ...
    'BinCounts', prob_slow, ...
    'FaceColor', 'k', 'FaceAlpha', 0.6, ...
    'Normalization', 'probability');
rlim([0 0.25]);
title('Polar Histogram - Slow');

% ==== Tile 3: Polar histogram - Fast ====
nexttile(4);
polarhistogram('BinEdges', bin_edges_rad, ...
    'BinCounts', prob_fast, ...
    'FaceColor', 'r', 'FaceAlpha', 0.6, ...
    'Normalization', 'probability');
rlim([0 0.25]);
title('Polar Histogram - Fast');

% Save the figure
cd(folder.summary)
plotname = 'hd_slowfast_wtna';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

% Save vectorized version
cd(folder.vector)
set(gcf, 'renderer', 'Painters');
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.');

%% end
disp('ALL ANALYSES COMPLETE.')
end

