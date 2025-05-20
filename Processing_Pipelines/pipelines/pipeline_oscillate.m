% pipeline_oscillate
%
% Pipeline Function
% Pulls all processed files from ALL flies in a given oscillatory experiment
% and performs necessary analyses and plots accordingly. Can be used for both
% behavior-only experiments and ephys experiments.
%
% INPUTS
% exptFolder - overarching experiment folder
% trialTypes - list of trial conditions (e.g., different oscillation speeds)
%
% The function pools data across all flies, performs analysis on spike rate,
% target motion, cross-correlations, pursuit performance, and tracking indices.
% Results are saved as plots, including raster, summary figures, and comparisons
% across trial types.
%
% 06/21/2022 - MC adapted from visual pursuit pipeline
% 07/08/2024 - MC cleaned and simplified
% 07/16/2024 - MC added pursuit index
%
function pipeline_oscillate(exptFolder,trialTypes)
%% initialize
disp('STARTING ANALYSES FOR POOLED PURSUIT EXPERIMENT...')
close all

% load in processing settings
settings = processSettings();

% set filename info and create necessary directories
filebase = strrep(exptFolder,' ','_');
% generate folder structures as needed
folder = generateFolders(exptFolder);

% number of target speeds used by this battery
nSpeeds = length(trialTypes);

%% set plotting variables

% plot values
sr_limit = [0 150];
sr_limitA = [0 90];
sr_limitD2 = [-70 70];
sr_limitD = [0 50];

ps_limit = [-35 35];

fwd_limit = [0 12];
ang_limit = [-250 250];
sid_limit = [-7 7];
fwd_limit2 = [0 8];
ang_limit2 = [-200 200];
sid_limit2 = [-2.5 2.5];
fwd_limita = [-100 100];
ang_limita = [-800 800];
sid_limita = [-100 100];

xc_lim = [-400 400];
xc_range = [0 1];

% set common variables
cm = [0.7 0.7 0.7 ; flip(colormap(cool(nSpeeds-1)))]; %color map
close all


%% load in and pool all trials from each experiment folder
disp('Loading in datasets...')
% pull file info
cd(folder.int)
allFiles = dir('*int.mat');
nFlies = length(allFiles);

% initialize data storage arrays
nFliesThresh = nFlies; %number of flies used in behavior analysis
nt_t = 0; %counter
allForward = [];
allSideway = [];
allAngular = [];
allPanelPs = [];
allSpikeRt = [];
normSweep = [];

for nt = 1:nFlies
    % load this trial
    disp(['Loading fly ' num2str(nt) '/' num2str(nFlies)])
    thisTrial = allFiles(nt).name;
    thisFly = thisTrial(6:16);
    cd(folder.int); load(thisTrial)

    % pool this trial data
    allForward = [allForward,int_forward];
    allSideway = [allSideway,int_sideway];
    allAngular = [allAngular,int_angular];
    allPanelPs = [allPanelPs,int_panelps];
    allSpikeRt = [allSpikeRt,int_spikert];

    % determine if this fly exhibitted sufficient walking behavior
    flyRunTime(nt,1) = (sum(int_forward>settings.runThreshE,'all')/length(int_time))*60;

    % for each trial type
    for s = 1:nSpeeds
        % pull this fly data for ONLY this trial type
        thisForward = int_forward(:,:,s);
        thisSideway = int_sideway(:,:,s);
        thisAngular = int_angular(:,:,s);
        thisPanelPs = int_panelps(:,:,s);
        thisSpikeRt = int_spikert(:,:,s);

        % analyze visual tuning for this trial type
        [thisSweep,thisMean_all] = osc_v_output(thisPanelPs,thisForward,thisSpikeRt,-1); %all
        [~,thisMean_high] = osc_v_output(thisPanelPs,thisForward,thisSpikeRt,settings.runThreshE); %run only
        [~,thisMean_rest] = osc_v_output(thisPanelPs,thisForward,thisSpikeRt,0); %rest only
        sweepPos{s} = thisSweep; %store sweeps
        % fetch size of arrays
        checkSweep(1) = size(thisMean_all,1);
        checkSweep(2) = size(thisMean_high,1);
        checkSweep(3) = size(thisMean_rest,1);
        if nt==1 % ensure all arrays are the same length
            normSweep(s) = checkSweep(1);
        end
        thisMean_all = thisMean_all(1:normSweep(s),:);
        thisMean_high = thisMean_high(1:normSweep(s),:);
        thisMean_rest = thisMean_rest(1:normSweep(s),:);
        % store visual tuning analysis from this trial type
        srvpos_all{nt,s} = thisMean_all;
        srvpos_high{nt,s} = thisMean_high;
        srvpos_rest{nt,s} = thisMean_rest;

        % only include flies that showed sufficient running in behavior analyses
        if flyRunTime(nt,1)>settings.minRunTime
            if s==1
                nt_t = nt_t+1; %update counter
            end
            % analyze pursuit performance for this trial type
            [fidelity,vigor,w_sr,w_time] = pursuit_performance(thisPanelPs,thisAngular,thisSpikeRt,int_time);
            allFidelity{nt_t,s} = fidelity;
            allVigor{nt_t,s} = vigor;
            allWinSR{nt_t,s} = w_sr;

            % analyze directional velocity for this trial type
            [~,fwdMean] = osc_v_output(thisPanelPs,thisForward,thisForward,settings.runThreshE); %run only
            [~,angMean] = osc_v_output(thisPanelPs,thisForward,thisAngular,settings.runThreshE); %run only
            [~,sidMean] = osc_v_output(thisPanelPs,thisForward,thisSideway,settings.runThreshE); %run only
            oscDur = size(angMean,1);
            oscDur2 = round(oscDur/2);
            % adjust for individual fly biases by centering
            ang_bias = 0; (max(angMean)+min(angMean))/2;
            sid_bias = 0; (max(sidMean)+min(sidMean))/2;
            % ensure all arrays are the same length
            fwdMean = fwdMean(1:normSweep(s),:);
            angMean = angMean(1:normSweep(s),:);
            sidMean = sidMean(1:normSweep(s),:);
            % store directional velocity analysis for this trial type
            posvfwd{nt_t,s} = fwdMean;
            posvang{nt_t,s} = angMean-ang_bias;
            posvsid{nt_t,s} = sidMean-sid_bias;
            % store R+L directional velocity analysis for this trial type
            posvfwdRL{nt_t,s} = mean([fwdMean circshift(fwdMean,oscDur2)],2,'omitnan');
            posvangRL{nt_t,s} = mean([angMean -circshift(angMean,oscDur2)],2,'omitnan');
            posvsidRL{nt_t,s} = mean([sidMean -circshift(sidMean,oscDur2)],2,'omitnan');

            % analyze behavior tuning for this trial type
            [binVel] = spikert_binvelocity(thisForward,thisAngular,thisSideway,thisSpikeRt,int_time,0); % velocity w/o lag
            [binVelL] = spikert_binvelocity(thisForward,thisAngular,thisSideway,thisSpikeRt,int_time,1); %velocity w/ lag
            [binAcc] = spikert_binacceleration(thisForward,thisAngular,thisSideway,thisSpikeRt,int_time,1); %acceleration
            % store behavior tuning analysis from this trial type
            srvfwd_vel{nt_t,s} = binVel.fwdMean';
            srvang_vel{nt_t,s} = binVel.angMean';
            srvsid_vel{nt_t,s} = binVel.sidMean';
            srvfwd_velL{nt_t,s} = binVelL.fwdMean';
            srvang_velL{nt_t,s} = binVelL.angMean';
            srvsid_velL{nt_t,s} = binVelL.sidMean';
            srvfwd_acc{nt_t,s} = binAcc.fwdMean';
            srvang_acc{nt_t,s} = binAcc.angMean';
            srvsid_acc{nt_t,s} = binAcc.sidMean';

            % run xcorr analysis for this trial type
            cd(folder.xcorr)
            thisXCorrFile = [thisFly '_' num2str(s) '_xc.mat'];
            if exist(thisXCorrFile,'file')
                % load in previously run xcorr analysis
                load(thisXCorrFile)
            else
                % run modified xcorr analysis
                [r_val,lag_t] = spikert_xcorr(thisSpikeRt,thisForward,thisAngular,thisSideway,int_time);
                % save
                save(thisXCorrFile,'r_val','lag_t','-v7.3');
            end
            % find xcorr peaks for each result
            [~,f_locs] = findpeaks(r_val.fwd,'MinPeakProminence',settings.minXCorrProm,'SortStr','descend');
            [~,a_locs] = findpeaks(r_val.ang,'MinPeakProminence',settings.minXCorrProm,'SortStr','descend');
            [~,s_locs] = findpeaks(r_val.sid,'MinPeakProminence',settings.minXCorrProm,'SortStr','descend');
            % store xcorr analysis and peaks (largest)
            % omit trials in which the xcorr peak too small or absent
            if ~isempty(f_locs)
                r_val_fwd(:,nt_t,s) = r_val.fwd;
                peak_fwd(nt_t,s) = lag_t(f_locs(1));
            else
                r_val_fwd(:,nt_t,s) = nan(length(r_val.fwd),1);
                peak_fwd(nt_t,s) = nan;
            end
            if ~isempty(a_locs)
                r_val_ang(:,nt_t,s) = r_val.ang;
                peak_ang(nt_t,s) = lag_t(a_locs(1));
            else
                r_val_ang(:,nt_t,s) = nan(length(r_val.ang),1);
                peak_ang(nt_t,s) = nan;
            end
            if ~isempty(s_locs)
                r_val_sid(:,nt_t,s) = r_val.sid;
                peak_sid(nt_t,s) = lag_t(s_locs(1));
            else
                r_val_sid(:,nt_t,s) = nan(length(r_val.sid),1);
                peak_sid(nt_t,s) = nan;
            end
            % run autocorr analysis for this trial type
            % thisAutoCorrFile = [thisFly '_' num2str(s) '_ac.mat'];
            % if exist(thisAutoCorrFile,'file')
            %     % load in previously run xcorr analysis
            %     load(thisAutoCorrFile)
            % else
            %     % run modified xcorr analysis
            %     [r_val,lag_t] = auto_xcorr(thisForward,thisAngular,thisSideway,int_time);
            %     % save
            %     save(thisAutoCorrFile,'r_val','lag_t','-v7.3');
            % end
            % % find xcorr peaks for each result
            % [~,f_locs] = findpeaks(r_val.fwd,'MinPeakProminence',settings.minXCorrProm,'SortStr','descend');
            % [~,a_locs] = findpeaks(r_val.ang,'MinPeakProminence',settings.minXCorrProm,'SortStr','descend');
            % [~,s_locs] = findpeaks(r_val.sid,'MinPeakProminence',settings.minXCorrProm,'SortStr','descend');
            % % store xcorr analysis and peaks (largest)
            % % omit trials in which the xcorr peak too small or absent
            % if ~isempty(f_locs)
            %     ar_val_fwd(:,nt_t,s) = r_val.fwd;
            % else
            %     ar_val_fwd(:,nt_t,s) = nan(length(r_val.fwd),1);
            % end
            % if ~isempty(a_locs)
            %     ar_val_ang(:,nt_t,s) = r_val.ang;
            % else
            %     ar_val_ang(:,nt_t,s) = nan(length(r_val.ang),1);
            % end
            % if ~isempty(s_locs)
            %     ar_val_sid(:,nt_t,s) = r_val.sid;
            % else
            %     ar_val_sid(:,nt_t,s) = nan(length(r_val.sid),1);
            % end
        else
            if s==1
                disp([thisFly ' omitted from behavior analyses.'])
                nFliesThresh = nFliesThresh-1;
            end
        end
    end
end
% store velocity bins
binFwd_v = binVel.fwdBin;
binAng_v = binVel.angBin;
binSid_v = binVel.sidBin;
% store acceleration bins
binFwd_a = binAcc.fwdBin;
binAng_a = binAcc.angBin;
binSid_a = binAcc.sidBin;

%% by trial: plot spike rate vs target motion
disp('By trial type: analyzing spike rate vs target motion...')

% for each trial type
for s=1:nSpeeds
    % initialize
    figure; set(gcf,'Position',[100 100 1000 500])
    tiledlayout(1,3,"TileSpacing","compact")
    t = int_time(1:length(sweepPos{s}));
    cT = mean(t);
    sweepLengths(s) = length(t);

    % for each run
    for r = 1:3
        switch r
            case 1 %all
                thisTrialset = cat(2,srvpos_all{:,s});
            case 2 %low
                thisTrialset = cat(2,srvpos_rest{:,s});
            case 3 %high
                thisTrialset = cat(2,srvpos_high{:,s});
        end
        thisSweepMean = mean(thisTrialset,2,'omitnan');

        % plot trials and mean for this trial type
        nexttile
        % bar position
        fs = gca; yyaxis right
        plot(t,sweepPos{s},'k','LineWidth',2)
        fs.YAxis(1).Color = 'k'; fs.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k'); xline(cT,'Color','k'); cT = mean(t);
        % firing rate
        yyaxis left
        plot(t,thisTrialset,'Color', settings.trialColor,'LineStyle','-','LineWidth',settings.lwTri,'Marker','none');hold on
        plot(t,thisSweepMean,'Color', settings.spkColor,'LineStyle','-','LineWidth',settings.lwAvg,'Marker','none'); hold on
        % adjust axes
        axis tight; ylim(sr_limit); ylabel([settings.behaviorGroup{r} ' ' settings.spkLabel]); xlabel('Time (s)')
        hold off

        %store group means for further analysis
        switch r
            case 1
                means_srvpos_all{s} = thisSweepMean;
                sem_srvpos_all{s} = std(thisTrialset,0,2,'omitnan')./sqrt(nFlies);
            case 2
                means_srvpos_rest{s} = thisSweepMean;
                sem_srvpos_rest{s} = std(thisTrialset,0,2,'omitnan')./sqrt(nFlies);
            case 3
                means_srvpos_high{s} = thisSweepMean;
                sem_srvpos_high{s} = std(thisTrialset,0,2,'omitnan')./sqrt(nFlies);
        end
    end

    thisletter = settings.letters(s);
    sgtitle([strrep(filebase,'_',' ') ' ' trialTypes{s}])
    % save plot
    cd(folder.summary)
    plotname = ['srvpos_' thisletter];
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

end
disp('Complete.')


%% comparison: plot spike rate vs target motion
disp('Comparison: analyzing spike rate vs target motion...')
[maxLength,maxIdx] = max(sweepLengths); % reshape everything to longest sweep size
x = 1:maxLength;
% initialize storage variables
group_trophs = [];
group_peaksr = [];
group_peakIdx = [];
group_peakpos = [];

% for all movement conditions
for r = 1
    %pull run info
    thisName = settings.behaviorGroup{r};
    % initialize
    figure; set(gcf,'Position',[100 100 1000 500])

    % for each trial type
    for s = 1:nSpeeds
        switch r
            case 1 %all
                thisGroup = cat(2,srvpos_all{:,s});
                thisMean = means_srvpos_all{s};
                thisSEM = sem_srvpos_all{s};
            case 2 %quiescent
                thisGroup = cat(2,srvpos_rest{:,s});
                thisMean = means_srvpos_rest{s};
                thisSEM = sem_srvpos_rest{s};
            case 3 %high
                thisGroup = cat(2,srvpos_high{:,s});
                thisMean = means_srvpos_high{s};
                thisSEM = sem_srvpos_high{s};
        end
        % pull this sweep
        thisSweep = sweepPos{s};
        % pull indices when target was in ipsilateral hemifield ONLY
        ipsiIdx = find(thisSweep > 0);
        thisIpsiSweep = thisSweep(ipsiIdx);

        % resize
        resizeMean = imresize(thisMean,[maxLength 1], 'nearest');
        resizeSEM = imresize(thisSEM,[maxLength 1], 'nearest');

        % pull troph, peak rate, peak position
        % first pull trophs
        group_trophs(:,s,r) = min(thisGroup)';
        trophMean = mean(group_trophs(:,s,r),'omitnan');
        trophSEM = std(group_trophs(:,s,r),1,'omitnan')/sqrt(nFlies);
        % second pull peaks
        [thispeak,thisidx] = max(thisGroup(ipsiIdx,:));
        group_peaksr(:,s,r) = thispeak';
        group_peakIdx(:,s,r) = thisidx';
        peaksrMean = mean(group_peaksr(:,s,r),'omitnan');
        peaksrSEM = std(group_peaksr(:,s,r),1,'omitnan')/sqrt(nFlies);
        % third pull corresponding peak position
        group_peakpos(:,s,r) = thisIpsiSweep(group_peakIdx(:,s,r));
        peakposMean = mean(group_peakpos(:,s,r),'omitnan');
        peakposSEM = std(group_peakpos(:,s,r),1,'omitnan')/sqrt(nFlies);

        % plot mean firing rates for each condition
        subplot(1,5,[1 2])
        % if first run, plot bar position below
        if s ==1
            fc = gca;
            % bar position
            yyaxis right
            plot(x,sweepPos{maxIdx},'k','LineWidth',2)
            axis tight
            fc.YAxis(1).Color = 'k';
            fc.YAxis(2).Visible = 'off';
            % add reference lines
            yline(0,'Color','k')
            cT = mean(x);
            xline(cT,'Color','k')
        end
        yyaxis left; hold on
        se(1) = patch([x'; flipud(x')],[(resizeMean-resizeSEM); flipud((resizeMean+resizeSEM))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
        plot(x,resizeMean,'-','LineWidth',2,'Color',cm(s,:))
        axis tight; ylim(sr_limitA); ylabel(['Average Pursuit ' settings.spkLabel]); xticklabels([])

        % plot troph firing rate for each condition
        % add a bit of jitter such that points are not overlapping
        jitter = 0.1;
        tt_j = s + (-jitter+(jitter*2)*rand(1,nFlies));
        subplot(1,5,3); hold on
        plot(tt_j,group_trophs(:,s,r),'.','Color', [0.8 0.8 0.8])
        errorbar(s,trophMean,trophSEM,'o','Color',cm(s,:))
        xlim([0 nSpeeds+1]); ylim(sr_limit)
        xticks(1:nSpeeds); xticklabels(trialTypes); ylabel(['Min ' settings.spkLabel])

        % plot peak firing rate for each condition
        subplot(1,5,4); hold on
        plot(tt_j,group_peaksr(:,s,r),'.','Color', [0.8 0.8 0.8])
        errorbar(s,peaksrMean,peaksrSEM,'o','Color',cm(s,:))
        xlim([0 nSpeeds+1]); ylim(sr_limit)
        xticks(1:nSpeeds); xticklabels(trialTypes); ylabel(['Peak ' settings.spkLabel])

        % plot position peak for each condition
        subplot(1,5,5); hold on
        plot(tt_j,group_peakpos(:,s,r),'.','Color', [0.8 0.8 0.8])
        errorbar(s,peakposMean,peakposSEM,'o','Color',cm(s,:))
        yline(0); xlim([0 nSpeeds+1]); ylim(ps_limit)
        xticks(1:nSpeeds); xticklabels(trialTypes); ylabel('Peak Firing Position (deg)')

    end
    sgtitle([strrep(filebase,'_',' ') ' ' thisName])
    % save plot
    cd(folder.summary)
    plotname = ['srvpos_across_' thisName];
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

    disp('Complete.')

    % Determine if difference is sig
    fetchP1peaks = group_peaksr(:,2:end);
    [num_flies, num_speeds] = size(fetchP1peaks);

    % Create a factor for object speed
    object_speed = repmat(1:num_speeds, num_flies, 1); % Speed indices (1 to num_speeds)
    object_speed = object_speed(:); % Column vector

    % Create a factor for fly
    fly = repmat((1:num_flies)', 1, num_speeds); % Fly indices (1 to num_flies)
    fly = fly(:); % Column vector

    % Reshape the peak firing rate data to a column
    peak_firing_rates = fetchP1peaks(:);

    % Specify folder to save results
    output_folder = folder.summary;
    if ~isfolder(output_folder)
        mkdir(output_folder); % Create the folder if it doesn't exist
    end

    object_speed = categorical(object_speed);
    fly = categorical(fly);

    % Define string labels for trial types (object speeds)
    subsetTrialTypes = {'Speed15', 'Speed25', 'Speed35', 'Speed55', 'Speed75'};

    % Get dimensions
    [numFlies, numSpeeds] = size(fetchP1peaks);

    % Reshape data for long format
    response = fetchP1peaks(:);  % Column of peak firing rates

    % Create factors
    speed = repmat(subsetTrialTypes, numFlies, 1);
    speed = categorical(speed(:));  % Object speed as categorical

    fly = repelem((1:numFlies)', numSpeeds);
    fly = categorical(fly);  % Each fly ID, repeated for each speed

    % Build table for LME
    T_lme = table(response, speed, fly, ...
        'VariableNames', {'PeakFR', 'Speed', 'Fly'});

    % Fit linear mixed-effects model
    lme = fitlme(T_lme, 'PeakFR ~ Speed + (1|Fly)');

    % Show ANOVA table for fixed effects
    anovaTbl = anova(lme);
    disp(anovaTbl);

    % Assuming group_peaksr is a matrix where rows = flies and columns = speeds
    [num_flies, num_speeds] = size(fetchP1peaks);

    % Define x-axis values (object speed conditions)
    object_speeds = 1:num_speeds;

    % Create the figure
    figure;
    hold on;

    % Plot each fly as a grey line
    for f = 1:num_flies
        plot(object_speeds, fetchP1peaks(f, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    end

    % Calculate the median firing rate for each speed
    median_rates = median(fetchP1peaks, 1);

    % Plot the median as black underscore markers
    plot(object_speeds, median_rates, '_', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2);

    % Customize plot
    xlabel('Object Speed Conditions');
    ylabel('Peak Firing Rate');
    title('Peak Firing Rate Across Object Speeds');
    set(gca, 'XTick', object_speeds); % Ensure x-axis matches the speed indices
    grid on;
    box on;

    % Save the figure (optional)
    % saveas(gcf, fullfile(folder.summary, 'peak_firing_rate_plot.png'));

    hold off;
    sgtitle([strrep(filebase,'_',' ') ' ' thisName])
    % save plot
    cd(folder.summary)
    plotname = ['srvpos_across2_' thisName];
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

end


%% diference: plot R-L spike rate vs target motion
disp('Difference: analyzing expected R-L spike rate vs target motion...')
[maxLength,maxIdx] = max(sweepLengths); % reshape everything to longest sweep size
x = 1:maxLength;

% for all movement conditions
for r = 1:3
    %pull run info
    thisName = settings.behaviorGroup{r};

    % initialize
    figure; set(gcf,'Position',[100 100 400 600])
    subx = 1; suby = 3;

    % for each trial type
    for s = 1:nSpeeds
        switch r
            case 1 %all
                thisSRSet = cat(2,srvpos_all{:,s});
            case 2 %quiescent
                thisSRSet = cat(2,srvpos_rest{:,s});
            case 3 %high
                thisSRSet = cat(2,srvpos_high{:,s});
        end
        % calculate group mean for each data subset
        thisIpsiAll = thisSRSet;
        thisCntrAll = circshift(thisIpsiAll,floor(length(thisIpsiAll)/2));
        % pull mean and resize
        meanIpsiOnly = imresize(mean(thisIpsiAll,2,'omitnan'),[maxLength 1], 'nearest');
        meanCntrOnly = imresize(mean(thisCntrAll,2,'omitnan')*-1,[maxLength 1], 'nearest');
        thisDiffMean = imresize(mean(thisIpsiAll-thisCntrAll,2,'omitnan'),[maxLength 1], 'nearest');
        thisDiffSEM = imresize(std(thisIpsiAll-thisCntrAll,0,2,'omitnan')./sqrt(nFlies),[maxLength 1], 'nearest');
        thisSweep = sweepPos{s};

        % plot ipsi sr only
        subplot(suby,subx,1)
        plot(x,meanIpsiOnly,'Color',cm(s,:)); hold on
        plot(x,meanCntrOnly,'Color',cm(s,:)); hold on
        ylim([-sr_limitA(2) sr_limitA(2)])
        axis tight
        yline(0)
        ylabel('R & L')
        xticklabels([])

        % plot expected R-L
        subplot(suby,subx,[2 3])
        % if first run, plot bar position below
        if s ==1
            fc = gca;
            % bar position
            yyaxis right
            plot(x,sweepPos{maxIdx},'k','LineWidth',2)
            axis tight
            fc.YAxis(1).Color = 'k';
            fc.YAxis(2).Visible = 'off';
            % add reference lines
            yline(0,'Color','k')
            cT = mean(x);
            xline(cT,'Color','k')
        end
        yyaxis left
        se(1) = patch([x'; flipud(x')],[(thisDiffMean-thisDiffSEM); flipud((thisDiffMean+thisDiffSEM))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
        hold on
        plot(x,thisDiffMean,'-','LineWidth',2,'Color',cm(s,:)); hold on
        axis tight
        ylim(sr_limitD2)
        ylabel(['R-L ' settings.spkLabel])
        xticklabels([])

    end

    sgtitle([strrep(filebase,'_',' ') ' ' thisName])
    % save plot
    cd(folder.summary)
    plotname = ['srvpos_difference_' thisName];
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');
    disp('Complete.')
end


%% by trial: plot pursuit behavior vs target motion
disp('By trial type: analyzing pursuit behavior across conditions...')

% for each trial type
for s = 1:nSpeeds
    % fetch data
    groupFwd = cat(2,posvfwd{:,s});
    groupAng = cat(2,posvang{:,s});
    groupSid = cat(2,posvsid{:,s});
    nFliesActual = sum(~isnan(groupFwd(1,:)));
    % calculate mean for this trial type
    mean_pursuit_fwd{s} = mean(groupFwd,2,'omitnan');
    mean_pursuit_ang{s} = mean(groupAng,2,'omitnan');
    mean_pursuit_sid{s} = mean(groupSid,2,'omitnan');
    sweepLength(s) = length(mean_pursuit_fwd{s});
    % calculate SEM for pursuit of this trial type
    sem_pursuit_fwd{s} = std(groupFwd,0,2,'omitnan')/sqrt(nFliesActual);
    sem_pursuit_ang{s} = std(groupAng,0,2,'omitnan')/sqrt(nFliesActual);
    sem_pursuit_sid{s} = std(groupSid,0,2,'omitnan')/sqrt(nFliesActual);

    % initialize
    figure; set(gcf,'Position',[100 100 1500 500])
    t = int_time(1:length(sweepPos{s}));

    % plot ALL data for this trial type
    % forward
    subplot(1,3,1)
    fs = gca;
    % bar position
    yyaxis right
    plot(t,sweepPos{s},'k','LineWidth',2)
    fs.YAxis(1).Color = 'k';
    fs.YAxis(2).Visible = 'off';
    % add reference lines
    yline(0,'Color','k')
    cT = mean(t);
    xline(cT,'Color','k')
    % firing rate
    yyaxis left
    plot(t,cat(2,posvfwd{:,s}),'Color', settings.trialColor,'LineStyle','-','LineWidth',settings.lwTri,'Marker','none');hold on
    plot(t,mean_pursuit_fwd{s},'Color', settings.velColor{1},'LineStyle','-','LineWidth',settings.lwAvg,'Marker','none'); hold on
    % adjust axes
    axis tight
    ylim(fwd_limit)
    ylabel(settings.velLabel{1})
    hold off

    % angular
    subplot(1,3,2)
    fs = gca;
    % bar position
    yyaxis right
    plot(t,sweepPos{s},'k','LineWidth',2)
    fs.YAxis(1).Color = 'k';
    fs.YAxis(2).Visible = 'off';
    % add reference lines
    yline(0,'Color','k')
    cT = mean(t);
    xline(cT,'Color','k')
    % firing rate
    yyaxis left
    plot(t,cat(2,posvang{:,s}),'Color', settings.trialColor,'LineStyle','-','LineWidth',settings.lwTri,'Marker','none');hold on
    plot(t,mean_pursuit_ang{s},'Color', settings.velColor{2},'LineStyle','-','LineWidth',settings.lwAvg,'Marker','none'); hold on
    % adjust axes
    axis tight
    ylim(ang_limit)
    ylabel(settings.velLabel{2})
    hold off

    % sideways
    subplot(1,3,3)
    fs = gca;
    % bar position
    yyaxis right
    plot(t,sweepPos{s},'k','LineWidth',2)
    fs.YAxis(1).Color = 'k';
    fs.YAxis(2).Color = 'k';
    ylabel('Object Position (deg)')
    % add reference lines
    yline(0,'Color','k')
    cT = mean(t);
    xline(cT,'Color','k')
    % firing rate
    yyaxis left
    plot(t,cat(2,posvsid{:,s}),'Color', settings.trialColor,'LineStyle','-','LineWidth',settings.lwTri,'Marker','none');hold on
    plot(t,mean_pursuit_sid{s},'Color', settings.velColor{3},'LineStyle','-','LineWidth',settings.lwAvg,'Marker','none'); hold on
    % adjust axes
    axis tight
    ylim(sid_limit)
    ylabel(settings.velLabel{3})
    hold off

    thisletter = settings.letters(s);
    sgtitle([strrep(filebase,'_',' ') ' ' trialTypes{s} ' (n = ' num2str(nFliesActual) '/' num2str(nFlies) ')'])
    % save plot
    cd(folder.summary)
    plotname = ['velvpos_' thisletter];
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');

end
disp('Complete.')

%% comparison: plot pursuit behavior vs target motion
disp('Comparison: analyzing pursuit behavior across conditions...')

% initialize
figure; set(gcf,'Position',[100 100 1500 600])
[maxLength,maxIdx] = max(sweepLengths); % reshape everything to longest sweep size
x = 1:maxLength;
xval = 1:nSpeeds;

for s=1:nSpeeds
    % pull individual trials for each directional velocity
    groupFwd = cell2mat(posvfwd(:,s)');
    groupAng = cell2mat(posvang(:,s)');
    groupSid = cell2mat(posvsid(:,s)');
    nFliesActual = sum(~isnan(groupFwd(1,:)));
    % calculate mean and sem
    meanFwd = mean(groupFwd,2,'omitnan');
    meanAng = mean(groupAng,2,'omitnan');
    meanSid = mean(groupSid,2,'omitnan');
    semFwd = std(groupFwd,0,2,'omitnan')./sqrt(nFliesActual);
    semAng = std(groupAng,0,2,'omitnan')./sqrt(nFliesActual);
    semSid = std(groupSid,0,2,'omitnan')./sqrt(nFliesActual);

    % resize to plot together
    meanFwdrs = imresize(meanFwd,[maxLength 1], 'nearest');
    meanAngrs = imresize(meanAng,[maxLength 1], 'nearest');
    meanSidrs = imresize(meanSid,[maxLength 1], 'nearest');
    semFwdrs = imresize(semFwd,[maxLength 1], 'nearest');
    semAngrs = imresize(semAng,[maxLength 1], 'nearest');
    semSidrs = imresize(semSid,[maxLength 1], 'nearest');

    % pull max values for each directional velocity
    group_max_forward = max(groupFwd);
    group_max_angular = max(groupAng);
    group_max_sideway = max(groupSid);
    % calculate mean and SEM for max values
    mean_max_forward = mean(group_max_forward,'omitnan');
    mean_max_angular = mean(group_max_angular,'omitnan');
    mean_max_sideway = mean(group_max_sideway,'omitnan');
    sem_max_forward = std(group_max_forward,'omitnan')/sqrt(nFliesActual);
    sem_max_angular = std(group_max_angular,'omitnan')/sqrt(nFliesActual);
    sem_max_sideway = std(group_max_sideway,'omitnan')/sqrt(nFliesActual);

    % forward
    subplot(3,4,[1 5 9])
    if s ==1
        fc = gca;
        % bar position
        yyaxis right
        plot(x,sweepPos{maxIdx},'k','LineWidth',2);
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
        cT = mean(x);
        xline(cT,'Color','k')
    end
    yyaxis left
    r(1) = patch([x'; flipud(x')],[(meanFwdrs-semFwdrs); flipud((meanFwdrs+semFwdrs))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(x,meanFwdrs,'LineStyle','-','LineWidth',2,'Marker','none','Color',cm(s,:)); hold on
    % adjust axes
    axis tight
    ylim(fwd_limit)
    ylabel(['Average Pursuit ' settings.velLabel{1}])
    hold on
    % angular
    subplot(3,4,[2 6 10])
    if s ==1
        fc = gca;
        % bar position
        yyaxis right
        plot(x,sweepPos{maxIdx},'k','LineWidth',2)
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
        cT = mean(x);
        xline(cT,'Color','k')
    end
    yyaxis left
    r(1) = patch([x'; flipud(x')],[(meanAngrs-semAngrs); flipud((meanAngrs+semAngrs))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(x,meanAngrs,'LineStyle','-','LineWidth',2,'Marker','none','Color',cm(s,:)); hold on
    % adjust axes
    axis tight
    ylim(ang_limit)
    ylabel(['Average Pursuit ' settings.velLabel{2}])
    hold on
    % sideway
    subplot(3,4,[3 7 11])
    if s == 1
        fc = gca;
        % bar position
        yyaxis right
        plot(x,sweepPos{maxIdx},'k','LineWidth',2)
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
        cT = mean(x);
        xline(cT,'Color','k')
    end
    yyaxis left
    r(1) = patch([x'; flipud(x')],[(meanSidrs-semSidrs); flipud((meanSidrs+semSidrs))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(x,meanSidrs,'LineStyle','-','LineWidth',2,'Marker','none','Color',cm(s,:)); hold on
    % adjust axes
    axis tight
    ylim(sid_limit)
    ylabel(['Average Pursuit ' settings.velLabel{3}])
    hold on

    % add a bit of jitter such that points are not overlapping
    jitter = 0.1;
    tt_j = s + (-jitter+(jitter*2)*rand(1,nFliesThresh));

    % max forward
    subplot(3,4,4)
    plot(tt_j,group_max_forward,'.','Color', [0.8 0.8 0.8]); hold on
    errorbar(xval(s),mean_max_forward,sem_max_forward,'o','Color',cm(s,:))
    ylim(fwd_limit)
    xlim([0 nSpeeds+1])
    xticks(xval)
    xticklabels(trialTypes)
    ylabel('Max Forward Velocity')
    hold on
    % max angular
    subplot(3,4,8)
    plot(tt_j,group_max_angular,'.','Color', [0.8 0.8 0.8]); hold on
    errorbar(xval(s),mean_max_angular,sem_max_angular,'o','Color',cm(s,:))
    ylim([0 ang_limit(2)])
    xlim([0 nSpeeds+1])
    xticks(xval)
    xticklabels(trialTypes)
    ylabel('Max Angular Velocity')
    hold on
    % max sideway
    subplot(3,4,12)
    plot(tt_j,group_max_sideway,'.','Color', [0.8 0.8 0.8]); hold on
    errorbar(xval(s),mean_max_sideway,sem_max_sideway,'o','Color',cm(s,:))
    ylim([0 sid_limit(2)])
    xlim([0 nSpeeds+1])
    xticks(xval)
    xticklabels(trialTypes)
    ylabel('Max Sideway Velocity')
    hold on
end

sgtitle([strrep(filebase,'_',' ') ' (n = ' num2str(nFliesThresh) ')'])
% save plot
cd(folder.summary)
plotname = 'velvpos_across';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vector plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');


%% comparison: plot pursuit behavior vs target motion R+L
disp('Comparison: analyzing pursuit behavior across conditions R+L...')

% initialize
figure; set(gcf,'Position',[100 100 1500 600])
[maxLength,maxIdx] = max(sweepLengths); % reshape everything to longest sweep size
x = 1:maxLength;
xval = 1:nSpeeds;

for s=1:nSpeeds
    % pull individual trials for each directional velocity
    groupFwd = cell2mat(posvfwdRL(:,s)');
    groupAng = cell2mat(posvangRL(:,s)');
    groupSid = cell2mat(posvsidRL(:,s)');
    nFliesActual = sum(~isnan(groupFwd(1,:)));
    % calculate mean and sem
    meanFwd = mean(groupFwd,2,'omitnan');
    meanAng = mean(groupAng,2,'omitnan');
    meanSid = mean(groupSid,2,'omitnan');
    semFwd = std(groupFwd,0,2,'omitnan')./sqrt(nFliesActual);
    semAng = std(groupAng,0,2,'omitnan')./sqrt(nFliesActual);
    semSid = std(groupSid,0,2,'omitnan')./sqrt(nFliesActual);

    % resize to plot together
    meanFwdrs = imresize(meanFwd,[maxLength 1], 'nearest');
    meanAngrs = imresize(meanAng,[maxLength 1], 'nearest');
    meanSidrs = imresize(meanSid,[maxLength 1], 'nearest');
    semFwdrs = imresize(semFwd,[maxLength 1], 'nearest');
    semAngrs = imresize(semAng,[maxLength 1], 'nearest');
    semSidrs = imresize(semSid,[maxLength 1], 'nearest');

    % pull max values for each directional velocity
    group_max_forward = max(groupFwd);
    group_max_angular = max(groupAng);
    group_max_sideway = max(groupSid);
    % calculate mean and SEM for max values
    mean_max_forward = mean(group_max_forward,'omitnan');
    mean_max_angular = mean(group_max_angular,'omitnan');
    mean_max_sideway = mean(group_max_sideway,'omitnan');
    sem_max_forward = std(group_max_forward,'omitnan')/sqrt(nFliesActual);
    sem_max_angular = std(group_max_angular,'omitnan')/sqrt(nFliesActual);
    sem_max_sideway = std(group_max_sideway,'omitnan')/sqrt(nFliesActual);

    % forward
    subplot(3,4,[1 5 9])
    if s ==1
        fc = gca;
        % bar position
        yyaxis right
        plot(x,sweepPos{maxIdx},'k','LineWidth',2);
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
        cT = mean(x);
        xline(cT,'Color','k')
    end
    yyaxis left
    r(1) = patch([x'; flipud(x')],[(meanFwdrs-semFwdrs); flipud((meanFwdrs+semFwdrs))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(x,meanFwdrs,'LineStyle','-','LineWidth',2,'Marker','none','Color',cm(s,:)); hold on
    % adjust axes
    axis tight
    ylim(fwd_limit)
    ylabel(['Average Pursuit ' settings.velLabel{1}])
    hold on
    % angular
    subplot(3,4,[2 6 10])
    if s ==1
        fc = gca;
        % bar position
        yyaxis right
        plot(x,sweepPos{maxIdx},'k','LineWidth',2)
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
        cT = mean(x);
        xline(cT,'Color','k')
    end
    yyaxis left
    r(1) = patch([x'; flipud(x')],[(meanAngrs-semAngrs); flipud((meanAngrs+semAngrs))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(x,meanAngrs,'LineStyle','-','LineWidth',2,'Marker','none','Color',cm(s,:)); hold on
    % adjust axes
    axis tight
    ylim(ang_limit)
    ylabel(['Average Pursuit ' settings.velLabel{2}])
    hold on
    % sideway
    subplot(3,4,[3 7 11])
    if s == 1
        fc = gca;
        % bar position
        yyaxis right
        plot(x,sweepPos{maxIdx},'k','LineWidth',2)
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
        cT = mean(x);
        xline(cT,'Color','k')
    end
    yyaxis left
    r(1) = patch([x'; flipud(x')],[(meanSidrs-semSidrs); flipud((meanSidrs+semSidrs))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(x,meanSidrs,'LineStyle','-','LineWidth',2,'Marker','none','Color',cm(s,:)); hold on
    % adjust axes
    axis tight
    ylim(sid_limit)
    ylabel(['Average Pursuit ' settings.velLabel{3}])
    hold on

    % add a bit of jitter such that points are not overlapping
    jitter = 0.1;
    tt_j = s + (-jitter+(jitter*2)*rand(1,nFliesThresh));

    % max forward
    subplot(3,4,4)
    plot(tt_j,group_max_forward,'.','Color', [0.8 0.8 0.8]); hold on
    errorbar(xval(s),mean_max_forward,sem_max_forward,'o','Color',cm(s,:))
    ylim(fwd_limit)
    xlim([0 nSpeeds+1])
    xticks(xval)
    xticklabels(trialTypes)
    ylabel('Max Forward Velocity')
    hold on
    % max angular
    subplot(3,4,8)
    plot(tt_j,group_max_angular,'.','Color', [0.8 0.8 0.8]); hold on
    errorbar(xval(s),mean_max_angular,sem_max_angular,'o','Color',cm(s,:))
    ylim([0 ang_limit(2)])
    xlim([0 nSpeeds+1])
    xticks(xval)
    xticklabels(trialTypes)
    ylabel('Max Angular Velocity')
    hold on
    % max sideway
    subplot(3,4,12)
    plot(tt_j,group_max_sideway,'.','Color', [0.8 0.8 0.8]); hold on
    errorbar(xval(s),mean_max_sideway,sem_max_sideway,'o','Color',cm(s,:))
    ylim([0 sid_limit(2)])
    xlim([0 nSpeeds+1])
    xticks(xval)
    xticklabels(trialTypes)
    ylabel('Max Sideway Velocity')
    hold on
end

sgtitle([strrep(filebase,'_',' ') ' R+L (n = ' num2str(nFliesThresh) ')'])
% save plot
cd(folder.summary)
plotname = 'velvpos_acrossRL';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vector plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');


%% by trial: plot spike rate versus directional velocity
cd(folder.summary)
disp('By trial type: plotting spike rate vs behavior...')

for s = 1:nSpeeds
    % initialize
    figure; set(gcf,'Position',[100 100 1000 500])

    % fetch data
    groupFwd = cat(2,srvfwd_vel{:,s});
    groupAng = cat(2,srvang_vel{:,s});
    groupSid = cat(2,srvsid_vel{:,s});
    % calculate means
    velFwdMean(:,s) = mean(groupFwd,2,'omitnan');
    velAngMean(:,s) = mean(groupAng,2,'omitnan');
    velSidMean(:,s) = mean(groupSid,2,'omitnan');
    % calculate sem
    velFwdSEM(:,s) = std(groupFwd,0,2,'omitnan')/sqrt(nFliesThresh);
    velAngSEM(:,s) = std(groupAng,0,2,'omitnan')/sqrt(nFliesThresh);
    velSidSEM(:,s) = std(groupSid,0,2,'omitnan')/sqrt(nFliesThresh);

    % plot forward
    subplot(1,3,1)
    plot(binFwd_v,cat(2,srvfwd_vel{:,s}),'-','Color',settings.trialColor,'LineWidth',settings.lwTri)
    hold on
    plot(binFwd_v,velFwdMean(:,s),'Color',settings.velColor{1},'LineWidth',settings.lwAvg)
    ylim(sr_limit)
    xlim(fwd_limit2)
    xlabel(settings.velLabel{1})
    % plot angular
    subplot(1,3,2)
    plot(binAng_v,cat(2,srvang_vel{:,s}),'-','Color',settings.trialColor,'LineWidth',settings.lwTri)
    hold on
    plot(binAng_v,velAngMean(:,s),'Color',settings.velColor{2},'LineWidth',settings.lwAvg)
    ylim(sr_limit)
    xlim(ang_limit2)
    xline(0)
    xlabel(settings.velLabel{2})
    % plot sideway
    subplot(1,3,3)
    plot(binSid_v,cat(2,srvsid_vel{:,s}),'-','Color',settings.trialColor,'LineWidth',settings.lwTri)
    hold on
    plot(binSid_v,velSidMean(:,s),'Color',settings.velColor{3},'LineWidth',settings.lwAvg)
    ylim(sr_limit)
    xlim(sid_limit2)
    xline(0)
    xlabel(settings.velLabel{3})

    thisletter = settings.letters(s);
    sgtitle([strrep(filebase,'_',' ') ' ' trialTypes{s} ' v Velocity (n = ' num2str(nFliesThresh) ')'])
    % save plot
    cd(folder.summary)
    plotname = ['srvvel_' thisletter];
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vector plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');
end
disp('Complete.')



%% comparison: plot spike rate versus directional velocity
disp('Comparison: plotting spike rate vs behavior...')
% initialize
figure; set(gcf,'Position',[100 100 1000 500])
% check for nans
velFwdMean(isnan(velFwdMean)) = 0;
velAngMean(isnan(velAngMean)) = 0;
velSidMean(isnan(velSidMean)) = 0;
velFwdSEM(isnan(velFwdSEM)) = 0;
velAngSEM(isnan(velAngSEM)) = 0;
velSidSEM(isnan(velSidSEM)) = 0;

for s = 1:nSpeeds
    % plot forward
    subplot(1,3,1)
    r(1) = patch([binFwd_v'; flipud(binFwd_v')],[(velFwdMean(:,s)-velFwdSEM(:,s)); flipud((velFwdMean(:,s)+velFwdSEM(:,s)))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(binFwd_v,velFwdMean(:,s),'Color',cm(s,:),'LineWidth',settings.lwAvg)
    ylim(sr_limitA)
    xlim(fwd_limit2)
    xlabel(settings.velLabel{1})
    ylabel(['Mean ' settings.spkLabel])
    % plot angular
    subplot(1,3,2)
    r(1) = patch([binAng_v'; flipud(binAng_v')],[(velAngMean(:,s)-velAngSEM(:,s)); flipud((velAngMean(:,s)+velAngSEM(:,s)))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(binAng_v,velAngMean(:,s),'Color',cm(s,:),'LineWidth',settings.lwAvg)
    ylim(sr_limitA)
    xlim(ang_limit2)
    xline(0)
    xlabel(settings.velLabel{2})
    % plot sideway
    subplot(1,3,3)
    r(1) = patch([binSid_v'; flipud(binSid_v')],[(velSidMean(:,s)-velSidSEM(:,s)); flipud((velSidMean(:,s)+velSidSEM(:,s)))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(binSid_v,velSidMean(:,s),'Color',cm(s,:),'LineWidth',settings.lwAvg)
    ylim(sr_limitA)
    xlim(sid_limit2)
    xline(0)
    xlabel(settings.velLabel{3})
end

sgtitle([strrep(filebase,'_',' ') ' v Velocity (n = ' num2str(nFliesThresh) ')'])
% save plot
cd(folder.summary)
plotname = 'srvvel_across';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vector plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

disp('Complete.')


%% comparison: plot expected R-L spike rate versus directional velocity
disp('Comparison: plotting expected R-L spike rate vs behavior...')
% initialize
figure; set(gcf,'Position',[100 100 700 500])

% find ipsi vs contra indices
ipsiAngIdx = find(binAng_v>=0);
cntrAngIdx = find(binAng_v<=0);
ipsiSidIdx = find(binSid_v>=0);
cntrSidIdx = find(binSid_v<=0);

angBinD = binAng_v(ipsiAngIdx);
sidBinD = binSid_v(ipsiSidIdx);

for s = 1:nSpeeds
    % pull all data for this trial type
    groupAng = cat(2,srvang_vel{:,s});
    groupSid = cat(2,srvsid_vel{:,s});
    % calculate R-L difference
    thisAngDiff = groupAng(ipsiAngIdx,:)-flip(groupAng(cntrAngIdx,:));
    thisSidDiff = groupSid(ipsiSidIdx,:)-flip(groupSid(cntrSidIdx,:));

    % calculate R-L diff means
    thisAngDiffMean(:,s) = mean(thisAngDiff,2,'omitnan');
    thisSidDiffMean(:,s) = mean(thisSidDiff,2,'omitnan');
    % calculate sem
    thisAngDiffSEM(:,s) = std(thisAngDiff,0,2,'omitnan')/sqrt(nFliesThresh);
    thisSidDiffSEM(:,s) = std(thisSidDiff,0,2,'omitnan')/sqrt(nFliesThresh);
    % check for nans
    thisAngDiffMean(isnan(thisAngDiffMean)) = 0;
    thisSidDiffMean(isnan(thisSidDiffMean)) = 0;
    thisAngDiffSEM(isnan(thisAngDiffSEM)) = 0;
    thisSidDiffSEM(isnan(thisSidDiffSEM)) = 0;

    % plot angular
    subplot(1,2,1)
    r(1) = patch([angBinD'; flipud(angBinD')],[(thisAngDiffMean(:,s)-thisAngDiffSEM(:,s)); flipud((thisAngDiffMean(:,s)+thisAngDiffSEM(:,s)))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(angBinD,thisAngDiffMean(:,s),'Color',cm(s,:),'LineWidth',settings.lwAvg)
    ylim(sr_limitD)
    xlim(ang_limit2)
    xline(0)
    xlabel(settings.velLabel{2})
    if s == 1
        ylabel(['R-L ' settings.spkLabel])
    end
    % plot sideway
    subplot(1,2,2)
    r(1) = patch([sidBinD'; flipud(sidBinD')],[(thisSidDiffMean(:,s)-thisSidDiffSEM(:,s)); flipud((thisSidDiffMean(:,s)+thisSidDiffSEM(:,s)))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(sidBinD,thisSidDiffMean(:,s),'Color',cm(s,:),'LineWidth',settings.lwAvg)
    ylim(sr_limitD)
    xlim(sid_limit2)
    xline(0)
    xlabel(settings.velLabel{3})
end

sgtitle([strrep(filebase,'_',' ') ' R-L Difference' '(n = ' num2str(nFliesThresh) ')'])
% save plot
cd(folder.summary)
plotname = 'srvvel_difference_across';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vector plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

disp('Complete.')


%% by trial: plot spike rate versus directional velocity w/lag
cd(folder.summary)
disp('By trial type: plotting spike rate vs behavior w/lag...')

for s = 1:nSpeeds
    % initialize
    figure; set(gcf,'Position',[100 100 1000 500])

    % fetch data
    groupFwd = cat(2,srvfwd_velL{:,s});
    groupAng = cat(2,srvang_velL{:,s});
    groupSid = cat(2,srvsid_velL{:,s});
    % calculate means
    velFwdMeanL(:,s) = mean(groupFwd,2,'omitnan');
    velAngMeanL(:,s) = mean(groupAng,2,'omitnan');
    velSidMeanL(:,s) = mean(groupSid,2,'omitnan');
    % calculate sem
    velFwdSEML(:,s) = std(groupFwd,0,2,'omitnan')/sqrt(nFliesThresh);
    velAngSEML(:,s) = std(groupAng,0,2,'omitnan')/sqrt(nFliesThresh);
    velSidSEML(:,s) = std(groupSid,0,2,'omitnan')/sqrt(nFliesThresh);

    % plot forward
    subplot(1,3,1)
    plot(binFwd_v,cat(2,srvfwd_velL{:,s}),'-','Color',settings.trialColor,'LineWidth',settings.lwTri)
    hold on
    plot(binFwd_v,velFwdMeanL(:,s),'Color',settings.velColor{1},'LineWidth',settings.lwAvg)
    ylim(sr_limit)
    xlim(fwd_limit2)
    xlabel(settings.velLabel{1})
    % plot angular
    subplot(1,3,2)
    plot(binAng_v,cat(2,srvang_velL{:,s}),'-','Color',settings.trialColor,'LineWidth',settings.lwTri)
    hold on
    plot(binAng_v,velAngMeanL(:,s),'Color',settings.velColor{2},'LineWidth',settings.lwAvg)
    ylim(sr_limit)
    xlim(ang_limit2)
    xline(0)
    xlabel(settings.velLabel{2})
    % plot sideway
    subplot(1,3,3)
    plot(binSid_v,cat(2,srvsid_velL{:,s}),'-','Color',settings.trialColor,'LineWidth',settings.lwTri)
    hold on
    plot(binSid_v,velSidMeanL(:,s),'Color',settings.velColor{3},'LineWidth',settings.lwAvg)
    ylim(sr_limit)
    xlim(sid_limit2)
    xline(0)
    xlabel(settings.velLabel{3})

    thisletter = settings.letters(s);
    sgtitle([strrep(filebase,'_',' ') ' ' trialTypes{s} ' v Velocity w/lag (n = ' num2str(nFliesThresh) ')'])
    % save plot
    cd(folder.summary)
    plotname = ['srvvel_lag_' thisletter];
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vector plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');
end
disp('Complete.')



%% comparison: plot spike rate versus directional velocity w/lag
disp('Comparison: plotting spike rate vs behavior w/lag...')
% initialize
figure; set(gcf,'Position',[100 100 1000 500])
% check for nans
velFwdMeanL(isnan(velFwdMeanL)) = 0;
velAngMeanL(isnan(velAngMeanL)) = 0;
velSidMeanL(isnan(velSidMeanL)) = 0;
velFwdSEML(isnan(velFwdSEML)) = 0;
velAngSEML(isnan(velAngSEML)) = 0;
velSidSEML(isnan(velSidSEML)) = 0;

for s = 1:nSpeeds
    % plot forward
    subplot(1,3,1)
    r(1) = patch([binFwd_v'; flipud(binFwd_v')],[(velFwdMeanL(:,s)-velFwdSEML(:,s)); flipud((velFwdMeanL(:,s)+velFwdSEML(:,s)))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(binFwd_v,velFwdMeanL(:,s),'Color',cm(s,:),'LineWidth',settings.lwAvg)
    ylim(sr_limitA)
    xlim(fwd_limit2)
    xlabel(settings.velLabel{1})
    ylabel(['Mean ' settings.spkLabel])
    % plot angular
    subplot(1,3,2)
    r(1) = patch([binAng_v'; flipud(binAng_v')],[(velAngMeanL(:,s)-velAngSEML(:,s)); flipud((velAngMeanL(:,s)+velAngSEML(:,s)))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(binAng_v,velAngMeanL(:,s),'Color',cm(s,:),'LineWidth',settings.lwAvg)
    ylim(sr_limitA)
    xlim(ang_limit2)
    xline(0)
    xlabel(settings.velLabel{2})
    % plot sideway
    subplot(1,3,3)
    r(1) = patch([binSid_v'; flipud(binSid_v')],[(velSidMeanL(:,s)-velSidSEM(:,s)); flipud((velSidMeanL(:,s)+velSidSEML(:,s)))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(binSid_v,velSidMeanL(:,s),'Color',cm(s,:),'LineWidth',settings.lwAvg)
    ylim(sr_limitA)
    xlim(sid_limit2)
    xline(0)
    xlabel(settings.velLabel{3})
end

sgtitle([strrep(filebase,'_',' ') ' v Velocity w/lag (n = ' num2str(nFliesThresh) ')'])
% save plot
cd(folder.summary)
plotname = 'srvvel_lag_across';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vector plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

disp('Complete.')


%% by trial: plot spike rate versus directional acceleration
cd(folder.summary)
disp('By trial type: plotting spike rate vs acceleration...')
for s = 1:nSpeeds
    % initialize
    figure; set(gcf,'Position',[100 100 1000 500])

    % fetch data
    groupFwd = cat(2,srvfwd_acc{:,s});
    groupAng = cat(2,srvang_acc{:,s});
    groupSid = cat(2,srvsid_acc{:,s});
    % calculate means
    accFwdMean(:,s) = mean(groupFwd,2,'omitnan');
    accAngMean(:,s) = mean(groupAng,2,'omitnan');
    accSidMean(:,s) = mean(groupSid,2,'omitnan');
    % calculate sem
    accFwdSEM(:,s) = std(groupFwd,0,2,'omitnan')/sqrt(nFliesThresh);
    accAngSEM(:,s) = std(groupAng,0,2,'omitnan')/sqrt(nFliesThresh);
    accSidSEM(:,s) = std(groupSid,0,2,'omitnan')/sqrt(nFliesThresh);

    % plot forward
    subplot(1,3,1)
    plot(binFwd_a,cat(2,srvfwd_acc{:,s}),'-','Color',settings.trialColor,'LineWidth',settings.lwTri)
    hold on
    plot(binFwd_a,accFwdMean(:,s),'Color',settings.velColor{1},'LineWidth',settings.lwAvg)
    ylim(sr_limit)
    xlim(fwd_limita)
    xline(0)
    xlabel(settings.accLabel{1})
    % plot angular
    subplot(1,3,2)
    plot(binAng_a,cat(2,srvang_acc{:,s}),'-','Color',settings.trialColor,'LineWidth',settings.lwTri)
    hold on
    plot(binAng_a,accAngMean(:,s),'Color',settings.velColor{2},'LineWidth',settings.lwAvg)
    ylim(sr_limit)
    xlim(ang_limita)
    xline(0)
    xlabel(settings.accLabel{2})
    % plot sideway
    subplot(1,3,3)
    plot(binSid_a,cat(2,srvsid_acc{:,s}),'-','Color',settings.trialColor,'LineWidth',settings.lwTri)
    hold on
    plot(binSid_a,accSidMean(:,s),'Color',settings.velColor{3},'LineWidth',settings.lwAvg)
    ylim(sr_limit)
    xlim(sid_limita)
    xline(0)
    xlabel(settings.accLabel{3})

    thisletter = settings.letters(s);
    sgtitle([strrep(filebase,'_',' ') ' ' trialTypes{s} ' v Acceleration (n = ' num2str(nFliesThresh) ')'])
    % save plot
    cd(folder.summary)
    plotname = ['srvacc_' thisletter];
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vector plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');
end
disp('Complete.')



%% comparison: plot spike rate versus directional acceleration
disp('Comparison: plotting spike rate vs acceleration...')
% initialize
figure; set(gcf,'Position',[100 100 1000 500])
% check for nans
velFwdMean(isnan(velFwdMean)) = 0;
velAngMean(isnan(velAngMean)) = 0;
velSidMean(isnan(velSidMean)) = 0;
velFwdSEM(isnan(velFwdSEM)) = 0;
velAngSEM(isnan(velAngSEM)) = 0;
velSidSEM(isnan(velSidSEM)) = 0;

for s = 1:nSpeeds
    % plot forward
    subplot(1,3,1)
    r(1) = patch([binFwd_a'; flipud(binFwd_a')],[(accFwdMean(:,s)-accFwdSEM(:,s)); flipud((accFwdMean(:,s)+accFwdSEM(:,s)))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(binFwd_a,accFwdMean(:,s),'Color',cm(s,:),'LineWidth',settings.lwAvg)
    ylim(sr_limitA)
    xlim(fwd_limita)
    xline(0)
    xlabel(settings.accLabel{1})
    ylabel(['Mean ' settings.spkLabel])
    % plot angular
    subplot(1,3,2)
    r(1) = patch([binAng_a'; flipud(binAng_a')],[(accAngMean(:,s)-accAngSEM(:,s)); flipud((accAngMean(:,s)+accAngSEM(:,s)))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(binAng_a,accAngMean(:,s),'Color',cm(s,:),'LineWidth',settings.lwAvg)
    ylim(sr_limitA)
    xlim(ang_limita)
    xline(0)
    xlabel(settings.accLabel{2})
    % plot sideway
    subplot(1,3,3)
    r(1) = patch([binSid_a'; flipud(binSid_a')],[(accSidMean(:,s)-accSidSEM(:,s)); flipud((accSidMean(:,s)+accSidSEM(:,s)))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',cm(s,:));
    hold on
    plot(binSid_a,accSidMean(:,s),'Color',cm(s,:),'LineWidth',settings.lwAvg)
    ylim(sr_limitA)
    xlim(sid_limita)
    xline(0)
    xlabel(settings.accLabel{3})
end

sgtitle([strrep(filebase,'_',' ') ' v Acceleration (n = ' num2str(nFliesThresh) ')'])
% save plot
cd(folder.summary)
plotname = 'srvaccel_across';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vector plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

disp('Complete.')


%% determine lag and cross-correlation for each trial type
disp('Analyzing cross-correlation for firing vs pursuit behavior across conditions...')
cd(folder.summary)

% initialize
figure; set(gcf,'Position',[100 100 1200 800])
tiledlayout(3,nSpeeds+1,'TileSpacing','compact')

% for each directional velocity
for v = 1:3
    % initialize
    peak_mean = [];
    peak_sem = [];
    switch v
        case 1
            r_val = r_val_fwd;
            r_peaks = peak_fwd;
            thisV = 'Forward Rval';
        case 2
            r_val = r_val_ang;
            r_peaks = peak_ang;
            thisV = 'Angular Rval';
        case 3
            r_val = r_val_sid;
            r_peaks = peak_sid;
            thisV = 'Sideways Rval';
    end

    % for each speed
    for s = 1:nSpeeds
        % calculate mean
        r_mean = mean(r_val(:,:,s),2,'omitnan');
        % calculate peak lag mean and sem
        peak_mean(s) = mean(r_peaks(:,s),'omitnan');
        peak_sem(s) = std(r_peaks(:,s),'omitnan')./sqrt(nFliesThresh);

        % plot cross correlations
        nexttile; hold on
        plot(lag_t,r_val(:,:,s),':','LineWidth',0.5,'Color',settings.trialColor)
        plot(lag_t,r_mean,'LineWidth',1.5,'Color',settings.velColor{v})
        xline(0); ylim(xc_range); xlabel('Lag (msec)');
        if s==1
            ylabel(thisV)
        end
        if v==1
            title(trialTypes{s})
        end
    end
    % plot peak lags
    nexttile; hold on
    plot(1:nSpeeds,r_peaks,'.','Color',settings.trialColor)
    errorbar(1:nSpeeds,peak_mean,peak_sem,'o','Color',settings.velColor{v})
    axis padded; yline(0); ylim(xc_lim); xticks(1:nSpeeds)
    xticklabels(trialTypes); ylabel('Peak Lag (msec)')
    % add labels above each point
    for s = 1:nSpeeds
        text(s, 250, num2str(peak_mean(s)),'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom','FontSize',8,'Rotation',90);
    end
end

sgtitle([strrep(filebase,'_',' ') ' XCorr' '(n = ' num2str(nFliesThresh) ')'])
% save plot
cd(folder.summary)
plotname = 'osc_xcorr';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');


%% determine auto-correlation for each trial type
% disp('Analyzing auto-correlation for behavior across conditions...')
% cd(folder.summary)
%
% % initialize
% figure; set(gcf,'Position',[100 100 1050 800])
% tiledlayout(3,nSpeeds,'TileSpacing','compact')
%
% % for each directional velocity
% for v = 1:3
%     % initialize
%     peak_mean = [];
%     peak_sem = [];
%     switch v
%         case 1
%             r_val = ar_val_fwd;
%             thisV = 'Forward Rval';
%         case 2
%             r_val = ar_val_ang;
%             thisV = 'Angular Rval';
%         case 3
%             r_val = ar_val_sid;
%             thisV = 'Sideways Rval';
%     end
%
%     % for each speed
%     for s = 1:nSpeeds
%         % calculate mean
%         r_mean = mean(r_val(:,:,s),2,'omitnan');
%
%         % plot cross correlations
%         nexttile; hold on
%         plot(lag_t,r_val(:,:,s),':','LineWidth',0.5,'Color',settings.trialColor)
%         plot(lag_t,r_mean,'LineWidth',1.5,'Color',settings.velColor{v})
%         xline(0); ylim(xc_range); xlabel('Lag (msec)');
%         if s==1
%             ylabel(thisV)
%         end
%         if v==1
%             title(trialTypes{s})
%         end
%     end
% end
%
% sgtitle([strrep(filebase,'_',' ') ' AutoCorr' '(n = ' num2str(nFliesThresh) ')'])
% % save plot
% cd(folder.summary)
% plotname = 'osc_autocorr';
% saveas(gcf,[plotname '.png']);
% copyfile([plotname '.png'], folder.dropbox,'f');
% % save vectorized plot
% cd(folder.vector)
% set(gcf,'renderer','Painters')
% saveas(gcf, [plotname '.svg'])
% copyfile([plotname '.svg'], folder.dropbox,'f');


%% analyze pursuit performance/tracking as a function of time
disp('Analyzing pursuit performance over time...')

% initialize
normVigor = cell(size(allVigor));
trackingIdx = cell(size(allVigor));

% for each fly
for nt = 1:nFliesThresh
    %normalize vigor to that fly's maximum
    % fetch vigor from this fly
    thisVigor = allVigor(nt,:);
    % determine max
    maxVigor = max(cellfun(@(x) max(x(:)), thisVigor));
    normVigor(nt,:) = cellfun(@(x) x / maxVigor, thisVigor, 'UniformOutput', false);

    % calculate tracking index
    % for each trial condition
    for s = 1:nSpeeds
        trackingIdx{nt,s} = allFidelity{nt,s} .* normVigor{nt,s};

        % calculate median and sem for tracking
        trackingMedian{nt,s} = median(trackingIdx{nt,s},1,'omitnan');
        trackingStd{nt,s} = std(trackingIdx{nt,s},0,1,'omitnan');

        % calculate median and sem for tracking
        trackingPeak{nt,s} = max(trackingIdx{nt,s});
    end
    % store number of trials for this fly
    nTrials(nt,1) = size(trackingIdx{nt,s},2);
end

% plot tracking index across all flies for first trial only
% initialize
figure; set(gcf,'Position',[50 50 1800 900])
tiledlayout(nSpeeds,1,'TileSpacing','compact')

% for each speed
for s = 1:nSpeeds
    thistracking = [];
    % for each fly
    for nt = 1:nFliesThresh
        % fetch tracking for this speed
        thistracking(:,nt) = trackingIdx{nt,s}(:,1);
    end
    nexttile; hold on
    % plot tracking per fly
    plot(w_time,thistracking,'Color',settings.trialColor,'LineWidth',settings.lwTri)
    % plot mean across flies
    trackingMean = mean(thistracking,2);
    plot(w_time,trackingMean,'Color',settings.velColor{2},'LineWidth',2)
    axis padded; yline(0); ylim([-0.4 1])
    title(trialTypes{s}); xlabel('Time (msec.)'); ylabel('Tracking Idx')
end
sgtitle([strrep(filebase,'_',' ') ' Tracking for First Trial' ' (n = ' num2str(nFliesThresh) ')'])
% save plot
cd(folder.summary)
plotname = 'osc_tracking_firsttrial';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

% plot median tracking index across all flies and all trials
% initialize
figure; set(gcf,'Position',[100 100 1000 900])
tiledlayout(nSpeeds,1,'TileSpacing','compact')

% for each speed
for s = 1:nSpeeds
    nexttile; hold on
    % for each fly
    for nt = 1:nFliesThresh
        % fetch this mean, std, and number of trials
        nX = length(trackingMedian{nt,s});
        thisMedian = trackingMedian{nt,s};
        thisSTD = trackingStd{nt,s};
        % plot median
        errorbar(1:nX, thisMedian, thisSTD, 'Marker', 'none', 'CapSize', 0)
    end
    axis padded; yline(0); ylim([-0.2 1])
    title(trialTypes{s}); xlabel('Trial N'); ylabel('Tracking Idx')
end
sgtitle([strrep(filebase,'_',' ') ' Median Tracking v Trials' ' (n = ' num2str(nFliesThresh) ')'])
% save plot
cd(folder.summary)
plotname = 'osc_tracking_alltrials';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

% plot peak tracking index across all flies and all trials
% initialize
figure; set(gcf,'Position',[100 100 1000 900])
tiledlayout(nSpeeds,1,'TileSpacing','compact')
maxTrials = max(nTrials);
% for each speed
for s = 1:nSpeeds
    thisPeak = nan(maxTrials,nFliesThresh);
    % for each fly
    for nt = 1:nFliesThresh
        % fetch this peak, std, and number of trials
        nX = length(trackingPeak{nt,s});
        thisPeak(1:nX,nt) = trackingPeak{nt,s};

    end
    nexttile; hold on
    % plot peak per fly
    plot(1:maxTrials,thisPeak,'-','Color',settings.trialColor,'LineWidth',settings.lwTri)
    % calculate peak mean and sem across flies
    thisPeakMean = mean(thisPeak,2,'omitnan');
    thisPeakSEM = std(thisPeak,0,2,'omitnan')./sqrt(nFliesThresh);
    errorbar(1:maxTrials,thisPeakMean,thisPeakSEM,'Marker', 'none', 'CapSize', 0,'Color',settings.velColor{2})
    axis padded; yline(0); ylim([-0.2 1])
    title(trialTypes{s}); xlabel('Trial N'); ylabel('Peak Tracking Idx')
end
sgtitle([strrep(filebase,'_',' ') ' Peak Tracking v Trials' ' (n = ' num2str(nFliesThresh) ')'])
% save plot
cd(folder.summary)
plotname = 'osc_tracking_peaktrials';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');


%% analyze pursuit performance/tracking vs spike rate
disp('Analyzing pursuit performance vs spike rate...')

% initialize
idx_limit = [0 0.6];

srvidx = [];

% for each trial type
for s = 1:nSpeeds
    % for each fly
    for nt = 1:nFliesThresh
        thisIndex = trackingIdx{nt,s};
        thisSpikert  = allWinSR{nt,s};
        % generate average spikerate, binned by tracking index
        this_srvidx = spikert_binindex(thisSpikert,thisIndex);
        % store binned spikerate
        srvidx(:,nt,s) = this_srvidx(:,2);
    end
end
% store index bins
iBin = this_srvidx(:,1);

% calculate mean and sem
srvidx_mean = mean(srvidx,2,'omitnan');
srvidx_sem = std(srvidx,0,2,'omitnan')./sqrt(nFliesThresh);

% generate plot
% initialize
figure; set(gcf,'Position',[100 100 1800 900])
tiledlayout(2,nSpeeds,'TileSpacing','compact')
% plot trials and mean
for s = 1:nSpeeds
    nexttile; hold on
    % plot trials and mean
    plot(iBin,srvidx(:,:,s),'Color',settings.trialColor,'LineWidth',settings.lwTri)
    plot(iBin,srvidx_mean(:,:,s),'Color',settings.spkColor,'LineWidth',2)
    yline(0); xlim(idx_limit); ylim(sr_limit)
    xlabel('Pursuit Idx'); title(trialTypes{s})
    if s==1
        ylabel(['Mean ' settings.spkLabel])
    end
end

% plot mean and sem
for s = 1:nSpeeds
    nexttile; hold on
    % plot trials and mean
    dataAvail = ~isnan(srvidx_mean(:,:,s));
    r(1) = patch([iBin(dataAvail); flipud(iBin(dataAvail))],[(srvidx_mean(dataAvail,:,s)-srvidx_sem(dataAvail,:,s)); flipud((srvidx_mean(dataAvail,:,s)+srvidx_sem(dataAvail,:,s)))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none','FaceColor',settings.spkColor);
    plot(iBin,srvidx_mean(:,:,s),'Color',settings.spkColor,'LineWidth',2)
    yline(0); xlim(idx_limit); ylim(sr_limitA)
    xlabel('Pursuit Idx'); title(trialTypes{s})
    if s==1
        ylabel(['Mean ' settings.spkLabel])
    end
end

sgtitle([strrep(filebase,'_',' ') ' Tracking v Firing Rate' ' (n = ' num2str(nFliesThresh) ')'])
% save plot
cd(folder.summary)
plotname = 'srvtracking';
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

