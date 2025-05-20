% pipeline_kir_openloop
%
% Pipeline Function
% Pulls all processed files from ALL flies in a given experiment, performs
% necessary analyses and plots accordingly. This function compares the 
% behavior of flies with or without KIR perturbation across multiple 
% experimental conditions.
%
% INPUTS
% exptFolder - overarching experiment folder
%
% The function pools and analyzes data from oscillatory, motion pulse, 
% and optomotor experiments, generating behavior and velocity plots, 
% histograms, and summary figures. Results are saved in the form of 
% visualizations and statistical comparisons.
%
% 04/09/2024 - MC adapted from kir menotaxis pipeline
%
function pipeline_kir_openloop(exptFolder)
%% initialize
disp('STARTING ANALYSES FOR POOLED KIR PURSUIT...')
close all

% load in processing settings
settings = processSettings();


% set filename info and create necessary directories
filebase = strrep(exptFolder,' ','_');
% generate folder structures as needed
folder = generateFolders(exptFolder);

cd(folder.int)
% find all files for each experiment set
oscFiles = dir('*osc_int.mat');
mopFiles = dir('*mopulse_int.mat');
optFiles = dir('*opto_int.mat');
nFlies = length(optFiles);

%% set key variables

exptNames = {'osc';'mopulse';'opto'}; %experiment trial sets
hLabel = 'Normalized Probability';
hNorm = 1; %1 to normalize, 0 for counts
nExpt = length(exptNames);
nSweep = 13;

% set offset for overlapping data
o = 0.15;
offset = [-o; 0; o]; %used to offset overlapping plots

%% load in and pool oscillatory data
disp('Loading in and analyzing oscillatory datasets...')
nKIR = 0; nWT = 0; nNA = 0; nR = 1;
for e = 1:nFlies
    disp(['Processing fly ' num2str(e) '/' num2str(nFlies) '...'])
    % load this trial
    cd(folder.int)
    thisTrial = oscFiles(e).name;
    load(thisTrial)

    % find timepoints where the fly was running and for how long in total
    [run_forward,run_angular,~,~] = pursuitFinder(int_forward,int_angular,int_sideway,0,int_time,settings.runThreshB);
    int_time(end+1) = 60;
    thisRunTime = sum(int_time(sum(~isnan(run_forward))+1)); %sec
    thisRunSpeed = mean(run_forward,'all','omitnan'); %avg run speed
    thisTurnSpeed = mean(abs(run_angular),'all','omitnan'); %avg turn speed

    % analyze directional velocity during oscillation
    [sweepPos,meanFwd] = osc_v_output(int_panelps,int_forward,int_forward,settings.runThreshB,0,int_time,settings);
    [~,meanAng] = osc_v_output(int_panelps,int_forward,int_angular,settings.runThreshB,0,int_time,settings);
    [~,meanSid] = osc_v_output(int_panelps,int_forward,int_sideway,settings.runThreshB,0,int_time,settings);
    oscDur = size(meanAng,1);
    RLshift = round(oscDur/2);
    % combine R+L by rotating L
    meanFwdRL = mean([meanFwd circshift(meanFwd,RLshift)],2,'omitnan');
    meanAngRL = mean([meanAng -circshift(meanAng,RLshift)],2,'omitnan');
    meanSidRL = mean([meanSid -circshift(meanSid,RLshift)],2,'omitnan');

    % analyze directional acceleration during oscillation
    [accel] = velocity2acceleration(int_forward,int_angular,int_sideway,int_time(1:end-1)');
    [~,meanAccFwd] = osc_v_output(int_panelps,int_forward,accel.forward,settings.runThreshB,0,int_time,settings);
    [~,meanAccAng] = osc_v_output(int_panelps,int_forward,accel.angular,settings.runThreshB,0,int_time,settings);
    [~,meanAccSid] = osc_v_output(int_panelps,int_forward,accel.sideway,settings.runThreshB,0,int_time,settings);
    % combine R+L by rotating L
    meanAccFwdRL = mean([meanAccFwd circshift(meanAccFwd,RLshift)],2,'omitnan');
    meanAccAngRL = mean([meanAccAng -circshift(meanAccAng,RLshift)],2,'omitnan');
    meanAccSidRL = mean([meanAccSid -circshift(meanAccSid,RLshift)],2,'omitnan');

    % analyze velocity distributions
    [fwdHist,angHist,sidHist] = velocity_histogram(int_forward,int_angular,int_sideway,hNorm);

    % pool this trial data according to the assigned trial condition
    if contains(thisTrial,'KIR') %kir perturbation flies
        nKIR = nKIR+1;
        kirRunTime(nR,nKIR) = thisRunTime;
        kirRunSpeed(nR,nKIR) = thisRunSpeed;
        kirTurnSpeed(nR,nKIR) = thisTurnSpeed;

        kir_oscFwd(:,nKIR) = meanFwd;
        kir_oscAng(:,nKIR) = meanAng;
        kir_oscSid(:,nKIR) = meanSid;
        kir_oscFwdRL(:,nKIR) = meanFwdRL;
        kir_oscAngRL(:,nKIR) = meanAngRL;
        kir_oscSidRL(:,nKIR) = meanSidRL;

        kir_oscAccFwd(:,nKIR) = meanAccFwd;
        kir_oscAccAng(:,nKIR) = meanAccAng;
        kir_oscAccSid(:,nKIR) = meanAccSid;
        kir_oscAccFwdRL(:,nKIR) = meanAccFwdRL;
        kir_oscAccAngRL(:,nKIR) = meanAccAngRL;
        kir_oscAccSidRL(:,nKIR) = meanAccSidRL;

        kir_oscFwdHist(:,nKIR) = fwdHist(:,2);
        kir_oscAngHist(:,nKIR) = angHist(:,2);
        kir_oscSidHist(:,nKIR) = sidHist(:,2);
    elseif contains(thisTrial,'WT') %wildtype control flies
        nWT = nWT+1;
        wtRunTime(nR,nWT) = thisRunTime;
        wtRunSpeed(nR,nWT) = thisRunSpeed;
        wtTurnSpeed(nR,nWT) = thisTurnSpeed;

        wt_oscFwd(:,nWT) = meanFwd;
        wt_oscAng(:,nWT) = meanAng;
        wt_oscSid(:,nWT) = meanSid;
        wt_oscFwdRL(:,nWT) = meanFwdRL;
        wt_oscAngRL(:,nWT) = meanAngRL;
        wt_oscSidRL(:,nWT) = meanSidRL;

        wt_oscAccFwd(:,nWT) = meanAccFwd;
        wt_oscAccAng(:,nWT) = meanAccAng;
        wt_oscAccSid(:,nWT) = meanAccSid;
        wt_oscAccFwdRL(:,nWT) = meanAccFwdRL;
        wt_oscAccAngRL(:,nWT) = meanAccAngRL;
        wt_oscAccSidRL(:,nWT) = meanAccSidRL;

        wt_oscFwdHist(:,nWT) = fwdHist(:,2);
        wt_oscAngHist(:,nWT) = angHist(:,2);
        wt_oscSidHist(:,nWT) = sidHist(:,2);
    elseif contains(thisTrial,'NA') %na control flies
        nNA=nNA+1;
        naRunTime(nR,nNA) = thisRunTime;
        naRunSpeed(nR,nNA) = thisRunSpeed;
        naTurnSpeed(nR,nNA) = thisTurnSpeed;

        na_oscFwd(:,nNA) = meanFwd;
        na_oscAng(:,nNA) = meanAng;
        na_oscSid(:,nNA) = meanSid;
        na_oscFwdRL(:,nNA) = meanFwdRL;
        na_oscAngRL(:,nNA) = meanAngRL;
        na_oscSidRL(:,nNA) = meanSidRL;

        na_oscAccFwd(:,nNA) = meanAccFwd;
        na_oscAccAng(:,nNA) = meanAccAng;
        na_oscAccSid(:,nNA) = meanAccSid;
        na_oscAccFwdRL(:,nNA) = meanAccFwdRL;
        na_oscAccAngRL(:,nNA) = meanAccAngRL;
        na_oscAccSidRL(:,nNA) = meanAccSidRL;

        na_oscFwdHist(:,nNA) = fwdHist(:,2);
        na_oscAngHist(:,nNA) = angHist(:,2);
        na_oscSidHist(:,nNA) = sidHist(:,2);
    end
end
% pull stimulus position and time
oscSweepPos = sweepPos;
t_osc = int_time(1:length(oscSweepPos))*1000; %msec
disp('Complete.')

%% load in and pool motion pulse data
disp('Loading in and analyzing motion pulse datasets...')
nKIR = 0; nWT = 0; nNA = 0; nR = 2;
for e = 1:nFlies
    disp(['Processing fly ' num2str(e) '/' num2str(nFlies) '...'])
    % load this trial
    cd(folder.int)
    thisTrial = mopFiles(e).name;
    load(thisTrial)

    % find timepoints where the fly was running and for how long in total
    [run_forward,run_angular,~,~] = pursuitFinder(int_forward,int_angular,int_sideway,0,int_time,settings.runThreshB);
    int_time(end+1) = 60;
    thisRunTime = sum(int_time(sum(~isnan(run_forward))+1)); %sec
    thisRunSpeed = mean(run_forward,'all','omitnan'); %avg run speed
    thisTurnSpeed = mean(abs(run_angular),'all','omitnan'); %avg turn speed

    % analyze directional velocity during pulses
    [~, angMean] = pulse_v_output(int_panelps,int_forward,int_angular,int_time,1,1,nSweep,settings.runThreshB);
    [~, angMeanThresh] = pulse_v_onlyturns(int_panelps,int_forward,int_angular,int_time,1,1,nSweep,settings.runThreshB);
    pD = size(angMean.varOutR,1);

    % analyze velocity distributions
    [fwdHist,angHist,sidHist] = velocity_histogram(int_forward,int_angular,int_sideway,hNorm);

    % pool this trial data according to the assigned trial condition
    if contains(thisTrial,'KIR') %kir perturbation flies
        nKIR = nKIR+1;
        kirRunTime(nR,nKIR) = thisRunTime;
        kirRunSpeed(nR,nKIR) = thisRunSpeed;
        kirTurnSpeed(nR,nKIR) = thisTurnSpeed;

        kir_mopR(1:pD,nKIR,:) = angMean.varOutR;
        kir_mopL(1:pD,nKIR,:) = angMean.varOutL;
        kir_mopRL(1:pD,nKIR,:) = angMean.varOutRL;
        kir_mopRLthresh(1:pD,nKIR,:) = angMeanThresh.varOutRL;

        kir_mopFwdHist(:,nKIR) = fwdHist(:,2);
        kir_mopAngHist(:,nKIR) = angHist(:,2);
        kir_mopSidHist(:,nKIR) = sidHist(:,2);
    elseif contains(thisTrial,'WT') %wildtype control flies
        nWT = nWT+1;
        wtRunTime(nR,nWT) = thisRunTime;
        wtRunSpeed(nR,nWT) = thisRunSpeed;
        wtTurnSpeed(nR,nWT) = thisTurnSpeed;

        wt_mopR(1:pD,nWT,:) = angMean.varOutR;
        wt_mopL(1:pD,nWT,:) = angMean.varOutL;
        wt_mopRL(1:pD,nWT,:) = angMean.varOutRL;
        wt_mopRLthresh(1:pD,nWT,:) = angMeanThresh.varOutRL;

        wt_mopFwdHist(:,nWT) = fwdHist(:,2);
        wt_mopAngHist(:,nWT) = angHist(:,2);
        wt_mopSidHist(:,nWT) = sidHist(:,2);
    elseif contains(thisTrial,'NA') %na control flies
        nNA=nNA+1;
        naRunTime(nR,nNA) = thisRunTime;
        naRunSpeed(nR,nNA) = thisRunSpeed;
        naTurnSpeed(nR,nNA) = thisTurnSpeed;

        na_mopR(1:pD,nNA,:) = angMean.varOutR;
        na_mopL(1:pD,nNA,:) = angMean.varOutL;
        na_mopRL(1:pD,nNA,:) = angMean.varOutRL;
        na_mopRLthresh(1:pD,nNA,:) = angMeanThresh.varOutRL;

        na_mopFwdHist(:,nNA) = fwdHist(:,2);
        na_mopAngHist(:,nNA) = angHist(:,2);
        na_mopSidHist(:,nNA) = sidHist(:,2);
    end
end
% store panel data
nSweep = size(angMean.panelpsR,3);
pos_mopR = reshape(angMean.panelpsR,[],nSweep);
pos_mopL = reshape(angMean.panelpsL,[],nSweep);

% fetch time
t_mop = int_time(1:length(pos_mopL))*1000; %msec
disp('Complete.')

%% load in and pool optomotor data
disp('Loading in and analyzing optomotor datasets...')
nKIR = 0; nWT = 0; nNA = 0; nR = 3;
for e = 1:nFlies
    disp(['Processing fly ' num2str(e) '/' num2str(nFlies) '...'])
    % load this trial
    cd(folder.int)
    thisTrial = optFiles(e).name;
    load(thisTrial)

    % find timepoints where the fly was running and for how long in total
    [run_forward,run_angular,~,~] = pursuitFinder(int_forward,int_angular,int_sideway,0,int_time,settings.runThreshB);
    int_time(end+1) = 60;
    thisRunTime = sum(int_time(sum(~isnan(run_forward))+1)); %sec
    thisRunSpeed = mean(run_forward,'all','omitnan'); %avg run speed
    thisTurnSpeed = mean(abs(run_angular),'all','omitnan'); %avg turn speed

    % analyze directional velocity during optomotor task
    [mean_angularR,mean_angularL,panelWin] = optomotor_v_velocity(int_panelps, int_forward,int_angular,int_time,settings.runThreshB,0);
    mean_angularRL = mean([mean_angularR -mean_angularL],2);

    % analyze velocity distributions
    [fwdHist,angHist,sidHist] = velocity_histogram(int_forward,int_angular,int_sideway,hNorm);

    % pool this trial data according to the assigned trial condition
    if contains(thisTrial,'KIR') %kir perturbation flies
        nKIR = nKIR+1;
        kirRunTime(nR,nKIR) = thisRunTime;
        kirRunSpeed(nR,nKIR) = thisRunSpeed;
        kirTurnSpeed(nR,nKIR) = thisTurnSpeed;

        kir_optoR(:,nKIR) = mean_angularR;
        kir_optoL(:,nKIR) = mean_angularL;
        kir_optoRL(:,nKIR) = mean_angularRL;

        kir_optoFwdHist(:,nKIR) = fwdHist(:,2);
        kir_optoAngHist(:,nKIR) = angHist(:,2);
        kir_optoSidHist(:,nKIR) = sidHist(:,2);
    elseif contains(thisTrial,'WT') %wildtype control flies
        nWT = nWT+1;
        wtRunTime(nR,nWT) = thisRunTime;
        wtRunSpeed(nR,nWT) = thisRunSpeed;
        wtTurnSpeed(nR,nWT) = thisTurnSpeed;

        wt_optoR(:,nWT) = mean_angularR;
        wt_optoL(:,nWT) = mean_angularL;
        wt_optoRL(:,nWT) = mean_angularRL;

        wt_optoFwdHist(:,nWT) = fwdHist(:,2);
        wt_optoAngHist(:,nWT) = angHist(:,2);
        wt_optoSidHist(:,nWT) = sidHist(:,2);
    elseif contains(thisTrial,'NA') %na control flies
        nNA=nNA+1;
        naRunTime(nR,nNA) = thisRunTime;
        naRunSpeed(nR,nNA) = thisRunSpeed;
        naTurnSpeed(nR,nNA) = thisTurnSpeed;

        na_optoR(:,nNA) = mean_angularR;
        na_optoL(:,nNA) = mean_angularL;
        na_optoRL(:,nNA) = mean_angularRL;

        na_optoFwdHist(:,nNA) = fwdHist(:,2);
        na_optoAngHist(:,nNA) = angHist(:,2);
        na_optoSidHist(:,nNA) = sidHist(:,2);
    end
end
% store time
t_opto = int_time(1:length(mean_angularR));
% fetch opto start/stop
sweepIdx = find(~isnan(panelWin(:,1)));
sweepStartStop = [t_opto(sweepIdx(1));t_opto(sweepIdx(end))];

disp('Complete.')

%% compare basic behavior parameters
disp('Comparing basic behavior parameters...')
% initialize
nPlot = 3;
figure; set(gcf,'Position',[100 100 400 600])
tiledlayout(nPlot,3,'TileSpacing','compact')

for prm = 1:nPlot
    % for each experiment set
    for g = 1:nExpt
        % pull parameter to compare
        switch prm
            case 1 %run time
                yname = 'Run Time (min)';
                yrange = [0 15];
                thisKiR = kirRunTime(g,:)./60;
                thisWT = wtRunTime(g,:)./60;
                thisNA = naRunTime(g,:)./60;
            case 2 %mean forward speed
                yname = 'Mean Forward (mm/s)';
                yrange = [0 18];
                thisKiR = kirRunSpeed(g,:);
                thisWT = wtRunSpeed(g,:);
                thisNA = naRunSpeed(g,:);
            case 3 %mean angular speed
                yname = 'Mean Angular (deg/s)';
                yrange = [0 150];
                thisKiR = kirTurnSpeed(g,:);
                thisWT = wtTurnSpeed(g,:);
                thisNA = naTurnSpeed(g,:);
        end

        % calculate means and SEM
        thisKiRMean = mean(thisKiR);
        thisKiRSEM = std(thisKiR)./sqrt(nKIR);
        thisWTMean = mean(thisWT);
        thisWTSEM = std(thisWT)./sqrt(nWT);
        thisNAMean = mean(thisNA);
        thisNASEM = std(thisNA)./sqrt(nNA);

        % generate plot
        nexttile; hold on
        plot(1,thisKiR,'.','Color',settings.trialColor)
        errorbar(1,thisKiRMean,thisKiRSEM,'o','Color',settings.geneColor{1},'LineWidth',1)
        plot(2,thisWT,'.','Color',settings.trialColor)
        errorbar(2,thisWTMean,thisWTSEM,'o','Color',settings.geneColor{2},'LineWidth',1)
        plot(3,thisNA,'.','Color',settings.trialColor)
        errorbar(3,thisNAMean,thisNASEM,'o','Color',settings.geneColor{3},'LineWidth',1)

        axis padded
        xticks([1 2 3])
        xticklabels(settings.geneLabel)
        ylim(yrange)
        if g==1
            ylabel(yname)
        end
        if prm==1
            title(exptNames{g})
        end
    end
end
% save plot
cd(folder.summary)
plotname = 'basics_summary';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

disp('Done.')

%% oscillatory pursuit: plot directonal velocity trials+means separately
disp('Comparing pursuit behavior for oscillatory stimulus...')
% plot means and trials
dir_vel = {"fwd","ang","side"};
for e = 1:3
    % plot directional velocity separately
    figure; set(gcf,'Position',[100 100 1200 600])
    tiledlayout(1,3,'TileSpacing','compact')
    % select datasets
    switch e
        case 1 %fwd
            thisKIR = kir_oscFwd;
            thisWT = wt_oscFwd;
            thisNA = na_oscFwd;
            yname = "Forward Velocity (mm/s)";
            yrange = [0 18];
        case 2 %ang
            thisKIR = kir_oscAng;
            thisWT = wt_oscAng;
            thisNA = na_oscAng;
            yname = "Angular Velocity (deg/s)";
            yrange = [-300 300];
        case 3 %sid
            thisKIR = kir_oscSid;
            thisWT = wt_oscSid;
            thisNA = na_oscSid;
            yname = "Sideways Velocity (mm/s)";
            yrange = [-9.5 9.5];
    end
    % calculate mean for each genotype
    kirMean = mean(thisKIR,2);
    wtMean = mean(thisWT,2);
    naMean = mean(thisNA,2);
    
    % generate plot
    nexttile; hold on
    % add target position reference
    fc = gca;
    yyaxis right
    plot(t_osc,oscSweepPos,'k')
    axis tight
    fc.YAxis(1).Color = 'k';
    fc.YAxis(2).Visible = 'off';
    % add reference lines
    yline(0,'Color','k')
    yyaxis left
    plot(t_osc,thisKIR,'-','Color',settings.trialColor)
    plot(t_osc,kirMean,'-','Color', settings.geneColor{1},'LineWidth',1.5)
    axis tight
    ylim(yrange)
    yline(0)
    ylabel(yname)
    title(settings.geneLabel{1})

    nexttile; hold on
    % add target position reference
    fc = gca;
    yyaxis right
    plot(t_osc,oscSweepPos,'k')
    axis tight
    fc.YAxis(1).Color = 'k';
    fc.YAxis(2).Visible = 'off';
    % add reference lines
    yline(0,'Color','k')
    yyaxis left
    plot(t_osc,thisWT,'-','Color',settings.trialColor)
    plot(t_osc,wtMean,'-','Color', settings.geneColor{2},'LineWidth',1.5)
    axis tight
    ylim(yrange)
    yline(0)
    title(settings.geneLabel{2})

    nexttile; hold on
    % add target position reference
    fc = gca;
    yyaxis right
    plot(t_osc,oscSweepPos,'k')
    axis tight
    fc.YAxis(1).Color = 'k';
    fc.YAxis(2).Visible = 'off';
    % add reference lines
    yline(0,'Color','k')
    yyaxis left
    plot(t_osc,thisNA,'-','Color',settings.trialColor)
    plot(t_osc,naMean,'-','Color', settings.geneColor{3},'LineWidth',1.5)
    axis tight
    ylim(yrange)
    yline(0)
    title(settings.geneLabel{3})
    xlabel('Time (msec.)')

    % save plot
    cd(folder.summary)
    plotname = ['oscpursuit_' char(dir_vel{e}) '_trials'];
    saveas(gcf,[plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox,'f');
end


%% oscillatory pursuit: plot directonal velocity means together
disp('Comparing pursuit behavior for oscillatory stimulus...')

% initialize
figure; set(gcf,'Position',[100 100 1200 600])
nMax = max([nKIR nWT nNA]);
maxFwd = nan(nMax,3);
maxAng = nan(nMax,3);
maxSid = nan(nMax,3);

limit_fwd = [0 18];
limit_ang = [-175 175];
limit_sid = [-6 6];

for e = 1:3
    % select dataset
    switch e
        case 1 %kir
            thisFwd = kir_oscFwd;
            thisAng = kir_oscAng;
            thisSid = kir_oscSid;
            thisN = nKIR;
        case 2 %wt
            thisFwd = wt_oscFwd;
            thisAng = wt_oscAng;
            thisSid = wt_oscSid;
            thisN = nWT;
        case 3 %na
            thisFwd = na_oscFwd;
            thisAng = na_oscAng;
            thisSid = na_oscSid;
            thisN = nNA;
    end
    % calculate mean for each directional velocity
    fwdMean = mean(thisFwd,2);
    angMean = mean(thisAng,2);
    sidMean = mean(thisSid,2);
    fwdSEM = std(thisFwd,[],2)./sqrt(thisN);
    angSEM = std(thisAng,[],2)./sqrt(thisN);
    sidSEM = std(thisSid,[],2)./sqrt(thisN);

    % store max for plotting later
    maxFwd(1:thisN,e) = max(thisFwd)';
    maxAng(1:thisN,e) = max(abs(thisAng))';
    maxSid(1:thisN,e) = max(abs(thisSid))';

    % generate plot
    subplot(1,3,1); hold on
    if e==1
        % add target position reference
        fc = gca;
        yyaxis right
        plot(t_osc,oscSweepPos,'k')
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
    end
    yyaxis left
    plot(t_osc,fwdMean,'-','Color', settings.geneColor{e},'LineWidth',1.5)
    sp1 = patch([t_osc'; flipud(t_osc')],[fwdMean-fwdSEM; flipud(fwdMean+fwdSEM)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{e};
    ylabel('Forward Velocity (mm/s)')
    axis tight
    ylim(limit_fwd)

    subplot(1,3,2); hold on
    if e==1
        % add target position reference
        fc = gca;
        yyaxis right
        plot(t_osc,oscSweepPos,'k')
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
    end
    yyaxis left
    plot(t_osc,angMean,'-','Color', settings.geneColor{e},'LineWidth',1.5)
    sp1 = patch([t_osc'; flipud(t_osc')],[angMean-angSEM; flipud(angMean+angSEM)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{e};
    ylabel('Angular Velocity (deg/s)')
    axis tight
    ylim(limit_ang)

    subplot(1,3,3); hold on
    if e==1
        % add target position reference
        fc = gca;
        yyaxis right
        plot(t_osc,oscSweepPos,'k')
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
    end
    yyaxis left
    plot(t_osc,sidMean,'-','Color', settings.geneColor{e},'LineWidth',1.5)
    sp1 = patch([t_osc'; flipud(t_osc')],[sidMean-sidSEM; flipud(sidMean+sidSEM)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{e};
    ylabel('Sideways Velocity (mm/s)')
    axis tight
    xlabel('Time (msec.)')
    ylim(limit_sid)
end
% save plot
sgtitle('Velocity Comparison')
cd(folder.summary)
plotname = 'oscpursuit_summary';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

%% oscillatory pursuit: plot directonal velocity means together R+L
disp('Comparing R+L pursuit behavior for oscillatory stimulus...')

% initialize
figure; set(gcf,'Position',[100 100 1200 600])
nMax = max([nKIR nWT nNA]);
maxFwd = nan(nMax,3);
maxAng = nan(nMax,3);
maxSid = nan(nMax,3);

limit_fwd = [0 18];
limit_ang = [-175 175];
limit_sid = [-6 6];

% Get middle 50% indices
len = length(t_osc);
q1 = round(len * 0.25);
q3 = round(len * 0.75);
midIdx = q1:q3;

for e = 1:3
    % select dataset
    switch e
        case 1 %kir
            thisFwd = kir_oscFwdRL;
            thisAng = kir_oscAngRL;
            thisSid = kir_oscSidRL;
            thisN = nKIR;
        case 2 %wt
            thisFwd = wt_oscFwdRL;
            thisAng = wt_oscAngRL;
            thisSid = wt_oscSidRL;
            thisN = nWT;
        case 3 %na
            thisFwd = na_oscFwdRL;
            thisAng = na_oscAngRL;
            thisSid = na_oscSidRL;
            thisN = nNA;
    end
    % calculate mean for each directional velocity
    fwdMean = mean(thisFwd,2);
    angMean = mean(thisAng,2);
    sidMean = mean(thisSid,2);
    fwdSEM = std(thisFwd,[],2)./sqrt(thisN);
    angSEM = std(thisAng,[],2)./sqrt(thisN);
    sidSEM = std(thisSid,[],2)./sqrt(thisN);

    % store max for plotting later
    maxFwd(1:thisN,e) = max(thisFwd)';
    maxAng(1:thisN,e) = max(abs(thisAng))';
    maxSid(1:thisN,e) = max(abs(thisSid))';

    % generate plot
    subplot(1,3,1); hold on
    if e==1
        % add target position reference
        fc = gca;
        yyaxis right
        plot(t_osc,oscSweepPos,'k')
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
    end
    yyaxis left
    plot(t_osc,fwdMean,'-','Color', settings.geneColor{e},'LineWidth',1.5)
    sp1 = patch([t_osc'; flipud(t_osc')],[fwdMean-fwdSEM; flipud(fwdMean+fwdSEM)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{e};
    ylabel('Forward Velocity (mm/s)')
    axis tight
    ylim(limit_fwd)
    xlim([t_osc(q1) t_osc(q3)])

    subplot(1,3,2); hold on
    if e==1
        % add target position reference
        fc = gca;
        yyaxis right
        plot(t_osc,oscSweepPos,'k')
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
    end
    yyaxis left
    plot(t_osc,angMean,'-','Color', settings.geneColor{e},'LineWidth',1.5)
    sp1 = patch([t_osc'; flipud(t_osc')],[angMean-angSEM; flipud(angMean+angSEM)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{e};
    ylabel('Angular Velocity (deg/s)')
    axis tight
    xlim([t_osc(q1) t_osc(q3)])
    ylim(limit_ang)

    subplot(1,3,3); hold on
    if e==1
        % add target position reference
        fc = gca;
        yyaxis right
        plot(t_osc,oscSweepPos,'k')
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
    end
    yyaxis left
    plot(t_osc,sidMean,'-','Color', settings.geneColor{e},'LineWidth',1.5)
    sp1 = patch([t_osc'; flipud(t_osc')],[sidMean-sidSEM; flipud(sidMean+sidSEM)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{e};
    ylabel('Sideways Velocity (mm/s)')
    axis tight
    xlabel('Time (msec.)')
    xlim([t_osc(q1) t_osc(q3)])
    ylim(limit_sid)
end
% save plot
sgtitle('Velocity Comparison (R+L)')
cd(folder.summary)
plotname = 'oscpursuit_RLsummary';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

% Fit oscillatory data
plotOscSlopes(kir_oscAngRL, wt_oscAngRL, na_oscAngRL, oscSweepPos,folder)

cd(folder.summary)
plotname = 'oscpursuit_RLfit';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

run_forward_velocity_anova(oscSweepPos, kir_oscFwdRL, wt_oscFwdRL, na_oscFwdRL)
cd(folder.summary)
plotname = 'oscpursuit_changefwd';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

%% oscillatory pursuit: plot directonal acceleration means together
disp('Comparing pursuit behavior for oscillatory stimulus...')

% initialize
figure; set(gcf,'Position',[100 100 1200 600])
nMax = max([nKIR nWT nNA]);
maxFwd = nan(nMax,3);
maxAng = nan(nMax,3);
maxSid = nan(nMax,3);

limit_fwd = [-15 15];
limit_ang = [-600 600];
limit_sid = [-25 25];

for e = 1:3
    % select dataset
    switch e
        case 1 %kir
            thisFwd = kir_oscAccFwdRL;
            thisAng = kir_oscAccAngRL;
            thisSid = kir_oscAccSidRL;
            thisN = nKIR;
        case 2 %wt
            thisFwd = wt_oscAccFwdRL;
            thisAng = wt_oscAccAngRL;
            thisSid = wt_oscAccSidRL;
            thisN = nWT;
        case 3 %na
            thisFwd = na_oscAccFwdRL;
            thisAng = na_oscAccAngRL;
            thisSid = na_oscAccSidRL;
            thisN = nNA;
    end
    % calculate mean for each directional velocity
    fwdMean = mean(thisFwd,2);
    angMean = mean(thisAng,2);
    sidMean = mean(thisSid,2);
    fwdSEM = std(thisFwd,[],2)./sqrt(thisN);
    angSEM = std(thisAng,[],2)./sqrt(thisN);
    sidSEM = std(thisSid,[],2)./sqrt(thisN);

    % store max for plotting later
    maxFwd(1:thisN,e) = max(thisFwd)';
    maxAng(1:thisN,e) = max(abs(thisAng))';
    maxSid(1:thisN,e) = max(abs(thisSid))';

    % generate plot
    subplot(1,3,1); hold on
    if e==1
        % add target position reference
        fc = gca;
        yyaxis right
        plot(t_osc,oscSweepPos,'k')
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
    end
    yyaxis left
    plot(t_osc,fwdMean,'-','Color', settings.geneColor{e},'LineWidth',1.5)
    sp1 = patch([t_osc'; flipud(t_osc')],[fwdMean-fwdSEM; flipud(fwdMean+fwdSEM)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{e};
    ylabel('Forward Acceleration (mm/s^2)')
    axis tight
    ylim(limit_fwd)

    subplot(1,3,2); hold on
    if e==1
        % add target position reference
        fc = gca;
        yyaxis right
        plot(t_osc,oscSweepPos,'k')
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
    end
    yyaxis left
    plot(t_osc,angMean,'-','Color', settings.geneColor{e},'LineWidth',1.5)
    sp1 = patch([t_osc'; flipud(t_osc')],[angMean-angSEM; flipud(angMean+angSEM)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{e};
    ylabel('Angular Acceleration (deg/s^2)')
    axis tight
    ylim(limit_ang)

    subplot(1,3,3); hold on
    if e==1
        % add target position reference
        fc = gca;
        yyaxis right
        plot(t_osc,oscSweepPos,'k')
        axis tight
        fc.YAxis(1).Color = 'k';
        fc.YAxis(2).Visible = 'off';
        % add reference lines
        yline(0,'Color','k')
    end
    yyaxis left
    plot(t_osc,sidMean,'-','Color', settings.geneColor{e},'LineWidth',1.5)
    sp1 = patch([t_osc'; flipud(t_osc')],[sidMean-sidSEM; flipud(sidMean+sidSEM)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{e};
    ylabel('Sideways Acceleration (mm/s^2)')
    axis tight
    xlabel('Time (msec.)')
    ylim(limit_sid)
end
% save plot
sgtitle('Acceleration Comparison (R+L)')
cd(folder.summary)
plotname = 'oscpursuit_accel_summary';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');


%% oscillatory pursuit: summary stats
disp('Comparing behavior stats for oscillatory stimulus...')
nPlot = 3;
figure; set(gcf,'Position',[100 100 600 300])
tiledlayout(1,nPlot,'TileSpacing','compact')
for prm = 1:nPlot
    switch prm
        case 1 %fwd max
            yname = 'Max Forward Speed (mm/s)';
            yrange = [0 40];
            thisData = maxFwd;
        case 2 %angular max
            yname = 'Max Angular Speed (deg/s)';
            yrange = [0 800];
            thisData = maxAng;
        case 3 %sideways max
            yname = 'Max Sideways Speed (mm/s)';
            yrange = [0 40];
            thisData = maxSid;
    end

    % calculate means and SEM
    thisMean = mean(thisData,1,'omitnan');
    thisSEM = std(thisData,1,'omitnan')./sqrt(sum(~isnan(thisData)));

    % generate plot
    nexttile; hold on
    for x = 1:3
        plot(x,thisData(:,x),'.','Color',settings.trialColor)
        errorbar(x,thisMean(x),thisSEM(x),'o','Color',settings.geneColor{x},'LineWidth',1)
    end
    axis padded
    xticks([1 2 3])
    xticklabels(settings.geneLabel)
    ylim(yrange)
    ylabel(yname)
end

% save plot
cd(folder.summary)
plotname = 'basics_oscpursuit';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

%% oscillatory pursuit: histogram distribution
figure; set(gcf,'Position',[100 100 1800 800])
tiledlayout(3,4,'TileSpacing','compact')
fBin = fwdHist(:,1);
aBin = angHist(:,1);
sBin = sidHist(:,1);

for v = 1:3
    switch v
        case 1 %forward
            hBin = fBin;
            kirVel = kir_oscFwdHist;
            wtVel = wt_oscFwdHist;
            naVel = na_oscFwdHist;
            xvar = 'Binned Forward Velocity (mm/s)';
        case 2 %angular
            hBin = aBin;
            kirVel = kir_oscAngHist;
            wtVel = wt_oscAngHist;
            naVel = na_oscAngHist;
            xvar = 'Binned Angular Velocity (deg/s)';
        case 3 %sideways
            hBin = sBin;
            kirVel = kir_oscSidHist;
            wtVel = wt_oscSidHist;
            naVel = na_oscSidHist;
            xvar = 'Binned Sideways Velocity (mm/s)';
    end
    % calculate mean and sem
    kirVel_mean = mean(kirVel,2);
    wtVel_mean = mean(wtVel,2);
    naVel_mean = mean(naVel,2);
    kirVel_sem = std(kirVel,[],2)./sqrt(nKIR);
    wtVel_sem = std(wtVel,[],2)./sqrt(nKIR);
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
sgtitle('Oscillatory Pursuit - Velocity Distribution')
cd(folder.summary)
plotname = 'hist_oscpursuit';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

%% motion pulse: plot directional velocity trials+means separately
disp('Analyzing pursuit behavior during motion pulse experiment...')
% variables
ang_range = [-350 350]; %velocity range
tmax = max(t_mop); %time range

% for each genotype condition
for g = 1:3
    % fetch corresponding data
    switch g
        case 1 %kir
            thisG = settings.geneLabel{g};
            mopR = kir_mopR;
            mopL = kir_mopL;
        case 2 %wt
            thisG = settings.geneLabel{g};
            mopR = wt_mopR;
            mopL = wt_mopL;
        case 3 %na
            thisG = settings.geneLabel{g};
            mopR = na_mopR;
            mopL = na_mopL;
    end

    % for each sweep direction (left/right)
    for d = 1:2
        % initialize
        figure; set(gcf,'Position',[100 100 1800 900])
        tiledlayout(2,nSweep,'TileSpacing','compact')

        % plot panel position
        for ns = 1:nSweep
            nexttile
            if d==1 %right
                plot(t_mop,pos_mopR(:,ns),'-','Color','k','Linewidth',1.5)
            else %left
                plot(t_mop,pos_mopL(:,ns),'-','Color','k','Linewidth',1.5)
            end
            xlim([0 tmax])
            ylim([-160 160])
            yline(0)
            if ns==1
                ylabel('Object Position (deg)')
            end
        end
        % plot angular velocity
        for ns = 1:nSweep
            % fetch data for this sweep
            if d==1
                thisTitle = "Rightward Sweeps";
                thisSave = "R";
                thisMOP = mopR(:,:,ns);
            else
                thisTitle = "Leftward Sweeps";
                thisSave = "L";
                thisMOP = mopL(:,:,ns);
            end

            % calculate mean and SEM
            meanMOP = mean(thisMOP,2,'omitnan');

            % generate plot
            nexttile; hold on
            plot(t_mop,thisMOP,'-','Color',settings.trialColor)
            plot(t_mop,meanMOP,'-','Color',settings.geneColor{g},'Linewidth',1.5)
            xlim([0 tmax])
            ylim(ang_range)
            yline(0)
            if ns==1
                ylabel('Angular Velocity (deg/s)')
            end
        end
        sgtitle([thisG thisTitle])
        xlabel('Time (msec.)')

        % save plot
        cd(folder.summary)
        plotname = append(thisG,'_mopulse_', thisSave,'_trials');
        saveas(gcf,append(plotname,'.png'));
        copyfile(append(plotname,'.png'), folder.dropbox,'f');
        % save vectorized plot
        cd(folder.vector)
        set(gcf,'renderer','Painters')
        saveas(gcf, append(plotname,'.svg'))
        copyfile(append(plotname,'.svg'), folder.dropbox,'f');

    end
end


%% motion pulse: plot directional velocity means together
disp('Analyzing pursuit behavior during motion pulse experiment...')
% variables
ang_range = [-250 250]; %velocity range
tmax = max(t_mop); %time range
idxSweepStart = find(~isnan(pos_mopR(:,1)),1,'first');
peakKIR = []; peakWT = []; peakNA = [];
peakControls =[];

% for each sweep direction (left/right)
for d = 1:2
    % initialize
    figure; set(gcf,'Position',[100 100 1800 900])
    tiledlayout(2,nSweep,'TileSpacing','compact')

    % plot panel position
    for ns = 1:nSweep
        nexttile
        if d==1 %right
            plot(t_mop,pos_mopR(:,ns),'-','Color','k','Linewidth',1.5)
        else %left
            plot(t_mop,pos_mopL(:,ns),'-','Color','k','Linewidth',1.5)
        end
        xlim([0 tmax])
        ylim([-160 160])
        yline(0)
        if ns==1
            ylabel('Object Position (deg)')
        end
    end
    % plot angular velocity
    for ns = 1:nSweep
        % fetch data for this sweep
        if d==1
            thisTitle = "Rightward Sweeps";
            thisSave = "R";
            pk = 6;
            thisKIR = kir_mopR(:,:,ns);
            thisWT = wt_mopR(:,:,ns);
            thisNA = na_mopR(:,:,ns);
        else
            thisTitle = "Leftward Sweeps";
            thisSave = "L";
            pk = 9;
            thisKIR = kir_mopL(:,:,ns);
            thisWT = wt_mopL(:,:,ns);
            thisNA = na_mopL(:,:,ns);
        end

        % calculate mean and SEM
        meanKIR = mean(thisKIR,2,'omitnan');
        meanWT = mean(thisWT,2,'omitnan');
        meanNA = mean(thisNA,2,'omitnan');
        semKIR = std(thisKIR,[],2,'omitnan')./sqrt(nKIR);
        semWT = std(thisWT,[],2,'omitnan')./sqrt(nWT);
        semNA = std(thisNA,[],2,'omitnan')./sqrt(nNA);

        % pull turn peak from pulse onset to window end
        % pull left for sweeps on left and right for sweeps on right
        if ns<pk
            peakKIR(ns,:) = min(thisKIR(idxSweepStart:end,:));
            peakWT(ns,:) = min(thisWT(idxSweepStart:end,:));
            peakNA(ns,:) = min(thisNA(idxSweepStart:end,:));
        else
            peakKIR(ns,:) = max(thisKIR(idxSweepStart:end,:));
            peakWT(ns,:) = max(thisWT(idxSweepStart:end,:));
            peakNA(ns,:) = max(thisNA(idxSweepStart:end,:));
        end
        

        % generate plot
        nexttile; hold on
        plot(t_mop,meanKIR,'-','Color',settings.geneColor{1},'Linewidth',1.5)
        sp1 = patch([t_mop'; flipud(t_mop')],[meanKIR-semKIR; flipud(meanKIR+semKIR)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
        sp1.FaceColor = settings.geneColor{1};
        plot(t_mop,meanWT,'-','Color',settings.geneColor{2},'Linewidth',1.5)
        sp2 = patch([t_mop'; flipud(t_mop')],[meanWT-semWT; flipud(meanWT+semWT)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
        sp2.FaceColor = settings.geneColor{2};
        plot(t_mop,meanNA,'-','Color',settings.geneColor{3},'Linewidth',1.5)
        sp3 = patch([t_mop'; flipud(t_mop')],[meanNA-semNA; flipud(meanNA+semNA)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
        sp3.FaceColor = settings.geneColor{3};
        xlim([0 tmax])
        ylim(ang_range)
        yline(0)
        if ns==1
            ylabel('Angular Velocity (deg/s)')
        end
    end
    sgtitle(thisTitle)
    xlabel('Time (msec.)')

    % save plot
    cd(folder.summary)
    plotname = append('mopulse_', thisSave,'_summary');
    saveas(gcf,append(plotname,'.png'));
    copyfile(append(plotname,'.png'), folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, append(plotname,'.svg'))
    copyfile(append(plotname,'.svg'), folder.dropbox,'f');


    % plot a quick summary
    figure; set(gcf,'Position',[100 100 900 900]); hold on
    x = 1:nSweep;
    % plot trial data
    plot(x+offset(1),peakKIR,'.','Color',settings.trialColor)
    plot(x+offset(2),peakWT,'.','Color',settings.trialColor)
    plot(x+offset(3),peakNA,'.','Color',settings.trialColor)
    % plot mean +/- SEM
    errorbar(x+offset(1),mean(peakKIR,2,'omitnan'),std(peakKIR,[],2,'omitnan')./sqrt(nKIR),'o','Color',settings.geneColor{1},'Linewidth',1)
    errorbar(x+offset(2),mean(peakWT,2,'omitnan'),std(peakWT,[],2,'omitnan')./sqrt(nWT),'o','Color',settings.geneColor{2},'Linewidth',1)
    errorbar(x+offset(3),mean(peakNA,2,'omitnan'),std(peakNA,[],2,'omitnan')./sqrt(nNA),'o','Color',settings.geneColor{3},'Linewidth',1)
    axis padded
    yline(0);
    ylim([-350 350])
    ylabel('Peak Angular Velocity (deg/s)')
    xticks(x)
    sgtitle(thisTitle)
    xlabel('Sweep Order')

    % save plot
    cd(folder.summary)
    plotname = append('mopulse_', thisSave,'_plot');
    saveas(gcf,append(plotname,'.png'));
    copyfile(append(plotname,'.png'), folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, append(plotname,'.svg'))
    copyfile(append(plotname,'.svg'), folder.dropbox,'f');

    % store controls together
    if d==1 %right
        peakControls(:,:,d) = [peakWT(:,:) peakNA(:,:)];
    elseif d==2 %left
        peakControls(:,:,d) = [peakWT(:,:) peakNA(:,:)];
    end

end

%% Prepare data for repeated-measures ANOVA using fitlme
peakControlsRL = [];
peakControlsRL(:,:,1) = peakControls(:,:,1) - flip(peakControls(:,:,2),1);
peakControlsRL(:,:,2) = peakControls(:,:,2) - flip(peakControls(:,:,1),1);
posR = [-123.75, -101.25, -78.75, -56.25, -33.75, -11.25, 11.25, 33.75, 56.25, 78.75, 101.25, 123.75, 146.25];
posL = [-146.25, -123.75, -101.25, -78.75, -56.25, -33.75, -11.25, 11.25, 33.75, 56.25, 78.75, 101.25, 123.75];

% Remove columns (flies) containing NaNs
validCols = ~any(any(isnan(peakControls), 1), 3); % Find valid flies (no NaNs across all rows and directions)
filteredControls = peakControlsRL(:, validCols, :);

%% Define overlapping sweep positions
% These are the positions that exist in both posR and posL
overlapPos = intersect(posR, posL); % Gives symmetric central range

% Find indices of overlapping positions for both directions
[~, idxR] = ismember(overlapPos, posR);
[~, idxL] = ismember(overlapPos, posL);

% Filter data to just overlapping sweep positions
filteredControls_overlap = filteredControls(idxR, :, :); % Apply to both directions

% Calculate mean and SEM for plotting
meanRight = mean(filteredControls_overlap(:,:,1), 2, 'omitnan');
semRight = std(filteredControls_overlap(:,:,1), 0, 2, 'omitnan') ./ sqrt(sum(~isnan(filteredControls_overlap(:,:,1)), 2));

meanLeft = mean(filteredControls_overlap(:,:,2), 2, 'omitnan');
semLeft = std(filteredControls_overlap(:,:,2), 0, 2, 'omitnan') ./ sqrt(sum(~isnan(filteredControls_overlap(:,:,2)), 2));

% Prepare data for repeated-measures ANOVA using fitlme

% Dimensions: [positions x animals x directions]
numPositions = size(filteredControls_overlap, 1);
numAnimals   = size(filteredControls_overlap, 2);

% Flatten data
data = filteredControls_overlap(:);

% Build label vectors
[positions, animals, directions] = ndgrid(1:numPositions, 1:numAnimals, 1:2);
SweepPos = categorical(positions(:));
Direction = categorical(directions(:));
FlyID = categorical(animals(:));

% Combine into table
T = table(data, SweepPos, Direction, FlyID, ...
    'VariableNames', {'Turning', 'SweepPos', 'Direction', 'FlyID'});

% Remove NaNs
T = T(~isnan(T.Turning), :);

% Fit LME with interaction
lme = fitlme(T, 'Turning ~ SweepPos*Direction + (1|FlyID)');

% Extract ANOVA p-values
a = anova(lme);
pSweep = a.pValue(strcmp(a.Term, 'SweepPos'));
pDir   = a.pValue(strcmp(a.Term, 'Direction'));
pInt   = a.pValue(strcmp(a.Term, 'SweepPos:Direction'));

% Plotting

figure;
hold on;

% X-axis values: overlapping physical positions
x_overlap = overlapPos;

% Plot rightward and leftward sweeps
errorbar(x_overlap, meanRight, semRight, '-', 'Color', 'b', 'CapSize', 0, 'LineWidth', 1, 'DisplayName', 'Rightward sweeps');
errorbar(x_overlap, meanLeft, semLeft, '-', 'Color', 'k', 'CapSize', 0, 'LineWidth', 1, 'DisplayName', 'Leftward sweeps');

% Plot formatting
xlabel('Sweep Position');
ylabel('Peak Turning Response');
title(['Peak Turning Response by Sweep Position (n = ' num2str(numAnimals) ')']);
legend('Location', 'best');
xticks([-150:30:150]);
axis padded;
grid on;
yline(0);

% Add p-values to bottom-right of plot
text(0.75, 0.1, sprintf('p_{Sweep} = %.3g', pSweep), 'Units', 'normalized', 'HorizontalAlignment', 'right');
text(0.75, 0.05, sprintf('p_{Direction} = %.3g', pDir), 'Units', 'normalized', 'HorizontalAlignment', 'right');
text(0.75, 0.00, sprintf('p_{Interaction} = %.3g', pInt), 'Units', 'normalized', 'HorizontalAlignment', 'right');

hold off;

% Save plot
cd(folder.summary)
plotname = append('mopulse_motiondir_plot');
saveas(gcf, append(plotname, '.png'));
copyfile(append(plotname, '.png'), folder.dropbox, 'f');

cd(folder.vector)
set(gcf, 'renderer', 'Painters');
saveas(gcf, append(plotname, '.svg'));
copyfile(append(plotname, '.svg'), folder.dropbox, 'f');



%% motion pulse: plot directional velocity means together pool L+R
disp('Analyzing pursuit behavior during motion pulse experiment...')
% variables
ang_range = [-250 250]; %velocity range
tmax = max(t_mop); %time range
idxSweepStart = find(~isnan(pos_mopR(:,1)),1,'first');
firstSweep = 1; %first sweep to plot
sweepSelect = firstSweep:nSweep;

% initialize
figure; set(gcf,'Position',[100 100 1800 900])
tiledlayout(2,length(sweepSelect),'TileSpacing','compact')
peakKIR = []; peakWT = []; peakNA = [];
peakIdxKIR = []; peakIdxWT = []; peakIdxGFP = [];

% plot panel position
for ns = sweepSelect
    nexttile
    plot(t_mop,pos_mopR(:,ns),'-','Color','k','Linewidth',1.5)
    xlim([0 tmax])
    ylim([-160 160])
    yline(0)
    if ns==firstSweep
        ylabel('Object Position (deg)')
    end
end

% plot angular velocity
s = 0;
for ns = sweepSelect
    % fetch data for this sweep
    thisTitle = "Combined RL Sweeps";
    thisSave = "RL";
    thisKIR = kir_mopRL(:,:,ns);
    thisWT = wt_mopRL(:,:,ns);
    thisNA = na_mopRL(:,:,ns);

    % calculate mean and SEM
    meanKIR = mean(thisKIR,2,'omitnan');
    meanWT = mean(thisWT,2,'omitnan');
    meanNA = mean(thisNA,2,'omitnan');
    semKIR = std(thisKIR,[],2,'omitnan')./sqrt(nKIR);
    semWT = std(thisWT,[],2,'omitnan')./sqrt(nWT);
    semNA = std(thisNA,[],2,'omitnan')./sqrt(nNA);

    % pull turn peak from pulse onset to window end
    s = s+1;
    if s>=floor(nSweep/2)
        [peakKIR(s,:),peakIdxKIR(s,:)] = max(thisKIR(idxSweepStart:end,:));
        [peakWT(s,:),peakIdxWT(s,:)] = max(thisWT(idxSweepStart:end,:));
        [peakNA(s,:),peakIdxGFP(s,:)] = max(thisNA(idxSweepStart:end,:));
    else
        [peakKIR(s,:),peakIdxKIR(s,:)] = min(thisKIR(idxSweepStart:end,:));
        [peakWT(s,:),peakIdxWT(s,:)] = min(thisWT(idxSweepStart:end,:));
        [peakNA(s,:),peakIdxGFP(s,:)] = min(thisNA(idxSweepStart:end,:));
    end

    % generate plot
    nexttile; hold on
    plot(t_mop,meanKIR,'-','Color',settings.geneColor{1},'Linewidth',1.5)
    sp1 = patch([t_mop'; flipud(t_mop')],[meanKIR-semKIR; flipud(meanKIR+semKIR)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{1};
    plot(t_mop,meanWT,'-','Color',settings.geneColor{2},'Linewidth',1.5)
    sp2 = patch([t_mop'; flipud(t_mop')],[meanWT-semWT; flipud(meanWT+semWT)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp2.FaceColor = settings.geneColor{2};
    plot(t_mop,meanNA,'-','Color',settings.geneColor{3},'Linewidth',1.5)
    sp3 = patch([t_mop'; flipud(t_mop')],[meanNA-semNA; flipud(meanNA+semNA)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp3.FaceColor = settings.geneColor{3};
    xlim([0 tmax])
    ylim(ang_range)
    yline(0)
    if ns==firstSweep
        ylabel('Angular Velocity (deg/s)')
    end
end
sgtitle(thisTitle)
xlabel('Time (msec.)')

% save plot
cd(folder.summary)
plotname = append('mopulse_', thisSave,'_summary');
saveas(gcf,append(plotname,'.png'));
copyfile(append(plotname,'.png'), folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, append(plotname,'.svg'))
copyfile(append(plotname,'.svg'), folder.dropbox,'f');


% plot a quick summary
figure; set(gcf,'Position',[100 100 1200 900]);tiledlayout(1,2)
% plot peak angular velocity
nexttile; hold on
% plot trial data
plot(sweepSelect+offset(1),peakKIR,'.','Color',settings.trialColor)
plot(sweepSelect+offset(2),peakWT,'.','Color',settings.trialColor)
plot(sweepSelect+offset(3),peakNA,'.','Color',settings.trialColor)
% plot mean +/- SEM
errorbar(sweepSelect+offset(1),mean(peakKIR,2,'omitnan'),std(peakKIR,[],2,'omitnan')./sqrt(nKIR),'o','Color',settings.geneColor{1},'Linewidth',1)
errorbar(sweepSelect+offset(2),mean(peakWT,2,'omitnan'),std(peakWT,[],2,'omitnan')./sqrt(nWT),'o','Color',settings.geneColor{2},'Linewidth',1)
errorbar(sweepSelect+offset(3),mean(peakNA,2,'omitnan'),std(peakNA,[],2,'omitnan')./sqrt(nNA),'o','Color',settings.geneColor{3},'Linewidth',1)
axis padded
xticks(sweepSelect); xline(5.5); yline(0); ylim([-350 350])
xlabel('Sweep Order'); ylabel('Peak Angular Velocity (deg/s)'); sgtitle(thisTitle)

% plot time of peak angular velocity
% convert from idx to time
peaktKIR = t_mop(peakIdxKIR);
peaktWT = t_mop(peakIdxWT);
peaktGFP = t_mop(peakIdxGFP);
% remove empties
peaktKIR(peaktKIR==0) = nan;
peaktWT(peaktWT==0) = nan;
peaktGFP(peaktGFP==0) = nan;

nexttile; hold on
% plot trial data
plot(sweepSelect+offset(1),peaktKIR,'.','Color',settings.trialColor)
plot(sweepSelect+offset(2),peaktWT,'.','Color',settings.trialColor)
plot(sweepSelect+offset(3),peaktGFP,'.','Color',settings.trialColor)
% plot mean +/- SEM
errorbar(sweepSelect+offset(1),mean(peaktKIR,2,'omitnan'),std(peaktKIR,[],2,'omitnan')./sqrt(nKIR),'o','Color',settings.geneColor{1},'Linewidth',1)
errorbar(sweepSelect+offset(2),mean(peaktWT,2,'omitnan'),std(peaktWT,[],2,'omitnan')./sqrt(nWT),'o','Color',settings.geneColor{2},'Linewidth',1)
errorbar(sweepSelect+offset(3),mean(peaktGFP,2,'omitnan'),std(peaktGFP,[],2,'omitnan')./sqrt(nNA),'o','Color',settings.geneColor{3},'Linewidth',1)
axis padded
xticks(sweepSelect); xline(5.5); yline(0); %ylim([-120 120])
xlabel('Sweep Order'); ylabel('Peak Time (msec)'); sgtitle(thisTitle)

% save plot
cd(folder.summary)
plotname = append('mopulse_', thisSave,'_plot');
saveas(gcf,append(plotname,'.png'));
copyfile(append(plotname,'.png'), folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, append(plotname,'.svg'))
copyfile(append(plotname,'.svg'), folder.dropbox,'f');

% Save turn response for behavior
cd(folder.compare);
saveFilename = 'HeadClosedBehavior.mat';
save(saveFilename, 'peakKIR', 'peakWT', 'peakNA');

%% motion pulse: plot directional velocity means together pool L+R
disp('Analyzing pursuit behavior during motion pulse experiment...')
% variables
ang_range = [-250 250]; %velocity range
tmax = max(t_mop); %time range
idxSweepStart = find(~isnan(pos_mopR(:,1)),1,'first');
firstSweep = 1; %first sweep to plot
sweepSelect = firstSweep:nSweep;

% initialize
figure; set(gcf,'Position',[100 100 1800 900])
tiledlayout(2,length(sweepSelect),'TileSpacing','compact')
peakKIR = []; peakCTRL = [];
peakIdxKIR = []; peakIdxCTRL = [];

% plot panel position
for ns = sweepSelect
    nexttile
    plot(t_mop,pos_mopR(:,ns),'-','Color','k','Linewidth',1.5)
    xlim([0 tmax])
    ylim([-160 160])
    yline(0)
    if ns==firstSweep
        ylabel('Object Position (deg)')
    end
end

% plot angular velocity
s = 0;
for ns = sweepSelect
    % fetch data for this sweep
    thisTitle = "Combined RL Sweeps";
    thisSave = "RL";
    thisKIR = kir_mopRL(:,:,ns);
    thisWT = wt_mopRL(:,:,ns);
    thisNA = na_mopRL(:,:,ns);

    % combine WT and NA data
    thisCTRL = cat(2, thisWT, thisNA);

    % calculate mean and SEM
    meanKIR = mean(thisKIR,2,'omitnan');
    semKIR = std(thisKIR,[],2,'omitnan')./sqrt(nKIR);
    meanCTRL = mean(thisCTRL,2,'omitnan');
    semCTRL = std(thisCTRL,[],2,'omitnan')./sqrt(size(thisCTRL,2));

    % pull turn peak from pulse onset to window end
    s = s+1;
    if s>=floor(nSweep/2)
        [peakKIR(s,:),peakIdxKIR(s,:)] = max(thisKIR(idxSweepStart:end,:));
        [peakCTRL(s,:),peakIdxCTRL(s,:)] = max(thisCTRL(idxSweepStart:end,:));
    else
        [peakKIR(s,:),peakIdxKIR(s,:)] = min(thisKIR(idxSweepStart:end,:));
        [peakCTRL(s,:),peakIdxCTRL(s,:)] = min(thisCTRL(idxSweepStart:end,:));
    end

    % generate plot
    nexttile; hold on
    plot(t_mop,meanKIR,'-','Color',settings.geneColor{1},'Linewidth',1.5)
    sp1 = patch([t_mop'; flipud(t_mop')],[meanKIR-semKIR; flipud(meanKIR+semKIR)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{1};
    
    plot(t_mop,meanCTRL,'-','Color',[0.2 0.2 0.2],'Linewidth',1.5) % pooled control
    sp2 = patch([t_mop'; flipud(t_mop')],[meanCTRL-semCTRL; flipud(meanCTRL+semCTRL)], 'k', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp2.FaceColor = [0.2 0.2 0.2];
    
    xlim([0 tmax])
    ylim(ang_range)
    yline(0)
    if ns==firstSweep
        ylabel('Angular Velocity (deg/s)')
    end
end
sgtitle(thisTitle)
xlabel('Time (msec.)')

% save plot
cd(folder.summary)
plotname = append('mopulse_', thisSave,'_combinedsummary');
saveas(gcf,append(plotname,'.png'));
copyfile(append(plotname,'.png'), folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, append(plotname,'.svg'))
copyfile(append(plotname,'.svg'), folder.dropbox,'f');

%% motion pulse: plot directional velocity means together pool L+R with threshold
disp('Analyzing pursuit behavior during motion pulse experiment...')
% variables
ang_range = [-250 250]; %velocity range
tmax = max(t_mop); %time range
idxSweepStart = find(~isnan(pos_mopR(:,1)),1,'first');
firstSweep = 1; %first sweep to plot
sweepSelect = firstSweep:nSweep;

% initialize
figure; set(gcf,'Position',[100 100 1800 900])
tiledlayout(2,length(sweepSelect),'TileSpacing','compact')
peakKIR = []; peakWT = []; peakNA = [];
peakIdxKIR = []; peakIdxWT = []; peakIdxGFP = [];

% plot panel position
for ns = sweepSelect
    nexttile
    plot(t_mop,pos_mopR(:,ns),'-','Color','k','Linewidth',1.5)
    xlim([0 tmax])
    ylim([-160 160])
    yline(0)
    if ns==firstSweep
        ylabel('Object Position (deg)')
    end
end

% plot angular velocity
s = 0;
for ns = sweepSelect
    % fetch data for this sweep
    thisTitle = "Combined RL Sweeps";
    thisSave = "RL";
    thisKIR = kir_mopRLthresh(:,:,ns);
    thisWT = wt_mopRLthresh(:,:,ns);
    thisNA = na_mopRLthresh(:,:,ns);

    % calculate mean and SEM
    meanKIR = mean(thisKIR,2,'omitnan');
    meanWT = mean(thisWT,2,'omitnan');
    meanNA = mean(thisNA,2,'omitnan');
    semKIR = std(thisKIR,[],2,'omitnan')./sqrt(nKIR);
    semWT = std(thisWT,[],2,'omitnan')./sqrt(nWT);
    semNA = std(thisNA,[],2,'omitnan')./sqrt(nNA);

    % pull turn peak from pulse onset to window end
    s = s+1;
    if s>=floor(nSweep/2)
        [peakKIR(s,:),peakIdxKIR(s,:)] = max(thisKIR(idxSweepStart:end,:));
        [peakWT(s,:),peakIdxWT(s,:)] = max(thisWT(idxSweepStart:end,:));
        [peakNA(s,:),peakIdxGFP(s,:)] = max(thisNA(idxSweepStart:end,:));
    else
        [peakKIR(s,:),peakIdxKIR(s,:)] = min(thisKIR(idxSweepStart:end,:));
        [peakWT(s,:),peakIdxWT(s,:)] = min(thisWT(idxSweepStart:end,:));
        [peakNA(s,:),peakIdxGFP(s,:)] = min(thisNA(idxSweepStart:end,:));
    end

    % generate plot
    nexttile; hold on
    plot(t_mop,meanKIR,'-','Color',settings.geneColor{1},'Linewidth',1.5)
    sp1 = patch([t_mop'; flipud(t_mop')],[meanKIR-semKIR; flipud(meanKIR+semKIR)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{1};
    plot(t_mop,meanWT,'-','Color',settings.geneColor{2},'Linewidth',1.5)
    sp2 = patch([t_mop'; flipud(t_mop')],[meanWT-semWT; flipud(meanWT+semWT)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp2.FaceColor = settings.geneColor{2};
    plot(t_mop,meanNA,'-','Color',settings.geneColor{3},'Linewidth',1.5)
    sp3 = patch([t_mop'; flipud(t_mop')],[meanNA-semNA; flipud(meanNA+semNA)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp3.FaceColor = settings.geneColor{3};
    xlim([0 tmax])
    ylim(ang_range)
    yline(0)
    if ns==firstSweep
        ylabel('Angular Velocity (deg/s)')
    end
end
sgtitle(thisTitle)
xlabel('Time (msec.)')

% save plot
cd(folder.summary)
plotname = append('mopulse_', thisSave,'_threshsummary');
saveas(gcf,append(plotname,'.png'));
copyfile(append(plotname,'.png'), folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, append(plotname,'.svg'))
copyfile(append(plotname,'.svg'), folder.dropbox,'f');


% plot a quick summary
figure; set(gcf,'Position',[100 100 1200 900]);tiledlayout(1,2)
% plot peak angular velocity
nexttile; hold on
% plot trial data
plot(sweepSelect+offset(1),peakKIR,'.','Color',settings.trialColor)
plot(sweepSelect+offset(2),peakWT,'.','Color',settings.trialColor)
plot(sweepSelect+offset(3),peakNA,'.','Color',settings.trialColor)
% plot mean +/- SEM
errorbar(sweepSelect+offset(1),mean(peakKIR,2,'omitnan'),std(peakKIR,[],2,'omitnan')./sqrt(nKIR),'o','Color',settings.geneColor{1},'Linewidth',1)
errorbar(sweepSelect+offset(2),mean(peakWT,2,'omitnan'),std(peakWT,[],2,'omitnan')./sqrt(nWT),'o','Color',settings.geneColor{2},'Linewidth',1)
errorbar(sweepSelect+offset(3),mean(peakNA,2,'omitnan'),std(peakNA,[],2,'omitnan')./sqrt(nNA),'o','Color',settings.geneColor{3},'Linewidth',1)
axis padded
xticks(sweepSelect); xline(5.5); yline(0); ylim([-350 350])
xlabel('Sweep Order'); ylabel('Peak Angular Velocity (deg/s)'); sgtitle(thisTitle)

% plot time of peak angular velocity
% convert from idx to time
peaktKIR = t_mop(peakIdxKIR);
peaktWT = t_mop(peakIdxWT);
peaktGFP = t_mop(peakIdxGFP);
% remove empties
peaktKIR(peaktKIR==0) = nan;
peaktWT(peaktWT==0) = nan;
peaktGFP(peaktGFP==0) = nan;

nexttile; hold on
% plot trial data
plot(sweepSelect+offset(1),peaktKIR,'.','Color',settings.trialColor)
plot(sweepSelect+offset(2),peaktWT,'.','Color',settings.trialColor)
plot(sweepSelect+offset(3),peaktGFP,'.','Color',settings.trialColor)
% plot mean +/- SEM
errorbar(sweepSelect+offset(1),mean(peaktKIR,2,'omitnan'),std(peaktKIR,[],2,'omitnan')./sqrt(nKIR),'o','Color',settings.geneColor{1},'Linewidth',1)
errorbar(sweepSelect+offset(2),mean(peaktWT,2,'omitnan'),std(peaktWT,[],2,'omitnan')./sqrt(nWT),'o','Color',settings.geneColor{2},'Linewidth',1)
errorbar(sweepSelect+offset(3),mean(peaktGFP,2,'omitnan'),std(peaktGFP,[],2,'omitnan')./sqrt(nNA),'o','Color',settings.geneColor{3},'Linewidth',1)
axis padded
xticks(sweepSelect); xline(5.5); yline(0); %ylim([-120 120])
xlabel('Sweep Order'); ylabel('Peak Time (msec)'); sgtitle(thisTitle)

% save plot
cd(folder.summary)
plotname = append('mopulse_', thisSave,'_threshplot');
saveas(gcf,append(plotname,'.png'));
copyfile(append(plotname,'.png'), folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, append(plotname,'.svg'))
copyfile(append(plotname,'.svg'), folder.dropbox,'f');

% Save turn response for behavior
cd(folder.compare);
saveFilename = 'HeadClosedBehaviorThresh.mat';
save(saveFilename, 'peakKIR', 'peakWT', 'peakNA','wt_mopRLthresh','na_mopRLthresh');

%% motion pulse: histogram distribution
figure; set(gcf,'Position',[100 100 1800 800])
tiledlayout(3,4,'TileSpacing','compact')
fBin = fwdHist(:,1);
aBin = angHist(:,1);
sBin = sidHist(:,1);

for v = 1:3
    switch v
        case 1 %forward
            hBin = fBin;
            kirVel = kir_mopFwdHist;
            wtVel = wt_mopFwdHist;
            naVel = na_mopFwdHist;
            xvar = 'Binned Forward Velocity (mm/s)';
        case 2 %angular
            hBin = aBin;
            kirVel = kir_mopAngHist;
            wtVel = wt_mopAngHist;
            naVel = na_mopAngHist;
            xvar = 'Binned Angular Velocity (deg/s)';
        case 3 %sideways
            hBin = sBin;
            kirVel = kir_mopSidHist;
            wtVel = wt_mopSidHist;
            naVel = na_mopSidHist;
            xvar = 'Binned Sideways Velocity (mm/s)';
    end
    % calculate mean and sem
    kirVel_mean = mean(kirVel,2);
    wtVel_mean = mean(wtVel,2);
    naVel_mean = mean(naVel,2);
    kirVel_sem = std(kirVel,[],2)./sqrt(nKIR);
    wtVel_sem = std(wtVel,[],2)./sqrt(nKIR);
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
sgtitle('Motion Pulse - Velocity Distribution')
cd(folder.summary)
plotname = 'hist_moppursuit';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

%% optomotor: plut directional velocity trials+means separately
disp('Plotting optomotor control...')

ang_range = [-150 150];
for g = 1:3
    % initialize
    figure; set(gcf,'Position',[100 100 500 700])
    tiledlayout(1,2,'TileSpacing','compact')

    % fetch dataset
    switch g
        case 1
            thisG = settings.geneLabel{g};
            optoR = kir_optoR;
            optoL = kir_optoL;
        case 2
            thisG = settings.geneLabel{g};
            optoR = wt_optoR;
            optoL = wt_optoL;
        case 3
            thisG = settings.geneLabel{g};
            optoR = na_optoR;
            optoL = na_optoL;
    end
    % calculate mean
    mean_optoR = mean(optoR,2,'omitnan');
    mean_optoL = mean(optoL,2,'omitnan');

    % generate plot
    nexttile; hold on
    plot(t_opto,optoR,'-','Color',settings.trialColor)
    plot(t_opto,mean_optoR,'-','Color',settings.geneColor{g},'Linewidth',1.5)
    axis tight
    ylim(ang_range); yline(0); xline(sweepStartStop,':')
    ylabel('Angular Velocity (deg/s)')

    nexttile; hold on
    plot(t_opto,optoL,'-','Color',settings.trialColor)
    plot(t_opto,mean_optoL,'-','Color',settings.geneColor{g},'Linewidth',1.5)
    axis tight
    ylim(ang_range); yline(0); xline(sweepStartStop,':')
    xlabel('Time (msec.)')

    sgtitle(join([thisG 'optomotor']," "))
    % save plot
    cd(folder.summary)
    plotname = append(thisG,'_optomotor','_trials');
    saveas(gcf,append(plotname,'.png'));
    copyfile(append(plotname,'.png'), folder.dropbox,'f');
    % save vectorized plot
    cd(folder.vector)
    set(gcf,'renderer','Painters')
    saveas(gcf, append(plotname,'.svg'))
    copyfile(append(plotname,'.svg'), folder.dropbox,'f');
end

disp('Done.')

%% optomotor: plut directional velocity means together
disp('Plotting optomotor control...')

% initialize
figure; set(gcf,'Position',[100 100 750 700])
tiledlayout(1,3,'TileSpacing','compact')
ang_range = [-70 70];

for s = 1:3
    % fetch dataset
    switch s
        case 1 %right
            thisName = 'CW';
            optoKIR = kir_optoR;
            optoWT = wt_optoR;
            optoNA = na_optoR;
        case 2 %left
            thisName = 'CCW';
            optoKIR = kir_optoL;
            optoWT = wt_optoL;
            optoNA = na_optoL;
        case 3 %right+left
            thisName = 'CW+CCW';
            optoKIR = kir_optoRL;
            optoWT = wt_optoRL;
            optoNA = na_optoRL;
    end
    % calculate mean and SEM
    mean_optoKIR = mean(optoKIR,2,'omitnan');
    mean_optoWT = mean(optoWT,2,'omitnan');
    mean_optoNA = mean(optoNA,2,'omitnan');

    sem_optoKIR = std(optoKIR,[],2,'omitnan')./sqrt(nKIR);
    sem_optoWT = std(optoWT,[],2,'omitnan')./sqrt(nWT);
    sem_optoNA = std(optoNA,[],2,'omitnan')./sqrt(nNA);

    % generate plot
    nexttile; hold on
    plot(t_opto,mean_optoKIR,'-','Color',settings.geneColor{1},'Linewidth',1.5)
    sp1 = patch([t_opto'; flipud(t_opto')],[mean_optoKIR-sem_optoKIR; flipud(mean_optoKIR+sem_optoKIR)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{1};
    plot(t_opto,mean_optoWT,'-','Color',settings.geneColor{2},'Linewidth',1.5)
    sp1 = patch([t_opto'; flipud(t_opto')],[mean_optoWT-sem_optoWT; flipud(mean_optoWT+sem_optoWT)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{2};
    plot(t_opto,mean_optoNA,'-','Color',settings.geneColor{3},'Linewidth',1.5)
    sp1 = patch([t_opto'; flipud(t_opto')],[mean_optoNA-sem_optoNA; flipud(mean_optoNA+sem_optoNA)], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
    sp1.FaceColor = settings.geneColor{3};
    axis tight
    ylim(ang_range); yline(0); xline(sweepStartStop,':')
    ylabel('Angular Velocity (deg/s)'); title(thisName)
end
xlabel('Time (s)')
sgtitle("optomotor summary")

% save plot
cd(folder.summary)
plotname = append('optomotor','_summary');
saveas(gcf,append(plotname,'.png'));
copyfile(append(plotname,'.png'), folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, append(plotname,'.svg'))
copyfile(append(plotname,'.svg'), folder.dropbox,'f');

plotOptoMaxes(kir_optoRL, wt_optoRL, na_optoRL, folder)
% save plot
cd(folder.summary)
plotname = append('optomotor','_maxsummary');
saveas(gcf,append(plotname,'.png'));
copyfile(append(plotname,'.png'), folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, append(plotname,'.svg'))
copyfile(append(plotname,'.svg'), folder.dropbox,'f');

disp('Done.')

%% optomotor: histogram distribution
figure; set(gcf,'Position',[100 100 1800 800])
tiledlayout(3,4,'TileSpacing','compact')
fBin = fwdHist(:,1);
aBin = angHist(:,1);
sBin = sidHist(:,1);

for v = 1:3
    switch v
        case 1 %forward
            hBin = fBin;
            kirVel = kir_optoFwdHist;
            wtVel = wt_optoFwdHist;
            naVel = na_mopFwdHist;
            xvar = 'Binned Forward Velocity (mm/s)';
        case 2 %angular
            hBin = aBin;
            kirVel = kir_optoAngHist;
            wtVel = wt_optoAngHist;
            naVel = na_optoAngHist;
            xvar = 'Binned Angular Velocity (deg/s)';
        case 3 %sideways
            hBin = sBin;
            kirVel = kir_optoSidHist;
            wtVel = wt_optoSidHist;
            naVel = na_optoSidHist;
            xvar = 'Binned Sideways Velocity (mm/s)';
    end
    % calculate mean and sem
    kirVel_mean = mean(kirVel,2);
    wtVel_mean = mean(wtVel,2);
    naVel_mean = mean(naVel,2);
    kirVel_sem = std(kirVel,[],2)./sqrt(nKIR);
    wtVel_sem = std(wtVel,[],2)./sqrt(nKIR);
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
sgtitle('Optomotor - Velocity Distribution')
cd(folder.summary)
plotname = 'hist_opto';
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

