% pipeline_kir_menotaxis
%
% Pipeline Function
% Pulls all processed files from ALL flies in a given experiment, performs
% necessary analyses and plots accordingly. The experiment involves
% menotaxis behavior with jumps, comparing flies with and without KIR perturbation.
%
% INPUTS
% exptFolder - overarching experiment folder
%
% The function pools data across all flies, analyzes heading direction,
% jump responses, and basic behavior parameters (running time, speed, etc.).
% Results are saved as plots, including histograms, comparison plots, 
% and summary figures.
%
% 02/27/2024 - MC created
% 07/10/2024 - MC cleaned up
%
function pipeline_kir_menotaxis(exptFolder)
%% initialize 
disp('STARTING ANALYSES FOR POOLED KIR MENOTAXIS...')
close all

% set filename info and create necessary directories
filebase = strrep(exptFolder,' ','_');
folder = generateFolders(exptFolder);

% load processing settings
settings = processSettings();

%% set key variables

jumps_simp = [30, 60, 180]; %deg
nCon = length(settings.geneLabel);
minRunTime = 6*60; %s, min time spend running to be considered "running well"

% set offset for overlapping data
o = 0.1;
offset = [-o; 0; o]; %used to offset overlapping plots

%% load in and pool all trials from each experiment folder
disp('Loading in datasets...')

% find all files in this directory
% pull all file info
cd(folder.int)
allFiles = dir('*int.mat');
nFlies = length(allFiles);

% initialize data storage arrays
nKIR = 0;
nWT = 0;
nNA = 0;

for e = 1:nFlies
    disp(['Processing fly ' num2str(e) '/' num2str(nFlies) '...'])
    % load this trial
    thisTrial = allFiles(e).name;
    load(thisTrial)

    % find timepoints where the fly was running and for how long in total
    [run_forward,run_angular,~,run_panelps] = pursuitFinder(int_forward,int_angular,int_sideway,int_panelps,int_time,settings.runThreshB);
    thisN = size(int_forward,2);
    thisRunTime(e) = sum(int_time(sum(~isnan(run_forward))+1)); %sec
    thisRunSpeed = mean(run_forward,'all','omitnan'); %avg run speed
    thisTurnSpeed = mean(abs(run_angular),'all','omitnan'); %avg turn speed

    % ANALYSIS #1: if fly ran well, infer menotaxis goal HD direction
    if thisRunTime(e)>minRunTime
        ph = polarhistogram(deg2rad(run_panelps),settings.HDBins,'Normalization','probability');
        [thisGoal_val,hd_idx] = max(ph.Values); %find most prominent heading
        thisGoal_hd = (ph.BinEdges(hd_idx)+ph.BinEdges(hd_idx+1))/2; %assign heading bin
        
        % analyze turning v object position relative to goal
        [posvang, posvangRL, posBins] = setpoint_errorvturn(run_panelps, run_angular, int_time, settings, 1, 0);

        close all
    else
        thisGoal_val = nan;
        thisGoal_hd = nan;
        posvang=nan;
        povangRL=nan;
    end

    % ANALYSIS #2: perform jump analysis
    [jumpData,restoreFT] = barjump_analysis(int_panelps,int_jumptrg,int_forward,int_angular,0,int_time,0);

    % pool this trial data according to the assigned trial condition
    % but only if the fly made a sufficient number of corrections
    if ~sum((restoreFT.trialsWithCorrect)<2)
        if contains(thisTrial,'KIR') %kir perturbation flies
            nKIR = nKIR+1;
            kirTrialN(nKIR) = thisN;
            kirRunTime(nKIR) = thisRunTime(e);
            kirRunSpeed(nKIR) = thisRunSpeed;
            kirTurnSpeed(nKIR) = thisTurnSpeed;
            kirGoalHD(nKIR,:) = [thisGoal_hd,thisGoal_val];

            kirJumpHD(:,:,nKIR) = jumpData.panelps;
            kirJumpHD_STD(:,:,nKIR) = jumpData.panelpsstd;
            kirJumpHDpool(:,:,nKIR) = jumpData.panelpsPool;
            kirJumpAng(:,:,nKIR) = jumpData.angular;
            kirJumpAng_STD(:,:,nKIR) = jumpData.angularstd;
            kirJumpAVG(nKIR,:) = restoreFT.meanPool;
            kirJumpSTD(nKIR,:) = restoreFT.stdPool;
            kirJumpVARpre(nKIR,:) = restoreFT.prevarPool;
            kirJumpVAR(nKIR,:) = restoreFT.varPool;
            kirJumpNCor(nKIR,:) = restoreFT.trialsWithCorrect;
            kirJumpNoCor(nKIR,:) = restoreFT.trialsWithoutCorrect;
            kirJumpPER(nKIR,:) = restoreFT.percentPool;
            kirJumpFwdMax(nKIR,:) = jumpData.fwdMax;
            kirJumpAngMax(nKIR,:) = jumpData.angMax;

            kirEVT(:,nKIR) = posvangRL;

        elseif contains(thisTrial,'WT') %wildtype control flies
            nWT = nWT+1;
            wtTrialN(nWT) = thisN;
            wtRunTime(nWT) = thisRunTime(e);
            wtRunSpeed(nWT) = thisRunSpeed;
            wtTurnSpeed(nWT) = thisTurnSpeed;
            wtGoalHD(nWT,:) = [thisGoal_hd,thisGoal_val];

            wtJumpHD(:,:,nWT) = jumpData.panelps;
            wtJumpHD_STD(:,:,nWT) = jumpData.panelpsstd;
            wtJumpHDpool(:,:,nWT) = jumpData.panelpsPool;
            wtJumpAng(:,:,nWT) = jumpData.angular;
            wtJumpAng_STD(:,:,nWT) = jumpData.angularstd;
            wtJumpAVG(nWT,:) = restoreFT.meanPool;
            wtJumpSTD(nWT,:) = restoreFT.stdPool;
            wtJumpVARpre(nWT,:) = restoreFT.prevarPool;
            wtJumpVAR(nWT,:) = restoreFT.varPool;
            wtJumpNCor(nWT,:) = restoreFT.trialsWithCorrect;
            wtJumpNoCor(nWT,:) = restoreFT.trialsWithoutCorrect;
            wtJumpPER(nWT,:) = restoreFT.percentPool;
            wtJumpFwdMax(nWT,:) = jumpData.fwdMax;
            wtJumpAngMax(nWT,:) = jumpData.angMax;

            wtEVT(:,nWT) = posvangRL;

        elseif contains(thisTrial,'GFP') %na control flies
            nNA=nNA+1;
            naTrialN(nNA) = thisN;
            naRunTime(nNA) = thisRunTime(e);
            naRunSpeed(nNA) = thisRunSpeed;
            naTurnSpeed(nNA) = thisTurnSpeed;
            naGoalHD(nNA,:) = [thisGoal_hd,thisGoal_val];

            naJumpHD(:,:,nNA) = jumpData.panelps;
            naJumpHD_STD(:,:,nNA) = jumpData.panelpsstd;
            naJumpHDpool(:,:,nNA) = jumpData.panelpsPool;
            naJumpAng(:,:,nNA) = jumpData.angular;
            naJumpAng_STD(:,:,nNA) = jumpData.angularstd;
            naJumpAVG(nNA,:) = restoreFT.meanPool;
            naJumpSTD(nNA,:) = restoreFT.stdPool;
            naJumpVARpre(nNA,:) = restoreFT.prevarPool;
            naJumpVAR(nNA,:) = restoreFT.varPool;
            naJumpNCor(nNA,:) = restoreFT.trialsWithCorrect;
            naJumpNoCor(nNA,:) = restoreFT.trialsWithoutCorrect;
            naJumpPER(nNA,:) = restoreFT.percentPool;
            naJumpFwdMax(nNA,:) = jumpData.fwdMax;
            naJumpAngMax(nNA,:) = jumpData.angMax;

            naEVT(:,nNA) = posvangRL;

        end
    else
        disp([thisTrial ' omitted due to insufficient data'])
    end
end
% pull jump window time
timeJump = jumpData.time;

disp(['All datasets loaded: ' 'KIR = ' num2str(nKIR) ', WT = ' num2str(nWT) ', NA = ' num2str(nNA)])


%% plot goal HD across settings.geneLabel
disp('Comparing goal HD across settings.geneLabel...')

% initialize
figure; set(gcf,'Position',[100 100 1400 400])
tiledlayout(1,3,'TileSpacing','compact')

% plot KIR HD polar
nexttile
for n = 1:nKIR
    polarplot([kirGoalHD(n,1);kirGoalHD(n,1)],[0;kirGoalHD(n,2)])
    hold on
end
title(settings.geneLabel{1})

% plot WT HD polar
nexttile
for n = 1:nWT
    polarplot([wtGoalHD(n,1);wtGoalHD(n,1)],[0;wtGoalHD(n,2)])
    hold on
end
title(settings.geneLabel{2})

% plot NA HD polar
nexttile
for n = 1:nNA
    polarplot([naGoalHD(n,1);naGoalHD(n,1)],[0;naGoalHD(n,2)])
    hold on
end
title(settings.geneLabel{3})

% save plot
cd(folder.summary)
plotname = 'summary_goalHD';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

disp('Complete.')

% plot heading summary
% initialize
figure; set(gcf,'Position',[100 100 600 400])
tiledlayout(1,2,'TileSpacing','compact')

% vector direction
nexttile
for c = 1:3
    switch c
        case 1 %kir
            allStrengths = abs(rad2deg(kirGoalHD(:,1)));
            avgStrengths = mean(allStrengths);
            semStrengths = std(allStrengths)./sqrt(nKIR);
        case 2 %wt
            allStrengths = abs(rad2deg(wtGoalHD(:,1)));
            avgStrengths = mean(allStrengths);
            semStrengths = std(allStrengths)./sqrt(nWT);
        case 3 %na
            allStrengths = abs(rad2deg(naGoalHD(:,1)));
            avgStrengths = mean(allStrengths);
            semStrengths = std(allStrengths)./sqrt(nNA);
    end
    plot(c,allStrengths,'.','Color', settings.trialColor)
    hold on
    errorbar(c,avgStrengths,semStrengths,'o','Color', settings.geneColor{c},'LineWidth',1.5)
end
axis padded
xticks([1 2 3])
xticklabels(settings.geneLabel)
%ylim([0 .15])
ylabel('goal direction (deg)')
xlabel('genotype')

% vector strength
nexttile
for c = 1:3
    switch c
        case 1 %kir
            allStrengths = kirGoalHD(:,2);
            avgStrengths = mean(allStrengths);
            semStrengths = std(allStrengths)./sqrt(nKIR);
        case 2 %wt
            allStrengths = wtGoalHD(:,2);
            avgStrengths = mean(allStrengths);
            semStrengths = std(allStrengths)./sqrt(nWT);
        case 3 %na
            allStrengths = naGoalHD(:,2);
            avgStrengths = mean(allStrengths);
            semStrengths = std(allStrengths)./sqrt(nNA);
    end
    plot(c,allStrengths,'.','Color', settings.trialColor)
    hold on
    errorbar(c,avgStrengths,semStrengths,'o','Color', settings.geneColor{c},'LineWidth',1.5)
end
axis padded
xticks([1 2 3])
xticklabels(settings.geneLabel)
ylim([0 .15])
ylabel('goal probability / goal vector strength')
xlabel('genotype')

% save plot
cd(folder.summary)
plotname = 'summary2_goalHD';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

disp('Complete.')

%% plot jump response across settings.geneLabel
disp('Comparing jump response across settings.geneLabel...')

% initialize
figure; set(gcf,'Position',[100 100 1800 850])
tiledlayout(3,6,'TileSpacing','compact')
t_range = [-2 20]; %time window
p_range = [-180 180]; %panels range

% for each fly genotype condition
for c = 1:3
    % for each jump options
    for j = 1:6
        % select the dataset and calculate corresponding mean
        switch c
            case 1
                thisJumpData = reshape(kirJumpHD(:,j,:),[],nKIR);
                thisJumpMean = mean(thisJumpData,2,'omitnan');
            case 2
                thisJumpData = reshape(wtJumpHD(:,j,:),[],nWT);
                thisJumpMean = mean(thisJumpData,2,'omitnan');
            case 3
                thisJumpData = reshape(naJumpHD(:,j,:),[],nNA);
                thisJumpMean = mean(thisJumpData,2,'omitnan');
        end
        
        % plot
        nexttile
        plot(timeJump,thisJumpData,'Color', [0.4 0.4 0.4])
        hold on
        plot(timeJump,thisJumpMean,'Color',settings.geneColor{c},'LineWidth',2)
        xline(0);yline(0)
        xlim(t_range)
        ylim(p_range)

        if j==1
            ylabel([settings.geneLabel{c} ' obj pos (deg)'])
        end
        if c==1
            title([num2str(settings.jumps_inuse(j)) ' jump'])
        elseif c==3
            xlabel('time (s)')
        end
    end
end

% save plot
cd(folder.summary)
plotname = 'summary_jumpHD';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');


% initialize
figure; set(gcf,'Position',[100 100 1800 600])
tiledlayout(1,6,'TileSpacing','compact')
p_range = [-180 180]; %panels range

% offset correction
[jumpIdx] = fetchTimeIdx(timeJump,-0.01);

% for each jump options
for j = 1:6
    % select the dataset and calculate corresponding mean
    nexttile
    % plot kir
    thisJumpData = reshape(kirJumpHD(:,j,:),[],nKIR);
    thisJumpMean = mean(thisJumpData,2,'omitnan');
    thisJumpMean = thisJumpMean - thisJumpMean(jumpIdx);
    plot(timeJump,thisJumpMean,'Color',settings.geneColor{1},'LineWidth',1.5)
    hold on
    % plot wt control
    thisJumpData = reshape(wtJumpHD(:,j,:),[],nWT);
    thisJumpMean = mean(thisJumpData,2,'omitnan');
    thisJumpMean = thisJumpMean - thisJumpMean(jumpIdx);
    plot(timeJump,thisJumpMean,'Color',settings.geneColor{2},'LineWidth',1.5)
    hold on
    % plot na control
    thisJumpData = reshape(naJumpHD(:,j,:),[],nNA);
    thisJumpMean = mean(thisJumpData,2,'omitnan');
    thisJumpMean = thisJumpMean - thisJumpMean(jumpIdx);
    plot(timeJump,thisJumpMean,'Color',settings.geneColor{3},'LineWidth',1.5)

    xline(0);yline(0)
    xlim(t_range)
    ylim(p_range)

    title([num2str(settings.jumps_inuse(j)) ' jump'])
    xlabel('time (s)')
end
legend(settings.geneLabel)


% save plot
cd(folder.summary)
plotname = 'summary_jumpHD_mean';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');


%% plot jump response across settings.geneLabel pooled
disp('Comparing jump response across settings.geneLabel pooled...')

% initialize
figure; set(gcf,'Position',[100 100 1000 850])
tiledlayout(3,3,'TileSpacing','compact')
t_range = [-2 20]; %time window
p_range = [-180 180]; %panels range

% for each fly genotype condition
for c = 1:3
    % for each jump options
    for j = 1:3
        % select the dataset and calculate corresponding mean
        switch c
            case 1
                thisJumpData = reshape(kirJumpHDpool(:,j,:),[],nKIR);
                thisJumpMean = mean(thisJumpData,2,'omitnan');
            case 2
                thisJumpData = reshape(wtJumpHDpool(:,j,:),[],nWT);
                thisJumpMean = mean(thisJumpData,2,'omitnan');
            case 3
                thisJumpData = reshape(naJumpHDpool(:,j,:),[],nNA);
                thisJumpMean = mean(thisJumpData,2,'omitnan');
        end
        
        % plot
        nexttile
        plot(timeJump,thisJumpData,'Color', [0.4 0.4 0.4])
        hold on
        plot(timeJump,thisJumpMean,'Color',settings.geneColor{c},'LineWidth',2)
        xline(0);yline(0)
        xlim(t_range)
        ylim(p_range)

        if j==1
            ylabel([settings.geneLabel{c} ' obj pos (deg)'])
        end
        if c==1
            title([num2str(jumps_simp(j)) ' jump'])
        elseif c==3
            xlabel('time (s)')
        end
    end
end

% save plot
cd(folder.summary)
plotname = 'summary_jumpHDpool';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');


% initialize
figure; set(gcf,'Position',[100 100 1000 600])
tiledlayout(1,3,'TileSpacing','compact')
p_range = [-180 180]; %panels range

% for each jump options
for j = 1:3
    % select the dataset and calculate corresponding mean
    nexttile
    % plot kir
    thisJumpData = reshape(kirJumpHDpool(:,j,:),[],nKIR);
    thisJumpMean = mean(thisJumpData,2,'omitnan');
    thisJumpMean = thisJumpMean - thisJumpMean(jumpIdx);
    plot(timeJump,thisJumpMean,'Color',settings.geneColor{1},'LineWidth',1.5)
    hold on
    % plot wt control
    thisJumpData = reshape(wtJumpHDpool(:,j,:),[],nWT);
    thisJumpMean = mean(thisJumpData,2,'omitnan');
    thisJumpMean = thisJumpMean - thisJumpMean(jumpIdx);
    plot(timeJump,thisJumpMean,'Color',settings.geneColor{2},'LineWidth',1.5)
    hold on
    % plot na control
    thisJumpData = reshape(naJumpHDpool(:,j,:),[],nNA);
    thisJumpMean = mean(thisJumpData,2,'omitnan');
    thisJumpMean = thisJumpMean - thisJumpMean(jumpIdx);
    plot(timeJump,thisJumpMean,'Color',settings.geneColor{3},'LineWidth',1.5)

    xline(0);yline(0)
    xlim(t_range)
    ylim(p_range)

    title([num2str(jumps_simp(j)) ' jump'])
    xlabel('time (s)')
end
legend(settings.geneLabel)


% save plot
cd(folder.summary)
plotname = 'summary_jumpHDpool_mean';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

%% compare basic behavior parameters

% initialize
figure; set(gcf,'Position',[100 100 400 900])
tiledlayout(3,2,'TileSpacing','compact')

% compare number of trials collected
nexttile
for c = 1:3
    % for each condition
    switch c
        case 1 %kir
            thisN = kirTrialN;
            thisMean = mean(thisN);
            thisSEM = std(thisN)./sqrt(nKIR);
        case 2 %wt
            thisN = wtTrialN;
            thisMean = mean(thisN);
            thisSEM = std(thisN)./sqrt(nWT);
        case 3 %na
            thisN = wtTrialN;
            thisMean = mean(thisN);
            thisSEM = std(thisN)./sqrt(nNA);
    end

    % plot
    hold on
    plot(c, thisN,'.','Color',settings.trialColor)
    errorbar(c, thisMean,thisSEM,'o','Color', settings.geneColor{c},'LineWidth', 1.5)
end
axis padded
xticks([1 2 3])
xticklabels(settings.geneLabel)
ylabel('total trials run')
xlabel('genotype')

% compare run times
nexttile
for c = 1:3
    % for each condition
    switch c
        case 1 %kir
            thisRun = kirRunTime./60;
            thisMean = mean(thisRun);
            thisSEM = std(thisRun)./sqrt(nKIR);
        case 2 %wt
            thisRun = wtRunTime./60;
            thisMean = mean(thisRun);
            thisSEM = std(thisRun)./sqrt(nWT);
        case 3 %na
            thisRun = naRunTime./60;
            thisMean = mean(thisRun);
            thisSEM = std(thisRun)./sqrt(nNA);
    end

    % plot
    hold on
    plot(c, thisRun,'.','Color',settings.trialColor)
    errorbar(c, thisMean,thisSEM,'o','Color', settings.geneColor{c},'LineWidth', 1.5)
end
axis padded
xticks([1 2 3])
xticklabels(settings.geneLabel)
ylabel('total run time ALL (min)')
xlabel('genotype')

% compare run speeds
nexttile
for c = 1:3
    % for each condition
    switch c
        case 1 %kir
            thisSpeed = kirRunSpeed;
            thisMean = mean(thisSpeed);
            thisSEM = std(thisSpeed)./sqrt(nKIR);
        case 2 %wt
            thisSpeed = wtRunSpeed;
            thisMean = mean(thisSpeed);
            thisSEM = std(thisSpeed)./sqrt(nWT);
        case 3 %na
            thisSpeed = naRunSpeed;
            thisMean = mean(thisSpeed);
            thisSEM = std(thisSpeed)./sqrt(nNA);
    end

    % plot
    hold on
    plot(c, thisSpeed,'.','Color',settings.trialColor)
    errorbar(c, thisMean,thisSEM,'o','Color', settings.geneColor{c},'LineWidth', 1.5)
end
axis padded
xticks([1 2 3])
xticklabels(settings.geneLabel)
ylim([0 20])
ylabel('mean forward velocity ALL (mm/s)')
xlabel('genotype')

% compare turn speeds
nexttile
for c = 1:3
    % for each condition
    switch c
        case 1 %kir
            thisSpeed = kirTurnSpeed;
            thisMean = mean(thisSpeed);
            thisSEM = std(thisSpeed)./sqrt(nKIR);
        case 2 %wt
            thisSpeed = wtTurnSpeed;
            thisMean = mean(thisSpeed);
            thisSEM = std(thisSpeed)./sqrt(nWT);
        case 3 %na
            thisSpeed = naTurnSpeed;
            thisMean = mean(thisSpeed);
            thisSEM = std(thisSpeed)./sqrt(nNA);
    end

    % plot
    hold on
    plot(c, thisSpeed,'.','Color',settings.trialColor)
    errorbar(c, thisMean,thisSEM,'o','Color', settings.geneColor{c},'LineWidth', 1.5)
end
axis padded
xticks([1 2 3])
xticklabels(settings.geneLabel)
ylim([0 100])
ylabel('mean turn speed ALL (deg/s)')
xlabel('genotype')

% compare max forward velocities
nexttile
for c = 1:3
    % for each condition
    switch c
        case 1 %kir
            thisMaxGroup = median(kirJumpFwdMax,2);
            thisMean = mean(thisMaxGroup);
            thisSEM = std(thisMaxGroup)./sqrt(nKIR);
        case 2 %wt
            thisMaxGroup = median(wtJumpFwdMax,2);
            thisMean = mean(thisMaxGroup);
            thisSEM = std(thisMaxGroup)./sqrt(nWT);
        case 3 %na
            thisMaxGroup = median(naJumpFwdMax,2);
            thisMean = mean(thisMaxGroup);
            thisSEM = std(thisMaxGroup)./sqrt(nNA);
    end

    % plot
    hold on
    plot(c, thisMaxGroup,'.','Color',settings.trialColor)
    errorbar(c, thisMean,thisSEM,'o','Color', settings.geneColor{c},'LineWidth', 1.5)
end
axis padded
xticks([1 2 3])
xticklabels(settings.geneLabel)
ylim([0 40])
ylabel('max forward velocity during jump (mm/s)')
xlabel('genotype')

% compare max angular velocities
nexttile
for c = 1:3
    % for each condition
    switch c
        case 1 %kir
            thisMaxGroup = median(kirJumpAngMax,2);
            thisMean = mean(thisMaxGroup);
            thisSEM = std(thisMaxGroup)./sqrt(nKIR);
        case 2 %wt
            thisMaxGroup = median(wtJumpAngMax,2);
            thisMean = mean(thisMaxGroup);
            thisSEM = std(thisMaxGroup)./sqrt(nWT);
        case 3 %na
            thisMaxGroup = median(naJumpAngMax,2);
            thisMean = mean(thisMaxGroup);
            thisSEM = std(thisMaxGroup)./sqrt(nNA);
    end

    % plot
    hold on
    plot(c, thisMaxGroup,'.','Color',settings.trialColor)
    errorbar(c, thisMean,thisSEM,'o','Color', settings.geneColor{c},'LineWidth', 1.5)
end
axis padded
xticks([1 2 3])
xticklabels(settings.geneLabel)
ylim([100 500])
ylabel('max angular speed during jump (deg/s)')
xlabel('genotype')

% save plot
cd(folder.summary)
plotname = 'summary_basiccompare';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

%% compare jump performance

% initialize
nV = 3;
figure; set(gcf,'Position',[100 100 1000 600])
tiledlayout(1,3,'TileSpacing','compact')
% set offset for overlapping data

% compare N jumps NOT corrected
nexttile
for c = 1:3
    % for each condition
    switch c
        case 1 %kir
            thisJUMP = reshape(kirJumpNoCor,[],3);
            thisMean = mean(thisJUMP,1,'omitnan');
            thisSEM = std(thisJUMP,0,1,'omitnan')./sqrt(nKIR);
        case 2 %wt
            thisJUMP = reshape(wtJumpNoCor,[],3);
            thisMean = mean(thisJUMP,1,'omitnan');
            thisSEM = std(thisJUMP,0,1,'omitnan')./sqrt(nWT);
        case 3 %na
            thisJUMP = reshape(naJumpNoCor,[],3);
            thisMean = mean(thisJUMP,1,'omitnan');
            thisSEM = std(thisJUMP,0,1,'omitnan')./sqrt(nNA);
    end

    % plot
    hold on
    x = (1:nV)+offset(c);
    plot(x, thisJUMP',':.','Color',settings.trialColor)
    errorbar(x, thisMean,thisSEM,'o','Color', settings.geneColor{c},'LineWidth', 1.5)
end
axis padded
ylim([0 25])
xticks(1:6)
xticklabels(num2str(jumps_simp'))
ylabel('N jumps not corrected for')
xlabel('jump size (deg)')

% compare N jumps corrected
nexttile
for c = 1:3
    % for each condition
    switch c
        case 1 %kir
            thisJUMP = reshape(kirJumpNCor,[],3);
            thisMean = mean(thisJUMP,1,'omitnan');
            thisSEM = std(thisJUMP,0,1,'omitnan')./sqrt(nKIR);
        case 2 %wt
            thisJUMP = reshape(wtJumpNCor,[],3);
            thisMean = mean(thisJUMP,1,'omitnan');
            thisSEM = std(thisJUMP,0,1,'omitnan')./sqrt(nWT);
        case 3 %na
            thisJUMP = reshape(naJumpNCor,[],3);
            thisMean = mean(thisJUMP,1,'omitnan');
            thisSEM = std(thisJUMP,0,1,'omitnan')./sqrt(nNA);
    end

    % plot
    hold on
    x = (1:nV)+offset(c);
    plot(x, thisJUMP',':.','Color',settings.trialColor)
    errorbar(x, thisMean,thisSEM,'o','Color', settings.geneColor{c},'LineWidth', 1.5)
end
axis padded
ylim([0 25])
xticks(1:6)
xticklabels(num2str(jumps_simp'))
ylabel('N jumps not corrected for')
xlabel('jump size (deg)')

% compare correction percentage
nexttile
for c = 1:3
    % for each condition
    switch c
        case 1 %kir
            thisPER = kirJumpPER*100;
            thisMean = mean(thisPER,1,'omitnan');
            thisSEM = std(thisPER,0,1,'omitnan')./sqrt(nKIR);
        case 2 %wt
            thisPER = wtJumpPER*100;
            thisMean = mean(thisPER,1,'omitnan');
            thisSEM = std(thisPER,0,1,'omitnan')./sqrt(nWT);
        case 3 %na
            thisPER = wtJumpPER*100;
            thisMean = mean(thisPER,1,'omitnan');
            thisSEM = std(thisPER,0,1,'omitnan')./sqrt(nNA);
    end

    % plot
    hold on
    x = (1:nV)+offset(c);
    plot(x, thisPER',':.','Color',settings.trialColor)
    errorbar(x, thisMean,thisSEM,'o','Color', settings.geneColor{c},'LineWidth', 1.5)
end
axis padded
ylim([0 100])
xticks(1:6)
xticklabels(num2str(jumps_simp'))
ylabel('median jump correction rate (%)')
xlabel('jump size (deg)')

% save plot
cd(folder.summary)
plotname = 'summary_jumpcompare1';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

% initialize
nV = 3;
figure; set(gcf,'Position',[100 100 1000 600])
tiledlayout(1,3,'TileSpacing','compact')

% compare correction time
nexttile
for c = 1:3
    % for each condition
    switch c
        case 1 %kir
            thisAVG = kirJumpAVG;
            thisMean = mean(thisAVG,1,'omitnan');
            thisSEM = std(thisAVG,0,1,'omitnan')./sqrt(nKIR);
        case 2 %wt
            thisAVG = wtJumpAVG;
            thisMean = mean(thisAVG,1,'omitnan');
            thisSEM = std(thisAVG,0,1,'omitnan')./sqrt(nWT);
        case 3 %na
            thisAVG = naJumpAVG;
            thisMean = mean(thisAVG,1,'omitnan');
            thisSEM = std(thisAVG,0,1,'omitnan')./sqrt(nNA);
    end

    % plot
    hold on
    x = (1:nV)+offset(c);
    plot(x, thisAVG',':.','Color',settings.trialColor)
    errorbar(x, thisMean,thisSEM,'o','Color', settings.geneColor{c},'LineWidth', 1.5)
end
axis padded
xticks(1:6)
xticklabels(num2str(jumps_simp'))
ylabel('median jump correction time (s)')
xlabel('jump size (deg)')

% compare correction variance
nexttile
for c = 1:3
    % for each condition
    switch c
        case 1 %kir
            thisVAR = kirJumpVARpre;
            thisMean = mean(thisVAR,1,'omitnan');
            thisSEM = std(thisVAR,0,1,'omitnan')./sqrt(nKIR);
        case 2 %wt
            thisVAR = wtJumpVARpre;
            thisMean = mean(thisVAR,1,'omitnan');
            thisSEM = std(thisVAR,0,1,'omitnan')./sqrt(nWT);
        case 3 %na
            thisVAR = naJumpVARpre;
            thisMean = mean(thisVAR,1,'omitnan');
            thisSEM = std(thisVAR,0,1,'omitnan')./sqrt(nNA);
    end

    % plot
    x = (1:nV)+offset(c);
    plot(x, thisVAR',':.','Color',settings.trialColor)
    hold on
    errorbar(x, thisMean,thisSEM,'o','Color', settings.geneColor{c},'LineWidth', 1.5)
end
axis padded
xticks(1:6)
xticklabels(num2str(jumps_simp'))
ylim([0 3500])
ylabel('median pre-jump variance')
xlabel('jump size (deg)')

nexttile
for c = 1:3
    % for each condition
    switch c
        case 1 %kir
            thisVAR = kirJumpVAR;
            thisMean = mean(thisVAR,1,'omitnan');
            thisSEM = std(thisVAR,0,1,'omitnan')./sqrt(nKIR);
        case 2 %wt
            thisVAR = wtJumpVAR;
            thisMean = mean(thisVAR,1,'omitnan');
            thisSEM = std(thisVAR,0,1,'omitnan')./sqrt(nWT);
        case 3 %na
            thisVAR = naJumpVAR;
            thisMean = mean(thisVAR,1,'omitnan');
            thisSEM = std(thisVAR,0,1,'omitnan')./sqrt(nNA);
    end

    % plot
    hold on
    x = (1:nV)+offset(c);
    plot(x, thisVAR,':.','Color',settings.trialColor)
    errorbar(x, thisMean,thisSEM,'o','Color', settings.geneColor{c},'LineWidth', 1.5)
end
axis padded
xticks(1:6)
xticklabels(num2str(jumps_simp'))
ylim([0 3500])
ylabel('median post-jump variance')
xlabel('jump size (deg)')

% save plot
cd(folder.summary)
plotname = 'summary_jumpcompare2';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
% save vectorized plot
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');

disp('Complete.')

%% plot error v turning
% CREATED: 06/17/2024 - MC

% Set up tiled layout
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
titles = {'WT', 'kirEV', 'naEVT'};
colors = {[0 0 0], [0 0.3 0.7], [0.9 0.2 0]}; % Colors for genotype means

% Plot WT
nexttile;
hold on;
plot(posBins, wtEVT, 'Color', [0.8 0.8 0.8]); % Individual flies in grey
plot(posBins, mean(wtEVT, 2, 'omitnan'), 'Color', colors{1}, 'LineWidth', 1.5); % Mean
title(titles{1});
xlabel('Error from Goal HD (deg)');
ylabel('Rotational Velocity (deg/s)');
ylim([-20 20])
hold off;

% Plot kirEV
nexttile;
hold on;
plot(posBins, kirEVT, 'Color', [0.8 0.8 0.8]); % Individual flies in grey
plot(posBins, mean(kirEVT, 2, 'omitnan'), 'Color', colors{2}, 'LineWidth', 1.5); % Mean
title(titles{2});
xlabel('Error from Goal HD (deg)');
ylim([-20 20])
hold off;

% Plot naEVT
nexttile;
hold on;
plot(posBins, naEVT, 'Color', [0.8 0.8 0.8]); % Individual flies in grey
plot(posBins, mean(naEVT, 2, 'omitnan'), 'Color', colors{3}, 'LineWidth', 1.5); % Mean
title(titles{3});
xlabel('Error from Goal HD (deg)');
ylim([-20 20])
hold off;

% Final formatting
title(t, 'Rotational Velocity vs Object Position');

%% end
disp('ALL ANALYSES COMPLETE.')
end

