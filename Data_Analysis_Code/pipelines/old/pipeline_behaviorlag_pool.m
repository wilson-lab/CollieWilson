% pipeline_behaviorlag_pool
%
% Pipeline Function
% Pulls all processed files from ALL flies in a given experiment, each fly
% is then analyzed to determine the best fit lag time between spike rate
% and each directional velocity. Flies are then pooled at the end to
% estimate the best lag time across flies.
%
% INPUTS
% exptName - name of experiment selection
% exptFilePath - directory for all flies to be included in analysis
%
% 12/08/2022 - MC original
%

function pipeline_behaviorlag_pool(exptName,exptFilePath)
%% initialize
disp('STARTING LAG ANALYSES FOR SELECTED EXPERIMENT...')
close all

% set folder names
mainFolder = ['D:\' exptName]; %master folder
cd(mainFolder)
saveFolder = [mainFolder '\laganalysis']; %folder for saving individual data
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder)
end
summaryFolder = [mainFolder '\summary']; %folder for saving summary data
if ~exist(summaryFolder, 'dir')
    mkdir(summaryFolder)
end
dropboxFolder = ['C:\Users\wilson\Dropbox (HMS)\Data\' exptName '\laganalysis']; %dropbox
if ~exist(dropboxFolder, 'dir')
    mkdir(dropboxFolder)
end

% set variables
totFlies = length(exptFilePath); %number of flies
runThresh = 5; %mm/s
minTimeSpentRunning = 30; %sec, must have run for at least this
onesided=0;


%% pull data and calculate rsquared for each experiment

% initialize
fwd_r=[];
ang_r=[];
sid_r=[];
nChase = 0;

% load in cell tracker
cd([mainFolder '\interpolated'])
load('trackLR.mat')
trackLR=flip(trackLR);

% for each experiment
for f=1:totFlies
    disp(['Analyzing lag data: ' num2str(f) '/' num2str(totFlies)])
    %initialize
    allForward = [];
    allSideway = [];
    allAngular = [];
    allSpikeRt = [];

    % jump to this fly and find all relevant experiment files
    cd(cell2mat(exptFilePath{f}))
    load('metaDat.mat') %pull expt meta
    filebase = [exptInfo.dateDir '_' exptInfo.flyDir ]; %set file name
    allFiles = dir('*pro.mat'); % pull all processed trials
    allFiles(find(contains(string({allFiles.name}), 'Acclimate'))) = []; %remove acclimation trials

    % load in this experiment
    for e = 1:length(allFiles)
        % load in this file
        trialName = allFiles(e).name;
        load(trialName)
        % pool data
        allForward(:,e) = exptData.forwardVelocity;
        allAngular(:,e) = exptData.angularVelocity;
        allSideway(:,e) = exptData.sidewaysVelocity;

        allSpikeRt(:,e) = exptData.spikeRate;
    end
    expttime = exptData.t;

    % normalize left/right
    if trackLR(f)==1
        vCorrect = -1;
    else
        vCorrect = 1;
    end
    allAngular=allAngular.*vCorrect;
    allSideway=allSideway.*vCorrect;
    % if contra speed should be ignored and only look at ipsi
    if onesided
        allForward(allForward<0)=NaN;
        allAngular(allAngular<0)=NaN;
        allSideway(allSideway<0)=NaN;
    end

    % run linear model w/all
    [fwd_model,ang_model,sid_model] = spikerate_v_behaviorlag(allForward,allAngular,allSideway,allSpikeRt,expttime);
    % pull rsquared array for each directional velocity
    fwd_r(:,f,1) = fwd_model.rsquared;
    ang_r(:,f,1) = ang_model.rsquared;
    sid_r(:,f,1) = sid_model.rsquared;

    % save plot
    sgtitle(strrep(filebase,'_',' '))
    cd(saveFolder)
    plotname = ['lag_analysis_' filebase '.png'];
    saveas(gcf,plotname);
    copyfile(plotname, dropboxFolder,'f');

    % threshold for pursuit
    [purForward,purAngular,purSideways,purSpikeRt] = pursuitFinder(allForward,allAngular,allSideway,allSpikeRt,expttime,runThresh);
    % if the fly engaged in pursuit, pull pursuit behavior for analysis
    timeFlySpentRunning = (sum(purForward>runThresh,'all')/length(expttime))*60;
    if timeFlySpentRunning>minTimeSpentRunning
        nChase=nChase+1; %update counter
        % run linear model w/threshold
        [fwd_model,ang_model,sid_model] = spikerate_v_behaviorlag(purForward,purAngular,purSideways,purSpikeRt,expttime);
        fwd_r(:,f,2) = fwd_model.rsquared;
        ang_r(:,f,2) = ang_model.rsquared;
        sid_r(:,f,2) = sid_model.rsquared;

        % save plot
        sgtitle([strrep(filebase,'_',' ') ' pursuit only'])
        cd(saveFolder)
        plotname = ['lag_analysis_thresh_' filebase '.png'];
        saveas(gcf,plotname);
        copyfile(plotname, dropboxFolder,'f');
    else
        % do not include
        fwd_r(:,f,2) = NaN;
        ang_r(:,f,2) = NaN;
        sid_r(:,f,2) = NaN;
    end
end
% pull shift array
shift = fwd_model.shift;
disp('Complete.')


%% find mean lag

% for each curve...
for f=1:totFlies
    % find peak for each velocity
    [f_maxr,f_pk] = findpeaks(fwd_r(:,f,1),'NPeaks',1);
    [a_maxr,a_pk] = findpeaks(ang_r(:,f,1),'NPeaks',1);
    [s_maxr,s_pk] = findpeaks(sid_r(:,f,1),'NPeaks',1);

    % if there is no peak, convert to nan
    if isempty(f_pk)
        f_pk=NaN;
        f_maxr=NaN;
    else
        f_pk=shift(f_pk); %convert to delay time
    end
    if isempty(a_pk)
        a_pk=NaN;
        a_maxr=NaN;
    else
        a_pk=shift(a_pk); %convert to delay time
    end
    if isempty(s_pk)
        s_pk=NaN;
        s_maxr=NaN;
    else
        s_pk=shift(s_pk); %convert to delay time
    end

    % store peak values
    f_pk_t(f)=f_pk;
    f_maxr_all(f)=f_maxr;
    a_pk_t(f)=a_pk;
    a_maxr_all(f)=a_maxr;
    s_pk_t(f)=s_pk;
    s_maxr_all(f)=s_maxr;
end
% repeat with pursuit data
% for each curve...
for f=1:totFlies
    % find peak for each velocity
    [f_maxr,f_pk] = findpeaks(fwd_r(:,f,2),'NPeaks',1);
    [a_maxr,a_pk] = findpeaks(ang_r(:,f,2),'NPeaks',1);
    [s_maxr,s_pk] = findpeaks(sid_r(:,f,2),'NPeaks',1);

    % if there is no peak, convert to nan
    if isempty(f_pk)
        f_pk=NaN;
        f_maxr=NaN;
    else
        f_pk=shift(f_pk); %convert to delay time
    end
    if isempty(a_pk)
        a_pk=NaN;
        a_maxr=NaN;
    else
        a_pk=shift(a_pk); %convert to delay time
    end
    if isempty(s_pk)
        s_pk=NaN;
        s_maxr=NaN;
    else
        s_pk=shift(s_pk); %convert to delay time
    end

    % store peak values
    fp_pk_t(f)=f_pk;
    fp_maxr_all(f)=f_maxr;
    ap_pk_t(f)=a_pk;
    ap_maxr_all(f)=a_maxr;
    sp_pk_t(f)=s_pk;
    sp_maxr_all(f)=s_maxr;
end

% calculate mean "best lag" time for each velocity
f_meanlag = round(mean(f_pk_t,'omitnan'));
a_meanlag = round(mean(a_pk_t,'omitnan'));
s_meanlag = round(mean(s_pk_t,'omitnan'));

fp_meanlag = round(mean(fp_pk_t,'omitnan'));
ap_meanlag = round(mean(ap_pk_t,'omitnan'));
sp_meanlag = round(mean(sp_pk_t,'omitnan'));

%% plot cross correlation curves
disp('Plotting all flies...')

% plot all data
% initialize
figure; set(gcf,'Position',[100 100 1200 500])
colorselect = {'#D95319';'#0072BD';'#7E2F8E'}; %velocity colors
ylabels = {'forward'; 'angular'; 'sideway'}; %velocity names
lw = 1; %linewidth

% plot forward data
subplot(1,3,1)
xline(f_meanlag); hold on
plot(shift,fwd_r(:,:,1),'Color',colorselect{1},'LineWidth',lw)
legend([num2str(f_meanlag) ' msec.'])
ylabel([ylabels{1} ' rsquared'])
% plot angular data
subplot(1,3,2)
xline(a_meanlag); hold on
plot(shift,ang_r(:,:,1),'Color',colorselect{2},'LineWidth',lw)
legend([num2str(a_meanlag) ' msec.'])
ylabel([ylabels{2} ' rsquared'])
xlabel('shift (msec.)')
% plot sideway data
subplot(1,3,3)
xline(s_meanlag); hold on
plot(shift,sid_r(:,:,1),'Color',colorselect{3},'LineWidth',lw)
legend([num2str(s_meanlag) ' msec.'])
ylabel([ylabels{3} ' rsquared'])

% save plot
sgtitle(exptName)
cd(summaryFolder)
plotname = ['summary_lag_analysis_' strrep(exptName,' ','_') '.png'];
saveas(gcf,plotname);
copyfile(plotname, dropboxFolder,'f');


% plot threshold data
% initialize
figure; set(gcf,'Position',[100 100 1200 500])

% plot forward data
subplot(1,3,1)
xline(fp_meanlag); hold on
plot(shift,fwd_r(:,:,2),'Color',colorselect{1},'LineWidth',lw)
legend([num2str(fp_meanlag) ' msec.'])
ylabel([ylabels{1} ' rsquared'])
% plot angular data
subplot(1,3,2)
xline(ap_meanlag); hold on
plot(shift,ang_r(:,:,2),'Color',colorselect{2},'LineWidth',lw)
ylabel([ylabels{2} ' rsquared'])
legend([num2str(ap_meanlag) ' msec.'])
xlabel('shift (msec.)')
% plot sideway data
subplot(1,3,3)
xline(sp_meanlag); hold on
plot(shift,sid_r(:,:,2),'Color',colorselect{3},'LineWidth',lw)
legend([num2str(sp_meanlag) ' msec.'])
ylabel([ylabels{3} ' rsquared'])

% save plot
sgtitle([exptName ' pursuit only'])
cd(summaryFolder)
plotname = ['summary_lag_analysis_thresh_' strrep(exptName,' ','_') '.png'];
saveas(gcf,plotname);
copyfile(plotname, dropboxFolder,'f');

disp('Complete.')

%% plot peaks across flies

% initialize
disp('Plotting r-values across flies...')
figure; set(gcf,'Position',[100 100 1000 500])
ticklabels = {'fwd'; 'ang'; 'side'}; %velocity names

for f=1:totFlies
    % pull data for this fly
    thisFly = [f_maxr_all(f) a_maxr_all(f) s_maxr_all(f)];
    subplot(1,2,1)
    hold on
    plot([1 2 3],thisFly,'-o')
    ylabel('rsquared values')
    xlabel(['all behavior (n=' num2str(totFlies) ')'])
    xticks([1 2 3])
    xticklabels(ticklabels)
    hold off

    % pull data for this fly during pursuit
    thisFlyP = [fp_maxr_all(f) ap_maxr_all(f) sp_maxr_all(f)];
    subplot(1,2,2)
    hold on
    plot([1 2 3],thisFlyP,'-o')
    xlabel(['pursuit behavior only (n=' num2str(nChase) ')'])
    xticks([1 2 3])
    xticklabels(ticklabels)
    hold off
end

% save plot
sgtitle(exptName)
cd(summaryFolder)
plotname = ['summary_lag_acrossflies_' strrep(exptName,' ','_') '.png'];
saveas(gcf,plotname);
copyfile(plotname, dropboxFolder,'f');

disp('Complete.')

%% end
disp('ALL ANALYSES COMPLETE.')
end

