% pipeline_menotaxis_single
%
% Pipeline Function
% Pulls all processed files from a single electrophysiology experiment in
% which fly was in closed loop with a bright object that jumped to the left
% or right every 60 seconds
%
% INPUTS
% exptFolder- overarching experiment folder
% thisFly   - file location for data from just this fly
% trackLR   - which side for AOTU019 w/ 1=L and 2=R
% runThresh - forward speed
%
% 11/14/2023 - MC adapted from ihold pipeline
%

function pipeline_menotaxis_single(exptFolder,thisFly,trackLR,runThresh)
%% initialize
disp('STARTING ANALYSES FOR SELECTED MENOTAXIS EXPERIMENT...')
close all

% set folders
menotaxFolder = char(thisFly{1});
cd(menotaxFolder)

% pull all file info
allFiles = dir('*pro.mat');
% remove any acclimate trials
allFiles(find(contains(string({allFiles.name}), 'Acclimate'))) = [];
% pull expt meta
load('metaDat.mat')

% set filename info and create necessary directories
filebase = [exptInfo.dateDir '_' exptInfo.flyDir ];
mainfolder = ['E:\' exptFolder];
cd(mainfolder)
intFolder = [mainfolder '\interpolated']; %for saving interpolated data
if ~exist(intFolder, 'dir')
    mkdir(intFolder)
end
plotFolder = [mainfolder '\plot']; %for saving plots
if ~exist(plotFolder, 'dir')
    mkdir(plotFolder)
end
dropboxFolder = ['C:\Users\wilson\Dropbox (HMS)\Data\' exptFolder '\individual']; %for saving data to dropbox
if ~exist(dropboxFolder, 'dir')
    mkdir(dropboxFolder)
end

[~, ~, settings] = ephysSettings();


%% load in dataset if not already processed previously

% if the interpolated data set already exists, load in
if exist([intFolder '\' filebase '_int.mat'],'file')
    disp('Existing dataset loading...')
    cd(intFolder)
    load([filebase '_int.mat'])
    disp('Complete.')

    %else data does not exist so load, downsample, and save
else
    disp('Loading in dataset...')
    cd(menotaxFolder)

    % initialize data storage arrays
    allForward = [];
    allSideway = [];
    allAngular = [];
    allSpikeRt = [];
    allPanelPs = [];

    % if trackLR == 'R'
    %     vnorm = 1; %no flip
    % else
    %     vnorm = -1; %flip
    % end

    % pull data by trial type
    for e = 1:length(allFiles)

        thisTrial = allFiles(e).name;
        % load in this trial data
        load(thisTrial)

        % pool this trial data
        allForward(:,e) = exptData.forwardVelocity;
        allSideway(:,e) = exptData.sidewaysVelocity;
        allAngular(:,e) = exptData.angularVelocity;

        allPanelPs(:,e) = exptData.g4displayXPos;
        allJumpTrig(:,e) = exptData.pythonJumpTrig;
        
        % convolve spikerate
        ksize= 1000; %kernel size
        [~,~,spikeRate] = convolveSpikeRate_input(settings,exptData,ksize,'gaussian');
        allSpikeRt(:,e) = spikeRate;

    end
    % pull time
    expttime = exptData.t;

    %% interpolate (downsample) dataset
    % optional, but dramatically increases analysis time
    disp('Interpolating dataset...')
    cd(intFolder)

    % set downsampling parameters
    curSamp = size(allSpikeRt,1); %total number of current sample points
    newSR = 10; %new sample rate (must be shorter than shortest pixel dwell time)

    % downsample dataset
    int_forward = interp1((1:curSamp),allForward,(1:newSR:curSamp),'linear');
    int_angular = interp1((1:curSamp),allAngular,(1:newSR:curSamp),'linear');
    int_sideway = interp1((1:curSamp),allSideway,(1:newSR:curSamp),'linear');
    int_panelps = interp1((1:curSamp),allPanelPs,(1:newSR:curSamp),'nearest');
    int_jumptrg = interp1((1:curSamp),allJumpTrig,(1:newSR:curSamp),'nearest');
    int_spikert = interp1((1:curSamp),allSpikeRt,(1:newSR:curSamp),'linear');

    int_time = interp1((1:curSamp),expttime,(1:newSR:curSamp),'linear');

    % save interpolated dataset w/ephys
    save([filebase '_int.mat'], 'int_forward','int_angular','int_sideway','int_panelps','int_jumptrg','int_spikert','int_time','-v7.3');

    disp('Dataset processed and saved.')
end


%% visualize behavior

%initialize
nTrial = size(int_forward,2);
figure; set(gcf,'Position',[100 100 1700 900])
% plot a max of 10 trials
if nTrial>10
    n1 = nTrial-9;
    nmax = 10;
else
    n1 = 1;
    nmax = nTrial;
end
tiledlayout(nmax,1,'TileSpacing','compact')

% for each trial
for e = n1:nTrial
    nexttile
    colororder({'#77AC30','#D95319'})
    yyaxis left
    plot(int_time,int_panelps(:,e),'Color','#77AC30')
    yline(0,':')
    ylabel('HD')
    axis tight
    yyaxis right
    plot(int_time,int_forward(:,e),'Color','#D95319')
    ylabel('fwd')
    axis tight
end
sgtitle(strrep(filebase,'_','/'))
xlabel('time (sec)')

% save plot
cd(plotFolder)
plotname = ['raw_data_' filebase '.png'];
saveas(gcf,plotname);
copyfile(plotname, dropboxFolder,'f');


%% plot spikerate vs velocity
disp('Analyzing spikerate vs velocity...')

% plot pursuit behavior
summaryData = spikerate_v_velocityPlot(int_forward,int_angular,int_sideway,int_time,int_spikert,runThresh,1,1);

% save plot
cd(plotFolder)
sgtitle([strrep(filebase,'_','/') ' ' exptFolder])
plotname = ['spikert_v_velocity_' filebase '.png'];
saveas(gcf,plotname);
copyfile(plotname, dropboxFolder,'f');

disp('Complete.')

%% find trials when the fly was actually running well
[run_forward,~,~,run_panelps] = pursuitFinder(int_forward,int_angular,int_sideway,int_panelps,int_time,runThresh);
[~,~,~,run_spikert] = pursuitFinder(int_forward,int_angular,int_sideway,int_spikert,int_time,runThresh);

% min time spent running per trial
minTimeSpentRunning = 6*60; %sec
% determine if this trial type met minimum requirements to be considered "running"
timeSpentRunning = int_time(sum(~isnan(run_forward))+1); %sec
runLog = (timeSpentRunning>minTimeSpentRunning);


%% plot heading distribution and corresponding vector strength
disp('Analyzing heading distribution...')

% set the number of heading bins
headingBins = 25;

% run for all data
figure; set(gcf,'Position',[100 100 450 800])
tiledlayout(2,1,'TileSpacing','compact')
nexttile
h1 = polarhistogram(deg2rad(int_panelps),headingBins,'FaceColor',"#77AC30",'FaceAlpha',0.5,'Normalization','probability');
[h1_max,h1_imax] = max(h1.Values); %find most prominent heading
h1_i = (h1.BinEdges(h1_imax)+h1.BinEdges(h1_imax+1))/2; %find prominent heading bin center
hold on
polarplot([h1_i;h1_i],[0;h1_max],'Color',"#7E2F8E",'LineWidth',2)

sgtitle(strrep(filebase,'_','/'))
title('all behavior')

% run for only running data
if sum(runLog)
    nexttile
    h2 = polarhistogram(deg2rad(run_panelps),headingBins,'FaceColor',"#77AC30",'FaceAlpha',0.5,'Normalization','probability');
    [h2_max,h2_imax] = max(h2.Values); %find most prominent heading
    h2_i = (h2.BinEdges(h2_imax)+h2.BinEdges(h2_imax+1))/2; %find prominent heading bin center
    hold on
    polarplot([h2_i;h2_i],[0;h2_max],'Color',"#7E2F8E",'LineWidth',2)

    title('run behavior only')

    goalHD = round(rad2deg(h2_i));
    disp(['Goal HD = ' num2str(goalHD) 'deg'])
end

% save plot
cd(plotFolder)
plotname = ['heading_dist_' filebase '.png'];
saveas(gcf,plotname);
copyfile(plotname, dropboxFolder,'f');

disp('Complete.')


%% plot spikerate vs heading
disp('Analyzing spikerate vs heading...')

% run analysis
[~] = spikerate_v_heading(int_panelps,int_spikert,1);
sgtitle([strrep(filebase,'_','/') ' all'])
% save plot
cd(plotFolder)
plotname = ['spikerate_v_heading_' filebase '_all' '.png'];
saveas(gcf,plotname);
copyfile(plotname, dropboxFolder,'f');

% run for only running data
if sum(runLog)
    % run analysis
    [~] = spikerate_v_heading(run_panelps,run_spikert,1);
    sgtitle([strrep(filebase,'_','/') ' all'])
    % save plot
    cd(plotFolder)
    plotname = ['spikerate_v_heading_' filebase '_run' '.png'];
    saveas(gcf,plotname);
    copyfile(plotname, dropboxFolder,'f');
end

disp('Complete.')

%% plot behavior and spikerate changes at each bar jump

% plot behavior and spikerate around ALL jumps
disp('Analyzing bar jumps...')
[~] = barjump_nocorrect_analysis(int_panelps,int_jumptrg,int_forward,int_angular,int_spikert,int_time,15,1);
sgtitle(strrep(filebase,'_','/'))
% save plot
cd(plotFolder)
plotname = ['jump_aligned_all_' filebase '.png'];
saveas(gcf,plotname);
copyfile(plotname, dropboxFolder,'f');

% plot behavior and spikerate around corrected jumps
[~,restoreFT] = barjump_analysis(int_panelps,int_jumptrg,int_forward,int_angular,int_spikert,int_time,15,1);
sgtitle(strrep(filebase,'_','/'))
% save plot
cd(plotFolder)
plotname = ['jump_aligned_' filebase '.png'];
saveas(gcf,plotname);
copyfile(plotname, dropboxFolder,'f');

% plot restore features
figure; set(gcf,'Position',[100 100 450 800])
tiledlayout(2,1,'TileSpacing','compact')
nexttile
% plot individual trials
plot(1:6,restoreFT.all,'Color',[0.8 0.8 0.8],'Marker','o','LineStyle','none')
hold on
% plot mean +/- std
errorbar(1:6,restoreFT.mean,restoreFT.std,'k','Marker','o','LineStyle','none')
axis padded
xlabel('jump size (deg)')
ylabel('correction time (s)')
xticklabels({'+30','-30','+60','-60','+180','-180'})
% plot number of corrected vs uncorrected trials
nexttile
bh = bar(1:6,[restoreFT.trialsWithoutCorrect;restoreFT.trialsWithCorrect],'stacked','FaceColor','Flat');
bh(1).CData = [0.8 0.8 0.8];
bh(2).CData = "#77AC30";
axis padded
xlabel('jump size (deg)')
ylabel('correction rate')
xticklabels({'+30','-30','+60','-60','+180','-180'})
legend('Uncorrected','Corrected', 'Location', 'southoutside')
sgtitle(strrep(filebase,'_','/'))
% save plot
cd(plotFolder)
plotname = ['jump_corrected_' filebase '.png'];
saveas(gcf,plotname);
copyfile(plotname, dropboxFolder,'f');


%% end
disp('ALL ANALYSES FOR THIS BATTERY ARE COMPLETE.')
end

