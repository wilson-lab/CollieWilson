% preprocess_kir_menotaxis
%
% Pipeline Function
% Pulls all processed files from a single behavior experiment in
% which the fly was in closed loop with a bright object that jumped to the left
% or right every 60 seconds. The function processes and analyzes the data
% to assess behavior and bar-jump responses for flies with or without KIR expression.
%
% INPUTS
% exptFolder      - overarching experiment folder
% thisFly         - file location for data from just this fly
% condition       - specify which genotype group (GFP, WT, or KIR)
% settings.runThreshB - forward speed threshold to determine "running"
%
% The function loads, processes, and interpolates data, performs behavior
% analysis, including heading distributions, bar-jump responses, and 
% running behavior. It also visualizes the results and saves both raw and 
% processed data.
%
% 02/27/2024 - MC adapted from menotaxis ephys pipeline
%
function preprocess_kir_menotaxis(exptFolder,thisFly,condition)
%% initialize

% set folders
kirFolder = char(thisFly{1});
cd(kirFolder)

% pull all file info
allFiles = dir('*pro.mat');
% remove any acclimate trials
allFiles(find(contains(string({allFiles.name}), 'acclimate'))) = [];
nTrial = length(allFiles);
% pull expt meta
load('metaDat.mat')

% set filename info and create necessary directories
filebase = [condition '_' exptInfo.dateDir '_' exptInfo.flyDir];
folder = generateFolders(exptFolder);
folder.plot = [folder.main '\plot']; %for saving plots
if ~exist(folder.plot, 'dir')
    mkdir(folder.plot)
end

% load processing settings
settings = processSettings();

%% load in dataset if not already processed previously

% if the interpolated data set already exists, load in
if exist([folder.int '\' filebase '_int.mat'],'file')
    disp('Existing dataset loading...')
    cd(folder.int)
    load([filebase '_int.mat'])
    disp('Complete.')

    %else data does not exist so load, downsample, and save
else
    disp('Loading in dataset...')
    cd(kirFolder)

    % initialize data storage arrays
    allForward = [];
    allSideway = [];
    allAngular = [];
    allPanelPs = [];
    allJumpTrg = [];

    % pull data by trial type
    for e = 1:nTrial

        thisTrial = allFiles(e).name;
        % load in this trial data
        load(thisTrial)

        % pool this trial data
        allForward(:,e) = exptData.forwardVelocity;
        allSideway(:,e) = exptData.sidewaysVelocity;
        allAngular(:,e) = exptData.angularVelocity;

        allPanelPs(:,e) = exptData.g4displayXPos-180;
        allJumpTrg(:,e) = exptData.pythonJumpTrig;

    end
    % pull time
    expttime = exptData.t;

    %% interpolate (downsample) dataset
    % optional, but dramatically increases analysis time
    disp('Interpolating dataset...')
    cd(folder.int)

    % set downsampling parameters
    curSamp = size(allPanelPs,1); %total number of current sample points
    newSR = 30; %new sample rate (must be shorter than shortest pixel dwell time)

    % downsample dataset
    int_forward = allForward;
    int_angular = allAngular;
    int_sideway = allSideway;
    int_panelps = interp1((1:curSamp),allPanelPs,(1:newSR:curSamp),'nearest');
    int_jumptrg = interp1((1:curSamp),allJumpTrg,(1:newSR:curSamp),'nearest');

    int_time = interp1((1:curSamp),expttime,(1:newSR:curSamp),'linear');

    % save interpolated dataset w/ephys
    save([filebase '_int.mat'], 'int_forward','int_angular','int_sideway','int_panelps','int_jumptrg','int_time','-v7.3');

    disp('Dataset processed and saved.')
end


%% visualize behavior and threshold

%initialize
figure; set(gcf,'Position',[100 100 1700 900])
tiledlayout(10,1,'TileSpacing','compact')
% plot a max of 10 trials
if nTrial>10
    n1 = nTrial-9;
else
    n1 = 1;
end

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
cd(folder.plot)
plotname = ['raw_data_' filebase '.png'];
saveas(gcf,plotname);
copyfile(plotname, folder.dropbox,'f');

%% find trials when the fly was actually running well
[run_forward,~,~,run_panelps] = pursuitFinder(int_forward,int_angular,int_sideway,int_panelps,int_time,settings.runThreshB);

% min time spent running per trial
minTimeSpentRunning = 6*60; %sec
% determine if this trial type met minimum requirements to be considered "running"
timeSpentRunning = int_time(sum(~isnan(run_forward))+1); %sec
runLog = (timeSpentRunning>minTimeSpentRunning);


%% plot heading distribution and corresponding vector strength
disp('Analyzing heading distribution...')

% run for all data
figure; set(gcf,'Position',[100 100 450 800])
tiledlayout(2,1,'TileSpacing','compact')
nexttile
h1 = polarhistogram(deg2rad(int_panelps),settings.HDBins,'FaceColor',settings.HDColor,'FaceAlpha',0.5,'Normalization','probability');
[h1_max,h1_imax] = max(h1.Values); %find most prominent heading
h1_i = (h1.BinEdges(h1_imax)+h1.BinEdges(h1_imax+1))/2; %find prominent heading bin center
hold on
polarplot([h1_i;h1_i],[0;h1_max],'Color',"#7E2F8E",'LineWidth',2)

sgtitle(strrep(filebase,'_','/'))
title('all behavior')

% run for only running data
if sum(runLog)
    nexttile
    h2 = polarhistogram(deg2rad(run_panelps),settings.HDBins,'FaceColor',settings.HDColor,'FaceAlpha',0.5,'Normalization','probability');
    [h2_max,h2_imax] = max(h2.Values); %find most prominent heading
    h2_i = (h2.BinEdges(h2_imax)+h2.BinEdges(h2_imax+1))/2; %find prominent heading bin center
    hold on
    polarplot([h2_i;h2_i],[0;h2_max],'Color',"#7E2F8E",'LineWidth',2)

    title('run behavior only')

    goalHD = round(rad2deg(h2_i));
    disp(['Goal HD = ' num2str(goalHD) 'deg'])
end

% save plot
cd(folder.plot)
plotname = ['heading_dist_' filebase '.png'];
saveas(gcf,plotname);
copyfile(plotname, folder.dropbox,'f');

disp('Complete.')

%% plot behavior at each bar jump

% plot behavior around jumps
disp('Analyzing bar jumps...')
[~,restoreFT] = barjump_analysis(int_panelps,int_jumptrg,int_forward,int_angular,0,int_time,1);
sgtitle(strrep(filebase,'_','/'))
% save plot
cd(folder.plot)
plotname = ['jump_aligned_' filebase '.png'];
saveas(gcf,plotname);
copyfile(plotname, folder.dropbox,'f');

% plot restore features
figure; set(gcf,'Position',[100 100 450 800])
tiledlayout(2,1,'TileSpacing','compact')
nexttile
% plot individual trials
plot(1:6,restoreFT.all,'Color',settings.trialColor,'Marker','o','LineStyle','none')
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
bh(1).CData = settings.trialColor;
bh(2).CData = "#77AC30";
axis padded
xlabel('jump size (deg)')
ylabel('correction rate')
xticklabels({'+30','-30','+60','-60','+180','-180'})
legend('Uncorrected','Corrected', 'Location', 'southoutside')
sgtitle(strrep(filebase,'_','/'))
% save plot
cd(folder.plot)
plotname = ['jump_corrected_' filebase '.png'];
saveas(gcf,plotname);
copyfile(plotname, folder.dropbox,'f');

end

