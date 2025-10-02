% preprocess_motionpulse
%
% Preprocess Function
% Pulls all processed files from a single electrophysiology motion pulse
% battery experiment. First organizes the data, then downsamples, and
% plots relevant analyses, including spike rate vs. velocity relationships.
%
% INPUTS
% exptFolder - overarching experiment folder
% thisFly    - file location for data from just this fly
% trialTypes - trial types denoted by motion speed
% trackLR    - which side for AOTU019 w/ 1=L and 2=R
%
% The function processes, downsample, and analyzes data. It accounts for
% side normalization, performs liquid junction potential corrections, and
% generates heatmaps of spike rate vs. directional velocity.
%
% 05/01/2023 - MC created from single battery pipeline
% 07/08/2024 - MC cleaned up and simplified
% 07/18/2024 - MC added voltage processing
%
function preprocess_motionpulse(exptFolder,thisFly,trialTypes,trackLR)
%% initialize
close all

% set folders
dataFolder = char(thisFly{1});
cd(dataFolder)

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
vectorFolder = [plotFolder '\vector']; %for saving plots
if ~exist(vectorFolder, 'dir')
    mkdir(vectorFolder)
end
dropboxFolder = ['C:\Users\wilson\Dropbox (HMS)\Data\' exptFolder '\individual']; %for saving data to dropbox
if ~exist(dropboxFolder, 'dir')
    mkdir(dropboxFolder)
end

% trial types denoted via lettering
nTypes = length(trialTypes);
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
motionCheck = ~contains(trialTypes{1},'0dps');

% load ephys settings
[~,~, settings] = ephysSettings();

% set liquid junction potential subtraction for voltage traces
% source: Gouwens, N. W. & Wilson, R.I. Signal propagation in Drosophila central neurons. J.Neurosci. 29, 6239-6249
ljPotential = 13; %mV

%% load in dataset if not already processed previously
% if the interpolated data set already exists, load in
if exist([intFolder '\' filebase '_int.mat'],'file')
    disp('Data alreaded preprocessed.')
    %else data does not exist so load, downsample, and save
else
    disp('Loading in dataset by trial type...')
    cd(dataFolder)

    % initialize data storage arrays
    allForward = [];
    allSideway = [];
    allAngular = [];
    allPanelPs = [];
    allSpikeRt = [];
    allVoltage = [];
    c = 1;% counter

    % pull data by trial type
    for e = 1:length(allFiles)

        thisTrial = allFiles(e).name;
        % load in this trial data
        load(thisTrial)
        % if first, initialize
        if e==1
            % check if ephys data present
            ephysCheck = isfield(exptData,'spikeRate');
            % check if iinj was delivered (only w/ some)
            iinjCheck = isfield(exptData,'iInj');
            % determine arena midpoint
            midPos = (88 - (exptMeta.objSize/2 - 1))/192 *360; %center position, in degrees
            % normalize directions, such that + always ipsi and - always contra
            % aka flip turn and panel data for cells on the left
            if trackLR=='R'
                flp = 1; %do not flip
            else
                flp = -1; %flip
            end
        end
        % in just a few experiments, iinj was added for some trials
        % ensure those trials are omitted
        if iinjCheck
            iinjThisTrial = mean(exptData.iInj)>10;
        else
            iinjThisTrial = 0;
        end

        if ~iinjThisTrial
            % post-process panel data
            panelData = exptData.g4displayXPos - midPos; %center
            hiddenPos = (184 - (exptMeta.objSize/2))/192 * 360;
            panelData(exptData.g4displayXPos>hiddenPos) = nan; %remove hidden pos
            panelData(abs(diff(panelData))>2) = nan; %remove acquisition jumps
            panelData(isoutlier(panelData)) = nan; %remove any remaining outliers

            % pool this trial data
            allPanelPs(:,c) = panelData.*flp;
            allForward(:,c) = exptData.forwardVelocity;
            allSideway(:,c) = exptData.sidewaysVelocity.*flp;
            allAngular(:,c) = exptData.angularVelocity.*flp;
            if ephysCheck
                % convolve firing rate
                [~,~,spikeRate] = convolveSpikeRate_input(settings,exptData,exptFolder,'gaussian');
                allSpikeRt(:,c) = spikeRate;
                % pool voltage
                allVoltage(:,c) = exptData.scaledVoltage - ljPotential;
            end
            c=c+1; %update counter
        end
    end
    % pull time
    expttime = exptData.t;


    %% interpolate (downsample) dataset
    % optional, but dramatically increases analysis time
    disp('Interpolating dataset...')
    cd(intFolder)

    % set downsampling parameters
    curSamp = size(allPanelPs,1); %total number of current sample points
    newSR = 30; %new sample rate (must be shorter than shortest pixel dwell time)

    % downsample dataset
    int_panelps = round(interp1((1:curSamp),allPanelPs,(1:newSR:curSamp),'linear'),2);
    int_forward = allForward;
    int_angular = allAngular;
    int_sideway = allSideway;
    if isfield(exptData,'spikeRate')
        int_spikert = interp1((1:curSamp),allSpikeRt,(1:newSR:curSamp),'linear');
        int_voltage = interp1((1:curSamp),allVoltage,(1:newSR:curSamp),'linear');
    end
    int_time = interp1((1:curSamp),expttime,(1:newSR:curSamp),'linear');

    % save interpolated dataset w/ephys
    if ephysCheck
        save([filebase '_int.mat'], 'int_panelps', 'int_forward','int_angular','int_sideway','int_spikert','int_voltage','int_time','-v7.3');
    else
        save([filebase '_int.mat'], 'int_panelps', 'int_forward','int_angular','int_sideway','int_time','-v7.3');
    end

    disp('Dataset processed and saved.')

    %% plot spikerate vs directional velocity
    disp('Generating heatmap...')

    % without lag
    [~] = spikert_heatmapvelocity(int_forward,int_angular,int_sideway,int_spikert,int_time,0,1);
    sgtitle(strrep(filebase,'_','/'))
    % save plot
    cd(plotFolder)
    plotname = ['spikert_v_vel_heat_' filebase '.png'];
    saveas(gcf,plotname);
    copyfile(plotname, dropboxFolder,'f');
    % save vectorized plot
    cd(vectorFolder)
    plotname = ['spikert_v_vel_heat_' filebase '.svg'];
    set(gcf,'renderer','Painters')
    saveas(gcf, plotname)
    copyfile(plotname, dropboxFolder,'f');

    % with lag
    [~] = spikert_heatmapvelocity(int_forward,int_angular,int_sideway,int_spikert,int_time,1,1);
    sgtitle([strrep(filebase,'_','/') ' w/lag'])
    % save plot
    cd(plotFolder)
    plotname = ['spikert_v_vel_heatlag_' filebase '.png'];
    saveas(gcf,plotname);
    copyfile(plotname, dropboxFolder,'f');
    % save vectorized plot
    cd(vectorFolder)
    plotname = ['spikert_v_vel_heat_' filebase '.svg'];
    set(gcf,'renderer','Painters')
    saveas(gcf, plotname)
    copyfile(plotname, dropboxFolder,'f');

    disp('Complete.')
end

end

