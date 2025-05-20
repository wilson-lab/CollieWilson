% preprocess_oscillate
%
% Pipeline Function
% Pulls all processed files from a single electrophysiology battery
% experiment. First organizes the data based on trial and trial type, then
% performs a series of analyses and plots.
%
% INPUTS
% exptFolder - overarching experiment folder
% thisFly    - file location for data from just this fly
% trialTypes - list of trial types (e.g., different oscillation speeds)
% trackLR    - which side for AOTU019 w/ 1 = L and 2 = R
%
% This function loads, processes, and interpolates data, dividing it by
% trial type and side. It also normalizes directions and saves the
% processed dataset for further analysis.
%
% 02/03/2023 - MC created from visual pursuit pipeline
% 07/08/2024 - MC cleaned and simplified
%
function preprocess_oscillate(exptFolder,thisFly,trialTypes,trackLR)
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

% plotting
srlimit = 20;

% min time spent running
minTimeSpentRunning = 10; %sec

% trial types denoted via lettering
nSpeeds = length(trialTypes);
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

[~, ~, settings] = ephysSettings();


%% load in dataset if not already processed previously

% if the interpolated data set already exists, load in
if exist([intFolder '\' filebase '_int.mat'],'file')
    disp('Existing dataset loading...')
    cd(intFolder)
    load([filebase '_int.mat'])
    ephysCheck = exist('int_spikert');
    disp('Complete.')

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

    % pull data by trial type
    for e = 1:length(allFiles)

        thisTrial = allFiles(e).name;
        % pull this trial info
        thisN = str2double(thisTrial(30:31)); %pull trial number
        if nSpeeds>1
            thisA = find(ismember(letters,thisTrial)); %pull trial type based on letter denoted
        else
            thisA = 1;
        end
        % load in this trial data
        load(thisTrial)
        % if first time, check if ephys data present
        if e==1
            ephysCheck = isfield(exptData,'spikeRate');
        end
        % normalize directions, such that + always ipsi and - always contra
        % aka flip turn and panel data for cells on the left
        if trackLR=='R'
            flp = 1; %do not flip
        else
            flp = -1; %flip
        end

        % pool this trial data
        allPanelPs(:,thisN,thisA) = exptData.g4displayXPos.*flp;
        allForward(:,thisN,thisA) = exptData.forwardVelocity;
        allSideway(:,thisN,thisA) = exptData.sidewaysVelocity.*flp;
        allAngular(:,thisN,thisA) = exptData.angularVelocity.*flp;
        if ephysCheck
            [~,~,spikeRate] = convolveSpikeRate_input(settings,exptData,exptFolder,'gaussian');
            allSpikeRt(:,thisN,thisA) = spikeRate;
        end

    end
    % pull time
    expttime = exptData.t;

    % center panel data at zero
    % for each trial type
    for tt = 1:nSpeeds
        thisMidpoint = min(allPanelPs(1e3:end,1,tt)) + (max(allPanelPs(1e3:end,1,tt))-min(allPanelPs(1e3:end,1,tt)))/2;
        allPanelPs_c(:,:,tt) = allPanelPs(:,:,tt)-thisMidpoint;
    end

    %% interpolate (downsample) dataset
    % optional, but dramatically increases analysis time
    disp('Interpolating dataset...')
    cd(intFolder)

    % set downsampling parameters
    curSamp = size(allPanelPs,1); %total number of current sample points
    newSR = 30; %new sample rate (must be shorter than shortest pixel dwell time)

    % downsample dataset
    int_panelps = round(interp1((1:curSamp),allPanelPs_c,(1:newSR:curSamp),'linear'),2);
    int_forward = allForward;
    int_angular = allAngular;
    int_sideway = allSideway;
    if isfield(exptData,'spikeRate')
        int_spikert = interp1((1:curSamp),allSpikeRt,(1:newSR:curSamp),'linear');
    end
    int_time = interp1((1:curSamp),expttime,(1:newSR:curSamp),'linear');

    % save interpolated dataset w/ephys
    if ephysCheck
        save([filebase '_int.mat'], 'int_panelps', 'int_forward','int_angular','int_sideway','int_spikert','int_time','-v7.3');
    else
        save([filebase '_int.mat'], 'int_panelps', 'int_forward','int_angular','int_sideway','int_time','-v7.3');
    end

    disp('Dataset processed and saved.')
end
end

