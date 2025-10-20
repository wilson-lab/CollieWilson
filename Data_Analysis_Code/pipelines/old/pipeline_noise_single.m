% pipeline_noise_single
%
% Pipeline Function
% Pulls all processed files from a single behavior battery
% experiment. First organizes the data based on trial and trial type, then
% performs a series of analyses and plots.
%
% INPUTS
% exptFolder- overarching experiment folder
% thisFly   - file location for data from just this fly
% condition - specify which genotype group (GFP, WT, or KIR)
% trialTypes - notation for each trial type in the battery
%
% 06/07/2024 - MC created from visual pursuit pipeline
%

function pipeline_noise_single(exptFolder,thisFly,condition,trialTypes)
%% initialize
disp('STARTING ANALYSES FOR SELECTED PURSUIT EXPT...')
close all

% set folders
dataFolder = char(thisFly{1});
cd(dataFolder)

% pull all file info
allFiles = dir('*pro.mat');
% remove any acclimate trials
allFiles(find(contains(string({allFiles.name}), 'acclimate'))) = [];
% pull expt meta
load('metaDat.mat')

% set filename info and create necessary directories
filebase = [condition '_' exptInfo.dateDir '_' exptInfo.flyDir ];
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
    disp('Loading in dataset by trial type...')
    cd(dataFolder)

    % initialize data storage arrays
    allForward = [];
    allSideway = [];
    allAngular = [];
    allPanelPs = [];

    % pull data by trial type
    for e = 1:length(allFiles)

        thisTrial = allFiles(e).name;
        % pull this trial info
        thisN = str2double(thisTrial(30:31)); %pull trial number
        if nTypes>1
            thisA = find(ismember(letters,thisTrial)); %pull trial type based on letter denoted
        else
            thisA = 1;
        end
        % load in this trial data
        load(thisTrial)
        % pool this trial data
        allPanelPs(:,thisN,thisA) = exptData.g4displayXPos;
        allForward(:,thisN,thisA) = exptData.forwardVelocity;
        allSideway(:,thisN,thisA) = exptData.sidewaysVelocity;
        allAngular(:,thisN,thisA) = exptData.angularVelocity;

    end
    % pull time
    expttime = exptData.t;

    % center panel data at zero
    % for each trial type
    for tt = 1:nTypes
        thisMidpoint = min(allPanelPs(1e3:end,1,tt)) + (max(allPanelPs(1e3:end,1,tt))-min(allPanelPs(1e3:end,1,tt)))/2;
        allPanelPs_c(:,:,tt) = allPanelPs(:,:,tt)-thisMidpoint;
    end

    %% interpolate (downsample) dataset
    % optional, but dramatically increases analysis time
    disp('Interpolating dataset...')
    cd(intFolder)

    % set downsampling parameters
    curSamp = size(allPanelPs,1); %total number of current sample points
    newSR = 10; %new sample rate (must be shorter than shortest pixel dwell time)

    % downsample dataset
    int_panelps = round(interp1((1:curSamp),allPanelPs_c,(1:newSR:curSamp),'linear'),2);
    int_forward = interp1((1:curSamp),allForward,(1:newSR:curSamp),'linear');
    int_angular = interp1((1:curSamp),allAngular,(1:newSR:curSamp),'linear');
    int_sideway = interp1((1:curSamp),allSideway,(1:newSR:curSamp),'linear');
    int_time = interp1((1:curSamp),expttime,(1:newSR:curSamp),'linear');

    % save interpolated dataset
    save([filebase '_int.mat'], 'int_panelps', 'int_forward','int_angular','int_sideway','int_time','-v7.3');

    disp('Dataset processed and saved.')
end

%% end
disp('ALL ANALYSES FOR THIS BATTERY ARE COMPLETE.')
end

