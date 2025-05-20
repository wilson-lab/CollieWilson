% preprocess_p1background
%
% Pipeline Function
% Pulls all processed files from a single electrophysiology experiment in
% which no visual stimulus was given, and P1 was only activated for portions
% of the trial to investigate activity correlates with behavior.
%
% INPUTS
% exptFolder - overarching experiment folder
% thisFly    - file location for data from just this fly
% trackLR    - which side for AOTU019 w/ 1 = L and 2 = R
%
% The function processes electrophysiology data, normalizes for side 
% differences, interpolates the dataset for downsampling, divides it based on 
% P1 light on vs off, and saves the results.
%
% 02/03/2023 - MC created from visual pursuit pipeline
% 07/08/2024 - MC cleaned and simplified
%
function preprocess_p1background(exptFolder,thisFly,trackLR)
%% initialize
close all

% set folders
pulseFolder = char(thisFly{1});
cd(pulseFolder)

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

[~, ~, settings] = ephysSettings();
settingsP = processSettings();
ljPotential = settingsP.ljp;

%% load in dataset if not already processed previously

% if the interpolated data set already exists, load in
if exist([intFolder '\' filebase '_pulse.mat'],'file')
    disp('Existing pulse dataset.')
    cd(intFolder)
    %load([filebase '_pulse.mat'])
    %disp('Complete.')

    %else data does not exist so load, downsample, and save
else
    disp('Loading in pulse dataset...')
    cd(pulseFolder)

    % pull all file info
    pulseFiles = dir('*pro.mat');
    % remove any acclimate trials
    pulseFiles(find(contains(string({pulseFiles.name}), 'Acclimate'))) = [];

    % initialize data storage arrays
    pulseForward = [];
    pulseSideway = [];
    pulseAngular = [];
    pulseSpikeRt = [];
    pulseVoltage = [];

    % normalize directions, such that + always ipsi and - always contra
    % aka flip turn and panel data for cells on the left
    if trackLR=='R'
        flp = 1; %do not flip
    else
        flp = -1; %flip
    end

    % pull data by trial type
    for e = 1:length(pulseFiles)

        thisTrial = pulseFiles(e).name;
        % load in this trial data
        load(thisTrial)

        % pool pulse data
        pulseForward(:,e) = exptData.forwardVelocity;
        pulseSideway(:,e) = exptData.sidewaysVelocity.*flp;
        pulseAngular(:,e) = exptData.angularVelocity.*flp;

        [~,~,spikeRate] = convolveSpikeRate_input(settings,exptData,exptFolder,'gaussian');
        pulseSpikeRt(:,e) = spikeRate;
        pulseVoltage(:,e) = exptData.scaledVoltage - ljPotential;

    end
    % pull opto stim
    pulseOptoStm = exptData.optoStim;
    % pull time
    expttime = exptData.t;

    %% interpolate (downsample) dataset
    % optional, but dramatically increases analysis time
    disp('Interpolating dataset...')
    cd(intFolder)

    % set downsampling parameters
    curSamp = size(pulseOptoStm,1); %total number of current sample points
    newSR = 30; %new sample rate

    % downsample dataset
    int_forward_full = pulseForward;
    int_angular_full = pulseAngular;
    int_sideway_full = pulseSideway;
    int_spikert_full = interp1((1:curSamp),pulseSpikeRt,(1:newSR:curSamp),'linear');
    int_voltage_full = interp1((1:curSamp),pulseVoltage,(1:newSR:curSamp),'linear');
    int_optostm = interp1((1:curSamp),pulseOptoStm,(1:newSR:curSamp),'nearest');
    int_time = interp1((1:curSamp),expttime,(1:newSR:curSamp),'linear');

    %% divide dataset based on light on vs off

    % initialize
    int_forward = [];
    int_angular = [];
    int_sideway = [];
    int_spikert = [];

    % find stimulus onset/offsets
    stimonset = find(diff(int_optostm)>0)+1;
    restonset = [1 find(diff(int_optostm)<0)+1];
    stimduration = stimonset(1)-restonset(1)-1;

    % set buffer to exclude
    bufferSize = 0.250; %s, size buffer att onset/offset
    [bufferIdx] = fetchTimeIdx(int_time,bufferSize); %find nearest index
    
    % for each stim on/off
    onIdx = [];
    offIdx = [];
    for s = 1:length(stimonset)
        offIdx = [offIdx restonset(s)+bufferIdx:restonset(s)+stimduration];
        onIdx = [onIdx stimonset(s)+bufferIdx: stimonset(s)+stimduration];
    end

    % pull on points
    int_forward(:,:,1) = int_forward_full(onIdx,:);
    int_angular(:,:,1) = int_angular_full(onIdx,:);
    int_sideway(:,:,1) = int_sideway_full(onIdx,:);
    int_spikert(:,:,1) = int_spikert_full(onIdx,:);
    int_voltage(:,:,1) = int_voltage_full(onIdx,:);
    % pull off points
    int_forward(:,:,2) = int_forward_full(offIdx,:);
    int_angular(:,:,2) = int_angular_full(offIdx,:);
    int_sideway(:,:,2) = int_sideway_full(offIdx,:);
    int_spikert(:,:,2) = int_spikert_full(offIdx,:);
    int_voltage(:,:,2) = int_voltage_full(offIdx,:);

    % save interpolated dataset
    save([filebase '_pulse.mat'], 'int_forward','int_angular','int_sideway','int_spikert','int_voltage','int_time','-v7.3');

    disp('Dataset processed and saved.')
end

%% end
end

