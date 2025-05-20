% preprocess_iinj_hold
%
% Pipeline Function
% Pulls all processed files from a single electrophysiology experiment in
% which held current pulses were delivered with or without P1-activation.
%
% INPUTS
% exptFolder - overarching experiment folder
% thisFly    - file location for data from just this fly
% trackLR    - which side for AOTU019 w/ 1 = L and 2 = R
%
% The function processes the electrophysiology data, including spike rate, 
% current injection, and voltage, interpolates the dataset, and extracts 
% timepoints where current pulses were delivered. The processed data is saved 
% for further analysis.
%
% 08/31/2023 - MC adapted from iinj pipeline
% 07/08/2024 - MC cleaned up and simplified
%
function preprocess_iinj_hold(exptFolder,thisFly,trackLR)
%% initialize

% set folders
iinjFolder = char(thisFly{1});
cd(iinjFolder)

% pull all file info
allFiles = dir('*pro.mat');
% remove any acclimate trials
allFiles(find(contains(string({allFiles.name}), 'Acclimate'))) = [];
% pull expt meta
load('metaDat.mat')

% set filename info and create necessary directories
filebase = [exptInfo.dateDir '_' exptInfo.flyDir ];
folder = generateFolders(exptFolder);
folder.main = ['E:\' exptFolder];
cd(folder.main)
folder.int = [folder.main '\interpolated']; %for saving interpolated data
if ~exist(folder.int, 'dir')
    mkdir(folder.int)
end

% load relevant settings
settings = processSettings();
[~, ~, e_settings] = ephysSettings();


%% load in dataset if not already processed previously
% if the interpolated data set already exists, load in
if exist([folder.int '\' filebase '_int.mat'],'file')
    disp('Existing dataset loading...')
    cd(folder.int)
    load([filebase '_int.mat'])
    disp('Complete.')

    %else data does not exist so load, downsample, and save
else
    disp('Loading in dataset by trial type...')
    cd(iinjFolder)

    % initialize data storage arrays
    allForward = [];
    allSideway = [];
    allAngular = [];
    allSpikeRt = [];
    allIInject = [];

    if trackLR == 'R'
        vnorm = 1; %no flip
    else
        vnorm = -1; %flip
    end

    % pull data by trial type
    for e = 1:length(allFiles)

        thisTrial = allFiles(e).name;
        % pull this trial info
        thisN = str2double(thisTrial(30:31)); %pull trial number
        thisSub = find(ismember(settings.letters,thisTrial)); %pull trial type based on letter denoted
        % load in this trial data
        load(thisTrial)

        % pool this trial data
        allForward(:,thisN,thisSub) = exptData.forwardVelocity;
        allSideway(:,thisN,thisSub) = exptData.sidewaysVelocity.*vnorm;
        allAngular(:,thisN,thisSub) = exptData.angularVelocity.*vnorm;
        allVoltage(:,thisN,thisSub) = exptData.scaledVoltage - settings.ljp;
        allIInject(:,thisN,thisSub) = exptData.iInj;

        [~,~,spikeRate] = convolveSpikeRate_input(e_settings,exptData,exptFolder,'gaussian');
        allSpikeRt(:,thisN,thisSub) = spikeRate;

    end
    % pull time
    expttime = exptData.t;

    %% interpolate (downsample) dataset
    % optional, but dramatically increases analysis time
    disp('Interpolating dataset...')
    cd(folder.int)

    % set downsampling parameters
    curSamp = size(allSpikeRt,1); %total number of current sample points
    newSR = 30; %new sample rate (must be shorter than shortest pixel dwell time)

    % downsample dataset
    int_forward = allForward;
    int_angular = allAngular;
    int_sideway = allSideway;
    int_iinject = interp1((1:curSamp),allIInject,(1:newSR:curSamp),'nearest');
    int_voltage = interp1((1:curSamp),allVoltage,(1:newSR:curSamp),'linear');
    int_spikert = interp1((1:curSamp),allSpikeRt,(1:newSR:curSamp),'linear');

    int_time = interp1((1:curSamp),expttime,(1:newSR:curSamp),'linear');

    %% pull only timepoints where ihold was delivered

    % find trials where current was delivered
    iinj_trials = find(sum(ischange(int_iinject(:,1,:))));

    % find start/stop indices for when current was delivered
    iinj_idx = find(ischange(int_iinject(:,1,iinj_trials(1))));
    iinj_rng = iinj_idx(1):iinj_idx(2)-1;

    % pull data accordingly
    int_forward = int_forward(iinj_rng,:,:);
    int_angular = int_angular(iinj_rng,:,:);
    int_sideway = int_sideway(iinj_rng,:,:);
    int_iinject = int_iinject(iinj_rng,:,:);
    int_voltage = int_voltage(iinj_rng,:,:);
    int_spikert = int_spikert(iinj_rng,:,:);

    int_time = int_time(1:(length(iinj_rng)));

    % save interpolated dataset w/ephys
    save([filebase '_int.mat'], 'int_forward','int_angular','int_sideway','int_iinject','int_spikert','int_voltage','int_time','-v7.3');

    disp('Dataset processed and saved.')
end
end

