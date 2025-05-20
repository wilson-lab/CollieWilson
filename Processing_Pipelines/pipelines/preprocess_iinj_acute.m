% preprocess_iinj_acute
%
% Pipeline Function
% Pulls all processed files from a single electrophysiology experiment in
% which current pulses were delivered with or without P1-activation.
%
% INPUTS
% exptFolder - overarching experiment folder
% thisFly    - file location for data from just this fly
% trackLR    - which side for AOTU019 w/ 1 = L and 2 = R
%
% The function processes the electrophysiology data, including spike rate,
% current injection, and voltage traces. It interpolates the dataset,
% downsamples it, and saves the processed data for further analysis.
%
% 07/10/2023 - MC adapted from no stimulus pipeline
% 06/03/2024 - MC switched SR kernel to causal
% 07/08/2024 - MC cleaned up and simplified
%
function preprocess_iinj_acute(exptFolder,thisFly,trackLR)
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

% trial types denoted via lettering
nTypes = 2;

% load processing settings
settings = processSettings();
[~, ~, esettings] = ephysSettings();


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
    allVoltage = [];
    allIInject = [];
    trialN = 1; %initialize trial counter

    if trackLR == 'R'
        vnorm = 1; %no flip
    else
        vnorm = -1; %flip
    end

    % pull data by trial type
    for e = 1:length(allFiles)

        thisTrial = allFiles(e).name;
        % pull this trial info
        thisSub = find(ismember(settings.letters,thisTrial)); %pull trial type based on letter denoted
        % load in this trial data
        load(thisTrial)

        % pool this trial data
        allForward(:,trialN,thisSub) = exptData.forwardVelocity;
        allSideway(:,trialN,thisSub) = exptData.sidewaysVelocity.*vnorm;
        allAngular(:,trialN,thisSub) = exptData.angularVelocity.*vnorm;
        allVoltage(:,trialN,thisSub) = exptData.scaledVoltage - settings.ljp;
        allIInject(:,trialN,thisSub) = exptData.iInj;

        [~,~,spikeRate] = convolveSpikeRate_input(esettings,exptData,exptFolder,'causal');
        allSpikeRt(:,trialN,thisSub) = spikeRate;

        % update counter as needed
        if ~rem(e,2)
            trialN = trialN+1;
        end

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
    int_voltage = interp1((1:curSamp),allVoltage,(1:newSR:curSamp),'linear');
    int_iinject = interp1((1:curSamp),allIInject,(1:newSR:curSamp),'nearest');
    int_spikert = interp1((1:curSamp),allSpikeRt,(1:newSR:curSamp),'linear');

    int_time = interp1((1:curSamp),expttime,(1:newSR:curSamp),'linear');

    % save interpolated dataset w/ephys
    save([filebase '_int.mat'], 'int_forward','int_angular','int_sideway','int_voltage','int_iinject','int_spikert','int_time','-v7.3');

    disp('Dataset processed and saved.')
end
end

