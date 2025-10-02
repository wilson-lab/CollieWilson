% preprocess_battery
%
% Pipeline Function
%
% INPUTS
% exptFolder - overarching experiment folder
% thisFly    - file location for data from just this fly
% condition  - specify which genotype group (GFP, WT, or KIR)
%
% The function organizes and downsamples data, interpolates, estimates
% panel velocity, and performs behavior analyses such as fixation and
% direction change performance. Plots are generated and saved for each step.
%
% 07/24/2025 - MC adapted from kir openloop pipeline
%
function preprocess_aIPg(exptFolder,thisFly,trialTypes)
%% initialize
close all

% set filename info and create necessary directories
folder = generateFolders(exptFolder);
folder.plot = folder.main;
folder.plot = [folder.main '\plot']; %for saving individual plots
if ~exist(folder.plot, 'dir')
    mkdir(folder.plot)
end
kirFolder = char(thisFly{1});
cd(kirFolder)
% pull all file info
allFiles = dir('*pro.mat');
% remove any acclimate trials
allFiles(contains(string({allFiles.name}), 'Acclimate')) = [];
nTrial = 10; %first 10
% pull expt meta
load('metaDat.mat')
filebase = strrep(char(extractBetween(char(thisFly{1}),'aIPg\','_processed')),'\','_');

% load processing settings
settings = processSettings();
nObj = length(trialTypes); % fetch number of visual objects
newSR = 30; %new sample rate (must be shorter than shortest pixel dwell time)

%% load in dataset if not already processed previously

% if the interpolated data set already exists, load in
if exist([folder.int '\' filebase '_int.mat'],'file')
    disp('Dataset alread pre-processed.')
else
    disp('Loading in dataset...')
    cd(kirFolder)

    % initialize data storage arrays
    allForward = [];
    allSideway = [];
    allAngular = [];
    allPanelPs = [];
    allJumpTrg = [];

    % pull data for each trial
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

    %% organize dataset by condition
    % conditions repeated every multiple of nObj
    % set z dimension to correspond to each gain condition
    for g = 1:nObj
        orgForward(:,:,g) = allForward(:,g:nObj:nTrial);
        orgSideway(:,:,g) = allSideway(:,g:nObj:nTrial);
        orgAngular(:,:,g) = allAngular(:,g:nObj:nTrial);
        orgPanelPs(:,:,g) = allPanelPs(:,g:nObj:nTrial);
        orgJumpTrg(:,:,g) = allJumpTrg(:,g:nObj:nTrial);
    end

    %% interpolate (downsample) dataset
    % optional, but dramatically increases analysis time
    disp('Interpolating dataset...')
    cd(folder.int)

    % set downsampling parameters
    curSamp = size(orgPanelPs,1); %total number of current sample points

    % downsample dataset
    int_forward = orgForward;
    int_angular = orgAngular;
    int_sideway = orgSideway;
    int_panelps = interp1((1:curSamp),orgPanelPs,(1:newSR:curSamp),'nearest');
    int_jumptrg = interp1((1:curSamp),orgJumpTrg,(1:newSR:curSamp),'nearest');
    int_time = interp1((1:curSamp),expttime,(1:newSR:curSamp),'linear');

    % estimate panel velocity
    int_panelvel = computePanelVelocity(int_panelps,int_time);

    % save interpolated dataset w/ephys
    save([filebase '_int.mat'], 'int_forward','int_angular','int_sideway','int_panelps','int_panelvel','int_jumptrg','int_time','-v7.3');

    disp('Dataset processed and saved.')

    %% plot fixation performance
    thisFixation = fixationFinder(int_panelps,int_forward,int_time,1);

    % Find all axes in the current figure
    axList = findall(gcf, 'type', 'axes');
    axList = flip(axList);

    % Loop through and update ylabel
    n = 1;
    for g = 1:nObj
        ylabel(axList(n), trialTypes{g});
        n=n+nObj;
    end

    % save plot
    cd(folder.plot)
    plotname = ['fixation_' filebase];
    saveas(gcf,[plotname '.png']);

    disp('Complete.')
    accFound = 0;

    if ~exist([folder.accl '\' filebase '_acc.mat'], 'file')
        % Check if the raw acclimation data exists in kirFolder
        if exist(fullfile(kirFolder, 'preExptAcclimate_1_pro.mat'), 'file')
            disp('Loading pre-experiment acclimation data from kirFolder...');
            load(fullfile(kirFolder, 'preExptAcclimate_1_pro.mat'));
            accFound = 1;

            % If not in kirFolder, check the folder with cell01_processed
        else
            % Define the alternate folder path by replacing any 'cellXX_processed' with 'cell01_processed'
            altFolder = regexprep(kirFolder, 'cell\d+_processed', 'cell01_processed');
            if exist(fullfile(altFolder, 'preExptAcclimate_1_pro.mat'), 'file')
                disp('Loading pre-experiment acclimation data from cell01_processed folder...');
                load(fullfile(altFolder, 'preExptAcclimate_1_pro.mat'));
                accFound = 1;
            else
                disp('No pre-experiment acclimation file found in either cell01_processed or kirFolder.');
                accFound = 0;
            end
        end

        % Proceed with data extraction and downsampling if data was loaded
        if accFound
            % Extract variables
            accl_forwardVelocity = exptData.forwardVelocity;
            accl_angularVelocity = exptData.angularVelocity;
            accl_panelPs = exptData.g4displayXPos - 180;
            accl_time = exptData.t;

            % Downsample data
            curSamp = length(accl_panelPs); % Total number of current sample points
            int_accl_forwardVelocity = accl_forwardVelocity;
            int_accl_angularVelocity = accl_angularVelocity;
            int_accl_panelPs = interp1((1:curSamp), accl_panelPs, (1:newSR:curSamp), 'nearest');
            int_accl_time = interp1((1:curSamp), accl_time, (1:newSR:curSamp), 'linear');

            % Save interpolated data
            save([folder.accl '\' filebase '_acc.mat'], 'int_accl_forwardVelocity', 'int_accl_angularVelocity', 'int_accl_panelPs', 'int_accl_time', '-v7.3');
            disp('Acclimation data processed and saved.');
        end
    end
end
end

