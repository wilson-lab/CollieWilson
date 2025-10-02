% preprocess_kir_closedloop
%
% Pipeline Function
% Pulls all processed files from a single behavior experiment in which 
% the fly was presented with a target moving in closed-loop at various gains.
% The function processes and analyzes the data, including fixation 
% performance and direction change behavior.
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
% 08/21/2024 - MC adapted from kir openloop pipeline
%
function preprocess_kir_closedloop(exptFolder,thisFly,condition)
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
nTrial = length(allFiles);
% pull expt meta
load('metaDat.mat')
filebase = [condition '_' strrep(char(extractBetween(char(thisFly{1}),'Pursuit\','_processed')),'\','_')];

% load processing settings
settings = processSettings();
if exptFolder == 'AOTU025 KIR'
    nGain = 1;
else
    nGain = length(settings.pursuitGain); % fetch number of gain settings
end

newSR = 30; %new sample rate (must be shorter than shortest pixel dwell time)

%% load in dataset if not already processed previously

% if the interpolated data set already exists, load in
if exist([folder.int '\' filebase '_int.mat'],'file')
    disp('Data already pre-processed.')
    %load([folder.int '\' filebase '_int.mat'])
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
    % conditions repeated every multiple of nGain
    % set z dimension to correspond to each gain condition

    % note that low gain conditions (60 and 85x) conditions were remove in
    % the second half of data collection, so they are omitted
    checkGain = contains(thisTrial,'G','IgnoreCase',false);
    if checkGain
        % low gain conditions included
         for g = 3:nGain+2
            orgForward(:,:,g-2) = allForward(:,g:nGain+2:nTrial);
            orgSideway(:,:,g-2) = allSideway(:,g:nGain+2:nTrial);
            orgAngular(:,:,g-2) = allAngular(:,g:nGain+2:nTrial);
            orgPanelPs(:,:,g-2) = allPanelPs(:,g:nGain+2:nTrial);
            orgJumpTrg(:,:,g-2) = allJumpTrg(:,g:nGain+2:nTrial);
         end
    else
        % no low gain conditions
        for g = 1:nGain
            orgForward(:,:,g) = allForward(:,g:nGain:nTrial);
            orgSideway(:,:,g) = allSideway(:,g:nGain:nTrial);
            orgAngular(:,:,g) = allAngular(:,g:nGain:nTrial);
            orgPanelPs(:,:,g) = allPanelPs(:,g:nGain:nTrial);
            orgJumpTrg(:,:,g) = allJumpTrg(:,g:nGain:nTrial);
        end
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
    fix_panelps = thisFixation.panelps_run;
    fix_angular = int_angular;
    fix_angular(~thisFixation.idx_run) = nan;

    thisFixationT = reshape(sum(int_time(sum(~isnan(fix_panelps))+1)),1,nGain); %time spent fixating, s

    % save plot
    cd(folder.plot)
    plotname = ['fixation_' filebase];
    saveas(gcf,[plotname '.png']);

    disp('Complete.')
end

%% plot direction change performance
% if all(thisFixationT>settings.minFixationTime)
%     [~, ~, ~, ~] = setpoint_dirchange(fix_panelps, fix_angular, int_jumptrg, int_time, 1);
%     sgtitle(strrep(filebase,'_',' '))
%     % save plot
%     cd(folder.plot)
%     plotname = ['dirchange_' filebase];
%     saveas(gcf,[plotname '.png']);
% end

%% plot peak change performance
% if all(thisFixationT>settings.minFixationTime)
%     [~, ~] = setpoint_max_to_max(fix_panelps, fix_angular, int_jumptrg, int_time, 1);
%     sgtitle(strrep(filebase,'_',' '))
%     % save plot
%     cd(folder.plot)
%     plotname = ['pkchange_' filebase];
%     saveas(gcf,[plotname '.png']);
% end

%% Load and process pre-experiment acclimation data
% Check if the interpolated acclimation data already exists
if ~exist([folder.accl '\' filebase '_acc.mat'], 'file')
    % Check if the raw acclimation data exists in kirFolder
    if exist(fullfile(kirFolder, 'preExptAcclimate_1_pro.mat'), 'file')
        disp('Loading pre-experiment acclimation data from kirFolder...');
        load(fullfile(kirFolder, 'preExptAcclimate_1_pro.mat'));
    
    % If not in kirFolder, check the folder with cell01_processed
    else
        % Define the alternate folder path by replacing any 'cellXX_processed' with 'cell01_processed'
        altFolder = regexprep(kirFolder, 'cell\d+_processed', 'cell01_processed');
        if exist(fullfile(altFolder, 'preExptAcclimate_1_pro.mat'), 'file')
            disp('Loading pre-experiment acclimation data from cell01_processed folder...');
            load(fullfile(altFolder, 'preExptAcclimate_1_pro.mat'));
        else
            disp('No pre-experiment acclimation file found in either cell01_processed or kirFolder.');
        end
    end
    
    % Proceed with data extraction and downsampling if data was loaded
    if exist('exptData', 'var')
        % Extract variables
        accl_forwardVelocity = exptData.forwardVelocity;
        accl_panelPs = exptData.g4displayXPos - 180;
        accl_time = exptData.t;

        % Downsample data
        curSamp = length(accl_panelPs); % Total number of current sample points
        int_accl_forwardVelocity = accl_forwardVelocity;
        int_accl_panelPs = interp1((1:curSamp), accl_panelPs, (1:newSR:curSamp), 'nearest');
        int_accl_time = interp1((1:curSamp), accl_time, (1:newSR:curSamp), 'linear');

        % Save interpolated data
        save([folder.accl '\' filebase '_acc.mat'], 'int_accl_forwardVelocity', 'int_accl_panelPs', 'int_accl_time', '-v7.3');
        disp('Acclimation data processed and saved.');
    end
end


end

