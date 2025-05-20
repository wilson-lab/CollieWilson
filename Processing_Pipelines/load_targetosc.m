% load_targetosc
% Master script for processing electrophysiology experiments in which 
% the fly was presented with a battery of oscillatory target stimuli, 
% varying the speed of the target (15-75 deg/s). This script handles the 
% initialization, organization, and preprocessing for individual experiments, 
% as well as running the pooled analysis pipeline.
%
% Inputs:
%   - userClear: Optional input to clear existing data folders before processing. 
%                If 'y', all relevant folders are deleted and recreated, which is 
%                recommended after any changes in the analysis pipeline.
%
% Process Overview:
%   1. Loads all included experiments from a designated Excel tracker file 
%      based on oscillatory target trials.
%   2. Sets file paths for each experiment, including both main and stim test 
%      data (if available).
%   3. Executes preprocessing for individual experiments (via preprocess_oscillate).
%   4. Runs the pooled analysis pipeline for all included experiments (via pipeline_oscillate).
%
% Outputs:
%   - Processed data is saved into the specified experiment folder.
%
function load_targetosc(userClear)
%% initialize
% pick which battery analysis to run through pipeline
exptFolder = 'AOTU019 Oscillate';
trialTypes = {'no stim','15deg/sec','25deg/sec','35deg/sec','55deg/sec','75deg/sec'};
exc_name = 'AOTU019_Tracker.xlsx'; %excel file containing meta data
exc_sheet = 1; %sheet for this experiment

% load in all experiments of interest
% designate folder to save all processed data in
dataFolder = 'E:\';
dropboxFolder = 'C:\Users\wilson\Dropbox (HMS)\Data';
cd(dataFolder)

% asks if user would like to clear the experiment folder
% note: this is recommended any time the analysis pipeline is modified
if contains(userClear,'y')
    %clear dropbox folder
    cd(dropboxFolder)
    if exist(exptFolder, 'dir')
        rmdir(fullfile(dropboxFolder,exptFolder),'s') %remove
    end
    mkdir(exptFolder) %remake

    %clear data folder
    cd(dataFolder)
    if exist(exptFolder, 'dir')
        rmdir(fullfile(dataFolder,exptFolder),'s') %remove
    end
    mkdir(exptFolder) %remake
else
    cd(dataFolder) %else, jump to data folder
end
% load in tracker for this experiment set
tracker = readtable(exc_name,'Sheet',exc_sheet); %load all expts
expt_tracker = tracker(tracker.Include==1,:); %select only included expts

% for each included expt, pull the full file path
for et = 1:height(expt_tracker)
    day = expt_tracker.Date(et);
    fly = sprintf('fly%02d',expt_tracker.Fly(et));

    % load in main data set
    cell = sprintf('cell%02d',expt_tracker.Cell(et));
    exptFilePath{et,1} = fullfile('E:',day,fly,[cell '_processed']);

    % load in stim test data set, if available
    if expt_tracker.TestCell(et)
        testcell = sprintf('cell%02d',expt_tracker.TestCell(et));
        exptFilePath{et,2} = fullfile('E:',day,fly,[testcell '_processed']);
    else
        exptFilePath{et,2} = 0;
    end
    trackLR(et) = expt_tracker.Side{et};
end
nExpt = height(exptFilePath);

%% run battery analysis for EACH INDIVIDUAL experiment
% for each included experiment
% parfor nf = 1:nExpt
%     thisFly = exptFilePath(nf,:);
%     preprocess_oscillate(exptFolder,thisFly,trialTypes,trackLR(nf))
% end
% filename = fullfile('E:',exptFolder,'interpolated','trackLR.mat');
% save(filename,'trackLR')


%% run battery analysis for POOLED experiments
% for all included experiments
pipeline_oscillate(exptFolder,trialTypes)
end
