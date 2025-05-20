% load_targetpulse
% Master script for processing electrophysiology experiments in which 
% the fly was presented with a battery of target pulses, varying the 
% speed of the target (0, 25, or 75 deg/s). This script handles the 
% initialization, organization, and pre-processing for individual experiments, 
% as well as running the pooled analysis pipeline.
%
% Inputs:
%   - thisRun: Specifies the experiment set to run, including target motion 
%              and control conditions. Options include '19_motion', '19_none', 
%              '19_stationary', '25_motion', '25_none', 'a02_motion', and 'a02_none'.
%   - userClear: Optional input to clear existing data folders before processing. 
%                If 'y', all relevant folders are deleted and recreated, recommended 
%                after changes in the analysis pipeline.
%
% Process Overview:
%   1. Loads included experiments from a designated Excel tracker file based on the 
%      selected battery analysis (thisRun).
%   2. Sets file paths and prepares data for each included experiment.
%   3. Executes preprocessing for individual experiments (via preprocess_motionpulse).
%   4. Runs the pooled analysis pipeline for all included experiments (via pipeline_motion_pulse).
%
% Outputs:
%   - Processed data is saved into the specified experiment folder.
%
function load_targetpulse(thisRun,userClear)
%% initialize
% pick which battery analysis to run through pipeline based on input
switch thisRun
    case '19_motion'
        exptFolder = 'AOTU019 Motion Pulse';
        dataFolder = 'E:';
        exc_name = 'AOTU019_Tracker.xlsx'; %excel file containing meta data
        exc_sheet = 2; %sheet for this experiment
        trialTypes = {'75deg/sec','25deg/sec'};
    case '19_none'
        exptFolder = 'AOTU019 Motion Pulse No P1';
        dataFolder = 'E:';
        exc_name = 'AOTU019_Tracker.xlsx';
        exc_sheet = 3;
        trialTypes = {'75deg/sec','25deg/sec'};
    case '19_stationary'
        exptFolder = 'AOTU019 Stationary Pulse';
        dataFolder = 'E:';
        exc_name = 'AOTU019_Tracker.xlsx';
        exc_sheet = 4;
        trialTypes = {'0dps'};
        
    case '25_motion'
        exptFolder = 'AOTU025 Motion Pulse';
        dataFolder = 'E:';
        exc_name = 'AOTU025_Tracker.xlsx'; %excel file containing meta data
        exc_sheet = 1; %sheet for this experiment
        trialTypes = {'75deg/sec','25deg/sec'};
    case '25_none'
        exptFolder = 'AOTU025 Motion Pulse No P1';
        dataFolder = 'E:';
        exc_name = 'AOTU025_Tracker.xlsx'; %excel file containing meta data
        exc_sheet = 2; %sheet for this experiment
        trialTypes = {'75deg/sec','25deg/sec'};
    case '25_stationary'
        exptFolder = 'AOTU025 Stationary Pulse';
        dataFolder = 'E:';
        exc_name = 'AOTU025_Tracker.xlsx';
        exc_sheet = 5;
        trialTypes = {'0dps'};

    case 'a02_motion'
        exptFolder = 'DNa02 Motion Pulse';
        dataFolder = 'E:';
        exc_name = 'DNa02_Tracker.xlsx'; %excel file containing meta data
        exc_sheet = 1; %sheet for this experiment
        trialTypes = {'25deg/sec'};
    case 'a02_none'
        exptFolder = 'DNa02 Motion Pulse No P1';
        dataFolder = 'E:';
        exc_name = 'DNa02_Tracker.xlsx'; %excel file containing meta data
        exc_sheet = 2; %sheet for this experiment
        trialTypes = {'25deg/sec'};
end

% load in all experiments of interest
% designate folder to save all processed data in
saveFolder = 'E:\';
dropboxFolder = 'C:\Users\wilson\Dropbox (HMS)\Data';
cd(saveFolder)

% asks if user would like to clear the experiment folder
% note: this is recommended any time the analysis pipeline is modified
if contains(userClear,'y')
    %clear data folder
    cd(saveFolder)
    if exist(exptFolder, 'dir')
        rmdir(fullfile(saveFolder,exptFolder),'s') %remove
    end
    mkdir(exptFolder) %remake
else
    cd(saveFolder) %else, jump to data folder
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
    exptFilePath{et,1} = fullfile(dataFolder,day,fly,[cell '_processed']);
    trackLR(et) = expt_tracker.Side{et};
end
nExpt = height(exptFilePath);


%% run battery analysis for EACH INDIVIDUAL experiment
% pre-process motion pulse data
parfor nf = 1:nExpt
    thisFly = exptFilePath(nf,:);
    disp(['pre-processing ' num2str(nf) '/' num2str(nExpt)])
    preprocess_motionpulse(exptFolder,thisFly,trialTypes,trackLR(nf))
end
filename = fullfile('E:',exptFolder,'interpolated','trackLR.mat');
save(filename,'trackLR')


%% run battery analysis for ALL experiments together
% run for all included experiments
pipeline_motion_pulse(exptFolder,trialTypes)

end