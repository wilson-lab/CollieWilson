% load_kir_openloop
% Master script for processing behavioral genetics experiments in which 
% the fly was presented with three different experimental conditions: 
% (1) an oscillating target, (2) a motion pulse target, and (3) an optomotor assay.
% This script handles the initialization, data loading, and preprocessing 
% of individual experiments, as well as the pooled analysis pipeline.
%
% Inputs:
%   - userClear: Optional input to clear existing data folders before processing. 
%                If 'y', all relevant folders are deleted and recreated, which is 
%                recommended after changes in the analysis pipeline.
%
% Process Overview:
%   1. Loads all experiments from a designated Excel tracker file, selecting 
%      only included experiments for processing.
%   2. For each experiment, sets file paths for the three conditions: oscillating target, 
%      motion pulse target, and optomotor assay.
%   3. Executes preprocessing for individual experiments (via preprocess_kir_openloop).
%   4. Runs the pooled analysis pipeline for all experiments (via pipeline_kir_openloop).
%
% Outputs:
%   - Processed data is saved into the specified experiment folder.
%
function load_kir_openloop(userClear)
%% load in all experiments of interest

% designate folder to save all processed data in
exptFolder = 'AOTU019 KIR Openloop';
exc_name = 'KIR_Behavior_Tracker.xlsx'; %excel file containing meta data
exc_sheet = 1; %sheet for this experiment

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
nExpt = size(expt_tracker,1);

% for each included expt, pull the full file path
for et = 1:height(expt_tracker)
    day = expt_tracker.Date(et);
    fly = sprintf('fly%02d',expt_tracker.Fly(et));
    oscFile = sprintf('cell%02d',expt_tracker.Oscillate_File(et));
    mopFile = sprintf('cell%02d',expt_tracker.MotionPulse_File(et));
    optFile = sprintf('cell%02d',expt_tracker.Optomotor_File(et));
    % set corresponding filepaths
    exptFilePath{et,1} = fullfile('G:','Behavior Rig','Pursuit',day,fly,[oscFile '_processed']);
    exptFilePath{et,2} = fullfile('G:','Behavior Rig','Pursuit',day,fly,[mopFile '_processed']);
    exptFilePath{et,3} = fullfile('G:','Behavior Rig','Pursuit',day,fly,[optFile '_processed']);

    genotype(et,1) = expt_tracker.Genotype(et);
end

%% run analysis for EACH INDIVIDUAL experiment

% for each included experiment
% parfor nf = 1:nExpt
%     thisFly = exptFilePath(nf,:);
%     thisCondition = genotype{nf};
%     disp(['PROCESSING INDIVIDUAL KIR ' num2str(nf) '/' num2str(nExpt)])
%     preprocess_kir_openloop(exptFolder,thisFly,thisCondition)
% end
% disp('COMPLETE.')

%% run analysis for POOLED experiment
% for all included experiments
pipeline_kir_openloop(exptFolder)
end
