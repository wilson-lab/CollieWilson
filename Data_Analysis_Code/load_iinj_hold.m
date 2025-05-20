% load_iinj_ihold
% Master script for processing electrophysiology experiments in which 
% the fly walked in darkness, and the cell was either recorded normally or 
% held below spike threshold via current injection. This script manages the 
% initialization, data loading, and preprocessing for individual experiments, 
% as well as the pooled analysis pipeline.
%
% Inputs:
%   - userClear: Optional input to clear existing data folders before processing. 
%                If 'y', all relevant folders are deleted and recreated, which is 
%                recommended after changes in the analysis pipeline.
%
% Process Overview:
%   1. Loads all experiments from a designated Excel tracker file, selecting 
%      only included experiments for processing.
%   2. Sets file paths for each experiment and organizes data for analysis.
%   3. Executes preprocessing for individual experiments (via preprocess_iinj_hold).
%   4. Runs the pooled analysis pipeline for all experiments (via pipeline_iinj_hold).
%
% Outputs:
%   - Processed data is saved into the specified experiment folder.
%
function load_iinj_hold(userClear)
%% load in all experiments of interest
% pick which battery analysis to run through pipeline
exptFolder = 'AOTU019 IInj Hold';
exc_name = 'AOTU019_Tracker.xlsx'; %excel file containing meta data
exc_sheet = 7; %sheet for this experiment

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
nExpt = size(expt_tracker,1);

% for each included expt, pull the full file path
for et = 1:height(expt_tracker)
    day = expt_tracker.Date(et);
    fly = sprintf('fly%02d',expt_tracker.Fly(et));
    cell = sprintf('cell%02d',expt_tracker.Cell(et));
    exptFilePath{et,1} = fullfile('E:',day,fly,[cell '_processed']);
    trackLR(et) = expt_tracker.Side{et};
end
filename = fullfile('E:',exptFolder,'trackLR.mat');
save(filename,'trackLR')

%% run analysis for EACH INDIVIDUAL experiment

% for each included experiment
% parfor nf = 1:nExpt
%     thisFly = exptFilePath(nf,:);
%     preprocess_iinj_hold(exptFolder,thisFly,trackLR(nf))
% end
% filename = fullfile('E:',exptFolder,'interpolated','trackLR.mat');
% save(filename,'trackLR')

%% run analysis for POOLED experiment
% for all included experiments
pipeline_iinj_hold(exptFolder)
end
