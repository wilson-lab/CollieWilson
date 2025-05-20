% load_kir_menotaxis
% Master script for processing behavioral genetics experiments in which 
% the fly was presented with a bright target in closed-loop, periodically 
% jumping by 30, 60, or 180 degrees. This script manages the initialization, 
% data loading, and preprocessing for individual fly experiments, as well as 
% the pooled analysis pipeline.
%
% Inputs:
%   - userClear: Optional input to clear existing data folders before processing. 
%                If 'y', all relevant folders are deleted and recreated, which is 
%                recommended after changes in the analysis pipeline.
%
% Process Overview:
%   1. Loads all experiments of interest from a designated Excel tracker file, 
%      selecting only included experiments for processing.
%   2. Sets file paths for each experiment and organizes data for further analysis.
%   3. Executes preprocessing for individual experiments (via preprocess_kir_menotaxis).
%   4. Runs the pooled analysis pipeline for all included experiments (via pipeline_kir_menotaxis).
%
% Outputs:
%   - Processed data is saved into the specified experiment folder.
%
function load_kir_menotaxis(userClear)
%% load in all experiments of interest

% designate folder to save all processed data in
exptFolder = 'AOTU019 KIR Menotaxis';
exc_name = 'KIR_Behavior_Tracker.xlsx'; %excel file containing meta data
exc_sheet = 4; %sheet for this experiment

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
    cellname = sprintf('cell%02d',expt_tracker.Cell(et));
    % set corresponding filepaths
    exptFilePath{et,1} = fullfile('G:','Behavior Rig','Menotaxis',day,fly,[cellname '_processed']);

    genotype(et,1) = expt_tracker.Genotype(et);
end

%% run analysis for EACH INDIVIDUAL experiment
% for each included experiment
% parfor nf = 1:nExpt
%     thisFly = exptFilePath(nf,:);
%     thisGene = genotype{nf};
%     preprocess_kir_menotaxis(exptFolder,thisFly,thisGene)
% end

%% run analysis for POOLED experiment
pipeline_kir_menotaxis(exptFolder)
end
