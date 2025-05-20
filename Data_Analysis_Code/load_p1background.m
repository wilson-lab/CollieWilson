% load_p1background
% Master script for processing electrophysiology experiments in which 
% the fly was presented with a blank environment, and optogenetic activation 
% was delivered in alternating 15-second intervals. This script manages 
% the initialization, organization, and preprocessing of individual fly 
% experiments, as well as the pooled analysis pipeline.
%
% Inputs:
%   - thisRun: Specifies the experiment set to run, including '19' for 
%              AOTU019, '25' for AOTU025, and 'a02' for DNa02.
%   - userClear: Optional input to clear existing data folders before 
%                processing. If 'y', all relevant folders are deleted and 
%                recreated, recommended after changes in the analysis pipeline.
%
% Process Overview:
%   1. Loads all included experiments from a designated Excel tracker file 
%      based on the selected experiment type (thisRun).
%   2. Sets file paths for each experiment and loads the main data set.
%   3. Executes preprocessing for individual experiments (via preprocess_p1background).
%   4. Runs the pooled analysis pipeline for all included experiments (via pipeline_p1background).
%
% Outputs:
%   - Processed data is saved into the specified experiment folder.
%
function load_p1background(thisRun,userClear)
%% initialize

% pick which battery analysis to run through pipeline
switch thisRun
    case '19'
        exptFolder = 'AOTU019 Background P1';
        dataFolder = 'E:';
        exc_name = 'AOTU019_Tracker.xlsx'; %excel file containing meta data
        exc_sheet = 5; %sheet for this experiment
    case '25'
        exptFolder = 'AOTU025 Background P1';
        dataFolder = 'E:';
        exc_name = 'AOTU025_Tracker.xlsx'; %excel file containing meta data
        exc_sheet = 3; %sheet for this experiment
    case 'a02'
        exptFolder = 'DNa02 Background P1';
        dataFolder = 'E:';
        exc_name = 'DNa02_Tracker.xlsx'; %excel file containing meta data
        exc_sheet = 3; %sheet for this experiment
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

for nf = 1:nExpt
    thisFly = exptFilePath(nf,:);
    preprocess_p1background(exptFolder,thisFly,trackLR(nf))
end
filename = fullfile('E:',exptFolder,'interpolated','trackLR.mat');
save(filename,'trackLR')


%% run battery analysis for ALL experiments together
% for all included experiments
pipeline_p1background(exptFolder)

end