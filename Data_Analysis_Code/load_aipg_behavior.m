% load_targetbattery_behavior
% Master script for processing behavioral experiments.
%
% Inputs:
%   - userClear: Optional input to clear existing experiment data folders
%                before proceeding. If 'y', all relevant folders will be
%                deleted and recreated. This is recommended when changes
%                have been made to the analysis pipeline.
%
% Process Overview:
%   1. Loads all included experiments from an Excel tracker file.
%   2. Prepares individual data for each experiment by setting file paths.
%   3. Executes preprocessing for each experiment.
%   4. Runs pooled analysis across all experiments.
%
% Outputs:
%   - Processed data is saved into the designated experiment folder.
%
function load_aipg_behavior(thisRun,userClear)
%% load in all experiments of interest

dataFolder = 'E:\';
dropboxFolder = 'C:\Users\wilson\Dropbox (HMS)\Data';
cd(dataFolder)

% designate folder to save all processed data in
switch thisRun
    case 'aIPg'
        exptFolder = 'aIPg Pursuit Behavior';
        exc_name = 'KIR_Behavior_Tracker.xlsx'; %excel file containing meta data
        exc_sheet = 8; %sheet for this experiment

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
            dataFile = sprintf('cell%02d',expt_tracker.Folder(et));
            % set corresponding filepaths
            exptFilePath{et,1} = fullfile('G:','Behavior Rig','aIPg',day,fly,[dataFile '_processed']);
        end

        %% run analysis for EACH INDIVIDUAL experiment

        %for each included experiment
        trialTypes = {'na'};
        for nf = 1:nExpt
            thisFly = exptFilePath{nf,:};
            disp(['PROCESSING ' num2str(nf) '/' num2str(nExpt)])
            preprocess_aIPg(exptFolder,thisFly,trialTypes)
        end
        disp('COMPLETE.')

    case 'P1'
        % Data pre-processed elsewhere, just running same analyses for ease
        exptFolder = 'P1 Pursuit Behavior';
        trialTypes = {'na'};
end

%% run analysis for POOLED experiment
% for all included experiments
pipeline_aipg_behavior(exptFolder,trialTypes)

end
