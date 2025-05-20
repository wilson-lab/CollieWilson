% preprocess_kir_openloop
%
% Pipeline Function
% Pulls all processed files from a single behavior experiment in
% which the fly was presented with 3 different trial sets for assessing the
% effect of KIR on pursuit behavior.
%
% INPUTS
% exptFolder - overarching experiment folder
% thisFly    - file location for data from just this fly
% condition  - specifies which genotype group (GFP, WT, or KIR)
%
% The function processes data from three experimental conditions (oscillatory, 
% motion pulse, and optomotor) for each fly, downsamples the dataset, and 
% saves the processed data. It handles various preprocessing tasks such as 
% normalizing directions and handling trial-specific conditions.
%
% 04/09/2024 - MC adapted from kir menotaxis pipeline
%
function preprocess_kir_openloop(exptFolder,thisFly,condition)
%% initialize
disp('STARTING ANALYSES FOR SELECTED KIR PURSUIT EXPERIMENT...')
close all

% set filename info and create necessary directories
mainfolder = ['E:\' exptFolder];
cd(mainfolder)
intFolder = [mainfolder '\interpolated']; %for saving interpolated data
if ~exist(intFolder, 'dir')
    mkdir(intFolder)
end
plotFolder = [mainfolder '\plot']; %for saving plots
if ~exist(plotFolder, 'dir')
    mkdir(plotFolder)
end
dropboxFolder = ['C:\Users\wilson\Dropbox (HMS)\Data\' exptFolder '\individual']; %for saving data to dropbox
if ~exist(dropboxFolder, 'dir')
    mkdir(dropboxFolder)
end

[~, ~, settings] = ephysSettings();

% set base folders and file
flyFolder = extractBefore(char(thisFly{1}),'\cell');
filebase = [condition '_' strrep(extractAfter(flyFolder,'Rig\'),'\','_')];
cd(flyFolder)

% set trial conditions
exptNames = {"osc";"mopulse";"opto"};
nCon =  size(exptNames,1);

%% load in dataset if not already processed previously
% if the interpolated data set already exists, load in
% for each experiment trial set, load in the dataset, interpolate, and save
for ts = 1:nCon
    % pull trial set info
    thisSet = char(thisFly{ts});
    thisExptName = char(exptNames{ts});
    intFilebase = [filebase '_' thisExptName  '_int.mat'];
    if exist([intFolder '\' intFilebase],'file')
        disp('Data alreaded preprocessed.')
        cd(intFolder)

    else
        %else data does not exist so load, downsample, and save
        disp('Loading in dataset...')
        cd(thisSet)

        % find all files in this directory
        allFiles = dir('*pro.mat');
        % remove any acclimate trials
        allFiles(contains(string({allFiles.name}), 'acclimate')) = [];
        nTrial = length(allFiles);

        % initialize data storage arrays
        allForward = [];
        allSideway = [];
        allAngular = [];
        allPanelPs = [];

        % pull data by trial type
        for e = 1:nTrial
            thisTrial = allFiles(e).name;
            % load in this trial data
            load(thisTrial)

            % pool this trial data
            allForward(:,e) = exptData.forwardVelocity;
            allSideway(:,e) = exptData.sidewaysVelocity;
            allAngular(:,e) = exptData.angularVelocity;

            thisMidpoint = (88 - (exptMeta.objSize/2 - 1))/192 *360; %center position, in degrees
            allPanelPs(:,e) = exptData.g4displayXPos-thisMidpoint;

        end
        % pull time
        expttime = exptData.t;

        % additional post-processing (if necessary)
        switch ts
            case 1 %for oscillatory experiments
            case 2 %for motion pulse experiments
                hiddenPos = 160;
                allPanelPs(allPanelPs>hiddenPos) = nan; %remove hidden pos
                allPanelPs(abs(diff(allPanelPs))>2) = nan; %remove acquisition jumps
                allPanelPs(isoutlier(allPanelPs)) = nan; %remove any remaining outliers
            case 3 %for optomotor experiments
                rndPanelPs = round(allPanelPs);
                % remove points where optomotor sweep was stationary
                allPanelPs(rndPanelPs==max(rndPanelPs)) = nan;
                allPanelPs(rndPanelPs==min(rndPanelPs)) = nan;
        end

        %% interpolate (downsample) dataset
        % optional, but dramatically increases analysis time
        disp('Interpolating dataset...')
        cd(intFolder)

        % set downsampling parameters
        curSamp = size(allPanelPs,1); %total number of current sample points
        newSR = 30; %new sample rate (must be shorter than shortest pixel dwell time)

        % downsample dataset
        int_forward = allForward;
        int_angular = allAngular;
        int_sideway = allSideway;
        int_panelps = interp1((1:curSamp),allPanelPs,(1:newSR:curSamp),'nearest');
        int_time = interp1((1:curSamp),expttime,(1:newSR:curSamp),'linear');

        % save interpolated dataset w/ephys
        save(intFilebase, 'int_forward','int_angular','int_sideway','int_panelps','int_time','-v7.3');

        disp('Dataset processed and saved.')
    end
end
end

