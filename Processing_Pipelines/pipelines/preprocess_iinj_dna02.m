% preprocess_iinj_dna02
%
% Preprocess Function
% Pulls all processed files from a single electrophysiology motion pulse
% battery experiment. First organizes the data, then downsamples it, and
% plots relevant analyses including spike rate and voltage.
%
% INPUTS
% exptFolder - overarching experiment folder
% thisFly    - file location for data from just this fly
% trackLR    - which side for DNa02 neuron, 1 = L and 2 = R
%
% The function processes, interpolates, and downsamples data based on current
% injection conditions. It handles voltage traces, firing rates, and behavioral 
% variables, while normalizing directions and excluding hidden positions.
%
% 05/01/2023 - MC created from single battery pipeline
% 07/08/2024 - MC cleaned up and simplified
% 07/18/2024 - MC added voltage processing
%
function preprocess_iinj_dna02(exptFolder,thisFly,trackLR)
%% initialize
close all

% set folders
dataFolder = char(thisFly{1});
cd(dataFolder)

% pull all file info
allFiles = dir('*pro.mat');
% remove any acclimate trials
allFiles(find(contains(string({allFiles.name}), 'Acclimate'))) = [];
% pull expt meta
load('metaDat.mat')

% set filename info and create necessary directories
filebase = [exptInfo.dateDir '_' exptInfo.flyDir ];
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
vectorFolder = [plotFolder '\vector']; %for saving plots
if ~exist(vectorFolder, 'dir')
    mkdir(vectorFolder)
end
dropboxFolder = ['C:\Users\wilson\Dropbox (HMS)\Data\' exptFolder '\individual']; %for saving data to dropbox
if ~exist(dropboxFolder, 'dir')
    mkdir(dropboxFolder)
end

% load ephys settings
[~,~, settings] = ephysSettings();

% set liquid junction potential subtraction for voltage traces
% source: Gouwens, N. W. & Wilson, R.I. Signal propagation in Drosophila central neurons. J.Neurosci. 29, 6239-6249
ljPotential = 13; %mV

%% load in dataset if not already processed previously
% if the interpolated data set already exists, load in
if exist([intFolder '\' filebase '_int.mat'],'file')
    disp('Data alreaded preprocessed.')

    %else data does not exist so load, downsample, and save
else
    disp('Loading in dataset by trial type...')
    cd(dataFolder)

    % initialize data storage arrays
    allForward = [];
    allSideway = [];
    allAngular = [];
    allPanelPs = [];
    allSpikeRt = [];
    allVoltage = [];

    % pull data by trial type
    for e = 1:length(allFiles)

        thisTrial = allFiles(e).name;
        % load in this trial data
        load(thisTrial)
        % if first, initialize
        if e==1
            % check if ephys data present
            ephysCheck = isfield(exptData,'spikeRate');
            % determine arena midpoint
            midPos = (88 - (exptMeta.objSize/2 - 1))/192 *360; %center position, in degrees
        end
        
        % post-process panel data
        panelData = exptData.g4displayXPos - midPos; %center
        hiddenPos = (184 - (exptMeta.objSize/2))/192 * 360;
        panelData(exptData.g4displayXPos>hiddenPos) = nan; %remove hidden pos
        panelData(abs(diff(panelData))>2) = nan; %remove acquisition jumps
        panelData(isoutlier(panelData)) = nan; %remove any remaining outliers

        % for each trial, check if current was delivered
        if mean(exptData.iInj)<50 % no current added
            c = 1;
        else %depolarizing current added
            c = 2;
        end
        
        % normalize directions, such that + always ipsi and - always contra
        % aka flip turn and panel data for cells on the left
        if trackLR=='R'
            flp = 1; %do not flip
        else
            flp = -1; %flip
        end
        % pool this trial data
        allPanelPs(:,e,c) = panelData.*flp;
        allForward(:,e,c) = exptData.forwardVelocity;
        allSideway(:,e,c) = exptData.sidewaysVelocity.*flp;
        allAngular(:,e,c) = exptData.angularVelocity.*flp;
        if ephysCheck
            % convolve firing rate
            [~,~,spikeRate] = convolveSpikeRate_input(settings,exptData,exptFolder,'gaussian');
            allSpikeRt(:,e,c) = spikeRate;
            % pool voltage
            allVoltage(:,e,c) = exptData.scaledVoltage - ljPotential;
        end

    end
    % pull time
    expttime = exptData.t;


    %% interpolate (downsample) dataset
    % optional, but dramatically increases analysis time
    disp('Interpolating dataset...')
    cd(intFolder)

    % set downsampling parameters
    curSamp = size(allPanelPs,1); %total number of current sample points
    newSR = 30; %new sample rate (must be shorter than shortest pixel dwell time)

    % downsample dataset
    int_panelps = round(interp1((1:curSamp),allPanelPs,(1:newSR:curSamp),'linear'),2);
    int_forward = allForward;
    int_angular = allAngular;
    int_sideway = allSideway;
    if isfield(exptData,'spikeRate')
        int_spikert = interp1((1:curSamp),allSpikeRt,(1:newSR:curSamp),'linear');
        int_voltage = interp1((1:curSamp),allVoltage,(1:newSR:curSamp),'linear');
    end
    int_time = interp1((1:curSamp),expttime,(1:newSR:curSamp),'linear');

    % save interpolated dataset w/ephys
    if ephysCheck
        save([filebase '_int.mat'], 'int_panelps', 'int_forward','int_angular','int_sideway','int_spikert','int_voltage','int_time','-v7.3');
    else
        save([filebase '_int.mat'], 'int_panelps', 'int_forward','int_angular','int_sideway','int_time','-v7.3');
    end

    disp('Dataset processed and saved.')
end

end

