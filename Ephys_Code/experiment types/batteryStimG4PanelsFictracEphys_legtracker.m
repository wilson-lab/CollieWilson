% recordFictracEphysIInj.m
%
% Trial Type Function for recording electrophysiological (ephys) channels 
% while performing current injections. This function captures data from 
% both the FicTrac system and ephys channels, and triggers leg camera frames 
% during the trial.
%
% INPUTS:
% - settings  : Struct containing electrophysiological setup settings, 
%               typically obtained from the ephysSettings() function.
% - duration  : Duration of the trial (in seconds).
%
% OUTPUTS:
% - rawData   : Matrix of raw data measured by the DAQ, where each column 
%               corresponds to data from a different channel.
% - inputParams : Struct containing parameters specific to this experiment 
%                 type, including experimental conditions and channel settings.
% - rawOutput  : Matrix of raw output sent by the DAQ, where each column 
%                corresponds to a different channel.
%
% ADAPTED: 10/30/2023 by MC
%
function [rawData, inputParams, rawOutput] = batteryStimG4PanelsFictracEphys_legtracker(settings,duration,pattN,funcN,stim)

%% INITIALIZE DAQ
% Initialize global variable for raw data collection
global daqData
% Initialize global variable for knowing when to prompt for
%  starting acquisition on leg camera (i.e. only first time leg video
%  is acquired for cell)

% binary variable for whether DAQ should be stopped; used in nested
%  function queueLegTrig()
acqStopBin = 0; % starts as 0 (no)

% start index (in scans) of output
% will count up every time DataRequired event is called
whichOutScan = 1; % start at 1

% EXPERIMENT-SPECIFIC PARAMETERS
inputParams.exptCond = 'legvidG4PanelsFictracEphys'; % name of trial type
% leg tracking camera frame rate - make sure it's a whole number of
%  DAQ scans
legCamFrameRate = 50; % in Hz
legCamFrameRateScans = round(settings.bob.sampRate / legCamFrameRate);
inputParams.legCamFrameRate = settings.bob.sampRate / ...
    legCamFrameRateScans;

% FicTrac camera frame rate - make sure it's a whole number of DAQ
%  scans
ftCamFrameRate = 150; % in Hz
ftCamFrameRateScans = round(settings.bob.sampRate / ftCamFrameRate);
inputParams.ftCamFrameRate = settings.bob.sampRate / ...
    ftCamFrameRateScans;

% which input and output data streams used in this experiment
inputParams.aInCh = {'ampI', 'amp10Vm', 'ampScaledOut', ...
    'ampMode', 'ampGain', 'ampFreq', ...
    'ficTracHeading', 'ficTracIntX', 'ficTracIntY', 'g4panelXPosition'};
inputParams.aOutCh = {'optoExtCmd'};
inputParams.dInCh = {'legCamFrames'};
inputParams.dOutCh = {'legCamFrameStartTrig'};

% index for output channels: analog before digital
optoTrigChInd = 1;
legTrigChInd = 2;
totNumOutCh = 2; % total number of output channels

% initialize DAQ, including channels
[userDAQ, ~, ~, ~, ~] = initUserDAQ(settings, ...
    inputParams.aInCh, inputParams.aOutCh, inputParams.dInCh, ...
    inputParams.dOutCh);
% make DAQ acquisition continuous (runs until stopped)
userDAQ.IsContinuous = true;

% pre-allocate global variables
% maximum number of scans in acquisition, with 5 second buffer
maxNumScans = round((duration + 5) * userDAQ.Rate);
% number of channels data is being acquired on
numInCh = length(inputParams.aInCh) + length(inputParams.dInCh);
% pre-allocate global variable for data acquired
daqData = zeros(maxNumScans, numInCh);
% number of channels data is being output on
numOutCh = length(inputParams.aOutCh) + length(inputParams.dOutCh);
% pre-allocate variable for data output
daqOutput = zeros(maxNumScans, numOutCh);


%% SET OPTO SETTINGS (IF SELECTED)
% generate output matrix for opto stimulation
if stim
    optoOut = ones(duration*settings.bob.sampRate,1)*5; %generate stim array, 5V output
else
    optoOut = zeros(duration*settings.bob.sampRate,1); %generate empty stim array, 0V output
end


%% SET EXPERIMENT TIMING

% delay start of acquisition on other hardware by 0.5 sec to ensure
%  user DAQ starts first and captures everything
startDelay = 0.5; % in seconds
startDelayScans = round(startDelay * userDAQ.Rate); % in scans
% save actual into inputParams
inputParams.startDelay = startDelayScans / userDAQ.Rate;

% timing for triggers to leg camera
% amount of data in seconds to queue each time DataRequired event is
%  fired
camTrigQueuedLen = 1;
% queue in scans
camTrigQueuedScans = round(camTrigQueuedLen * userDAQ.Rate);
% save actual queued length into inputParams
inputParams.legTrigQueuedLen = camTrigQueuedScans / userDAQ.Rate;

% DataRequired event fires whenever queued data falls below threshold -
%  use default here of 0.5 sec
userDAQ.NotifyWhenScansQueuedBelow = 0.5 * userDAQ.Rate;
inputParams.legTrigQueuedBelow = userDAQ.NotifyWhenScansQueuedBelow ...
    / userDAQ.Rate;

% how often to fire DataAvailable event for background acquisition (0.5
%  sec)
dataAvailExceeds = 0.5; % in seconds
dataAvailExceedsScans = round(dataAvailExceeds * userDAQ.Rate); % in scans
% save actual into inputParams
inputParams.dataAvailExceeds = dataAvailExceedsScans / userDAQ.Rate;
% set value on DAQ
userDAQ.NotifyWhenDataAvailableExceeds = dataAvailExceedsScans;

% delay end of acquisition on user DAQ by leg trigger queue length * 2
%  plus threshold of DataRequired event to ensure user DAQ captures
%  everything
inputParams.endDelay = inputParams.legTrigQueuedLen * 2 + ...
    inputParams.legTrigQueuedBelow;

% save initial experiment duration here into inputParams
inputParams.initialExptDuration = duration;


% QUEUE INITIAL OUTPUT - leg camera frame triggers
% number of scans to queue initially - delay + initial bout of data
numScans = startDelayScans + camTrigQueuedScans;

outputInit = zeros(numScans, totNumOutCh);
legCamTrigInit = zeros(numScans, 1);
camStartInd = startDelayScans + 1;

% current injection output: delay and then start protocol
outputInit((startDelayScans+1):numScans, optoTrigChInd) = ...
    optoOut(1:(numScans-startDelayScans));

% generate trigger pattern for leg camera - starts with 0 for start
%  delay time, then a 1 at leg camera frame rate
legCamTrigInit(camStartInd:legCamFrameRateScans:end) = 1;


% add camera trigger patterns to outputInit
outputInit(:, legTrigChInd) = legCamTrigInit;

% queue output on DAQ
userDAQ.queueOutputData(outputInit);

% save queued output into daqOutput
lengthOut = size(outputInit,1);
daqOutput(whichOutScan:(whichOutScan + lengthOut - 1),:) = outputInit;
% update whichOutScan for next iteration
whichOutScan = whichOutScan + lengthOut;

% generate leg camera trigger pattern (to be queued after intial set)
legCamTrig = zeros(camTrigQueuedScans, 1);
% trigger pattern for leg camera
legCamTrig(1:legCamFrameRateScans:end) = 1;

% output matrix preallocate
outputMatrix = zeros(camTrigQueuedScans, totNumOutCh);

% output matrix of all zeros, for end
outputMatrixEnd = zeros(camTrigQueuedScans, totNumOutCh);

% nested function for queuing more output (camera trigger and current
%  injection; called by event listener for DataRequired
    function queueOut(src, event)
        if ~acqStopBin
            % add leg camera triggers to output matrix
            outputMatrix(:, legTrigChInd) = legCamTrig;
            % grab next set of output from current injection
            optoStartInd = whichOutScan - startDelayScans;
            optoEndInd = whichOutScan + camTrigQueuedScans - ...
                startDelayScans - 1;
            outputMatrix(:, optoTrigChInd) = optoOut(optoStartInd:optoEndInd);

            queueOutputData(src, outputMatrix);

            % save queued output into daqOutput
            lenOut = size(outputMatrix, 1);
            daqOutput(whichOutScan:(whichOutScan + lenOut - 1),:) = ...
                outputMatrix;
        else % when DAQ acquisition is stopped
            queueOutputData(src, outputMatrixEnd);

            % save queued output into daqOutput
            lenOut = size(outputMatrixEnd, 1);
            daqOutput(whichOutScan:(whichOutScan + lenOut - 1),:) = ...
                outputMatrixEnd;
        end
        % update whichOutScan for next iteration
        whichOutScan = whichOutScan + lenOut;
    end

% create listeners for DataAvailable and DataRequired events
dataAvailLh = addlistener(userDAQ, 'DataAvailable', @collectData);
dataReqLh = addlistener(userDAQ, 'DataRequired', @queueOut);

% first time leg video is acquired for cell, set up folder for saving leg video
if ~exist('rawLegVid','dir')
    % make folder for raw images
    mkdir('rawLegVid');

end
% raw leg video full path
legVidFileName = sprintf('%s%srawLegVid%slegVid', pwd, filesep, ...
    filesep);
% copy path to clipboard
clipboard('copy', legVidFileName);
% prompt user to copy path into spinview
fprintf(['Leg Video Acquisition. \n'...
    'Press RECORD button and paste directory from system clipboard '...
    'into the *Filename* section. \n Set *Image Format* to Tiff ' ...
    'and *Compression Method* to Rle. Then press Start Recording.' ...
    '\n Make sure camera is acquiring (green play button). \n']);

% prompt user for current number of leg vid frames grabbed
prompt = 'Enter current number of leg video frames grabbed: ';
inputParams.startLegVidNum = str2double(input(prompt, 's'));

% get time stamp of approximate experiment start
inputParams.startTimeStamp = datestr(now, 'HH:MM:SS');
fprintf('Start time: %s \n', inputParams.startTimeStamp);
disp('Starting legvidFictracvidEphysIInj acquisition');


%% SET PANEL PARAMETERS

disp('Initializing G4 panels...');
% panel settings
mode = 1; %pos change func

% pull settings and connect
Panel_com('change_root_directory', 'C:\Users\wilson\Dropbox (HMS)\MATLAB\G4_Display_Tools\PControl_Matlab\Experiment');

% load panel parameters
Panel_com('set_control_mode', mode);
Panel_com('set_pattern_id', pattN);
Panel_com('set_pattern_func_id', funcN);


%% ACQUIRE IN BACKGROUND
disp(['[' datestr(now,'HH:MM') '] Beginning: alternating background w/ leg tracking'])
Panel_com('start_display', duration); %sec
userDAQ.startBackground();

% total number of scans (less than maxNumScans, which has 5 sec buffer
totNumScans = round(duration * userDAQ.Rate);
% loop that allows acquisition to stop smoothly
while 1
    % if any key on keyboard is pressed, initialize stopping of
    %  acquisition; or if specified duration of acquisition is reached
    if (userDAQ.ScansAcquired > totNumScans)
        % binary, will stop triggering leg camera through
        %  queueLegTrig function
        acqStopBin = 1;
        disp('Stopping acquisition, please wait');
        % pause to give DAQ enough time to stop triggering leg camera
        %  and acquire all data
        pause(inputParams.endDelay);
        break; % stop loop
    end
    % this loop doesn't need to go that quickly to register keyboard
    %  presses
    pause(0.2);
end

% once looping stops, stop acquisition
userDAQ.stop();
disp('Acquisition stopped');
fprintf('End time: %s \n', datestr(now, 'HH:MM:SS'));

% save actual experiment duration into inputParams
inputParams.actualExptDuration = userDAQ.ScansAcquired / userDAQ.Rate;

% only keep data and output up until point when acquisition stopped
daqData = daqData(1:userDAQ.ScansAcquired, :);
daqOutput = daqOutput(1:userDAQ.ScansAcquired, :);

% display number of leg video frames triggered
numLegVidTrigs = sum(daqOutput(:,legTrigChInd));
fprintf('%d leg video frames triggered.\n', numLegVidTrigs);

% prompt user for current number of leg video frames grabbed
prompt = 'Enter current number of leg video frames grabbed: ';
inputParams.endLegVidNum = str2double(input(prompt, 's'));

% display total number of leg video frames acquired
numLegVidAcq = inputParams.endLegVidNum - inputParams.startLegVidNum;
fprintf('%d leg video frames grabbed. \n', numLegVidAcq);

% save global variables into variables returned by this function
rawData = daqData;
rawOutput = daqOutput;

% delete global variable
clear global daqData

% delete listeners
delete(dataAvailLh);
delete(dataReqLh);
end