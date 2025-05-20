% recordFictracOnly.m
%
% Trial Type Function for recording only FicTrac channels during an experiment. 
% This function focuses solely on capturing the motion data from the FicTrac system 
% without additional output channels or stimuli.
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
% - rawOutput  : Empty matrix for this trial type, included to maintain 
%                consistency with the trial type function format.
%
% ORIGINAL: 01/30/2021 by HY
% UPDATED: 04/01/2021 by MC
%
function [rawData, inputParams, rawOutput] = recordFictracOnly(settings, ...
    duration)

% EXPERIMENT-SPECIFIC PARAMETERS
inputParams.exptCond = 'fictracOnly'; % name of trial type

% which input and output data streams used in this experiment
inputParams.aInCh = {'ficTracHeading', 'ficTracIntX', 'ficTracIntY'};
inputParams.aOutCh = {};
inputParams.dInCh = {'ficTracCamFrames'};
inputParams.dOutCh = {};

% output matrix - empty for this trial type (there are no output
%  channels initialized for DAQ; no current injection, triggers, etc)
% placeholder for different trial types
rawOutput = [];

% save trial duration here into inputParams
inputParams.trialDuration = duration;

% initialize DAQ, including channels
[userDAQ, ~, ~, ~, ~] = initUserDAQ(settings, ...
    inputParams.aInCh, inputParams.aOutCh, inputParams.dInCh, ...
    inputParams.dOutCh);

% set duration of acquisition
userDAQ.DurationInSeconds = duration;

% get time stamp of approximate experiment start
inputParams.startTimeStamp = datestr(now, 'HH:MM:SS');

disp(['[' datestr(now,'HH:MM') '] Beginning behavior only'])
% acquire data (in foreground)
rawData = userDAQ.startForeground();

disp('Finished behavior trial.');
end