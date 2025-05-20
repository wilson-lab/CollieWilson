% measureSpikeRate.m
%
% Function to record a brief 10-second trace and report the spike rate. 
% This function is helpful to run at the start of experiments to ensure 
% the proper amount of current is being injected.
%
% INPUTS:
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
% ORIGINAL: 12/08/2021 by MC
%
function measureSpikeRate(duration)

    % load settings
    [~, ~ , settings] = ephysSettings();
    
    % which input and output data streams used in this experiment
    inputParams.exptCond = 'ephys'; % name of trial type
    inputParams.aInCh = {'ampI', 'amp10Vm', 'ampScaledOut','ampMode', 'ampGain', 'ampFreq'};
    inputParams.aOutCh = {};
    inputParams.dInCh = {};
    inputParams.dOutCh = {};
    rawOutput = [];
    inputParams.trialDuration = duration; 

    % initialize DAQ, including channels
    [userDAQ, ~, ~, ~, ~] = initUserDAQ(settings, inputParams.aInCh, inputParams.aOutCh, inputParams.dInCh, inputParams.dOutCh);
    
    % set duration of acquisition
    userDAQ.DurationInSeconds = duration;
    % get time stamp of approximate experiment start
    inputParams.startTimeStamp = datestr(now, 'HH:MM:SS');
    
    
    disp('Acquiring spike rate recording.');
    % acquire data (in foreground)
    rawData = userDAQ.startForeground();
    disp('Spike rate recording acquired.');
    
    % process
    [daqData, daqOutput, daqTime] = preprocessUserDaq(inputParams, rawData, rawOutput, settings);
    [exptData, ~] = processExptData(daqData, daqOutput, daqTime, inputParams, settings);
    
    % find spike rate
    [~,~,spikeRate] = findSpikeRate(settings,exptData,7,0);
    vc_avgspikerate = nanmean(spikeRate);
    disp(['Spike rate = ' num2str(round(vc_avgspikerate)) '/s']);


end