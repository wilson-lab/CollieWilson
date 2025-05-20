% recordEphysIInj.m
%
% Trial Type Function for recording electrophysiological (ephys) channels 
% while performing current injections. This function captures ephys data 
% during the trial and executes the specified current injection protocol.
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
% CREATED: 07/17/2020 by HHY
% UPDATED: 04/02/2021 by MC
%
function [rawData, inputParams, rawOutput] = recordEphysIInj(settings, duration)

    % EXPERIMENT-SPECIFIC PARAMETERS
    inputParams.exptCond = 'ephysIInj'; % name of trial type
    
    % which input and output data streams used in this experiment
    inputParams.aInCh = {'ampI', 'amp10Vm', 'ampScaledOut', ...
        'ampMode', 'ampGain', 'ampFreq'};
    inputParams.aOutCh = {'ampExtCmdIn'};
    inputParams.dInCh = {};
    inputParams.dOutCh = {};
    
    % save trial duration here into inputParams
    inputParams.trialDuration = duration; 

    % initialize DAQ, including channels
    [userDAQ, ~, ~, ~, ~] = initUserDAQ(settings, ...
        inputParams.aInCh, inputParams.aOutCh, inputParams.dInCh, ...
        inputParams.dOutCh);
    
    % path to current injection protocol functions
    iPath = settings.VOut.iInjFnDir;
    
    % prompt user to enter function call to current injection function
        % prompt user to select an experiment
    iInjSelected = 0;
    disp('Select a current injection protocol');
    while ~iInjSelected
        iInjTypeFileName = uigetfile('*.m', ...
            'Select a current injection protocol', iPath);
        % if user cancels or selects valid file
        if (iInjTypeFileName == 0)
            disp('Selection cancelled');
            iInjSelected = 1; % end loop
        elseif (contains(iInjTypeFileName, '.m'))
            disp(['Protocol: ' iInjTypeFileName]);
            iInjSelected = 1; % end loop
        else
            disp('Select a current injection .m file or cancel');
            iInjSelected = 0;
        end
    end

    % if user cancels at this point 
    if (iInjTypeFileName == 0)
        % throw error message; ends run of this function
        error('No current injection protocol was run. Ending ephysIInj()');
    end
    
    % convert selected experiment file into function handle
    % get name without .m
    iInjTypeName = extractBefore(iInjTypeFileName, '.');
    iInjFn = str2func(iInjTypeName);

    % run current injection function to get output vector
    [iInjOut, iInjParams] = iInjFn(settings, duration); 
    
    % save info into returned variables
    rawOutput = iInjOut; % output commanded into rawOutput
    % record current injection function name
    inputParams.iInjProtocol = iInjTypeName; 
    inputParams.iInjParams = iInjParams; % current injection parameters
    
    % queue current injection output
    userDAQ.queueOutputData(iInjOut);
    
    % get time stamp of approximate experiment start
    inputParams.startTimeStamp = datestr(now, 'HH:MM:SS');
    
    disp(['[' datestr(now,'HH:MM') '] Beginning: ephys w/ ' iInjTypeName])
    % acquire data (in foreground)
    rawData = userDAQ.startForeground();
    
    % to stop it from presenting non-zero values if current injection
    %  protocol ends on non-zero value
    userDAQ.outputSingleScan(0);

    disp('Finished i-inj trial.');
end