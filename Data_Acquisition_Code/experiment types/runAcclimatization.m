% runAcclimatization.m
%
% Trial Type Function for displaying a blank G4 pattern/function and recording 
% data from G4 panels and FicTrac channels. This function delivers a specified 
% duration of light pulse at the start of each trial.
%
% INPUTS:
% - settings  : Struct containing electrophysiological setup settings, 
%               typically obtained from the ephysSettings() function.
% - duration  : Duration of the acclimatization period (in seconds).
% - loop      : Specifies the trial type; options are:
%               0 - open loop, 1 - closed loop with specific settings,
%               2 - alternative closed loop settings.
%
% OUTPUTS:
% - rawData   : Matrix of raw data measured by the DAQ, where each column 
%               corresponds to data from a different channel.
% - inputParams : Struct containing parameters specific to this experiment 
%                 type, including experimental conditions and channel settings.
% - rawOutput  : Empty matrix for this trial type as there is no output data.
%
% CREATED: 03/22/2022 by MC
% UPDATED: 09/14/2022 by MC, G4 panels now controlled through DAC instead of log.
%
function [rawData, inputParams, rawOutput] = runAcclimatization(settings,duration,loop)

%% INITIALIZE DAQ
inputParams.exptCond = 'G4PanelsFictrac'; % name of trial type

% which input and output data streams used in this experiment
inputParams.aInCh = {'ficTracHeading', 'ficTracIntX', 'ficTracIntY', 'g4panelXPosition'};
inputParams.aOutCh = {};
inputParams.dInCh = {};
inputParams.dOutCh = {};

% output matrix - empty for this trial type (there are no output)
rawOutput = [];

% initialize DAQ, including channels
[userDAQ, ~, ~, ~, ~] = initUserDAQ(settings, inputParams.aInCh, inputParams.aOutCh,...
    inputParams.dInCh, inputParams.dOutCh);


%% SET PANEL PARAMETERS

disp('Initializing G4 panels...');
% panel settings
switch loop
    case 0
        pattN = 1; %blank background
        funcN = 1; %hold
        mode = 1; %pos change func

        % pull settings and connect
        userSettings %load settings
        Panel_com('change_root_directory', exp_path);
        load([exp_path '\currentExp'])
        expt_t = duration;

        %store pattern parameters
        load(['patt_lookup_' sprintf('%04d', pattN)])
        inputParams.pattern_name = patlookup.fullname;
        inputParams.objectSize = patlookup.size;
        inputParams.objectGS = patlookup.objgs;
        inputParams.backgroundGS = patlookup.bckgs;
        %store function parameters
        load(['func_lookup_' sprintf('%04d', funcN)])
        inputParams.function_name = funlookup.name;
        inputParams.sweepRange = funlookup.sweepRange;
        inputParams.sweepRate = funlookup.sweepRange;
        inputParams.sweepDur = (funlookup.sweepRange/funlookup.sweepRate)*2;

        % load panel parameters
        Panel_com('set_control_mode', mode);
        Panel_com('set_pattern_id', pattN);
        Panel_com('set_pattern_func_id', funcN);

    case 1
        % panel settings
        pattN = 2; %6px dark bar only

        %mode = 4; %closed loop - frame rate
        mode = 7; %closed loop - frame index

        panel_gain = 80;

        % pull settings and connect
        userSettings %load settings
        Panel_com('change_root_directory', exp_path);
        load([exp_path '\currentExp'])
        expt_t = duration;

        %store pattern parameters
        load(['patt_lookup_' sprintf('%04d', pattN)])
        inputParams.pattern_name = patlookup.name;
        inputParams.objectSize = patlookup.size;
        inputParams.objectGS = patlookup.objgs;
        inputParams.backgroundGS = patlookup.bckgs;

        % load panel parameters
        Panel_com('set_control_mode', mode);
        Panel_com('set_gain_bias', [panel_gain 0]);
        Panel_com('set_pattern_id', pattN);
    case 2
        % panel settings
        pattN = 2; %6px dark bar only
        %pattN = 6; %6px bright bar only

        %mode = 4; %closed loop - frame rate
        mode = 7; %closed loop - frame index

        panel_gain = 80;

        % pull settings and connect
        userSettings %load settings
        Panel_com('change_root_directory', exp_path);
        load([exp_path '\currentExp'])
        expt_t = duration;

        %store pattern parameters
        load(['patt_lookup_' sprintf('%04d', pattN)])
        inputParams.pattern_name = patlookup.name;
        inputParams.objectSize = patlookup.size;
        inputParams.objectGS = patlookup.objgs;
        inputParams.backgroundGS = patlookup.bckgs;

        % load panel parameters
        Panel_com('set_control_mode', mode);
        Panel_com('set_gain_bias', [panel_gain 0]);
        Panel_com('set_pattern_id', pattN);
end





%% ACQUIRE
disp(['[' datestr(now,'HH:MM') '] Beginning acclimation period: ' patlookup.name ' w/ ' 'closed-loop'])
inputParams.trialDuration = expt_t;
userDAQ.DurationInSeconds = expt_t;

% start experiment
pause(0.1) %slight delay
inputParams.startTimeStamp = datestr(now, 'HH:MM:SS');
Panel_com('start_display', expt_t); %sec
rawData = userDAQ.startForeground(); % acquire data (in foreground)
% stop experiment
pause(0.1)

disp('Finished acclimatization block.');
end