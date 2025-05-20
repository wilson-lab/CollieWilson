% ephysSettings.m
%
% Function that returns all the CONSTANTS for an electrophysiology
%  experiment
%
% OUTPUT:
%   dataDir - name of base data directory
%   exptFnDir - name of folder containing all experiment functions
%   settings - struct of settings
%
% CREATED:
% 01/30/2021 - HY
% 02/02/2021 - MC adapted
% 05/12/2021 - MC added panels, removed leg cam
% 08/19/2021 - MC changed dataDir
% 09/14/2022 - MC g4 panels now through DAC instead of log
% 10/06/2022 - MC added opto stim output
%

function [dataDir, exptFnDir, settings] = ephysSettings()    
    %% Data folder
    % Determine which computer this code is running on
    comptype = computer; % get the string describing the computer type
    PC_STRING = 'PCWIN64';  % string for PC on 2P rig
    
    % path to folder for saving all experiment data
    dataDir = 'D:\';
    % path to folder containing all experiment defining functions
    mfile_path = mfilename('fullpath');
    [root_path, ~, ~] = fileparts(mfile_path);
    exptFnDir = fullfile(root_path, 'experiment types');
    % add experiment function path
    addpath(exptFnDir);


    %% DAQ settings

    % Devices
    settings.devVendor = 'ni';
    settings.bob.devID = 'Dev1';
    settings.temp.devID = 'Dev1';

    % Sampling Rate
    settings.bob.sampRate  = 20e3; % 20 kHz, same as Yvette
     % Sampling rate for measuring access resistance
     % NI 6361 DAQ max for 6 ephys channels 166 kHz (1MS/s)
     % settings.bob.sampRateRacc = 10e4; % 100 kHz 
    

    %% Break out box, channel assignments
    % which analog input channels are used
    settings.bob.aInChUsed  = [0:12];
    
    % which digital input channels are used (matrix each row is channel, column
    %  1 is port number, column 2 is line number)
    settings.bob.dInChUsed = [0 5];
    
    % which analog output channels are used
    settings.bob.aOutChUsed = [0 1];
    % which digital output channels are used (notation as above)
    settings.bob.dOutChUsed = [0 4];
    
    %% Break out box, channel decode
    % to decode which column in raw data output from data acquisition
    %  corresponds to what information; ordered by order channels will be added
    %  to DAQ session; add analog before digital
    settings.bob.aInChAssign = {'ampI', 'amp10Vm', ...
        'ampScaledOut', 'ampMode', 'ampGain', 'ampFreq', ...
        'g3panelXPosition', 'g3panelYPosition', ...
        'ficTracHeading', 'ficTracIntX', 'ficTracIntY',...
        'g4panelXPosition','pythonJumpTrig'};
    settings.bob.dInChAssign = {'legCamFrames'};
    settings.bob.inChAssign = [settings.bob.aInChAssign ...
        settings.bob.dInChAssign];
    % output channel assignments (notation like input, above)
    settings.bob.aOutChAssign = {'ampExtCmdIn','optoExtCmd'};
    settings.bob.dOutChAssign = {'legCamFrameStartTrig'};
    settings.bob.outChAssign = [settings.bob.aOutChAssign ...
        settings.bob.dOutChAssign];

    %% Channel settings
    MV_PER_V = 1000; % millivolts per volts

    % Analog input channel settings
    % analog input type 'SingleEnded' as opposed to 'Differential' (no
    %  comparisons across 2 BNCs of break out box; switch on SE)
    settings.bob.aiMeasType = 'Voltage';
    settings.bob.aiInType = 'SingleEnded'; 
    % voltage range - for channels, in order in aInChAssign
    settings.bob.aiRange = repmat([-10 10],length(settings.bob.aInChAssign),1);

    % Digital input channel settings
    % digital input type - 'InputOnly', not 'Bidirectional'
    settings.bob.diType = 'InputOnly'; 
    
    % Analog output channel settings
    % analog output can be 'Voltage' or 'Current'
    settings.bob.aoMeasType = 'Voltage';
    % analog output as 'SingleEnded'
    settings.bob.aoOutType = 'SingleEnded';
    % voltage range - for channels, in order of aOutChAssign
    settings.bob.aoRange = repmat([-5 5],length(settings.bob.aOutChAssign),1);
    
    % Digital output channel settings
    % digital output type - 'OutputOnly', not 'Bidirectional'
    settings.bob.doType = 'OutputOnly';
    
    % Thermocouple (USB-TC01) only channel
    settings.temp.aiChUsed = 'ai0';
    settings.temp.aiMeasType = 'Thermocouple';
    settings.temp.tcType = 'J'; % thermocouple type
    
    
    % Static conversion factors on ephys channels
    settings.amp.beta = 1; % beta value for Axopatch 200B, whole cell
    
    % Current (beta mV/pA)
    settings.I.sigCondGain = 1; % signal conditioner not currently in use
    settings.I.sigCondFreq = nan; % LP filter, in kHz, sig cond not in use
    settings.I.ampGain = 1; % has option in back for this to be 100
    % conversion from V reading from DAQ to pA measured
    settings.I.softGain = MV_PER_V / (settings.I.sigCondGain * ...
        settings.amp.beta * settings.I.ampGain);
    
    % Voltage (10 Vm)
    settings.Vm.sigCondGain = 1; % signal conditioner not currently in use
    settings.Vm.sigCondFreq = nan; % LP filter, in kHz, sig cond not in use
    settings.Vm.ampGain = 10; % amp gain
    % conversion from V reading from DAQ to mV measured
    settings.Vm.softGain = MV_PER_V / (settings.Vm.sigCondGain * ...
        settings.Vm.ampGain);
    
    % path to folder containing all current functions
    settings.VOut.iInjFnDir = [exptFnDir '\Current Injection'];
    
    % Voltage Output (V Clamp 20 mV/V, I Clamp 2/beta nA/V)
    settings.VOut.vDivGain = 1; % currently, no voltage divider
    settings.VOut.ampVCmdGain = 20; % 20 mV/V
    settings.VOut.VConvFactor = 1 / (settings.VOut.ampVCmdGain * ...
        settings.VOut.vDivGain);
    settings.VOut.zeroV = -0.012880558953762; % mV, measured 7/23/20 - HHY
    settings.VOut.ampICmdGain = 2000 / settings.amp.beta; % 2000/beta pA/V
    settings.VOut.IConvFactor = 1 / (settings.VOut.ampICmdGain * ...
        settings.VOut.vDivGain);
    settings.VOut.zeroI = -0.965853722646187; % pA, measured 7/23/20 - HHY
    
    
    % Some static parameters for testing pipette/seal/access resistances
    % duration of recording for pipette/seal resistances measurements in
    %  pre-expt routine
    settings.sealTestDur = 2; % in sec
    
end
