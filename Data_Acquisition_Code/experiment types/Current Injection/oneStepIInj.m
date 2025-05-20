% oneStepIInj.m
%
% Current Injection Function. Injects square current steps of user 
%  specified amplitudes of user specified durations. Durations can be 
%  different among steps. Has option to randomize amplitudes or not. User 
%  specified step duration and time between steps. User specified space 
%  amplitude.
% Returns vector of to feed to output of DAQ (i.e. appropriately scaled to
%  be read directly by amplifier)
% Starts with space
% User specifies actual amplitudes in pA
%
% INPUTS:
%   settings - struct returned by ephysSettings()
%   durScans - duration of trial in scans
%
% OUTPUTS:
%   iInjOut - col vector of current injection output, of length durScans
%   iInjParams -struct with all user specified parameter values
%
% Created: 7/23/20 - HHY
%
% Updated: 4/04/21 - MC
%

function [iInjOut, iInjParams] = oneStepIInj(settings, duration)

    % option 1 - set parameters
    iInjParams.stepAmp = 50;  % step amplitude (pA)
    iInjParams.stepDur = 5;   % step duration (s)
    iInjParams.spaceAmp = 0;   % amplitude b/n steps (pA)
    iInjParams.spaceDur = 20;  % duration b/n steps (s)
    
%     % option 2 - user enters parameters
%     % prompt user for input parameters, as dialog box
%     inputParams = {'Step Amplitudes (pA):', 'Step Durations (s):', ...
%         'Space Amplitude (pA):', 'Space Duration (s):', ...
%         'Randomize steps? y/n'};
%     dlgTitle = 'Enter parameter values';
%     dlgDims = [1 35]; % dimensions of dialog box input fields
%     
%     % dialog box
%     dlgAns = inputdlg(inputParams, dlgTitle, dlgDims);
%     
%     % convert user input into actual variables
%     stepAmps = str2double(dlgAns{1});
%     stepDurs = str2double(dlgAns{2});
%     spaceAmp = str2double(dlgAns{3});
%     spaceDur = str2double(dlgAns{4});
%     randomize = dlgAns{5};
    
    % convert user input into correct units for output (amplitude in volts,
    %  duration in scans);
    % not included: compensate for non-zero output from DAQ when zero commanded
    stepAmpV = (iInjParams.stepAmp - settings.VOut.zeroI) * settings.VOut.IConvFactor;
    stepDurSR = round(iInjParams.stepDur * settings.bob.sampRate);
    
    spaceAmpV = (iInjParams.spaceAmp - settings.VOut.zeroI)) * settings.VOut.IConvFactor;
    spaceDurSR = round(iInjParams.spaceDur * settings.bob.sampRate);
    
    % each stimulus as a column, starts with space, ends with step
    spaceMatrix = ones(spaceDurSR, 1) * spaceAmpV;
    stepMatrix = ones(stepDurSR, 1) .* stepAmpV;
    stimMatrix = [spaceMatrix; stepMatrix];
    
    % generate stimulus output
    % reshape stim matrix to single column vector
    oneRepStim = reshape(stimMatrix, numel(stimMatrix), 1);
    
    % add spacer at the end and output
    iInjOut = [oneRepStim ; spaceMatrix(:,1)];

end