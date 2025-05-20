% multiStepIInj.m
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
% Updated: 4/04/21 - MC
%

function [iInjOut, iInjParams] = multiStepIInj(settings, duration)

    % option 1 - set parameters
    iInjParams.stepAmp = 100;  % step amplitude (pA)
    iInjParams.stepDur = 0.5;   % step duration (s)
    iInjParams.spaceAmp = 0;   % amplitude b/n steps (pA)
    iInjParams.spaceDur = 9.5;  % duration b/n steps (s)
    
    % convert user input into correct units for output (amplitude in volts,
    %  duration in scans);
    % included: compensate for non-zero output from DAQ when zero commanded
    stepAmpV = (iInjParams.stepAmp - settings.VOut.zeroI) * settings.VOut.IConvFactor;
    stepDurSR = round(iInjParams.stepDur * settings.bob.sampRate);
    
    spaceAmpV = (iInjParams.spaceAmp- settings.VOut.zeroI) * settings.VOut.IConvFactor;
    spaceDurSR = round(iInjParams.spaceDur * settings.bob.sampRate);
    
    % number of scans to be run, round down based on trial duration
    iInjParams.numScans = floor(duration/(iInjParams.stepDur+iInjParams.spaceDur));
    
    % each stimulus as a column, starts with space, ends with step
    spaceMatrix = ones(spaceDurSR, iInjParams.numScans) * spaceAmpV;
    stepMatrix = ones(stepDurSR, iInjParams.numScans) .* stepAmpV;
    stimMatrix = [spaceMatrix; stepMatrix];
    
    % generate stimulus output
    % reshape stim matrix to single column vector
    oneRepStim = reshape(stimMatrix, numel(stimMatrix), 1);
    
    % add spacer at the end and output
    iInjOut = [oneRepStim ; spaceMatrix(:,1)];

end