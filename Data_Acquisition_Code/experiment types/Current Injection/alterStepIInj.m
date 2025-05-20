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
%          7/04/23 - MC alternating depol/hyperpol
%

function [iInjOut, iInjParams] = alterStepIInj(settings, duration)

    % set parameters
    iInjParams.stepUpAmp = 100;  % depolarizing step amplitude (pA)\
    iInjParams.stepDownAmp = 100;  % hyperpolarizing step amplitude (pA)
    iInjParams.spaceAmp = 0;   % amplitude b/n steps (pA)
    iInjParams.stepDur = 1;   % step duration (s)
    iInjParams.spaceDur = 5;  % duration b/n steps (s)
    
    % convert user input into correct units for output (amplitude in volts, duration in scans);
    % included: compensate for non-zero output from DAQ when zero commanded
    stepUpAmpV = (iInjParams.stepUpAmp - settings.VOut.zeroI) * settings.VOut.IConvFactor;
    stepDownAmpV = (iInjParams.stepDownAmp - settings.VOut.zeroI) * settings.VOut.IConvFactor;
    stepDurSR = round(iInjParams.stepDur * settings.bob.sampRate);
    
    spaceAmpV = (iInjParams.spaceAmp- settings.VOut.zeroI) * settings.VOut.IConvFactor;
    spaceDurSR = round(iInjParams.spaceDur * settings.bob.sampRate);
    
    % number of scans to be run, round down based on trial duration
    iInjParams.numScans = floor(duration/(iInjParams.stepDur+iInjParams.spaceDur));
    numUpScans = ceil(iInjParams.numScans/2); %n depol pulses
    numDownScans = floor(iInjParams.numScans/2); %n hyperpol pulses
    
    % randomize step up/down delivery
    stimSelect = [ones(numUpScans,1) ; zeros(numDownScans,1)]; %initialize
    stimSelect = stimSelect(randperm(length(stimSelect))); %randomize
    
    % set space matrix
    spaceMatrix = ones(spaceDurSR, 1) * spaceAmpV;
    % set step matrix, randomizing delivery
    stimMatrix = [spaceMatrix ;spaceMatrix]; %initialize
    for s = 1:iInjParams.numScans
        switch stimSelect(s)
            case 1 %depol pulse
                stepMatrix = ones(stepDurSR, 1) .* stepUpAmpV;
            case 0 %hyperpol pulse
                stepMatrix = ones(stepDurSR, 1) .* stepDownAmpV;
        end
        % store stim matrix
        stimMatrix = [stimMatrix; stepMatrix; spaceMatrix];
    end
    
    % add spacer at the end and output
    iInjOut = stimMatrix;

end