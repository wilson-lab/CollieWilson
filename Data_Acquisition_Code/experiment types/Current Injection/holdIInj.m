% multiStepIInj.m
%
% Current Injection Function. Injects constant current for specified
% duration
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
% Created: 07/23/20 - HHY
% Updated: 08/24/23 - MC generated from alternating step protocol
%

function [iInjOut, iInjParams] = holdIInj(settings, holdAmp, duration)

    % set parameters
    iInjParams.stepDownAmp = holdAmp;  % hyperpolarizing step amplitude (pA)
    iInjParams.spaceAmp = 0;   % amplitude b/n steps (pA)
    iInjParams.stepDur = duration;   % step duration (s)
    iInjParams.spaceDur = 5;  % duration b/n steps (s)
    
    % convert user input into correct units for output (amplitude in volts, duration in scans);
    % included: compensate for non-zero output from DAQ when zero commanded
    stepDownAmpV = (iInjParams.stepDownAmp - settings.VOut.zeroI) * settings.VOut.IConvFactor;
    stepDurSR = round(iInjParams.stepDur * settings.bob.sampRate);
    
    spaceAmpV = (iInjParams.spaceAmp- settings.VOut.zeroI) * settings.VOut.IConvFactor;
    spaceDurSR = round(iInjParams.spaceDur * settings.bob.sampRate);
    
    % set space matrix
    spaceMatrix = ones(spaceDurSR, 1) * spaceAmpV;
    % set hold matrix
    holdMatrix = ones(stepDurSR, 1) .* stepDownAmpV;

    % store stim matrix
    iInjOut = [spaceMatrix ;spaceMatrix; holdMatrix; spaceMatrix];

end