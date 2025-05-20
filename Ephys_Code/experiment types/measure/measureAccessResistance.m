% measureAccessResistance.m
%
% Function that acquires a trial of voltage clamp data with the seal test
%  on and uses it to calculate the access resistance, input resistance, 
%  and the holding current.
% Saves the recording in the present working directory
%
% INPUT:
%   settings - struct of ephys setup settings, from ephysSettings()
%
% OUTPUT:
%   holdingCurrent - average current being injected to hold voltage (pA)
%   accessResistance - deltaV / delta transient current (MOhms)
%   inputResistance - deltaV / delta steady state current (MOhms) 
%
% CREATED: 2/27/20 - HHY
% UPDATED: 02/09/2021 - MC
%
function [holdingCurrent, accessResistance, inputResistance] = ...
    measureAccessResistance(settings)

    % some constants
    startSteadyState = 2/3; % steady state as 2/3 from end of pulse
    endSteadyState = 1/3; % to 1/3 from end of pulse
    VOLTAGE_STEP_AMP = 5; %mV  (seal test from the amplifier)
    V_PER_mV = 1e-3; % V /1000 mV
    A_PER_pA = 1e-12; % 1e-12 A / 1 pA
    MOhm_PER_Ohm = 1e-6; % 1 MOhm / 1e6 Ohm
    
    
    % acquire recording 
    % use duration specified in ephys settings
    duration = settings.sealTestDur; 
    [rawData, inputParams, rawOutput] = recordEphysOnly(settings, duration);
    
    % process recording snippet
    [daqData, daqOutput, daqTime] = preprocessUserDaq(inputParams, ...
        rawData, rawOutput, settings);
    [ephysData, ephysMeta] = preprocessEphysData(daqData, daqOutput, ...
        daqTime, inputParams, settings);

    % holding current as mean current during trial
    holdingCurrent = mean(ephysData.current);
    
    % logical array for when voltage is above the mean
    highV = ephysData.voltage > mean(ephysData.voltage);
    
    % steps up from 0 to 1 (0.8 as buffer) to find all starts of pulse; +1
    %  to get start index
    pulseStarts = find(diff(highV) > 0.8) + 1;
    % steps down from 1 to 0 (-0.8 as buffer) to find all ends of pulse;
    pulseEnds = find(diff(highV) < -0.8); 
    
    % clip ends that don't have corresponding start (at beginning of trace)
    if pulseEnds(1) < pulseStarts(1)
        pulseEnds = pulseEnds(2:end);
    end
    % clip starts that don't have corresponding end (at end of trace)
    if pulseStarts(end) > pulseEnds(end)
        pulseStarts = pulseStarts(1:end-1);
    end
    
    % trough - when pulse off; calculate start, end
    troughStarts = pulseEnds + 1;
    troughEnds = pulseStarts - 1;
    
    % for troughs, clip ends that don't have a corresponding start
    if troughEnds(1) < troughStarts(1)
        troughEnds = troughEnds(2:end);
    end
    % for troughs, clip starts that don't have corresponding end
    if troughStarts(end) > troughEnds(end)
        troughStarts = troughStarts(1:end-1);
    end
    
    % pulse steady state start and ind 
    pulseSSStarts = round(pulseEnds - ((pulseEnds - pulseStarts) * ...
        startSteadyState));
    pulseSSEnds = round(pulseEnds - ((pulseEnds - pulseStarts) * ...
        endSteadyState));
    
    % trough steady state start and end ind
    troughSSStarts = round(troughEnds - ((troughEnds - troughStarts) * ...
        startSteadyState));
    troughSSEnds = round(troughEnds - ((troughEnds - troughStarts) * ...
        endSteadyState));
    
    % find steadyState current during pulses, troughs using only points b/w
    %  startSteadyState and endSteadyState 
    % also find peak current during pulses, only between pulse start and
    %  startSteadyState
    pulseSSCurrents = zeros(size(pulseSSStarts));
    peakCurrents = pulseSSCurrents;
    troughSSCurrents = zeros(size(troughSSStarts));

    % loop through all pulses
    for i = 1:length(pulseSSStarts)
        pulseSSCurrents(i) = mean(...
            ephysData.current(pulseSSStarts(i):pulseSSEnds(i)));
        peakCurrents(i) = max(ephysData.current(...
            pulseStarts(i):pulseSSStarts(i)));
    end
    % loop through all troughs
    for i = 1:length(troughSSStarts)
        troughSSCurrents(i) = mean(...
            ephysData.current(troughSSStarts(i):troughSSEnds(i)));
    end

    % get mean steady state current, peak current - relative to mean steady
    %  state current during trough
    meanPeakCurrent = mean(peakCurrents) - mean(troughSSCurrents);
    meanSSCurrent = mean(pulseSSCurrents) - mean(troughSSCurrents);
    
    % solve for access resistance in MOhm (R = V/I) using peak current
    accessResistance = ((VOLTAGE_STEP_AMP * V_PER_mV) / ...
        (meanPeakCurrent * A_PER_pA)) * MOhm_PER_Ohm;
    
    % solve for input resistance in MOhm using steady state current
    inputResistance = ((VOLTAGE_STEP_AMP * V_PER_mV) / ...
        (meanSSCurrent * A_PER_pA)) * MOhm_PER_Ohm;

    % save trace in current directory (preExptTrials)
    save('accessResistance.mat', 'ephysData', 'ephysMeta', '-v7.3');

end