% measurePipetteResistance.m
%
% Function to calculate pipette or seal resistance by acquiring voltage
%  and current data when seal test is on.
% Will run for duration set by settings.sealTestDur
% Saves actual recordings in a folder in the pwd called preExptTrials
%
% INPUTS:
%   settings - struct of ephys setup settings, from ephysSettings()
%   type - string indicating whether we're measuring pipette or seal
%       resistance (options are either 'pipette' or 'seal')
%
% OUTPUTS:
%   pipetteResistance - pipette/seal resistance in MOhms
%
% CREATED: 11/4/19
% UPDATED
%   11/4/19 - HHY
%   2/27/20 - HHY
%   3/24/20 - HHY - debugging error with how troughs extracted
%

function pipetteResistance = measurePipetteResistance(settings, type)

    % some constants
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
    
    % compute resistance
    % borrows method from Yvette's code but averages over all pulses of
    %  seal test instead of just using the first
    
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
    
    % midpoint indices of pulse on
    pulseMids = round(pulseEnds - ((pulseEnds - pulseStarts)/2));
    
    % trough - when pulse off; calculate start, end, midpoint indicies
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
    
    % midpoint indices of pulse off
    troughMids = round(troughEnds - ((troughEnds - troughStarts)/2));
    
    % find current and voltage during pulses and troughs, use only the last
    %  half of each pulse/trough (from midpoint to endpoint)
    pulseCurrents = zeros(size(pulseStarts));
    pulseVoltages = pulseCurrents;
    troughCurrents = zeros(size(troughStarts));
    troughVoltages = troughCurrents;
    
    % loop through all pulses
    for i = 1:length(pulseStarts)
        pulseCurrents(i) = mean(...
            ephysData.current(pulseMids(i):pulseEnds(i)));
        pulseVoltages(i) = mean(...
            ephysData.voltage(pulseMids(i):pulseEnds(i)));
    end
    % loop through all troughs
    for i = 1:length(troughStarts)
        troughCurrents(i) = mean(...
            ephysData.current(troughMids(i):troughEnds(i)));
        troughVoltages(i) = mean(...
            ephysData.voltage(troughMids(i):troughEnds(i)));
    end
    
    % get mean current and voltage during pulse and trough
    meanPulseCurrent = mean(pulseCurrents);
    meanPulseVoltage = mean(pulseVoltages);
    meanTroughCurrent = mean(troughCurrents);
    meanTroughVoltage = mean(troughVoltages);
        
    % voltage and current differences b/w pulse and trough
    voltageDiff = meanPulseVoltage - meanTroughVoltage;
    currentDiff = meanPulseCurrent - meanTroughCurrent;
    
    % calculate pipette resistance using R = V/I, in MOhms
    pipetteResistance = ((voltageDiff * V_PER_mV) / ...
        (currentDiff * A_PER_pA)) * MOhm_PER_Ohm;
    
    % save trial in current directory (preExptTrials folder) in file named
    %  for type
    % note: overrides previous measurement for same cell (e.g. 1st pipette
    %  failed, try again with new pipette)
    
    switch type
        case 'pipette'
            filename = 'pipetteResistance.mat';
        case 'seal'
            filename = 'sealResistance.mat';
    end
    
    % save actual processed recording
    save(filename, 'ephysData', 'ephysMeta', '-v7.3');
    
end