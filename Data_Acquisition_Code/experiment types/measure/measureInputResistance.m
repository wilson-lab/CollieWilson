% measureInputResistance.m
%
% Function that measures input resistance between trials by applying 
% a current step and recording the resulting voltage change.
%
% INPUT:
% - settings  : Struct containing electrophysiological setup settings, 
%               typically obtained from the ephysSettings() function.
%
% OUTPUT:
% - inputResistance : Calculated input resistance (in MOhms) based on the 
%                     change in voltage and current during the current 
%                     injection protocol.
%
% ORIGINAL: 04/06/2021 by MC
%
function inputResistance = measureInputResistance(settings)
    
    %% Measure input resistance    
    
    % which input and output data streams used in this experiment
    inputParams.aInCh = {'ampI', 'amp10Vm', 'ampScaledOut', 'ampMode', 'ampGain', 'ampFreq'};
    inputParams.aOutCh = {'ampExtCmdIn'};
    inputParams.dInCh = {};
    inputParams.dOutCh = {};

    % initialize DAQ, including channels
    [userDAQ, ~, ~, ~, ~] = initUserDAQ(settings, inputParams.aInCh, ...
        inputParams.aOutCh, inputParams.dInCh, inputParams.dOutCh);
    
    
    % set output parameters
    spaceDur = 1;     % spacer duration (s)
    stepDur = 2;   % step duration (s)
    totDur = spaceDur+stepDur;
    
    stepAmp = -5;  % step amplitude (pA)
    spaceAmp = 0;   % amplitude b/n steps (pA)

    
    % convert user input into correct units for output (amplitude in volts,
    %  duration in scans);
    %  compensate for non-zero output from DAQ when zero commanded
    spaceDurSR = round(spaceDur * settings.bob.sampRate);
    stepDurSR = round(stepDur * settings.bob.sampRate);

    spaceAmpV = (spaceAmp - settings.VOut.zeroI) * settings.VOut.IConvFactor;
    stepAmpV = (stepAmp - settings.VOut.zeroI) * settings.VOut.IConvFactor;
    
    % each stimulus as a column, starts with space, ends with step
    spaceMatrix = ones(spaceDurSR, 1) * spaceAmpV;
    stepMatrix = ones(stepDurSR, 1) .* stepAmpV;
    % generate stimulus output
    rawOutput = [spaceMatrix; stepMatrix];
    
    % queue current injection output
    userDAQ.queueOutputData(rawOutput);
    
    disp('Acquiring ephys trial for input resistance measurement');
    % acquire data (in foreground)
    rawData = userDAQ.startForeground();
    
    % to stop it from presenting non-zero values if current injection
    %  protocol ends on non-zero value
    userDAQ.outputSingleScan(0);
    
    %% Calculate input resistance
    
    % conversions
    V_PER_mV = 1e-3; % V /1000 mV
    A_PER_pA = 1e-12; % 1e-12 A / 1 pA
    MOhm_PER_Ohm = 1e-6; % 1 MOhm / 1e6 Ohm
    
    % decode telegraphed output
    gain = decodeTelegraphedOutput(rawData(1:10,5), 'gain');
        
    % quickly process the data
    softGain = 1000/gain;
    voltage = softGain .* rawData(:,3);
    current = settings.I.softGain .* rawData(:,1);
    Q = length(voltage)/totDur;
    
    % calculate change in voltage/current by subtracting end of each step
    voltageDiff = mean(voltage(Q*2:end)) - mean(voltage(5000:Q));
    currentDiff = mean(current(Q*2:end)) - mean(current(5000:Q));
    
    % calculate input resistance accordingly
    inputResistance = ((voltageDiff* V_PER_mV)/(currentDiff*A_PER_pA)) * MOhm_PER_Ohm;
    


