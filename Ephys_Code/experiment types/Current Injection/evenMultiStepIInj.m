% evenMultiStepIInj.m
%
% Current Injection Function. Injects square current steps of user 
%  specified number of amplitudes between max and min user specified 
%  values. Has option to randomize amplitudes or not. User specified step
%  duration and time between steps. User specified space amplitude.
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
% CREATED: 7/23/20 - HHY
%
% UPDATED:
%   7/24/20 - HHY
%

function [iInjOut, iInjParams] = evenMultiStepIInj(settings, durScans)

    % prompt user for input parameters, as dialog box
    inputParams = {'Min Step Amplitude (pA):', ...
        'Max Step Amplitude (pA):', 'Num Steps:', 'Step Duration (s):', ...
        'Space Amplitude (pA):', 'Space Duration (s):', ...
        'Randomize steps? y/n'};
    dlgTitle = 'Enter parameter values';
    dlgDims = [1 35]; % dimensions of dialog box input fields
    
    % dialog box
    dlgAns = inputdlg(inputParams, dlgTitle, dlgDims);
    
    % convert user input into actual variables
    minStepAmp = str2double(dlgAns{1});
    maxStepAmp = str2double(dlgAns{2});
    numSteps = str2double(dlgAns{3});
    stepDur = str2double(dlgAns{4});
    spaceAmp = str2double(dlgAns{5});
    spaceDur = str2double(dlgAns{6});
    randomize = dlgAns{7};
    
    % standarize randomize; defaults to no if input isn't y or Y
    if strcmpi(randomize,'y')
        randomize = 'y';
    else
        randomize = 'n';
    end
    
    % get amplitudes of steps to present 
    allStepAmps = linspace(minStepAmp, maxStepAmp, numSteps);
    
    % convert user input into correct units for output (amplitude in volts,
    %  duration in scans); compensate for non-zero output from DAQ when
    %  zero commanded
    allStepAmpsV = (allStepAmps - settings.VOut.zeroI) .* ...
        settings.VOut.IConvFactor;
    stepDurScans = round(stepDur * settings.bob.sampRate);
    spaceAmpV = (spaceAmp - settings.VOut.zeroI) * ...
        settings.VOut.IConvFactor;
    spaceDurScans = round(spaceDur * settings.bob.sampRate);
    
    % save user input into parameters struct (convert duration to actual
    %  duration delivered, if rounded)
    iInjParams.allStepAmps = allStepAmps;
    iInjParams.stepDur = stepDurScans / settings.bob.sampRate;
    iInjParams.spaceAmp = spaceAmp;
    iInjParams.spaceDur = spaceDurScans / settings.bob.sampRate;
    iInjParams.randomize = randomize;
    
    % each stimulus as a column, in order from min to max; starts with
    %  space, ends with step
    spaceMatrix = ones(spaceDurScans, numSteps) * spaceAmpV;
    stepMatrix = ones(stepDurScans, numSteps) .* allStepAmpsV;
    stimMatrix = [spaceMatrix; stepMatrix];
      
    % generate stimulus, depends on whether or not to randomize
    if strcmpi(randomize, 'y')
        
        % number of full repeats (one step and space)
        numReps = floor(durScans / (spaceDurScans + stepDurScans));

        % number of scans in remainder
        numScansLeft = mod(durScans, (spaceDurScans + stepDurScans));
    
        % initialize iInjOut
        iInjOut = ones(durScans, 1);
        
        % generate random order of steps, index numbers into columns
        stepOrderInd = randi(numSteps, numReps+1, 1); 
       
        % loop and index into stimMatrix to pull appropriate stimuli
        for i = 1:numReps
            startInd = 1 + (spaceDurScans + stepDurScans) * (i - 1);
            endInd = (spaceDurScans + stepDurScans) * i;
            
            iInjOut(startInd:endInd) = stimMatrix(:, stepOrderInd(i));        
        end
        
        % if there is a remainder, add in appropriate part of last stimulus
        if numScansLeft
            iInjOut((end-numScansLeft+1):end) = ...
                stimMatrix(1:numScansLeft, stepOrderInd(end));
        end

    else % in order from min to max
        % number of full repeats of stimulus, all steps
        numFullReps = floor(durScans / ...
            ((spaceDurScans + stepDurScans)*numSteps));
        numScansLeft = mod(durScans, ...
            ((spaceDurScans + stepDurScans)*numSteps));
        
        % reshape stim matrix to single column vector
        oneRepStim = reshape(stimMatrix, numel(stimMatrix), 1);
        
        % generate output
        iInjOut = repmat(oneRepStim, numFullReps, 1);
        
        % if there is a remainder, concatenate onto end
        if numScansLeft
            iInjOut = [iInjOut; oneRepStim(1:numScansLeft)];
        end  
    end
end