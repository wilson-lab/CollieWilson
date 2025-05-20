% preExptRoutine.m
%
% Function that measures pipette resistance, seal resistance, cell attached
%  spikes, access resistance, and resting Vm as the early part of the 
%  patching process. Gives option to run or not run each of these.
%
% INPUT:
%   settings - struct of ephys setup settings, from ephysSettings()
%
% OUTPUT:
%   preExptData - struct of pre-experimental data, with fields:
%       pipetteResistance
%       sealResistance
%       initialHoldingCurrent
%       initialAccessResistance
%       initialInputResistance
%       initialRestingVoltage
%       internal
%
% Adapted:  01/30/2021 - HY
% Updated:  04/02/2021 - MC added cell attached and i=0 processing/plotting
%           08/16/2021 - MC added mode reminders
%           11/01/2021 - MC added cell-attached/i=0 plot save
%           12/07/2021 - MC added spikerate

function preExptData = preExptRoutine(settings)

    preExptFolderName = 'preExptTrials';
    
    % make preExptTrials folder if it does't already exist
    if ~isfolder(preExptFolderName)
        mkdir(preExptFolderName);
    end
    
    % go to preExptFolder
    cd(preExptFolderName);
    
    preExptPath = pwd;
    
    % set sample durations
    duration = 20;
    
    %% Ask about pipette internal
    intSln = input('\nWhich internal? ', 's');
    preExptData.internal = intSln;
    
    %% Measure pipette resistance loops until 'no' selected
    contA = input('\nMeasure pipette resistance? (y/n) ','s');
    
    while 1
        if strcmpi(contA,'y') || strcmpi(contA,'') % 'y' or enter
            type = 'pipette';
            preExptData.pipetteResistance = measurePipetteResistance(settings, type);
            
            printVariable(preExptData.pipetteResistance,'Pipette Resistance', ' MOhms');
            
            contA = input(...
                '\nMeasure pipette resistance AGAIN? (y/n) ', 's');
            if  strcmpi(contA,'n')
                break;
            end
        else
            break
        end
    end

    %% Measure seal resistance
    contAns = input('\nMeasure seal resistance? (y/n) ','s');
    
    if ~strcmpi(contAns,'n')
        type = 'seal';
        preExptData.sealResistance = measurePipetteResistance(settings, type);

        % function returns MOhms, divide by 1000 to report GOhms
        printVariable(preExptData.sealResistance/1000 , 'Seal Resistance', ' GOhms');
    end

    %% Measure voltage trial to look at cell attached spikes
    contAns = input('\nRun cell-attached trial? (y/n)','s');
    
    if ~strcmpi(contAns,'n')
        
        [rawData, inputParams, rawOutput] = recordEphysOnly(settings, duration);
        
        % process recording snippet
        [daqData, daqOutput, daqTime] = preprocessUserDaq(inputParams, rawData, rawOutput, settings);
        [exptData, exptMeta] = processExptData(daqData, daqOutput, daqTime, inputParams, settings);
        % low pass filter
        %exptData.scaledCurrent = lowpass(exptData.scaledCurrent,1e3,settings.bob.sampRate);
        
        % find spike rate
        [~,~,spikeRate] = convolveSpikeRate(settings,exptData,'gaussian');
        vc_avgspikerate = nanmean(spikeRate);
        disp(['Spike Rate: ' num2str(round(vc_avgspikerate)) '/s']);
        % Save cellAttachedTrial trial in preExptTrials folder
        save('cellAttachedTrial.mat', 'exptData', 'exptMeta','inputParams','spikeRate', '-v7.3');
        
        % save plot
        sgtitle([datestr(now,'yyyy mm dd') ' cell attached'])
        saveas(gcf,[datestr(now,'yyyy_mm_dd') '_cellattached_plot.png']);
        
    end


    %% Measure access and input resistance and holding current
    contAns = input('\nMeasure access resistance? (y/n) ','s');
    
    if ~strcmpi(contAns,'n')
        [preExptData.initialHoldingCurrent, ...
            preExptData.initialAccessResistance, ...
            preExptData.initialInputResistance] = ...
            measureAccessResistance(settings);

        printVariable(preExptData.initialHoldingCurrent, 'Holding Current', ' pA');
        printVariable(preExptData.initialAccessResistance, 'Access Resistance', ' MOhms');
        printVariable(preExptData.initialInputResistance, 'Input Resistance', ' MOhms');

    end

    %% Measure resting voltage (I = 0)
    contAns = input('\nRun I=0 trial? (y/n)','s');
    
    if ~strcmpi(contAns,'n')
        
        % acquire trial
        [rawData, inputParams, rawOutput] = recordEphysOnly(settings,duration);
        
        % process recording snippet
        [daqData, daqOutput, daqTime] = preprocessUserDaq(inputParams, rawData, rawOutput, settings);
        [exptData, exptMeta] = processExptData(daqData, daqOutput, daqTime, inputParams, settings);

        % find resting voltage
        preExptData.initialRestingVoltage = mean(exptData.voltage);
        printVariable(preExptData.initialRestingVoltage,'Resting Voltage', 'mV');
        
        % find spike rate
        [~,~,spikeRate] = convolveSpikeRate(settings,exptData,'gaussian');
        vc_avgspikerate = nanmean(spikeRate);
        disp(['Spike Rate: ' num2str(round(vc_avgspikerate)) '/s']);

        % Save resting voltage trial
        save('restingVoltageTrial.mat', 'exptData', 'exptMeta', 'inputParams', 'spikeRate','-v7.3');
        
        sgtitle([datestr(now,'yyyy mm dd') ' i=0'])
        saveas(gcf,[datestr(now,'yyyy_mm_dd') '_restingV_plot.png']);
        
    end


    %% Check if preExptData was created
    if ~exist( 'preExptData', 'var')
        disp('WARNING: preExptData varible is empty!!');
        % create an empty struct as a place holder
        preExptData = struct;
    end
    
    cd .. % go back to previous folder, not preExptFolder
end

%% Helper Functions
function printVariable(value, label, unit)
    fprintf(['\n' label, ' = ', num2str(value), unit]);
end


