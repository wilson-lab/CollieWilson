% load_kir_validation
% Master script for processing electrophysiology experiments comparing 
% basic electrophysiological parameters for cells expressing or not expressing 
% KIR2.1. This script manages the initialization, data loading, and generation 
% of summary plots for the validation of KIR expression in AOTU019 neurons.
%
% Process Overview:
%   1. Loads all experiments of interest from a designated Excel tracker file.
%   2. Extracts the relevant data for each fly, including holding current, 
%      access resistance, input resistance, resting voltage, and firing rate.
%   3. Plots the basic electrophysiological parameters for KIR2.1-expressing 
%      versus non-expressing cells (GFP controls).
%   4. Saves the final plot as both PNG and vectorized SVG for publication-quality figures.
%
% Outputs:
%   - A summary plot is generated comparing key electrophysiological parameters 
%     between KIR2.1-expressing cells and controls. The plot is saved in the 
%     specified experiment folder.
%
function load_kir_validation()
%% load in all experiments of interest
% clean up
clear
close all

% designate folder to save all processed data in
exptFolder = 'AOTU019 KIR Openloop';
exc_name = 'KIR_Behavior_Tracker.xlsx'; %excel file containing meta data
exc_sheet = 3; %sheet for this experiment

dataFolder = 'E:\';
cd(dataFolder)
exptFolder = fullfile(dataFolder,exptFolder,'validation');
mkdir(exptFolder)

% load in tracker for this experiment set
tracker = readtable(exc_name,'Sheet',exc_sheet); %load all expts
expt_tracker = tracker(tracker.Include==1,:); %select only included expts
nFly = size(expt_tracker,1);

%% generate summary plot

% initialize
close all
figure;  set(gcf,'Position',[100 100 500 900])
sx = 3;
sy = 3;

ljp = 13; %mv, liquid junction potential

for f = 1:nFly
    % fetch genotype and set plot color
    thisGenotype = tracker.Genotype(f);
    if strcmp(thisGenotype,'KIR')
        x = 2;
        pointcolor = "#77AC30"; %green
    else
        x = 1;
        pointcolor = [0.2 0.2 0.2]; %dark grey
    end
    % fetch filepaths
    thisfly = fullfile(dataFolder,tracker.Date(f),sprintf('fly%02d',expt_tracker.Fly(f)),'cell01');
    preData = fullfile(thisfly,'preExptData.mat');
    i0Data = fullfile(thisfly,'preExptTrials','restingVoltageTrial.mat');
    % load data from this fly
    load(preData{1})
    load(i0Data{1})
    
    % plot break in data
    subplot(sx,sy,1); hold on
    scatter(x,preExptData.initialHoldingCurrent,'filled','MarkerEdgeColor','none','MarkerFaceColor',pointcolor)
    xlim([0 3]); xticks([1,2]); xticklabels({'GFP','KIR'});
    yline(0); ylabel('VClamp (@-40mV) Holding Current (pA)')

    subplot(sx,sy,2); hold on
    scatter(x,preExptData.initialAccessResistance,'filled','MarkerEdgeColor','none','MarkerFaceColor',pointcolor)
    xlim([0 3]); xticks([1,2]); xticklabels({'GFP','KIR'});
    ylim([0 26]); ylabel('Access Resistance (mOhm)')

    subplot(sx,sy,3); hold on
    scatter(x,preExptData.initialInputResistance,'filled','MarkerEdgeColor','none','MarkerFaceColor',pointcolor)
    xlim([0 3]); xticks([1,2]); xticklabels({'GFP','KIR'});
    ylim([50 90]); ylabel('Input Resistance (mOhm)')
    
    % plot i=0 data
    subplot(sx,sy,4); hold on
    restingVoltage = preExptData.initialRestingVoltage-ljp;
    scatter(x,restingVoltage,'filled','MarkerEdgeColor','none','MarkerFaceColor',pointcolor)
    xlim([0 3]); xticks([1,2]); xticklabels({'GFP','KIR'});
    ylim([-100 -40]); ylabel('i=0 Resting Voltage (mV)')

    subplot(sx,sy,5); hold on
    restingFiringRate = median(spikeRate);
    scatter(x,restingFiringRate,'filled','MarkerEdgeColor','none','MarkerFaceColor',pointcolor)
    xlim([0 3]); xticks([1,2]); xticklabels({'GFP','KIR'});
    ylim([0 3]); ylabel('i=0 Firing Rate (spikes/s)')

    % plot iclamp data
    if ~isnan(tracker.vThresh(f))
            subplot(sx,sy,7); hold on
            spikeVoltage = tracker.vThresh(f)-ljp;
            scatter(x,spikeVoltage,'filled','MarkerEdgeColor','none','MarkerFaceColor',pointcolor)
            xlim([0 3]); xticks([1,2]); xticklabels({'GFP','KIR'});
            yline(0); ylim([-100 200]); ylabel('Spike Threshold Voltage (mV)')

            subplot(sx,sy,8); hold on
            scatter(x,tracker.iThresh(f),'filled','MarkerEdgeColor','none','MarkerFaceColor',pointcolor)
            xlim([0 3]); xticks([1,2]); xticklabels({'GFP','KIR'});
            yline(0); ylim([-100 2000]); ylabel('Spike Threshold Current (pA)')
    else
        subplot(sx,sy,7); hold on
        spikeVoltage = 200;
        scatter(x,spikeVoltage,'x','MarkerEdgeColor',pointcolor)

        subplot(sx,sy,8); hold on
        scatter(x,2000,'x','MarkerEdgeColor',pointcolor)
    end
end

% save plot
sgtitle('AOTU019 KIR Validation')
cd(exptFolder)
plotname = 'validate_kirAOTU019';
saveas(gcf,[plotname '.png'])
% save vectorized plot
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
end