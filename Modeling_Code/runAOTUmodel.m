% runAOTUmodel
%
% Main script for running multiple iterations of a simplified AOTU (anterior optic tubercle) model 
% to evaluate the impact of varying parameters such as synapse weight, synapse sign (inhibitory vs excitatory),
% noise levels, open-loop vs closed-loop conditions, and feedforward modulation on steering behavior.
%
% Created: 10/05/2024 - MC
% Updated: 11/01/2024 - MC
%
%% Initialize workspace and variables
% Clear all variables and close any open figures to ensure a clean workspace
clear
close all

% Initialize folder paths, plot settings, and run-specific parameters by calling the modelSettings function.
[folder, plotSettings, runSettings] = modelSettings();

% Change the directory to the script's file path to access data files.
cd(folder.filePath)

% Load the receptive field (RF) data from the 'Pursuit_RFs.mat' file, which contains information about
% the RFs of AOTU neurons for visual pursuit behavior.
load("Pursuit_RFs.mat");

% Process the raw RF data using the 'processRFdata' function to generate predicted receptive fields (predicted_RF).
% These predicted RFs will be used in subsequent model performance comparisons.
predicted_RF = processRFdata(Pursuit_RFs);

% Plot RFs and ELU for reference
plotModelSettings(predicted_RF, runSettings)
% Save the plot
saveas(gcf, fullfile(folder.settings, ['Model_Settings' '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, ['Model_Settings' '.svg']));

%% Run AOTU019 direction selectivity comparison model
restricted_RF = predicted_RF;
restricted_RF{:,1} = zeros(height(restricted_RF), 1);

modelPerformance_trajectory(restricted_RF, 'dirselective');
modelPerformance_step(restricted_RF, 'dirselective')

%% Run AOTU019 elimination comparison model
modelPerformance(predicted_RF, 'strength');  % Specify 'dirselective' as the comparison type.