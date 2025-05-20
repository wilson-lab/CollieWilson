% reprocessExpt.m
%
% Function that retrieves all raw data files from a specified experiment 
% folder and re-processes them. This function is useful for applying 
% changes to the post-experiment routine across multiple files when needed.
%
% INPUTS:
% - exptFolder : Location of the raw data files to be reprocessed.
%
% ORIGINAL: 11/11/2021 by MC
%
function reprocessExpt(exptFolder)

% load ephys settings
[~, ~, settings] = ephysSettings();

disp('Re-processing data...')
% pull all raw files from the experiment folder of interest
cd(exptFolder)
close all
rawFiles = dir('*raw.mat');
rawFiles(contains(string({rawFiles.name}), 'Acclimate')) = [];
rawFiles(contains(string({rawFiles.name}), 'acclimate')) = [];

% for each raw file
for e = 1:length(rawFiles)
    % ensure folder is correct
    cd(exptFolder)
    % load in the file
    load(rawFiles(e).name)
    
    % re-run post experiment routine
    [~] = postExptRoutine(inputParams, rawData, rawOutput, settings);
    
end

copyfile('metaDat.mat', [exptFolder '_processed'],'f');

disp('Data re-processed successfully')
end

