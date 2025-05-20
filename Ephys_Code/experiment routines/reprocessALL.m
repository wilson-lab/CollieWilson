%% one script to clear them all
% initialize
clear
close all

%% one script to find them
% locate all dates in this data directory
dataFolder = 'G:\Behavior Rig\Gain\';
cd(dataFolder)
allFolders = dir('20*');
nFolders = length(allFolders);

%% one script to bring them all, and in a series of for loops re-process them
parfor f = 1:nFolders
    disp(['PROCESSING ' num2str(f) '/' num2str(nFolders) '...'])
    thisDate = fullfile(dataFolder,allFolders(f).name);
    cd(thisDate)
    allFlies = dir('fly*');
    nFlies = length(allFlies);
    for e = 1:nFlies
        thisFly = fullfile(thisDate,allFlies(e).name);
        cd(thisFly)
        allCells = dir('cell*');
        allCells(contains(string({allCells.name}), 'processed')) = [];
        nCells = length(allCells);
        for c = 1:nCells
            thisDataset = fullfile(thisFly,allCells(c).name);
            reprocessExpt(thisDataset)
        end

    end
    
end
disp('Complete!')